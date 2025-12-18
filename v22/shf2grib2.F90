PROGRAM CONVERT_SHF_TO_GRIB2

! Program reads in input a shf file: short model (Moloch) history file
! and writes in output data in grib2 format using the ECMWF GRIB_API library
! version for .shf file with multiple instants

!USE grib_api
!USE GRIB2_CODING_DATA
USE eccodes

IMPLICIT NONE

INTEGER :: INPUT_UNIT=54, IGRIBSAMPLE1, IGRIBSAMPLE, IGRIBOUT, IGRIBCLONE, IERR
CHARACTER(LEN=30) :: INPUT_FILE="moloch.shf"
CHARACTER :: MODEL_NAME*10

INTEGER, DIMENSION(50) :: NFDR
REAL, DIMENSION(200) :: PDR

INTEGER :: NX, NY, NLEV, NSTEP_PERIOD, NDAY_FC, NHOUR_FC, NMIN_FC, NDAY, NDM, IIMON, &
 I, J, K, IFIELD, N_VERT_COORD_PAR, NMEM, IMEM, nmsgg2
REAL :: PI, DZ, DTSTEP, DELTA_FORECAST, DELTA_PERIOD
INTEGER, DIMENSION(12) ::IMON=(/31,28,31,30,31,30,31,31,30,31,30,31/)
REAL*4, DIMENSION(:), ALLOCATABLE :: VERT_COORD_PAR

! Complete pass and name of grib sample file

 OPEN (INPUT_UNIT, FILE="grib_sample.inp", STATUS="OLD")
 READ (INPUT_UNIT, GRIB_SAMPLE)
 CLOSE (INPUT_UNIT)

! For ensemble forecasts only

! PRINT *
! PRINT *,' Enter total number of ensemble forecast members'
! READ (*,*) NMEM
! PRINT *,' You entered ',NMEM,' of ensemble forecast members'
!
! PRINT *
! PRINT *,' Enter index of current ensemble forecast member'
! READ (*,*) IMEM
! PRINT *,' You entered ',IMEM,' index of current ensemble forecast member'

! Opening of input file

 OPEN (INPUT_UNIT, FILE=INPUT_FILE, FORM="UNFORMATTED", STATUS="OLD")

! Reading of input file header

 READ (INPUT_UNIT) NFDR
 READ (INPUT_UNIT) PDR

! Grid parameters

 DATA(1) % NX = NFDR(2)
 DATA(1) % NY = NFDR(3)
 NX = DATA(1) % NX
 NY = DATA(1) % NY
 NLEV = NFDR(4)
 DATA(1) % X0 = PDR(39)
 DATA(1) % Y0 = PDR(38)
 DATA(1) % DX = PDR(2)
 DATA(1) % DY = PDR(1)
 DATA(1) % X00 = PDR(5)
 DATA(1) % Y00 = PDR(4) + DATA(1) % DY*0.5

! Initial date and time

 DATA(1) % IDATE0(1) = NFDR(5)
 DATA(1) % IDATE0(2) = NFDR(6)
 DATA(1) % IDATE0(3) = NFDR(7)
 DATA(1) % IDATE0(4) = NFDR(8)
 DATA(1) % IDATE0(5) = NFDR(9)
 NDAY_FC  = NFDR(10)
 NHOUR_FC = NFDR(11)
 NMIN_FC  = NFDR(12)
 DTSTEP = PDR(3)
 NSTEP_PERIOD = NFDR(18)

! Length of forecast term and period (for statistical data fields)

 DELTA_FORECAST = FLOAT(NDAY_FC)*86400. + FLOAT(NHOUR_FC)*3600. + FLOAT(NMIN_FC)*60.
 DELTA_PERIOD = FLOAT(NSTEP_PERIOD)*DTSTEP

 IF (DELTA_FORECAST < 300.) THEN
   DELTA_FORECAST=0.
   NMIN_FC=0
 ENDIF

! Date and time of current forecast term

 DATA(1) % IDATEC(:)=DATA(1) % IDATE0(:)
 DATA(1) % IDATEC(3)=DATA(1) % IDATEC(3)+NDAY_FC
 DATA(1) % IDATEC(4)=DATA(1) % IDATEC(4)+NHOUR_FC
 DATA(1) % IDATEC(5)=DATA(1) % IDATEC(5)+NMIN_FC
 IF (MOD(DATA(1) % IDATEC(1),4)==0) THEN
   IMON(2)=29
 ELSE
   IMON(2)=28
 ENDIF
 NDAY=NDAY_FC
 NDM=0
 DO IIMON=1,50
   NDAY=NDAY-IMON(DATA(1) % IDATEC(2))
   IF (NDAY<0) THEN
     EXIT
   ELSE
     NDM=NDM+IMON(DATA(1) % IDATEC(2))
     DATA(1) % IDATEC(2)=DATA(1) % IDATEC(2)+1
     IF (DATA(1) % IDATEC(2)>12) THEN
       DATA(1) % IDATEC(2)=DATA(1) % IDATEC(2)-12
       DATA(1) % IDATEC(1)=DATA(1) % IDATEC(1)+1
       IF (MOD(DATA(1) % IDATEC(1),4)==0) THEN
         IMON(2)=29
       ELSE
         IMON(2)=28
       ENDIF
     ENDIF
   ENDIF
 ENDDO
 DATA(1) % IDATEC(3)=DATA(1) % IDATEC(3)-NDM
 IF (DATA(1) % IDATEC(5)>=60) THEN
   DATA(1) % IDATEC(5)=DATA(1) % IDATEC(5)-60
   DATA(1) % IDATEC(4)=DATA(1) % IDATEC(4)+1
 ENDIF
 IF (DATA(1) % IDATEC(4)>=24) THEN
   DATA(1) % IDATEC(4)=DATA(1) % IDATEC(4)-24
   DATA(1) % IDATEC(3)=DATA(1) % IDATEC(3)+1
 ENDIF
 IF (DATA(1) % IDATEC(3)>IMON(DATA(1) % IDATEC(2))) THEN
   DATA(1) % IDATEC(3)=DATA(1) % IDATEC(3)-IMON(DATA(1) % IDATEC(2))
   DATA(1) % IDATEC(2)=DATA(1) % IDATEC(2)+1
   IF (DATA(1) % IDATEC(2)>12) THEN
     DATA(1) % IDATEC(2)=DATA(1) % IDATEC(2)-12
     DATA(1) % IDATEC(1)=DATA(1) % IDATEC(1)+1
   ENDIF
 ENDIF

 READ (INPUT_UNIT) DATA(1) % GRIB2_DESCRIPT

 REWIND (INPUT_UNIT)

! loop over instants...

100 continue

 READ (INPUT_UNIT, iostat=ierr) NFDR
 IF (IERR /= 0) stop 0
 nmsgg2 = nfdr(19) !! number of expected grib messages
 READ (INPUT_UNIT) PDR

 NDAY_FC  = NFDR(10)
 NHOUR_FC = NFDR(11)
 NMIN_FC  = NFDR(12)

! Date and time of current forecast term

 DATA(1) % IDATEC(:)=DATA(1) % IDATE0(:)
 DATA(1) % IDATEC(3)=DATA(1) % IDATEC(3)+NDAY_FC
 DATA(1) % IDATEC(4)=DATA(1) % IDATEC(4)+NHOUR_FC
 DATA(1) % IDATEC(5)=DATA(1) % IDATEC(5)+NMIN_FC
 IF (MOD(DATA(1) % IDATEC(1),4)==0) THEN
   IMON(2)=29
 ELSE
   IMON(2)=28
 ENDIF
 NDAY=NDAY_FC
 NDM=0
 DO IIMON=1,50
   NDAY=NDAY-IMON(DATA(1) % IDATEC(2))
   IF (NDAY<0) THEN
     EXIT
   ELSE
     NDM=NDM+IMON(DATA(1) % IDATEC(2))
     DATA(1) % IDATEC(2)=DATA(1) % IDATEC(2)+1
     IF (DATA(1) % IDATEC(2)>12) THEN
       DATA(1) % IDATEC(2)=DATA(1) % IDATEC(2)-12
       DATA(1) % IDATEC(1)=DATA(1) % IDATEC(1)+1
       IF (MOD(DATA(1) % IDATEC(1),4)==0) THEN
         IMON(2)=29
       ELSE
         IMON(2)=28
       ENDIF
     ENDIF
   ENDIF
 ENDDO
 DATA(1) % IDATEC(3)=DATA(1) % IDATEC(3)-NDM
 IF (DATA(1) % IDATEC(5)>=60) THEN
   DATA(1) % IDATEC(5)=DATA(1) % IDATEC(5)-60
   DATA(1) % IDATEC(4)=DATA(1) % IDATEC(4)+1
 ENDIF
 IF (DATA(1) % IDATEC(4)>=24) THEN
   DATA(1) % IDATEC(4)=DATA(1) % IDATEC(4)-24
   DATA(1) % IDATEC(3)=DATA(1) % IDATEC(3)+1
 ENDIF
 IF (DATA(1) % IDATEC(3)>IMON(DATA(1) % IDATEC(2))) THEN
   DATA(1) % IDATEC(3)=DATA(1) % IDATEC(3)-IMON(DATA(1) % IDATEC(2))
   DATA(1) % IDATEC(2)=DATA(1) % IDATEC(2)+1
   IF (DATA(1) % IDATEC(2)>12) THEN
     DATA(1) % IDATEC(2)=DATA(1) % IDATEC(2)-12
     DATA(1) % IDATEC(1)=DATA(1) % IDATEC(1)+1
   ENDIF
 ENDIF

! Definition of forecast and period length in defined time unit

 IF (DATA(1) % GRIB2_DESCRIPT(4) == 0) THEN ! Forecast time step and Statistical period in minutes
   DATA(1) % IFORECAST = NINT(DELTA_FORECAST/60.)
   DATA(1) % IPERIOD = NINT(DELTA_PERIOD/60.)
 ENDIF

 IF (DATA(1) % GRIB2_DESCRIPT(4) == 1) THEN ! Forecast time step and Statistical period in hours
   DATA(1) % IFORECAST = NINT(DELTA_FORECAST/3600.)
   DATA(1) % IPERIOD = NINT(DELTA_PERIOD/3600.)
 ENDIF

 IF (DATA(1) % GRIB2_DESCRIPT(4) == 2) THEN ! Forecast time step and Statistical period in days
   DATA(1) % IFORECAST = NINT(DELTA_FORECAST/86400.)
   DATA(1) % IPERIOD = NINT(DELTA_PERIOD/86400.)
 ENDIF

! Model vertical coordinate parameters and model name

 IF (DATA(1) % GRIB2_DESCRIPT(1) == 1) THEN ! Model code, 1 - Bolam
   N_VERT_COORD_PAR = NLEV+2
   ALLOCATE (VERT_COORD_PAR (N_VERT_COORD_PAR) )
   VERT_COORD_PAR(1)=PDR(40)*0.5
   DO K=2,NLEV
     VERT_COORD_PAR(K)=(PDR(39+K)+PDR(39+K-1))*0.5 ! sigma of integer levels
   ENDDO
   VERT_COORD_PAR(NLEV+1) = PDR(37) ! Alfa
   VERT_COORD_PAR(NLEV+2) = PDR(36) ! p0
   MODEL_NAME="bolam"
 ENDIF

 IF (DATA(1) % GRIB2_DESCRIPT(1) == 2) THEN ! Model code, 2 - Moloch
   N_VERT_COORD_PAR = NLEV+1
   ALLOCATE (VERT_COORD_PAR (N_VERT_COORD_PAR) )
   DZ = PDR(40)/FLOAT(NLEV) ! h/float(nlev)
   DO K=1,NLEV
     VERT_COORD_PAR(K) = 1.-(FLOAT(K-1)*DZ+DZ*0.5)/PDR(40) ! fz = 1.-zita/h
   ENDDO
   VERT_COORD_PAR(NLEV+1) = PDR(40) ! h
   MODEL_NAME="moloch"
 ENDIF

 IF (DATA(1) % GRIB2_DESCRIPT(1) == 3) THEN ! Model code, 3 - Globo
   N_VERT_COORD_PAR = NLEV+2
   ALLOCATE (VERT_COORD_PAR (N_VERT_COORD_PAR) )
   DO K=1,NLEV
     VERT_COORD_PAR(K)=(PDR(14+K)+PDR(14+K+1))*0.5 ! sigma of integer levels
   ENDDO
   VERT_COORD_PAR(NLEV+1) = PDR(2) ! Alfa
   VERT_COORD_PAR(NLEV+2) = PDR(1) ! p0
   MODEL_NAME="globo"
 ENDIF

! Reading data fields of input file

 DO IFIELD = 1, nmsgg2

   READ (INPUT_UNIT) DATA(IFIELD) % GRIB2_DESCRIPT

   ALLOCATE (DATA(IFIELD) % FIELD(NX,NY))
   DO J=1,NY
     READ (INPUT_UNIT) (DATA(IFIELD) % FIELD(I,J), I=1,NX)
   ENDDO

! Grid parameters

   DATA(IFIELD) % NX = DATA(1) % NX
   DATA(IFIELD) % NY = DATA(1) % NY
   DATA(IFIELD) % X0 = DATA(1) % X0
   DATA(IFIELD) % Y0 = DATA(1) % Y0
   DATA(IFIELD) % DX = DATA(1) % DX
   DATA(IFIELD) % DY = DATA(1) % DY
   DATA(IFIELD) % X00 = DATA(1) % X00
   DATA(IFIELD) % Y00 = DATA(1) % Y00

! Initial date and time

   DATA(IFIELD) % IDATE0(:) = DATA(1) % IDATE0(:)

! Length of forecast term and period (for statistical data fields)

   DATA(IFIELD) % IFORECAST = DATA(1) % IFORECAST
   DATA(IFIELD) % IPERIOD = DATA(1) % IPERIOD

! Date and time for current forecast term

   DATA(IFIELD) % IDATEC(:) = DATA(1) % IDATEC(:)

! If level type is Hybrid level, then vertical coordinate parameters are defined

   IF (DATA(IFIELD) % GRIB2_DESCRIPT(10) == 105 ) THEN ! level (layer) type, 105 - Hybrid level

     DATA(IFIELD) % N_VERT_COORD_PAR = N_VERT_COORD_PAR
     ALLOCATE ( DATA(IFIELD) % VERT_COORD_PAR (N_VERT_COORD_PAR) )
     DATA(IFIELD) % VERT_COORD_PAR(:) = VERT_COORD_PAR(:)

   ENDIF

 ENDDO

 IF (ALLOCATED(DATA(IFIELD) % FIELD)) DEALLOCATE (DATA(IFIELD) % FIELD)
 nfield = nmsgg2

! Output file name

 WRITE (OUTPUT_FILE_NAME,'(2A,I4.4,3I2.2,A,I3.3,2I2.2,A)') &
 TRIM(MODEL_NAME),"_",DATA(1) % IDATE0(1:4),"_",NDAY_FC,NHOUR_FC,NMIN_FC,".grib2.shf"

! Coding read data in grib2 format

 CALL WRITE_GRIB2_DATA

! Deallocation of all allocated arrays

 DO IFIELD = 1, nmsgg2
   IF (ALLOCATED(DATA(IFIELD) % FIELD) ) DEALLOCATE (DATA(IFIELD) % FIELD)
 ENDDO

 DEALLOCATE (VERT_COORD_PAR)

 go to 100

 STOP
 END
