!-----------------------------------------------------------------------------

SUBROUTINE WRITE_GRIB2_HORIZONTAL_GRID_STATIC_DATA(IMODEL,NLON,NLAT,DX,DY,X0,Y0,X1,Y1,NPOINT_JUMP, &
 IDATE0,IPERIOD_INP,IPERIOD_ACCUM,IDATEC,ZOROGR,FMASK,FLAG_LAT_LON)

! Procedure that prepares fields of model orography height and model land-sea mask
! for coding in grib2 format,
! with the option of grid point rarefied.

! Uses module grib2_coding_data in write_grib2_data.F90

USE GRIB2_CODING_DATA

IMPLICIT NONE

INTEGER :: IMODEL, NLON, NLAT, NPOINT_JUMP, NX, NY
REAL :: DX, DY, X0, Y0, X1, Y1
INTEGER, DIMENSION(5) :: IDATE0, IDATEC
INTEGER, DIMENSION(3) :: IPERIOD_INP, IPERIOD_ACCUM
REAL, DIMENSION(NLON,NLAT) :: ZOROGR, FMASK, LAT, LON
CHARACTER (LEN=10) :: MODEL_NAME
INTEGER :: IFIELD, IX1, IX2, IY1, IY2, FLAG_LAT_LON

! FLAG_LAT_LON 0 or 1: write field of Geographical latitude (°N) and
! Geographical longitude (°E) of rotated grid points

! Name of ouput file

  IF (IMODEL == 1) MODEL_NAME="bolam"
  IF (IMODEL == 2) MODEL_NAME="moloch"
  IF (IMODEL == 3) MODEL_NAME="globo"
  IF (IMODEL == 12) MODEL_NAME="blended"

  WRITE (OUTPUT_FILE_NAME,'(2A)') TRIM(MODEL_NAME),"_domain_orogr_lsm.grib2"

  NFIELD = 2
  IF (FLAG_LAT_LON == 1) NFIELD=NFIELD+2

! Cleaning of possible previous allocations

  DO IFIELD=1,NFIELD0
    IF (ALLOCATED(DATA(IFIELD) % FIELD) ) DEALLOCATE (DATA(IFIELD) % FIELD)
  ENDDO
  DO IFIELD=1,NFIELD
    DATA(IFIELD) % GRIB2_DESCRIPT(:) = IVALMISS
  ENDDO

! Description of data for grib2 coding (INTEGER, DIMENSION(50) :: GRIB2_DESCRIPT=IVALMISS)
! see module_write_grib2_data.F90

! Content of GRIB2_DESCRPT array - see bottom of file

! Model index: 1 - Bolam, 2 - Moloch, 3 - Globo

  DATA(1:NFIELD) % GRIB2_DESCRIPT(1) = IMODEL

! Status and type of data

  DATA(1:NFIELD) % GRIB2_DESCRIPT(6) = 0  ! operational products
  DATA(1:NFIELD) % GRIB2_DESCRIPT(7) = 1  ! forecast
  DATA(1:NFIELD) % GRIB2_DESCRIPT(30) = 0 ! bit-map absent

! Grid parameters

  DATA(1:NFIELD) % GRIB2_DESCRIPT(2) = 1 ! grid template index - horizontal grid
  NX = NLON/NPOINT_JUMP
  NY = NLAT/NPOINT_JUMP

  DATA(1:NFIELD) % NX = NX
  DATA(1:NFIELD) % NY = NY
  DATA(1:NFIELD) % X0 = X0
  DATA(1:NFIELD) % Y0 = Y0
  DATA(1:NFIELD) % DX = DX*FLOAT(NPOINT_JUMP)
  DATA(1:NFIELD) % DY = DY*FLOAT(NPOINT_JUMP)
  DATA(1:NFIELD) % X00 = X1
  DATA(1:NFIELD) % Y00 = Y1

! Initial date and time

  DO IFIELD=1,NFIELD

    DATA(IFIELD) % IDATE0(1:5) = IDATE0(1:5)

! Date and time of current forecast term

    DATA(IFIELD) % IDATEC(1:5) = IDATEC(1:5)

  ENDDO

! Definition of time unit, forecast and period length in defined time unit

  DATA(1:NFIELD) % GRIB2_DESCRIPT(4) = 1 ! Time unit for Forecast time step and Statistical period : 0 minute, 1 hour, 2 day
  DATA(1:NFIELD) % IFORECAST = IPERIOD_INP(1)*24+IPERIOD_INP(2)   ! Forecast time step
  DATA(1:NFIELD) % IPERIOD = IPERIOD_ACCUM(1)*24+IPERIOD_ACCUM(2) ! Statistical period

! Product template index

  DATA(1:NFIELD) % GRIB2_DESCRIPT(3) = 0 ! instant

! Level/layer parameters (for forecast satellite products parameters of satellite platform and sensor)

  DO IFIELD=1,NFIELD
    DATA(IFIELD) % GRIB2_DESCRIPT(10:15) = IVALMISS
  ENDDO

  DO IFIELD=1,NFIELD
    ALLOCATE(DATA(IFIELD) % FIELD(NX,NY))
  ENDDO

! Definition of data fields and data parameters in grib2 terms


  DATA(1:NFIELD) % GRIB2_DESCRIPT(10) = 1 ! Ground or water surface

! Orography

  IFIELD = 1
  DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2 ! Discipline: Land surface products
  DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0 ! Category: Vegetation/Biomass
  DATA(IFIELD) % GRIB2_DESCRIPT(22) = 7 ! Parameter: Model terrain height (m)
  IY2=0
  DO IY1=1,NLAT,NPOINT_JUMP
  IY2=IY2+1
  IX2=0
  DO IX1=1,NLON,NPOINT_JUMP
  IX2=IX2+1
    DATA(IFIELD) % FIELD(IX2,IY2)=ZOROGR(IX1,IY1)
  ENDDO
  ENDDO

! Land-Sea Mask

  IFIELD = 2
  DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2 ! Discipline: Land surface products
  DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0 ! Category: Vegetation/Biomass
  DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0 ! Parameter: Land cover (0=land, 1=sea) (proportion)
  IY2=0
  DO IY1=1,NLAT,NPOINT_JUMP
  IY2=IY2+1
  IX2=0
  DO IX1=1,NLON,NPOINT_JUMP
  IX2=IX2+1
    DATA(IFIELD) % FIELD(IX2,IY2)=1.-FMASK(IX1,IY1)
  ENDDO
  ENDDO

  IF (FLAG_LAT_LON == 1) THEN

    CALL ROT_GRID (X0, Y0, X1, Y1, DX, DY, LON, LAT, NLON, NLAT)

! Geographical latitude of grid points
! (see WMO table
! http://www.wmo.int/pages/prog/www/WMOCodes/WMO306_vI2/LatestVERSION/WMO306_vI2_GRIB2_CodeFlag_en.pdf)

    IFIELD = 3
    DATA(IFIELD) % GRIB2_DESCRIPT(20) =   0 ! Discipline: Meteorological products
    DATA(IFIELD) % GRIB2_DESCRIPT(21) = 191 ! Category:  miscellaneous
    DATA(IFIELD) % GRIB2_DESCRIPT(22) =   1 ! Parameter: Geographical latitude (°N)
    IY2=0
    DO IY1=1,NLAT,NPOINT_JUMP
    IY2=IY2+1
    IX2=0
    DO IX1=1,NLON,NPOINT_JUMP
    IX2=IX2+1
      DATA(IFIELD) % FIELD(IX2,IY2)=LAT(IX1,IY1)
    ENDDO
    ENDDO

! Geographical latitude of grid points
! (see WMO table
! http://www.wmo.int/pages/prog/www/WMOCodes/WMO306_vI2/LatestVERSION/WMO306_vI2_GRIB2_CodeFlag_en.pdf)

    IFIELD = 4
    DATA(IFIELD) % GRIB2_DESCRIPT(20) =   0 ! Discipline: Meteorological products
    DATA(IFIELD) % GRIB2_DESCRIPT(21) = 191 ! Category:  miscellaneous
    DATA(IFIELD) % GRIB2_DESCRIPT(22) =   2 ! Parameter: Geographical longitude (°E)
    IY2=0
    DO IY1=1,NLAT,NPOINT_JUMP
    IY2=IY2+1
    IX2=0
    DO IX1=1,NLON,NPOINT_JUMP
    IX2=IX2+1
      DATA(IFIELD) % FIELD(IX2,IY2)=LON(IX1,IY1)
    ENDDO
    ENDDO

  ENDIF

! End of preparation of output data for coding in grib2 format

! Write output data in grib2 format

  CALL WRITE_GRIB2_DATA

RETURN
END

!-----------------------------------------------------------------------------

#ifdef rttov

SUBROUTINE WRITE_GRIB2_HORIZONTAL_GRID_RTTOV_DATA(IMODEL, NLON, NLAT, NPOINT_JUMP,&
 X0, Y0, X1, Y1, DX, DY, &
 IDATE0, IPERIOD_INP, IPERIOD_ACCUM, IDATEC,&
 NCHAN, RTTOV_SAT_SERIES, RTTOV_SAT_ID, RTTOV_SAT_SENSOR,&
 SENSOR_CHAN_ID, SENSOR_CHAN_CW,&
 RADIANCE, RADIANCE_BT, RADIANCE_CLEAR, RADIANCE_CLEAR_BT, EMIS_RTTOV)

! Procedure that prepares fields of RTTOV model products for coding in grib2 format,
! with the option of grid point rarefied.

! Uses module grib2_coding_data in write_grib2_data.F90

USE GRIB2_CODING_DATA

IMPLICIT NONE

INTEGER :: IMODEL, NLON, NLAT, NPOINT_JUMP, NCHAN, RTTOV_SAT_SERIES, RTTOV_SAT_ID, RTTOV_SAT_SENSOR
REAL :: DX, DY, X0, Y0, X1, Y1
INTEGER, DIMENSION(5) :: IDATE0, IDATEC
INTEGER, DIMENSION(3) :: IPERIOD_INP, IPERIOD_ACCUM
CHARACTER (LEN=10) :: MODEL_NAME

INTEGER, DIMENSION(NCHAN) :: SENSOR_CHAN_ID
REAL*8, DIMENSION(NCHAN) :: SENSOR_CHAN_CW

REAL, DIMENSION(NLON,NLAT,NCHAN) :: RADIANCE, RADIANCE_BT, RADIANCE_CLEAR, RADIANCE_CLEAR_BT, EMIS_RTTOV

INTEGER :: NPARAM = 2
INTEGER :: NX, NY, IFIELD, ICHAN, IPARAM, IX1, IX2, IY1, IY2
INTEGER :: WMO_SAT_NUM, WMO_SAT_SER, WMO_INST_ID

INCLUDE 'convert_sat_number.F90'
INCLUDE 'convert_sat_series.F90'
INCLUDE 'convert_sat_sensor.F90'

! Name of ouput file

  IF (IMODEL == 1) MODEL_NAME="bolam"
  IF (IMODEL == 2) MODEL_NAME="moloch"
  IF (IMODEL == 3) MODEL_NAME="globo"

  WRITE (OUTPUT_FILE_NAME,'(2A,I4.4,3I2.2,A,I3.3,2I2.2,A)') TRIM(MODEL_NAME),"_rttov_prod_",IDATE0(1:4),"_",IPERIOD_INP(1:3),".grib2"

! Conversion of RTTOV satellite index into the satellite number in WMO classification
! accordingly to BUFR Code table 0 01 007 - missing value is 999

  IF (RTTOV_SAT_SERIES >=1 .AND. RTTOV_SAT_SERIES <= 29 &
 .AND. RTTOV_SAT_ID >= 1 .AND. RTTOV_SAT_ID <= 20) THEN
    WMO_SAT_NUM = SAT_NUM_RTTOV_TO_WMO( RTTOV_SAT_ID, RTTOV_SAT_SERIES )
  ELSE
    WMO_SAT_NUM = 999
  ENDIF

! Conversion of RTTOV satellite index into the satellite series in WMO classification
! accordingly to BUFR Code table 0 02 020, missing value is 511

  IF (RTTOV_SAT_SERIES >=1 .AND. RTTOV_SAT_SERIES <= 29 ) THEN
    WMO_SAT_SER = SAT_SER_RTTOV_TO_WMO( RTTOV_SAT_SERIES )
  ELSE
    WMO_SAT_SER = 511
  ENDIF

! Conversion of RTTOV instrument (sensor) index into satellite instrument index
! in WMO classification accordingly to BUFR Code table 0 02 019 - missing value is 2047

  IF (RTTOV_SAT_SENSOR >=0 .AND. RTTOV_SAT_SENSOR <= 53 ) THEN
    WMO_INST_ID = INST_ID_RTTOV_TO_WMO( RTTOV_SAT_SENSOR )
  ELSE
    WMO_INST_ID = 2047
  ENDIF

! Calculation of number of data fields and set of common parameters to all data description arrays

  NFIELD = NPARAM*NCHAN

! Cleaning of possible previous allocations

  DO IFIELD=1,NFIELD0
    IF (ALLOCATED(DATA(IFIELD) % FIELD) ) DEALLOCATE (DATA(IFIELD) % FIELD)
  ENDDO
  DO IFIELD=1,NFIELD
    DATA(IFIELD) % GRIB2_DESCRIPT(:) = IVALMISS
  ENDDO

! Description of data for grib2 coding (INTEGER, DIMENSION(50) :: GRIB2_DESCRIPT=IVALMISS),
! see module_write_grib2_data.F90

! Content of GRIB2_DESCRPT array see bottom of file

! Model index: 1 - Bolam, 2 - Moloch, 3 - Globo

  DATA(1:NFIELD) % GRIB2_DESCRIPT(1) = IMODEL

! Status and type of data

  DATA(1:NFIELD) % GRIB2_DESCRIPT(6) = 0  ! operational products
  DATA(1:NFIELD) % GRIB2_DESCRIPT(7) = 1  ! forecast
  DATA(1:NFIELD) % GRIB2_DESCRIPT(30) = 0 ! bit-map absent

! Grid parameters

  DATA(1:NFIELD) % GRIB2_DESCRIPT(2) = 1 ! grid template index - horizontal grid
  NX = NLON/NPOINT_JUMP
  NY = NLAT/NPOINT_JUMP

  DATA(1:NFIELD) % NX = NX
  DATA(1:NFIELD) % NY = NY
  DATA(1:NFIELD) % X0 = X0
  DATA(1:NFIELD) % Y0 = Y0
  DATA(1:NFIELD) % DX = DX*FLOAT(NPOINT_JUMP)
  DATA(1:NFIELD) % DY = DY*FLOAT(NPOINT_JUMP)
  DATA(1:NFIELD) % X00 = X1
  DATA(1:NFIELD) % Y00 = Y1

  DO IFIELD=1,NFIELD

! Initial date and time

    DATA(IFIELD) % IDATE0(1:5) = IDATE0(1:5)

! Date and time of current forecast term

    DATA(IFIELD) % IDATEC(1:5) = IDATEC(1:5)

  ENDDO

! Definition of time unit, forecast and period length in defined time unit

  DATA(1:NFIELD) % GRIB2_DESCRIPT(4) = 1 ! Time unit for Forecast time step and Statistical period : 0 minute, 1 hour, 2 day
  DATA(1:NFIELD) % IFORECAST = IPERIOD_INP(1)*24+IPERIOD_INP(2)   ! Forecast time step
  DATA(1:NFIELD) % IPERIOD = IPERIOD_ACCUM(1)*24+IPERIOD_ACCUM(2) ! Statistical period

! Product template index

  DATA(1:NFIELD) % GRIB2_DESCRIPT(3) = 32 ! forecast satellite

! Level/layer parameters (for forecast satellite products parameters of satellite platform and sensor)

  DO IFIELD=1,NFIELD
    DATA(IFIELD) % GRIB2_DESCRIPT(10:15) = IVALMISS
  ENDDO

! Definition of data fields and data parameters in grib2 terms

  DO IFIELD=1,NFIELD
    ALLOCATE(DATA(IFIELD) % FIELD(NX,NY))
  ENDDO

! Satellite and sensor parameters

  DATA(1:NFIELD) % GRIB2_DESCRIPT(10) = & ! level (layer) type (for forecast satellite products code of satellite platform and sensor)
 WMO_SAT_NUM*10000000 + WMO_SAT_SER*10000 + WMO_INST_ID

  DATA(1:NFIELD) % GRIB2_DESCRIPT(20) = 3 ! Discipline: Space products
  DATA(1:NFIELD) % GRIB2_DESCRIPT(21) = 1 ! Category: Quantitative products

  IFIELD = 0

  DO ICHAN=1,NCHAN

    DO IPARAM=1,NPARAM

      IFIELD = IFIELD+1

      DATA(IFIELD) % GRIB2_DESCRIPT(11) = NINT(SENSOR_CHAN_CW(ICHAN)*1.D3) ! Scaled value of channel central wave number
      DATA(IFIELD) % GRIB2_DESCRIPT(12) = 3                                ! Scaled factor of channel central wave number

      IF (IPARAM == 1) THEN ! Radiance
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 16 ! Parameter: Cloudy radiance (with respect to wave number) (W m^-1 sr^-1)
        IY2=0
        DO IY1=1,NLAT,NPOINT_JUMP
          IY2=IY2+1
          IX2=0
          DO IX1=1,NLON,NPOINT_JUMP
            IX2=IX2+1
            DATA(IFIELD) % FIELD(IX2,IY2)=RADIANCE(IX1,IY1,ICHAN)*1.E-5 ! (mW cm)/(ster m) -> (W)/(ster m)
          ENDDO
        ENDDO
      ENDIF

      IF (IPARAM == 2) THEN ! Brightness temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 14 ! Parameter: Cloudy brightness temperature (K)
        IY2=0
        DO IY1=1,NLAT,NPOINT_JUMP
          IY2=IY2+1
          IX2=0
          DO IX1=1,NLON,NPOINT_JUMP
            IX2=IX2+1
            DATA(IFIELD) % FIELD(IX2,IY2)=RADIANCE_BT(IX1,IY1,ICHAN)
          ENDDO
        ENDDO
      ENDIF

!      IF (IPARAM == 3) THEN ! Clear sky radiance
!        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 17 ! Parameter: Clear-sky radiance (with respect to wave number) (W m^-1 sr^-1)
!        IY2=0
!        DO IY1=1,NLAT,NPOINT_JUMP
!          IY2=IY2+1
!          IX2=0
!          DO IX1=1,NLON,NPOINT_JUMP
!            IX2=IX2+1
!            DATA(IFIELD) % FIELD(IX2,IY2)=RADIANCE_CLEAR(IX1,IY1,ICHAN)*1.E-5 ! (mW cm)/(ster m) -> (W)/(ster m)
!          ENDDO
!        ENDDO
!      ENDIF
!
!      IF (IPARAM == 4) THEN ! Clear sky brightness temperature
!        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 15 ! Parameter: Clear-sky brightness temperature (K)
!        IY2=0
!        DO IY1=1,NLAT,NPOINT_JUMP
!          IY2=IY2+1
!          IX2=0
!          DO IX1=1,NLON,NPOINT_JUMP
!            IX2=IX2+1
!            DATA(IFIELD) % FIELD(IX2,IY2)=RADIANCE_CLEAR_BT(IX1,IY1,ICHAN)
!          ENDDO
!        ENDDO
!      ENDIF
!
!      IF (IPARAM == 5) THEN ! Emissivity
!        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 201 ! Parameter: our local use: emissivity
!        IY2=0
!        DO IY1=1,NLAT,NPOINT_JUMP
!          IY2=IY2+1
!          IX2=0
!          DO IX1=1,NLON,NPOINT_JUMP
!            IX2=IX2+1
!            DATA(IFIELD) % FIELD(IX2,IY2)=EMIS_RTTOV(IX1,IY1,ICHAN)
!          ENDDO
!        ENDDO
!      ENDIF

    ENDDO

  ENDDO

! End of preparation of output data for coding in grib2 format

! Write output data in grib2 format

  CALL WRITE_GRIB2_DATA

RETURN
END

#endif

!-----------------------------------------------------------------------------

! Content of GRIB2_DESCRPT array:

!  1 - model index (1 Bolam, 2 Moloch, 3 Globo)
!  2 - grid template index (1 horizontal grid, 1000 vertical cross-section)
!  3 - product template index (0 instant, 8 statistical, 32 forecast satellite, 1 instant individual ensemble, 11 statistical individual ensemble)
!  4 - time unit index (0- minute, 1 hour, 2 day, 3 month, 4 year, 13 second)
!  5 - statistical elaboration type (for statistical products only) : 0 average, 1 accumulation, 2 maximum, 3 minimum
!  6 - production status of data (0 Operational products, 2 Research products, 3 Re-analysis products, 7 Sub-seasonal to seasonal prediction S2S)
!  7 - type of data (0 Analysis, 1 Forecast, 2 Analysis and forecast, 3 Control forecast, 4 Perturbed forecast)
!  8 - indicator of unit of time for the increment between successive fields used for statistical elaboration
! 10 - level (layer) type (for forecast satellite products code of satellite platform and sensor)
! 11 - first scaled value of level (layer)
! 12 - scale of first value of level (layer)
! 13 - second scaled value of level (layer)
! 14 - scale of second value of level (layer)
! 15 - level type of second level for a layer
! 20 - product discipline
! 21 - product category
! 22 - product parameter
! 30 - flag of bit-map presence: 0 not present, 1 present (use VALMISS)
! 31 - member number of Perturbed forecast
! 32 - member index of Perturbed forecast
