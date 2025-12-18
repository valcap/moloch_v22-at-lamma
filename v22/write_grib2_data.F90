!------------------------------------------------------------------------------
! Last update 17/02/2016
!------------------------------------------------------------------------------

MODULE GRIB2_CODING_DATA

IMPLICIT NONE

! For writing GRIB file

  CHARACTER (LEN=200) :: GRIBSAMPLENAME
  NAMELIST /GRIB_SAMPLE/ GRIBSAMPLENAME
  CHARACTER (LEN=50) :: OUTPUT_FILE_NAME

  INTEGER :: IGRIB_EDITION=2, ICENTRE_CODE=80, ICENTRE_SUBCODE=102, IBITNUM=16 ! CNR-ISAC codes

! Definition of data type to pass data to the procedure of grib2 format coding
! (write_grib2_data.F90)

  INTEGER, PARAMETER :: NFIELD0=500
  INTEGER :: NFIELD
  INTEGER, PARAMETER :: IVALMISS=-9999, IVALMISS1=255, IVALMISS2=0
  REAL*4, PARAMETER ::  VALMISS=-9999.

  TYPE DATA_FOR_CODING

    INTEGER :: NX=IVALMISS, NY=IVALMISS, IFORECAST=IVALMISS, IPERIOD=IVALMISS, &
 N_VERT_COORD_PAR=IVALMISS, IND_VERT_COORD=IVALMISS
    REAL*4 :: X0=VALMISS, Y0=VALMISS, X00=VALMISS, Y00=VALMISS, DX=VALMISS, DY=VALMISS, &
 LON_INI, LON_FIN, LAT_INI, LAT_FIN
    REAL*4, DIMENSION(:), ALLOCATABLE :: VERT_COORD_PAR
    INTEGER, DIMENSION(5) :: IDATE0=IVALMISS, IDATEC=IVALMISS

    INTEGER, DIMENSION(50) :: GRIB2_DESCRIPT=IVALMISS1
    REAL*4, DIMENSION(:,:), ALLOCATABLE :: FIELD

! Content of GRIB2_DESCRIPT array:

!  1  - model index (1 Bolam, 2 Moloch, 3 Globo)
!  2  - grid template index (1 horizontal grid, 1000 vertical cross-section)
!  3  - product template index (0 instant, 8 statistical, 32 forecast satellite, 1 instant individual ensemble, 11 statistical individual ensemble)
!  4  - time unit index (0- minute, 1 - hour, 2 - day, 3 - month, 4 - year, 13 - second)
!  5  - statistical elaboration type (for statistical prodicts only) : 0 average, 1 accumulation, 2 maximum, 3 minimum
!  6  - production status of data (0 Operational products, 2 Research products, 3 Re-analysis products, 7 Sub-seasonal to seasonal prediction project (S2S) test)
!  7  - type of data (0 Analysis, 1 Forecast, 2 Analysis and forecast, 3 Control forecast, 4 Perturbed forecast)
!  8  - indicator of unit of time for the increment between the successive fields used for statistical elaboration
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

  END TYPE DATA_FOR_CODING

  TYPE (DATA_FOR_CODING), DIMENSION(NFIELD0) :: DATA

  CONTAINS

  SUBROUTINE READ_GRIB_SAMPLE_NAME

! Complete pass and name of grib sample file

    OPEN (11, FILE="grib_sample.inp", STATUS="OLD")
    READ (11, GRIB_SAMPLE)
    CLOSE (11)

  END SUBROUTINE READ_GRIB_SAMPLE_NAME

END MODULE GRIB2_CODING_DATA

!------------------------------------------------------------------------------

SUBROUTINE WRITE_GRIB2_DATA

! Procedure of coding and writing of meteorological parameter fields in GRIB2 format
! using the ECMWF ECCODES library

! AUTHOR: Oxana Drofa, ISAC-CNR, ITALY (o.drofa@isac.cnr.it)

!------------------------------------------------------------------------------

USE eccodes
USE GRIB2_CODING_DATA

IMPLICIT NONE

INTEGER :: IND_LINE_TYPE=0 ! For space vertical cross-section (Grib template 1000)
                           ! Horizontal line type = 0 (Rhumb) or 1 (Great circle)

INTEGER :: IGRIBSAMPLE1, IGRIBSAMPLE, IGRIBOUT, IGRIBCLONE, IERR, I, III, J, K, IFIELD, NSIZE, &
 NX, NY, IND_TEMPLATE_GRID, IND_TEMPLATE_PROD, &
 IDATEINI, ITIMEINI, LEV_TYPE, WMO_SAT_NUM, WMO_SAT_SER, WMO_INST_ID
REAL :: PI, Y0, DX, DY, XFIN, YFIN

REAL*4, DIMENSION(:), ALLOCATABLE :: VALUE

 PI = ABS(ACOS(-1.))

! Complete pass and name of grib sample file

 CALL READ_GRIB_SAMPLE_NAME

 IF (NFIELD==0) THEN
  PRINT *,'There are not defined files for coding of GRIB format data in file ',TRIM(OUTPUT_FILE_NAME)
  GOTO 105
 ENDIF

!------------------------------------------------------------------------------

! A new grib message is loaded from an existing sample

 III=0
 DO I=200,1,-1
   IF (GRIBSAMPLENAME(I:I).EQ.' ') III=I
 ENDDO

 CALL GRIB_OPEN_FILE(IGRIBSAMPLE1, GRIBSAMPLENAME(1:III-1), 'r')

 CALL GRIB_NEW_FROM_FILE (IGRIBSAMPLE1, IGRIBSAMPLE, IERR)

 IF (IERR/=0) THEN
   PRINT *,'Not found the sample of GRIB file: ',GRIBSAMPLENAME(1:III-1)
   PRINT *,'Procedure writing output file in GRIB format interrupts'
   GOTO 105
 ENDIF

! Open output grib file

  III=0
  DO I=50,1,-1
   IF (OUTPUT_FILE_NAME(I:I).EQ.' ') III=I
  ENDDO

 CALL GRIB_OPEN_FILE (IGRIBOUT, OUTPUT_FILE_NAME(1:III-1), 'w')

! Set of header parameters in Sample Grib Message

! General parameters

 CALL GRIB_SET (IGRIBSAMPLE, 'editionNumber', IGRIB_EDITION)                             ! Format Grib Number
 CALL GRIB_SET (IGRIBSAMPLE, 'originatingCentre', ICENTRE_CODE)                          ! Code of Originating Centre
 CALL GRIB_SET (IGRIBSAMPLE, 'subCentre', ICENTRE_SUBCODE)                               ! Code of Originating SubCentre
 CALL GRIB_SET (IGRIBSAMPLE, 'generatingProcessIdentifier', DATA(1) % GRIB2_DESCRIPT(1)) ! Code of generating Process (Model):
         ! 1 - Bolam, 2 - Moloch, 3 -Globo
!! CALL GRIB_SET (IGRIBSAMPLE, 'dataRepresentationTemplateNumber', 50002) ! Compression mode
 CALL GRIB_SET (IGRIBSAMPLE, 'bitsPerValue', IBITNUM) ! Precision
!! CALL GRIB_SET (IGRIBSAMPLE, 'secondOrderFlags', 128)
!! CALL GRIB_SET (IGRIBSAMPLE, 'packingType', 'grid_second_order') ! Compression mode

! Production status of data

! 0 Operational products, 2 Research products, 3 Re-analysis products, 7 Sub-seasonal to seasonal prediction project (S2S) test
 CALL GRIB_SET (IGRIBSAMPLE, 'productionStatusOfProcessedData', DATA(1) % GRIB2_DESCRIPT(6))

! Type of data
! 0 Analysis, 1 Forecast, 2 Analysis and forecast, 3 Control forecast, 4 Perturbed forecast
 CALL GRIB_SET (IGRIBSAMPLE, 'typeOfProcessedData', DATA(1) % GRIB2_DESCRIPT(7))

 LOOP_FIELDS: DO IFIELD=1,NFIELD

! Create Clone grib message form Sample grib message
! Clone grib message will be used to create output grib message

   IGRIBCLONE=0
   CALL GRIB_CLONE (IGRIBSAMPLE, IGRIBCLONE)

! Set TemplateNumber for grid type and grid parameters

   IND_TEMPLATE_GRID = DATA(IFIELD) %GRIB2_DESCRIPT(2) ! grid template index (0 - horizontal grid, 1000 - vertical cross-section)

   CALL GRIB_SET (IGRIBCLONE, 'gridDefinitionTemplateNumber', IND_TEMPLATE_GRID)

   NX = DATA(IFIELD) % NX
   NY = DATA(IFIELD) % NY

   DX = DATA(IFIELD) % DX
   DY = DATA(IFIELD) % DY

   IF (IND_TEMPLATE_GRID == 1) THEN ! horizontal grid

     IF (ABS(DATA(IFIELD) % X0)>1.E-3.OR.ABS(DATA(IFIELD) % Y0)>1.E-3) THEN
       CALL GRIB_SET (IGRIBCLONE, 'typeOfGrid','rotated_ll') ! Grid type is "rotated lat-lon"
       CALL GRIB_SET (IGRIBCLONE, 'longitudeOfSouthernPoleInDegrees', DATA(IFIELD) % X0)
! South Pole coordinates (of rotation)

       Y0 = DATA(IFIELD) % Y0
       Y0 = ASIN(-COS(Y0*PI/180.))*180./PI
       CALL GRIB_SET (IGRIBCLONE, 'latitudeOfSouthernPoleInDegrees', Y0)

     ELSE

       CALL GRIB_SET (IGRIBCLONE, 'typeOfGrid','regular_ll') ! Grid type is "regular lat-lon (not rotated)"

     ENDIF

     CALL GRIB_SET (IGRIBCLONE, 'numberOfPointsAlongAParallel', NX)
     CALL GRIB_SET (IGRIBCLONE, 'numberOfPointsAlongAMeridian', NY)
     CALL GRIB_SET (IGRIBCLONE, 'longitudeOfFirstGridPointInDegrees', DATA(IFIELD) % X00)
     CALL GRIB_SET (IGRIBCLONE, 'latitudeOfFirstGridPointInDegrees', DATA(IFIELD) % Y00)
     XFIN = DATA(IFIELD) % X00 + DX * FLOAT(NX-1)
     YFIN = DATA(IFIELD) % Y00 + DY * FLOAT(NY-1)
     CALL GRIB_SET (IGRIBCLONE, 'longitudeOfLastGridPointInDegrees', XFIN)
     CALL GRIB_SET (IGRIBCLONE, 'latitudeOfLastGridPointInDegrees', YFIN)
     CALL GRIB_SET (IGRIBCLONE, 'iDirectionIncrementInDegrees', DX)
     IF (DY >= 0.) THEN
       CALL GRIB_SET (IGRIBCLONE, 'jDirectionIncrementInDegrees', DY)
     ELSE
       CALL GRIB_SET (IGRIBCLONE, 'jDirectionIncrementInDegrees', -DY)
     ENDIF
   ENDIF

   IF (IND_TEMPLATE_GRID == 1000) THEN ! space vertical cross-section
     CALL GRIB_SET (IGRIBCLONE, 'numberOfHorizontalPoints', NX)
     CALL GRIB_SET (IGRIBCLONE, 'typeOfHorizontalLine', IND_LINE_TYPE)
     CALL GRIB_SET (IGRIBCLONE, 'latitudeOfFirstGridPoint', INT(DATA(IFIELD) % LAT_INI*1.E6))
     CALL GRIB_SET (IGRIBCLONE, 'latitudeOfLastGridPoint',  INT(DATA(IFIELD) % LAT_FIN*1.E6))
     IF (DATA(IFIELD) % LON_INI >= 0.) THEN
       CALL GRIB_SET (IGRIBCLONE, 'longitudeOfFirstGridPoint', INT(DATA(IFIELD) % LON_INI*1.E6))
     ELSE
       CALL GRIB_SET (IGRIBCLONE, 'longitudeOfFirstGridPoint', INT((DATA(IFIELD) % LON_INI+360.)*1.E6))
     ENDIF
     IF (DATA(IFIELD) % LON_FIN >= 0.) THEN
       CALL GRIB_SET (IGRIBCLONE, 'longitudeOfLastGridPoint', INT(DATA(IFIELD) % LON_FIN*1.E6))
     ELSE
       CALL GRIB_SET (IGRIBCLONE, 'longitudeOfLastGridPoint', INT((DATA(IFIELD) % LON_FIN+360.)*1.E6))
     ENDIF
     CALL GRIB_SET (IGRIBCLONE, 'numberOfVerticalPoints', NY)
     CALL GRIB_SET (IGRIBCLONE, 'NC', NY)
     CALL GRIB_SET (IGRIBCLONE, 'meaningOfVerticalCoordinate',  DATA(IFIELD) % IND_VERT_COORD) ! Type for vertical coordinate (3.15.table):
!                100 - Pressure (Pa)
!                102 - Altitude above mean sea level (m)
     CALL GRIB_SET (IGRIBCLONE, 'verticalCoordinate', 0) ! Explicit coordinate values set
!!     CALL GRIB_SET (IGRIBCLONE, '????', DATA(IFIELD) % VERT_COORD_PAR(:)) ! In eccodes not defined (definitions/grib2/template.3.1000.def)
!!   This problem is avoided by using pv in section 4 (template 0)
   ENDIF ! Grib template type

   CALL GRIB_SET (IGRIBCLONE, 'iScansNegatively', 0) ! Scanning for X XINI < XFIN
   CALL GRIB_SET (IGRIBCLONE, 'jScansPositively', 1) ! Scanning for Y YINI < YFIN
   IF (IND_TEMPLATE_GRID == 1.AND.DY < 0.) CALL GRIB_SET (IGRIBCLONE, 'jScansPositively', 0) ! Scanning for Y YINI > YFIN
   CALL GRIB_SET (IGRIBCLONE, 'jPointsAreConsecutive', 0) ! Scanning order: first - X, second - Y

! Analysis date and time

   IDATEINI=DATA(IFIELD) % IDATE0(1)*10000 + DATA(IFIELD) % IDATE0(2)*100 + DATA(IFIELD) % IDATE0(3)
   ITIMEINI=DATA(IFIELD) % IDATE0(4)*100 + DATA(IFIELD) % IDATE0(5)
   CALL GRIB_SET (IGRIBCLONE, 'dataDate', IDATEINI)
   CALL GRIB_SET (IGRIBCLONE, 'dataTime', ITIMEINI)

! Set TemplateNumber for product type and period of forecast definition

   IND_TEMPLATE_PROD = DATA(IFIELD) % GRIB2_DESCRIPT(3) ! product template index :
!                         0 - instant, 8 - statistical, 32 - forecast satellite

   CALL GRIB_SET (IGRIBCLONE, 'productDefinitionTemplateNumber', IND_TEMPLATE_PROD)
   CALL GRIB_SET (IGRIBCLONE, 'typeOfProcessedData', DATA(IFIELD) % GRIB2_DESCRIPT(7))
   CALL GRIB_SET (IGRIBCLONE, 'indicatorOfUnitOfTimeRange', DATA(IFIELD) % GRIB2_DESCRIPT(4)) ! time unit index (0- minute, 1 - hour, 2 - day, 3 - month, 4 - year, 13 - second)
   CALL GRIB_SET (IGRIBCLONE, 'forecastTime', DATA(IFIELD) % IFORECAST ) ! Forecast time step in defined time unit

! For ensemble forecasts

   IF (IND_TEMPLATE_PROD == 1.OR.IND_TEMPLATE_PROD == 11) THEN
     CALL GRIB_SET (IGRIBCLONE, 'perturbationNumber', DATA(IFIELD) % GRIB2_DESCRIPT(31))
     CALL GRIB_SET (IGRIBCLONE, 'numberOfForecastsInEnsemble', DATA(IFIELD) % GRIB2_DESCRIPT(32))
   ENDIF

! Statistical forecast variable

   IF (IND_TEMPLATE_PROD == 8.OR.IND_TEMPLATE_PROD == 11) THEN
     CALL GRIB_SET (IGRIBCLONE, 'forecastTime', MAX( DATA(IFIELD) % IFORECAST - DATA(IFIELD) % IPERIOD, 0 )) ! DATA(IFIELD) % IPERIOD - Statistical period in defined time unit
     CALL GRIB_SET (IGRIBCLONE, 'yearOfEndOfOverallTimeInterval', DATA(IFIELD) % IDATEC(1))
     CALL GRIB_SET (IGRIBCLONE, 'monthOfEndOfOverallTimeInterval', DATA(IFIELD) % IDATEC(2))
     CALL GRIB_SET (IGRIBCLONE, 'dayOfEndOfOverallTimeInterval', DATA(IFIELD) % IDATEC(3))
     CALL GRIB_SET (IGRIBCLONE, 'hourOfEndOfOverallTimeInterval', DATA(IFIELD) % IDATEC(4))
     CALL GRIB_SET (IGRIBCLONE, 'minuteOfEndOfOverallTimeInterval', DATA(IFIELD) % IDATEC(5))
     CALL GRIB_SET (IGRIBCLONE, 'secondOfEndOfOverallTimeInterval', 0)
     CALL GRIB_SET (IGRIBCLONE, 'typeOfStatisticalProcessing', DATA(IFIELD) % GRIB2_DESCRIPT(5)) ! Statistical elaboration type : 0 - average, 1 - accumulation, 2 - maximum, 3 - minimum
     CALL GRIB_SET (IGRIBCLONE, 'typeOfTimeIncrement', 2)
     CALL GRIB_SET (IGRIBCLONE, 'indicatorOfUnitForTimeRange', DATA(IFIELD) % GRIB2_DESCRIPT(4)) ! time unit index (0- minute, 1 - hour, 2 - day, 3 - month, 4 - year, 13 - second)
     IF (DATA(IFIELD) % IFORECAST -  DATA(IFIELD) % IPERIOD >= 0) THEN
       CALL GRIB_SET (IGRIBCLONE, 'lengthOfTimeRange',  DATA(IFIELD) % IPERIOD)
     ELSE
       CALL GRIB_SET (IGRIBCLONE, 'lengthOfTimeRange', 0)
     ENDIF
   ENDIF

   IF (IND_TEMPLATE_GRID /= 1000) THEN ! horizontal grid

! Set Level parameters or Satellite instrument parameters

     LEV_TYPE = DATA(IFIELD) % GRIB2_DESCRIPT(10) ! level (layer) type (for forecast satellite products code of satellite platform and sensor)

     IF (IND_TEMPLATE_PROD /= 32 ) THEN ! No satellite products

       IF (LEV_TYPE /= IVALMISS) THEN
         CALL GRIB_SET (IGRIBCLONE, 'typeOfFirstFixedSurface', LEV_TYPE)
       ELSE
         CALL GRIB_SET_MISSING (IGRIBCLONE, 'typeOfFirstFixedSurface')
       ENDIF
       IF (DATA(IFIELD) % GRIB2_DESCRIPT(11) /= IVALMISS) THEN
         CALL GRIB_SET (IGRIBCLONE, 'scaledValueOfFirstFixedSurface', DATA(IFIELD) % GRIB2_DESCRIPT(11))
       ELSE
         CALL GRIB_SET_MISSING (IGRIBCLONE, 'scaledValueOfFirstFixedSurface')
       ENDIF
       IF (DATA(IFIELD) % GRIB2_DESCRIPT(12) /= IVALMISS) THEN
         CALL GRIB_SET (IGRIBCLONE, 'scaleFactorOfFirstFixedSurface', DATA(IFIELD) % GRIB2_DESCRIPT(12))
       ELSE
         CALL GRIB_SET_MISSING (IGRIBCLONE, 'scaleFactorOfFirstFixedSurface')
       ENDIF
       IF (DATA(IFIELD) % GRIB2_DESCRIPT(13) /= IVALMISS) THEN
         CALL GRIB_SET (IGRIBCLONE, 'scaledValueOfSecondFixedSurface', DATA(IFIELD) % GRIB2_DESCRIPT(13))
       ELSE
         CALL GRIB_SET_MISSING (IGRIBCLONE, 'scaledValueOfSecondFixedSurface')
       ENDIF
       IF (DATA(IFIELD) % GRIB2_DESCRIPT(14) /= IVALMISS) THEN
         CALL GRIB_SET (IGRIBCLONE, 'scaleFactorOfSecondFixedSurface', DATA(IFIELD) % GRIB2_DESCRIPT(14))
       ELSE
         CALL GRIB_SET_MISSING (IGRIBCLONE, 'scaleFactorOfSecondFixedSurface')
       ENDIF
       IF (DATA(IFIELD) % GRIB2_DESCRIPT(15) /= IVALMISS) THEN
         CALL GRIB_SET (IGRIBCLONE, 'typeOfSecondFixedSurface', DATA(IFIELD) % GRIB2_DESCRIPT(15))
       ELSE
         CALL GRIB_SET (IGRIBCLONE, 'typeOfSecondFixedSurface', IVALMISS1)
       ENDIF

       IF (LEV_TYPE == 104.OR.LEV_TYPE == 105) THEN ! Sigma level or Hybrid level
         NSIZE = DATA(IFIELD) % N_VERT_COORD_PAR
         CALL GRIB_SET (IGRIBCLONE, 'PVPresent', 1)
         CALL GRIB_SET (IGRIBCLONE, 'numberOfVerticalCoordinateValues', NSIZE)
         ALLOCATE(VALUE(NSIZE))
         VALUE(:) = DATA(IFIELD) % VERT_COORD_PAR(:)
         CALL GRIB_SET (IGRIBCLONE, 'pv', VALUE)
         DEALLOCATE(VALUE)
       ELSE
         CALL GRIB_SET (IGRIBCLONE, 'PVPresent', 0)
         CALL GRIB_SET (IGRIBCLONE, 'numberOfVerticalCoordinateValues', 0)
       ENDIF

     ELSE ! Satellite simulated products

       WMO_SAT_NUM = LEV_TYPE / 10000000
       WMO_SAT_SER=(LEV_TYPE - WMO_SAT_NUM*10000000) / 10000
       WMO_INST_ID=LEV_TYPE - WMO_SAT_NUM*10000000 - WMO_SAT_SER*10000
       CALL GRIB_SET (IGRIBCLONE, 'NB', 1)
       CALL GRIB_SET (IGRIBCLONE, 'satelliteNumber', WMO_SAT_NUM)
       CALL GRIB_SET (IGRIBCLONE, 'satelliteSeries', WMO_SAT_SER)
       CALL GRIB_SET (IGRIBCLONE, 'instrumentType', WMO_INST_ID)
       CALL GRIB_SET (IGRIBCLONE, 'scaledValueOfCentralWaveNumber', DATA(IFIELD) % GRIB2_DESCRIPT(11))
       CALL GRIB_SET (IGRIBCLONE, 'scaleFactorOfCentralWaveNumber', DATA(IFIELD) % GRIB2_DESCRIPT(12))

     ENDIF

   ELSE ! Vertical cross-section

! Vertical coordinate definition (to avoid the lack in section 3 template 1000)

     IF (NY > 1) THEN
       CALL GRIB_SET (IGRIBCLONE, 'PVPresent', 1)
       CALL GRIB_SET (IGRIBCLONE, 'numberOfVerticalCoordinateValues', NY)
       CALL GRIB_SET (IGRIBCLONE, 'pv', DATA(IFIELD) % VERT_COORD_PAR(1:NY))
     ENDIF

   ENDIF ! Grid template type

! Set Parameter Index (Code Table 4.2)

   CALL GRIB_SET (IGRIBCLONE, 'discipline', DATA(IFIELD) % GRIB2_DESCRIPT(20))        ! product discipline
   CALL GRIB_SET (IGRIBCLONE, 'parameterCategory', DATA(IFIELD) % GRIB2_DESCRIPT(21)) ! product category
   CALL GRIB_SET (IGRIBCLONE, 'parameterNumber', DATA(IFIELD) % GRIB2_DESCRIPT(22))   ! product parameter
!print *,'discipline', DATA(IFIELD) % GRIB2_DESCRIPT(20)
!print *,'parameterCategory', DATA(IFIELD) % GRIB2_DESCRIPT(21)
!print *,'parameterNumber', DATA(IFIELD) % GRIB2_DESCRIPT(22)

! Set 2D Data Field

   IF (DATA(IFIELD) % GRIB2_DESCRIPT(30) == 1) THEN ! bit-map present
     CALL GRIB_SET (IGRIBCLONE, 'missingValue', VALMISS)
     CALL GRIB_SET (IGRIBCLONE,  'bitmapPresent', 1)
   ELSE
     CALL GRIB_SET (IGRIBCLONE,  'bitmapPresent', 0) ! bit-map absent
   ENDIF

   NSIZE= NX * NY
   ALLOCATE(VALUE(NSIZE))
   DO J=1,NY
   DO I=1,NX
     III=(J-1)*NX+I
     VALUE(III) = DATA(IFIELD) % FIELD(I,J)
   ENDDO
   ENDDO
   CALL GRIB_SET (IGRIBCLONE, 'values', VALUE)
   DEALLOCATE(VALUE)

! Write grib message

   CALL GRIB_WRITE (IGRIBCLONE, IGRIBOUT)
   CALL GRIB_RELEASE (IGRIBCLONE)

 ENDDO LOOP_FIELDS

! Close output grib file

 CALL GRIB_CLOSE_FILE (IGRIBOUT)

 PRINT *,NFIELD,' 2D data fields are wtitten in GRIB2 format file ',OUTPUT_FILE_NAME

! Deallocate data fields

 DO IFIELD=1,NFIELD
   IF (ALLOCATED(DATA(IFIELD) % VERT_COORD_PAR)) DEALLOCATE(DATA(IFIELD) % VERT_COORD_PAR)
   IF (ALLOCATED(DATA(IFIELD) % FIELD)) DEALLOCATE(DATA(IFIELD) % FIELD)
 ENDDO
 NFIELD=0

105 RETURN
END SUBROUTINE WRITE_GRIB2_DATA

!------------------------------------------------------------------------------

! Discipline, category and parameter codes - see in ../eccodes-[version_number]/definitions/grib2/tables/2

! Level type codes (../eccodes-[version_number]/definitions/grib2/tables/0/4.5.table):

! Ground or water surface                        001
! Isobaric surface (Pa)                          100
! Mean sea level                                 101
! Specified height level above mean sea lev. (m) 102
! Specified height level above ground (m)        103
! Sigma level (sigma value)                      104
! Hybrid level                                   105
! Depth below land surface (m)                   106
! Depth below sea level (m)                      160
! pt Isentropic (theta) level Potential Temp.(K) 107
! Tropopause                                     007
! Maximum wind level                             006
! For RTTOV simulated data = WMO_SAT_NUM*10000000 + WMO_SAT_SER*10000 + WMO_INST_ID (BUFR Code tables: 0 01 007, 0 02 020, 0 02 019)

!------------------------------------------------------------------------------