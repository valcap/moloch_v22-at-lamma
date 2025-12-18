! Last update 25/01/2019

subroutine read_grib2_data(ist,filename,nfile,ini_flag,flag_cut_paste,icentre_code,isubcentre_code,imodel_code, &
    iana,jana,nlev_atm_inp,nlev_atm_inp_max,nlev_soil_inp,nlev_soil_inp_max,                                    &
    x0a,y0a,x1a,y1a,dxa,dya,idate,iperiod,lev_list,lev_list_soil,level_type,                                    &
    npar3d,npar2d,npar3d_soil,field3d,field2d,field3d_soil,alev,blev,val_missing)

! Author: Oxana Drofa, ISAC-CNR, Bologna, Italy (o.drofa@isac.cnr.it)

! Procedure for reading and decoding of meteorol. parameter fields in GRIB2 format
! using the ECMWF GRIB_API library

! Version: for both IFS and GFS model data

! 18/08/2017: after 19 Jul. 2017 the new parameter "LANDN" (land-sea mask, nearest neighbour)
! of GFS data (index 0200218) is the only land-sea mask consistent with surface data - before,
! the "old" parameter "LAND" (index 0200000) was the consistent land-sea mask to be used.
! To maintain backward compatibility, the existence of both is checked here.
! Index 3 is attributed to LAND, index 46 to LANDN.

! 27/03/2017: automatic convers. from geop. height to geopotential introduced.
!------------------------------------------------------------------------------

!USE grib_api
USE eccodes

implicit none

integer :: ist, nfile, level_type, flag_cut_paste
character, dimension(nfile) :: filename*80
logical :: ini_flag
integer :: icentre_code, isubcentre_code, imodel_code, &
           iana, jana, nlev_atm_inp, nlev_atm_inp_max, nlev_soil_inp, nlev_soil_inp_max, &
           idate(5), iperiod(4), npar3d, npar3d_soil, npar2d
real :: x0a, y0a, x1a, y1a, dxa, dya, lev_list(nlev_atm_inp_max), lev_list_soil(nlev_soil_inp_max), &
        alev(nlev_atm_inp_max+1), blev(nlev_atm_inp_max+1), val_missing
real, dimension(iana,jana,nlev_atm_inp_max,npar3d) :: field3d
real, dimension(iana,jana,npar2d) :: field2d
real, dimension(iana,jana,nlev_soil_inp_max,npar3d_soil) :: field3d_soil

! For reading GRIB2 file

type grib_field
  integer :: index, lev_type, lev_top(2), lev_bot(2), nxgrib, nygrib, idate(5), iperiod(4), lev_tot
  real :: x00grib, y00grib, x0grib, y0grib, xngrib, yngrib, dxgrib, dygrib
  real, dimension(:,:), allocatable :: field_r
  real, dimension(:), allocatable :: alev, blev
end type grib_field

integer*4 :: ifile=11, igrib, iif, ierr
integer*4 :: jcount, ncount, nsize, param_discipl, param_categ, param_ind, &
             iscan, jscan, jscancons, iini, ifin, di, jini, jfin, dj, order(2), klev
real*4, dimension(:), allocatable :: value, vert_coord
real :: zzz

! For data processing

real, dimension(iana,jana,nlev_atm_inp_max) :: work_3d
real, dimension(iana,jana,nlev_soil_inp_max) :: work_3d_soil
integer, dimension(npar3d) :: nlev
integer, dimension(npar3d_soil) :: nlevg
integer, dimension(300,npar3d) ::  levels
real, dimension(30,npar3d_soil) ::  levels_soil
integer, dimension(nlev_atm_inp_max) :: lev, lev1
integer, dimension(nlev_soil_inp_max) :: lev_soil, lev1_soil
integer, dimension(nlev_atm_inp_max,npar3d) :: ifield3d
integer, dimension(npar2d) :: ifield2d
integer, dimension(nlev_soil_inp_max,npar3d_soil) :: ifield3d_soil
integer :: iflag, ihrfc, idayfc, ihrfc2, idepth, itop, lev_max
integer :: iii, i, j, ij, ir, jr, k, k1, nk, knum, first, ind1(1), ind, flag_lapse_rate

real*4 :: pi, zs1, zs2, zs3, zs4, zw2_q,  zw3_q, zw2_t, zw3_t

! For reading GRIB2 data

type (grib_field), dimension(1000) :: read_2d

!------------------------------------------------------------------------------
! 3D fields in input, GRIB2 codes:
!                                     discipline category number
! field3d(:,:,:,1) = Geopotential         00         03      005
! field3d(:,:,:,2) = Temperature          00         00      000
! field3d(:,:,:,3) = U-velocity           00         02      002
! field3d(:,:,:,4) = V-velocity           00         02      003
! field3d(:,:,:,5) = Specific Humidity    00         01      000
! field3d(:,:,:,6) = Relative Humidity    00         01      001
! field3d(:,:,:,7) = Cloud water m.r.     00         01      022
! field3d(:,:,:,8) = Cloud liquid water   00         01      083
! field3d(:,:,:,9) = Cloud ice water      00         01      084
! field3d(:,:,:,10) = Vertical velocity   00         02      009
! field3d(:,:,:,11) = Pressure (Pa)       00         03      000
! field3d(:,:,:,12) = Geometric height(m) 00         03      006
!------------------------------------------------------------------------------
! 3D fields at soil vertical levels in input, GRIB2 codes:
!                                         discipline category number
! field3d_soil(:,:,:,1) = Soil temperature         02         00      002
! field3d_soil(:,:,:,2) = Volumetric soil wetness  02         00      009 (0200192)
! field3d_soil(:,:,:,3) = Liquid moisture content  02         00      198
! field3d_soil(:,:,:,4) = Frozen moisture content  02         00      199
! field3d_soil(:,:,:,5) = Sea temperature          10         03      000
! field3d_soil(:,:,:,6) = Sea ice temperature      10         02      192
! field3d_soil(:,:,:,7) = Max soil water content  (Proportion = m**3/m**3)         02         03      009
! field3d_soil(:,:,:,8) = Min soil water content  (Proportion = m**3/m**3)         xx         xx      xxx
!------------------------------------------------------------------------------
! 2D fields in input, GRIB2 codes:
!                                                               discipline category par.number
! field2d(:,:, 1) = Geopot. height at the surface                  00         03      005
! or
! field2d(:,:, 1) = Geopotential at the surface (surf.)            00         03      004
! field2d(:,:, 2) = Geopot. height at the surf. (mod.lev. 1)       00         03      005
! or
! field2d(:,:, 2) = Geopotential at the surf. (mod.lev. 1 )        00         03      004
! field2d(:,:, 3) = Land/sea mask (GFS: valid before 19.07.2017)   02         00      000
! field2d(:,:, 4) = Surface pressure (Pa)                          00         03      000
! or
! field2d(:,:, 4) = Logarithm of surface pressure                  00         03      025
! field2d(:,:, 5) = Land surface skin temperature                  00         00      000 (00017)
! field2d(:,:, 6) = Total precipitation (kg/m/m)                   00         01      008
! field2d(:,:, 7) = Convective precipitation (kg/m/m)              00         01      010
! field2d(:,:, 8) = Total snow precipitation (m---> mm)            00         01      029
! field2d(:,:, 9) = Terrain height (m)                             02         00      007
! field2d(:,:,14) = Snow depth (equiv. water layer) (mm)           00         01      013 (01060)
! field2d(:,:,15) = Mean sea level pressure                        00         03      001 (03000)
! field2d(:,:,16) = Cloud cover (total)                            00         06      001
! field2d(:,:,17) = U wind component at 10 m                       00         02      002
! field2d(:,:,18) = V wind component at 10 m                       00         02      003
! field2d(:,:,19) = 2 metre temperature                            00         00      000
! field2d(:,:,20) = 2 metre dew point temperature                  00         00      006
! field2d(:,:,21) = Soil type                                      02         03      000
! field2d(:,:,22) = Fraction of sea ice                            10         02      000
! field2d(:,:,23) = Total column integrated water (kg/m/m)         00         01      051
! field2d(:,:,24) = Total column int. condensate (kg/m/m)          00         06      020
! field2d(:,:,25) = Land surface skin humidity (kg/kg)             00         01      000
! field2d(:,:,26) = Sea ice thickness (m)                          10         02      001
! field2d(:,:,27) = Lapse rate (K m-1)                             00         00      008
! field2d(:,:,28) = Lapse rate (K m-1) - 2                         00         00      008
! field2d(:,:,30) = Relative humidity at 2 m                       00         01      001
! field2d(:,:,31) = Low cloud cover (%)                            00         06      003
! field2d(:,:,32) = Medium cloud cover (%)                         00         06      004
! field2d(:,:,33) = High cloud cover (%)                           00         06      005
! field2d(:,:,34) = Average sensible heat flux (W/m/m)             00         00      011
! field2d(:,:,35) = Average latent heat flux (W/m/m)               00         00      010
! field2d(:,:,36) = Aver. downward Short-wave radiat. flux (W/m/m) 00         04      007
! field2d(:,:,37) = Average long-wave radiation flux (W/m/m)       00         05      000
! field2d(:,:,38) = Temperature at 80 m above ground (K)           00         00      000
! field2d(:,:,39) = Specific humidity at 80 m above ground (kg/kg) 00         01      000
! field2d(:,:,40) = Pressure at 80 m above ground (Pa)             00         03      000
! field2d(:,:,41) = U wind component at 80 m above ground (m/s)    00         02      002
! field2d(:,:,42) = V wind component at 80 m above ground (m/s)    00         02      003
! field2d(:,:,45) = Soil porosity (Proportion = m**3/m**3)         02         03      009
! field2d(:,:,46) = Land/sea mask (GFS, valid after 19.07.2017)    02         00      218
!------------------------------------------------------------------------------
! Level types, grib2 codes:

! Ground or water surface                  001
! Isobaric surface (Pa)                    100
! Mean sea level                           101
! Specified height level above m.sea.l.(m) 102
! Specified height level above ground (m)  103
! Sigma level (sigma value)                104
! Hybrid level                             105
! Depth below land surface (m)             106
! Depth below sea level (m)                160
!------------------------------------------------------------------------------

 pi = abs(acos(-1.))

! Reading of GRIB2 files

  iflag=0 ! for hybrid level type
  jcount = 0
  levels(:,:)=0
  levels_soil(:,:)=0

! Loop for GRIB2 files

 fileloop: do iif=1,nfile

! Open file

  call grib_open_file(ifile, trim(filename(iif)),"r")

! Turn on support for multi fields messages */

  call grib_multi_support_on()

! Turn off support for multi fields messages */
!  call grib_multi_support_off()

! Loop on all the messages in a opened file

50  continue

! Read next message into IGRIB

  call grib_new_from_file(ifile,igrib)
  if (igrib == -1 )  then
    goto 60
  endif

! GRIB field's parameters

  JCOUNT = JCOUNT + 1

  if (jcount ==1 ) then
    call grib_get(igrib,'originatingCentre', icentre_code)          ! Code of originating centre
    call grib_get(igrib,'subCentre', isubcentre_code)               ! Code of originating sub-centre
    call grib_get(igrib,'generatingProcessIdentifier', imodel_code) ! Code of generating process (model)
  endif

! Set value (-9999) for missing value in reading message

  call grib_set(igrib,'missingValue',val_missing)

! Parameter index:

  call grib_get(igrib,'discipline',param_discipl)
  call grib_get(igrib,'parameterCategory',param_categ)
  call grib_get(igrib,'parameterNumber',param_ind)
  read_2d(jcount) % index = param_discipl*100000 + param_categ*1000 + param_ind

! Level parameters:

  call grib_get(igrib,'typeOfFirstFixedSurface', read_2d(jcount) % lev_type)
  call grib_get(igrib,'scaledValueOfFirstFixedSurface', read_2d(jcount) % lev_top(1))
  call grib_get(igrib,'scaleFactorOfFirstFixedSurface', read_2d(jcount) % lev_top(2))
  call grib_get(igrib,'scaledValueOfSecondFixedSurface', read_2d(jcount) % lev_bot(1))
  call grib_get(igrib,'scaleFactorOfSecondFixedSurface', read_2d(jcount) % lev_bot(2))
  call grib_get(igrib,'numberOfCoordinatesValues', read_2d(jcount) % lev_tot)
  read_2d(jcount) % lev_tot = read_2d(jcount) % lev_tot /2+1

! Analysis date and time

  call grib_get(igrib,'year',   read_2d(jcount) % idate(1) )
  call grib_get(igrib,'month',  read_2d(jcount) % idate(2) )
  call grib_get(igrib,'day',    read_2d(jcount) % idate(3) )
  call grib_get(igrib,'hour',   read_2d(jcount) % idate(4) )
  call grib_get(igrib,'minute', read_2d(jcount) % idate(5) )

! Forecast range

  call grib_get(igrib,'indicatorOfUnitOfTimeRange', read_2d(jcount) % iperiod(1) )
  call grib_get(igrib,'forecastTime', read_2d(jcount) % iperiod(2) )
  call grib_get(igrib,'productDefinitionTemplateNumber', read_2d(jcount) % iperiod(4) )
  read_2d(jcount) % iperiod(3)=0
  if (read_2d(jcount) % iperiod(4)==8) then
    call grib_get(igrib,'lengthOfTimeRange', read_2d(jcount) % iperiod(3) )
    read_2d(jcount) % iperiod(2)=read_2d(jcount) % iperiod(2)+read_2d(jcount) % iperiod(3)
  endif

!print *,'idate0', read_2d(jcount) % idate(1:5)
!print *,'iperiod', read_2d(jcount) % iperiod(1:4)

! Grid parameters

  call grib_get(igrib,'numberOfPointsAlongAParallel',read_2d(jcount) % nxgrib)
  call grib_get(igrib,'numberOfPointsAlongAMeridian',read_2d(jcount) % nygrib)
  call grib_get(igrib,'longitudeOfSouthernPoleInDegrees',read_2d(jcount) % x00grib,ierr)
  if (ierr /= 0) read_2d(jcount) % x00grib=val_missing
  call grib_get(igrib,'latitudeOfSouthernPoleInDegrees',READ_2D(JCOUNT) % Y00GRIB,IERR)
  if (ierr /= 0) read_2d(jcount) % y00grib=val_missing
  call grib_get(igrib,'longitudeOfFirstGridPointInDegrees',read_2d(jcount) % x0grib)
  call grib_get(igrib,'latitudeOfFirstGridPointInDegrees',read_2d(jcount) % y0grib)
  call grib_get(igrib,'longitudeOfLastGridPointInDegrees',read_2d(jcount) % xngrib)
  call grib_get(igrib,'latitudeOfLastGridPointInDegrees',read_2d(jcount) % yngrib)
  call grib_get(igrib,'iDirectionIncrementInDegrees',read_2d(jcount) % dxgrib)
  call grib_get(igrib,'jDirectionIncrementInDegrees',read_2d(jcount) % dygrib)
  if (read_2d(jcount) % dxgrib < 0.) &
       read_2d(jcount) % dxgrib=abs(read_2d(jcount) % xngrib-read_2d(jcount) % x0grib)/float(read_2d(jcount) % nxgrib-1)
  if (read_2d(jcount) % dygrib < 0.) &
      read_2d(jcount) % dygrib=abs(read_2d(jcount) % yngrib-read_2d(jcount) % &
      y0grib)/float(read_2d(jcount) % nygrib-1)

  if (ini_flag) then
    iana= read_2d(jcount) % nxgrib
    jana= read_2d(jcount) % nygrib
    x1a= read_2d(jcount) % x0grib
    if(read_2d(jcount) % y0grib < read_2d(jcount) % yngrib) then
    y1a= read_2d(jcount) % y0grib
    else
    y1a= read_2d(jcount) % yngrib
    endif
    dxa= read_2d(jcount) % dxgrib
    dya= read_2d(jcount) % dygrib
    if (int(read_2d(jcount) % x00grib)==int(val_missing).and.int(read_2d(jcount) % &
        y00grib)==int(val_missing)) then  ! non rotated grid
      x0a=0.
      y0a=0.
!      if (x1a <= 180..and.x1a > 0.) x1a = x1a-360.
!      if (x1a > -180..and.x1a < 0.) x1a = x1a+360.
      if (x1a > 0..and.x1a+dxa*float(iana-1) >= 360.) x1a = x1a-360.
    else

! Rotated grid: conversion of South Pole coordinates (of rotation) into
! coordinates of intersection point of zero latitude and zero longitude

      x0a=read_2d(jcount) % x00grib
      y0a=acos(-sin((pi/180.)*read_2d(jcount) % y00grib))
      y0a=y0a*180./pi
      if (x1a > 180.) x1a=x1a-360.
      print *,'Input data grid is rotated;'
      print*, 'coordinates of rotation centre:', x0a, y0a
    endif
    goto 105
  endif

! Dinamic allocation array for reading array data

  call grib_get_size(igrib,'values',nsize)
  allocate(value(nsize))

! Get values array

  call grib_get(igrib,'values',value)

! Dinamic arrays allocation

  allocate ( read_2d(jcount)%field_r (read_2d(jcount)%nxgrib,read_2d(jcount)%nygrib) )

! For hybrid levels

  if (icentre_code == 98 ) then ! ECMWF
    if (read_2d(jcount) % lev_type==105.and.iflag==0) then

     iflag=1

     nsize=(read_2d(jcount)%lev_tot-1)*2
!!!     allocate ( read_2d(jcount)%alev(read_2d(jcount)%lev_tot+1) )
!!!     allocate ( read_2d(jcount)%blev(read_2d(jcount)%lev_tot+1) )

     allocate ( vert_coord(nsize) )

     call grib_get(igrib,'pv',vert_coord)

     do klev=1,nsize/2
!!!       read_2d(jcount)%alev(klev)=vert_coord(klev)
!!!       read_2d(jcount)%blev(klev)=vert_coord(klev+nsize/2)
      alev(klev)=vert_coord(klev)
      blev(klev)=vert_coord(klev+nsize/2)
     enddo

     deallocate ( vert_coord )

    endif
  endif

! Storing of read array

! Scanning mode flags

  call grib_get(igrib,'iScansNegatively',iscan)
  call grib_get(igrib,'jScansPositively',jscan)
  call grib_get(igrib,'jPointsAreConsecutive',jscancons)

  if (iscan == 0) then
    iini = 1
    ifin = read_2d(jcount)%nxgrib
    di = 1
  else
    iini = read_2d(jcount)%nxgrib
    ifin = 1
    di = -1
    zzz=read_2d(jcount)%x0grib
    read_2d(jcount)%x0grib=read_2d(jcount)%xngrib
    read_2d(jcount)%xngrib=zzz
  endif

  if (jscan == 0) then
    jini = read_2d(jcount)%nygrib
    jfin = 1
    dj = -1
    zzz=read_2d(jcount)%y0grib
    read_2d(jcount)%y0grib=read_2d(jcount)%yngrib
    read_2d(jcount)%yngrib=zzz
  else
    jini = 1
    jfin = read_2d(jcount)%nygrib
    dj = 1
  endif

  if ( jscancons == 0) then
    order = (/1,2/)
  else
    order = (/2,1/)
  endif

  read_2d(jcount)%field_r(iini:ifin:di,jini:jfin:dj) =  &
          reshape(value, (/read_2d(jcount)%nxgrib,read_2d(jcount)%nygrib/), order=order)

! For global data grid (when input longitude of grid points range from 0 to 360
! deg. but local grid longitude ranges from -180 to 180 deg.)

!  if (read_2d(jcount) % x0grib >= 0..and.x1a < 0..and. &
!  (read_2d(jcount) % xngrib+read_2d(jcount) % dxgrib) >= 359.9) then

  if (flag_cut_paste /= 0) then
    read_2d(jcount) % x0grib = read_2d(jcount) % x0grib -180.
    read_2d(jcount) % xngrib = read_2d(jcount) % xngrib -180.
    do j=1,read_2d(jcount) % nygrib
      value(1 : read_2d(jcount) % nxgrib)=read_2d(jcount) % field_r(1 : read_2d(jcount) % nxgrib, j)
      read_2d(jcount) % field_r(1 : int(read_2d(jcount) % nxgrib/2), j)= &
              value(int(read_2d(jcount) % nxgrib/2)+1 : read_2d(jcount) % nxgrib)
      read_2d(jcount) % field_r(int(read_2d(jcount) % nxgrib/2)+1 : read_2d(jcount) % nxgrib, j)= &
              value(1 : int(read_2d(jcount) % nxgrib/2))
    enddo
  endif

! End of grib message elaboration

 deallocate(value)
 call grib_release(igrib)

  goto 50

! End data extraction

60 call grib_close_file(ifile)

 enddo fileloop

! End of reading of input grib2 files

 ncount=jcount

! Extraction of meteorological parameter fields from 2d-fields

! Indexes of meteorological parameters in FIELD3D and FIELD2d array

  field3d(:,:,:,:)      = val_missing
  field2d(:,:,:)        = val_missing
  field3d_soil(:,:,:,:) = val_missing

  nlev(:)=0
  nlevg(:)=0
  ihrfc=0

  iperiod(3:4)=0
  first=0
  flag_lapse_rate=0

  do jcount=1,ncount

  if (jcount == 1) then

! Set date and time

   idate(:)=read_2d(jcount) % idate(:)
   iperiod(1:2)=read_2d(jcount) % iperiod(1:2)

! If time unit is hour:

   if (read_2d(jcount) % iperiod(1) == 1) then
      ihrfc=read_2d(jcount) % iperiod(2) ! (0 is indicator of analysis)
      idayfc=ihrfc/24
      ihrfc2=ihrfc-24*idayfc
   endif
   write(*,*) "Date (YYYY MM DD) and time (HH MM, analysis time) of the input fields:"
   print '(i11, i3.2, i3.2, i14.2, i3.2)', idate(:)
   print '(a44, i4.2)', "Validity time if forecast (00 if analysis):", ihrfc

  else ! jcount == 1

   if (read_2d(jcount) % index/=0200000.and.read_2d(jcount) % index/=0200007.and. &
read_2d(jcount) % index/=0203000.and. &
(read_2d(jcount) % lev_type==1.and.read_2d(jcount) % index/=0003005).and. &
(read_2d(jcount) % lev_type==105.and.read_2d(jcount) % index/=0003005)) then ! non constant field
   if (any(idate/=read_2d(jcount) % idate).or. &
       any(iperiod(1:2)/=read_2d(jcount) % iperiod(1:2))) then
     write (*,*) 'Caution: no coincidence in current date'
     write (*,*) idate(:)
     write (*,*) read_2d(jcount) % idate(:)
     write (*,*) iperiod(1:2)
     write (*,*) read_2d(jcount) % iperiod(1:2)
     write (*,'(a17,i5.5)') 'Parameter index',read_2d(jcount) % index,' message ',jcount
!     stop
   endif
   endif

  endif ! jcount == 1

  if (read_2d(jcount) % iperiod(3)/=0) iperiod(3)=read_2d(jcount) % iperiod(3)

! Parameters at atmospheric vertical levels (caution: level type is important!)

  if (read_2d(jcount) % lev_type == 100.or. read_2d(jcount) % lev_type == 105) then

    if (ist==1.and.first==0) then
      first=1
      if (read_2d(jcount) % lev_type == 100) then
        level_type=1 ! isobaric levels
      else
        level_type=2 ! hybrid or sigma levels
      endif
    endif

    iii=0

!   Geopotential

    if (read_2d(jcount) % index == 0003005.or.read_2d(jcount) % index == 0003004) iii=1

!   Temperature

    if (read_2d(jcount) % index == 0000000) iii=2

!   U wind component

    if (read_2d(jcount) % index == 0002002) iii=3

!   V wind component

    if (read_2d(jcount) % index == 0002003) iii=4

!   Specific humidity

    if (read_2d(jcount) % index == 0001000) iii=5

!   Relative humidity

    if (read_2d(jcount) % index == 0001001) iii=6

! Cloud total (water+ice) content

    if (read_2d(jcount) % index == 0001022) iii=7

! Cloud liquid water content

    if (read_2d(jcount) % index == 0001083) iii=8

! Cloud ice water content

    if (read_2d(jcount) % index == 0001084) iii=9

! Vertical velocity

    if (read_2d(jcount) % index == 0002009) iii=10

! Pressure (Pa)

    if (read_2d(jcount) % index == 0003000) iii=11

! Geometric height (m)

    if (read_2d(jcount) % index == 0003006) iii=12

    if (iii /= 0) then

      nlev(iii)=nlev(iii)+1
      levels(nlev(iii),iii) =  &
         nint( float(read_2d(jcount) % lev_top(1))*(0.1**read_2d(jcount) % lev_top(2)) )

      if (any( levels(nlev(iii),iii) == levels(1:nlev(iii)-1,iii) )) then

        ! This field has already been stored

        levels(nlev(iii),iii)=0
        nlev(iii)=nlev(iii)-1
      else

        ! New field

        klev=nlev(iii)
        if (read_2d(jcount) % lev_type == 105) then ! analysis model levels
          if (read_2d(jcount) % lev_tot > nlev_atm_inp_max) then
            write (*,*) 'Caution: the no. of vert. levels in input data exceeds maximum declared:', &
                 read_2d(jcount) % lev_tot,nlev_atm_inp_max,' for parameter with index',read_2d(jcount) % index
!            stop
          endif
        else
          if (klev > nlev_atm_inp_max) then
            write (*,*) 'Caution: the no. of atm. levels in input data exceeds maximum declared', &
                 nlev_atm_inp_max,' for parameter with index',read_2d(jcount) % index
!            stop
          endif
        endif
        field3d(:,:,klev,iii)=read_2d(jcount)%field_r(:,:)
        if (iii == 1.and.read_2d(jcount) % index == 0003005) then ! Conversion from geop. height (gpm) to geopotential (m**2 s**-2)
          field3d(:,:,klev,iii)=field3d(:,:,klev,iii)*9.8
        endif
      endif

!      if (read_2d(jcount) % lev_type == 100) iflag=iflag+1 ! pressure standard levels

    endif

   endif ! atmospheric vertical levels

!   Surface parameters at the 1-st hybrid model level

   if (read_2d(jcount) % lev_type == 105.and.read_2d(jcount) % lev_top(1) == 1) then

     iii=0

!     Geopotential at the surface (at level 1)

      if (read_2d(jcount) % index == 0003005.or.read_2d(jcount) % index == 0003004) iii=2

!     Natural log of surface pressure

      if (read_2d(jcount) % index == 0003025) iii=4

      if (iii/=0) then
        field2d(:,:,iii)=read_2d(jcount)%field_r(:,:)
        if (iii == 2.and.read_2d(jcount) % index == 0003005) then ! Conversion from geop. height (gpm) to geopotential (m**2 s**-2)
          field2d(:,:,iii)=field2d(:,:,iii)*9.8
        endif
      endif

  endif

! Parameters at soil levels or sea levels (level type is important!)

  if (read_2d(jcount) % lev_type == 106.or.read_2d(jcount) % lev_type == 160) then

    iii=0

    if (read_2d(jcount) % lev_type == 106) then ! depth below land surface (m)

!   Soil temperature

      if (read_2d(jcount) % index == 0000000.or.read_2d(jcount) % index == 0200002) iii=1

!   Soil volumetric water content

      if (read_2d(jcount) % index == 0200009.or.read_2d(jcount) % index == 0200192) iii=2

!   Liquid moisture content (kg m**-2)

      if (read_2d(jcount) % index == 0200198) iii=3

!   Frozen moisture content (kg m**-2)

      if (read_2d(jcount) % index == 0200199) iii=4

    endif

    if (read_2d(jcount) % lev_type == 160) then ! depth below sea level (m)

!   Sea water temperature

      if (read_2d(jcount) % index == 1003000) iii=5

!   Sea ice temperature

      if (read_2d(jcount) % index == 1002192) iii=6

    endif

    if (iii /= 0) then

      nlevg(iii)=nlevg(iii)+1

      if (read_2d(jcount) % lev_bot(1)<=0) then  ! level
        levels_soil(nlevg(iii),iii) =  &
                 float(read_2d(jcount) % lev_top(1))*(0.1**read_2d(jcount) % lev_top(2))
      else ! layer
        levels_soil(nlevg(iii),iii)= 0.5* &
                (float(read_2d(jcount) % lev_top(1))*(0.1**read_2d(jcount) % lev_top(2))+ &
                 float(read_2d(jcount) % lev_bot(1))*(0.1**read_2d(jcount) % lev_bot(2)))
      endif

      if (any( levels_soil(nlevg(iii),iii) == levels_soil(1:nlevg(iii)-1,iii) )) then

        ! This field has already been stored

        levels_soil(nlevg(iii),iii)=0
        nlevg(iii)=nlevg(iii)-1
      else

        ! New field

        klev=nlevg(iii)
        if (klev > nlev_soil_inp_max) then
          write (*,*) 'Caution: the no. of soil levels in input data exceeds maximum declared', &
                 nlev_soil_inp_max,' for parameter with index',read_2d(jcount) % index
!          stop
        endif
        field3d_soil(:,:,klev,iii)=read_2d(jcount)%field_r(:,:)
      endif

    endif

  endif ! soil or sea levels

! Parameters at ground and sea surface

  iii=0

  if (read_2d(jcount) % lev_type == 001) then

!   Geopotential at the surface

    if (read_2d(jcount) % index == 0003005.or.read_2d(jcount) % index == 0003004) iii=1

!   Land/sea mask (GFS: valid before 19.07.2017)

    if (read_2d(jcount) % index == 0200000) iii=3

!   Land/sea mask (GFS, valid after 19.07.2017)

    if (read_2d(jcount) % index == 0200218) iii=46

!   Surface pressure (Pa)

    if (read_2d(jcount) % index == 0003000) iii=4

!   Skin temperature

    if (read_2d(jcount) % index == 0000000.or.read_2d(jcount) % index == 0000017) iii=5

!   Total precipitation (kg/m/m)

    if (read_2d(jcount) % index == 0001008) iii=6

!   Convective precipitation (kg/m/m)

    if (read_2d(jcount) % index == 0001010) iii=7

!   Total snow precipitation (m) ---> Total snow precipitation (mm)

    if (read_2d(jcount) % index == 0001029) iii=8

!   Terrain height (m)

    if (read_2d(jcount) % index == 0200007) iii=9

!   Snow depth

    if (read_2d(jcount) % index == 0001013.or.read_2d(jcount) % index == 0001060) iii=14

!   Soil type

    if (read_2d(jcount) % index == 0203000) iii=21

!   Sea ice cover (fraction)

    if (read_2d(jcount) % index == 1002000) iii=22

!   Total column integrated water (kg/m/m)

    if (read_2d(jcount) % index == 0001051) iii=23

!   Total column-integrated condensate (kg/m/m)

    if (read_2d(jcount) % index == 0006020) iii=24

!   Skin humidity

    if (read_2d(jcount) % index == 0001000) iii=25

!   Sea ice thickness (m)

    if (read_2d(jcount) % index == 1002001) iii=26

!   Soil porosity (proportion) = maximum specific volumetric water soil content (m**3/m**3)

!    Lapse rate (K m-1)

     if (read_2d(jcount) % index == 0000008) then
       if (flag_lapse_rate == 0) then
         iii=27
         flag_lapse_rate=1
       else
         iii=28
       endif
     endif

!   Low cloud cover

    if (read_2d(jcount) % index == 0006003) iii=31

!   Medium cloud cover

    if (read_2d(jcount) % index == 0006004) iii=32

!   High cloud cover

    if (read_2d(jcount) % index == 0006005) iii=33

! Average sensible heat flux (W/m/m)

    if (read_2d(jcount) % index == 0000011) iii=34

! Average latent heat flux (W/m/m)

    if (read_2d(jcount) % index == 0000010) iii=35

! Average downward Short-wave radiation flux (W/m/m)

    if (read_2d(jcount) % index == 0004007) iii=36

! Average long-wave radiation flux (W/m/m)

    if (read_2d(jcount) % index == 0005000) iii=37

! Soil porosity (proportion = m**3/m**3)

    if (read_2d(jcount) % index == 0203009) iii=45

  endif ! level type is "land surface"

! Mean sea level pressure

  if (read_2d(jcount) % lev_type == 1.and.read_2d(jcount) % index == 0003001) iii=15
  if (read_2d(jcount) % lev_type == 101.and.(read_2d(jcount) % index == 0003000.or.read_2d(jcount) % index == 0003001)) iii=15

! Total cloud cover

  if ((read_2d(jcount) % lev_type == 1.or.read_2d(jcount) % lev_type == 244).and. &
          read_2d(jcount) % index == 0006001) iii=16

! Specified levels above ground

  if (read_2d(jcount) % lev_type == 103) then

!  10 m

   if (nint(float(read_2d(jcount) %lev_top(1))*(0.1**read_2d(jcount) %lev_top(2))) == 10) then

!    U wind component

     if (read_2d(jcount) % index == 0002002) iii=17

!    V wind component

     if (read_2d(jcount) % index == 0002003) iii=18

   endif

!  2 m

   if (nint(float(read_2d(jcount) %lev_top(1))*(0.1**read_2d(jcount) %lev_top(2))) == 2) then

!    Temperature

     if (read_2d(jcount) % index == 0000000) iii=19

!    Dew point temperature

     if (read_2d(jcount) % index == 0000006) iii=20

!    Relative humidity

     if (read_2d(jcount) % index == 0001001) iii=30

   endif

!  80 m (GFS)

   if (nint(float(read_2d(jcount) %lev_top(1))*(0.1**read_2d(jcount) %lev_top(2))) == 80) then

!    Temperature

     if (read_2d(jcount) % index == 0000000) iii=38

!    Specific humidity

     if (read_2d(jcount) % index == 0001000) iii=39

!    Pressure

     if (read_2d(jcount) % index == 0003000) iii=40

!    U wind component

     if (read_2d(jcount) % index == 0002002) iii=41

!    V wind component

     if (read_2d(jcount) % index == 0002003) iii=42

   endif

  endif

  if (iii/=0) field2d(:,:,iii)=read_2d(jcount)%field_r(:,:)
  if (iii == 1.and.read_2d(jcount) % index == 0003005) then ! Conversion from geop. height (gpm) to geopotential (m**2 s**-2)
    field2d(:,:,iii)=field2d(:,:,iii)*9.8
  endif
  if (iii == 8) field2d(:,:,iii)=field2d(:,:,iii)*1.e3 ! Total snow precipitation (m) ---> Total snow precipitation (mm)

  enddo  ! jcount=1,ncount

  do jcount=1,ncount
    deallocate(read_2d(jcount)%field_r)
!!!    deallocate(read_2d(jcount)%alev)
!!!    deallocate(read_2d(jcount)%blev)
  enddo

! Vertical level ordering for 3D field parameters

 if (nlev_atm_inp_max>0) then

! First check

  iii=0
  do klev=1,npar3d
    if (nlev(klev).ne.nlev(2).and.nlev(klev).ne.0) iii=iii+1
  enddo

!  print*, "npar3d =", npar3d

! Check (1st instant only) on the number of levels for each atmospheric variable

  if (iii > 0.and.ist <= 2) then ! Check only for 1st and 2nd instants (analysis an one fc)
    write(*,*) "Subroutine read_grib2_data:"
    write(*,*) "the no. of atmospheric levels is not the same for all variables;"
    write(*,*) "the no. of levels for each variable is printed below."
    write(*,'(a)') " Variable order:  GPH,    T,    U,    V,    Q,   RH, CTWC, CLWC, CIWC,    W,    P,    Z"
    write(*,'(a,20i6)') "               ", nlev(1:12)
  endif

! Second check

  if (nlev(2)>nlev_atm_inp_max) then
    write (*,*) 'Caution: the no. of vert. levels in input data exceeds maximum declared:', &
          nlev(2),nlev_atm_inp_max
!    stop
  endif

!  nlev_atm_inp=nlev(2)
  nlev_atm_inp=maxval(nlev(:))
  ind1=maxloc(nlev(:))
  ind=ind1(1)

 else

  nlev_atm_inp=0
  lev_list(:)=val_missing
  field3d(:,:,:,:)=val_missing

 endif

 if (nlev_atm_inp>0) then

! Level ordering

  do k=1,nlev_atm_inp
  lev1(k)=levels(k,ind)
  enddo

  nk=nlev_atm_inp
 100 continue
  lev_max=-100000
  do k=1,nk
  if (lev1(k) > lev_max) then
    lev_max=lev1(k)
    knum=k
  endif
  enddo
  lev(nk)=lev_max
  do k=knum,nk-1
  lev1(k)=lev1(k+1)
  enddo
  nk=nk-1
  if (nk >= 1) goto 100

! Values at vertical levels read from input file

  do k=1,nlev_atm_inp
   lev_list(k)=float(lev(k))
  enddo

! Levels ordering of all 3D fields

 do iii=1,npar3d

  work_3d(:,:,:)=val_missing
  do k=1,nlev_atm_inp
    do k1=1,nlev_atm_inp
      if (levels(k1,iii) == lev(k)) then
       work_3d(:,:,k)=field3d(:,:,k1,iii)
      endif
    enddo
  enddo
  field3d(:,:,:,iii)=work_3d(:,:,:)

  enddo

 endif ! nlev_atm_inp > 0

! Vertical level ordering for 3D soil field parameters

 if (nlev_soil_inp_max>0) then

! First check

  iii=0
  do klev=1,npar3d_soil
    if (nlevg(klev).ne.nlevg(1).and.nlevg(klev).ne.0) iii=iii+1
  enddo

! Check (1st instant only) on the number of levels for each atmospheric variable

  if (iii > 0.and.ist == 1) then  ! Check only for 1st instant
    write(*,*) "Subroutine read_grib2_data:"
    write(*,*) "The no. of soil levels is not the same for all variables;"
    write(*,*) "the no. of soil levels for each variable is printed below."
    write(*,*) "Variable order: TSoil  QSoil QLiqSoil QIceSoil SeaWatT SeaIceT"
    write(*,'(a,20i8)') "            ",nlevg(1:6)
!    stop
  endif

! Second check

  if (nlevg(1)>nlev_soil_inp_max) then
    write (*,*) 'Caution: the no. of soil levels in input data exceeds maximum declared:', &
           nlevg(1),nlev_soil_inp_max
!    stop
  endif

  nlev_soil_inp=maxval(nlevg(:))
  ind1=maxloc(nlevg(:))
  ind=ind1(1)

 else

  nlev_soil_inp=0
  lev_list_soil(:)=val_missing
  field3d_soil(:,:,:,:)=val_missing

 endif

 if (nlev_soil_inp>0) then

! Level ordering

  do k=1,nlev_soil_inp
  lev1_soil(k)=nint(levels_soil(k,ind)*1.e3)
  enddo

  nk=nlev_soil_inp
 101 continue
  lev_max=-100000
  do k=1,nk
  if (lev1_soil(k) > lev_max) then
    lev_max=lev1_soil(k)
    knum=k
  endif
  enddo
  lev_soil(nk)=lev_max
  do k=knum,nk-1
  lev1_soil(k)=lev1_soil(k+1)
  enddo
  nk=nk-1
  if (nk >= 1) goto 101

! Values at vertical levels read from input file

  do k=1,nlev_soil_inp
   lev_list_soil(k)=float(lev_soil(k))*1.e-3
  enddo

! Level ordering of all 3D fields

 do iii=1,npar3d_soil

  work_3d_soil(:,:,:)=val_missing
  do k=1,nlev_soil_inp
    do k1=1,nlev_soil_inp
      if (nint(levels_soil(k1,iii)*1.e3) == lev_soil(k)) work_3d_soil(:,:,k)=field3d_soil(:,:,k1,iii)
    enddo
  enddo
  field3d_soil(:,:,:,iii)=work_3d_soil(:,:,:)

  enddo

 endif  ! nlev_soil_inp > 0

! Check for the presence of data

  ifield3d(:,:)=0
  ifield2d(:)=0
  ifield3d_soil(:,:)=0

  do iii=1,npar3d
  do k=1,nlev_atm_inp
  do j=1,jana,jana
  do i=1,iana,iana
    if (int(field3d(i,j,k,iii)) == int(val_missing)) ifield3d(k,iii)=ifield3d(k,iii)+1
  enddo
  enddo
  enddo
  enddo

  do iii=1,npar2d
  do j=1,jana,jana
  do i=1,iana,iana
    if (int(field2d(i,j,iii)) == int(val_missing)) ifield2d(iii)=ifield2d(iii)+1
  enddo
  enddo
  enddo

  do iii=1,npar3d_soil
  do k=1,nlev_soil_inp
  do j=1,jana,jana
  do i=1,iana,iana
    if (int(field3d_soil(i,j,k,iii)) == int(val_missing)) ifield3d_soil(k,iii)=ifield3d_soil(k,iii)+1
  enddo
  enddo
  enddo
  enddo

! For optional check

!  do iii=1,npar3d
!  do k=1,nlev_atm_inp
!    if (ifield3d(k,iii) > 0) write (*,'(a27,i3,a50,i3,i8)')  &
!       "Caution: 3D Field (index=",III,") not available in the input GRIB file at level",K,LEV(K)
!  enddo
!  enddo
!  do iii=1,npar2d
!    if (ifield2d(iii) > 0) write (*,'(a27,i3,a40)')  &
!        "Caution: 2D Field (index=",III,") not available in the input GRIB file"
!  enddo
!  do iii=1,npar3d_soil
!  do k=1,nlev_soil_inp
!    if (ifield3d_soil(k,iii) > 0) write (*,'(a27,i3,a50,i3,i8)')  &
!        "Caution: 3D Field Soil (index=",III,") not available in the input GRIB file at level", &
!         k,lev_list_soil(k)
!  enddo
!  enddo

! end of degrib part

105 return
end
!=======================================================================
