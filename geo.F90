    program geo

! Last update 06/02/2019  V19

! Defines all the geographical fields (orography, sea, lakes), landuse, soil
! types and soil physical properties in the domains of Globo, Bolam or Moloch.
! Grid parameters in file geo.inc (except for Globo)
! Output is written to geo.bin (to be read by preprocessings and models).

! Compilation:
!   ncargf77 -D ncarg (-D globo) geo.F90 : to check geography
!   ifort geo.F90                        : to create file geo.bin for Bolam/Moloch
!   ifort -D globo geo.F90               : to create file geo.bin for Globo

#ifdef globo
!    parameter (nlon=1538, nlat=1058, npolcap=12)   ! 19 km globo
    parameter (nlon=2050, nlat=1442, npolcap=16)   ! 14 km globo
    parameter (dlon=360./float(nlon-2), dlat=180./float(nlat-1), x0d=0., y0d=0., alon0=0., alat0=-90.)
#else
    include 'geo.inc'  ! contains parameters: dlon, dlat, x0d, y0d, alon0, alat0
#endif

    real, dimension(nlon,nlat) :: orog, orogstd, fmask, alont, alatt, flake
    real suolo(nlon,nlat,14), vegeta(nlon,nlat,13), soilvegpar(nlon,nlat,21), vect(14)
    real(kind=8) :: zfac, zx0, zy0, zlatt, zlatv, zlont, zzlatt, zargt, zaarg
    integer mty(nlon,nlat)

#ifdef globo
    print*, "GLOBO model"
    do jlat = 1, nlat
    do jlon = 1, nlon
    alont(jlon,jlat) = alon0 + (jlon-1)*dlon
    alatt(jlon,jlat) = alat0 + (jlat-1)*dlat
    enddo
    enddo
#else
    zfac = dabs(dacos(-1.d0))/180.d0
    zx0  = dble(x0d)*zfac
    zy0  = dble(y0d)*zfac
    do jlat = 1, nlat
    zlatv = (dble(alat0)+                (jlat-1)*dble(dlat))*zfac
    zlatt = (dble(alat0)+dble(dlat)/2.d0+(jlat-1)*dble(dlat))*zfac
    do jlon = 1, nlon
    zlont  = (dble(alon0)+(jlon-1)*dble(dlon))*zfac
    zzlatt = 1.d0/zfac*dasin( dcos(zy0)*dsin(zlatt) + dsin(zy0)*dcos(zlatt)*dcos(zlont))
    zargt  = -dsin(zlatt)*dsin(zy0)+dcos(zy0)*dcos(zlatt)*dcos(zlont)
    zaarg  = zargt/dcos(zfac*zzlatt)
    alatt(jlon,jlat) = zzlatt
    if (zaarg.lt.-1..and.zaarg.gt.-1.00001) zaarg = -1.d0
    if (zaarg.gt. 1..and.zaarg.lt. 1.00001) zaarg =  1.d0
      if (zlont.lt.0.d0) then
      alont(jlon,jlat) = 1.d0/zfac*(zx0-dacos(zaarg))
      else
      alont(jlon,jlat) = 1.d0/zfac*(zx0+dacos(zaarg))
      endif
    enddo
    enddo
#endif

! Grid coordinates

    print*, "Grid extremes (deg., T-points):"
    print*, ' min. alont,    max. alont'
    print*, minval(alont), maxval(alont)
    print*, ' min. alatt,    max. alatt'
    print*, minval(alatt), maxval(alatt)

! Orogvegdef defines orography, orog. std., fmask, flake and soil-vegetation types
! (in vegeta)

    call orogvegdef (orog, orogstd, fmask, flake, vegeta,        &
                  alont, alatt, dlon, dlat, x0d, y0d, nlon, nlat)
    print*, 'Orog, orogstd, fmask, flake, vegeta defined in orogvegdef.'

#ifdef ncarg
! Routine for plot (ncar-graphics)

!    dis = 250. 
!    ris = 1755.
    dis = 250. 
    ris = 1500.4
    call bolcolor (orog,fmask,alont,alatt,nlon,nlat,dlat,ris,dis)
#else
#ifdef globo
    call polavert (orog   , nlon, nlat, npolcap)
    call polavert (orogstd, nlon, nlat, npolcap)
    orogstd = max (orogstd, 0.)
#endif

    call defsoil (suolo, soilvegpar, vegeta, orogstd, fmask, flake, alont, alatt, nlon, nlat)
    print*, 'Suolo defined in defsoil.'

#ifdef globo
    call polavert (soilvegpar(:,:,15), nlon, nlat, npolcap)
    call polavert (soilvegpar(:,:,16), nlon, nlat, npolcap)
    soilvegpar(:,:,15) = max (soilvegpar(:,:,15), 1.e-8)   ! rgm
    soilvegpar(:,:,16) = max (soilvegpar(:,:,16), 1.e-8)   ! rgq
#endif

!    mty = maxloc(suolo,3)
!    write(99) float(mty(:,:))
!    write(99) fmask

    open (10, file='geo.bin', form='unformatted', status='unknown')

!  alat0 = latitude of first T point!

    write(10) nlon, nlat, x0d, y0d, dlon, dlat, alon0, alat0
    do jlat = 1, nlat
    write(10) (orog(jlon,jlat), jlon=1,nlon)
    enddo
    do jlat = 1, nlat
    write(10) (fmask(jlon,jlat), jlon=1,nlon)
    enddo
    do k = 1, 14
    do jlat = 1, nlat
    write (10) (suolo(jlon,jlat,k), jlon=1,nlon)
    enddo
    enddo
    do k = 1, 13
    do jlat = 1, nlat
    write (10) (vegeta(jlon,jlat,k), jlon=1,nlon)
    enddo
    enddo
    do k = 1, 21
    do jlat = 1, nlat
    write (10) (soilvegpar(jlon,jlat,k), jlon=1,nlon)
    enddo
    enddo
    close (10)
    print*, 'File geo.bin written.'
#endif

    stop
    end
!##################################################################################################
    subroutine orogvegdef (orog, orogstd, fmask, flake, vegeta,        &
                        alont, alatt, dlon, dlat, x0d, y0d, nlon, nlat)

! Defines orography, orog. std., fmask, flake and soil-vegetation types (in vegeta)

! Landuse (soil types, vegetation):

!-1  Permanent ice shelf (redefined as bare ground)
! 0  Water (sea, lakes/rivers - put in fsea, flake)
! 1  Evergreen needleleaf forest
! 2  Evergreen broadleaf forest (does not appear in this dataset)
! 3  Deciduous needleleaf forest (does not appear in this dataset)
! 4  Deciduous broadleaf forest
! 5  Mixed cover
! 6  Woodland
! 7  Wooded grassland
! 8  Closed shrubland
! 9  Open shrubland
! 10 Grassland
! 11 Cropland
! 12 Bare ground
! 13 Urban and built-up

    parameter (imoro=360*120, jmoro=180*120)
    integer*1, dimension(:,:), allocatable :: ivegbase
    integer*2, dimension(:,:), allocatable :: iorogbase
    real, dimension(nlon,nlat) :: orog, fmask, orogstd, grad, alont, alatt, wrk, fsea, flake
    real vegeta(nlon,nlat,13)

    i1 = nint(.5*dlon*120./max(cos(3.1415927/180.*y0d),.17))
    j1 = nint(.5*dlat*120.)
    i1 = max(i1,1)
    j1 = max(j1,1)
    print*, 'Original pixels (lon, lat) contained in one gridbox:'
    print*, i1, j1

    dist02 = 8.e-6 + (dlon**2 + dlat**2)/5.5   ! squared Gaussian half-width
    ntot = (2*j1+1)*(2*i1+1)
    write(*,'(a, f9.5)') ' Gaussian half-width (deg) =', sqrt(dist02)

!********************************************************************************
! Orography and fmask
!********************************************************************************

    orog    = 0.
    orogstd = 0.
    flake   = 0.
    fsea    = 0.
    fmask   = 0.
    vegeta  = 0.

    allocate (iorogbase(imoro,jmoro))
    open (17, file='orography120.bin', status='unknown', form='unformatted')
    do j = 1, jmoro
    read (17) iorogbase(:,j)
    enddo
    close (17)
    allocate (ivegbase(imoro,jmoro))
    open (17, file='landuse120.bin', status='unknown', form='unformatted')
    do j = 1, jmoro
    read (17) ivegbase(:,j)
    enddo
    close (17)

    do jlat = 1, nlat
    do jlon = 1, nlon
    ztot = 0.
    zalont = alont(jlon,jlat)
    if (zalont.ge.180.) zalont = zalont-360.
    i0 = int((zalont+180.)*120.)           +1   ! pixel containing the longitude alont
    j0 = int((90.-alatt(jlon,jlat))*120.)  +1   ! pixel containing the latitude alatt
!    zalat0 = 90.-float(j0-1)/120. - 1./240.    ! latitude of j-th pixel
!    zalon0 = float(i0-1)/120. -180. + 1./240.  ! longitude of i-th pixel

    do j = j0-2*j1, j0+2*j1
    jj = j
    if (jj.gt.jmoro) jj = jmoro
    if (jj.le.1    ) jj = 1
    do i = i0-2*i1, i0+2*i1
    ii = i
    if (ii.gt.imoro) ii = ii - imoro
    if (ii.lt.1    ) ii = ii + imoro
    dalati = 90.-float(jj-1)/120. - 1./240.   - alatt(jlon,jlat)
    daloni = float(ii-1)/120. -180. + 1./240. - zalont
    if (daloni.gt. 180.) daloni = daloni-360.
    if (daloni.lt.-180.) daloni = daloni+360.
    dist2 = daloni**2 + dalati**2
    factij = exp(-dist2/dist02) ! gaussian weight
    orog(jlon,jlat) = orog(jlon,jlat) + float(iorogbase(ii,jj))*factij
    ztot = ztot + factij
    enddo
    enddo

    orog(jlon,jlat) = orog(jlon,jlat)/ztot
    enddo
    enddo

! Filtering of orography (before computing the orography standard deviation)

    wei = .55/(1. + 7.*.5*(dlon+dlat))
!    print*, "2nd order filter weight =", wei
    do nit = 1,4
    call filt2 (orog, nlon, nlat, wei)
    enddo

    do jlat = 1, nlat
    do jlon = 1, nlon
    ztot = 0.
    zalont = alont(jlon,jlat)
    if (zalont.ge.180.) zalont = zalont-360.
    i0 = int((zalont+180.)*120.)           +1   ! pixel containing the longitude alont
    j0 = int((90.-alatt(jlon,jlat))*120.)  +1   ! pixel containing the latitude alatt

    do j = j0-2*j1, j0+2*j1
    jj = j
    if (jj.gt.jmoro) jj = jmoro
    if (jj.le.1    ) jj = 1
    do i = i0-2*i1, i0+2*i1
    ii = i
    if (ii.gt.imoro) ii = ii - imoro
    if (ii.lt.1    ) ii = ii + imoro
    dalati = 90.-float(jj-1)/120. - 1./240.   - alatt(jlon,jlat)
    daloni = float(ii-1)/120. -180. + 1./240. - zalont
    if (daloni.gt. 180.) daloni = daloni-360.
    if (daloni.lt.-180.) daloni = daloni+360.
    dist2 = daloni**2 + dalati**2
    factij = exp(-dist2/dist02) ! gaussian weight
    orogstd(jlon,jlat) = orogstd(jlon,jlat) + (orog(jlon,jlat)-float(iorogbase(ii,jj)))**2*factij
    ztot = ztot + factij
    enddo
    enddo
    orogstd(jlon,jlat) = sqrt(orogstd(jlon,jlat)/ztot) ! orography standard deviation

!********************************************************************************
! Landuse (soil type and vegetation)
!********************************************************************************

    do j = j0-j1, j0+j1
    jj = j
    if (jj.gt.jmoro) jj = jmoro
    if (jj.le.1    ) jj = 1
    do i = i0-i1, i0+i1
    ii = i
    if (ii.gt.imoro) ii = ii - imoro
    if (ii.lt.1    ) ii = ii + imoro

    if (ivegbase(ii,jj).eq.0) then
    if (iorogbase(ii,jj).eq.0.or.iorogbase(ii,jj).eq.1) then 
    fsea(jlon,jlat) = fsea(jlon,jlat) + 1.
    else
    flake(jlon,jlat) = flake(jlon,jlat) + 1.
    endif
    endif

      if (ivegbase(ii,jj).eq.1 ) vegeta(jlon,jlat,1 ) = vegeta(jlon,jlat,1 ) + 1.
      if (ivegbase(ii,jj).eq.2 ) vegeta(jlon,jlat,2 ) = vegeta(jlon,jlat,2 ) + 1.
      if (ivegbase(ii,jj).eq.3 ) vegeta(jlon,jlat,3 ) = vegeta(jlon,jlat,3 ) + 1.
      if (ivegbase(ii,jj).eq.4 ) vegeta(jlon,jlat,4 ) = vegeta(jlon,jlat,4 ) + 1.
      if (ivegbase(ii,jj).eq.5 ) vegeta(jlon,jlat,5 ) = vegeta(jlon,jlat,5 ) + 1.
      if (ivegbase(ii,jj).eq.6 ) vegeta(jlon,jlat,6 ) = vegeta(jlon,jlat,6 ) + 1.
      if (ivegbase(ii,jj).eq.7 ) vegeta(jlon,jlat,7 ) = vegeta(jlon,jlat,7 ) + 1.
      if (ivegbase(ii,jj).eq.8 ) vegeta(jlon,jlat,8 ) = vegeta(jlon,jlat,8 ) + 1.
      if (ivegbase(ii,jj).eq.9 ) vegeta(jlon,jlat,9 ) = vegeta(jlon,jlat,9 ) + 1.
      if (ivegbase(ii,jj).eq.10) vegeta(jlon,jlat,10) = vegeta(jlon,jlat,10) + 1.
      if (ivegbase(ii,jj).eq.11) vegeta(jlon,jlat,11) = vegeta(jlon,jlat,11) + 1.
      if (ivegbase(ii,jj).eq.12) vegeta(jlon,jlat,12) = vegeta(jlon,jlat,12) + 1.
      if (ivegbase(ii,jj).eq.-1) vegeta(jlon,jlat,12) = vegeta(jlon,jlat,12) + 1.
      if (ivegbase(ii,jj).eq.13) vegeta(jlon,jlat,13) = vegeta(jlon,jlat,13) + 1.
    enddo
    enddo
    flake  (jlon,jlat)   = flake(jlon,jlat)/float(ntot)
    fsea   (jlon,jlat)   = fsea (jlon,jlat)/float(ntot)
    vegeta (jlon,jlat,:) = vegeta(jlon,jlat,:)/float(ntot)

    enddo
    enddo

!********************************************************************************
    deallocate (iorogbase)
    deallocate (ivegbase)
!********************************************************************************

! Filtering of orography standard deviation (same as for orography)

    do nit = 1,4
    call filt2 (orogstd, nlon, nlat, wei)
    enddo

    grad = 0.
    do jlat = 2, nlat-1
    do jlon = 2, nlon-1
    grad(jlon,jlat) = sqrt((orog(jlon+1,jlat)-orog(jlon-1,jlat))**2+(orog(jlon,jlat+1)-orog(jlon,jlat-1))**2)
    orogstd(jlon,jlat) = orogstd(jlon,jlat)*(1.-flake(jlon,jlat)-fsea(jlon,jlat))
    enddo
    enddo

#ifdef globo
    grad(:,1) = grad(:,2)
    grad(:,nlat) = grad(:,nlat-1)
    grad(1 ,:) = grad(nlon-1,:)
    grad(nlon,:) = grad(2,:)
    orogstd(:,1)    = orogstd(:,2)
    orogstd(:,nlat) = orogstd(:,nlat-1)
    orogstd(1 ,:)   = orogstd(nlon-1,:)
    orogstd(nlon,:) = orogstd(2,:)
#endif

! "Envelope orography" using modified std. deviation

    zcovar = .45
    do jlat = 1, nlat
    do jlon = 1, nlon
    std2 = orogstd(jlon,jlat)/(0.0008*grad(jlon,jlat)+1.)
    orog(jlon,jlat) = orog(jlon,jlat) + zcovar*std2
    if(orog(jlon,jlat).gt.1.)                                                       &
    orog(jlon,jlat) = orog(jlon,jlat)*(1. + 0.25/(700. + orog(jlon,jlat))*std2)
    enddo
    enddo

! Reset orography over seas and lakes

    do jlat = 1, nlat
    do jlon = 1, nlon
    fsea(jlon,jlat) = fsea(jlon,jlat)**1.7        ! ad hoc correction to enlarge islands
    orog(jlon,jlat) = orog(jlon,jlat)*(1.-fsea(jlon,jlat))
    if (fsea(jlon,jlat).gt.0.5) orog(jlon,jlat) = 0.
    enddo
    enddo

! The following defines flatter lakes (at high res. only)

    if (dlat.lt.0.04) then
    call seatemp (orog, flake, wrk, nlon, nlat, 4, 1, .9)
    do jlat = 1,nlat
    do jlon = 1,nlon
    if (flake(jlon,jlat).gt.0.5) orog(jlon,jlat) = wrk(jlon,jlat)
    enddo
    enddo
    endif

    fmask = fsea + flake
    fmask = min (fmask, 1.)
    fmask = max (fmask, 0.)

!  Consistency of vegeta and fmask:
!  sum(vegeta(1:13))+fmask = 1. is imposed

    do jlat = 1, nlat
    do jlon = 1, nlon
    zs = sum(vegeta(jlon,jlat,1:13))
    if (zs.eq.0.) then
    vegeta(jlon,jlat,12) = 1.-fmask(jlon,jlat) ! put into bare soil
    else
    zc = (1.-fmask(jlon,jlat))/zs
    vegeta(jlon,jlat,:) = vegeta(jlon,jlat,:)*zc
    endif
    enddo
    enddo

#ifdef globo
    orog (nlon,:) = orog (2,     :)
    orog (1,   :) = orog (nlon-1,:)
    orog (1:nlon,nlat) = 0.
    fmask(nlon,:) = fmask(2,     :)
    fmask(1,   :) = fmask(nlon-1,:)
    fmask(1:nlon,1) = 0.
    fmask(1:nlon,nlat) = 1.
    vegeta(nlon,:,: ) = vegeta(2,     :,:)
    vegeta(1,   :,: ) = vegeta(nlon-1,:,:)
    vegeta(:,   1,1:13) = 0.
    vegeta(:,   1,12  ) = 1.
    vegeta(:,nlat,1:13) = 0.
#endif

    return
    end subroutine orogvegdef
!##################################################################################################
  subroutine seatemp (tana,zlsm,sst,iana,jana,niter,nfilt,w)

! Defines SST extending open sea values towards and beyond coasts

! 1) Definition of a logical vector initially TRUE only at "open sea" points
!    (defined as sea points not adjacent to any land point)
! 2) Expansion is repeated NITER times:
!    for any point with logical flag FALSE (land or coastal)
!    temperature is re-defined as average with only adjacent points having flag TRUE,
!    if they exist, and in this case the flag is set to TRUE.
!    In this way the temp. of "open sea" is gradually extended to coastal sea points
!    and also to land sea points near the coast (depending on the no. of iterations NITER).
! 3) A 4-point filter is optionally applied NFILT times with weight W
!   (W: 0.5-1; if W=1, no filtering)

  real  tana(iana,jana),zlsm(iana,jana)
  real  sst(iana,jana)
  real  twork(iana,jana)
  logical ilwork(iana,jana),ilsst(iana,jana)

! Identification of open sea points

  do j=1,jana
  do i=1,iana
  twork(i,j) = tana(i,j)

    if (zlsm(i,j).lt.0.5) then        ! land
    ilwork(i,j) = .false.
    else                              ! sea
    im1 = max(i-1,1)
    ip1 = min(i+1,iana)
    jm1 = max(j-1,1)
    jp1 = min(j+1,jana)

      if(zlsm(i,jp1).lt.0.5.or.zlsm(i,jm1).lt.0.5.or.zlsm(ip1,j).lt.0.5.or.zlsm(im1,j).lt.0.5) then
      ilwork(i,j) = .false.
      else
      ilwork(i,j) = .true.
      endif

    endif
  enddo
  enddo

! Correction of the flag on lateral boundaries (it is assumed that a sea point is a coastal point)

 do i=1,iana
   ilwork(i,1) = .false.
   ilwork(i,jana) = .false.
 enddo
 do j=1,jana
   ilwork(1,j) = .false.
   ilwork(iana,j) = .false.
 enddo

! Expansion of the sea temperature field towards the land

do 300 n=1,niter
  do j=1,jana
    do i=1,iana
      if (ilwork(i,j)) then
        sst(i,j) = twork(i,j)
        ilsst(i,j) = .true.
      else
        im1 = max(i-1,1)
        ip1 = min(i+1,iana)
        jm1 = max(j-1,1)
        jp1 = min(j+1,jana)

        nsi = 0
        tsi = 0.

        do je = jm1,jp1
        do ie = im1,ip1
          if (ilwork(ie,je)) then
          nsi = nsi+1
          tsi = tsi+twork(ie,je)
          endif
        enddo
        enddo

        if (nsi.gt.0) then
        sst(i,j) = tsi/nsi
        ilsst(i,j) = .true.
        else
        sst(i,j) = twork(i,j)
        ilsst(i,j) = .false.
        endif
      endif
   enddo
 enddo

 do j=1,jana
 do i=1,iana
 twork(i,j) = sst(i,j)
 ilwork(i,j) = ilsst(i,j)
 enddo
 enddo

300 continue

! Filtering

  do n=1,nfilt
   do j=2,jana-1
   do i=2,iana-1
   t4med=(twork(i-1,j)+twork(i+1,j)+twork(i,j-1)+twork(i,j+1))/4.
   sst(i,j) = w*twork(i,j) + (1.-w)*t4med
   enddo
   enddo

   do j=2,jana-1
   do i=2,iana-1
   twork(i,j) = sst(i,j)
   enddo
   enddo
  enddo

  return
  end subroutine seatemp
!##################################################################################################
   subroutine filt2 (p, im, jm, anu)

   real p(im,jm), p2(im,0:jm+1)

!  2-grid interval filter over the whole domain

   do jlat = 1, jm
   do jlon = 2, im-1
   p2(jlon,jlat) = .25*(p(jlon-1,jlat)+p(jlon+1,jlat))+.5*p(jlon,jlat)
   enddo
#ifdef globo
   p2(1 ,jlat) = p2(im-1,jlat)
   p2(im,jlat) = p2(2,jlat)
#else
   p2(1 ,jlat) = p2(2,jlat)
   p2(im,jlat) = p2(im-1,jlat)
#endif
   enddo

#ifdef globo
   do jlon = 1, im
   jlon1 = jlon + im/2
   if (jlon1 > im) jlon1 = jlon1-im
   p2(jlon,jm+1) = p2(jlon1,jm-1)
   p2(jlon,0   ) = p2(jlon1,2   )
   enddo
#else
   do jlon = 1, im
   p2(jlon,jm+1) = p2(jlon,jm)
   p2(jlon,0   ) = p2(jlon,1 )
   enddo
#endif

   do jlat = 1, jm
   p(1:im,jlat) = (1.-anu)*p(1:im,jlat) + anu*    &
                  (.25*p2(1:im,jlat+1)+.25*p2(1:im,jlat-1)+.5*p2(1:im,jlat))
   enddo

   return
   end subroutine filt2
!##################################################################################################
   subroutine defsoil (suolo, soilvegpar, vegeta, orogstd, fmask, flake, alont, alatt, nlon, nlat)

   parameter (imsol=360*12, jmsol=180*12)
   integer*1 isoilbase(imsol,jmsol,14)
   real, dimension(nlon,nlat) :: orogstd, fmask, flake, alont, alatt
   real suolo(nlon,nlat,14), vegeta(nlon,nlat,13), soilvegpar(nlon,nlat,21)
   real, dimension(14) :: psig0, rogdr, cgdry, ykghy, qgmin, qgmax, qgwilt, &
                          dryalb, wetalb, dryemis1, dryemis2, wetemis1, wetemis2
   real, dimension(13) :: vegalb, vegemis1, vegemis2, root, vegroughm, vegroughq

! Soil properties:
! 1  Sand
! 2  Loamy sand
! 3  Sandy loam
! 4  Silty loam
! 5  Loam
! 6  Sandy clay loam
! 7  Silty clay loam
! 8  Clay loam
! 9  Sandy clay
! 10 Silty clay
! 11 Clay
! 12 Peat
! 13 Stone
! 14 Glaciers

! Saturated moisture potential of soil (at qg=qgmax) (m)

      data psig0/-0.121, -0.090, -0.218, -0.786, -0.478, -0.299, -0.356, &
                 -0.630, -0.153, -0.490, -0.405, -0.356, -0.100, -5.171/

! Dry soil and ice density (kg/m**3)
! rogdr(1), rogdr(10), rogdr(11) are defined to obtain the values rog*cg of Pielke

      data rogdr/1.55e3, 1.57e3, 1.51e3, 1.37e3, 1.46e3, 1.55e3, 1.40e3, &
                 1.40e3, 1.53e3, 1.57e3, 1.60e3, 0.40e3, 1.80e3, 0.90e3/

! Dry soil heat capacity (J/kg/K)

      data cgdry/0.948e3, 0.898e3, 0.887e3, 0.927e3, 0.826e3, 0.761e3, 0.951e3, &
                 0.878e3, 0.771e3, 0.732e3, 0.681e3, 2.100e3, 0.950e3, 2.093e3/

! Saturated hydraulic conductivity (m/s) (revised values, feb. 2014)

      data ykghy/2.7e-5, 1.6e-5, 7.0e-6, 1.4e-6, 1.4e-6, 1.3e-6, 7.0e-7, &
                 7.5e-7, 4.9e-7, 4.7e-7, 4.8e-7, 1.4e-6, 1.0e-6, 0.0e0/

! Minimum volumetric water content of the soil (m**3/m**3)

      data qgmin/0.015, 0.023, 0.050, 0.090, 0.060, 0.085, 0.113, &
                 0.137, 0.134, 0.140, 0.145, 0.188, 0.010, 0.999/

! Maximum volumetric water content of the soil (m**3/m**3) (qgmax*row is soil porosity)

      data qgmax/0.395, 0.410, 0.435, 0.485, 0.451, 0.420, 0.477, &
                 0.476, 0.426, 0.492, 0.482, 0.600, 0.400, 1.000/

! Volumetric water content of the soil at wilting point (m**3/m**3)

      data qgwilt/0.0677, 0.0750, 0.1142, 0.1794, 0.1547, 0.1749, 0.2181, &
                  0.2498, 0.2193, 0.2832, 0.2864, 0.2744, 0.0200, 0.0000/

! Dry soil (qg=qgmin) emissivity in broadband window

      data dryemis1/0.940, 0.950, 0.960, 0.965, 0.975, 0.965, 0.970, &
                    0.960, 0.950, 0.955, 0.940, 0.965, 0.950, 0.997/

! Dry soil (qg=qgmin) emissivity in 8-12 micron window

      data dryemis2/0.860, 0.880, 0.950, 0.950, 0.960, 0.955, 0.965, &
                    0.945, 0.940, 0.945, 0.932, 0.960, 0.945, 0.997/

! Wet soil (qg=qgmax) emissivity in broadband window

      data wetemis1/0.970, 0.972, 0.985, 0.985, 0.995, 0.985, 0.990, &
                    0.985, 0.975, 0.978, 0.968, 0.988, 0.965, 0.997/

! Wet soil (qg=qgmax) emissivity in 8-12 micron window

      data wetemis2/0.890, 0.930, 0.955, 0.970, 0.985, 0.975, 0.980, &
                    0.965, 0.960, 0.965, 0.950, 0.982, 0.958, 0.997/

! Dry soil (qg=qgmin) albedo

      data dryalb/0.36, 0.30, 0.25, 0.15, 0.10, 0.18, 0.17, &
                  0.22, 0.25, 0.28, 0.35, 0.15, 0.18, 0.55/

! Wet soil (qg=qgmax) albedo

      data wetalb/0.22, 0.20, 0.18, 0.08, 0.05, 0.10, 0.13, &
                  0.15, 0.20, 0.15, 0.20, 0.05, 0.10, 0.55/

! Vegeta:
! 1  Evergreen needleleaf forest
! 2  Evergreen broadleaf forest
! 3  Deciduous needleleaf forest
! 4  Deciduous broadleaf forest
! 5  Mixed cover
! 6  Woodland
! 7  Wooded grassland
! 8  Closed shrubland
! 9  Open shrubland
! 10 Grassland
! 11 Cropland
! 12 Bare ground (including permanent ice shelf)
! 13 Urban and built-up

! Depth of plant's root (m)

      data root/1.00, 1.25, 1.00, 1.25, 0.80, 0.90, 0.70, &
                0.50, 0.50, 0.20, 0.40, 0.00, 0.00/

! Land roughness for momentum

      data vegroughm/1.65, 1.70, 1.65, 1.70, 1.55, 1.60, 0.90, &
                     0.21, 0.12, 0.14, 0.22, 0.05, 1.20/

! Land roughness for heat and moisture

      data vegroughq/0.90, 1.00, 0.90, 0.90, 0.80, 0.80, 0.50, &
                     0.11, 0.07, 0.08, 0.12, 0.03, 0.60/

! Vegetation emissivity in broadband window

      data vegemis1/0.997, 0.995, 0.991, 0.992, 0.990, 0.987, 0.984, &
                    0.986, 0.982, 0.985, 0.983, 0.963, 0.986/

! Vegetation emissivity in 8-12 micron window

      data vegemis2/0.993, 0.985, 0.983, 0.988, 0.982, 0.975, 0.975, &
                    0.975, 0.970, 0.978, 0.983, 0.958, 0.975/

! Vegetation albedo

      data vegalb/0.14, 0.13, 0.15, 0.15, 0.13, 0.15, 0.18, &
                  0.19, 0.20, 0.19, 0.19, 0.20, 0.15/

!********************************************************************************
! Fraction of soil types:
! 1 - sum(suolo(1:14)) + fmask = 1. is imposed
! 2 - if glaciers execeed bare soil, extra ice is converted into stone
!********************************************************************************

    open (17, file='soil12.bin', status='unknown', form='unformatted')
    do js = 1, 14
    do j = 1, jmsol
    read (17) isoilbase(:,j,js)
    enddo
    enddo
    close (17)

    suolo = 0.

    do jlat = 1, nlat
    do jlon = 1, nlon

    i0 = int((alont(jlon,jlat)+180.)*12.) +1   ! pixel containing the longitude alont
    j0 = int((90.-alatt(jlon,jlat))*12.)  +1   ! pixel containing the latitude alatt
    suolo(jlon,jlat,:) = float(isoilbase(i0,j0,:))*.01

    if (alatt(jlon,jlat).lt.-60.) then ! Antarctica
    suolo(jlon,jlat,:) = 0.
    suolo(jlon,jlat,14) = 1.-fmask(jlon,jlat)
    else
    sums = sum (suolo(jlon,jlat,1:14))
      if (sums.gt.0.) then
      eps = (1.-fmask(jlon,jlat))/sums
      suolo(jlon,jlat,1:14) = suolo(jlon,jlat,1:14)*eps
       if (suolo(jlon,jlat,14).gt.vegeta(jlon,jlat,12)) then
       drock = suolo(jlon,jlat,14) - vegeta(jlon,jlat,12)
       suolo(jlon,jlat,14) = suolo(jlon,jlat,14) - drock
       suolo(jlon,jlat,13) = suolo(jlon,jlat,13) + drock
       endif
      else
      suolo(jlon,jlat,5) = 1.-fmask(jlon,jlat)  ! defines soil type where not present: loam
      endif
    endif
    suolo(jlon,jlat,:) = max (suolo(jlon,jlat,:), 0.0)

    enddo
    enddo

#ifdef globo
    suolo(nlon,:,:) = suolo(2,     :,:)
    suolo(1,   :,:) = suolo(nlon-1,:,:)
    suolo(:,   1,1:13) = 0.
    suolo(:,   1,14  ) = 1. ! 100% glacier
    suolo(:,nlat,1:14) = 0.
#endif

!    print*, 'Suolo has been defined'

!********************************************************************************
! Definition of weighted soil and vegetation parameters (soilvegpar)
!********************************************************************************

    soilvegpar = 0.

    do jlat = 1, nlat
    do jlon = 1, nlon

    zsum = sum (suolo(jlon,jlat,1:13)) ! except glaciers
    zqgwilt = 0.
    if (zsum.gt.1.e-6) then
    do js = 1, 13
    soilvegpar(jlon,jlat,1 ) = soilvegpar(jlon,jlat,1 ) + psig0 (js)*suolo(jlon,jlat,js)
    soilvegpar(jlon,jlat,2 ) = soilvegpar(jlon,jlat,2 ) + rogdr (js)*suolo(jlon,jlat,js)*cgdry(js)
    soilvegpar(jlon,jlat,3 ) = soilvegpar(jlon,jlat,3 ) + ykghy (js)*suolo(jlon,jlat,js)
    soilvegpar(jlon,jlat,4 ) = soilvegpar(jlon,jlat,4 ) + qgmin (js)*suolo(jlon,jlat,js)
    soilvegpar(jlon,jlat,5 ) = soilvegpar(jlon,jlat,5 ) + qgmax (js)*suolo(jlon,jlat,js)
    zqgwilt                  = zqgwilt                  + qgwilt(js)*suolo(jlon,jlat,js)
    enddo
    soilvegpar(jlon,jlat,1:5) = soilvegpar(jlon,jlat,1:5)/zsum
    zqgwilt = zqgwilt/zsum

!  Exponent 'b' in soilvegpar(6) (from analytic expression)

    zrw = zqgwilt/soilvegpar(jlon,jlat,5)
    soilvegpar(jlon,jlat,6) = -log(-153./soilvegpar(jlon,jlat,1))/log(zrw) ! -153 = psig0*zrw**(-b)
     if(fmask(jlon,jlat).ge.0.5) then
     soilvegpar(jlon,jlat,4 ) = soilvegpar(jlon,jlat,4 )*(1.-fmask(jlon,jlat)) + fmask(jlon,jlat)
     soilvegpar(jlon,jlat,5 ) = soilvegpar(jlon,jlat,5 )*(1.-fmask(jlon,jlat)) + fmask(jlon,jlat)
     endif
    else
    soilvegpar(jlon,jlat,4 ) = 1.
    soilvegpar(jlon,jlat,5 ) = 1.
    endif

!  Albedo and emissivity over land

    zsum = sum (suolo(jlon,jlat,1:13))  ! except glaciers
    if (zsum.gt.1.e-6) then
    do js = 1, 13
    soilvegpar(jlon,jlat,7 ) = soilvegpar(jlon,jlat,7 ) + dryalb  (js)*suolo(jlon,jlat,js)
    soilvegpar(jlon,jlat,8 ) = soilvegpar(jlon,jlat,8 ) + wetalb  (js)*suolo(jlon,jlat,js)
    soilvegpar(jlon,jlat,9 ) = soilvegpar(jlon,jlat,9 ) + dryemis1(js)*suolo(jlon,jlat,js)
    soilvegpar(jlon,jlat,10) = soilvegpar(jlon,jlat,10) + dryemis2(js)*suolo(jlon,jlat,js)
    soilvegpar(jlon,jlat,11) = soilvegpar(jlon,jlat,11) + wetemis1(js)*suolo(jlon,jlat,js)
    soilvegpar(jlon,jlat,12) = soilvegpar(jlon,jlat,12) + wetemis2(js)*suolo(jlon,jlat,js)
    enddo
    soilvegpar(jlon,jlat,7:12) = soilvegpar(jlon,jlat,7:12)/zsum
    else ! for smoothing
    soilvegpar(jlon,jlat,7 ) = .26
    soilvegpar(jlon,jlat,8 ) = .20
    soilvegpar(jlon,jlat,9 ) = .960
    soilvegpar(jlon,jlat,10) = .940
    soilvegpar(jlon,jlat,11) = .980
    soilvegpar(jlon,jlat,12) = .960
    endif

!  Desert fraction in soilvegpar(13)

    if (alatt(jlon,jlat).gt.-50..and.alatt(jlon,jlat).lt.50.) then
     if (vegeta(jlon,jlat,9)+vegeta(jlon,jlat,12).gt..6) then  ! deserts in soilvegpar(13)
     if (suolo(jlon,jlat,14).lt.1.e-3) soilvegpar(jlon,jlat,13)=.4*vegeta(jlon,jlat,9)+.9*vegeta(jlon,jlat,12)
     soilvegpar(jlon,jlat,13) = soilvegpar(jlon,jlat,13)*(1.-.5*suolo(jlon,jlat,13))  ! reduction over rocks
     else
     soilvegpar(jlon,jlat,13) = 0.
     endif
    if (alatt(jlon,jlat).gt.-10..and.alatt(jlon,jlat).lt.10..and.alont(jlon,jlat).gt.90.) soilvegpar(jlon,jlat,13)=0.
    endif
    soilvegpar(jlon,jlat,13) = sqrt(soilvegpar(jlon,jlat,13))

!  Root depth in soilvegpar(14)

    zsum = sum (vegeta(jlon,jlat,1:11)) ! exept bare soil and urban
    if (zsum.gt.0.) then
    do jv = 1, 11
    soilvegpar(jlon,jlat,14) = soilvegpar(jlon,jlat,14) + root(jv)*vegeta(jlon,jlat,jv)
    enddo
    soilvegpar(jlon,jlat,14) = soilvegpar(jlon,jlat,14)/zsum
    endif

!  Momentum and temperature roughness length over land

    if (1.-fmask(jlon,jlat).gt.1.e-5) then  ! average roughness over all land fractions
    zrgm = 0.
    zrgq = 0.
    do jv = 1, 13
    zrgm = zrgm + vegroughm(jv)*vegeta(jlon,jlat,jv)
    zrgq = zrgq + vegroughq(jv)*vegeta(jlon,jlat,jv)
    enddo
    zrgm = zrgm/(1.-fmask(jlon,jlat))
    zrgq = zrgq/(1.-fmask(jlon,jlat))

    ordip = min (1.25e-5*(1.-fmask(jlon,jlat))*orogstd(jlon,jlat)**2, 6.)   ! Bolam and Moloch
    soilvegpar(jlon,jlat,15) = zrgm + ordip*.1     ! momentum roughness
    soilvegpar(jlon,jlat,16) = zrgq + ordip*.1     ! heat roughness
    else
    soilvegpar(jlon,jlat,15) = 1.e-5
    soilvegpar(jlon,jlat,16) = 1.e-8
    endif

!  Albedo and emissivity over land

    if (1.-fmask(jlon,jlat).gt.1.e-5) then  ! average over all land fractions
    do jv = 1, 13
    soilvegpar(jlon,jlat,17) = soilvegpar(jlon,jlat,17) + vegemis1(jv)*vegeta(jlon,jlat,jv)
    soilvegpar(jlon,jlat,18) = soilvegpar(jlon,jlat,18) + vegemis2(jv)*vegeta(jlon,jlat,jv)
    soilvegpar(jlon,jlat,19) = soilvegpar(jlon,jlat,19) + vegalb  (jv)*vegeta(jlon,jlat,jv)
    enddo
    soilvegpar(jlon,jlat,17) = soilvegpar(jlon,jlat,17)/(1.-fmask(jlon,jlat))
    soilvegpar(jlon,jlat,18) = soilvegpar(jlon,jlat,18)/(1.-fmask(jlon,jlat))
    soilvegpar(jlon,jlat,19) = soilvegpar(jlon,jlat,19)/(1.-fmask(jlon,jlat))
    else
    soilvegpar(jlon,jlat,17) = .990
    soilvegpar(jlon,jlat,18) = .980
    soilvegpar(jlon,jlat,19) = .15
    endif

    soilvegpar(jlon,jlat,20) = flake(jlon,jlat)
    soilvegpar(jlon,jlat,21) = max (orogstd(jlon,jlat)+vegeta(jlon,jlat,13)*15., 0.) ! increased over cities...

    enddo
    enddo

!  Smoothing of albedo and emissivity

    wei = .5
    call filt2 (soilvegpar(:,:,7 ), nlon, nlat, wei)
    call filt2 (soilvegpar(:,:,8 ), nlon, nlat, wei)
    call filt2 (soilvegpar(:,:,9 ), nlon, nlat, wei)
    call filt2 (soilvegpar(:,:,10), nlon, nlat, wei)
    call filt2 (soilvegpar(:,:,11), nlon, nlat, wei)
    call filt2 (soilvegpar(:,:,12), nlon, nlat, wei)
    call filt2 (soilvegpar(:,:,13), nlon, nlat, wei)
    call filt2 (soilvegpar(:,:,17), nlon, nlat, wei)
    call filt2 (soilvegpar(:,:,18), nlon, nlat, wei)
    call filt2 (soilvegpar(:,:,19), nlon, nlat, wei)

    print*, 'Soilvegpar defined in defsoil.'

    return
    end subroutine defsoil
!##################################################################################################
   subroutine polavert (p, nlon, nlat, npolcap)

!  Polar filter

   implicit none

   integer nlon, nlat, npolcap, jlon, jlat, nfilt, n
   real(4) pi, dlon, anu
   real(4) p(nlon,nlat), p1(nlon), zfiltt(nlat)

   pi = abs(acos(-1.))
   dlon = 360./float(nlon-2)

   do jlat = 1, npolcap+1

   call slofou1 (p(1:nlon,jlat), nlon, 3.14*(jlat-1))
   enddo
   do jlat = nlat-npolcap, nlat
   call slofou1 (p(1:nlon,jlat), nlon, 3.14*(nlat-jlat))
   enddo

   zfiltt = 0
   do jlat = 2, nlat/2
   zfiltt(jlat) = 4./(dlon*pi/180.)**2/(3.14*(jlat-1))**2
   zfiltt(jlat) = max (zfiltt(jlat)-.2, 0.)
   enddo
   do jlat = 1, nlat/2-1
   zfiltt(nlat-jlat) = zfiltt(jlat+1)
   enddo
   zfiltt(1:npolcap+1) = 0
   zfiltt(nlat-npolcap:nlat) = 0

   do jlat = 1, nlat
   nfilt = zfiltt(jlat) + .99999
   anu = 0.
   if (nfilt > 0) anu = zfiltt(jlat)/float(nfilt)

   do n = 1, nfilt
     do jlon = 2, nlon-1
     p1(jlon) = .25*(p(jlon-1,jlat)+p(jlon+1,jlat))-.5*p(jlon,jlat)
     enddo
   p1(1   ) = p1(nlon-1)
   p1(nlon) = p1(2     )
   p(1:nlon,jlat) = p(1:nlon,jlat) + anu*p1(1:nlon)
   enddo

   enddo

   return
   end subroutine polavert
!##################################################################################################
   subroutine slofou1 (p, nlon, nt)

   implicit none

   integer jlon, ntot, n, nfou, nlon
   real(4) nt, dlon
   real(4), dimension(nlon)   :: p, zps, psmoo
   real(4), dimension(nlon,nlon) :: snt, cst
   real(4), dimension(0:nlon) :: zpc
   real(8) zdpi, zxt

   dlon = 360./float(nlon-2)
   zdpi = dabs(dacos(-1.d0))
   nfou = (nlon-2)/2
   do n = 1, nfou
     do jlon = 1, nlon
     zxt = (jlon-1) * 360.d0/dfloat(nlon-2) * zdpi/180.d0
     snt(jlon,n) = dsin(zxt*n)
     cst(jlon,n) = dcos(zxt*n)
     enddo
   enddo

   ntot = nint(3.*nt)
   do n = 1, ntot
   psmoo(n) = exp(-float(n)**2/nt**2)*dlon/180.
   enddo

   zpc(0) = 0.
   do jlon = 2, nlon-1
   zpc(0) = zpc(0) + p(jlon)
   enddo
   zpc(0) = zpc(0)/float(nlon-2)

   do 10 n = 1, ntot
   zps(n) = 0.
   zpc(n) = 0.
     do jlon = 2 , nlon-1
     zps(n) = zps(n)+p(jlon)*snt(jlon,n)
     zpc(n) = zpc(n)+p(jlon)*cst(jlon,n)
     enddo
   zps(n) = zps(n)*psmoo(n)
   zpc(n) = zpc(n)*psmoo(n)
10 continue

   p(1:nlon) = zpc(0)
   do n = 1, ntot
   p(1:nlon) = p(1:nlon) + zps(n)*snt(1:nlon,n) + zpc(n)*cst(1:nlon,n)
   enddo

  return
  end subroutine slofou1
!#################################################################################################################
#ifdef ncarg
  subroutine bolcolor(a,lsm,lambda,theta,im,jm,dlatr,ris,dis)

! Plotta topografia (con param. di plot definiti sotto), lsm (linea di costa),
! meridiani e paralleli (derivata da bolcolor.f, semplificato)

      real a(im,jm),lsm(im,jm),lambda(im,jm),theta(im,jm)
      character*90 title
      character*4 llbs(21)
      dimension lind(100)

  parameter (lmap=10000000,lwrk=100000)
  external fill
  external cpdrpl
  real rwrk(lwrk),xcra(lwrk),ycra(lwrk)
  integer iwrk(lwrk),iama(lmap),iarea(1000),igrp(1000)

  nlonr = im
  nlatr = jm

  do i=1,100
  lind(i)=i
  enddo

! open gks - ncar graphics

  call gopks(6,idum)
  call gopwk(1,2,1)
  call gacwk(1)
  call gsclip(0)

! character data base selection for contour labels

  call gstxfp (-2,2) !-2=roman  2=max precisione

  ioffm=1
  ioffd=0

  xlt=.035
  ybt=.07  ! definisce l'altezza della base del grafico
  side=.81 ! lato del grafico

! RIS: reference isoline; DIS: intervallo tra isolinee;
! TITLE: titolo del grafico;
! NCOLOR: no. colori nella scala (deve essere pari);
! NOFFS: offset scala colori;
! INDCOL: indice per scelta scala colori (arg. della routine COLOR);
! ILAB1: per labels: se 2, limita cifre dec. nelle labels max. e min.;
!   se 3, una cifra dec. nei max. e min.; se 4, fa scrivere metri nella
!   scala colori orog.);
! ILAB2: se 0, no line labels ma con max. e min; altrimenti
!   2 o meglio 3: def. densità line labels;
! IMAX: se 2, mette max. e min. senza 'H' e 'L'; se 3 mette anche
!   'H' e 'L'; se 1, nessuna indic. max. e min.

! Qui def. i valori come in boldis.def (per orografia) - ris e dis sono passati in argom.

!  ris = 1755.
!  dis = 250.
!!!!!!!  title = "MODEL OROGRAPHY"
  title = "               "
!  ncolor = 16
  ncolor = 14
  offs =  -5.
  indcol = 2
  ilab1 =  4
  ilab2 =  0
  imax =   2

! COLORI

! Il no. colori va messo PARI e non comprende i due colori
! aggiuntivi di background settati nella subr. COLOR (nero e grigio)
! che sono di indici 1 e 2 (piu' il bianco di indice 0)
! i colori dipendono dal valore di ris def. per ogni quantita'

! OFFS: valore HUE del primo colore nella scala (arg. routine COLOR)
! INDCOL: indice per scelta scala colori (arg. routine COLOR)
!  0: default;
! -1: scala rovesciata, dal blu al rosso (per T, THETA, CAPE, PRECIP., SNOW, PV)
!  2: scala colori speciali per orografia
!  1: colori speciali  per nubi

  call gsfais (1)

!  coordinates of the plot in the page

  call cpsetr('VPL - Viewport Left',   xlt     )
  call cpsetr('VPR - Viewport Right',  xlt+side)
  call cpsetr('VPB - Viewport Bottom', ybt     )
  call cpsetr('VPT - Viewport Top',    ybt+side)

! notare che la seg. istruz. influisce sul no. di cifre di tutte le labels,
! comprese quelle della barra colori - controllare che non vi siano arrotondam.

  call cpseti('NSD', ILAB1)

!  levels are chosen by the user

  call cpseti('CLS - Contour Level Selection',0)
  call cpseti('NCL - Number of Contour Levels', ncolor)
  call cpsetc('ILT - INFORMATION LABEL TEXT', &
              '                                          ')
  call gsplci(1) ! colore per linee e numeri (1:nero)

  call color(ncolor,indcol,offs)

  call cprect(a,im,nlonr,nlatr,rwrk,lwrk,iwrk,lwrk)

  if(imax.eq.3) then ! high lows only on MSLP graph
    call cpsetc ('HLT','H''L') ! only L and H without values
    call cpsetr ('HLS' ,.020) ! size of highs and lows

! labels di max e min sulle seg. quantita'

  elseif(imax.eq.2) then
    call cpsetc ('HLT','$ZDV$') ! max. and min. values (but without H and L)
    call cpsetr ('HLS' ,.013) ! size of high and low values
  else
    call cpsetc ('HLT',' ') ! no high and low values nor labels H, L
  endif

  do i=1,ncolor

! con la scelta seg., la ris corrisponde alla meta' della scala col.
! tenendo conto anche che i primi due colori sono di background
! deve essere coerente con la def. livelli per barra laterale

  zlev=ris+(i-1-ncolor/2.)*dis
  call cpseti('PAI - Parameter Array Index',i)
  call cpsetr('CLV - Contour Level Value', zlev)
  call cpseti('CLU - CONTOUR LEVEL USE FLAGS',3)
  call cpseti('LLC - LINE LABEL COLOR INDEX',1)

! la seg. setta densità labels (2 o 3, meglio 3)
! elim. line labels su alcune quantita' senza elim.labels massimi e minimi

  if (ilab2.ne.0) then
    call cpseti('LLP - LINE LABEL POSITIONING FLAG',ILAB2)! label posit.
  else
    call cpseti('LLP',ILAB2) ! no line labels
  endif

  call cpsetr('LLS - Line Label Size',0.014)! size of std. label digits
  if(zlev.lt.0.0) then !dashed lines for neg. values
   call cpsetc ('CLD','$''''$''''$''''$''''$')
  endif
  enddo

! Initialize Areas

  call arinam(iama,lmap)

! Add contours to area map

  call cpclam(a,rwrk,iwrk,iama)
  call arscam(iama,xcra,ycra,lwrk,iarea,igrp,40,fill) ! fill colori per tutto tranne il vento
  call sfseti('TY - TYPE OF FILL',0) ! resets normal type of fills (colors)

! Add label boxes to area map
! (IF to exclude ALL labels for certain fields)

!!!!  call cplbam(a,rwrk,iwrk,iama)

! Draw all labels

  call cplbdr (a,rwrk,iwrk) !draw labels

! Draw contours, masking label boxes

!!!!  call cpcldm(a,rwrk,iwrk,iama,cpdrpl)

  ilab=1
  call setusv('LW',1300)
  ilab=0
  call gslwsc(2.) ! bold coastlines
  call gscr(1,1, 0.25, 0.25, 0.25) ! coastline color
  call onecont(lsm,im,nlonr,nlatr,0.5)
  call gscr(1,1,.0, 0.,.0) ! back to black line color
  call gslwsc(1.) ! back to thin lines

! automatizzazione interv. in gradi del plot di merid. e paralleli
! in funz. della larghezza in long. del dominio
! non tutte le longit. possono essere plottate con interv. di 5 gradi

  dint=5.
  elon=maxval(lambda) - minval(lambda)
  if(elon.lt.15.) dint=2.5
  if(elon.gt.90.) dint=10.
  alamin=-180.
  alamax=360.
  call cocont(lambda,im,nlonr,nlatr,alamin,alamax,dint)
  call cocont(theta,im,nlonr,nlatr,-90.,90.,dint)
!  call perim(1,nlonr-1,1,nlatr-1)
  call perim(1,1,1,1)
  call cprect(a,im,nlonr,nlatr,rwrk,lwrk,iwrk,lwrk) ! necessary to avoid problem with labels of label bar

  call set(0.,1.,0.,1.,0.,1.,0.,1.,1)

!  call pwritx(xlt,.97,title,40,25,0,-1)
  call pwritx(xlt,.94,title,40,25,0,-1)
!!!!  call pwritx(xlt,.05,'CNR-ISAC MODEL, BOLOGNA',23,12,0,-1)
  call pwritx(xlt,.05,'                       ',23,12,0,-1)

! Draw a label bar for the plot, relating colors to values.

   ibarlvl=ncolor
   ibarlvl2=ibarlvl/2
   do i=1,ibarlvl+1
    call cpsetr ('ZDV - Z DATA VALUE',ris+(I-1-IBARLVL2)*dis)
    call cpgetc ('ZDV - Z DATA VALUE',LLBS(I))
   enddo
   call lbseti ('CBL - COLOR OF BOX LINES',0) ! colore che separa le box
   call gsplci(17) !set black for bar numbers

! nell'arg. si parte da LIND(3) per evitare di mettere i due colori
! di background (nero e grigio) di indici 1 e 2

  call lblbar (1,.873,0.988,0.06,.9,ibarlvl,.3,1.,lind(3),0,llbs,ibarlvl+1,1)

! close gks - ncar graphics

  call gdawk(1)
  call gclwk(1)
  call gclks

      return
 1010 format(11f7.1)
 8000 format(1h+,6x,a40,i8,i6,i6)
 9100 format(20i4)
 9200 format(2i5,4f11.5)
      END

      subroutine fill (xcra,ycra,ncra,iarea,igrp,ngrps)

      real xcra(ncra),ycra(ncra)
      integer iarea(ngrps),igrp(ngrps)

! Get area identifiers for contour levels and vertical strips.

      ifill=0
      do 101 i=1,ngrps
      if (igrp(i).eq.3) ifill=iarea(i)
  101 continue
          if (ifill.gt.0) then
          call gsfaci(ifill+1)
          call gfa(ncra,xcra,ycra)
          endif
      return
      end

      SUBROUTINE COLOR(N,INDCOL,OFFS)

! INDCOL: indice introdotto per avere diverse scale di colori
! INDCOL = 0: default
! INDCOL =-1: scala colori colori rovesciata
! INDCOL = 1: scala colori spec. per nubi
! INDCOL = 2: scala colori spec. (e rovesciata) per orograf.

! notare che i colori di indici 1 e 2 sono il nero e grigio settati sotto
! mentre i rimanenti colori (in numero N) hanno indici da 2 a n+2!

! BACKGROUND COLOR WHITE

!      CALL GSCR(1,0,1.,1.,1.)
      CALL GSCR(1,0,1.,.96,.76)  ! cambiato in crema

! First foreground color is black

      CALL GSCR(1,1,0.,0.,0.)

! Second foreground color is gray

!      CALL GSCR(1,2,.75,.75,.75)
      CALL GSCR(1,2,.94,.94,.94)   ! grigio chiaro

! Choose other foreground colors spaced equally around the spectrum

      ICNT=0

      IF(INDCOL.EQ.1) GOTO 20

!      HUER = 360. ! usa tutta la scala colori, ma i col. estremi si assomigliano
      HUER = 300. ! rende piu' distinti i col. estremi

      HUES=HUER/N

      REDLN=OFFS ! def. offset scala colori (colore di partenza)

      DO 10, I=1,N
      IF(INDCOL.EQ.-1.OR.INDCOL.EQ.2) THEN
      II=N-I+1 ! caso scala colori rovesciata
      ELSE
      II=I ! caso normale
      ENDIF

      IF(INDCOL.EQ.-1) CALL GSCR(1, 90, 0., 0., 0.35) ! blu scuro, usato per streamlines

! introd. scala non lineare di colori (parte in sin)
! (comprime i verdi, amplia arancio e ciano, comprime bordi scala)

      XHUE=II*HUES+REDLN-7.*sin((II*HUES+REDLN)*3.14159*4./HUER)
      IF(XHUE.GT.360.) XHUE = XHUE-360.
      IF(XHUE.LT.0.) XHUE = XHUE+360.

! convert hue, lightness and saturation to RGB

      CALL HLSRGB(XHUE,50.,84.,RED,GREEN,BLUE) ! il 2o e 3o param. sono luminosità e saturazione

! Sort colors so that the redest is first, and violetest is last
! Si saltano i primi due indici gia' assegnati

        CALL GSCR(1,I+2,RED,GREEN,BLUE)
 10   CONTINUE

! definiz. scale speciali (anche solo parziali)

 20   continue

      if(INDCOL.eq.1) then

! le seg. def. le nubi partendo da una base azzurra (sereno)
! e aggiungendo bianco al crescere delle nubi
! si def. solo gli indici che servono per le nubi

      CALL GSCR(1, 3,  0.13, 0.36, 0.85)
      CALL GSCR(1, 4,  0.18, 0.40, 0.85)
      CALL GSCR(1, 5,  0.23, 0.43, 0.85)
      CALL GSCR(1, 6,  0.27, 0.46, 0.85)
      CALL GSCR(1, 7,  0.31, 0.49, 0.85)
      CALL GSCR(1, 8,  0.35, 0.52, 0.85)
      CALL GSCR(1, 9,  0.39, 0.55, 0.85)
      CALL GSCR(1, 10, 0.43, 0.58, 0.86)
      CALL GSCR(1, 11, 0.47, 0.61, 0.87)
      CALL GSCR(1, 12, 0.51, 0.64, 0.88)
      CALL GSCR(1, 13, 0.55, 0.67, 0.89)
      CALL GSCR(1, 14, 0.59, 0.70, 0.90)
      CALL GSCR(1, 15, 0.63, 0.73, 0.91)
      CALL GSCR(1, 16, 0.67, 0.76, 0.92)
      CALL GSCR(1, 17, 0.71, 0.79, 0.93)
      CALL GSCR(1, 18, 0.75, 0.82, 0.94)
      CALL GSCR(1, 19, 0.79, 0.85, 0.95)
      CALL GSCR(1, 20, 0.83, 0.88, 0.96)
      CALL GSCR(1, 21, 0.88, 0.91, 0.97)
      CALL GSCR(1, 22, 0.94, 0.95, 0.98)

      elseif(INDCOL.eq.2) then ! ridef. alcuni colori per orografia (mare, neve...)

      CALL GSCR(1,  3, 0.25, 0.55, 0.95)
      CALL GSCR(1, 12, 0.75, 0.35, 0.10)
      CALL GSCR(1, 13, 0.96, 0.50, 0.35)
      CALL GSCR(1, 14, 0.93, 0.75, 0.95)
      CALL GSCR(1, 15, 0.78, 0.78, 0.96)
      CALL GSCR(1, 16, 0.82, 0.85, 0.96)
      CALL GSCR(1, 17, 0.88, 0.90, 0.98)
      CALL GSCR(1, 18, 0.94, 0.96, 1.00)
      endif

      RETURN
      END

      subroutine onecont(a,im,m,n,cval)

! Plots a single contour line of value cval

      parameter(ndimr=100000,ndimi=100000)
      dimension a(im,n)
      dimension rwk(ndimr),iwk(ndimi)

! initialize plot

      call cprect(a,im,m,n,rwk,ndimr,iwk,ndimi)

! various flags

      call cpsetc('ILT',' ')
      call cpsetc('HIT',' ')
      call cpsetc('LOT',' ')

! contour level selection

      call cpseti('CLS',0)
      call cpseti('NCL',1)
      call cpseti('PAI',1)
      call cpsetr('CLV',cval)  ! Value of the contour line to be plotted
      call cpseti('CLU',1)

! draw contours

      call cpcldr(a,rwk,iwk)

      return
      end

      subroutine cocont(a,im,m,n,zmin,zmax,zint)

! Plots contour lines for values set in cclv

      parameter(ndimr=100000,ndimi=100000)
      dimension a(m,n)
      dimension rwk(ndimr),iwk(ndimi)

! initialize plot

      call cprect(a,im,m,n,rwk,ndimr,iwk,ndimi)

! various flags

      call cpsetc('ILT',' ')
      call cpsetc('HIT',' ')
      call cpsetc('LOT',' ')
      call cpseti('CLS',0)

! Definition of number of contour lines

      ncl= (zmax-zmin)/zint + 1

      call cpseti('NCL',ncl)

      do inn=1,ncl
      call cpseti('PAI',inn)
      cclv=zmin+(inn-1)*zint
      call cpsetr('CLV',cclv)
      call cpseti('CLU',1)

!      if(cclv.lt.0.)then                           ! dashed lines for negative values
!       call cpsetc('CLD','$''''''$''''''$''''''$') ! dashed lines pattern
!       call cpsetc('CLD','$''$''$''$''$''$')       ! dashed lines pattern
!      endif
! CLD is set as decimal value corresponding to a 16 bit binary no.
! (for ex. 1010101010101010=43690;

      call cpseti('CLD',43690)                           ! sets dashed lines pattern

      enddo

! draw contours

      call cpcldr(a,rwk,iwk)
      return
      end
#endif
