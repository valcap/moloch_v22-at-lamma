! Last update 26/01/2019

! Subroutines common to various bolam and moloch preprocessing and
! postprocessing codes.
!=======================================================================
subroutine qsat_tetens(t, p, eps, qs, qsw, qsi, es, esw, esi)

! Computes saturation specific humidity QS and water vapour sat. partial pressure ES
! following Tetens formula for saturation over liquid water and over ice,
! and formula derived from Mazin-Hargian Handbook for saturation over
! a mixture of liquid water and ice

 implicit none

 real :: t, p, eps, qs, qsw, qsi, es, esw, esi, t0=273.15, fracw, z

 esw = 611.*exp(17.40*(t-t0)/(t-33.65))
 esi = 611.*exp(22.45*(t-t0)/(t-0.75))

 z = eps*esw/(p-esw)
 qsw = z/(1.+z)
 z = eps*esi/(p-esi)
 qsi = z/(1.+z)

 fracw = 1.04979*(0.5 + 0.5*tanh((t-t0+9.)/6.))  ! Fraction of liquid water in mixed-phase clouds
 if(t.ge.t0) fracw = 1.
 if(t.lt.245.) fracw = 0.

 qs  = fracw*qsw+(1.-fracw)*qsi
 es  = fracw*esw+(1.-fracw)*esi

return
end subroutine qsat_tetens
!=======================================================================
subroutine td_tetens(t, p, q, eps, td)

! Computes dew point temperature using
! inverted Tetens formula for saturation over liquid water and over ice

 implicit none

 real :: t, p, q, eps, td, t0=273.15, ep

 ep = p*q/(eps+q*(1.-eps))

 if (t >= t0) then
   td = (t0-33.65/17.40*log(ep/611.))/(1.-1./17.40*log(ep/611.))
 else
   td = (t0- 0.75/22.45*log(ep/611.))/(1.-1./22.45*log(ep/611.))
 endif

return
end subroutine td_tetens
!=======================================================================
subroutine comp_esk(esat, qsat, t, p, iflag)

! Computes esat from temperature and qsat from absolute temperature and pressure
! IFLAG 1: esat and qsat with respect to water and ice, separately, depending if t>tzer or t<tzer
!       2: esat and qsat with an interpolation at t<tzer between water and ice
!       3: esat and qsat with respect to water also for t<tzer

  tzer = 273.15
  ezer = 611.
  cpv  = 1869.46
  cw   = 4186.8
  rd   = 287.05
  rv   = 461.51
  yliv = 2834170.5
  yliw = 333560.5
  eps  = rd/rv
  ci   = cw/2.
  ylwv = yliv-yliw
  ccw1 = (cpv-cw)/rv
  ccw2 = ylwv/tzer/rv-ccw1
  cci1 = (cpv-ci)/rv
  cci2 = yliv/tzer/rv-cci1

  zt0t = tzer/t
  if (zt0t.le.1.) then
   zesk = ezer*exp(-ccw1*alog(zt0t)+ccw2*(1.-zt0t))  ! Partial pressure over water
   else
    if (iflag.eq.1) then
    zesk = ezer*exp(-cci1*alog(zt0t)+cci2*(1.-zt0t)) ! Partial pressure over ice
    elseif (iflag.eq.2) then
    zratio = 1.04979*(0.5 + 0.5*tanh((t-tzer+9.)/6.))
    zesk = zratio*(ezer*exp(-ccw1*alog(zt0t)+ccw2*(1.-zt0t)))+  &
      (1.-zratio)*(ezer*exp(-cci1*alog(zt0t)+cci2*(1.-zt0t)))
    elseif (iflag.eq.3) then
    zesk = ezer*exp(-ccw1*alog(zt0t)+ccw2*(1.-zt0t))
    else
    print*, "IFLAG out of range in subroutine comp_esk", iflag
    stop
    endif
  endif

  esat=zesk
  qsat = zesk*eps/(p+zesk*(eps-1.))
  return
end subroutine comp_esk
!=======================================================================
subroutine qsat_entropy(t, p, qs, qsw, qsi, es, esw, esi, fracw)

! Computes saturation specific humidity QS and water vapour sat. partial pressure ES
! using an analytic formula (based entropy) for saturation over liquid water and ice,
! and a formula derived from Mazin-Hargian Handbook for saturation over
! a mixture liquid water and ice

IMPLICIT NONE

real :: t, p, qs, qsw, qsi, es, esw, esi, fracw, zt0t, z
real, parameter :: tzer=273.15, ezer=611., rd=287.05, rv=461.51, eps=rd/rv,       &
 cpv=1869.46, cw=4186.8, ci=cw/2., yliv=2834170.5, yliw=333560.5, ylwv=yliv-yliw, &
 ccw1=(cpv-cw)/rv, ccw2=ylwv/(tzer*rv)-ccw1, cci1=(cpv-ci)/rv, cci2=yliv/(tzer*rv)-cci1

 zt0t=tzer/t
 esw=ezer*exp(-ccw1*alog(zt0t)+ccw2*(1.-zt0t)) ! Partial pressure over water
 esi=ezer*exp(-cci1*alog(zt0t)+cci2*(1.-zt0t)) ! Partial pressure over ice

 z = eps*esw/(p-esw)
 qsw = z/(1.+z)
 z = eps*esi/(p-esi)
 qsi = z/(1.+z)

 fracw = 1.04979*(0.5 + 0.5*tanh((t-tzer+9.)/6.)) ! Fraction of liquid water in mixed-phase clouds
 if(t.ge.tzer) fracw = 1.
 if(t.lt.245.) fracw = 0.

 qs  = fracw*qsw+(1.-fracw)*qsi
 es  = fracw*esw+(1.-fracw)*esi

end subroutine qsat_entropy
!=======================================================================
subroutine controlq(q, t, p, eps)

! Checks that Q does not exceed (significantly) the saturation value, using
! a saturated vapour function, mixed for T<0, following Mazin-Hargian Handbook.
! An over-saturation allowed value is prescribed to depend on
! atmospheric pressure (larger in the mid-low troposphere), in order to mimic the
! presence of some cloud water/ice content
! To be used only if cloud water and ice are missing in input!

 implicit none

 real :: q, t, p, e, eps, argp, fp, rh, rhmax, qsat, qsatw, qsati, esat, esatw, esati

 argp=((p-650.e2)/300.e2)**2
 fp=exp(-argp)
 rhmax=1.+0.04*fp

 q = max(q,1.e-12)
 e = (q*p)/(eps+q-eps*q)

 call qsat_tetens(t, p, eps, qsat, qsatw, qsati, esat, esatw, esati)

 e = min(e, esat*rhmax)
 q = e*eps/(p+e*(eps-1.))

! Increase of q values close to saturation (the log argument is the
! factor by which q is increased for rh=1.)

 rh=e/esat
 q=q*exp(alog(1.032)*rh**12)

return
end subroutine controlq
!=======================================================================
subroutine rot_grid(x0,y0,x00,y00,dlon,dlat,alon,alat,nlon,nlat)

! Calculation of standard geographical coordinates (ALON, ALAT) of model
! (rotated) grid points (rotation centre: X0, Y0)
! X00, Y00 are the SW corner point coordinates in rotated coordinates
! Some computations require double precision

implicit none

real :: x0, y0, x00, y00, dlon, dlat
integer :: nlon, nlat, jlat, jlon
real, dimension(nlon,nlat) :: alat, alon
real*8 :: zfac, zx0, zy0, zlon, zlat, zzlat, zaarg, zarg

if (abs(x0)>0.01.or.abs(y0)>0.01) then

! Case of rotated grid

  zfac = dabs(dacos(-1.d0))/180.d0
  zx0  = dble(x0)*zfac
  zy0  = dble(y0)*zfac
  do jlat=1,nlat
    zlat=(dble(y00)+(jlat-1)*dble(dlat))*zfac
    do jlon=1,nlon
      zlon = (dble(x00)+(jlon-1)*dble(dlon))*zfac
      zzlat= 1.d0/zfac*dasin( dcos(zy0)*dsin(zlat) +  &
                       dsin(zy0)*dcos(zlat)*dcos(zlon) )
      zarg = -dsin(zlat)*dsin(zy0)+dcos(zy0)*dcos(zlat)*dcos(zlon)
      zaarg = zarg/dcos(zfac*zzlat)
      alat(jlon,jlat)=zzlat
      if(zaarg < -1.d0.and.zaarg > -1.00001d0) zaarg = -1.d0
      if(zaarg >  1.d0.and.zaarg <  1.00001d0) zaarg =  1.d0

      if (abs(abs(alat(jlon,jlat))-90.) > 1.e-7) then
        if (zlon < 0.d0) then
          alon(jlon,jlat) = 1.d0/zfac*(zx0-dacos(zaarg))
        else
          alon(jlon,jlat) = 1.d0/zfac*(zx0+dacos(zaarg))
        endif
      else
        alon(jlon,jlat) = 0.
      endif

      if (alon(jlon,jlat) >  180.) alon(jlon,jlat) = alon(jlon,jlat) - 360.
      if (alon(jlon,jlat) < -180.) alon(jlon,jlat) = alon(jlon,jlat) + 360.
    enddo
  enddo

else

! Case of non rotated grid

  do jlon=1,nlon
    alon(jlon,:)=x00+(jlon-1)*dlon
  enddo
  do jlat=1,nlat
    alat(:,jlat)=y00+(jlat-1)*dlat
  enddo

endif

return
end subroutine rot_grid
!=======================================================================
subroutine anti_rot_grid(x0,y0,alon,alat,xr,yr,nlon,nlat)

! Calculation of the rotated coordinates (put in XR, YR) with rotation centre X0, Y0
! (given in geographical coordinates) for input grid points given in
! geographical coordinates (defined in ALON, ALAT)
! Some computations require double precision

implicit none

real :: x0, y0
integer :: nlon, nlat, jlat, jlon
real, dimension(nlon,nlat) :: alat, alon, xr, yr
real*8 :: zfac, zx0, zy0, zx, zy, zlon, zlat, zz, pi

if (abs(x0)>0.01.or.abs(y0)>0.01) then

! Case of rotated grid

  pi  = dabs(dacos(-1.d0))
  zfac = pi/180.d0
  zx0  = dble(x0)*zfac
  zy0  = dble(y0)*zfac
  do jlat=1,nlat
  do jlon=1,nlon
    if (jlon==1.and.jlat==1) call sleep (0) ! fix for a problem of ifort 14 (but needs compil. with -O2 or less)
    zx = dble(alon(jlon,jlat))*zfac
    if (zx-zx0 > pi) zx = zx - 2.d0*pi
    if (zx-zx0 < -pi) zx = zx + 2.d0*pi
    zy = dble(alat(jlon,jlat))*zfac
    zlat = dasin( -dcos(zy)*dsin(zy0)*dcos(zx-zx0)+dsin(zy)*dcos(zy0) )
    zz = dsin(zy)-dcos(zy0)*dsin(zlat)
    zz = zz/(dsin(zy0)*dcos(zlat))
    if (zz < -1.d0.and.zz > -1.00001d0) zz = -1.d0
    if (zz > 1.d0.and.zz < 1.00001d0) zz = 1.d0
    if (zx < zx0) then
      zlon = -dacos(zz)
    else
      zlon = dacos(zz)
    endif
    zx = zlon/zfac
    zy = zlat/zfac
    xr(jlon,jlat) = sngl(zx)
    yr(jlon,jlat) = sngl(zy)
  enddo
  enddo

else

! Case of non rotated grid

  xr(:,:) = alon(:,:)
  yr(:,:) = alat(:,:)

endif

return
end subroutine anti_rot_grid
!=======================================================================
 subroutine rot_wind(upu,vpu,upv,vpv,xxu,xxv,yyu,yyv,upup,vpvp,nlon,nlat,x0,y0,x00,y00,dlon,dlat)

! Rotation of wind components on the rotated grid
! Some computations require double precision

 implicit none

 real*8 :: pi,fac,x0dr,y0dr,x00dr,y00dr,dlonr,dlatr,zdxxu,zdyyu,zdxur,zdyur,zdxxv,zdyyv,zdxvr,zdyvr

 real, dimension(nlon,nlat) :: upu, vpu, upv, vpv, xxu, xxv, yyu, yyv, upup, vpvp, upvp, vpup
 real(8), dimension(nlon) :: xur, xvr
 real(8), dimension(nlat) :: yur, yvr

 real :: x0, y0, x00, y00, dlon, dlat

 integer :: nlon, nlat, i, j, icount

! Winds before rotat.:  UPU, VPU: wind compon. at "u" points
!                       UPV, VPV: wind compon. at "v" points
! Winds after rotat.: UPUP, VPUP: wind compon. at "u" points
!                     UPVP, VPVP: wind compon. at "v" points
! NOTE: comput. of VPUP and UPVP only for verification - such quantities are not used

! Case of no rotation (Y0 = 0.)

 save icount
 if (abs(y0).lt.1.e-8) then
 icount=icount+1
 if(icount.le.2) write(*,*)"No rotation of grid coordinates!"
     do j=1,nlat
       do i=1,nlon
       upup(i,j) = upu(i,j)
       vpup(i,j) = vpu(i,j)
       upvp(i,j) = upv(i,j)
       vpvp(i,j) = vpv(i,j)
       enddo
     enddo
    return
 endif

 pi = dabs(dacos(-1.d0))
 fac= pi/180.d0

 x0dr=x0*fac
 y0dr=y0*fac
 x00dr=x00*fac
 y00dr=y00*fac
 dlonr=dlon*fac
 dlatr=dlat*fac

 do i=1,nlon
 xur(i)=x00dr+(i-1)*dlonr+dlonr/2.d0
 xvr(i)=x00dr+(i-1)*dlonr
 enddo

 do j=1,nlat
 yur(j)=y00dr+(j-1)*dlatr+dlatr/2.d0
 yvr(j)=y00dr+(j-1)*dlatr
 enddo

  do j=1,nlat
   do i=1,nlon

    zdxxu = dble(xxu(i,j))*fac-x0dr
    if (zdxxu.lt.-pi) zdxxu = zdxxu + 2.d0*pi
    if (zdxxu.gt. pi) zdxxu = zdxxu - 2.d0*pi

    zdyyu=yyu(i,j)
    zdxur=xur(i)
    zdyur=yur(j)

    zdxxv = dble(xxv(i,j))*fac-x0dr
    if (zdxxv.lt.-pi) zdxxv = zdxxv + 2.d0*pi
    if (zdxxv.gt. pi) zdxxv = zdxxv - 2.d0*pi

    zdyyv=yyv(i,j)
    zdxvr=xvr(i)
    zdyvr=yvr(j)

!  Distinction between cases in which longitude (in "true" coordinates) is 0 or 180.
!  Change of sign of u must be taken into account across the pole.

  if (abs(xur(i)).lt.1.d-4) then
  vpup(i,j) = vpu(i,j)
  upup(i,j) = upu(i,j)
  if (zdxxu.gt.pi/2.d0.or.zdxxu.lt.-pi/2.d0) upup(i,j) = -upu(i,j)
  if (zdxxu.gt.pi/2.d0.or.zdxxu.lt.-pi/2.d0) vpup(i,j) = -vpu(i,j)
  else
  vpup(i,j)=dsin(zdxxu)*dsin(y0dr)/dcos(zdyur)*upu(i,j)+(dcos(zdxxu) &
            *dsin(zdyyu*fac)*dsin(y0dr)/dcos(zdyur)+dcos(y0dr)*      &
            dcos(zdyyu*fac)/dcos(zdyur))*vpu(i,j)
  upup(i,j)=(dcos(y0dr)*dcos(zdyur)-dsin(y0dr)*dsin(zdyur)*dcos(zdxur))/(dsin(y0dr)* &
             dsin(zdxur))*vpup(i,j)-dcos(zdyyu*fac)/(dsin(y0dr)*dsin(zdxur))*vpu(i,j)
  endif

  if (abs(xvr(i)).lt.1.d-4) then
  vpvp(i,j)=vpv(i,j)
  upvp(i,j)=upv(i,j)
  if (zdxxv.gt.pi/2.d0.or.zdxxv.lt.-pi/2.d0) upvp(i,j)=-upv(i,j)
  if (zdxxv.gt.pi/2.d0.or.zdxxv.lt.-pi/2.d0) vpvp(i,j)=-vpv(i,j)
  else
  vpvp(i,j)=dsin(zdxxv)*dsin(y0dr)/dcos(zdyvr)*upv(i,j)+(dcos(zdxxv) &
            *dsin(zdyyv*fac)*dsin(y0dr)/dcos(zdyvr)+dcos(y0dr)*      &
            dcos(zdyyv*fac)/dcos(zdyvr))*vpv(i,j)
  upvp(i,j)=(dcos(y0dr)*dcos(zdyvr)-dsin(y0dr)*dsin(zdyvr)*dcos(zdxvr))/(dsin(y0dr)* &
            dsin(zdxvr))*vpvp(i,j)-dcos(zdyyv*fac)/(dsin(y0dr)*dsin(zdxvr))*vpv(i,j)
  endif
  enddo
 enddo

 return
end subroutine rot_wind
!=======================================================================
subroutine anti_rot_wind(x0,y0,alon,alat,alon_rot,alat_rot,u,v,uf,vf,nlon,nlat)

! Calculation of horizontal wind components in non rotated geography (UF, VF)
! for input wind components fields given in rotated geography (U, V)
! Some computations require double precision,
! but inaccuracies remain at poles and along the meridian of longitude x0

implicit none

real :: x0, y0
integer :: nlon, nlat, jlat, jlon
real, dimension(nlon,nlat) :: alat, alon, alon_rot, alat_rot, u, v, uf, vf
real*8 :: zfac, zx0, zy0, zx, zy, zxx, zyy, zu, zv, zzv, zmod, zphi, zpi

if (abs(x0) > 0.01.or.abs(y0) > 0.01) then

! Case of rotated grid

  zpi = dabs(dacos(-1.d0))
  zfac = zpi/180.d0
  zx0  = dble(x0)*zfac
  zy0  = dble(y0)*zfac
  do jlat=1,nlat
  do jlon=1,nlon
    zx=dble(alon(jlon,jlat))*zfac
    if (zx-zx0.le.-zpi) zx = zx +2.*zpi
    if (zx-zx0.ge. zpi) zx = zx -2.*zpi
    zy=dble(alat(jlon,jlat))*zfac
    zxx=dble(alon_rot(jlon,jlat))*zfac
    zyy=dble(alat_rot(jlon,jlat))*zfac
    zu=dble(u(jlon,jlat))
    zv=dble(v(jlon,jlat))

    if (dabs(dcos(zy))>1.d-2) then
      if (dabs(dsin(zx-zx0))>1.d-1) then
      zzv=-dsin(zy0)*dsin(zxx)*zu/dcos(zy)+zv*(dcos(zy0)*dcos(zyy)-dsin(zy0)*dsin(zyy)*dcos(zxx))/dcos(zy)
      zu=dcos(zyy)/dsin(zx-zx0)/dsin(zy0)*(zv-zzv/dcos(zyy)*(dcos(zx-zx0)*dsin(zy0)*dsin(zy)+dcos(zy0)*dcos(zy)))
      zv=zzv
      uf(jlon,jlat)=sngl(zu)
      vf(jlon,jlat)=sngl(zv)
      else
       if (zx-zx0.lt.-zpi/2.d0.or.zx-zx0.gt.zpi/2.d0) then
       uf(jlon,jlat)=-sngl(zu)
       vf(jlon,jlat)=-sngl(zv)
       else
       uf(jlon,jlat)=sngl(zu)
       vf(jlon,jlat)=sngl(zv)
       endif
      endif

    else
    zmod = dsqrt(zu**2+zv**2)
    zphi = zpi/2.d0+datan2(zv,zu)
    if(zphi.gt.zpi) zphi = zphi-2.d0*zpi
    if(zphi.lt.zpi) zphi = zphi+2.d0*zpi
    if (zy0.lt.0.d0) zx = -zx             ! south pole
    uf(jlon,jlat) = -zmod*dsin(zx-zphi)
    vf(jlon,jlat) = -zmod*dcos(zx-zphi)
    endif
  enddo
  enddo

else

! Case of non rotated grid

  uf(:,:)=u(:,:)
  vf(:,:)=v(:,:)

endif

return
end subroutine anti_rot_wind
!=======================================================================
subroutine antiwind(u,v,nx,ny,xdal,ydal,xgrid,ygrid,x0d,y0d)

!  Anti-rotates the wind vector (u,v) defined over a regular lat-lon grid.
!  The input is overwritten.
!  XDAL, YDAL are the coordinates of the lat-lon regular domain
!  XGRID, YGRID are the rotated coordinates of the regular lat-lon grid
!  (defined by subroutine antigrid) in degrees.

 real,dimension(nx,ny) :: u,v
 real*8,dimension(nx,ny) :: xgrid,ygrid
 real*8,dimension(nx) :: xdal(nx)
 real*8,dimension(ny) :: ydal(ny)
 real*8 :: x0d,y0d,x0,y0
 real*8 :: pi,fac,xp,yp,x,y

 pi=abs(acos(-1.))
 fac = pi/180.
 x0=x0d*fac
 y0=y0d*fac

 do i=1,nx
 do j=1,ny
 xp=xgrid(i,j)*fac
 yp=ygrid(i,j)*fac
 x=xdal(i)*fac
 y=ydal(j)*fac
   if (abs(x-x0).gt.1.e-3) then  ! With 10-4 or less, noise may appear at the central meridian
   vv = -sin(y0)*sin(xp)*u(i,j)/cos(y) + v(i,j)*                &
            (cos(y0)*cos(yp)-sin(y0)*sin(yp)*cos(xp))/cos(y)
   u(i,j) = cos(yp)/sin(x-x0)/sin(y0) * ( v(i,j) - vv/cos(yp)*  &
            (cos(x-x0)*sin(y0)*sin(y)+cos(y0)*cos(y)) )
   v(i,j) = vv
   endif
 enddo
 enddo

 return
end subroutine antiwind
!=======================================================================
subroutine near(x,n,xe,ne,iv)

! Finds point XE closest (less) to X, for any X
! Used for application to the spline interpolation

! Array IV(J) contains indexes of points with coordinates XE(J) nearest
! to points with coordinate X(J), in INCREASING order

  real x(n),xe(ne)
  integer iv(n)

  do j=1,n
    do je=1,ne
    if (xe(je).gt.x(j)) exit
    enddo
  iv(j)=je-1
  enddo

  return
end subroutine near
!=======================================================================
subroutine interp_spline_1d(f,x,n,fe,xe,ne,iv,alf,ex1,ex2)

! Interpolation (for 1D array) with cubic spline

! F:  output array, with dimension N, contains values of the interpolated values
! X:  input array, with dimension N, contains point coordinates of F
! FE: input array, with dimension NE, contains known values at prescribed point
! XE: input array, with dimension NE, contains coordinates of prescribed points
! (they must be in increasing order!)
! IV: input array, with dimension N, contains indexes of points with coordinates
! of X points nearest (less) to XE points (IV is defined by SUBROUTINE NEAR)
! ALF: tension coefficient 0-1: if ALF=1, then the interpolation is linear,
! if ALF=0, then the interpolation is cubic spline without tension
! EX1 and EX2 are parameters for possible extrapolation outside extremes of XE,
! at lower (EX1) and higher (EX2) extreme, respectively.
! If EX1 or EX2 =0, the extrapolation is constant; if =1, a linear extrapolation is applied.

implicit none

integer :: n, ne
real :: alf, ex1, ex2
real, dimension(n) :: f, x
real, dimension(ne) :: fe, xe
integer, dimension(n) :: iv

real, dimension(ne) :: fe1, fe4, xe1, xe4
integer :: k, k1, k2, i
real :: fm, fmm, fp, fpp, xm, xmm, xp, xpp, delx, delxp, delxm, delx1, delx2, &
        delxs, delx1s, delx2s, vm, vmm, vp, vpp, spl, flin

! First loop: extrapolation backward

 do k=1,n
   if (iv(k)>0) exit
   f(k)=f1(k)
 enddo
 k1=k

! Second loop: extrapolation forward

 do k=n,k1,-1
  if (iv(k)<ne) exit
  f(k)=f2(k)
 enddo
 k2=k

! Third loop: interpolation (in two separate phases)

! First phase: definit. of extreme values F1,F4 and X1,X4 for all possible values of IV(K)

 fe1(1)=2.*fe(1)-fe(2)
 xe1(1)=2.*xe(1)-xe(2)
 fe4(1)=fe(3)
 xe4(1)=xe(3)

 do i=2,ne-2
   fe1(i)=fe(i-1)
   xe1(i)=xe(i-1)
   fe4(i)=fe(i+2)
   xe4(i)=xe(i+2)
 enddo

 fe4(ne-1)=2.*fe(ne)-fe(ne-1)
 xe4(ne-1)=2.*xe(ne)-xe(ne-1)
 fe1(ne-1)=fe(ne-2)
 xe1(ne-1)=xe(ne-2)

! Second phase: interpolation into X points where interpolation is possible

 do k=k1,k2
   i=iv(k)
   fmm=fe1(i)
   fm=fe(i)
   fp=fe(i+1)
   fpp=fe4(i)

   xmm=xe1(i)
   xm=xe(i)
   xp=xe(i+1)
   xpp=xe4(i)

   delx=xp-xm
   delxp=xpp-xp
   delxm=xm-xmm
   delx1=x(k)-xm
   delx2=xp-x(k)
   delxs=delx**2
   delx1s=delx1**2
   delx2s=delx2**2

   vm=fm*(delx2/delx+delx1*delx2s/(delxs*delxm)-delx1s*delx2/((delx+delxp)*delxs))
   vp=fp*(delx1/delx+delx1s*delx2/(delxs*delxp)-delx2s*delx1/((delx+delxm)*delxs))
   vmm=fmm*delx1*delx2s/((delx+delxm)*delx*delxm)
   vpp=fpp*delx1s*delx2/((delx+delxp)*delx*delxp)
   spl=vm+vp-vmm-vpp
   flin=(fm*delx2+fp*delx1)/delx
   f(k)=alf*flin+(1.-alf)*spl

 enddo

return

contains

! Definition of functions

 real function f1(k)
 integer :: k
   f1=fe(1)+ex1*(fe(1)-fe(2))/(xe(1)-xe(2))*(x(k)-xe(1))
 end function f1
 real function f2(k)
 integer :: k
   f2=fe(ne)+ex2*(fe(ne)-fe(ne-1))/(xe(ne)-xe(ne-1))*(x(k)-xe(ne))
 end function f2

end subroutine interp_spline_1d
!=======================================================================
subroutine interp_spline_2d(a,nx,ny,xi,yi,xo,yo,ntot,f,alfa)

!  Horizontal interpolation between 2-D grids

!  Version to interpolate from a regular grid in input
!  with prescribed values defined in A(NX,NY) and with coordinates in XI and YI
! (but note that only XI(1), XI(2), YI(1) and YI(2) are used to recompute the
!  coordinates of the regualar grid) to points of coordinates XO,YO.
!  Output values defined in vector F(NTOT).
!  The second derivatives tend to zero near the boundaries
!  (equivalent to using linear extrap. out of the input grid).
!  Bicubic splines with tension are applied, combining linearly
!  pure splines and bilinear interpolation, using the input param. ALFA:
!  if ALFA = 1, then bilinear interp. is applied;
!  if ALFA=0, spline with no tension ("thin plate spline") is applied.
!  Points of the output grid must not be external to the area of the input grid.
!  It is assumed that fictitiuous (missing) values of the input grid A may exist,
!  being identified with values typically >= 1.E38 (for grid with only frames defined).
!  If such values are encountered, interpolation is not done and a conventional
!  value 1.E37 is attributed to F.

 real a(nx,ny),xi(nx),yi(ny),xo(ntot),yo(ntot),f(ntot)

 if(alfa.lt.0..or.alfa.gt.1.) then
   write(*,*) "Param. alfa out of range 0-1 in subroutine interp_spline_2d", alfa
   stop
 endif

 dx=xi(2)-xi(1)
 dy=yi(2)-yi(1)
 dx3=1./dx**3
 dy3=1./dy**3
 dxr=1./dx
 dyr=1./dy
 dxyr=dxr*dyr

! Case of linear interpolation only - it may be the case that spline
! interpolation is not possible because not all necessary points are defined
! (for example in a too narrow frame)

 if(alfa.gt.0.9999) then
  do n=1,ntot
   i=int((xo(n)-xi(1))*dxr+1.0)
   j=int((yo(n)-yi(1))*dyr+1.0)

! Check that output grid coordinates are internal

   if(i.lt.1.or.i.gt.(nx).or.j.lt.1.or.j.gt.(ny)) then
   write(*,*) "Point to be interpolated in subr. interp_spline_2d external to input grid:"
   print*, "  n,  i, j,     xo(n),    xi(1),    yo(n),   yi(1)"
   write(*,'(3i4,4f10.4)') n, i, j, xo(n), xi(1), yo(n), yi(1)
   print*, "Stop."
   stop
   endif

  f(n) = 0.
  enddo
  goto 200
 endif

! Interp. with splines

 do n=1,ntot

 i=int((xo(n)-xi(1))*dxr+1.0)
 j=int((yo(n)-yi(1))*dyr+1.0)

! Check that output grid coordinates are internal

 if(i.lt.1.or.i.gt.(nx).or.j.lt.1.or.j.gt.(ny)) then
   write(*,*) "Point to be interpolated in subr. interp_spline_2d external to input grid:"
   print*, "  n,  i, j,     xo(n),    xi(1),    yo(n),   yi(1)"
   write(*,'(3i4,4f10.4)') n, i, j, xo(n), xi(1), yo(n), yi(1)
   print*, "Stop."
   stop
 endif

 x1=xo(n)-xi(i)
 x2=dx-x1
 y1=yo(n)-yi(j)
 y2=dy-y1

 gx1=.5*x1*x2*x2*dx3
 gy1=.5*y1*y2*y2*dy3
 gx2=.5*x1*x1*x2*dx3
 gy2=.5*y1*y1*y2*dy3

! Vanishing of coefficients GX and GY in case extrapolation would be required
! (equivalent to implicit linear extrap., i.e diminishing the second derivative
! towards the border)

 if(i.eq.1)    gx1=0.
 if(j.eq.1)    gy1=0.
 if(i.eq.nx-1) gx2=0.
 if(j.eq.ny-1) gy2=0.
 if(i.eq.nx)   gx2=0.
 if(j.eq.ny)   gy2=0.

 ax = x2*dxr+2.*gx1-gx2
 ay = y2*dyr+2.*gy1-gy2
 bx = x1*dxr+2.*gx2-gx1
 by = y1*dyr+2.*gy2-gy1
 cx = -gx1
 cy = -gy1
 ex = -gx2
 ey = -gy2

! Reset of indices exceeding limits (null coefficients were attributed)

 im1 = max(i-1,1)
 jm1 = max(j-1,1)
 ip1 = min(i+1,nx)
 jp1 = min(j+1,ny)
 ip2 = min(i+2,nx)
 jp2 = min(j+2,ny)

! Check of possible fictitious values (indicating missing data)

 ztest= abs(a(im1,jm1))+abs(a(i,  jm1))+abs(a(ip1,jm1))+  &
        abs(a(ip2,jm1))+abs(a(im1,j  ))+abs(a(i,  j  ))+  &
        abs(a(ip1,j  ))+                                  &
        abs(a(ip2,j  ))+abs(a(im1,jp1))+abs(a(i,  jp1))+  &
        abs(a(ip1,jp1))+abs(a(ip2,jp1))+abs(a(im1,jp2))+  &
        abs(a(i,  jp2))+abs(a(ip1,jp2))+abs(a(ip2,jp2))

 if(ztest.gt.1.e14) then
   f(n) = 1.e37
   goto 100
 endif

! End of check

 f(n) =                                                        &
  cy*(cx*a(im1,jm1)+ax*a(i,jm1)+bx*a(ip1,jm1)+ex*a(ip2,jm1))+  &
  ay*(cx*a(im1,j  )+ax*a(i,j  )+bx*a(ip1,j  )+ex*a(ip2,j  ))+  &
  by*(cx*a(im1,jp1)+ax*a(i,jp1)+bx*a(ip1,jp1)+ex*a(ip2,jp1))+  &
  ey*(cx*a(im1,jp2)+ax*a(i,jp2)+bx*a(ip1,jp2)+ex*a(ip2,jp2))

100 continue    ! (for check)
  enddo

! Bilinear interp. and weighted mean using ALFA between the two interpolation algorithms

200 continue

 do n=1,ntot

! Check of "flag" values possibly set in the spline interpolation part
! (they are not changed and the linear interpolation part is avoided)

 if(f(n).gt.1.e36) goto 300

 i=int((xo(n)-xi(1))*dxr+1.0)
 j=int((yo(n)-yi(1))*dyr+1.0)

 ip1 = min(i+1,nx)
 jp1 = min(j+1,ny)

 x=(xo(n)-xi(i))*dxr
 y=(yo(n)-yi(j))*dyr
 f1=a(i,j)+x*(a(ip1,j)-a(i,j))
 f2=a(i,jp1)+x*(a(ip1,jp1)-a(i,jp1))

! Check of possible fictitiuous values encountered
! (needed only in the case ALFA=1 for wich the spline interp. is not done)

 ztest= abs(a(i,j))+abs(a(ip1,j))+abs(a(i,jp1))+abs(a(ip1,jp1))
 if(ztest.gt.1.e14) then
   f(n) = 1.e37
   goto 300
 endif

! End of check
! FLIN is the value obtained from the linear interpolation

 flin=f1+y*(f2-f1)

 if(alfa.gt.0.9999) then
  f(n) = flin
 else
  f(n)= (1.-alfa)*f(n)+alfa*flin
 endif
300 continue           ! (for check)
 enddo

 return
end subroutine interp_spline_2d
!=======================================================================
subroutine interp(alfa, ex1, ex2, npi, xi, g, x, f, nval)

! (Similar to INTERP_SPLINE_1D, but with different input, and does not require
!  call of SUBR. NEAR - moreover, it checks that input coordinates are monotonic).

!  Interpolates with splines with tension in one dimension.
!  The spline is defined imposing that the second derivative is the average
!  of second derivatives computed at the two adjacent points.
!  At interval extremes the second derivative is assumed null.
!  This subroutine also extrapolates out of the interval where the input funtion G is defined

!  INPUT:  function G defined at coordinates XI (CAUTION: can be changed by this subroutine)
!          G(1:NPI) values at irregular but strictly growing coordinates XI(1:NPI)
!  OUTPUT: F(1:NVAL) interpolated values at arbitrary coordinates X(1:NVAL)

!  ALFA: spline tension parameter, comprised between 0 and 1
!  If ALFA=1, pure linear interpolation; if ALFA=0, pure spline

!  EX1: param. determining extrapolation for X < XI(1)
!  EX2: param. determining extrapolation for X > XI(NPI)
!  If EX1=0 OR EX2=0, constant value extrapolation is used at corresponding extreme
!  If EX1=1 OR EX2=1, linear extrapolation is used at corresponding extreme
!  Intermediate values of EX1 AND EX2 give intermediate extrapolation values

  real, dimension(npi )  :: xi, g
  real, dimension(nval)  :: x,  f

  if(alfa.lt..0.or.alfa.gt.1) then
  print*, 'CAUTION: in INTERP, ALFA out of interval 0-1'
  endif
  if(ex1.lt..0.or.ex1.gt.1) then
  print*, 'CAUTION: in INTERP, EX1 out of interval 0-1'
  endif
  if(ex2.lt..0.or.ex2.gt.1) then
  print*, 'CAUTION: in INTERP, EX2 out of interval 0-1'
  endif

! Fix for the case in which coordinates of the input function are not strictly increasing
! Note that this changes the input coordinates in the calling programme

  do k=2,npi
   if(xi(k).le.xi(k-1)) then
   print*, "CAUTION: in INTERP, coordinates of input function changed because not monotonic!"
   exit
   endif
  enddo

  zeps=(xi(npi)-xi(1))*1.e-6   ! Small deviation used to set apart interlaced coordinates
200 do k=2,npi
     if(xi(k).le.xi(k-1)) then
     ximed=0.5*(xi(k)+xi(k-1))
     xi(k-1)=ximed-zeps
     xi(k)=ximed+zeps
     gmed=0.5*(g(k)+g(k-1))
     g(k-1)=gmed
     g(k)=gmed
     endif
    enddo

 do k=2,npi
  if(xi(k).le.xi(k-1)) then
  goto 200
  endif
 enddo

 do 100 jval =1, nval

!  2 cases of extrapolation

 if(x(jval).lt.xi(1)) then
 f(jval) = g(1) + ex1*(g(1)-g(2))/(xi(1)-xi(2)) * (x(jval)-xi(1))
 go to 100
 elseif (x(jval).ge.xi(npi)) then
 f(jval) = g(npi) + ex2*(g(npi)-g(npi-1))/(xi(npi)-xi(npi-1)) * (x(jval)-xi(npi))
 go to 100
 endif

 ir = 0

!  IR is a reference index determining the interpolation interval
!  The interpolation expression is applied also if X=XI(J)

 do j = 1, npi
 if (x(jval).ge.xi(j)) ir = ir + 1
 enddo

 if (ir.eq.1) then
 fmm = 2*g(1) - g(2)
 xmm = 2*xi(1) - xi(2)
 fpp = g(ir+2)
 xpp = xi(ir+2)
 elseif (ir.eq.(npi-1)) then
 fpp = 2*g(npi) - g(npi-1)
 xpp = 2*xi(npi) - xi(npi-1)
 fmm = g(ir-1)
 xmm = xi(ir-1)
 else
 fmm = g(ir-1)
 xmm = xi(ir-1)
 fpp = g(ir+2)
 xpp = xi(ir+2)
 endif

 fm     = g(ir)
 xm     = xi(ir)
 fp     = g(ir+1)
 xp     = xi(ir+1)
 delx   = xp - xm
 delxp  = xpp - xp
 delxm  = xm - xmm
 delx1  = x(jval) - xm
 delx2  = xp - x(jval)
 delxs  = delx**2
 delx1s = delx1**2
 delx2s = delx2**2

!  Spline contribution to interpolation

 spl = fm*(delx2/delx + delx1*delx2s/(delxs*delxm) - delx1s*     &
       delx2/((delx+delxp)*delxs)) + fp*(delx1/delx +            &
       delx1s*delx2/(delxs*delxp) - delx1*delx2s/((delx+delxm)*  &
       delxs)) - fmm * delx1*delx2s/((delx+delxm)*delx*delxm) -  &
       fpp * delx1s*delx2/((delx+delxp)*delx*delxp)

!  Linear interpolation contribution

 clin = (fm*delx2 + fp*delx1)/delx

!  Final interpolation combined using ALFA

 f(jval) = alfa*clin + (1.-alfa)*spl

 100  continue

 return
end subroutine interp
!=======================================================================
subroutine smooth(vt,vout,ni,nj,w,nsmooth)

!  Smoothing of a 2-D matrix VT(NI,NJ), based on a 5 point filter (laplacian).
!  Output in VOUT(NI,NJ) and also in VT (ATTENTION: VT is redefined)

!  It is assumed that grid distances in x and y are similar
!  The degree of filtering is proportional to W (from 0 to 1)
!  (better use multiple NSMOOTH steps rather than W larger than about 0.5)

 dimension vt(ni,nj), vout(ni,nj)

 do ns=1,nsmooth

   do j=2,nj-1
   do i=2,ni-1
   vout(i,j)=(1.-w)*vt(i,j)+.25*w*(vt(i+1,j)+vt(i-1,j)+vt(i,j+1)+vt(i,j-1))
   enddo
   enddo

!  Points at the borders (1-D filter)

   do i=2, ni-1
   vout(i,1)=(1.-.7*w)*vt(i,1)+.35*w*(vt(i+1,1)+vt(i-1,1))
   vout(i,nj)=(1.-.7*w)*vt(i,nj)+.35*w*(vt(i+1,nj)+vt(i-1,nj))
   enddo
   do j=2, nj-1
   vout(1,j)=(1.-.7*w)*vt(1,j)+.35*w*(vt(1,j+1)+vt(1,j-1))
   vout(ni,j)=(1.-.7*w)*vt(ni,j)+.35*w*(vt(ni,j+1)+vt(ni,j-1))
   enddo

!  Points at the corners (unchanged)

   vout(1,1)=vt(1,1)
   vout(1,nj)=vt(1,nj)
   vout(ni,nj)=vt(ni,nj)
   vout(ni,1)=vt(ni,1)

   do j=1,nj
   do i=1,ni
   vt(i,j) = vout(i,j)
   enddo
   enddo

  enddo

  return
end subroutine smooth
!=======================================================================
subroutine seatemp(tana,zlsm,sst,iana,jana,niter,nfilt,w)

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
 logical  ilwork(iana,jana),ilsst(iana,jana)

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
! ======================================================================
subroutine landtemp(tana,zlsm1,iana,jana,niter,nfilt,w)

! Defines soil quantities extending them from land to and beyond coasts

! Similar but inverted with respect to SEATEMP
! (note that input matrix TANA is overwritten, unlike in SEATEMP)

! 1) Definition of a logical vector initially TRUE only at "internal land" points
!    (defined as land points not adjacent to any sea point)
! 2) Expansion is repeated NITER times:
!    for any point with logical flag FALSE (sea or coastal)
!    temperature is redefined as average with only adjacent points having flag TRUE,
!    if they exist, and in this case the flag is set to TRUE.
! 3) A 4-point filter is optionally applied NFILT times with weight W
!   (W: 0.5-1; if W=1, no filtering)
! ----------------------------------------------------------------------

 real  tana(iana,jana),zlsm1(iana,jana),zlsm(iana,jana)
 real  aland(iana,jana)
 real  twork(iana,jana)
 logical ilwork(iana,jana),ilsst(iana,jana)

! Identification of internal land points
! LSM is inverted

  do j=1,jana
  do i=1,iana
  zlsm(i,j)=1.-zlsm1(i,j)
  enddo
  enddo

  do j=1,jana
    do i=1,iana
      twork(i,j) = tana(i,j)

      if (zlsm(i,j).lt.0.55) then                     ! sea
        ilwork(i,j) = .false.
      else                                            ! land
        ilwork(i,j) = .true.
      endif
   enddo
  enddo

! Correction of the flag on lateral boundaries (where it assumed that
! a sea point is a coastal point)

 do i=1,iana
   ilwork(i,1) = .false.
   ilwork(i,jana) = .false.
 enddo

 do j=1,jana
   ilwork(1,j) = .false.
   ilwork(iana,j) = .false.
 enddo

! Expansion of the land field towards the sea

 do 300 n=1,niter
   do j=1,jana
     do i=1,iana
       if (ilwork(i,j)) then
         aland(i,j) = twork(i,j)
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
           aland(i,j) = tsi/nsi
           ilsst(i,j) = .true.
         else
           aland(i,j) = twork(i,j)
           ilsst(i,j) = .false.
         endif
       endif
    enddo
  enddo

  do j=1,jana
    do i=1,iana
     twork(i,j) = aland(i,j)
     ilwork(i,j) = ilsst(i,j)
    enddo
  enddo

300 continue

! Filtering of field ALAND

  do n=1,nfilt

    do j=2,jana-1
      do i=2,iana-1
       t4med = (twork(i-1,j)+twork(i+1,j)+twork(i,j-1)+twork(i,j+1))/4.
       aland(i,j) = w*twork(i,j) + (1.-w)*t4med
      enddo
    enddo

    do j=2,jana-1
      do i=2,iana-1
       twork(i,j) = aland(i,j)
      enddo
    enddo
  enddo

! Input matrix overwritten

   do j=1,jana
     do i=1,iana
      tana(i,j) = aland(i,j)
     enddo
   enddo

 return
end subroutine landtemp
!=======================================================================
subroutine plotout(a,b,n,m,title,nhf)

! Creates a file containing 2D fields (matrix a) for plotting.
! Matrix must contain the fmask to be plotted superimposed on each field.
! Title of the graph and dimension of the matrix are written first.
! The title can be composed by two words (e.g. "word1 word2"): define the string as "word1_word2"

   real a(n,m)
   real b(n,m)
   character*30 title

   write(nhf,*) title
   write (nhf,*) n, m
   do j=1,m
   write (nhf,37) (b(i,j),i=1,n)
   enddo
   do j=1,m
   write (nhf,38) (a(i,j),i=1,n)
   enddo

37 format (20f6.3)
38 format (10e12.5)

   return
end subroutine plotout
!=======================================================================
subroutine outgraph(icentre_code,imodel_code,nx,ny,x0,y0,alon0,alat0,dx,dy,                &
           param_discipl,param_categ,param_ind,lev_type,lev1,lev2,idate0,iperiod,a,zf1,zf2)

  implicit none

  integer :: nx, ny, icentre_code, imodel_code, param_discipl, param_categ, param_ind, &
             lev_type, lev1, lev2, idate0(5), iperiod(3)
  real, dimension(nx,ny) :: a
  real :: x0, y0, alon0, alat0, dx, dy, zf1, zf2

! CNR-ISAC-BO: ICENTRE_CODE=80, ISUBCENTRE_CODE=102
! BOLAM: IMODEL_CODE=1, MOLOCH: IMODEL_CODE=2, GLOBO: IMODEL_CODE=3

!              PARAM_DISCIPL PARAM_CATEG PARAM_IND
!    T              0            0          0
!    Q              0            1          0

!                          LEV_TYPE
! Ground or water surface       1
! Hybrid level                105
! Isobaric surface (Pa)       100

  open (11,file='data_output.bin',status='unknown',form='unformatted',position='append')

  write (11) icentre_code,imodel_code
  write (11) nx,ny
  write (11) x0,y0,alon0,alat0,dx,dy
  write (11) param_discipl,param_categ,param_ind
  write (11) lev_type,lev1,lev2
  write (11) idate0(1:5),iperiod(1:3)
  write (11) (a(1:nx,1:ny)*zf1+zf2)

  close (11)

return
end subroutine outgraph
!=======================================================================
 subroutine redistr_snow(snow, topog, nx, ny, nit)

! Redistributes a snow height field (snow) as a function of topography height (topog),
! maintaining the average snow and increasing the snow over peaks at expenses of valleys
! with an iterative algorithm (nit: no. of steps) - boundary values remain unchanged

 real snow(nx,ny), topog(nx,ny), topincr(nx,ny), workt(nx,ny)

 topincr = 0.
 if(maxval(snow).lt.1.e-4) then
 print*, "No snow cover found (subr. redistr_snow)."
 return
 endif

 avsnowi = 0.
 do j = 2, ny-1
 do i = 2, nx-1
 avsnowi = avsnowi + snow(i,j)
 enddo
 enddo
 avsnowi = avsnowi/float((nx-2)*(ny-2))
! print*, "Initial average snow", avsnowi

 workt = topog

 do it = 1, nit

   do j = 2, ny-1
   do i = 2, nx-1

! 9-point weighted mean (area average mean) of snow and topography

   topav = 0.25*workt(i,j) +                                                         &
           0.125*(workt(i-1,j) + workt(i+1,j) + workt(i,j-1) + workt(i,j+1)) +       &
           0.0625*(workt(i-1,j-1) + workt(i-1,j+1) + workt(i+1,j-1) + workt(i+1,j+1))
           topincr(i,j) = workt(i,j) - topav
   enddo
   enddo

   workt = workt - topincr

   dincr = 0.00055
   do j = 2, ny-1
   do i = 2, nx-1
   snowav = 0.25*snow(i,j) +                                                 &
            0.125*(snow(i-1,j)+snow(i+1,j)+snow(i,j-1)+snow(i,j+1)) +        &
            0.0625*(snow(i-1,j-1)+snow(i-1,j+1)+snow(i+1,j-1)+snow(i+1,j+1))
   snow(i,j) = max(snow(i,j) + topincr(i,j)*dincr*snowav, 0.)
   enddo
   enddo

 enddo

 avsnowf = 0.
 do j = 2, ny-1
 do i = 2, nx-1
 avsnowf = avsnowf + snow(i,j)
 enddo
 enddo
 avsnowf = avsnowf/float((nx-2)*(ny-2))

! print*, "Final average snow", avsnowf
! print*, "avsnowi/avsnowf", avsnowi/avsnowf

 do j = 2, ny-1
 do i = 2, nx-1
 snow(i,j) = snow(i,j)*avsnowi/avsnowf  ! correction to preserve the snow average
 enddo
 enddo
 print*, "Snow redistributed over orography (subr. redistr_snow)."

 return
 end
!=======================================================================
      subroutine bextr(amatr, ni, nj, ex, klper)

!  amatr is input and output - defines boundary values in amatr as
!  extrapolated from interior, depending on parameter ex
!  0<=ex<=1; if ex=1: linear extrap.; if ex=0: constant value as interior
!  ni and nj are matrix dimensions

      dimension amatr(ni,nj)
      logical klper

      if(ex.lt.0..or.ex.gt.1.) then
      print*, 'Caution: param. ex in subr. bextr is out of interval 0-1!'
      endif
      do i = 2, ni-1
      amatr(i,1)  = (1.+ex)*amatr(i,2)    - ex*amatr(i,3)
      amatr(i,nj) = (1.+ex)*amatr(i,nj-1) - ex*amatr(i,nj-2)
      enddo

      if(klper) then
      do j = 1, nj
      amatr(1,j)  = amatr(ni-1,j)
      amatr(ni,j) = amatr(2,j)
      enddo
      else
      do j = 1, nj
      amatr(1,j)  = (1.+ex)*amatr(2,j)   - ex*amatr(3,j)
      amatr(ni,j) = (1.+ex)*amatr(ni-1,j) - ex*amatr(ni-2,j)
      enddo
      endif

      return
      end
! ======================================================================
subroutine interp_gauss (data_inp, data_out, alon_inp, alat_inp, alon_out, alat_out,&
np_inp, np_out, radius, par1, par2, valmiss)

implicit none

integer :: np_inp, np_out, i, npoint, np, ip
real, dimension(np_inp) :: data_inp, alon_inp, alat_inp, data_inp_val, alon_inp_val, alat_inp_val, data_inp_val_2
real, dimension(np_out) :: data_out, alon_out, alat_out
integer, parameter :: npoint0=10000
integer, dimension(npoint0) :: index_point
real, dimension(npoint0) :: dist_point, func_weight_point, weight_point
real :: radius, valmiss, pi, distfac, fac, par1, par2, par, zdist, zzz

pi=abs(acos(-1.))
distfac=6371.*pi/180.
fac=pi/180.
par=par1*(par2/radius)

 npoint=0
 do i=1,np_inp
   if (int(data_inp(i))/=int(valmiss)) then
     npoint=npoint+1
     data_inp_val(npoint)=data_inp(i)
     alon_inp_val(npoint)=alon_inp(i)
     alat_inp_val(npoint)=alat_inp(i)
   endif
 enddo

 do i=1,np_out
   np=0
   do ip=1,npoint
     zdist=sqrt(((alon_out(i)-alon_inp_val(ip))*cos(abs(alat_out(i)*fac))*distfac)**2+&
((alat_out(i)-alat_inp_val(ip))*distfac)**2)
     if (zdist<=radius) then
       np=np+1
       dist_point(np)=zdist
       func_weight_point(np)=exp(-(zdist*par)**2)
       data_inp_val_2(np)=data_inp_val(ip)
     endif
   enddo
   if (np>0) then
     zzz=0.
     do ip=1,np
       zzz=zzz+(func_weight_point(ip))
     enddo
     data_out(i)=0.
     do ip=1,np
       weight_point(ip)=func_weight_point(ip)/zzz
       data_out(i)=data_out(i)+data_inp_val_2(ip)*weight_point(ip)
     enddo
   else
     data_out(i)=valmiss
   endif
 enddo

return
end subroutine interp_gauss
! ======================================================================
subroutine interp_gauss_input_grid (data_inp, data_out, alon_inp, alat_inp, alon_out, alat_out, &
 ni_inp, nj_inp, np_out, nlev, npar, mask, semi_radius, gauss_par, valmiss)

implicit none

! Interpolation from grid input data (data_inp) with dimensions (ni_inp,nj_inp,nlev,npar)
! to output points (data_out) with dimensions (np_out,nlev,npar),
! using Gaussian approach to calculate weights of input grid points

! npar               : number of parameters (fields)
! nlev               : number of vertical levels
! ni_inp, nj_inp     : input grid point number along axes X and Y
! np_inp             : number of output points
! alon_inp, alat_inp : longitude and latitude of input data grid points
! alon_out, alat_out : longitude and latitude of output data points
! semi_radius        : semi-radius of Gaussing function (km)
! par                : a parameter of Gaussing function
! valmiss            : value used in the case of procedure failure in an output point

! Input :

integer :: ni_inp, nj_inp, np_out, nlev, npar
real, dimension(ni_inp,nj_inp,nlev,npar) :: data_inp
real, dimension(ni_inp,nj_inp) :: alon_inp, alat_inp
integer, dimension(ni_inp,nj_inp) :: mask
real, dimension(np_out) :: alon_out, alat_out
real :: semi_radius, gauss_par, valmiss

! Output :

real, dimension(np_out,nlev,npar) :: data_out

! Work :

integer, parameter :: npoint0=10000
real, dimension(npoint0,np_out) :: dist_point, func_weight_point, weight_point
integer, dimension(npoint0) :: npoint
integer, dimension(npoint0,np_out) :: i_inp, j_inp
real, dimension(np_out) :: func_weight_sum
integer :: i, ip, jp, ipoint, ipar, ilev
real :: pi, fac, distfac, radius, par, zdist

 pi=abs(acos(-1.))
 distfac=6371.*pi/180.
 fac=pi/180.
 radius=semi_radius*3.
 par=gauss_par/semi_radius

! Definition of gaussian weight of input data for output points

 do i=1,np_out ! Output data points

   npoint(i)=0
   do jp=1,nj_inp
   do ip=1,ni_inp
     zdist=sqrt(((alon_out(i)-alon_inp(ip,jp))*cos(abs(alat_out(i)*fac))*distfac)**2+ &
           ((alat_out(i)-alat_inp(ip,jp))*distfac)**2)
     if (zdist <= radius.and.mask(ip,jp) /= 0) then
       npoint(i)=npoint(i)+1
       ipoint=npoint(i)
       dist_point(ipoint,i)=zdist
       func_weight_point(ipoint,i)=exp(-(zdist*par)**2)
       i_inp(ipoint,i)=ip
       j_inp(ipoint,i)=jp
     endif
   enddo
   enddo

   func_weight_sum(i)=0.
   if (npoint(i) > 0) then
     do ipoint=1,npoint(i)
       func_weight_sum(i)=func_weight_sum(i)+func_weight_point(ipoint,i)
     enddo
     do ipoint=1,npoint(i)
       weight_point(ipoint,i)=func_weight_point(ipoint,i)/func_weight_sum(i)
     enddo
   endif

 enddo ! Output data points

! Definition of output data at all levels and for all parameters

 do ipar=1,npar
   do ilev=1,nlev
     do i=1,np_out

       if (npoint(i) > 0) then
         data_out(i,ilev,ipar)=0.
         do ipoint=1,npoint(i)
           ip=i_inp(ipoint,i)
           jp=j_inp(ipoint,i)
           data_out(i,ilev,ipar)=data_out(i,ilev,ipar)+data_inp(ip,jp,ilev,ipar)*weight_point(ipoint,i)
         enddo
       else
         data_out(i,ilev,ipar)=valmiss
       endif

     enddo
   enddo
 enddo

return
end subroutine interp_gauss_input_grid
