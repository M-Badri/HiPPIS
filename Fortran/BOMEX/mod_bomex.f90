module mod_bomex

use physcons, only : rair => con_rd, &
                     rh20 => con_rv, &
                     cp   => con_cp, &
                     cv   => con_cv, &
                     grav => con_g,  &
                     eps0 => con_eps

implicit none

integer, parameter :: r8 = SELECTED_REAL_KIND(8)

real(kind=r8), parameter :: p0 = 100000.0_r8
contains

subroutine bomex_ls_forcing(u,v,p,rho,th,qv,qc,qr, &
                         ug,vg,w,dthdt,dqvdt,z,dt,nz, bomex_type)

  implicit none

  integer,                      intent(in) :: nz  ! number of model levels
  real(kind=r8),                intent(in) :: dt
  real(kind=r8), dimension(nz), intent(in) :: w,     &  ! prescribed subsidence [m/s]
                                              ug,    &  ! prescribed zonal geostropic wind [K/s]
                                              vg,    &  ! prescribed meridional geostropic wind [K/s]
                                              dthdt, &  ! prescribed heating [K/s]
                                              dqvdt, &  ! prescribed drying  [kg/kg/s]
                                              z         ! prescribed drying  [kg/kg/s]
  character(len=16),            intent(in) :: bomex_type
  real(kind=r8), dimension(nz), intent(inout) :: u,    &  ! zonal velocity [m/s]
                                                 v,    &  ! meridional velocity [m/s]
                                                 p,    &  ! pressure [Pa]
                                                 rho,  &  ! density [kg/m/m/m]
                                                 th,   &  ! potential temperautre [K]
                                                 qv,   &  ! vapor mixing ratio [kg/kg]
                                                 qc,   &  ! cloud water mixing ratio [kg/kg]
                                                 qr       ! rain water mixing ratio [kg/kg]

  integer :: k, kp, km
  real(kind=r8) :: thv, dz

  integer  :: i
  real(kind=r8), dimension(nz) :: du, dv, dth, dqv, dqr, dqc
  real(kind=r8), dimension(nz) :: g1, g2, g3
  real(kind=r8), dimension(nz) :: fluxm, fluxp, u1, u_old, u_new, tmp2
  real(kind=r8), dimension(nz) :: u2, v2, th2, qv2, qr2, qc2
  real(kind=r8), dimension(nz) :: u3, v3, th3, qv3, qr3, qc3
  real(kind=r8), dimension(nz-1) :: fu, fv, fth, fqv, fqr, fqc
  real(kind=r8), dimension(nz) :: dudz, dvdz, dthdz, dqvdz, dqcdz, dqrdz
  real(kind=r8) :: tmp

  real(kind=r8), parameter :: f = 0.376e-4
  real(kind=r8) :: kappa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! bomex using linear and forward euler
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(bomex_type .eq. "linear")then
  !Upstream biased, first order differencing
  do k=1,nz
    kp=min(nz,k+1)
    km=min(nz-1,k)
    dz = z(kp)-z(km)
    du(k) = w(k)*(u(km) -u(kp) )/dz
    dv(k) = w(k)*(v(km) -v(kp) )/dz
    dth(k)= w(k)*(th(km)-th(kp))/dz
    dqv(k)= w(k)*(qv(km)-qv(kp))/dz
    dqc(k)= w(k)*(qc(km)-qc(kp))/dz
    dqr(k)= w(k)*(qr(km)-qr(kp))/dz
  end do
  do k=1,nz
    u(k)  = u(k) + dt*(f*(v(k)-vg(k)) + du(k))
    v(k)  = v(k) + dt*(-f*(u(k)-ug(k))+ dv(k))
    th(k) = th(k)+ dt*(dthdt(k)+dth(k))
    qv(k) = qv(k)+ dt*(dqvdt(k)+dqv(k))
    qc(k) = qc(k)+ dt*dqc(k)
    qr(k) = qr(k)+ dt*dqr(k)
    thv   = th(k)*(qv(k)+eps0)/(eps0*(1.0_r8+qv(k)))
    p(k)  = p0*(thv*rair*rho(k)/p0)**(cp/cv)
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! BOMEX advection using a weighted essential non-oscillatory (WENO) finite
!! differnce method with a third order Runge-kutta method. 
!! This implementation is based on the Chi-Wang Shu ``High-order 
!! Finite Difference and Finite Volume WENO Schemes and Discontinuous 
!! Galerkin Methods for CFD".
!! https://doi.org/10.1080/1061856031000104851
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  elseif(bomex_type .eq. "weno") then
  !!do k=1, nz
  !!  write(*,*), 'th(',k,')=', th(k)
  !!enddo
  do i=1,6
    !!if( i .eq. 1 .and. i .eq. 2)then
    if(i.eq.1) then
      u_new = u
      tmp2 = f*(v -vg)
    elseif(i.eq.2) then
      u_new = v
      tmp2 = -f*(u-ug)
    elseif(i.eq.3)then
      u_new = qv
      tmp2 = dqvdt
    elseif(i.eq.4) then
      u_new = qc
      tmp2 = 0.0
    elseif(i.eq.5) then
      u_new = qr
      tmp2 = 0.0
    elseif(i.eq.6) then
      u_new = th
      tmp2 = dthdt
    endif

    call computeflux(u_new, fluxm, fluxp, nz)
    u1=u_new
    u_old=u_new
    do k=4,nz-3
      if( w(k) >= 0 )then
        tmp = w(k)*(fluxm(k)-fluxm(k-1))/(z(k)-z(k-1))
      else
        tmp = -w(k)*(fluxp(k+1)-fluxp(k))/(z(k+1)-z(k))
      endif
      u1(k) = u_old(k) + dt*(tmp2(k) + tmp)
    enddo
    !!u_new = u1
    call computeflux(u1, fluxm, fluxp, nz)
    u2 = u1
    do k=4,nz-3
      if(w(k) >= 0)then
        tmp = w(k)*(fluxm(k)-fluxm(k-1))/(z(k)-z(k-1))
      else
        tmp = -w(k)*(fluxp(k+1)-fluxp(k))/(z(k+1)-z(k))
      endif
      u2(k) = 1.0/4.0*(3.0*u_old(k) + u1(k) + dt*(tmp2(k) + tmp))
    enddo

    call computeflux(u2, fluxm, fluxp, nz)
    u_new = u2
    do k=4,nz-3
      if(w(k) >= 0)then
        tmp = w(k)*(fluxm(k)-fluxm(k-1))/(z(k)-z(k-1))
      else
        tmp = -w(k)*(fluxp(k+1)-fluxp(k))/(z(k+1)-z(k))
      endif
      u_new(k) = 1.0/3.0*(u_old(k) + 2.0*u2(k) + 2.0*dt*(tmp2(k) + tmp))
    enddo

    do k=1,3
      kp=min(nz,k+1)
      km=min(nz-1,k)
      dz = z(kp)-z(km)
      du(k) = w(k)*( u_old(km)-u_old(kp) )/dz
      u_new(k)  = u_old(k) + dt*( tmp2(k) + du(k) )
    enddo
    do k=nz-2,nz
      kp=min(nz,k+1)
      km=min(nz-1,k)
      dz = z(kp)-z(km)
      du(k) = w(k)*( u_old(km)-u_old(kp) )/dz
      u_new(k)  = u_old(k) + dt*(tmp2(K) + du(k))
    enddo

    if(i.eq.1) then
      u = u_new
    elseif(i.eq.2) then
      v = u_new
    elseif(i.eq.3)then
      qv = u_new 
    elseif(i.eq.4) then
      qc = u_new
    elseif(i.eq.5) then
      qr = u_new
    elseif(i.eq.6) then
      th = u_new 
      
    endif
    !!endif
  enddo


  !! compute p(k)
  do k=1,nz
    !!write(*,*)  'k=',k, 'th(k)=', th(k), 'qv(k)=', qv(k), 'eps0=', eps0
    thv   = th(k)*(qv(k)+eps0)/(eps0*(1.0_r8+qv(k)))
    !!write(*,*) 'thv =', thv, 'rair=', rair, 'p0=', p0, 'cp', cp, 'cv=', cv, 'rho(k)=', rho(k)
    p(k)  = p0*(thv*rair*rho(k)/p0)**(cp/cv)
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! bomex using third order finite difference method with a third order
!! Runge Kutta method based on the work by Hundsdorfer et al. ``A 
!! positive finite-difference advection scheme".
!! https://doi.org/10.1006/jcph.1995.1042
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  elseif(bomex_type .eq. "cubic")then
  kappa = 1.0/3.0
  !Upstream biased, first order differencing
  do k=1,nz
    kp=min(nz,k+1)
    km=min(nz-1,k)
    dz = z(kp)-z(km)
    dudz(k) = (u(km) -u(kp) )/dz
    dvdz(k) = (v(km) -v(kp) )/dz
    dthdz(k)= (th(km)-th(kp))/dz
    dqvdz(k)= (qv(km)-qv(kp))/dz
    dqcdz(k)= (qc(km)-qc(kp))/dz
    dqrdz(k)= (qr(km)-qr(kp))/dz
  end do


 !!-- Compute the different flux functions --!!
 call flux(u, nz, fu, kappa)
   
  !!-- Compute spatial derivative of different functions --!!
  do k= 3, nz-2
    du(k) = -(fu(k)-fu(k-1)) /((z(k+1) - z(k-1))*0.5)
  enddo

  !!-- Compute time integration for different function --!!
  g1(1)  = f*(v(1)-vg(1)) + w(1)*dudz(1)
  u2(1)  = u(1) + dt*g1(1)
  g1(2)  = f*(v(2)-vg(2)) + w(2)*dudz(2)
  u2(2)  = u(2) + dt*g1(2)
  do k= 3, nz-2
    g1(k)  = f*(v(k)-vg(k)) + w(k)*du(k)
    u2(k)  = u(k) + dt*g1(k)
  enddo
  g1(nz-1)  = f*(v(nz-1)-vg(nz-1)) + w(nz-1)*dudz(nz-1)
  u2(nz-1)  = u(nz-1) + dt*g1(nz-1)
  g1(nz)  = f*(v(nz)-vg(nz)) + w(nz)*dudz(nz)
  u2(nz)  = u(nz) + dt*g1(nz)
  call flux(u2, nz, fu, kappa)
  do k= 3, nz-2
    du(k) = -(fu(k)-fu(k-1)) /((z(k+1) - z(k-1))*0.5)
  enddo

  g2(1)  = f*(v(1)-vg(1)) + w(1)*dudz(1)
  u3(1)  = u(1) + dt*0.25*(g1(1)+ g2(1))
  g2(2)  = f*(v(2)-vg(2)) + w(2)*dudz(2)
  u3(2)  = u(2) + dt*0.25*(g1(2) + g2(2) )
  do k= 3, nz-2
    g2(k)  = f*(v(k)-vg(k)) + w(k)*du(k)
    u3(k)  = u(k) + dt*0.25*(g1(k)+g2(k))
  enddo
  g2(nz-1)  = f*(v(nz-1)-vg(nz-1)) + w(nz-1)*dudz(nz-1)
  u3(nz-1)  = u(nz-1) + dt*0.25*(g1(nz-1)+g2(nz-1))
  g2(nz)  = f*(v(nz)-vg(nz)) + w(nz)*dudz(nz)
  u3(nz)  = u(nz) + dt*0.25*(g1(nz)+g2(nz))

  call flux(u3, nz, fu, kappa)
  do k= 3, nz-2
    du(k) = -(fu(k)-fu(k-1)) /((z(k+1) - z(k-1))*0.5)
  enddo
  g3(1)  = f*(v(1)-vg(1)) + w(1)*dudz(1)
  u(1)  = u(1) + dt*(1.0/6.0*g1(1)+ 1.0/6.0*g2(1) + 2.0/3.0*g3(1))
  g3(2)  = f*(v(2)-vg(2)) + w(2)*dudz(2)
  u(2)  = u(2) + dt*(1.0/6.0*g1(2) + 1.0/6.0*g2(2) + 2.0/3.0*g3(2) )
  do k= 3, nz-2
    g3(k)  = f*(v(k)-vg(k)) + w(k)*du(k)
    u(k)  = u(k) + dt*(1.0/6.0*g1(k)+1.0/6.0*g2(k)+2.0/3.0*g3(k))
  enddo
  g3(nz-1)  = f*(v(nz-1)-vg(nz-1)) + w(nz-1)*dudz(nz-1)
  u(nz-1)  = u(nz-1) + dt*(1.0/6.0*g1(nz-1)+1.0/6.0*g2(nz-1)+2.0/3.0*g3(nz-1))
  g1(nz)  = f*(v(nz)-vg(nz)) + w(nz)*dudz(nz)
  u3(nz)  = u(nz) + dt*(1.0/6.0*g1(nz)+1.0/6.0*g2(nz) + 2.0/3.0*g3(nz))

  call flux(v, nz, fv, kappa)    !! flux compuation
  do k= 3, nz-2
    dv(k) = -(fv(k)-fv(k-1)) /((z(k+1) - z(k-1))*0.5)
  enddo
  g1(1)  = -f*(u(1)-ug(1))+ w(1)*dvdz(1)
  v2(1)  = v(1) + dt*g1(1)
  g1(2)  = -f*(u(2)-ug(2))+ w(2)*dvdz(2)
  v2(2)  = v(2) + dt*g1(2)
  do k= 3, nz-2
    g1(k)  = -f*(u(k)-ug(k))+ w(k)*dv(k)
    v2(k)  = v(k) + dt*g1(k)
  enddo
  g1(nz-1)  = -f*(u(nz-1)-ug(nz-1))+ w(nz-1)*dvdz(nz-1)
  v2(nz-1)  = v(nz-1) + dt*g1(nz-1)
  g1(nz)  = -f*(u(nz)-ug(nz))+ w(nz)*dvdz(nz)
  v2(nz)  = v(nz) + dt*g1(nz)

  call flux(v2, nz, fv, kappa)    !! flux compuation
  do k= 3, nz-2
    dv(k) = -(fv(k)-fv(k-1)) /((z(k+1) - z(k-1))*0.5)
  enddo
  g2(1)  = -f*(u(1)-ug(1))+ w(1)*dvdz(1)
  v3(1)  = v(1) + dt*0.25*(g1(1) + g2(1))
  g2(2)  = -f*(u(2)-ug(2))+ w(2)*dvdz(2)
  v3(2)  = v(2) + dt*0.25*(g1(2) + g2(2))
  do k= 3, nz-2
    g2(k)  = -f*(u(k)-ug(k))+ w(k)*dv(k)
    v3(k)  = v(k) + dt*0.25*(g1(k)+g2(k))
  enddo
  g2(nz-1)  = -f*(u(nz-1)-ug(nz-1))+ w(nz-1)*dvdz(nz-1)
  v3(nz-1)  = v(nz-1) + dt*0.25*(g1(nz-1)+g2(nz-1))
  g2(nz)  = -f*(u(nz)-ug(nz))+ w(nz)*dvdz(nz)
  v3(nz)  = v(nz) + dt*0.25*(g1(nz)+g2(nz))
  
  call flux(v3, nz, fv, kappa)    !! flux compuation
  do k= 3, nz-2
    dv(k) = -(fv(k)-fv(k-1)) /((z(k+1) - z(k-1))*0.5)
  enddo
  g3(1)  = -f*(u(1)-ug(1))+ w(1)*dvdz(1)
  v(1)  = v(1) + dt*(1.0/6.0*g1(1) + 1.0/6.0*g2(1) +2.0/3.0*g3(1))
  g3(2)  = -f*(u(2)-ug(2))+ w(2)*dvdz(2)
  v(2)  = v(2) + dt*(1.0/6.0*g1(2) + 1.0/6.0*g2(2) + 2.0/3.0*g3(2))
  do k= 3, nz-2
    g3(k)  = -f*(u(k)-ug(k))+ w(k)*dv(k)
    v(k)  = v(k) + dt*(1.0/6.0*g1(k)+1.0/6.0*g2(k)+2.0/3.0*g3(k))
  enddo
  g3(nz-1)  = -f*(u(nz-1)-ug(nz-1))+ w(nz-1)*dvdz(nz-1)
  v(nz-1)  = v(nz-1) + dt*(1.0/6.0*g1(nz-1)+1.0/6.0*g2(nz-1)+2.0/3.0*g3(nz-1) )
  g3(nz)  = -f*(u(nz)-ug(nz))+ w(nz)*dvdz(nz)
  v(nz)  = v(nz) + dt*(1.0/6.0*g1(nz)+1.0/6.0*g2(nz)+2.0/3.0*g3(nz))
  call flux(qv, nz, fqv, kappa)
  do k= 3, nz-2
    dqv(k) = -(fqv(k)-fqv(k-1)) /((z(k+1) - z(k-1))*0.5)
  enddo
  g1(1) = dqvdt(1)+w(1)*dqvdz(1)
  qv2(1) = qv(1)+ dt*g1(1)
  g1(2) = dqvdt(2)+w(2)*dqvdz(2)
  qv2(2) = qv(2)+ dt*g1(2)
  do k= 3, nz-2
    g1(k) = dqvdt(k)+w(k)*dqv(k)
    qv2(k) = qv(k)+ dt*g1(k)
  enddo
  g1(nz-1) = dqvdt(nz-1)+w(nz-1)*dqvdz(nz-1)
  qv2(nz-1) = qv(nz-1)+ dt*g1(nz-1)
  g1(nz) = dqvdt(nz)+w(nz)*dqvdz(nz)
  qv2(nz) = qv(nz)+ dt*g1(nz)

  call flux(qv2, nz, fqv, kappa)
  do k= 3, nz-2
    dqv(k) = -(fqv(k)-fqv(k-1)) /((z(k+1) - z(k-1))*0.5)
  enddo
  g2(1) = dqvdt(1)+w(1)*dqvdz(1)
  qv3(1) = qv(1)+ dt*0.25*(g1(1)+g2(1))
  g2(2) = dqvdt(2)+w(2)*dqvdz(2)
  qv3(2) = qv(2)+ dt*0.25*(g1(2)+g2(2))
  do k= 3, nz-2
    g2(k) = dqvdt(k)+w(k)*dqv(k)
    qv3(k) = qv(k)+ dt*0.25*(g1(k)+g2(k))
  enddo
  g2(nz-1) = dqvdt(nz-1)+w(nz-1)*dqvdz(nz-1)
  qv3(nz-1) = qv(nz-1)+ dt*0.25*(g1(nz-1)+g2(nz-1))
  g2(nz) = dqvdt(nz)+w(nz)*dqvdz(nz)
  qv3(nz) = qv(nz)+ dt*0.25*(g1(nz)+g2(nz))

  call flux(qv3, nz, fqv, kappa)
  do k= 3, nz-2
    dqv(k) = -(fqv(k)-fqv(k-1)) /((z(k+1) - z(k-1))*0.5)
  enddo
  g3(1) = dqvdt(1)+w(1)*dqvdz(1)
  qv(1) = qv(1)+ dt*(1.0/6.0*g1(1)+1.0/6.0*g2(1) + 2.0/3.0*g3(1))
  g3(2) = dqvdt(2)+w(2)*dqvdz(2)
  qv(2) = qv(2)+ dt*(1.0/6.0*g1(2)+1.0/6.0*g2(2)+2.0/3.0*g3(2))
  do k= 3, nz-2
    g3(k) = dqvdt(k)+w(k)*dqv(k)
    qv(k) = qv(k)+ dt*(1.0/6.0*g1(k)+1.0/6.0*g2(k) + 2.0/3.0*g3(k))
  enddo
  g3(nz-1) = dqvdt(nz-1)+w(nz-1)*dqvdz(nz-1)
  qv(nz-1) = qv(nz-1)+ dt*(1.0/6.0*g1(nz-1)+1.0/6.0*g2(nz-1)+2.0/3.0*g3(nz-1))
  g3(nz) = dqvdt(nz)+w(nz)*dqvdz(nz)
  qv(nz) = qv(nz)+ dt*(1.0/6.0*g1(nz)+1.0/6.0*g2(nz)+2.0/3.0*g3(nz))


  call flux(qc, nz, fqc, kappa)
  do k= 3, nz-2
    dqc(k) = -(fqc(k)-fqc(k-1)) /((z(k+1) - z(k-1))*0.5)
  enddo
  g1(1) = w(1)*dqcdz(1)
  qc2(1) = qc(1)+ dt*g1(1)
  g1(2) = w(2)*dqcdz(2)
  qc2(2) = qc(2)+ dt*g1(2)
  do k= 3, nz-2
    g1(k) = w(k)*dqc(k)
    qc2(k) = qc(k)+ dt*g1(k)
  enddo
  g1(nz-1) = w(nz-1)*dqcdz(nz-1)
  qc2(nz-1) = qc(nz-1)+ dt*g1(nz-1)
  g1(nz) = w(nz)*dqcdz(nz)
  qc2(nz) = qc(nz)+ dt*g1(nz)

  call flux(qc2, nz, fqc, kappa)
  do k= 3, nz-2
    dqc(k) = -(fqc(k)-fqc(k-1)) /((z(k+1) - z(k-1))*0.5)
  enddo
  g2(1) = w(1)*dqcdz(1)
  qc3(1) = qc(1)+ dt*0.25*(g1(1)+g2(1))
  g2(2) = w(2)*dqcdz(2)
  qc3(2) = qc(2)+ dt*0.25*(g1(2)+g2(2))
  do k= 3, nz-2
    g2(k) = w(k)*dqc(k)
    qc3(k) = qc(k)+ dt*0.25*(g1(k)+g2(k))
  enddo
  g2(nz-1) = w(nz-1)*dqcdz(nz-1)
  qc3(nz-1) = qc(nz-1)+ dt*0.25*(g1(nz-1)+ g2(nz-1))
  g2(nz) = w(nz)*dqcdz(nz)
  qc3(nz) = qc(nz)+ dt*0.25*(g1(nz)+g2(nz))

  call flux(qc3, nz, fqc, kappa)
  do k= 3, nz-2
    dqc(k) = -(fqc(k)-fqc(k-1)) /((z(k+1) - z(k-1))*0.5)
  enddo
  g3(1) = w(1)*dqcdz(1)
  qc(1) = qc(1)+ dt*(1.0/6.0*g1(1)+1.0/6.0*g2(1)+2.0/3.0*g3(1))
  g3(2) = w(2)*dqcdz(2)
  qc(2) = qc(2)+ dt*(1.0/6.0*g1(2)+1.0/6.0*g2(2)+2.0/3.0*g3(2))
  do k= 3, nz-2
    g3(k) = w(k)*dqc(k)
    qc(k) = qc(k)+ dt*(1.0/6.0*g1(k)+1.0/6.0*g2(k)+2.0/3.0*g3(k))
  enddo
  g3(nz-1) = w(nz-1)*dqcdz(nz-1)
  qc(nz-1) = qc(nz-1)+ dt*(1.0/6.0*g1(nz-1)+ 1.0/6.0*g2(nz-1)+2.0/3.0*g3(nz-1))
  g3(nz) = w(nz)*dqcdz(nz)
  qc(nz) = qc(nz)+ dt*(1.0/6.0*g1(nz)+1.0/6.0*g2(nz)+2.0/3.0*g3(nz))


  call flux(qr, nz, fqr, kappa)
  do k= 3, nz-2
    dqr(k) = -(fqr(k)-fqr(k-1)) /((z(k+1) - z(k-1))*0.5)
  enddo
  g1(1) = w(1)*dqrdz(1)
  qr2(1) = qr(1)+ dt*g1(1)
  g1(2) = w(2)*dqrdz(2)
  qr2(2) = qr(2)+ dt*g1(2)
  do k= 3, nz-2
    g1(k) = w(k)*dqr(k)
    qr2(k) = qr(k)+ dt*g1(k)
  enddo
  g1(nz-1) = w(nz-1)*dqrdz(nz-1)
  qr2(nz-1) = qr(nz-1)+ dt*g1(nz-1)
  g1(nz) = w(nz)*dqrdz(nz)
  qr2(nz) = qr(nz)+ dt*g1(nz)

  call flux(qr2, nz, fqr, kappa)
  do k= 3, nz-2
    dqr(k) = -(fqr(k)-fqr(k-1)) /((z(k+1) - z(k-1))*0.5)
  enddo
  g2(1) = w(1)*dqrdz(1)
  qr3(1) = qr(1)+ dt*0.25*(g1(1)+g2(1))
  g2(2) = w(2)*dqrdz(2)
  qr3(2) = qr(2)+ dt*0.25*(g1(2)+g2(2))
  do k= 3, nz-2
    g2(k) = w(k)*dqr(k)
    qr3(k) = qr(k)+ dt*0.25*(g1(k)+g2(k))
  enddo
  g2(nz-1) = w(nz-1)*dqrdz(nz-1)
  qr3(nz-1) = qr(nz-1)+ dt*0.25*(g1(nz-1)+g2(nz-1))
  g2(nz) = w(nz)*dqrdz(nz)
  qr3(nz) = qr(nz)+ dt*0.25*(g1(nz)+g2(nz))

  call flux(qr3, nz, fqr, kappa)
  do k= 3, nz-2
    dqr(k) = -(fqr(k)-fqr(k-1)) /((z(k+1) - z(k-1))*0.5)
  enddo
  g3(1) = w(1)*dqrdz(1)
  qr(1) = qr(1)+ dt*(1.0/6.0*g1(1)+1.0/6.0*g2(1)+2.0/3.0*g3(1))
  g3(2) = w(2)*dqrdz(2)
  qr(2) = qr(2)+ dt*(1.0/6.0*g1(2)+1.0/6.0*g2(2)+2.0/3.0*g3(2))
  do k= 3, nz-2
    g3(k) = w(k)*dqr(k)
    qr(k) = qr(k)+ dt*(1.0/6.0*g1(k)+1.0/6.0*g2(k)+2.0/3.0*g3(k))
  enddo
  g3(nz-1) = w(nz-1)*dqrdz(nz-1)
  qr(nz-1) = qr(nz-1)+ dt*(1.0/6.0*g1(nz-1)+1.0/6.0*g2(nz-1)+2.0/3.0*g3(nz-1))
  g3(nz) = w(nz)*dqrdz(nz)
  qr(nz) = qr(nz)+ dt*(1.0/6.0*g1(nz)+1.0/6.0*g2(nz)+2.0/3.0*g3(nz))


  call flux(th, nz, fth, kappa)
  do k= 3, nz-2
    dth(k) = -(fth(k)-fth(k-1)) /((z(k+1) - z(k-1))*0.5)
  enddo
  g1(1) = dthdt(1)+w(1)*dthdz(1)
  th2(1) = th(1)+ dt*g1(1)
  g1(2) = dthdt(2)+w(2)*dthdz(2)
  th2(2) = th(2)+ dt*g1(2)
  do k= 3, nz-2
    g1(k) = dthdt(k)+w(k)*dth(k)
    th2(k) = th(k)+ dt*g1(k)
  enddo
  g1(nz-1) = dthdt(nz-1)+w(nz-1)*dthdz(nz-1)
  th2(nz-1) = th(nz-1)+ dt*g1(nz-1)
  g1(nz) = dthdt(nz)+w(nz)*dthdz(nz)
  th2(nz) = th(nz)+ dt*g1(nz)

  call flux(th2, nz, fth, kappa)
  do k= 3, nz-2
    dth(k) = -(fth(k)-fth(k-1)) /((z(k+1) - z(k-1))*0.5)
  enddo
  g2(1) = dthdt(1)+w(1)*dthdz(1)
  th3(1) = th(1)+ dt*0.25*(g1(1)+g2(1))
  g2(2) = dthdt(2)+w(2)*dthdz(2)
  th3(2) = th(2)+ dt*0.25*(g1(2)+g2(2))
  do k= 3, nz-2
    g2(k) = dthdt(k)+w(k)*dth(k)
    th3(k) = th(k)+ dt*0.25*(g1(k)+g2(k))
  enddo
  g2(nz-1) = dthdt(nz-1)+w(nz-1)*dthdz(nz-1)
  th3(nz-1) = th(nz-1)+ dt*0.25*(g1(nz-1)+g1(nz-1))
  g2(nz) = dthdt(nz)+w(nz)*dthdz(nz)
  th3(nz) = th(nz)+ dt*(g1(nz) + g2(nz))

  call flux(th3, nz, fth, kappa)
  do k= 3, nz-2
    dth(k) = -(fth(k)-fth(k-1)) /((z(k+1) - z(k-1))*0.5)
  enddo
  g3(1) = dthdt(1)+w(1)*dthdz(1)
  th(1) = th(1)+ dt*(1.0/6.0*g1(1)+1.0/6.0*g2(1)+2.0/3.0*g3(1))
  thv   = th(1)*(qv(1)+eps0)/(eps0*(1.0_r8+qv(1)))
  p(1)  = p0*(thv*rair*rho(1)/p0)**(cp/cv)
  g3(2) = dthdt(2)+w(2)*dthdz(2)
  th(2) = th(2)+ dt*(1.0/6.0*g1(2)+1.0/6.0*g2(2)+2.0/3.0*g3(2))
  thv   = th(2)*(qv(2)+eps0)/(eps0*(1.0_r8+qv(2)))
  p(2)  = p0*(thv*rair*rho(2)/p0)**(cp/cv)
  do k= 3, nz-2
    g3(k) = dthdt(k)+w(k)*dth(k)
    th(k) = th(k)+ dt*(1.0/6.0*g1(k)+1.0/6.0*g2(k)+2.0/3.0*g3(k))
    thv   = th(k)*(qv(k)+eps0)/(eps0*(1.0_r8+qv(k)))
    p(k)  = p0*(thv*rair*rho(k)/p0)**(cp/cv)
  enddo
  g3(nz-1) = dthdt(nz-1)+w(nz-1)*dthdz(nz-1)
  th(nz-1) = th(nz-1)+ dt*(1.0/6.0*g1(nz-1)+1.0/6.0*g1(nz-1)+2.0/3.0*g3(nz-1))
  thv   = th(nz-1)*(qv(nz-1)+eps0)/(eps0*(1.0_r8+qv(nz-1)))
  !print*, 'th', th
  p(nz-1)  = p0*(thv*rair*rho(nz-1)/p0)**(cp/cv)
  g3(nz) = dthdt(nz)+w(nz)*dthdz(nz)
  th(nz) = th(nz)+ dt*(1.0/6.0*g1(nz) + 1.0/6.0*g2(nz)+2.0/3.0*g3(nz))
  thv   = th(nz)*(qv(nz)+eps0)/(eps0*(1.0_r8+qv(nz)))
  p(nz)  = p0*(thv*rair*rho(nz)/p0)**(cp/cv)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! invalid bomex_type
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  else
    write(*,*) 'ERROR:Invalid bomex_type'
    write(*,*) 'bomex_type must be set to linear, cubic, or weno'
  endif
  
end subroutine bomex_ls_forcing

subroutine set_bomex_forcing(ug,vg,w,dthdt,dqvdt,heat,evap,stress,z,nz)

  implicit none

  integer,                      intent(in)  :: nz ! number of model levels
  real(kind=r8), dimension(nz), intent(in)  :: z  ! model levels [m]

  real(kind=r8), dimension(nz), intent(out) :: w,     &  ! prescribed subsidence [m/s]
                                               ug,    &  ! zonal geostropic wind [m/s]
                                               vg,    &  ! meridional geostropic wind [m/s]
                                               dthdt, &  ! prescribed heating [K/s]
                                               dqvdt     ! prescribed drying  [kg/kg/s]

  real(kind=r8), intent(out) :: heat,  &  ! prescribed sensible heat flux [K*m/s]
                                evap,  &  ! prescribed latent heat flux   [m/s]
                                stress      ! prescribed surface stress     [m*m/s/s]
  real(kind=r8) :: z0, z1, z2

  real(kind=r8), parameter :: w0= 0.0_r8, w1=-0.65_r8,w2=0.0_r8 ! [cm/s]
  real(kind=r8), parameter :: t0=-2.0_r8, t1=-2.0_r8, t2=0.0_r8 ! [K/day]
  real(kind=r8), parameter :: q0=-1.2_r8, q1=-1.2_r8, q2=0.0_r8 ! [1e8/s]
  !real(kind=r8), parameter :: w0= 0.0_r8, w1=-0.0_r8,w2=0.0_r8 ! [cm/s]
  !real(kind=r8), parameter :: t0=-0.0_r8, t1=-0.0_r8, t2=0.0_r8 ! [K/day]
  !real(kind=r8), parameter :: q0=-0.0_r8, q1=-0.0_r8, q2=0.0_r8 ! [1e8/s]

  real(kind=r8) :: zz, zfrac, u0, m
  integer :: k

  heat   = 8.0e-3  ! [K*m/s]
  evap   = 5.2e-5  ! [m/s]
  stress = 0.28**2 ! [m*m/s/s]

  z0 = 0.0_r8

  ! Set geostropic wind to match zonal wind above 700 m.
  m=(-4.61_r8+8.75_r8)/(3000.0_r8-700.0_r8)
  u0=-8.75_r8-m*700.0_r8
  do k=1,nz
    ug(k) = u0 + m*z(k)
    vg(k) = 0.0_r8
  end do

  ! Set prescribed subsidence
  z1=1500.0_r8 
  z2=2100.0_r8
  do k=1,nz
    zz = z(k)
    if(zz<=z1) then
      zfrac    = (zz-z0)/(z1-z0)
      w(k) =  w0 + (w1-w0)*zfrac
    elseif(zz<=z2 .and. zz > z1) then
      zfrac    = (zz-z1)/(z2-z1)
      w(k) =  w1 + (w2-w1)*zfrac
    else
      w(k) =  w2
    end if
    ! convert to [m/s]
    w(k) = w(k)/100.0_r8
  end do


  ! Set prescribed radiative cooling
  z1=1500.0_r8 
  z2=3000.0_r8
  do k=1,nz
    zz = z(k)
    if(zz<=z1) then
      zfrac    = (zz-z0)/(z1-z0)
      dthdt(k) =  t0 + (t1-t0)*zfrac
    elseif(zz<=z2 .and. zz > z1) then
      zfrac    = (zz-z1)/(z2-z1)
      dthdt(k) =  t1 + (t2-t1)*zfrac
    else
      dthdt(k) =  t2
    end if

    ! convert to [K/s]
    dthdt(k) = dthdt(k)/86400.0_r8

  end do

  z1=300.0_r8 
  z2=500.0_r8
  do k=1,nz
    zz = z(k)
    if(zz<=z1) then
      dqvdt(k) =  q0
    elseif(zz<=z2 .and. zz > z1) then
      zfrac    = (zz-z1)/(z2-z1)
      dqvdt(k) =  q1 + (q2-q1)*zfrac
    else
      dqvdt(k) =  q2
    end if

    !convert to mixing ratio tendency [1/s]
    dqvdt(k) = dqvdt(k)/1.0e8
    dqvdt(k) = dqvdt(k)/(1.0_r8-dqvdt(k))
  end do

end subroutine set_bomex_forcing

subroutine bomex_init(u,v,rho,theta,qv,qc,qr,p,z,ps,nz)

  implicit none

  integer, parameter :: r8 = SELECTED_REAL_KIND(8)

  integer,                      intent(in)  :: nz ! number of model levels
  real(kind=r8), dimension(nz), intent(in)  :: z  ! model levels [m]

  real(kind=r8),                intent(out) :: ps ! surface pressure [Pa]
  real(kind=r8), dimension(nz), intent(out) :: u,     &  ! zonal velocity [m/s]
                                               v,     &  ! meridional velocity [m/s]
                                               rho,   &  ! density [kg/m/m/m]
                                               theta, &  ! potential temperautre [K]
                                               qv,    &  ! vapor mixing ratio [kg/kg]
                                               qc,    &  ! cloud water mixing ratio [kg/kg]
                                               qr,    &  ! rain water mixing ratio [kg/kg]
                                               p         ! pressure [Pa]


  real(kind=r8), parameter :: z0=0._r8 , z1=520._r8, z2=1480._r8 , z3=2000._r8
  real(kind=r8), parameter :: qv0=17.0_r8, qv1=16.3_r8, qv2=10.7_r8, qv3=4.2_r8
  real(kind=r8), parameter :: pt0=298.7_r8, pt1=pt0, pt2=302.4_r8, pt3=308.2_r8

  integer :: k
  real(kind=r8) :: zz, zfrac, dz
  real(kind=r8) :: exner, theta_s, exner_s, qv_s

  ps   = 101500.0_r8

  do k=1,nz
    zz = z(k)

    v(k)  = 0.0_r8
    qc(k) = 0.0_r8
    qr(k) = 0.0_r8
    if(zz<700.0) then
      u(k) = -8.75
    else
      u(k) = -8.75+1.8e-3*(zz-700.0)
    end if

    ! qv is in specific humidity [g/kg]
    if(zz<z1) then
      zfrac    = (zz-z0)/(z1-z0)
      theta(k) = pt0+(pt1-pt0)*zfrac
      qv(k)    = qv0+(qv1-qv0)*zfrac
    elseif(zz>=z1 .and. zz<z2) then
      zfrac    = (zz-z1)/(z2-z1)
      theta(k) = pt1+(pt2-pt1)*zfrac
      qv(k)    = qv1+(qv2-qv1)*zfrac
    elseif(zz>=z2 .and. zz<z3) then
      zfrac    = (zz-z2)/(z3-z2)
      theta(k) = pt2+(pt3-pt2)*zfrac
      qv(k)    = qv2+(qv3-qv2)*zfrac
    else
      zfrac    = (zz-z3)
      theta(k) = pt3 + 3.65e-3_r8*zfrac
      qv(k)    = max(qv3 - 1.2e-3_r8*zfrac,qv(k-1))
    end if

    ! Convert qv to [kg/kg]
    qv(k) = qv(k)/1.0e3

    ! Convert qv to mixing ratio
    qv(k) = qv(k)/(1.0 - qv(k))

    ! Compute virtual potential temperature
    theta(k) = theta(k) !*(qv(k)+eps0)/(eps0*(1.0_r8+qv(k)))
  end do

  ! Set surface exner and pressure.
  ! Assume constant mixing ratio in bottom layer.
  exner_s = (ps/p0)**(rair/cp)
  qv_s = qv0/1e3
  qv_s = qv_s/(1.0 - qv_s)

  theta_s = pt0*(qv_s+eps0)/(eps0*(1.0_r8+qv_s))
  exner = exner_s - grav/cp*0.5*(1.0/theta(1) + 1/theta_s)*(z(1)-z0)

  ! Set bottom level variables
  k=1 
  p(k) = p0 * exner**(cp/rair)
  rho(k) = p(k)/(rair*exner*theta(k))

  ! Integrate hydrostatically to get pressure
  do k=2,nz
    dz = z(k) - z(k-1)
    exner = exner - grav/cp*0.5*(1.0/theta(k) + 1/theta(k-1))*dz
    p(k) = p0 * exner**(cp/rair)
    rho(k) = p(k)/(rair*exner*theta(k))
  end do

  ! Convert virtual potential temperature back to dry potential temperature
  do k=1,nz
   theta(k) = theta(k) !/((qv(k)+eps0)/(eps0*(1.0_r8+qv(k))))
  end do

end subroutine bomex_init

subroutine computeflux(f, fluxm, fluxp, n)
!!
!! Calculates upwind fith order WENO fluxes
!! f_{i+1/2}^{-} and f_{i+1/2}^{+}
!! Implementation based on Chi-Wang Shu ``High-order 
!! Finite Difference and Finite Volume WENO Schemes and Discontinuous 
!! Galerkin Methods for CFD".
!! https://doi.org/10.1080/1061856031000104851
!!
!! INPUT
!! f(n): data values to be used to calculate the fuxes
!! n: size of f and number of data values used.
!!
!! OUTPUT
!! fluxm(n): to hold the computed flux f_{i+1/2}^{-} 
!! fluxp(n): to hold the computed flux f_{i+1/2}^{+} 
!!
  integer, intent(in)                       :: n
  real(kind=8), intent(in)                  :: f(n)
  real(kind=8), intent(out)                 :: fluxm(n)
  real(kind=8), intent(out)                 :: fluxp(n)

  integer                                   :: i
  real(kind=8)                              :: eps0 
  real(kind=8)                              :: d0, d1, d2, dd0, dd1, dd2 
  real(kind=8)                              :: beta0(n), beta1(n), beta2(n) 
  real(kind=8)                              :: alpha0(n), alpha1(n), alpha2(n) 
  real(kind=8)                              :: w0(n), w1(n), w2(n) 
  real(kind=8)                              :: ww0(n), ww1(n), ww2(n) 
  real(kind=8)                              :: f0(n), f1(n), f2(n) 

  !!! Initialize variables
  !beta0=0.0
  !beta1=0.0
  !beta2=0.0
  !alpha0=0.0
  !alpha1=0.0
  !alpha2=0.0
  !w0=0.0
  !w1=0.0
  !w2=0.0
  !ww0=0.0
  !ww1=0.0
  !ww2=0.0
  !f0=0.0
  !f1=0.0
  !f2=0.0
  eps0 = 1.0e-16
  do i=3,n-2
    beta0(i) = 13.0/12.0*(f(i)  -2.0*f(i+1)+f(i+2))**2.0 + 1.0/4.0*(3.0*f(i)-4.0*f(i+1)+  f(i+2))**2.0 
    !
    beta1(i) = 13.0/12.0*(f(i-1)-2.0*f(i)  +f(i+1))**2.0 + 1.0/4.0*(f(i-1)         -  f(i+1))**2.0
    !
    beta2(i) = 13.0/12.0*(f(i-2)-2.0*f(i-1)+f(i))**2.0   + 1.0/4.0*(f(i-2)-4.0*f(i-1)+3.0*f(i))**2.0
    !write(*,*) 'i=', i,' beta0(i), beta1(i), beta2(i)'
    !write(*,'(3(1x,E30.15))')  beta0(i), beta1(i), beta2(i)
  enddo


  d0=3.0/10.0; d1=3.0/5.0; d2=1.0/10.0;
  do i=3,n-2
    alpha0(i) = d0 / (beta0(i)+eps0)**2.0
    alpha1(i) = d1 / (beta1(i)+eps0)**2.0
    alpha2(i) = d2 / (beta2(i)+eps0)**2.0
    !write(*,*) 'i=', i,' alpha0(i), alpha1(i), alpha2(i)'
    !write(*,'(3(1x,E30.15))') alpha0(i), alpha1(i), alpha2(i)
  enddo

  do i=3,n-2
    w0(i) =  alpha0(i) / (alpha0(i)+alpha1(i)+alpha2(i) );
    w1(i) =  alpha1(i) / (alpha0(i)+alpha1(i)+alpha2(i) );
    w2(i) =  alpha2(i) / (alpha0(i)+alpha1(i)+alpha2(i) );
  enddo

  dd0=1.0/10.0; dd1=3.0/5.0; dd2=3.0/10.0;
  do i=3,n-2
    alpha0(i) = dd0 / (beta0(i)+eps0)**2.0;
    alpha1(i) = dd1 / (beta1(i)+eps0)**2.0;
    alpha2(i) = dd2 / (beta2(i)+eps0)**2.0;
  enddo

  do i=3,n-2
    ww0(i) =  alpha0(i) / (alpha0(i)+alpha1(i)+alpha2(i));
    ww1(i) =  alpha1(i) / (alpha0(i)+alpha1(i)+alpha2(i));
    ww2(i) =  alpha2(i) / (alpha0(i)+alpha1(i)+alpha2(i));
  enddo


  do i=3,n-2
    f0(i) =  1.0/3.0*f(i)  +5.0/6.0*f(i+1)-1.0/6.0*f(i+2)  !! S0 = {x_{i}, x_{i+1}, x_{i+2}}
    f1(i) = -1.0/6.0*f(i-1)+5.0/6.0*f(i)  +1.0/3.0*f(i+1)  !! S1 = {x_{i-1}, x_{i}, x_{i+1}}
    f2(i) =  1.0/3.0*f(i-2)-7.0/6.0*f(i-1)+11.0/6.0*f(i)   !! S2 = {x_{i-2}, x_{i-1}, x_{i}}
  enddo
  
  do i=3,n-2
    !! f_{i+1/2}^{-}
    fluxm(i) = f0(i)*w0(i) + f1(i)*w1(i) + f2(i)*w2(i);

    !! f_{i+1/2}^{+}
    fluxp(i) = f0(i)*ww0(i) + f1(i)*ww1(i) + f2(i)*ww2(i);
  enddo


end subroutine

subroutine flux(f, n, ff, kappa)
!! This subroutine computes third order flux based on the work by
!! by Hundsdorfer et al. ``A positive finite-difference advection scheme".
!! https://doi.org/10.1006/jcph.1995.1042
!!
!! INPUT: 
!! f: data values to be used to compute fluxes
!! n: number of data points in f
!! kappa: 1 first order central, -1 second oder upwind, 1/3 third order upwind 
!!
!! OUTPUT:
!! ff: flux at cell centered 
  
  integer, intent(in)         :: n
  integer                     :: i
  real(kind=8), intent(in)    :: f(n), kappa
  real(kind=8), intent(out)   :: ff(n-1)
  real(kind=8)                :: lim(n-1), tmp, r(n-1), du0, du1


  !!!--  kappa = 1 second order central, kappa = -1 second order upwind
  !!!--  kappa = 1/3 third order upwind biased
  !!kappa = 1.0/3.0
  !!!-- compute the ratio
  !!!-- r_{i+1/2} = (f_{i+1}-f_{i}) / (f_{i}-f_{i-1})
  !do i=2, n-2
  !  r(i) = (f(i)-f(i+1))/ (f(i+1)-f(i+2) + 1e-30)
  !enddo
 
  !!!-- Compute the limter 
  !!! K(r) = (1-k)/2 + (1+k)/2*r 
  !!! phi_{i+1/2} = max( 0, min(2r, min(2, K(r_{i+1/2})) ) )
  !do i=2, n-2
  !   tmp = (1.0-kappa)/2.0 + (1.0+kappa)/2.0*r(i)
  !   lim(i) = max(0.0, min(2.0*r(i), min(2.0, tmp)))
  !   !tmp = 0.25 + 0.5 * r(i)
  !   !lim(i) = max(min(min(tmp, 4.0), min(tmp, 2*r(i))), 1e-8)
  !enddo

  do i=2, n-2
     du0 = f(i+1)-f(i+2)
     du1 = f(i)-f(i+1)
     call limiter(8, du0, du1, kappa, tmp)
     lim(i) = tmp
  enddo
  !!-- compute the 2 fluxes 
  !! ff_{i+1/2} = f_{i+1} + 0.5 phi_{i+1/2}(f_{i+1}-f_{i+2})
  do i=2, n-2
    ff(i) = f(i+1) + 0.5 * lim(i)*(f(i+1)-f(i+2))
    if(abs(ff(i)) < 1.0e-20 .and. ff(i) .ne. 0.0)then
      ff(i) = 0.0
    endif
  enddo

end subroutine

subroutine limiter( i, du0, du1, kappa, phi)
!!
!! limiter function
!! value depend on value of i
!! i = 1 first order  
!! i = 2 second order central
!! i = 3 van leer harnonic
!! i = 4 third order unlimited
!! i = 4 CCCT
!!

integer, intent(in)         :: i
real(kind=8), intent(in)    :: du0, du1, kappa
real(kind=8), intent(out)   :: phi
real(kind=8)                :: tol, r, tmp



tol = 1.0e-8
if(i==1) then
   phi = 0.0
elseif(i == 2) then
   phi = 1.0
endif 
if(i == 3) then
   r = du1 / ( du0 + 1.0e-26)
   phi = ( r + abs(r) ) / ( 1.0 + abs(r))
elseif (i == 4) then
   r = du1 / ( du0 + 1.0e-26)
   phi = 0.25 + 0.75*r
endif
if (i == 5) then
   r = du1 / ( du0 + 1.0e-26)
   phi = 0.25 + 0.75*r
    if(phi > 4.0) then
       phi = 4.0
   endif
   if (phi > (2.0*r)) then
      phi = 2.0*r
   endif
   if(phi < tol) then 
      phi = 0.0
   endif
 endif
if (i == 6) then
   r = du1 / ( du0 + 1.0e-26)
   phi = r;
    if (phi > 1.0) then
       phi = 1.0
   endif
   if (phi < tol) then 
      phi = 0.0
   endif
 endif
if (i == 7) then
   r = du1 / ( du0 + 1.0e-26)
   phi = 0.5 + 0.5*r
    if (phi > 2.0) then
       phi = 2.0
   endif
   if (phi > (2*r)) then
      phi = 2*r;
   endif
    if (phi > 2.0) then
       phi = 2.0
   endif
   if (phi < tol) then
      phi = 0.0
   endif
 endif
 if(i == 8)then
   r = du1 / ( du0 + 1.0e-26)
   tmp = (1.0-kappa)/2.0 + (1.0+kappa)/2.0*r
   phi = max(0.0, min(2.0*r, min(2.0, tmp)))
 endif

end subroutine  

end module mod_bomex
