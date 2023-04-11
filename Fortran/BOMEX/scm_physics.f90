module scm_physics

use machine,  only : r8 => kind_phys
use physcons, only : rair => con_rd, &
                     grav => con_g,  &
                     cp   => con_cp

implicit none

contains

subroutine run_scm_physics(u,v,th,qv,qc,qr,pr,pri,rho, &
                           heat,evap,stress,z,zi,delt,nz)

integer,                        intent(in)    :: nz
real(kind=r8), dimension(nz),   intent(inout) :: u, v, th, qv, qc, qr, pr, rho
real(kind=r8), dimension(nz),   intent(in)    :: z
real(kind=r8), dimension(nz+1), intent(in)    :: pri, zi
real(kind=r8),                  intent(in)    :: delt

real(kind=r8), intent(in)    :: evap, heat, stress

real(kind=r8), parameter :: p00=1.0e5

integer, parameter :: ix     = 1,    &
                      im     = 1,    &
                      ntrac  = 2,    &
                      ntcw   = ntrac

logical, parameter :: dspheat=.false., &
                      lprnt  =.false.

real(kind=r8), dimension(nz)  :: pk, t!, r

! Input variables to moninedmf
real(kind=r8), dimension(ix,nz)       :: u1, v1, t1, swh, hlw, del, prsl, &
                                         prslk, phil
real(kind=r8), dimension(ix,nz+1)     :: prsi,phii
real(kind=r8), dimension(ix,nz,ntrac) :: q1
real(kind=r8), dimension(ix)          :: xmu, psk, rbsoil, zorl, u10m, &
                                         v10m, tsea, qss, heat0, evap0, stress0, &
                                         spd1, fm, fh

integer,       dimension(ix)          :: kinver
integer                               :: ipr
real(kind=r8)                         :: xkzm_m,xkzm_h,xkzm_s 

! Input/Output variables to moninedmf
real(kind=r8), dimension(ix,nz)       :: du, dv, tau
real(kind=r8), dimension(ix,nz,ntrac) :: rtg

! Output variables to moninedmf
real(kind=r8), dimension(ix)          :: dusfc, dvsfc, dtsfc, dqsfc, &
                                         hpbl, hgamt, hgamq
real(kind=r8), dimension(ix,nz-1)     :: dkt!,dku
integer                               :: kpbl

integer :: k

real :: precl

!! Added so the input matches matches the dummy variables !!
real(kind=r8)    :: heat_in(ix), evap_in(ix), stress_in(ix)
integer    :: kpbl_in(ix)


! Set background diffusion to 0 for scm
kinver=0
xkzm_m=0.0_r8
xkzm_h=0.0_r8
xkzm_s=0.0_r8

spd1(ix)   = sqrt(u(1)**2+v(1)**2)
fm(ix)     = 1.0_r8
fh(ix)     = 1.0_r8
rbsoil(ix) = 0.0_r8
u10m(ix)   = 0.0_r8
v10m(ix)   = 0.0_r8
tsea(ix)   = 300.4_r8
qss(ix)    = 0.02245_r8

heat0(ix) = heat
evap0(ix) = evap
stress0(ix)= stress

do k=1,nz+1
  prsi(ix,k)  = pri(k)
  phii(ix,k)  = zi(k)*grav
end do

do k=1,nz
  pk(k) = (pr(k)/p00)**(rair/cp)
  t(k)  = th(k)*pk(k)
end do

do k=1,nz
  u1(ix,k)    = u(k)
  v1(ix,k)    = v(k)
  q1(ix,k,1)  = qv(k)
  q1(ix,k,2)  = qc(k)
  prsl(ix,k)  = pr(k)
  prslk(ix,k) = pk(k) 
  t1(ix,k)    = t(k)
  del(ix,k)   = prsi(ix,k) - prsi(ix,k+1)
  phil(ix,k)  = z(k)*grav
end do

psk(ix) = (1000.0/1000.0)**(rair/cp)

! Set radiative heating to 0
do k=1,nz
  swh(ix,k) = 0.0_r8
  hlw(ix,k) = 0.0_r8
end do
xmu(ix) = 0.0_r8

! Initialize tendencies to 0
do k=1,nz
  du(ix,k)    = 0.0_r8
  dv(ix,k)    = 0.0_r8
  tau(ix,k)   = 0.0_r8
  rtg(ix,k,1) = 0.0_r8
  rtg(ix,k,2) = 0.0_r8
end do

heat_in(ix) = heat
evap_in(ix) = evap
stress_in(ix) = stress
kpbl_in(ix) = kpbl
call moninedmf(ix,im,nz,ntrac,ntcw,dv,du,tau,rtg,            &
               u1,v1,t1,q1,swh,hlw,xmu,                      &
               psk,rbsoil,zorl,u10m,v10m,fm,fh,              &
               tsea,qss,heat_in,evap_in,stress_in,spd1,kpbl_in, &
               prsi,del,prsl,prslk,phii,phil,delt,dspheat,   &
               dusfc,dvsfc,dtsfc,dqsfc,hpbl,hgamt,hgamq,dkt, &
               kinver,xkzm_m,xkzm_h,xkzm_s,lprnt,ipr)

do k=1,nz
  !print '(I2,5(3x,E12.5))',k,du(ix,k),dv(ix,k),tau(ix,k),rtg(ix,k,1),rtg(ix,k,2)
  u(k) = u(k) + delt*du(ix,k)
  v(k) = v(k) + delt*dv(ix,k)
  t(k) = t(k) + delt*tau(ix,k)
  qv(k)= qv(k)+ delt*rtg(ix,k,1)
  qc(k)= qc(k)+ delt*rtg(ix,k,2)
end do

do k=1,nz
  !pr(k) = rho(k)*rair*t(k)
  !pk(k) = (pr(k)/p00)**(rair/cp)
  th(k) = t(k)/pk(k)
end do
call kessler(th,qv,qc,qr,rho,pk,delt,z,nz,precl)

end subroutine run_scm_physics

end module scm_physics
