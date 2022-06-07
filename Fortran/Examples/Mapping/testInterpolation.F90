program testInterpolation
!!
!!
!!
  implicit none
  integer 		:: nz(3)
  integer 		:: k

  nz = (/64, 127, 253/)
  do k=1,3
    call mapping(nz(k))
  enddo
   
end program testInterpolation



subroutine mapping(nz)
!!
!! This subrtouine is used to set up the mapping for the Runge and TWP-ICE examples.
!! The following files below are required for the experiment.
!! 'zd_qc_qv_pres_u_v_T_zp_127' and 'zd_qc_qv_pres_u_v_T_zp_253' are obtained by fitting
!! 'zd_qc_qv_pres_u_v_T_zp_64' using a radial basis function interpolation and the evaluating
!! the fitted function at the desired points.
!!
!! FILES
!! 'mapping_data/zd_qc_qv_pres_u_v_T_zp_64': obtained directly from TWP-ICE simulation at
!!   at t = XX s.
!! 'mapping_data/zd_qc_qv_pres_u_v_T_zp_127': obtained by adding at point at the center of each interval 
!! 'mapping_data/zd_qc_qv_pres_u_v_T_zp_253': obtained by adding 3 uniformly spaced points inside each 
!!   interval.
!!  
!! INPUT
!! nz: number of point to be used for the Runge and TWP-ICE examples 
!!
  use mod_adaptiveInterpolation

  implicit none

  integer, parameter    :: niter = 1          !! number of iteration
  integer               :: nz  				!! number of point used 64 127 253
  integer               :: d(15)			!! polynomial degree used 
  integer               :: i, j , k                     !! iteration ideces
  real(kind=8)          :: zd(nz),   zp(nz)             !! uniform and LGL mesh
  real(kind=8)          :: qcp(nz),   qcp2(nz),   qc(nz),   qc2(nz)
  real(kind=8)          :: qvp(nz),   qvp2(nz),   qv(nz),   qv2(nz)
  real(kind=8)          :: presp(nz), presp2(nz), pres(nz), pres2(nz)
  real(kind=8)          :: up(nz),    up2(nz),    u(nz),    u2(nz)                 !! data on uniform and LGL mesh
  real(kind=8)          :: vp(nz),    vp2(nz),    v(nz),    v2(nz)                 !! data on uniform and LGL mesh
  real(kind=8)          :: Tp(nz),    Tp2(nz),    T(nz),    T2(nz)                 !! data on uniform and LGL mesh
  character*32          :: name_u, name_v, name_T, name_pres, name_qc, name_qv

  !!!! Varibales needed for Runge and Heavisisde functions !!!!
  character*32          :: name_runge, name_heaviside
  real(kind=8)          :: zd_runge(nz),   zp_runge(nz)                 !! uniform and LGL mesh
  real(kind=8)          :: zd_heaviside(nz),   zp_heaviside(nz)                 !! uniform and LGL mesh
  real(kind=8)          :: rungep(nz), rungep2(nz), runge(nz), runge2(nz)                 !! data on uniform and LGL mesh
  real(kind=8)          :: heavisidep(nz), heavisidep2(nz), heaviside(nz), heaviside2(nz)                 !! data on uniform and LGL mesh
  real(kind=8)          :: dx, a_heaviside, b_heaviside, a_runge , b_runge 
    

  !!** Initialize data **!!
  qcp2=0.0; qvp2=0.0; up2=0.0; vp2=0.0; presp2=0.0; Tp2=0.0;

  !!** Read input data from file  **!!
  if(nz .eq. 64) then
  open(100, file='mapping_data/zd_qc_qv_pres_u_v_T_zp_64',  status='old')
  elseif(nz .eq. 127) then
  open(100, file='mapping_data/zd_qc_qv_pres_u_v_T_zp_127',  status='old')
  elseif(nz .eq. 253) then
  open(100, file='mapping_data/zd_qc_qv_pres_u_v_T_zp_253',  status='old')
  endif
  do i=1, nz
    !!read(100, *) zd(i), qc2(i), qv2(i), pres2(i), u2(i), v2(i), T2(i), &
    !!             zp(i), qcp2(i), qvp2(i), presp2(i), up2(i), vp2(i), Tp2(i)
    read(100, *) zd(i), qc2(i), qv2(i), pres2(i), u2(i), v2(i), T2(i), runge2(i), &
                 zp(i), qcp2(i), qvp2(i), presp2(i), up2(i), vp2(i), Tp2(i),  rungep2(i)

  enddo
  close(100) 

  !!** Initialize polynomial degree to be used **!!
  d = (/1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16/)

  !!** Initialize variables names **!!
  name_runge = 'runge'
  name_qc = 'qc'
  !!name_heaviside= 'heaviside'
  !!name_u = 'u'
  !!name_v = 'v'
  !!name_T = 'T'
  !!name_pres = 'pres'
  !!name_qv = 'qv'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Examples using runge function and heaviside function !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!! Map zd and zp  to xd_xx an xp_xx respectively 
  !!a_heaviside = -0.2
  !!b_heaviside = 0.2
  a_runge = -1.0
  b_runge = 1.0
  call scaleab(zd, zd_runge, nz, zd(1), zd(nz), a_runge, b_runge)
  call scaleab(zp, zp_runge, nz, zd(1), zd(nz), a_runge, b_runge)
  !!call scaleab(zd, zd_heaviside, nz, zd(1), zd(nz), a_heaviside, b_heaviside)
  !!call scaleab(zp, zp_heaviside, nz, zd(1), zd(nz), a_heaviside, b_heaviside)

  !! Evaluate coresponding values at x 
  dx = 0.01 !! dummy variables not used for function evaluation
  do i=1, nz
    call evalFun1D(1, zd_runge(i), runge2(i), dx)
    call evalFun1D(1, zp_runge(i), rungep2(i), dx)
    !!call evalFun1D(2, zd_heaviside(i), heaviside2(i), dx)
    !!call evalFun1D(2, zp_heaviside(i), heavisidep2(i), dx)
  enddo

  do i=1, 15
    runge = runge2
    rungep = rungep2
    qc=qc2
    qcp=qcp2
    !!heaviside = heaviside2
    !!heavisidep = heavisidep2
    write(*, *) '********** d= ', d(i), '**********'
    call massConservation2(nz, zd_runge, runge, zp_runge, rungep, niter, d(i), name_runge)
    call massConservation2(nz, zd, qc, zp, qcp, niter, d(i), name_qc)
    !!call massConservation2(nz, zd, runge, zp, rungep, niter, d(i), name_runge)
    !!call massConservation2(nz, zd_heaviside, heaviside, zp_heaviside, heavisidep, niter, d(i), name_heaviside)
  enddo


  !!do i=1, 15
  !!  qc=qc2
  !!  qv=qv2
  !!  pres=pres2
  !!  u=u2
  !!  v=v2
  !!  T=T2    
  !!  qcp=qcp2
  !!  qvp=qvp2
  !!  presp=presp2
  !!  up=up2
  !!  vp=vp2
  !!  Tp=Tp2
  !!  write(*, *) '********** d= ', d(i), '**********'
  !!  call massConservation2(nz, zd, u, zp, up, niter, d(i), name_u)
  !!  !!call massConservation2(nz, zd_runge, u, zp_runge, up, niter, d(i), name_u)
  !!  call massConservation2(nz, zd, v, zp, vp, niter, d(i), name_v)
  !!  call massConservation2(nz, zd, pres, zp, presp, niter, d(i), name_pres)
  !!  call massConservation2(nz, zd, qc, zp, qcp, niter, d(i), name_qc)
  !!  call massConservation2(nz, zd, qv, zp, qvp, niter, d(i), name_qv)
  !!  call massConservation2(nz, zd, T, zp, Tp, niter, d(i), name_T)
  !!enddo

end subroutine 

subroutine massConservation2(nz, zd, u, zp, u2, niter, dd, profile_name)
!!
!!
!!
  use mod_legendre
  use mod_adaptiveInterpolation

  implicit none

  integer, parameter    :: d = 4                !! degree
  integer, parameter    :: nzplot = 2000                  !! number of elements
  integer, intent(in)   :: nz               !! number of points
  integer, intent(in)   :: dd               !! number of points
  integer, intent(in)   :: niter           !! number of iteration

  real(kind=8), intent(in)  :: zd(nz), zp(nz)                 !! uniform and LGL mesh
  real(kind=8), intent(in)  :: u(nz), u2(nz)
  character*12, intent(in)  :: profile_name

  integer               :: ne                  !! number of elements
  real(kind=8)          :: upplot(nzplot), upplot_pchip(nzplot), upplot_makima(nzplot)                   !! number of elements
  real(kind=8)          :: udplot(nzplot), udplot_pchip(nzplot), udplot_makima(nzplot)                    !! number of elements
  real(kind=8)          :: zplot(nzplot), upplot_dbi(nzplot), udplot_dbi(nzplot)                   !! number of elements
  integer               :: fun                  !! determine which fucntion isused
  integer               :: i, j , k             !! iteration ideces
  integer               :: is, ie
  integer               :: limiter
  integer 		:: sten
  integer               :: iter
  real(kind=8)          :: up(nz), ud(nz)     !! data on uniform and LGL mesh
  real(kind=8)          :: up_pchip(nz), ud_pchip(nz)     !! data on uniform and LGL mesh
  real(kind=8)          :: up_dbi(nz), ud_dbi(nz)     !! data on uniform and LGL mesh
  real(kind=8)          :: up_makima(nz), ud_makima(nz)     !! data on uniform and LGL mesh
  real(kind=8)          :: udt(nz)
  real(kind=8)          :: tmp1(nz-1), tmp2(nz-1)
  integer               :: deg(nz-1), deg_dbi(nz-1)
  real(kind=8)          :: qx(d+1), qw(d+1)     !! reference quadrature points and weigths 
  real(kind=8)          :: ux(d+1)              !! reference points for uniform mesh 
  real(kind=8)          :: nw(d+1)
  real(kind=8)          :: psi(d+1, d+1), dpsi(d+1, d+1), eps0

  !!real(kind=8)          :: a, b                 !! interval of interest [a, b] 
  real(kind=8)          :: dz, xl, xr
  real(kind=8)          :: m_uniform((nz-1)/d), m_lgl((nz-1)/d)
  real(kind=8)          :: m_uniform_dbi((nz-1)/d), m_lgl_dbi((nz-1)/d)
  real(kind=8)          :: m_uniform_pchip((nz-1)/d), m_lgl_pchip((nz-1)/d)
  real(kind=8)          :: m_uniform_makima((nz-1)/d), m_lgl_makima((nz-1)/d)
  real(kind=8)          :: c(d+1)

  real(kind=8)			:: wk(nz*2), d_tmp(nz)
  real(kind=8)			:: fdl(nz)
  real(kind=8)			:: fdl2(2000)
  logical                       :: spline
  integer                       :: nwk, ierr
  integer                       :: fnumber
  character*80                  :: fname, tmp_str


  !!** To save data **!!
  real(kind=8)          :: ud_out(nz, niter+2), up_out(nz, niter+2)
  real(kind=8)          :: ud_pchip_out(nz, niter+2), up_pchip_out(nz, niter+2)
  real(kind=8)          :: ud_dbi_out(nz, niter+2), up_dbi_out(nz, niter+2)
  real(kind=8)          :: ud_makima_out(nz, niter+2), up_makima_out(nz, niter+2)
  real(kind=8)          :: s(nz-1), s_0, s_nz
  real(kind=8)          :: sp(nz-1), sp_0, sp_nz
  real(kind=8)          :: s_table(nz, niter+1), s_table_d(nz, niter+1)
  integer               :: extrema_interval_d(nz, niter), extrema_interval_p(nz, niter)
  integer               :: deg_ud_out(nz-1, niter+2), deg_up_out(nz-1, niter+2)
  integer               :: deg_ud_dbi_out(nz-1, niter+2), deg_up_dbi_out(nz-1, niter+2)
  real(kind=8)          :: u_table(nz-1, dd+1), x_table(nz-1, dd+1), lambda_table(nz-1, dd+1), sigma_table(nz-1, dd+1)
  real(kind=8)          :: up_table(nz-1, dd+1), xp_table(nz-1, dd+1), lambdap_table(nz-1, dd+1), sigmap_table(nz-1, dd+1)
  real(kind=8)          :: prod_sigma_table(nz-1, dd+1)
  real(kind=8)          :: prod_sigmap_table(nz-1, dd+1)
  real(kind=8)          :: uumin(nz-1), uumax(nz-1)


  spline = .false.
  nwk = nz*2
  eps0 = 0.01
  sten = 3
 
  !!** Initialize file number **!!
  if (nz .le. 1000) then
    fnumber = dd*1000 + nz
    write(*,*) 'The file number is ', fnumber
  else
    write(*,*) 'ERROR: file not set for input size nz lager than 999 '
    call exit(1)
  endif
 

  !!** Initialize variables **!!
  sp = 0.0; sp_0 =0.0; sp_nz=0;
  s =0.0; s_0=0.0;s_nz=0.0
  deg= 0; deg_dbi=0
  m_uniform = 0.0; m_uniform_pchip = 0.0; m_uniform_dbi = 0.0; m_uniform_makima = 0.0
  m_lgl= 0.0; m_lgl_pchip = 0.0; m_lgl_dbi = 0.0; m_lgl_makima = 0.0
  deg_ud_out = 0; deg_up_out=0; deg_ud_dbi_out=0; deg_up_dbi_out=0

  ud = u
  ud_pchip = u
  ud_dbi = u
  ud_makima = u

  !!** save Initial data **!!
  do i=1, nz
    ud_out(i,1) = zd(i)
    up_out(i,1) = zp(i)
    ud_pchip_out(i,1) = zd(i)
    up_pchip_out(i,1) = zp(i)
    ud_makima_out(i,1) = zd(i)
    up_makima_out(i,1) = zp(i)
    ud_dbi_out(i,1) = zd(i)
    up_dbi_out(i,1) = zp(i)
  enddo
  do i=1, nz-1
    deg_ud_out(i,1) = 0
    deg_up_out(i,1) = 0
    deg_ud_dbi_out(i,1) = 0
    deg_up_dbi_out(i,1) = 0
  enddo

  !!** Calculate the number of elements **!!
  ne = (nz-1) / d


  !!** Set limiter **!!
  limiter = 2

  !!** Set up mesh for plotting **!!
  dz = (zd(nz)-zd(1)) / 1999.0
  zplot(1) = zd(1)
  do i=2, nzplot-1
    zplot(i) = zplot(i-1)+dz
  enddo
  zplot(nzplot) = zd(nz)
  !!** LGL reference poinst for each element **!!
  call legendre_gauss_lobatto(d+1, qx, qw)           
  
  !!** Open file **!!
  write(tmp_str, '("integral", i5.5)')fnumber
  !!write(fname,*) trim(profile_name)//trim(tmp_str)
  fname = trim(profile_name)//trim(tmp_str)
  write(*,*) 'file to save integral is ', fname
   
  !!==open(100,file=fname, status='unknown')
  !!==do iter=1, niter+1
  iter=1
    !!** save data on physics grid **!!
    if(iter <= niter+1) then 
      do i=1, nz
       up_out(i,iter+1) = u2(i)
       up_pchip_out(i,iter+1) = u2(i)
       up_dbi_out(i,iter+1) = u2(i)
       up_makima_out(i,iter+1) = u2(i)
      enddo

      do i=1, nz-1
       deg_up_out(i,iter+1) = 0
       deg_up_dbi_out(i,iter+1) = 0
      enddo
    endif
   

    !!==!!** Calculate mass using quadrature rule **!!
    !!==is = 1
    !!==ie = 1
    !!==do i=1, ne
    !!==  is = ie
    !!==  ie = is + d
    !!==  m_lgl(i) = 0.0
    !!==  m_lgl_pchip(i) = 0.0
    !!==  m_lgl_dbi(i) = 0.0
    !!==  m_lgl_makima(i) = 0.0
    !!==  do j=0, d
    !!==    m_lgl(i) = m_lgl(i) + qw(j+1)*ud(is+j)
    !!==    m_lgl_pchip(i) = m_lgl_pchip(i) + qw(j+1)*ud_pchip(is+j)
    !!==    m_lgl_dbi(i) = m_lgl_dbi(i) + qw(j+1)*ud_dbi(is+j)
    !!==    m_lgl_makima(i) = m_lgl_makima(i) + qw(j+1)*ud_makima(is+j)
    !!==  enddo
    !!==  m_lgl(i) = (zd(ie)-zd(is)) /2.0 * m_lgl(i)
    !!==  m_lgl_pchip(i) = (zd(ie)-zd(is)) /2.0 * m_lgl_pchip(i)
    !!==  m_lgl_dbi(i) = (zd(ie)-zd(is)) /2.0 * m_lgl_dbi(i)
    !!==  m_lgl_makima(i) = (zd(ie)-zd(is)) /2.0 * m_lgl_makima(i)
    !!==enddo
    write(*,*) '------ iter =', iter, ' --------'
    !write(*,*) 'Mass LGL mesh  =', sum(m_lgl), sum(m_lgl_pchip), sum(m_lgl_dbi), sum(m_lgl_makima) 

    !!** Interpolate from Dynamics to Physics **!!
    !!call adaptiveInterpolation1D(zd, ud_dbi, nz, zp, up_dbi, nz, dd, 1, deg_dbi) 
    call adaptiveInterpolation1D(zd, ud_dbi, nz, zp, up_dbi, nz, dd, 1, deg_dbi) 
    !!call adaptiveInterpolation1D(zd, ud_dbi, nz, zplot, upplot_dbi, nzplot, dd, 1, deg_dbi) 
    
    !!call adaptiveInterpolation1D(zd, ud, nz, zp, up, nz, dd, 2, deg, &
    call adaptiveInterpolation1d(zd, ud, nz, zp, up, nz, dd, 2, deg, &
                                       sten, eps0, uumin, uumax, &
                                       x_table, u_table, lambda_table, &
                                       sigma_table, prod_sigma_table )
    !!call adaptiveInterpolation1d(zd, ud, nz, zplot, upplot, nzplot, dd, 2, deg)

    !call pchez(nz, zd, ud_pchip, d_tmp, spline, wk, nwk, ierr)
    !call pchev(nz, zd, ud_pchip, d_tmp, nz, zp, up_pchip, fdl, ierr)
    call pchez(nz, zd, ud_pchip, d_tmp, spline, wk, nwk, ierr)
    call pchev(nz, zd, ud_pchip, d_tmp, nz, zp, up_pchip, fdl, ierr)

    !!call pchez(nz, zd, ud_pchip, d_tmp, spline, wk, nwk, ierr)
    !!call pchev(nz, zd, ud_pchip, d_tmp, nzplot, zplot, upplot_pchip, fdl2, ierr)

    !!call makima(zd, ud_makima, nz, zp, up_makima, nz)
    call makima(zd, ud_makima, nz, zp, up_makima, nz)
    !!call makima(zd, ud_makima, nz, zplot, upplot_makima, nzplot)

    !!** save data that is on dynamics grid**!!
    if(iter <= niter+1) then 
      do i=1, nz
       ud_out(i,iter+1) = ud(i) 
       ud_pchip_out(i,iter+1) = ud_pchip(i)
       ud_dbi_out(i,iter+1) = ud_dbi(i)
       ud_makima_out(i,iter+1) = ud_makima(i)
      enddo
      !!
      do i=1, nz-1
       deg_ud_out(i,iter+1) = deg(i) 
       deg_ud_dbi_out(i,iter+1) = deg_dbi(i)
      enddo
    endif

    iter = 2
    !!** save data on physics grid **!!
    if(iter <= niter+1) then 
      do i=1, nz
       up_out(i,iter+1) = up(i)
       up_pchip_out(i,iter+1) = up_pchip(i)
       up_dbi_out(i,iter+1) = up_dbi(i)
       up_makima_out(i,iter+1) = up_makima(i)
      enddo

      !!do i=1, nz-1
      !! deg_up_out(i,iter+1) = deg(i)
      !! deg_up_dbi_out(i,iter+1) = deg_dbi(i)
      !!enddo
    endif
 
    !!==!!** Calculate mass using newton cotes rules **!!
    !!==is = 1
    !!==ie = 1
    !!==do i=1, ne
    !!==  is = ie
    !!==  ie = is + d
    !!==  !!** evaluate Lagrange basis functions at quadrature points **!!
    !!==  call lagrange_basis_eval(d+1, zp(is:ie), d+1, zd(is:ie), psi)

    !!==  !!** Calculate newton cote weights for each element **!!
    !!==  nw = 0.0
    !!==  do j=1, d+1 
    !!==    do k=1, d+1
    !!==      nw(j) = nw(j) + psi(j, k)* qw(k)
    !!==    enddo
    !!==  enddo
    !!==  nw = (zp(ie)-zp(is))/2.0 * nw
    !!==   
    !!==  !!** Calculate integral fo each element **!!
    !!==  m_uniform(i) = 0
    !!==  m_uniform_pchip(i) = 0
    !!==  m_uniform_dbi(i) = 0
    !!==  m_uniform_makima(i) = 0
    !!==  do j=0, d
    !!==   m_uniform(i) = m_uniform(i) + nw(j+1)*up(is+j)
    !!==   m_uniform_pchip(i) = m_uniform_pchip(i) + nw(j+1)*up_pchip(is+j)
    !!==   m_uniform_dbi(i) = m_uniform_dbi(i) + nw(j+1)*up_dbi(is+j)
    !!==  enddo 
    !!==enddo

    !!==!write(*,*) 'Mass uniform mesh  =', sum(m_uniform), sum(m_uniform_pchip), sum(m_uniform_dbi), sum(m_uniform_makima)
    !!==write(100,'(8(1x,E30.15))') sum(m_lgl_pchip), sum(m_lgl_makima), sum(m_lgl_dbi), sum(m_lgl), &
    !!==                           sum(m_uniform_pchip), sum(m_uniform_makima), sum(m_uniform_dbi), sum(m_uniform)  

    !!** Interpolate from Physics to Dynamics **!!
    call adaptiveInterpolation1D(zp, up_dbi, nz, zd(2:nz-1), ud_dbi(2:nz-1), nz-2, dd, 1, deg_dbi) 
    !!call adaptiveInterpolation1D(zp, up_dbi, nz, zplot, udplot_dbi, nzplot, dd, 1, deg_dbi) 
    !!call adaptiveInterpolation1D(zp, u2, nz, zd, ud_dbi, nz, dd, 1, deg_dbi) 

    !!call adaptiveInterpolation1d(zp, u2, nz, zd, ud, nz, dd, 2, deg, & 
    call adaptiveInterpolation1D(zp, up, nz, zd(2:nz-1), ud(2:nz-1), nz-2, dd, 2, deg, & 
                                       sten, eps0, uumin, uumax, &
                                       xp_table, up_table, lambdap_table,&
                                       sigmap_table, prod_sigmap_table )
    !!call adaptiveInterpolation1d(zp, up, nz, zplot, udplot, nzplot, dd, 2, deg) 

    call pchez(nz, zp, up_pchip, d_tmp, spline, wk, nwk, ierr)
    call pchev(nz, zp, up_pchip, d_tmp, nz-2, zd(2:nz-1), ud_pchip(2:nz-1), fdl, ierr)

    !!call pchez(nz, zp, up_pchip, d_tmp, spline, wk, nwk, ierr)
    !!call pchev(nz, zp, up_pchip, d_tmp, nzplot, zplot, udplot_pchip, fdl2, ierr)

    !!call makima(zp, up_makima, nz, zd, ud_makima, nz)
    call makima(zp, up_makima, nz, zd(2:nz-1), ud_makima(2:nz-1), nz-2)
    !!call makima(zp, up_makima, nz, zplot, udplot_makima, nzplot)

     iter = 2
    !!** save data on physics grid **!!
    if(iter <= niter+1) then 
      !!do i=1, nz
      !! up_out(i,iter+1) = up(i)
      !! up_pchip_out(i,iter+1) = up_pchip(i)
      !! up_dbi_out(i,iter+1) = up_dbi(i)
      !! up_makima_out(i,iter+1) = up_makima(i)
      !!enddo

      do i=1, nz-1
       deg_up_out(i,iter+1) = deg(i)
       deg_up_dbi_out(i,iter+1) = deg_dbi(i)
      enddo
    endif
  
    !!** save data that is on dynamics grid**!!
    if(iter <= niter+1) then 
      do i=1, nz
       ud_out(i,iter+1) = ud(i) 
       ud_pchip_out(i,iter+1) = ud_pchip(i)
       ud_dbi_out(i,iter+1) = ud_dbi(i)
       ud_makima_out(i,iter+1) = ud_makima(i)
      enddo
      !!
      do i=1, nz-1
       deg_ud_out(i,iter+1) = 0 
       deg_ud_dbi_out(i,iter+1) = 0
      enddo
    endif


   
  !!==enddo !! end of iter
  !!==close(100)


  !do i=1, nz+1
  !  write(*,'(4(1x,E30.15))') real(i), s_table(i, 1), s_table(i, 2), s_table(i, 1)
  !enddo
  !do i=1, niter
  !  iter = i
  !  write(*, '(57(1x, I5))') iter, (extrema_interval_d(j, i), j=11, nz-1)
  !enddo
  !!** write data to file **!!
    !!** Runge function **!
 
  !!** Write PPI  to different filea where
  !!   -- the first letter indicate the name of the profile
  !!   -- PPI, DBI, PCHIP, MAKIMA indicate the interpolation method
  !!   -- the 2 digit after the letter indicate the target polynomial degree
  !!   -- the last 3 digit indicate the total number point use for each method  **!!
  write(tmp_str, '("dPLOT", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nzplot
    write(100, '(9(1x,E30.15))') zplot(i), udplot_pchip(i), udplot_makima(i), udplot_dbi(i), &
                       udplot(i), upplot_pchip(i), upplot_makima(i), upplot_dbi(i), upplot(i)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname
 
  write(tmp_str, '("dPPI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(3(1x,E30.15))') (ud_out(i, j), j=1, 3)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname
  write(tmp_str, '("dErrPPI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(3(1x,E30.15))') (ud_out(i, j)-u(i), j=1, 3)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname
 
  write(tmp_str, '("dDEGPPI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz-1
    write(100, '(1002(1x,I2))') (deg_ud_out(i, j), j=1, niter+2)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname

  !!** To study divided differences **!!
  write(tmp_str, '("dXTABPPI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz-1
    write(100, '(17(1x,E30.15))') (x_table(i, j), j=1,dd+1)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname

  write(tmp_str, '("dUTABPPI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz-1
    write(100, '(17(1x,E30.15))') (u_table(i, j), j=1,dd+1)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname

  write(tmp_str, '("dLTABPPI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz-1
    write(100, '(17(1x,E30.15))') (lambda_table(i, j), j=1,dd+1)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname

  write(tmp_str, '("dSTABPPI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz-1
    write(100, '(17(1x,E30.15))') (sigma_table(i, j), j=1,dd+1)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname

  write(tmp_str, '("dErrTABPPI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz-1
    write(100, '(17(1x,E30.15))') (prod_sigma_table(i, j), j=1,dd+1)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname


  write(tmp_str, '("pPPI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.15))') (up_out(i, j), j=1, niter+2)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname

  write(tmp_str, '("pErrPPI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.15))') (up_out(i, j)-u2(i), j=1, niter+2)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname


  write(tmp_str, '("pDEGPPI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz-1
    write(100, '(1002(1x,I2))') (deg_up_out(i, j), j=1, niter+2)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname

  !!** To study divided differences **!!
  write(tmp_str, '("pXTABPPI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz-1
    write(100, '(17(1x,E30.15))') (xp_table(i, j), j=1,dd+1)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname
  write(tmp_str, '("pUTABPPI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz-1
    write(100, '(17(1x,E30.15))') (up_table(i, j), j=1,dd+1)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname

  write(tmp_str, '("pLTABPPI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz-1
    write(100, '(17(1x,E30.15))') (lambdap_table(i, j), j=1,dd+1)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname

  write(tmp_str, '("pSTABPPI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz-1
    write(100, '(17(1x,E30.15))') (sigmap_table(i, j), j=1,dd+1)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname

  write(tmp_str, '("pErrTABPPI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz-1
    write(100, '(17(1x,E30.15))') (prod_sigmap_table(i, j), j=1,dd+1)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname


  write(tmp_str, '("dDBI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.15))') (ud_dbi_out(i, j), j=1, niter+2)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname
  write(tmp_str, '("dErrDBI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.15))') (ud_dbi_out(i, j)-u(i), j=1, niter+2)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname
 
  write(tmp_str, '("dDEGDBI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz-1
    write(100, '(1002(1x,I2))') (deg_ud_dbi_out(i, j), j=1, niter+2)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname

  write(tmp_str, '("pDBI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.15))') (up_dbi_out(i, j), j=1, niter+2)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname

  write(tmp_str, '("pErrDBI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.15))') (up_dbi_out(i, j)-u2(i), j=1, niter+2)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname


  write(tmp_str, '("pDEGDBI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)
  open(102,file=fname, status='unknown')
  do i=1, nz-1
    write(102, '(1002(1x,I2))') (deg_up_dbi_out(i, j), j=1, niter+2)
  enddo
  close(102)
  write(*,*) 'Saved in file name ', fname

  write(tmp_str, '("dPCHIP", i5.5)') 3*1000+nz
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.15))') (ud_pchip_out(i, j), j=1, niter+2)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname

  write(tmp_str, '("dErrPCHIP", i5.5)') 3*1000+nz
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.15))') (ud_pchip_out(i, j), u(i), j=1, niter+2)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname
 
  write(tmp_str, '("pPCHIP", i5.5)') 3*1000+nz
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.15))') (up_pchip_out(i, j), j=1, niter+2)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname

  write(tmp_str, '("pErrPCHIP", i5.5)') 3*1000+nz
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.15))') (up_pchip_out(i, j)-u2(i), j=1, niter+2)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname
 
  write(tmp_str, '("dMAKIMA", i5.5)') 3*1000+nz
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.15))') (ud_makima_out(i, j), j=1, niter+2)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname

  write(tmp_str, '("dErrMAKIMA", i5.5)') 3*1000+nz
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.15))') (ud_makima_out(i, j)-u(i), j=1, niter+2)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname
 
  write(tmp_str, '("pMAKIMA", i5.5)') 3*1000+nz
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.15))') (up_makima_out(i, j), j=1, niter+2)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname

  write(tmp_str, '("pErrMAKIMA", i5.5)') 3*1000+nz
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.15))') (up_makima_out(i, j)-u2(i), j=1, niter+2)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname

  print*, 'DONE SAVING FILES'
end subroutine massConservation2


subroutine evalFun1D(fun, x, v, h)
!!
!! To evaluate different fuunction. fun determine
!! the function that will be evaluated at the points x
!!
  implicit none 

  integer, intent(in)           :: fun                  !! function type
  real(kind=8), intent(in)      :: x                    !! point
  real(kind=8), intent(in)      :: h                    !! elements size
  real(kind=8), intent(out)     :: v                    !! point
  real(kind=8)                  :: pi, k, t, delta,a,b  !! temporary variables
  integer                       :: i, j, ne           
  
  !!** intialize variables **!!
  k = 100
  pi = 4.0*atan(1.0) 
  t = x/h
  delta = 0.01
  !!** a and b are only used for fun =5. When using fun=4 the 
  !!   a and b must be the same as the ones Paperexample1D **!!
  a = -2.0
  b = 0.0
  ne = int((b-a)/h)
  

  !!!** 1D runge function **!!
  if(fun .eq. 1) then
    !!v = 0.1 / (0.1 + 25.0 * x * x)
    v = 1.0 / (1.0 + 25.0 * x * x)

  !!** heaviside function **!!
  else if(fun .eq. 2)then
    v = 1.0/(1.0 + exp(-2*k*x))

  !!** Gelb and Tanner function **!!
  else if(fun .eq. 3)then
    if(x < -0.5) then
      v = 1.0 + (2.0* exp(2.0 * pi *(x+1.0)) - 1.0 -exp(pi)) /(exp(pi)-1.0)
    else
      v = 1.0 - sin(2*pi*x / 3.0 + pi/3.0)
    endif

  !!!!** modified square function **!!
  else if(fun .eq. 4)then
    v = 1.0 - abs( 2.0/ pi * atan( sin(pi*t) /delta))


  !!** Modified tanh function **!!
  else if(fun .eq. 5)then
    k = 10.0
    if(a <= x .and. x <= a+h)then
      v = tanh(x*k)
    elseif(a+h<= x .and. x <= a+2*h)then
      v = 2*tanh(x*k)         - tanh((a+h)*k)
    endif
    do i=3, ne 
      !!** 
      if(a +(i-1)*h <= x .and. x <= a+i*h)then
        v = i*tanh(x*k) - tanh((a+h)*k)
        do j=2, i-1 
         v = v - tanh( (a+(j*h))*k )
        enddo 
      endif
    enddo
    v = 1.0 + v

  !!!!* sin function **!!
  else if(fun .eq. 6)then
    v = 1.0 + sin(x*pi) 
  endif

end subroutine evalFun1D


subroutine f ( n, t, y, ydot )

!*****************************************************************************80
!
!! F evaluates the right hand sides of the ODE's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) t
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) ydot(n)

  ydot(1) = y(2)
  ydot(2) = -y(1)

  return
end


