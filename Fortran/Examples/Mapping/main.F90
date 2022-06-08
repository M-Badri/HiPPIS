program main
!!
!! Driver to produce the mapping result based on the
!! Runge and the TWP-ICE.
!!
  implicit none
  integer 		:: nz(3)
  integer 		:: k

  nz = (/64, 127, 253/)
  do k=1,3
    call mapping(nz(k))
  enddo
   
end program 



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
  integer               :: d(3)			!! polynomial degree used 
  integer               :: i, j , k                     !! iteration ideces
  real(kind=8)          :: zd(nz),   zp(nz)             !! uniform and LGL mesh
  real(kind=8)          :: qcp(nz),   qcp2(nz),   qc(nz),   qc2(nz)

  character*32          :: name_runge, name_qc
  real(kind=8)          :: zd_runge(nz),   zp_runge(nz)                 !! uniform and LGL mesh
  real(kind=8)          :: rungep(nz), rungep2(nz), runge(nz), runge2(nz)                 !! data on uniform and LGL mesh
  real(kind=8)          :: dx, a_runge , b_runge 
    

  !!** Read input data from file  **!!
  if(nz .eq. 64) then
  open(100, file='mapping_data/zd_qc_zp_64',  status='old')
  elseif(nz .eq. 127) then
  open(100, file='mapping_data/zd_qc_zp_127',  status='old')
  elseif(nz .eq. 253) then
  open(100, file='mapping_data/zd_qc_zp_253',  status='old')
  endif
  do i=1, nz
    read(100, *) zd(i), qc2(i), zp(i), qcp2(i)
  enddo
  close(100) 

  !!** Initialize polynomial degree to be used **!!
  d = (/3, 5, 7/)

  !!** Initialize variables names **!!
  name_runge = 'runge'
  name_qc = 'qc'

  !!!! Map zd and zp  to xd_xx an xp_xx respectively 
  a_runge = -1.0
  b_runge = 1.0
  call scaleab(zd, zd_runge, nz, zd(1), zd(nz), a_runge, b_runge)
  call scaleab(zp, zp_runge, nz, zd(1), zd(nz), a_runge, b_runge)

  !! Evaluate coresponding values at x 
  dx = 0.01 !! dummy variables not used for function evaluation
  do i=1, nz
    call evalFun1D(1, zd_runge(i), runge2(i), dx)
    call evalFun1D(1, zp_runge(i), rungep2(i), dx)
  enddo

  do i=1, 3
    runge = runge2
    rungep = rungep2
    qc=qc2
    qcp=qcp2
    write(*, *) '********** d= ', d(i), '**********'
    call mapping2(nz, zd_runge, runge, zp_runge, rungep, niter, d(i), name_runge)
    call mapping2(nz, zd, qc, zp, qcp, niter, d(i), name_qc)
  enddo


end subroutine 

subroutine mapping2(nz, zd, u, zp, u2, niter, dd, profile_name)
!!
!!
!! Subroutine for mapping data form mesh points zd to zp and back to zp
!!
!! INPUT
!! nz: number of points
!! zd: first mesh points (dynamics mesh points)
!! u: data values associated with the first mesh
!! zp: fecond mesh points (physics mesh points)
!! u2: data values associated with the second mesh
!! dd: maximum degree used fr each interpolant
!! profile_name: profile name to be sused to save results
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

  iter=1
  !!** save data on physics grid **!!
  do i=1, nz
   up_out(i,iter+1) = u2(i)
   up_pchip_out(i,iter+1) = u2(i)
   up_dbi_out(i,iter+1) = u2(i)
   up_makima_out(i,iter+1) = u2(i)
  enddo
  !!
  do i=1, nz-1
   deg_up_out(i,iter+1) = 0
   deg_up_dbi_out(i,iter+1) = 0
  enddo
   

  !!** Mapping data values from zd (dynamics mesh) to zp (physics mesh) using DBI  **!!
  call adaptiveInterpolation1D(zd, ud_dbi, nz, zp, up_dbi, nz, dd, 1, deg_dbi) 
  
  !!** Mapping data values from zd (dynamics mesh) to zp (physics mesh) using PPI  **!!
  call adaptiveInterpolation1d(zd, ud, nz, zp, up, nz, dd, 2, deg, &
                                     sten, eps0, uumin, uumax, &
                                     x_table, u_table, lambda_table, &
                                     sigma_table, prod_sigma_table )

  !!** Mapping data values from zd (dynamics mesh) to zp (physics mesh) using PCHIP  **!!
  call pchez(nz, zd, ud_pchip, d_tmp, spline, wk, nwk, ierr)
  call pchev(nz, zd, ud_pchip, d_tmp, nz, zp, up_pchip, fdl, ierr)


  !!call makima(zd, ud_makima, nz, zp, up_makima, nz)

  !!** save data that is on dynamics grid**!!
  do i=1, nz
   ud_out(i,iter+1) = ud(i) 
   ud_pchip_out(i,iter+1) = ud_pchip(i)
   ud_dbi_out(i,iter+1) = ud_dbi(i)
   !!ud_makima_out(i,iter+1) = ud_makima(i)
  enddo
  !!
  do i=1, nz-1
   deg_ud_out(i,iter+1) = deg(i) 
   deg_ud_dbi_out(i,iter+1) = deg_dbi(i)
  enddo

  iter = 2
  !!** save data on physics grid **!!
  do i=1, nz
   up_out(i,iter+1) = up(i)
   up_pchip_out(i,iter+1) = up_pchip(i)
   up_dbi_out(i,iter+1) = up_dbi(i)
   !!up_makima_out(i,iter+1) = up_makima(i)
  enddo


  !!** Mapping data values from zp (physics mesh) to zd (dynamics mesh)  using DBI  **!!
  call adaptiveInterpolation1D(zp, up_dbi, nz, zd(2:nz-1), ud_dbi(2:nz-1), nz-2, dd, 1, deg_dbi) 

  !!** Mapping data values from zp (physics mesh) to zd (dynamics mesh)  using PPI  **!!
  call adaptiveInterpolation1D(zp, up, nz, zd(2:nz-1), ud(2:nz-1), nz-2, dd, 2, deg, & 
                                     sten, eps0, uumin, uumax, &
                                     xp_table, up_table, lambdap_table,&
                                     sigmap_table, prod_sigmap_table )

  !!** Mapping data values from zp (physics mesh) to zd (dynamics mesh)  using PCHIP  **!!
  call pchez(nz, zp, up_pchip, d_tmp, spline, wk, nwk, ierr)
  call pchev(nz, zp, up_pchip, d_tmp, nz-2, zd(2:nz-1), ud_pchip(2:nz-1), fdl, ierr)

  !!call makima(zp, up_makima, nz, zd(2:nz-1), ud_makima(2:nz-1), nz-2)

  iter = 2
  !!** save data on physics grid **!!
  do i=1, nz-1
    deg_up_out(i,iter+1) = deg(i)
    deg_up_dbi_out(i,iter+1) = deg_dbi(i)
  enddo
  
  !!** save data that is on dynamics grid **!!
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

 
  !!** Write PPI  to different filea where
  !!   -- the first letter indicate the name of the profile
  !!   -- PPI, DBI, PCHIP, MAKIMA indicate the interpolation method
  !!   -- the 2 digit after the letter indicate the target polynomial degree
  !!   -- the last 3 digit indicate the total number point use for each method  **!!
  !!!write(tmp_str, '("dPLOT", i5.5)')fnumber
  !!!fname = trim(profile_name)//trim(tmp_str)
  !!!open(100,file=fname, status='unknown')
  !!!do i=1, nzplot
  !!!  write(100, '(9(1x,E30.15))') zplot(i), udplot_pchip(i), udplot_makima(i), udplot_dbi(i), &
  !!!                     udplot(i), upplot_pchip(i), upplot_makima(i), upplot_dbi(i), upplot(i)
  !!!enddo
  !!!close(100)
  !!!write(*,*) 'Saved in file name ', fname
 
  write(tmp_str, '("dPPI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(3(1x,E30.15))') (ud_out(i, j), j=1, 3)
  enddo
  close(100)

  !!!write(*,*) 'Saved in file name ', fname
  !!!write(tmp_str, '("dErrPPI", i5.5)')fnumber
  !!!fname = trim(profile_name)//trim(tmp_str)
  !!!open(100,file=fname, status='unknown')
  !!!do i=1, nz
  !!!  write(100, '(3(1x,E30.15))') (ud_out(i, j)-u(i), j=1, 3)
  !!!enddo
  !!!close(100)
  !!!write(*,*) 'Saved in file name ', fname
 
  write(tmp_str, '("dDEGPPI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz-1
    write(100, '(1002(1x,I2))') (deg_ud_out(i, j), j=1, niter+2)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname

  !!!!!** To study divided differences **!!
  !!!write(tmp_str, '("dXTABPPI", i5.5)')fnumber
  !!!fname = trim(profile_name)//trim(tmp_str)
  !!!open(100,file=fname, status='unknown')
  !!!do i=1, nz-1
  !!!  write(100, '(17(1x,E30.15))') (x_table(i, j), j=1,dd+1)
  !!!enddo
  !!!close(100)
  !!!write(*,*) 'Saved in file name ', fname

  !!!write(tmp_str, '("dUTABPPI", i5.5)')fnumber
  !!!fname = trim(profile_name)//trim(tmp_str)
  !!!open(100,file=fname, status='unknown')
  !!!do i=1, nz-1
  !!!  write(100, '(17(1x,E30.15))') (u_table(i, j), j=1,dd+1)
  !!!enddo
  !!!close(100)
  !!!write(*,*) 'Saved in file name ', fname

  !!!write(tmp_str, '("dLTABPPI", i5.5)')fnumber
  !!!fname = trim(profile_name)//trim(tmp_str)
  !!!open(100,file=fname, status='unknown')
  !!!do i=1, nz-1
  !!!  write(100, '(17(1x,E30.15))') (lambda_table(i, j), j=1,dd+1)
  !!!enddo
  !!!close(100)
  !!!write(*,*) 'Saved in file name ', fname

  !!!write(tmp_str, '("dSTABPPI", i5.5)')fnumber
  !!!fname = trim(profile_name)//trim(tmp_str)
  !!!open(100,file=fname, status='unknown')
  !!!do i=1, nz-1
  !!!  write(100, '(17(1x,E30.15))') (sigma_table(i, j), j=1,dd+1)
  !!!enddo
  !!!close(100)
  !!!write(*,*) 'Saved in file name ', fname

  !!!write(tmp_str, '("dErrTABPPI", i5.5)')fnumber
  !!!fname = trim(profile_name)//trim(tmp_str)
  !!!open(100,file=fname, status='unknown')
  !!!do i=1, nz-1
  !!!  write(100, '(17(1x,E30.15))') (prod_sigma_table(i, j), j=1,dd+1)
  !!!enddo
  !!!close(100)
  !!!write(*,*) 'Saved in file name ', fname


  write(tmp_str, '("pPPI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.15))') (up_out(i, j), j=1, niter+2)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname

  !!!write(tmp_str, '("pErrPPI", i5.5)')fnumber
  !!!fname = trim(profile_name)//trim(tmp_str)
  !!!open(100,file=fname, status='unknown')
  !!!do i=1, nz
  !!!  write(100, '(1002(1x,E30.15))') (up_out(i, j)-u2(i), j=1, niter+2)
  !!!enddo
  !!!close(100)
  !!!write(*,*) 'Saved in file name ', fname


  write(tmp_str, '("pDEGPPI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz-1
    write(100, '(1002(1x,I2))') (deg_up_out(i, j), j=1, niter+2)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname

  !!!!!** To study divided differences **!!
  !!!write(tmp_str, '("pXTABPPI", i5.5)')fnumber
  !!!fname = trim(profile_name)//trim(tmp_str)
  !!!open(100,file=fname, status='unknown')
  !!!do i=1, nz-1
  !!!  write(100, '(17(1x,E30.15))') (xp_table(i, j), j=1,dd+1)
  !!!enddo
  !!!close(100)
  !!!write(*,*) 'Saved in file name ', fname
  !!!write(tmp_str, '("pUTABPPI", i5.5)')fnumber
  !!!fname = trim(profile_name)//trim(tmp_str)
  !!!open(100,file=fname, status='unknown')
  !!!do i=1, nz-1
  !!!  write(100, '(17(1x,E30.15))') (up_table(i, j), j=1,dd+1)
  !!!enddo
  !!!close(100)
  !!!write(*,*) 'Saved in file name ', fname

  !!!write(tmp_str, '("pLTABPPI", i5.5)')fnumber
  !!!fname = trim(profile_name)//trim(tmp_str)
  !!!open(100,file=fname, status='unknown')
  !!!do i=1, nz-1
  !!!  write(100, '(17(1x,E30.15))') (lambdap_table(i, j), j=1,dd+1)
  !!!enddo
  !!!close(100)
  !!!write(*,*) 'Saved in file name ', fname

  !!!write(tmp_str, '("pSTABPPI", i5.5)')fnumber
  !!!fname = trim(profile_name)//trim(tmp_str)
  !!!open(100,file=fname, status='unknown')
  !!!do i=1, nz-1
  !!!  write(100, '(17(1x,E30.15))') (sigmap_table(i, j), j=1,dd+1)
  !!!enddo
  !!!close(100)
  !!!write(*,*) 'Saved in file name ', fname

  !!!write(tmp_str, '("pErrTABPPI", i5.5)')fnumber
  !!!fname = trim(profile_name)//trim(tmp_str)
  !!!open(100,file=fname, status='unknown')
  !!!do i=1, nz-1
  !!!  write(100, '(17(1x,E30.15))') (prod_sigmap_table(i, j), j=1,dd+1)
  !!!enddo
  !!!close(100)
  !!!write(*,*) 'Saved in file name ', fname


  write(tmp_str, '("dDBI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.15))') (ud_dbi_out(i, j), j=1, niter+2)
  enddo
  close(100)

  !!!write(*,*) 'Saved in file name ', fname
  !!!write(tmp_str, '("dErrDBI", i5.5)')fnumber
  !!!fname = trim(profile_name)//trim(tmp_str)
  !!!open(100,file=fname, status='unknown')
  !!!do i=1, nz
  !!!  write(100, '(1002(1x,E30.15))') (ud_dbi_out(i, j)-u(i), j=1, niter+2)
  !!!enddo
  !!!close(100)
  !!!write(*,*) 'Saved in file name ', fname
 
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

  !!!write(tmp_str, '("pErrDBI", i5.5)')fnumber
  !!!fname = trim(profile_name)//trim(tmp_str)
  !!!open(100,file=fname, status='unknown')
  !!!do i=1, nz
  !!!  write(100, '(1002(1x,E30.15))') (up_dbi_out(i, j)-u2(i), j=1, niter+2)
  !!!enddo
  !!!close(100)
  !!!write(*,*) 'Saved in file name ', fname


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

  !!!write(tmp_str, '("dErrPCHIP", i5.5)') 3*1000+nz
  !!!fname = trim(profile_name)//trim(tmp_str)
  !!!open(100,file=fname, status='unknown')
  !!!do i=1, nz
  !!!  write(100, '(1002(1x,E30.15))') (ud_pchip_out(i, j), u(i), j=1, niter+2)
  !!!enddo
  !!!close(100)
  !!!write(*,*) 'Saved in file name ', fname
 
  write(tmp_str, '("pPCHIP", i5.5)') 3*1000+nz
  fname = trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.15))') (up_pchip_out(i, j), j=1, niter+2)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname

  !!!write(tmp_str, '("pErrPCHIP", i5.5)') 3*1000+nz
  !!!fname = trim(profile_name)//trim(tmp_str)
  !!!open(100,file=fname, status='unknown')
  !!!do i=1, nz
  !!!  write(100, '(1002(1x,E30.15))') (up_pchip_out(i, j)-u2(i), j=1, niter+2)
  !!!enddo
  !!!close(100)
  !!!write(*,*) 'Saved in file name ', fname
 
  !!!write(tmp_str, '("dMAKIMA", i5.5)') 3*1000+nz
  !!!fname = trim(profile_name)//trim(tmp_str)
  !!!open(100,file=fname, status='unknown')
  !!!do i=1, nz
  !!!  write(100, '(1002(1x,E30.15))') (ud_makima_out(i, j), j=1, niter+2)
  !!!enddo
  !!!close(100)
  !!!write(*,*) 'Saved in file name ', fname

  !!!write(tmp_str, '("dErrMAKIMA", i5.5)') 3*1000+nz
  !!!fname = trim(profile_name)//trim(tmp_str)
  !!!open(100,file=fname, status='unknown')
  !!!do i=1, nz
  !!!  write(100, '(1002(1x,E30.15))') (ud_makima_out(i, j)-u(i), j=1, niter+2)
  !!!enddo
  !!!close(100)
  !!!write(*,*) 'Saved in file name ', fname
 
  !!!write(tmp_str, '("pMAKIMA", i5.5)') 3*1000+nz
  !!!fname = trim(profile_name)//trim(tmp_str)
  !!!open(100,file=fname, status='unknown')
  !!!do i=1, nz
  !!!  write(100, '(1002(1x,E30.15))') (up_makima_out(i, j), j=1, niter+2)
  !!!enddo
  !!!close(100)
  !!!write(*,*) 'Saved in file name ', fname

  !!!write(tmp_str, '("pErrMAKIMA", i5.5)') 3*1000+nz
  !!!fname = trim(profile_name)//trim(tmp_str)
  !!!open(100,file=fname, status='unknown')
  !!!do i=1, nz
  !!!  write(100, '(1002(1x,E30.15))') (up_makima_out(i, j)-u2(i), j=1, niter+2)
  !!!enddo
  !!!close(100)
  !!!write(*,*) 'Saved in file name ', fname

  print*, 'DONE SAVING FILES'

end subroutine


subroutine evalFun1D(fun, x, v, h)
!!
!! To evaluate different function. fun determine
!! the function that will be evaluated at the points x
!!

  implicit none 

  integer, intent(in)                    :: fun                  !! function type
  real(kind=8), intent(in)               :: x                    !! point
  real(kind=8), intent(in), optional     :: h                    !! elements size
  real(kind=8), intent(out)              :: v                    !! point
  real(kind=8)                           :: pi, k, t, delta,a,b  !! temporary variables
  integer                                :: i, j, ne           
  
  
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


