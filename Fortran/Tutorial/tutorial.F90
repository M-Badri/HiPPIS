program tutorial
!!---------------------------------------------------------------------------------------------!
!! Tutorial showing how to use the 1D, 2D, and 3D DBI and PPI interpolation
!!---------------------------------------------------------------------------------------------!

  use mod_adaptiveInterpolation

  implicit none
 
  integer                               :: i, j, k
  integer                               :: sten
  integer                               :: d
  integer                               :: interpolation_type
  integer, parameter                    :: n=17
  integer, parameter                    :: m=50
  real(kind=dp)                          :: x(n), y(n), z(n)
  real(kind=dp)                          :: v(n), v2D(n,n), v3D(n,n,n)
  real(kind=dp)                          :: xout(m), yout(m), zout(m)
  real(kind=dp)                          :: vout_apprx(m), vout_apprx2D(m,m), vout_apprx3D(m,m,m)
  real(kind=dp)                          :: eps0, eps1 
  real(kind=dp)                          :: dx 
  real(kind=dp)                          :: pi 


  
  !-- 1D Tutorial -- !
  dx = 2.0_dp/real(n-1, kind=dp) 
  do i=1, n-1
    x(i) = -1.0_dp + real(i-1, kind=dp)*dx  ! input mesh points
    v(i) = 0.1_dp / (0.1_dp + 25_dp*x(i)**2)                   ! input data values
  enddo
  x(n) = 2.0_dp
  v(n) =  0.1_dp / (0.1_dp + 25_dp*x(n)**2); 
  dx = 2.0_dp/real(m-1, kind=dp) 
  do i=1, m-1
    xout(i) = -1.0_dp + real(i-1, kind=dp)*dx  ! output points
  enddo
  xout(m)= 2.0_dp
  
  d = 8                           ! target and maximum polynomial degree used for each interval
  interpolation_type = 2          ! 1 for DBI and 2 for PPI
  sten = 1_dp                     ! optional parameter to guide stencil selection 1, 2, and 3
  eps0 = 0.01_dp                  ! optional positive parameter to bound interpolant in PPI
  eps1 = 1.0_dp                   ! optional positive parameter to bound interpolant in PPI

  call adaptiveInterpolation1D(x, v, n, xout, vout_apprx, m, d, interpolation_type, sten, eps0, eps1 ) 
  

  !-- 2D Tutorial -- !
  y = x                             
  do j=1,n
    do i=1,n
      v2D(i,j) = 0.1_dp / (0.1_dp + 25_dp*(x(i)+y(j))**2)                   ! input data values
    enddo
  enddo
  yout = xout
  d = 8                             ! target and maximum polynomial degree used for each interval
  interpolation_type = 2            ! 1 for DBI and 2 for PPI
  sten = 1_dp                       ! optional parameter to guide stencil selection 1, 2, and 3
  eps0 = 0.01_dp                    ! optional positive parameter to bound interpolant in PPI
  eps1 = 1.0_dp                     ! optional positive parameter to bound interpolant in PPI

  call adaptiveInterpolation2D(x, y, n, n, v2D, xout, yout, m, m, vout_apprx2D, d, interpolation_type, sten, eps0, eps1 ) 
  

  !-- 3D Tutorial -- !
  z = x
  do k=1,n
    do j=1,n
      do i=1,n
        v3D(i,j,k) = 0.1_dp / (0.1_dp + 25_dp*(x(i)+y(j)+z(k))**2)                   ! input data values
      enddo
    enddo
  enddo
  zout =xout
  d = 8                             ! target and maximum polynomial degree used for each interval
  interpolation_type = 2            ! 1 for DBI and 2 for PPI
  sten = 1_dp                       ! optional parameter to guide stencil selection 1, 2, and 3
  eps0 = 0.01_dp                    ! optional positive parameter to bound interpolant in PPI
  eps1 = 1.0_dp                     ! optional positive parameter to bound interpolant in PPI

  call adaptiveInterpolation3D(x, y, z, n, n, n, v3D, xout, yout, zout, m, m, m, vout_apprx3D, &
                                         d, interpolation_type, sten, eps0, eps1 ); 
 
end program


