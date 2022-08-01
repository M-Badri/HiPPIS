program testing
!
!
!

  call test1()
  !call test2()
  !call test3()
  !call test4()
end program
!---------------------------------------------------------------------------------------------!
! test1
!---------------------------------------------------------------------------------------------!

subroutine test1()
!
!  test considering small intervals for input and output mesh
!

  use mod_adaptiveInterpolation

  implicit none
  
  integer, parameter		::n = 1e+6
  integer, parameter		::m = 1e+7
  integer       		::i, j
  integer       		::sten
  integer, parameter		::d =8
  real(kind=8), parameter	:: a = -1.0e-10;
  real(kind=8), parameter       :: b = 1.0e-10;
  real(kind=8), parameter       :: eps0=1.0;
  real(kind=8), parameter       :: eps1=1.0;
  real(kind=8)  	 	:: x(n), v1D(n) 
  real(kind=8)  	 	:: xout(m), v1Dout(m) 
  real(kind=8) 			:: dx

  dx = (b-a)/real(n-1, kind=8)
  x(1) = a
  do i=2,n-1
    x(i) = x(i-1)+dx
  enddo
  x(n) = b
  dx = (b-a)/real(m-1, kind=8)
  xout(1) = a
  do i=2, m-1
    xout(i) =xout(i-1)+dx
  enddo
  xout(m) = b
  do i=1, n
    v1D(i) = 0.1 /(0.1 + 25*(x(i)*1e+10)**2);
  enddo
  do  j=1,3
    sten = j
    call adaptiveInterpolation1D(x, v1D, n, xout, v1Dout, m, d, 2, sten, eps0, eps1) 
  enddo
end subroutine


!function test2()
!!
!!  test considering small example with negative values
!!
!
!  n = 20;
!  m = 1000;
!  a = -1.0;
!  b = 1.0;
!  x= linspace(a, b, n);
!  xout= linspace(a, b, m);
!  v1D = zeros(n,1);
!  v1Dout = zeros(m,1);
!  eps0= 1;
!  eps1= 1;
!  dxn = 1; ! dummy variable 
!  d = 8;
!  for j=1:3
!    sten = j
!    for i=1:n
!      v1D(i) = sin(x(i)*pi);
!    end
!    
!    v1Dout = adaptiveInterpolation1D(x, v1D, xout, d, 2, sten, eps0, eps1); 
!
!    figure
!    plot(x, v1D, '*', xout, v1Dout)
!    legend('data', 'apprx')
!    pause
!  end
!end 
!
!
!function test3()
!!
!!  test considering small constant fucntions
!!
!
!  n = 20;
!  m = 1000;
!  a = -1.0;
!  b = 1.0;
!  x= linspace(a, b, n);
!  xout= linspace(a, b, m);
!  v1D = zeros(n,1);
!  v1Dout = zeros(m,1);
!  eps0= 1;
!  eps1= 1;
!  dxn = 1; ! dummy variable 
!  d = 8;
!  for j=1:3
!    sten = j
!    for i=1:n
!      v1D(i) = 1.0;
!    end
!    v1Dout = adaptiveInterpolation1D(x, v1D, xout, d, 2, sten, eps0, eps1); 
!    figure
!    plot(x, v1D, '*', xout, v1Dout)
!    legend('data', 'apprx')
!    pause
!  end
!end 
!
!function test4()
!!
!!  test considering a linear function constant fucntions
!!
!
!  n = 20;
!  m = 1000;
!  a = -1.0;
!  b = 1.0;
!  x= linspace(a, b, n);
!  xout= linspace(a, b, m);
!  v1D = zeros(n,1);
!  v1Dout = zeros(m,1);
!  eps0= 1;
!  eps1= 1;
!  dxn = 1; ! dummy variable 
!  d = 8;
!  for j=1:3
!    sten = j
!    for i=1:n
!      v1D(i) = (-1.0/(b-a))*(x(i)+1) + 1;
!    end
!    v1Dout = adaptiveInterpolation1D(x, v1D, xout, d, 2, sten, eps0, eps1); 
!    figure
!    plot(x, v1D, '*', xout, v1Dout)
!    legend('data', 'apprx')
!    pause
!  end
!end 
!
!
!
!!---------------------------------------------------------------------------------------------!
!! Tutorial showing how to use the 1D, 2D, and 3D DBI and PPI interpolation
!!---------------------------------------------------------------------------------------------!
!
!  !use mod_adaptiveInterpolation
!
!  !implicit none
! 
!  !integer				:: i, j, k
!  !integer				:: sten
!  !integer				:: d
!  !integer				:: interpolation_type
!  !integer, parameter			:: n=17
!  !integer, parameter			:: m=100
!  !real(kind=8)				:: x(n), y(n), z(n)
!  !real(kind=8)				:: v(n), v2D(n,n), v3D(n,n,n)
!  !real(kind=8)				:: xout(m), yout(m), zout(m)
!  !real(kind=8)				:: vout_apprx(m), vout_apprx2D(m,m), vout_apprx3D(m,m,m)
!  !real(kind=8)				:: eps0, eps1 
!  !real(kind=8)				:: dx 
!  !real(kind=8)				:: pi 
!
!
!  !
!  !!-- 1D Tutorial -- !
!  !pi = atan(1.0)*4
!  !dx = 2.0*pi/real(n-1, kind=8) 
!  !do i=1, n-1
!  !  x(i) = -pi + real(i-1, kind=8)*dx  ! input mesh points
!  !  v(i) = sin(x(i))                   ! input data values
!  !enddo
!  !x(n) = pi
!  !v(n) = sin(x(n))
!  !dx = 2.0*pi/real(m-1, kind=8) 
!  !do i=1, m-1
!  !  xout(i) = -pi + real(i-1, kind=8)*dx  ! output points
!  !enddo
!  !xout(m)= pi
!  !
!  !d = 8                           ! target and maximum polynomial degree used for each interval
!  !interpolation_type = 2          ! 1 for DBI and 2 for PPI
!  !sten = 1                        ! optional parameter to guide stencil selection 1, 2, and 3
!  !eps0 = 0.01                     ! optional positive parameter to bound interpolant in PPI
!  !eps1 = 1.0                      ! optional positive parameter to bound interpolant in PPI
!
!  !call adaptiveInterpolation1D(x, v, n, xout, vout_apprx, m, d, interpolation_type, sten, eps0, eps1 ) 
!  !
!
!  !!-- 2D Tutorial -- !
!  !y = x                             
!  !do j=1,n
!  !  do i=1,n
!  !    v2D(i,j) = sin(x(i))*sin(y(j))  ! input data values
!  !  enddo
!  !enddo
!  !yout = xout
!  !d = 8                             ! target and maximum polynomial degree used for each interval
!  !interpolation_type = 2            ! 1 for DBI and 2 for PPI
!  !sten = 1                          ! optional parameter to guide stencil selection 1, 2, and 3
!  !eps0 = 0.01                       ! optional positive parameter to bound interpolant in PPI
!  !eps1 = 1.0                        ! optional positive parameter to bound interpolant in PPI
!
!  !call adaptiveInterpolation2D(x, y, n, n, v2D, xout, yout, m, m, vout_apprx2D, d, interpolation_type, sten, eps0, eps1 ) 
!  !
!
!  !!-- 3D Tutorial -- !
!  !z = x
!  !do k=1,n
!  !  do j=1,n
!  !    do i=1,n
!  !      v3D(i,j,k) = sin(x(i))*sin(y(j))*sin(z(k))  ! input data values
!  !    enddo
!  !  enddo
!  !enddo
!  !zout =xout
!  !d = 8                             ! target and maximum polynomial degree used for each interval
!  !interpolation_type = 2            ! 1 for DBI and 2 for PPI
!  !sten = 1                          ! optional parameter to guide stencil selection 1, 2, and 3
!  !eps0 = 0.01                       ! optional positive parameter to bound interpolant in PPI
!  !eps1 = 1.0                        ! optional positive parameter to bound interpolant in PPI
!
!  !call adaptiveInterpolation3D(x, y, z, n, n, n, v3D, xout, yout, zout, m, m, m, vout_apprx3D, &
!  !                                       d, interpolation_type, sten, eps0, eps1 ); 
! 
!
!
