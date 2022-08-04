program testing
!
!
!

  call test1()
  !call test2()
  !call test3()
  !call test4()
  call test6()
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
  
  integer(kind=8), parameter		::n = 1e+6
  integer(kind=8), parameter		::m = 1e+7
  integer(kind=8)       		::i, j
  integer       		::sten
  integer, parameter		::d =8
  real(kind=8), parameter	:: a = -1.0e-10;
  real(kind=8), parameter       :: b = 1.0e-10;
  real(kind=8), parameter       :: eps0=1.0;
  real(kind=8), parameter       :: eps1=1.0;
  real(kind=8)  	 	:: x(n), v1D(n) 
  real(kind=8)  	 	:: xout(m), v1Dout(m), v1Dout_vec(m) 
  real(kind=8) 			:: dx, v1Dout_true(m) 

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
    v1D(i) = 0.1 /(0.1 + 25*(x(i)*1.0e+10)**2);
  enddo
  do i=1, m
    v1Dout_true(i) = 0.1 /(0.1 + 25*(xout(i)*1.0e+10)**2);
  enddo
  do  j=1,3
    sten = j
    call adaptiveInterpolation1D(x, v1D, n, xout, v1Dout, m, d, 2, sten, eps0, eps1) 
    call adaptiveInterpolation1D_vec(x, v1D, n, xout, v1Dout_vec, m, d, 2, sten, eps0, eps1) 
    do i=1, m
      if(v1Dout(i) .ne. v1Dout_vec(i) .and. abs(v1Dout(i)-v1Dout_vec(i)) >1.0e-15 ) then
        write(*,*) 'FAILED:',i, abs(v1Dout_true(i)-v1Dout(i)), abs(v1Dout_true(i)-v1Dout_vec(i))
      endif
    enddo

  enddo
end subroutine


subroutine test2()
!
!  test considering small example with negative values
!

  use mod_adaptiveInterpolation

  implicit none

  integer, parameter            :: n = 20;
  integer, parameter            :: m = 1000;
  integer, parameter            :: d = 8;
  integer                       :: i, j
  integer                       :: sten
  real(kind=8), parameter       :: a = -1.0;
  real(kind=8), parameter       :: b = 1.0;
  real(kind=8)                  :: x(n), v1D(n)
  real(kind=8)                  :: xout(m), v1Dout(m), v1Dout_vec(m)
  real(kind=8)                  :: dx, pi
  real(kind=8)                  :: eps0, eps1
  
  eps0= 1;
  eps1= 1;
  pi = 4.0*atan(1.0)
  dx = (b-a) / real(n-1, kind=8)
  x(1) = a
  do i=2,n-1
    x(i) = x(i-1) + dx
    v1D(i) = sin(x(i)*pi);
  enddo
  x(n) = b
  dx = (b-a)/real(m-1, kind=8)
  xout(1) = a
  do i=2, m-1
    xout(i) =xout(i-1)+dx
  enddo
  xout(m) = b
  
  do  j=1,3
    sten = j
    call adaptiveInterpolation1D(x, v1D, n, xout, v1Dout, m, d, 2, sten, eps0, eps1) 
    call adaptiveInterpolation1D_vec(x, v1D, n, xout, v1Dout_vec, m, d, 2, sten, eps0, eps1) 
    do i=1, m
      if(v1Dout(i) .ne. v1Dout_vec(i) .and. abs(v1Dout(i)-v1Dout_vec(i)) >1.0e-10 ) then
        write(*,*) 'FAILED:',i, v1Dout(i), v1Dout_vec(i)
      endif
    enddo


  enddo
  
end subroutine 


subroutine test3()
!
!  test considering small constant fucntions
!

  use mod_adaptiveInterpolation

  implicit none

  integer, parameter            :: n = 20;
  integer, parameter            :: m = 1000;
  integer, parameter            :: d = 8;
  integer                       :: i, j
  integer                       :: sten
  real(kind=8), parameter       :: a = -1.0;
  real(kind=8), parameter       :: b = 1.0;
  real(kind=8)                  :: x(n), v1D(n)
  real(kind=8)                  :: xout(m), v1Dout(m), v1Dout_vec(m)
  real(kind=8)                  :: dx, pi
  real(kind=8)                  :: eps0, eps1
  
  eps0= 1;
  eps1= 1;
  pi = 4.0*atan(1.0)
  dx = (b-a) / real(n-1, kind=8)
  x(1) = a
  do i=2,n-1
    x(i) = x(i-1) + dx
    v1D(i) = 1.0
  enddo
  x(n) = b
  dx = (b-a)/real(m-1, kind=8)
  xout(1) = a
  do i=2, m-1
    xout(i) =xout(i-1)+dx
  enddo
  xout(m) = b
  
  do  j=1,3
    sten = j
    call adaptiveInterpolation1D(x, v1D, n, xout, v1Dout, m, d, 2, sten, eps0, eps1) 
    call adaptiveInterpolation1D_vec(x, v1D, n, xout, v1Dout_vec, m, d, 2, sten, eps0, eps1) 
    do i=1, m
      if(v1Dout(i) .ne. v1Dout_vec(i) .and. abs(v1Dout(i)-v1Dout_vec(i)) >1.0e-10 ) then
        write(*,*) 'FAILED:',i, v1Dout(i), v1Dout_vec(i)
      endif
    enddo


  enddo
 
end subroutine

subroutine test4()
!
!  test considering a linear function constant fucntions
!

  use mod_adaptiveInterpolation

  implicit none

  integer, parameter            :: n = 20;
  integer, parameter            :: m = 1000;
  integer, parameter            :: d = 8;
  integer                       :: i, j
  integer                       :: sten
  real(kind=8), parameter       :: a = -1.0;
  real(kind=8), parameter       :: b = 1.0;
  real(kind=8)                  :: x(n), v1D(n)
  real(kind=8)                  :: xout(m), v1Dout(m), v1Dout_vec(m)
  real(kind=8)                  :: dx, pi
  real(kind=8)                  :: eps0, eps1
  
  eps0= 1;
  eps1= 1;
  pi = 4.0*atan(1.0)
  dx = (b-a) / real(n-1, kind=8)
  x(1) = a
  do i=2,n-1
    x(i) = x(i-1) + dx
    v1D(i) = (-1.0/(b-a))*(x(i)+1) + 1;
  enddo
  x(n) = b
  dx = (b-a)/real(m-1, kind=8)
  xout(1) = a
  do i=2, m-1
    xout(i) =xout(i-1)+dx
  enddo
  xout(m) = b
  
  do  j=1,3
    sten = j
    call adaptiveInterpolation1D(x, v1D, n, xout, v1Dout, m, d, 2, sten, eps0, eps1) 
    call adaptiveInterpolation1D_vec(x, v1D, n, xout, v1Dout_vec, m, d, 2, sten, eps0, eps1) 
    do i=1, m
      if(v1Dout(i) .ne. v1Dout_vec(i) .and. abs(v1Dout(i)-v1Dout_vec(i)) >1.0e-10 ) then
        write(*,*) 'FAILED:',i, v1Dout(i), v1Dout_vec(i)
      endif
    enddo


  enddo
 
end subroutine 

subroutine test5()
!
!  test commapring L^{2}-norm error for st=1, 2 and 3 with Runge
!

  integer                               :: n(6), i

  do i=1, 6
    call test52(n(i))
  enddo
end subroutine

subroutine test52(n)
!
!
!
  use mod_adaptiveInterpolation
  
  implicit none
  integer, parameter                    :: m = 10000
  integer, parameter                    :: d=8
  integer, intent(in)                   :: n

  integer                               :: i, j
  integer                               :: sten
  real(kind=8), parameter               :: a = -1.0
  real(kind=8), parameter               :: b = 1.0
  real(kind=8)                          :: x(n), v1D(n)
  real(kind=8)                          :: xout(m), v1Dout(m)
  real(kind=8)                          :: v1Dout_true(m), v1Dout_vec(m)
  real(kind=8)                          :: dx, err(m)

  dx = (b-a) / real(n-1, kind=8)
  x(1) = a
  v1D(1) = 0.1/(0.1 + 25.0*x(1)*x(1));
  do i=2, n-1
    x(i) = x(i-1) + dx
    v1D(i) = 0.1/(0.1 + 25.0*x(i)*x(i));
  enddo
  x(n) = b
  v1D(n) = 0.1/(0.1 + 25.0*x(n)*x(n));

  dx = (b-a) / real(m-1, kind=8)
  xout(1) = a
  v1Dout_true(1) = 0.1/(0.1 + 25.0*xout(1)*xout(1));
  do i=2, m-1
    xout(i) = xout(i-1) + dx
    v1Dout_true(i) = 0.1/(0.1 + 25.0*xout(i)*xout(i));
  enddo
  xout(m) = b
  v1Dout_true(m) = 0.1/(0.1 + 25.0*xout(m)*xout(m));

  call adaptiveInterpolation1D(x, v1D, n, xout, v1Dout, m, d, 2, sten) 
  call adaptiveInterpolation1D_vec(x, v1D, n, xout, v1Dout_vec, m, d, 2) 
!
!
!    n = nn(i)                                               
!    x= linspace(a, b, n);
!    xout= linspace(a, b, m);
!    v1D = zeros(n,1);
!    v1Dout = zeros(m,1);
!    v1Dout_true = zeros(m,1);
!    eps0= 1;  ! set to default value
!    eps1= 0.01; ! set to default value 
!
!    d = 8;
!    fprintf('st      error \n')
!    for j=1:3
!      st = j;
!      for i=1:n
!        v1D(i) = 0.1/(0.1 + 25.0*x(i)*x(i));
!      end
!      for i=1:m
!        v1Dout_true(i) = 0.1/(0.1 + 25.0*xout(i)*xout(i));
!      end
!      v1Dout = adaptiveInterpolation1D(x, v1D, xout, d, 2, st); 
!      err = sqrt( trapz(xout, (v1Dout-v1Dout_true').^2));
!      fprintf('!d \t !.8E \n', j, err)
!    end
!  end
!  fprintf('No significant different between st=1, 2, and 3 \n')
!
end subroutine


subroutine test6()
!
!  test 2D examples with different number of points for each dimension
!

  use mod_adaptiveInterpolation

  integer, parameter            :: mx = 100;
  integer, parameter            :: my = 200;
  integer, parameter            :: nx = 17;
  integer, parameter            :: ny = 33;
  integer, parameter            :: d = 8;
  integer                       :: i, j, ii, jj
  integer                       :: sten
 
  real(kind=8), parameter       :: ax = -1.0;
  real(kind=8), parameter       :: bx = 1.0;
  real(kind=8), parameter       :: ay = -1.0;
  real(kind=8), parameter       :: by =  1.0;
  real(kind=8)                  :: x(nx), y(ny), v2D(nx, ny) 
  real(kind=8)                  :: xout(mx), yout(my)
  real(kind=8)                  :: v2Dout(mx, my),v2Dout_vec(mx, my) 
  real(kind=8), parameter       :: a = -1.0;
  real(kind=8), parameter       :: b = 1.0;
  real(kind=8)                  :: dx, dy, pi
  !real(kind=8)                  :: eps0, eps1


  dx = (bx-ax)/real(nx-1, kind=8)
  x(1) = ax
  do i=2,nx-1
    x(i) = x(i-1) + dx
  enddo
  x(nx) = bx
  !!
  dy = (by-ay)/real(nx-1, kind=8)
  y(1) = ay
  do i=2,ny-1
    y(i) = y(i-1) + dy
  enddo
  y(ny) = b
  !!
  dx = (bx-ax)/real(mx-1, kind=8)
  xout(1) = ax
  do i=2, mx-1
    xout(i) =xout(i-1)+dx
  enddo
  xout(mx) = bx
  !!
  dy = (by-ay)/real(my-1, kind=8)
  yout(1) = ay
  do i=2, my-1
    yout(i) =yout(i-1)+dy
  enddo
  yout(my) = by
  
  do j=1,ny
   do i=1,nx
     v2D(i,j) = 0.1/(0.1 + 25.0*(x(i)*x(i) + y(j)*y(j)))
   enddo
  enddo
  !
  do j=1,3
  sten = j
  call adaptiveInterpolation2D(x, y, nx, ny, v2D,  xout, yout, mx, my, v2Dout, d, 2, sten)
  call adaptiveInterpolation2D_vec(x, y, nx, ny, v2D,  xout, yout, mx, my, v2Dout_vec, d, 2, sten)

  do jj=1, my
    do ii=1,mx
      if(v2Dout(ii,jj) .ne. v2Dout_vec(ii,jj) .and. abs(v2Dout(ii,jj)-v2Dout_vec(ii,jj))> 1.0e-10) then
        write(*,*) 'FAILED: ', i, v2Dout(ii,jj), v2Dout_vec(ii,jj)
      endif
    enddo
  enddo
  enddo
end subroutine 

subroutine test7()
!
!  test 2D examples with different number of points for each dimension
!

  use mod_adaptiveInterpolation

  integer, parameter            :: mx = 100;
  integer, parameter            :: my = 200;
  integer, parameter            :: nx = 17;
  integer, parameter            :: ny = 33;
  integer, parameter            :: d = 8;
  integer                       :: i, j, ii, jj
  integer                       :: sten
 
  real(kind=8), parameter       :: ax = -1.0e-10;
  real(kind=8), parameter       :: bx = 1.0e-10;
  real(kind=8), parameter       :: ay = -1.0;
  real(kind=8), parameter       :: by =  1.0;
  real(kind=8)                  :: x(nx), y(ny), v2D(nx, ny) 
  real(kind=8)                  :: xout(mx), yout(my)
  real(kind=8)                  :: v2Dout(mx, my),v2Dout_vec(mx, my) 
  real(kind=8), parameter       :: a = -1.0;
  real(kind=8), parameter       :: b = 1.0;
  real(kind=8)                  :: dx, dy, pi
  !real(kind=8)                  :: eps0, eps1


  dx = (bx-ax)/real(nx-1, kind=8)
  x(1) = ax
  do i=2,nx-1
    x(i) = x(i-1) + dx
  enddo
  x(nx) = bx
  !!
  dy = (by-ay)/real(nx-1, kind=8)
  y(1) = ay
  do i=2,ny-1
    y(i) = y(i-1) + dy
  enddo
  y(ny) = b
  !!
  dx = (bx-ax)/real(mx-1, kind=8)
  xout(1) = ax
  do i=2, mx-1
    xout(i) =xout(i-1)+dx
  enddo
  xout(mx) = bx
  !!
  dy = (by-ay)/real(my-1, kind=8)
  yout(1) = ay
  do i=2, my-1
    yout(i) =yout(i-1)+dy
  enddo
  yout(my) = by
  
  do j=1,ny
    do i=1,nx
      v2D(i,j) = 0.1/(0.1 + 25.0*(x(i)*1.0e+10*x(i)*1.0e+10 + y(j)*y(j)));
    enddo
  enddo
  !
  do j=1,3
    sten = j
    call adaptiveInterpolation2D(x, y, nx, ny, v2D,  xout, yout, mx, my, v2Dout, d, 2, sten)
    call adaptiveInterpolation2D_vec(x, y, nx, ny, v2D,  xout, yout, mx, my, v2Dout_vec, d, 2, sten)

    do jj=1, my
      do ii=1,mx
        if(v2Dout(ii,jj) .ne. v2Dout_vec(ii,jj) .and. abs(v2Dout(ii,jj)-v2Dout_vec(ii,jj))> 1.0e-10) then
          write(*,*) 'FAILED: ', i, v2Dout(ii,jj), v2Dout_vec(ii,jj)
        endif
      enddo
    enddo
  enddo
end subroutine 


subroutine trapz(x, y, n, res)
!
! Integretion
!
  integer      :: n
  real(kind=8) :: x(n) 
  real(kind=8) :: y(n) 
  real(kind=8) :: res_even
  real(kind=8) :: res_odd
  real(kind=8) :: res

  res_even =0.0;
  res_odd = 0.0;
  do i=2,n-1
    if(mod(i,2) == 0 ) then
      res_even = res_even + y(i)
    else
      res_odd = res_odd + y(i)
    endif
    res = (x(2)-x(1))/3.0 * ( y(1) + 2.0 * res_even + 4.0* res_odd + y(n))
  enddo

end subroutine
