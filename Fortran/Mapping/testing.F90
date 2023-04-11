program testing
!
!
!

  call test1()
  call test2()
  call test3()
  call test4()
  call test5()
  call test6()
  call test7()
end program
!---------------------------------------------------------------------------------------------!
! test1
!---------------------------------------------------------------------------------------------!

subroutine test1()
!
!  test considering small intervals for input and output mesh
!  This test verifies that small values used for the input 
!  mesh points does not cause large oscillations or failure
!  of the DBI and PPI methods.
!

  use mod_adaptiveInterpolation

  implicit none
  
  integer, parameter            ::n = 1e+4
  integer, parameter            ::m = 1e+5
  integer                       ::i, j
  integer                       ::sten
  integer, parameter            ::d =8
  real(kind=8), parameter       :: a = -1.0e-10;
  real(kind=8), parameter       :: b = 1.0e-10;
  real(kind=8), parameter       :: eps0=1.0;
  real(kind=8), parameter       :: eps1=1.0;
  real(kind=8)                  :: x(n), v1D(n) 
  real(kind=8)                  :: xout(m), v1Dout(m), v1Dout_vec(m) 
  real(kind=8)                  :: dx, v1Dout_true(m) 
  logical                       :: check


  check = .true.
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
    v1Dout = 0.0
    v1Dout_vec = 0.0
    call adaptiveInterpolation1D(x, v1D, n, xout, v1Dout, m, d, 2, sten, eps0, eps1) 
    call adaptiveInterpolation1D_vec(x, v1D, n, xout, v1Dout_vec, m, d, 2, sten, eps0, eps1) 
    do i=1, m
      if( ( v1Dout(i) .ne. v1Dout_vec(i) .and. abs(v1Dout(i)-v1Dout_vec(i)) >1.0e-10 ) .or. &
         abs(v1Dout(i)-v1Dout_true(i))>1.0e-1 .or. abs(v1Dout(i)-v1Dout_true(i))>1.0e-1 ) then
        write(*,*) 'test1() FAILED: the difference betewn approximated and true solution at ',i
        write(*,*) abs(v1Dout_true(i)-v1Dout(i)), abs(v1Dout_true(i)-v1Dout_vec(i))
        write(*,*) v1Dout(i), v1Dout_vec(i), v1Dout_vec(i)
        check = .false.
      endif
    enddo
  enddo

  if(check .eqv. .true.) then
    write(*,*) '--- test01() Passed ---'
  endif
end subroutine


subroutine test2()
!
!  test considering example with negative values.
!  This verifies that the DBI and PPI can approximate
!  functions that are not positive.
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
  real(kind=8)                  :: xout(m), v1Dout(m)
  real(kind=8)                  :: v1Dout_vec(m), v1Dout_true(m)
  real(kind=8)                  :: dx, pi
  logical                       :: check


  check = .true.
  pi = 4.0*atan(1.0)
  dx = (b-a) / real(n-1, kind=8)
  x(1) = a
  do i=2,n-1
    x(i) = x(i-1) + dx
  enddo
  x(n) = b
  !!
  do i=1,n
    v1D(i) = sin(x(i)*pi);
  enddo
  dx = (b-a)/real(m-1, kind=8)
  xout(1) = a
  do i=2, m-1
    xout(i) =xout(i-1)+dx
  enddo
  xout(m) = b
  do i=1, m
   v1Dout_true(i)  = sin(xout(i)*pi);
  enddo
  
  do  j=1,3
    sten = j
    call adaptiveInterpolation1D(x, v1D, n, xout, v1Dout, m, d, 2, sten) 
    call adaptiveInterpolation1D_vec(x, v1D, n, xout, v1Dout_vec, m, d, 2, sten) 
    do i=1, m
      if( ( v1Dout(i) .ne. v1Dout_vec(i) .and. abs(v1Dout(i)-v1Dout_vec(i)) >1.0e-10 ) .or. &
         abs(v1Dout(i)-v1Dout_true(i))>1.0e-1 .or. abs(v1Dout(i)-v1Dout_true(i))>1.0e-1 ) then
        write(*,*) 'test2() FAILED: the difference betewn approximated and true solution at ',i
        write(*,*) abs(v1Dout_true(i)-v1Dout(i)), abs(v1Dout_true(i)-v1Dout_vec(i))
        write(*,*) v1Dout(i), v1Dout_vec(i), v1Dout_vec(i)
        check = .false.
      endif
    enddo
  enddo
  
  if(check .eqv. .true.) then
    write(*,*) '--- test02() Passed ---'
  endif

end subroutine 


subroutine test3()
!
!  test considering a constant fucntion.
!  This test verifies that the DBI and PPI methods
!  are able to approximate linear functions
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
  real(kind=8)                  :: xout(m), v1Dout(m)
  real(kind=8)                  :: v1Dout_vec(m), v1Dout_true(m)
  real(kind=8)                  :: dx, pi
  logical                       :: check
  
  check = .true.
  pi = 4.0*atan(1.0)
  dx = (b-a) / real(n-1, kind=8)
  x(1) = a
  do i=2,n-1
    x(i) = x(i-1) + dx
  enddo
  x(n) = b

  do i=1,n
    v1D(i) = 1.0
  enddo
  dx = (b-a)/real(m-1, kind=8)
  xout(1) = a
  do i=2, m-1
    xout(i) =xout(i-1)+dx
  enddo
  xout(m) = b
  do i=1, m
    v1Dout_true(i) = 1.0
  enddo
  do  j=1,3
    sten = j
    call adaptiveInterpolation1D(x, v1D, n, xout, v1Dout, m, d, 2, sten) 
    call adaptiveInterpolation1D_vec(x, v1D, n, xout, v1Dout_vec, m, d, 2, sten) 
    do i=1, m
      if( (v1Dout(i) .ne. v1Dout_vec(i) .and. abs(v1Dout(i)-v1Dout_vec(i)) >1.0e-10 ) .or. &
         abs(v1Dout(i)-v1Dout_true(i))>1.0e-10 .or. abs(v1Dout(i)-v1Dout_true(i))>1.0e-10 ) then
        write(*,*) 'test3() FAILED: the difference betewn approximated and true solution at ',i
        write(*,*) abs(v1Dout_true(i)-v1Dout(i)), abs(v1Dout_true(i)-v1Dout_vec(i))
        write(*,*) v1Dout(i), v1Dout_vec(i), v1Dout_vec(i)
        check = .false.
      endif
    enddo
  enddo
 
  if(check .eqv. .true.) then
    write(*,*) '--- test03() Passed ---'
  endif
end subroutine

subroutine test4()
!
!  test considering a constant fucntion.
!  This test verifies that the DBI and PPI methods
!  are able to approximate linear functions
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
  real(kind=8)                  :: xout(m), v1Dout(m)
  real(kind=8)                  :: v1Dout_vec(m), v1Dout_true(m)
  real(kind=8)                  :: dx, pi
  logical                       :: check
   
  check = .true. 
  pi = 4.0*atan(1.0)
  dx = (b-a) / real(n-1, kind=8)
  x(1) = a
  do i=2,n-1
    x(i) = x(i-1) + dx
  enddo
  do i=1,n
    v1D(i) = (-1.0/(b-a))*(x(i)+1) + 1;
  enddo
  x(n) = b
  dx = (b-a)/real(m-1, kind=8)
  xout(1) = a
  do i=2, m-1
    xout(i) =xout(i-1)+dx
  enddo
  xout(m) = b
  do i=1, m
    v1Dout_true(i) = (-1.0/(b-a))*(xout(i)+1) + 1;
  enddo
  
  do  j=1,3
    sten = j
    call adaptiveInterpolation1D(x, v1D, n, xout, v1Dout, m, d, 2, sten) 
    call adaptiveInterpolation1D_vec(x, v1D, n, xout, v1Dout_vec, m, d, 2, sten) 
    do i=1, m
      if( (v1Dout(i) .ne. v1Dout_vec(i) .and. abs(v1Dout(i)-v1Dout_vec(i)) >1.0e-10) .or. &
         abs(v1Dout(i)-v1Dout_true(i))>1.0e-10 .or. abs(v1Dout(i)-v1Dout_true(i))>1.0e-10 ) then
        write(*,*) 'test4() FAILED: the difference betewn approximated and true solution at ',i
        write(*,*) abs(v1Dout_true(i)-v1Dout(i)), abs(v1Dout_true(i)-v1Dout_vec(i))
        write(*,*) v1Dout(i), v1Dout_vec(i), v1Dout_vec(i)
        check = .false.
      endif
    enddo
  enddo
 
  if(check .eqv. .true.) then
    write(*,*) '--- test04() Passed ---'
  endif
end subroutine 

subroutine test5()
!
!  test commapring L^{2}-norm error for st=1, 2 and 3 with Runge
!

  integer                               :: n(5), i
  logical                               :: check

  n = (/17, 33, 65, 129, 257/)
  check = .true.
  do i=1,5 
    write(*,*) '---', n(i), ' ---'
    write(*,*) 'st      error unvectorized  error vecotorized ' 
    call test52(n(i), check)
  enddo

  if(check .eqv. .true.) then
    write(*,*) 'No significant different between st=1, 2, and 3'
    write(*,*) '--- test05() Passed ---'
  endif
end subroutine

subroutine test52(n, check)
!
! 
!
  use mod_adaptiveInterpolation
  
  implicit none
  integer, parameter                    :: m = 10000
  integer, parameter                    :: d=8
  integer, intent(in)                   :: n
  logical, intent(inout)                :: check

  integer                               :: i, j
  integer                               :: sten
  real(kind=8), parameter               :: a = -1.0
  real(kind=8), parameter               :: b = 1.0
  real(kind=8)                          :: x(n), v1D(n)
  real(kind=8)                          :: xout(m), v1Dout(m)
  real(kind=8)                          :: v1Dout_true(m), v1Dout_vec(m)
  real(kind=8)                          :: dx, err_v(m), err_v2(m)
  real(kind=8)                          :: tmp(3), tmp2(3)

  dx = (b-a) / real(n-1, kind=8)
  x(1) = a
  do i=2, n-1
    x(i) = x(i-1) + dx
  enddo
  x(n) = b
  do i=1, n
    v1D(i) = 0.1/(0.1 + 25.0*x(i)*x(i));
  enddo
  dx = (b-a) / real(m-1, kind=8)
  xout(1) = a
  do i=2, m-1
    xout(i) = xout(i-1) + dx
  enddo
  xout(m) = b
  do i=1, m
    v1Dout_true(i) = 0.1/(0.1 + 25.0*xout(i)*xout(i));
  enddo

  tmp = 20
  tmp2 = 20
  do j=1,3
    sten = j
    call adaptiveInterpolation1D(x, v1D, n, xout, v1Dout, m, d, 2, sten) 
    call adaptiveInterpolation1D_vec(x, v1D, n, xout, v1Dout_vec, m, d, 2, sten) 
    do i=1, m
      err_v(i) = abs(v1Dout(i)-v1Dout_true(i))**2
      err_v2(i) = abs(v1Dout_vec(i)-v1Dout_true(i))**2
    enddo
    
    call trapz(xout, err_v, m, tmp(j))
    call trapz(xout, err_v2, m, tmp2(j))
    tmp(j) = sqrt(tmp(j)) 
    tmp2(j) = sqrt(tmp2(j)) 
    write(*,*) j, tmp(j), tmp2(j)
  enddo

  if( abs(tmp(1)-tmp(2)) > 10 .or. abs(tmp(2)-tmp(3)) > 10 .or. abs(tmp(1)-tmp(3)) > 10 .or. & 
      abs(tmp2(1)-tmp2(2)) > 10 .or. abs(tmp2(2)-tmp2(3)) > 10 .or. abs(tmp2(1)-tmp2(3)) > 10) then
    check = .false.
    write(*,*) 'test5() FAILED'
  endif
end subroutine


subroutine test6()
!
!  test 2D examples with different number of points for each dimension.
!  This test verifies that the DBI and PPI works when for 2D examples 
!  with different number points in each dirrection i.e. nx not equal to ny 
!  and mx not eaqual to my 
!  

  use mod_adaptiveInterpolation

  integer, parameter            :: mx = 100;
  integer, parameter            :: my = 200;
  integer, parameter            :: nx = 33;
  integer, parameter            :: ny = 65;
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
  real(kind=8)                  :: v2Dout_true(mx, my) 
  real(kind=8)                  :: dx, dy, pi
  logical                       :: check

  check = .true.
  dx = (bx-ax)/real(nx-1, kind=8)
  x(1) = ax
  do i=2,nx-1
    x(i) = x(i-1) + dx
  enddo
  x(nx) = bx
  !!
  dy = (by-ay)/real(ny-1, kind=8)
  y(1) = ay
  do i=2,ny-1
    y(i) = y(i-1) + dy
  enddo
  y(ny) = by
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
  !!
  do j=1,my
    do i=1,mx
      v2Dout_true(i,j) = 0.1/(0.1 + 25.0*(xout(i)*xout(i) + yout(j)*yout(j)))
    enddo
  enddo
  !!
  do j=1,3
    sten = j
    call adaptiveInterpolation2D(x, y, nx, ny, v2D,  xout, yout, mx, my, v2Dout, d, 2, sten)
    call adaptiveInterpolation2D_vec(x, y, nx, ny, v2D,  xout, yout, mx, my, v2Dout_vec, d, 2, sten)
    do jj=1, my
      do ii=1, mx
        if( (v2Dout(ii,jj) .ne. v2Dout_vec(ii,jj) .and. abs(v2Dout(ii,jj)-v2Dout_vec(ii,jj)) >1.0e-10) .or. &
           abs(v2Dout(ii,jj)-v2Dout_true(ii,jj))>1.0e-1 .or. abs(v2Dout(ii,jj)-v2Dout_true(ii,jj))>1.0e-1 ) then
          write(*,*) 'test4() FAILED: the difference betewn approximated and true solution at ',ii,jj
          write(*,*) abs(v2Dout_true(ii,jj)-v2Dout(ii,jj)), abs(v2Dout_true(ii,jj)-v2Dout_vec(ii,jj))
          write(*,*) v2Dout(ii,jj), v2Dout_vec(ii,jj), v2Dout_vec(ii,jj)
          check = .false.
        endif
      enddo
    enddo
  enddo

  if(check .eqv. .true.) then
    write(*,*) '--- test06() Passed ---'
  endif
end subroutine 


subroutine test7()
!
!  test 2D examples with different number of points for each dimension.
!  with small input mesh values. This test verifies that the DBI and PPI 
!  works when for 2D examples with different number points in each dirrection 
!  i.e. nx not equal to ny and mx not eaqual to my 
!  

  use mod_adaptiveInterpolation

  integer, parameter            :: mx = 100;
  integer, parameter            :: my = 200;
  integer, parameter            :: nx = 33;
  integer, parameter            :: ny = 65;
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
  real(kind=8)                  :: v2Dout_true(mx, my) 
  real(kind=8)                  :: dx, dy, pi
  logical                       :: check

  check = .true.
  dx = (bx-ax)/real(nx-1, kind=8)
  x(1) = ax
  do i=2,nx-1
    x(i) = x(i-1) + dx
  enddo
  x(nx) = bx
  !!
  dy = (by-ay)/real(ny-1, kind=8)
  y(1) = ay
  do i=2,ny-1
    y(i) = y(i-1) + dy
  enddo
  y(ny) = by
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
  !!
  do j=1,my
    do i=1,mx
      v2Dout_true(i,j) = 0.1/(0.1 + 25.0*(xout(i)*1.0e+10*xout(i)*1.0e+10 + yout(j)*yout(j)));
    enddo
  enddo
  !!
  do j=1,3
    sten = j
    call adaptiveInterpolation2D(x, y, nx, ny, v2D,  xout, yout, mx, my, v2Dout, d, 2, sten)
    call adaptiveInterpolation2D_vec(x, y, nx, ny, v2D,  xout, yout, mx, my, v2Dout_vec, d, 2, sten)
    do jj=1, my
      do ii=1, mx
        if( (v2Dout(ii,jj) .ne. v2Dout_vec(ii,jj) .and. abs(v2Dout(ii,jj)-v2Dout_vec(ii,jj)) >1.0e-10) .or. &
           abs(v2Dout(ii,jj)-v2Dout_true(ii,jj))>1.0e-1 .or. abs(v2Dout(ii,jj)-v2Dout_true(ii,jj))>1.0e-1 ) then
          write(*,*) 'test4() FAILED: the difference betewn approximated and true solution at ',ii,jj
          write(*,*) abs(v2Dout_true(ii,jj)-v2Dout(ii,jj)), abs(v2Dout_true(ii,jj)-v2Dout_vec(ii,jj))
          write(*,*) v2Dout(ii,jj), v2Dout_vec(ii,jj), v2Dout_vec(ii,jj)
          check = .false.
        endif
      enddo
    enddo
  enddo

  if(check .eqv. .true.) then
    write(*,*) '--- test07() Passed ---'
  endif
end subroutine 

subroutine trapz(x, y, n, res)
!
! Integretion
!
  integer      :: n, i
  real(kind=8) :: x(n) 
  real(kind=8) :: y(n) 
  real(kind=8) :: res

  res = 0.0;
  do i=2, n
    res = res + (y(i)+y(i-1))/2.0 * (x(i)-x(i-1))
  enddo
end subroutine
