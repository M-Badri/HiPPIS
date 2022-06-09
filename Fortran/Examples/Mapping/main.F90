program main
!!
!! Driver to produce the mapping result based on the
!! Runge and the TWP-ICE.
!!
  implicit none
  integer 		:: nz(3)
  integer 		:: k

  !!call approximations1D()


  call approximations2D()

  !!nz = (/64, 127, 253/)
  !!do k=1,3
  !!  call mapping(nz(k))
  !!enddo
   
end program 

subroutine approximations1D()
!!
!! approximation1D is used to set up the 
!! diffferent configuration used to produces
!! the approximation results for the 1D functions
!! presented in the manuscript.
!!
  implicit none 


  integer                       :: n(5)  			!! total number points used 
  integer                       :: d(4)				!! target degree for each interpolant
  integer                       :: fun(3)			!! functions used
  integer                       :: i, ii, j, k, kk		
  integer                       :: sten(3)			!! stencil selection procedure
  integer, parameter            :: m = 10000			!! number of output points
  real(kind=8)                  :: a(3)				!! intervals left boundary
  real(kind=8)                  :: b(3)				!! intervals right boundary
  real(kind=8)                  :: eps0, eps1, eps_test(6)      !! parameters used to bound interpolants


  !!** Initialization **!!
  n = (/17, 33, 65, 129, 257/)                                              
  d = (/1, 4, 8, 16/)                                                       
  a = (/-1.0, -0.2, -1.0/)
  b = (/ 1.0,  0.2,  1.0/)
  eps_test = (/ 1.0,  0.1,  0.01, 0.001, 0.0001, 0.00 /)
  sten = (/1, 2, 3/)

  !!** modify eps0 and eps1 to change the bounds on the interpolant **!!
  eps0 = 0.01
  eps1 = 1.0

  !!** functions 1=Runge , 2= heaviside, 3=Gelb Tanner **!! 
  fun = (/1, 2, 3/)           

  !!** Used to evaluate different choices of eps0 **!!
  call testepsilon1D(sten(2), eps_test, eps1, d(3), n(1), a, b,  m)

  do ii=1, 3
    write(*,*) 'st=', sten(ii)
    do k=1, 3

      !!** Third order resulst using DBI, PPI, and PCHIP **!!
      print*, '*****  fun=', fun(k), '*****'
      do i=1, 5                                                             
        call test001(3, eps0, eps1, sten(ii), fun(k), n(i), a(k), b(k), m, 8)
      enddo

      !!** higher degree interpolation methods using DBI and PPI **!!
      do j=1, 4                                                             
        print*, '*****  d=', d(j), '*****'
        do i=1, 5                                                            
           call test001(d(j), eps0, eps1, sten(ii), fun(k), n(i), a(k), b(k), m, 8)
        enddo  
      enddo 
    enddo 
  enddo
 

end subroutine 


subroutine testepsilon1D(sten, eps0, eps1, d, n, a, b,  m)
!!
!! testepsilon1D aprroximates the Runge, smoothed Heaviside, and
!! Gelb and Tanner functions with different values of eps0 that
!! are used to bound the interpolant in the case of the PPI method. 
!!
!! INPUT
!! sten: stencil selction procedure (sten=1, sten=2, sten=3) 
!! eps0(6): array of values of eps0 
!! d:  traget polynomial degree for each interpolant
!! n: number of points
!! a(3): left boundaries
!! b(3): right boundaries
!! m: number of output points 
!!
  use mod_adaptiveInterpolation

  implicit none

  integer, intent(in)           :: n                    !! number of input points
  integer, intent(in)           :: m                    !! number of output points
  integer, intent(in)           :: d                    !! target interpolant degree
  integer, intent(in)           :: sten                 !! stencil selection procedure
  real(kind=8), intent(in)      :: a(3), b(3)           !! interval [a, b]
  real(kind=8), intent(in)      :: eps0(6), eps1        !! test values used for eps0

  integer                       :: degOut(n-1, 7)       !! degree used for each interval
  integer 		        :: i, j, k, fid
  real(kind=8)			:: x(n)                 !! uniform input mesh points  
  real(kind=8)			:: v1D(n)               !! input data values
  real(kind=8)			:: v1Dout(m, 9)         !! output values
  real(kind=8)			:: dxn, dxm             !! interval sizes 



  !!** Initialize parameters **!!
  degOut = 0
  v1Dout = 0.0
  do k=1, 3
    
    !!** calculates intreval sizes **!!
    dxn = (b(k)-a(k)) /real(n-1, kind=8)
    dxm = (b(k)-a(k)) /real(m-1, kind=8)
  

    !!** uniform mesh **!!
    do i=1,n
      x(i) = a(k) + real(i-1, kind=8)*dxn
    enddo

    !!** output mesh points **!
    do i=1,m
      v1Dout(i, 1) = a(k) + real(i-1, kind=8)*dxm
    enddo

    !!** Data values associated to input meshes **!!
    do i=1,n
      call evalFun1D(k, x(i), v1D(i), dxn)
    enddo

    !!** True solution **!!
    do i=1,m
      call evalFun1D(k, v1Dout(i, 1), v1Dout(i,2), dxn)
    enddo

    do i=1, 7
      if(i==7)then
      call adaptiveInterpolation1D(x, v1D, n, v1Dout(:,1), v1Dout(:,2+i), m, d, 1, sten, 0.0, 0.0,  degOut(1:n-1, i) ) 
      else
      call adaptiveInterpolation1D(x, v1D, n, v1Dout(:,1), v1Dout(:,2+i), m, d, 2, sten, eps0(i), eps1, degOut(1:n-1, i) ) 
      endif
    enddo

    !!** open file **!!
    fid = 10
    if( k == 1)then
      open(unit=fid, file='RungeEps', status='unknown')
    elseif( k == 2)then
      open(unit=fid, file='HeavisideEps', status='unknown')
    elseif( k == 3)then
      open(unit=fid, file='GelbTEps', status='unknown')
    endif
    !!** write to file **!!
    do i=1, m
      write(fid,'(9(3x,E30.15))') ( v1Dout(i, j), j=1, 9 )
    enddo
    !!** close file **!!
    close(fid)

    !!** open file **!!
    fid = 10
    if( k == 1)then
      open(unit=fid, file='RunegeEpsDeg', status='unknown')
    elseif( k == 2)then
      open(unit=fid, file='HeavisideEpsDeg', status='unknown')
    elseif( k == 3)then
      open(unit=fid, file='GelbTEpsDeg', status='unknown')
    endif

    !!** write to file **!!
    do i=1, n-1
      write(fid,'(7(1x, I2))') ( degOut(i, j), j=1, 7 )
    enddo
    !!** close file **!!
    close(fid)
  enddo

end subroutine testepsilon1D

subroutine test001(d, eps0, eps1, sten, fun, n, a, b, m, d_el)
!! 
!! test001 is used to approximate the Runge, smoothed Heaviside
!! and Gelbd and Tanner function using different interpolation 
!! methods.
!!
!! INPUT
!! d: maximum polynomial degree for each interval
!! eps:
!! 
  use mod_legendre
  use mod_adaptiveInterpolation

  implicit none

  integer, intent(in)           :: fun                  !! function type 
  integer, intent(in)           :: n                    !! number of input points
  integer, intent(in)           :: m                    !! number of output points
  integer, intent(in)           :: d                    !! target interpolant degree
  integer, intent(in)           :: sten
  real(kind=8), intent(in)      :: a                    !! left bounary
  real(kind=8), intent(in)      :: b                    !! right boundary
  real(kind=8), intent(in)      :: eps0			!! parameters used to bound interpolant in intervals with no hidden extrema  
  real(kind=8), intent(in)      :: eps1			!! parameters used to bound interpolant in intervals with hidden extrema  
  integer, intent(in)           :: d_el

  integer                       :: degOut(n-1, 2)      !! degree used for each interval
  integer                       :: degOut_lgl(n-1, 2)      !! degree used for each interval
  integer                       :: limiter             !!
  integer                       :: ne                   !! number of elments
  integer 		        :: i, j, k, fid, ierr, tmp_idx
  integer 		        :: is, ie, dd
  real(kind=8)			:: x(n), x_lgl(n)                     !! uniform and  LGL input mesh points  
  real(kind=8)			:: v1D(n), v1D_lgl(n)              !! input data values
  real(kind=8)			:: xout(m)                                      !! output points to be approximated 
  real(kind=8)			:: v1Dout(m), v1Dout_lgl(m)     !! approximated output values
  real(kind=8)			:: v1Dout_true(m)                               !! True values at output points
  real(kind=8)			:: dxn, dxm
  real(kind=8)			:: x_tmp(d_el+1), w_tmp(d_el+1), xl, xr
  character*16                  :: fun_name
  character*16                  :: fnumber
  character*16                  :: sst
  character*64                  :: fname

  !!** Local variables need for PCHIP **!!
  integer 		        :: nwk
  real(kind=8)			:: wk((n+1)*2), d_tmp(n+1)
  real(kind=8)			:: fdl(m)
  logical                       :: spline

  spline = .false.  !! needed for PCHIP
  nwk = (n+1)*2     !! needed for PCHIP

  write(fnumber, '("", i5.5)') d*1000+n

  !!** get function  name **!!
  if(fun ==1)then
    fun_name = "Runge"
  elseif(fun ==2)then
    fun_name = "Heaviside"
  elseif(fun ==3)then
    fun_name = "GelbT"
  else
    write(*,*) 'ERROR: Invalid fun =', fun
    write(*,*) 'Invalid function value the possible values are fun =1, 2, or 3'
    call exit(0)
  endif

  !!** get stencil selection procedure **!!
  if(sten ==1) then
    sst = "st=1"
  elseif(sten ==2) then
    sst = "st=2"
  elseif(sten ==3) then
    sst = "st=3"
  else
    write(*,*) 'ERROR: Invalid paparamter sten =', fun
    write(*,*) 'ERROR: Invalid paparamter st. The possible options are st=1, 2, or 3'
    call exit(0)
  endif
  

  !!** uniform mesh **!!
  dxn = (b-a) /real(n-1, kind=8)
  do i=1,n
    x(i) = a + real(i-1, kind=8)*dxn
  enddo

  dd = d_el
  ne = (n-1) / dd                        !! calculates the number of elements
  
  !!** LGL mesh **!!
  call legendre_gauss_lobatto(dd+1, x_tmp, w_tmp)           !! LGL nodes
  dxn = (b-a) / real(ne, kind = 8)                          !! calculates element size
  xl = a                                                    !! initialaze element left boundary 
  xr = a                                                    !! initialize element right boundary 
  is = 1
  ie = 1
  do i=1, ne
    xl = xr                                      !! update left boundary of element i
    xr = xl + dxn                                !! update right boun dary of element i
    is = ie
    ie = is + dd
    x_lgl(is:ie) = 0.50*( x_tmp* (xr-xl) + (xr+xl) )
  end do

  !!** output mesh points **!
  dxm = (b-a) /real(m-1, kind=8)
  do i=1,m
    xout(i) = a + real(i-1, kind=8)*dxm
  enddo

  !!** Data values associated to input meshes **!!
  dxn = (b-a)/real(ne, kind=8) !! dummy variable not used for calculations
  do i=1,n
    call evalFun1D(fun, x(i), v1D(i), dxn)
    call evalFun1D(fun, x_lgl(i), v1D_lgl(i), dxn)
  enddo

  !!** True solution **!!
  do i=1,m
    call evalFun1D(fun, xout(i), v1Dout_true(i), dxn)
  enddo
 
  !!if(n == 17 .or. n==33) then
  !!!!print*, '**** STD interpolation ****'
  !!call STDinterp(x, v1D, n, xout , v1Dout, m, d)
  !!call STDinterp(x_lgl, v1D_lgl, n, xout , v1Dout_lgl, m, d)

  !!!!** open file **!!
  !!fid = 10
  !!call openFile(fun, fid, 0, n, d)
  !!!!** write to file **!!
  !!do i=1, m
  !!  write(fid,'(5(3x,E30.15))') xout(i), v1Dout_true(i), v1Dout(i), v1Dout_lgl(i)
  !!enddo
  !!!!** close file **!!
  !!close(fid)
  !!endif

  !!** interpolation using PCHIP **!!
  if(d ==3 ) then
    call pchez(n, x, v1D, d_tmp, spline, wk, nwk, ierr)
    call pchev(n, x, v1D, d_tmp, m, xout, v1Dout, fdl, ierr)
 
    call pchez(n, x_lgl, v1D_lgl, d_tmp, spline, wk, nwk, ierr)
    call pchev(n, x_lgl, v1D_lgl, d_tmp, m, xout, v1Dout_lgl, fdl, ierr)

    !!** open file and write to file **!!
    fid = 10
    !!call openFile(fun, fid, d, n, d)
    fname = trim(fun_name)//trim("PCHIP")//trim(fnumber)
    open(unit=fid,file=fname, status='unknown')
    do i=1, m
      write(fid,'(4(3x,E30.15))') xout(i), v1Dout_true(i), v1Dout(i), v1Dout_lgl(i)
    enddo
    close(fid)
  endif


  !!** Interpolation using DBI **!!
  v1Dout =0.0
  call adaptiveInterpolation1D(x, v1D, n, xout, v1Dout, m, d, 1, sten, eps0, eps1, degOut(1:n-1, 1)) 
  call adaptiveInterpolation1D(x_lgl, v1D_lgl, n, xout, v1Dout_lgl, m, d, 1, sten, eps0, eps1, degOut_lgl(1:n-1, 1)) 

  !!** open file and write to file **!!
  fid = 10
  !!call openFile(fun, fid, 1, n, d)
  fname = trim(fun_name)//trim("DBI")//trim(fnumber)
  open(unit=fid,file=fname, status='unknown')
  do i=1, m
    write(fid,'(4(3x,E30.15))') xout(i), v1Dout_true(i), v1Dout(i), v1Dout_lgl(i)
  enddo
  close(fid)

  !!** Interpolation using PPI **!!
  call adaptiveInterpolation1D(x, v1D, n, xout, v1Dout, m, d, 2, sten, eps0, eps1, degOut(1:n-1, 2) ) 
  call adaptiveInterpolation1D(x_lgl, v1D_lgl, n, xout, v1Dout_lgl, m, d, 2, sten, eps0, eps1, degOut_lgl(1:n-1, 2) ) 

  !!** open file and write to file **!!
  fid = 10
  !!call openFile(fun, fid, 2, n, d)
  fname = trim(fun_name)//trim("PPI")//trim(fnumber)//trim(sst)
  open(unit=fid,file=fname, status='unknown')
  do i=1, m
    write(fid,'(4(3x,E30.15))') xout(i), v1Dout_true(i), v1Dout(i), v1Dout_lgl(i)
  enddo
  close(fid)

  !!** open file write degree used for each interpolant to file**!!
  fid = 10
  !!call openFile(fun, fid, 4, n, d)
  fname = trim(fun_name)//trim("Deg")//trim(fnumber)//trim(sst)
  open(unit=fid,file=fname, status='unknown')
  do i=1, n-1
    write(fid,'(4(1x,I2))') degOut(i, 1), degOut(i, 2), degOut_lgl(i, 1), degOut_lgl(i, 2)
  enddo
  close(fid)

end subroutine

!! 2D Examples
subroutine approximations2D()
!!
!!
!!
  implicit none

  integer                       :: nx(5), nx2(5)
  integer                       :: ny(5), ny2(5)
  integer                       :: d(4)
  integer                       :: fun(4)
  integer                       :: i, ii, j, k, kk
  integer                       :: sten(3)
  integer, parameter            :: m = 100!!1000
  real(kind=8)                  :: ax(4), bx(4)
  real(kind=8)                  :: ay(4), by(4)
  real(kind=8)                  :: eps0, eps1, eps_test(6)


  d = (/1, 4, 8, 16/)                                                       !! array with interpolants degrees
  nx = (/17, 33, 65, 129, 257/)                                                !! array with number of inputpoints
  ny = (/17, 33, 65, 129, 257/)                                                !! array with number of inputpoints
  nx2 = (/13, 25, 49, 97, 193/)                                                !! array with number of inputpoints
  ny2 = (/13, 25, 49, 97, 193/)                                                !! array with number of inputpoints

  eps_test = (/ 1.0,  0.1,  0.01, 0.001, 0.0001, 0.00 /)

  !!** set up interval x \in [ax(i), bx(i)] and y \in [ay(i), by(i)]**!! 
  ax = (/-1.0, -1.0, 0.0, -0.2 /)
  bx = (/ 1.0,  1.0, 2.0, 0.2 /)
  ay = (/-1.0, -1.0, 0.0, -0.2 /)
  by = (/ 1.0,  1.0, 1.0, 0.2 /)
  
  !!** function type 1=runge funtion , 2= heaviside, 3=Gelb Tanner **!! 
  fun = (/1, 2, 3, 4/)                                                     !! function type

  !!**
  sten = (/1, 2, 3/)
  eps0 = 0.01
  eps1 = 1.0
  call testepsilon2D(sten(2),eps0, eps1, d(3), nx(1), ny(1), ax, bx, ay, by, m)
  do ii=1, 3
    !!** comparing against PCHIP **!!
    do k=1,4 
      !!if(k .eq. 1 .or. k .eq. 4) then
        do i=1, 5
           !!** Perform interpolation and calculate different L2-error
           !    norms different methods **!!
           print*, 'nx=', nx(i), 'ny=', ny(i)
           print*, 'ax=', ax(k), 'bx=', bx(k)
           print*, 'ay=', ay(k), 'by=', by(k)
           call test002(3, eps0, eps1, sten(ii), fun(k), nx(i), ny(i), ax(k), bx(k), ay(k), by(k), m, 8)
        enddo


        print*, '*****  fun=', fun(k), '*****'
        do j=2, 4  
          print*, '*****  d=', d(j), '*****'
          do i=1,3! 5
             !!** Perform interpolation and calculate different L2-error
             !    norms different methods **!!
             print*, 'nx=', nx(i), 'ny=', ny(i)
             print*, 'ax=', ax(k), 'bx=', bx(k)
             print*, 'ay=', ay(k), 'by=', by(k)
             call test002(d(j), eps0, eps1, sten(ii), fun(k), nx(i), ny(i), ax(k), bx(k), ay(k), by(k), m, 8)
          enddo
        enddo
      !!endif
    enddo 
  enddo


end subroutine 



subroutine testepsilon2D(sten, eps0, eps1, d, nx, ny, ax, bx, ay, by, m)
!!
!!
!!

  use mod_legendre
  use mod_adaptiveInterpolation


  implicit none


  integer, intent(in)           :: nx                    !! number of input points
  integer, intent(in)           :: ny                    !! number of input points
  integer, intent(in)           :: m                    !! number of output points
  integer, intent(in)           :: d                    !! target interpolant degree
  integer, intent(in)           :: sten
  real(kind=8), intent(in)      :: ax(4), bx(4)                 !! interval [a, b]
  real(kind=8), intent(in)      :: ay(4), by(4)                 !! interval [a, b]
  real(kind=8), intent(in)      :: eps0(6), eps1

  integer 		        :: i, ii, kk, j, k, fid
  real(kind=8)			:: x(nx)                 !! input mesh points  
  real(kind=8)			:: y(ny)                 !! input mesh points  
  integer 			:: degx2(nx-1, ny)       !! input mesh points  
  integer 			:: degy2(ny-1, m)        !! input mesh points  
  real(kind=8)			:: v2D(nx, ny)           !! input data values
  real(kind=8)			:: xout(m)               !! output points to be approximated 
  real(kind=8)			:: yout(m)               !! output points to be approximated 
  real(kind=8)			:: v2Dout(m, m)          !! approximated output values
  real(kind=8)			:: v2Dout_true(m, m)       !! True values at output points
  real(kind=8)			:: v2D_tmp(m, ny)        !! True values at output points

  real(kind=8)			:: v2D_s(m*m, 10)        !! True values at output points
  real(kind=8)			:: dxn, dxm, dyn, dym
  real(kind=8)			:: h                    !! element spacing

  
  character*36                  :: filename

  do k=4, 4
    !!** calculates intreval sizes **!!
    dxn = (bx(k)-ax(k)) /real(nx-1, kind=8)
    dxm = (bx(k)-ax(k)) /real(m-1, kind=8)
    dyn = (by(k)-ay(k)) /real(ny-1, kind=8)
    dym = (by(k)-ay(k)) /real(m-1, kind=8)
  

     !!** uniform mesh **!!
     do i=1,nx
       x(i) = ax(k) + real(i-1, kind=8)*dxn
     enddo      
     do i=1,ny  
       y(i) = ay(k) + real(i-1, kind=8)*dyn
     enddo

    !!** output mesh points **!
    do i=1,m
      xout(i) = ax(k) + real(i-1, kind=8)*dxm
      yout(i) = ay(k) + real(i-1, kind=8)*dym
    enddo


    !!** only used in calculation inside of evalFun2D for fun == 4
    h = dxn

    !!** Data values associated to input meshes **!!
    do j=1, ny
      do i=1,nx
        call evalFun2D(k, x(i), y(j), v2D(i, j), h)
      enddo
    enddo

    !!** True solution **!!
    do j=1, m
      do i=1, m
        call evalFun2D(k, xout(i), yout(j), v2Dout_true(i, j), h)
      enddo
    enddo

    ii = 1
    do j=1, m
      do i=1, m
        v2D_s(ii, 1) = xout(i)
        v2D_s(ii, 2) = yout(j)
        v2D_s(ii, 3) = v2Dout_true(i, j)
        ii = ii+1
      enddo
    enddo
 
    do kk=1, 7
      !!**  Interpolation using Tensor product and PPI **!!
      v2Dout = 0.0
      v2D_tmp = 0.0
      do j=1, ny
        if(kk == 7)then
          call adaptiveInterpolation1D(x, v2D(:,j), nx, xout, v2D_tmp(:,j), m, d, 1, sten, 0.0, 0.0, degx2(:, j)) 
        else
          call adaptiveInterpolation1D(x, v2D(:,j), nx, xout, v2D_tmp(:,j), m, d, 2, sten, eps0(kk), eps1, degx2(:, j)) 
        endif
      enddo
     
      do i=1, m
        if(kk == 7)then
          call adaptiveInterpolation1D(y, v2D_tmp(i,:), ny, yout, v2Dout(i,:), m, d, 1, sten, 0.0, 0.0,  degy2(:, i)) 
        else
          call adaptiveInterpolation1D(y, v2D_tmp(i,:), ny, yout, v2Dout(i,:), m, d, 2, sten, eps0(kk), eps1, degy2(:, i)) 
        endif
      enddo
 
      ii = 1
      do j=1, m
        do i=1, m
          v2D_s(ii, kk+3) = v2Dout(i, j)
          ii = ii+1
        enddo
      enddo
      write(*,*) 'kk= ', kk, 'max =', maxval(v2D_s(:,kk+3))
    enddo

    !!** Open file **!! 
    fid = 10                                                      !! file ID
    if(k ==1 )then
      open(unit=fid, file='runge2DEps', status='unknown')
    elseif(k ==2 )then
      open(unit=fid, file='surface1Eps', status='unknown')
    elseif(k ==3 )then
      open(unit=fid, file='surface2Eps', status='unknown')
    elseif(k ==4 )then
      open(unit=fid, file='heaviside2DEps', status='unknown')
    endif
    !!** Write to open file **!!
    ii =1
    do j=1, m
      do i=1, m
        write(fid,'(10(3x,E30.15))') ( v2D_s(ii, kk), kk=1, 10 )
        ii = ii+1
      enddo
    enddo
    !!** close file **!!
    close(fid)
  enddo


end subroutine

subroutine test002(d, eps0, eps1, sten, fun, nx, ny, ax, bx, ay, by, m, d_el)
!!
!!
!!

  use mod_legendre
  use mod_adaptiveInterpolation


  implicit none


  integer, intent(in)           :: fun                  !! function type 
  integer, intent(in)           :: nx                    !! number of input points
  integer, intent(in)           :: ny                    !! number of input points
  integer, intent(in)           :: m                    !! number of output points
  integer, intent(in)           :: d                    !! target interpolant degree
  integer, intent(in)           :: d_el                    !! target interpolant degree
  integer, intent(in)           :: sten
  real(kind=8), intent(in)      :: ax, bx                 !! interval [a, b]
  real(kind=8), intent(in)      :: ay, by                 !! interval [a, b]
  real(kind=8), intent(in)      :: eps0, eps1

  integer                       :: limiter             !!
  integer 		        :: i, j, k, fid, ierr, tmp_idx
  integer 		        :: ii, jj
  integer 		        ::  nwk, seed
  integer 		        :: is, ie, nnx, nny, dd
  integer 		        :: nex, ney              !! number of elements in x and y directions respectively
  real(kind=8)			:: x(nx)                 !! input mesh points  
  real(kind=8)			:: x_lgl(nx)                 !! input mesh points  
  real(kind=8)			:: y(ny)                 !! input mesh points  
  real(kind=8)			:: y_lgl(ny)             !! input mesh points  
  integer 			:: degx(nx-1, ny)        !! input mesh points  
  integer 			:: degx_lgl(nx-1, ny)    !! input mesh points  
  integer 			:: degx2(nx-1, ny)       !! input mesh points  
  integer 			:: degx2_lgl(nx-1, ny)   !! input mesh points  
  integer 			:: degy(ny-1, m)         !! input mesh points  
  integer 			:: degy_lgl(ny-1, m)     !! input mesh points  
  integer 			:: degy2(ny-1, m)        !! input mesh points  
  integer 			:: degy2_lgl(ny-1, m)    !! input mesh points  
  real(kind=8)			:: v2D(nx, ny)           !! input data values
  real(kind=8)			:: v2D_lgl(nx, ny)          !! input data values
  real(kind=8)			:: xout(m)               !! output points to be approximated 
  real(kind=8)			:: yout(m)               !! output points to be approximated 
  real(kind=8)			:: v2Dout(m, m)          !! approximated output values
  real(kind=8)			:: v2Dout_lgl(m, m)         !! approximated output values
  real(kind=8)			:: v2Dout_true(m, m)       !! True values at output points
  real(kind=8)			:: v2D_tmp(m, ny)        !! True values at output points
  real(kind=8)			:: dxn, dxm, dyn, dym, err_L2, start_t, end_t
  real(kind=8)			:: x_tmp(d_el+1), w_tmp(d_el+1)
  real(kind=8)			:: xl, xr, yl, yr
  real(kind=8)			:: h                    !! element spacing


  !!real(kind=8)			:: kk, tmp
  real(kind=8)			:: wk((nx+1)*2), d_tmp(nx+1)
  real(kind=8)			:: wk2((ny+1)*2), d_tmp2(ny+1)
  real(kind=8)			:: fdl(m)
  logical                       :: spline

 
 
  character*16 			:: fnumber
  character*16 			:: sst
  character*16 			:: fun_name
  character*64                  :: fname

  write(fnumber, '("", i2.2, i3.3, "x", i3.3)')d, nx, ny

  !!** get function  name **!!
  if(fun ==1)then
    fun_name = "Runge2D"
  elseif(fun ==2)then
    fun_name = "T1"
  elseif(fun ==3)then
    fun_name = "T2"
  elseif(fun ==4)then
    fun_name = "Heaviside2D"
  else
    write(*,*) 'ERROR: Invalid fun =', fun
    write(*,*) 'Invalid function value the possible values are fun =1, 2, 3, or 4'
    call exit(0)
  endif

  !!** get stencil selection procedure **!!
  if(sten ==1) then
    sst = "st=1"
  elseif(sten ==2) then
    sst = "st=2"
  elseif(sten ==3) then
    sst = "st=3"
  else
    write(*,*) 'ERROR: Invalid paparamter sten =', fun
    write(*,*) 'ERROR: Invalid paparamter st. The possible options are st=1, 2, or 3'
    call exit(0)
  endif
 
  !!** Initialize variables **!!
  x = 0.0
  y = 0.0
  x_lgl = 0.0
  y_lgl = 0.0
  v2D = 0.0
  v2Dout = 0.0
  v2Dout_true = 0.0
  spline = .false.


  !!** calculates intreval sizes **!!
  dxn = (bx-ax) /real(nx-1, kind=8)
  dxm = (bx-ax) /real(m-1, kind=8)
  dyn = (by-ay) /real(ny-1, kind=8)
  dym = (by-ay) /real(m-1, kind=8)
  

  !!** unifnorm mesh **!!
  do i=1,nx
    x(i) = ax + real(i-1, kind=8)*dxn
  enddo
  do i=1,ny
    y(i) = ay + real(i-1, kind=8)*dyn
  enddo

  !!** number of elements **!!
  dd = d_el
  nex = (nx-1) / dd
  ney = (ny-1) / dd

  !!** Legendre gauss lobatto points **!!
  call legendre_gauss_lobatto(dd+1, x_tmp, w_tmp)           !! LGL mesh
  dxn = (bx-ax) / real(nex, kind = 8)               !! calculates element size
  xl = ax                                        !! initialaziation 
  xr = ax                                       !! initialization 
  is = 1
  ie = 1
  do i=1, nex
    xl = xr                                     !! left boundary of element i
    xr = xl + dxn                                !! right boun dary of element i
    is = ie
    ie = is + dd
    !!** maping from [-1,1] to [xl, xr] 
    x_lgl(is:ie) = x_tmp* (xr-xl)/2.0 + (xr+xl)/2.0
  end do
  dyn = (by-ay) / real(ney, kind = 8)               !! calculates element size
  yl = ay
  yr = ay
  is = 1
  ie = 1
  do i=1, ney
    yl = yr                                     !! left boundary of element i
    yr = yl + dyn                                !! right boun dary of element i
    is = ie
    ie = is + dd
    !!** maping from [-1,1] to [yl, yr] 
    y_lgl(is:ie) = x_tmp* (yr-yl)/2.0 + (yr+yl)/2.0
  enddo

  !!** output mesh points **!
  do i=1,m
    xout(i) = ax + real(i-1, kind=8)*dxm
    yout(i) = ay + real(i-1, kind=8)*dym
  enddo


  !!** only used in calculation inside of evalFun2D for fun == 4
  h = (bx-ax)/((nx-1)/d)

  !!** Data values associated to input meshes **!!
  do j=1, ny
    do i=1,nx
      call evalFun2D(fun, x(i), y(j), v2D(i, j), h)
      call evalFun2D(fun, x_lgl(i), y_lgl(j), v2D_lgl(i, j), h)
    enddo
  enddo

  !!** True solution **!!
  do j=1, m
    do i=1, m
      call evalFun2D(fun, xout(i), yout(j), v2Dout_true(i, j), h)
    enddo
  enddo
 
  !!**  Interpolation using Tensor product and PCHIP **!!
  if(d == 3) then 
    nwk = (nx+1)*2
    do j=1, ny
      call pchez(nx, x, v2D(:,j), d_tmp, spline, wk, nwk, ierr)
      call pchev(nx, x, v2D(:,j), d_tmp, m, xout, v2D_tmp(:, j), fdl, ierr)
    enddo
    nwk = (ny+1)*2
    do i=1, m
      call pchez(ny, y, v2D_tmp(i,:), d_tmp2, spline, wk2, nwk, ierr)
      call pchev(ny, y, v2D_tmp(i,:), d_tmp2, m, yout, v2Dout(i, :), fdl, ierr)
    enddo
    !!** interpolation using lgl **!!
    nwk = (nx+1)*2
    do j=1, ny
      call pchez(nx, x_lgl, v2D_lgl(:,j), d_tmp, spline, wk, nwk, ierr)
      call pchev(nx, x_lgl, v2D_lgl(:,j), d_tmp, m, xout, v2D_tmp(:, j), fdl, ierr)
    enddo
    nwk = (ny+1)*2
    do i=1, m
      call pchez(ny, y_lgl, v2D_tmp(i,:), d_tmp2, spline, wk2, nwk, ierr)
      call pchev(ny, y_lgl, v2D_tmp(i,:), d_tmp2, m, yout, v2Dout_lgl(i, :), fdl, ierr)
    enddo


    !!** Open file **!! 
    fid = 10                                                      !! file ID
    !!call openFile2D(fun, fid, 3, nx, ny, 3, 0)
    fname = trim(fun_name)//trim("PCHIP")//trim(fnumber)
    open(unit=fid,file=fname, status='unknown')
    !!** Write to open file **!!
    do j=1, m
      do i=1, m
        write(fid,'(5(3x,E30.15))') xout(i), yout(j), v2Dout_true(i, j), v2Dout(i, j), v2Dout_lgl(i, j)
      enddo
    enddo
    !!** close file **!!
    close(fid)
  endif
  !!**  Interpolation using Tensor product and DBI **!!
  v2Dout =0.0
  do j=1, ny
    call adaptiveInterpolation1D(x, v2D(:,j), nx, xout, v2D_tmp(:,j), m, d, 1, sten, eps0, eps1, degx(:, j)) 
  enddo
  do i=1, m
    call adaptiveInterpolation1D(y, v2D_tmp(i,:), ny, yout, v2Dout(i,:), m, d, 1, sten, eps0, eps1, degy(:, i)) 
  enddo

  !!**  Interpolation using Tensor product and DBI **!!
  v2Dout_lgl =0.0
  do j=1, ny
    call adaptiveInterpolation1D(x_lgl, v2D_lgl(:,j), nx, xout, v2D_tmp(:,j), m, d, 1, sten, eps0, eps1, degx_lgl(:, j)) 
  enddo
  do i=1, m
    call adaptiveInterpolation1D(y_lgl, v2D_tmp(i,:), ny, yout, v2Dout_lgl(i,:), m, d, 1, sten, eps0, eps1, degy_lgl(:, i)) 
  enddo
  
  !!** Open file **!! 
  fid = 10                                                      !! file ID
  !!call openFile2D(fun, fid, 1, nx, ny, d, 0)
  fname = trim(fun_name)//trim("DBI")//trim(fnumber)//trim(sst)
  open(unit=fid,file=fname, status='unknown')
  !!** Write to open file **!!
  do j=1, m
    do i=1, m
      write(fid,'(5(3x,E30.15))') xout(i), yout(j), v2Dout_true(i, j), v2Dout(i, j), v2Dout_lgl(i, j)
    enddo
  enddo
  !!** close file **!!
  close(fid)

  !!**  Interpolation using Tensor product and PPI **!!
  v2Dout = 0.0
  do j=1, ny
    call adaptiveInterpolation1D(x, v2D(:,j), nx, xout, v2D_tmp(:,j), m, d, 2, sten, eps0, eps1, degx2(:, j) ) 
  enddo

  do i=1, m
    call adaptiveInterpolation1D(y, v2D_tmp(i,:), ny, yout, v2Dout(i,:), m, d, 2, sten, eps0, eps1, degy2(:, i) ) 
  enddo
  !!**  Interpolation using Tensor product and PPI **!!
  v2Dout_lgl = 0.0
  do j=1, ny
    call adaptiveInterpolation1D(x_lgl, v2D_lgl(:,j), nx, xout, v2D_tmp(:,j), m, d, 2, sten, eps0, eps1, degx2_lgl(:, j) ) 
  enddo
  do i=1, m
    call adaptiveInterpolation1D(y_lgl, v2D_tmp(i,:), ny, yout, v2Dout_lgl(i,:), m, d, 2, sten, eps0, eps1, degy2_lgl(:, i) ) 
  enddo



  !!** Open file **!! 
  fid = 10                                                      !! file ID
  fname = trim(fun_name)//trim("PPI")//trim(fnumber)//trim(sst)
  open(unit=fid,file=fname, status='unknown')
  !!** Write to open file **!!
  do j=1, m
    do i=1, m
      write(fid,'(5(3x,E30.15))') xout(i), yout(j), v2Dout_true(i, j), v2Dout(i, j), v2Dout_lgl(i, j)
    enddo
  enddo
  !!** close file **!!
  close(fid)

end subroutine



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

!!  integer, parameter    :: niter = 1          !! number of iteration
  integer               :: nz  				!! number of point used 64 127 253
  integer               :: d(3)			!! polynomial degree used 
  integer               :: i, j , k                     !! iteration ideces
  real(kind=8)          :: zd(nz),   zp(nz)             !! uniform and LGL mesh
  real(kind=8)          :: qcp(nz),   qcp2(nz),   qc(nz),   qc2(nz)

  character*32          :: name_runge, name_qc
  character*32          :: sst
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

  do j=1,3
    do i=1, 3
      runge = runge2
      rungep = rungep2
      qc=qc2
      qcp=qcp2
      write(*, *) '********** d= ', d(i), '**********'
      call mapping2(nz, zd_runge, runge, zp_runge, rungep, d(i), j, name_runge)
      call mapping2(nz, zd, qc, zp, qcp, d(i), j, name_qc)
    enddo
  enddo


end subroutine 

subroutine mapping2(nz, zd, u, zp, u2, dd, st, profile_name)
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

  integer, intent(in)   :: nz               !! number of points
  integer, intent(in)   :: dd               !! number of points
  integer, intent(in)   :: st               !! number of points

  real(kind=8), intent(in)  :: zd(nz), zp(nz)                 !! uniform and LGL mesh
  real(kind=8), intent(in)  :: u(nz), u2(nz)
  character*12, intent(in)  :: profile_name
  character*12              :: sst 

  integer               :: fun                  !! determine which fucntion isused
  integer               :: i, j , k             !! iteration ideces
  integer               :: is, ie
  integer               :: limiter
  integer 		:: sten
  integer               :: iter
  integer               :: deg(nz-1), deg_dbi(nz-1)
  real(kind=8)          :: up(nz), ud(nz)     !! data on uniform and LGL mesh
  real(kind=8)          :: up_pchip(nz), ud_pchip(nz)     !! data on uniform and LGL mesh
  real(kind=8)          :: up_dbi(nz), ud_dbi(nz)     !! data on uniform and LGL mesh
  real(kind=8)          :: up_makima(nz), ud_makima(nz)     !! data on uniform and LGL mesh
  real(kind=8)          ::  eps0, eps1

  real(kind=8)          :: dz, xl, xr

  real(kind=8)			:: wk(nz*2), d_tmp(nz)
  real(kind=8)			:: fdl(nz)
  logical                       :: spline
  integer                       :: nwk, ierr
  integer                       :: fnumber
  character*80                  :: fname, tmp_str


  !!** To save data **!!
  real(kind=8)          :: ud_out(nz, 3), up_out(nz, 3)
  real(kind=8)          :: ud_pchip_out(nz, 3), up_pchip_out(nz, 3)
  real(kind=8)          :: ud_dbi_out(nz, 3), up_dbi_out(nz, 3)
  integer               :: deg_ud_out(nz-1, 3), deg_up_out(nz-1, 3)
  integer               :: deg_ud_dbi_out(nz-1, 3), deg_up_dbi_out(nz-1, 3)


  spline = .false.
  nwk = nz*2
  eps0 = 0.01
  eps1 = 1.00

  !!** set stencil type for file name **!!
  if(st==1) then
   sst = "st=1"
  elseif(st==2) then
   sst = "st=2"
  elseif(st==3) then
   sst = "st=3"
  endif
 
  !!** Initialize file number **!!
  if (nz .le. 1000) then
    fnumber = dd*1000 + nz
    write(*,*) 'The file number is ', fnumber
  else
    write(*,*) 'ERROR: file not set for input size nz lager than 999 '
    call exit(1)
  endif
 
  ud = u
  ud_pchip = u
  ud_dbi = u
  !!ud_makima = u

  !!** save Initial data **!!
  do i=1, nz
    ud_out(i,1) = zd(i)
    up_out(i,1) = zp(i)
    ud_pchip_out(i,1) = zd(i)
    up_pchip_out(i,1) = zp(i)
    ud_dbi_out(i,1) = zd(i)
    up_dbi_out(i,1) = zp(i)
  enddo
  do i=1, nz-1
    deg_ud_out(i,1) = 0
    deg_up_out(i,1) = 0
    deg_ud_dbi_out(i,1) = 0
    deg_up_dbi_out(i,1) = 0
  enddo

  !!** Set limiter **!!
  limiter = 2

  iter=1
  !!** save data on physics grid **!!
  do i=1, nz
   up_out(i,iter+1) = u2(i)
   up_pchip_out(i,iter+1) = u2(i)
   up_dbi_out(i,iter+1) = u2(i)
  enddo
  !!
  do i=1, nz-1
   deg_up_out(i,iter+1) = 0
   deg_up_dbi_out(i,iter+1) = 0
  enddo
   

  !!** Mapping data values from zd (dynamics mesh) to zp (physics mesh) using DBI  **!!
  call adaptiveInterpolation1D(zd, ud_dbi, nz, zp, up_dbi, nz, dd, 1, st, eps0, eps1,  deg_dbi) 
  
  !!** Mapping data values from zd (dynamics mesh) to zp (physics mesh) using PPI  **!!
  call adaptiveInterpolation1d(zd, ud, nz, zp, up, nz, dd, 2, st, eps0, eps1, deg)

  !!** Mapping data values from zd (dynamics mesh) to zp (physics mesh) using PCHIP  **!!
  call pchez(nz, zd, ud_pchip, d_tmp, spline, wk, nwk, ierr)
  call pchev(nz, zd, ud_pchip, d_tmp, nz, zp, up_pchip, fdl, ierr)


  !!** save data that is on dynamics grid**!!
  do i=1, nz
   ud_out(i,iter+1) = ud(i) 
   ud_pchip_out(i,iter+1) = ud_pchip(i)
   ud_dbi_out(i,iter+1) = ud_dbi(i)
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
  enddo


  !!** Mapping data values from zp (physics mesh) to zd (dynamics mesh)  using DBI  **!!
  call adaptiveInterpolation1D(zp, up_dbi, nz, zd(2:nz-1), ud_dbi(2:nz-1), nz-2, dd, 1, st, eps0, eps1, deg_dbi) 

  !!** Mapping data values from zp (physics mesh) to zd (dynamics mesh)  using PPI  **!!
  call adaptiveInterpolation1D(zp, up, nz, zd(2:nz-1), ud(2:nz-1), nz-2, dd, 2, st, eps0, eps1, deg)

  !!** Mapping data values from zp (physics mesh) to zd (dynamics mesh)  using PCHIP  **!!
  call pchez(nz, zp, up_pchip, d_tmp, spline, wk, nwk, ierr)
  call pchev(nz, zp, up_pchip, d_tmp, nz-2, zd(2:nz-1), ud_pchip(2:nz-1), fdl, ierr)


  iter = 2
  !!** save data on physics grid **!!
  do i=1, nz-1
    deg_up_out(i,iter+1) = deg(i)
    deg_up_dbi_out(i,iter+1) = deg_dbi(i)
  enddo
  
  !!** save data that is on dynamics grid **!!
  do i=1, nz
   ud_out(i,iter+1) = ud(i) 
   ud_pchip_out(i,iter+1) = ud_pchip(i)
   ud_dbi_out(i,iter+1) = ud_dbi(i)
  enddo
  !!
  do i=1, nz-1
   deg_ud_out(i,iter+1) = 0 
   deg_ud_dbi_out(i,iter+1) = 0
  enddo

 
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
  fname = trim(profile_name)//trim(tmp_str)//trim(sst)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(3(1x,E30.15))') (ud_out(i, j), j=1, 3)
  enddo
  close(100)


  write(tmp_str, '("dDEGPPI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)//trim(sst)
  open(100,file=fname, status='unknown')
  do i=1, nz-1
    write(100, '(1002(1x,I2))') (deg_ud_out(i, j), j=1, 3)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname

  write(tmp_str, '("pPPI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)//trim(sst)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.15))') (up_out(i, j), j=1, 3)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname


  write(tmp_str, '("pDEGPPI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)//trim(sst)
  open(100,file=fname, status='unknown')
  do i=1, nz-1
    write(100, '(1002(1x,I2))') (deg_up_out(i, j), j=1, 3)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname

  write(tmp_str, '("dDBI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)//trim(sst)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.15))') (ud_dbi_out(i, j), j=1, 3)
  enddo
  close(100)

  write(tmp_str, '("dDEGDBI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)//trim(sst)
  open(100,file=fname, status='unknown')
  do i=1, nz-1
    write(100, '(1002(1x,I2))') (deg_ud_dbi_out(i, j), j=1, 3)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname

  write(tmp_str, '("pDBI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)//trim(sst)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.15))') (up_dbi_out(i, j), j=1, 3)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname

  write(tmp_str, '("pDEGDBI", i5.5)')fnumber
  fname = trim(profile_name)//trim(tmp_str)//trim(sst)
  open(102,file=fname, status='unknown')
  do i=1, nz-1
    write(102, '(1002(1x,I2))') (deg_up_dbi_out(i, j), j=1, 3)
  enddo
  close(102)
  write(*,*) 'Saved in file name ', fname

  write(tmp_str, '("dPCHIP", i5.5)') 3*1000+nz
  fname = trim(profile_name)//trim(tmp_str)//trim(sst)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.15))') (ud_pchip_out(i, j), j=1, 3)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname


  write(tmp_str, '("pPCHIP", i5.5)') 3*1000+nz
  fname = trim(profile_name)//trim(tmp_str)//trim(sst)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.15))') (up_pchip_out(i, j), j=1, 3)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname


end subroutine


subroutine evalFun1D(fun, x, v, h)
!!
!! To evaluate different function. fun determine
!! the function that will be evaluated at the points x
!!

  implicit none 

  integer, intent(in)                    :: fun                  !! function type
  real(kind=8), intent(in)               :: x                    !! point
  real(kind=8), intent(in          )     :: h                    !! elements size
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

subroutine evalFun2D(fun, x, y, v, h)
!!
!! To evaluate different fuunction. fun determine
!! the function that will be evaluated at the points x
!!
  implicit none 

  integer, intent(in)           :: fun                  !! function type
  real(kind=8), intent(in)      :: x                    !! point
  real(kind=8), intent(in)      :: y                    !! point
  real(kind=8), intent(in)      :: h                    !! element spacing
  real(kind=8), intent(out)     :: v                    !! point
  real(kind=8)                  :: pi, k, delta                !! temporary variables
  
  !!** intialize variables **!!
  k = 100
  pi = 4.0*atan(1.0) 
  delta = 0.01
  

  !!** 1D runge function **!!
  if(fun .eq. 1) then
    v = 1.0 / ( 1.0 + 25.0 * ( x*x + y*y) )
  !!** **!!
  else if(fun .eq. 2)then
    if( (x-1.5)*(x-1.5) + (y-0.5)*(y-0.5) .le. 1.0/16.0 )then
      v = 2.0*0.5*( cos(8.0*atan(1.0)*sqrt((x-1.5)*(x-1.5) &
          + (y-0.5)*(y-0.5))))!!+1)
    else if(y-x .ge. 0.5)then
      v = 1.0
    else if(0.0 .le. y-x .and. y-x .le. 0.5)then
      v = 2.0*(y-x)
    else
      v = 0.0
    endif

  !!** **!!
  else if(fun .eq. 3)then
    v = max(0.0, sin(4.0*atan(1.0)*x)*sin(4.0*atan(1.0)*y) ) 

  !!** Smoothed Heaviside function
  else if(fun .eq. 4)then
    v = 1.0 / ( 1.0 + exp(-2*k*(x+y)*sqrt(2.0)/2.0) )
  endif

end subroutine evalFun2D

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


