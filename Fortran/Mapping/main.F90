program main
!!
!! Driver use to produce the approximation and mapping result based on the
!! different function and the TWP-ICE example.
!!
  implicit none
  integer 		:: nz(3)
  integer 		:: k

  call approximations1D()


  nz = (/64, 127, 253/)
  do k=1,3
    call mapping(nz(k))
  enddo
   
  call approximations2D()
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
  integer                       :: d(3)				!! target degree for each interpolant
  integer                       :: fun(3)			!! functions used
  integer                       :: i, ii, j, k, kk		
  integer                       :: sten 			!! stencil selection procedure
  integer, parameter            :: m = 10000			!! number of output points
  real(kind=8)                  :: a(3)				!! intervals left boundary
  real(kind=8)                  :: b(3)				!! intervals right boundary
  real(kind=8)                  :: eps0, eps1, eps_test(6)      !! parameters used to bound interpolants


  !!** Initialization **!!
  n = (/17, 33, 65, 129, 257/)                                              
  d = (/3, 4, 8/)                                                       
  a = (/-1.0, -0.2, -1.0/)
  b = (/ 1.0,  0.2,  1.0/)
  eps_test = (/ 1.0,  0.1,  0.01, 0.001, 0.0001, 0.00 /)
  sten = 3

  !!** modify eps0 and eps1 to change the bounds on the interpolant **!!
  eps0 = 0.01
  eps1 = 1.0

  !!** functions 1=Runge , 2= heaviside, 3=Gelb Tanner **!! 
  fun = (/1, 2, 3/)           

  !!** Used to evaluate different choices of eps0 **!!
  call testepsilon1D(2, eps_test, eps1, d(3), n(1), a, b,  m)

  do k=1,3

    !!** higher degree interpolation methods using DBI and PPI **!!
    do j=1, 3                                                             
      print*, '*****  d=', d(j), '*****'
      do i=1, 5                                                            
         call test001(d(j), eps0, eps1, sten, fun(k), n(i), a(k), b(k), m, 8)
      enddo  
    enddo 
  enddo 
 

end subroutine 


subroutine testepsilon1D(sten, eps0, eps1, d, n, a, b,  m)
!!
!! testepsilon1D aprroximates the Runge, smoothed Heaviside, and
!! Gelb and Tanner functions with different values of eps0 that
!! are used to bound the interpolant in the case of the PPI method. 
!! This function produces the results used to build the 1D figures 
!! In the manuscript.
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

  integer 		        :: i, j, k, fid
  real(kind=8)			:: x(n)                 !! uniform input mesh points  
  real(kind=8)			:: v1D(n)               !! input data values
  real(kind=8)			:: v1Dout(m, 9)         !! output values
  real(kind=8)			:: dxn, dxm             !! interval sizes 


  !!** Initialize parameters **!!
  do k=1, 3
    
    !!** calculates intreval sizes **!!
    dxn = (b(k)-a(k)) /real(n-1, kind=8)
    dxm = (b(k)-a(k)) /real(m-1, kind=8)
  

    !!** uniform mesh **!!
    do i=1,n-1
      x(i) = a(k) + real(i-1, kind=8)*dxn
    enddo
    x(n) = b(k)

    !!** output mesh points **!
    do i=1,m-1
      v1Dout(i, 1) = a(k) + real(i-1, kind=8)*dxm
    enddo
    v1Dout(m, 1) = b(k)
    

    !!** Data values associated to input meshes **!!
    do i=1,n
      call evalFun1D(k, x(i), v1D(i))
    enddo

    !!** True solution **!!
    do i=1,m
      call evalFun1D(k, v1Dout(i, 1), v1Dout(i,2))
    enddo

    do i=1, 7
      if(i==7)then
      call adaptiveInterpolation1D(x, v1D, n, v1Dout(:,1), v1Dout(:,2+i), m, d, 1, sten, eps1, eps1 ) 
      else
      call adaptiveInterpolation1D(x, v1D, n, v1Dout(:,1), v1Dout(:,2+i), m, d, 2, sten, eps0(i), eps1) 
      endif
    enddo

    !!** open file **!!
    fid = 10
    if( k == 1)then
      open(unit=fid, file='mapping_data/data/RungeEps', status='unknown')
    elseif( k == 2)then
      open(unit=fid, file='mapping_data/data/HeavisideEps', status='unknown')
    elseif( k == 3)then
      open(unit=fid, file='mapping_data/data/GelbTEps', status='unknown')
    endif
    !!** write to file **!!
    do i=1, m
      write(fid,'(9(3x,E30.15))') ( v1Dout(i, j), j=1, 9 )
    enddo
    !!** close file **!!
    close(fid)
  enddo

end subroutine testepsilon1D

subroutine test001(d, eps0, eps1, sten, fun, n, a, b, m, d_el)
!! 
!! test001 is used to approximate the Runge, smoothed Heaviside
!! and Gelbd and Tanner function using different interpolation 
!! methods. This function is used toproduce the 1D results presented
!! in the manuscript.
!!
!! INPUT
!! d: maximum polynomial degree for each interval
!! eps0: positive user-supplied value used to bound interpolant for 
!!       intervalswith no extrema.
!! eps1: positive user-supplied value used to bound interpolant for 
!!       intervals with extrema.
!! sten: user-supplied value used to indicate stencil selection process
!!       possible choices are sten=1, sten=2, sten=3.
!! fun: used to indicate function used
!! n: number of input points
!! m: number of output points
!! a: global interval left boundary
!! b: global right interval boundary
!! d_el: number of LGL points in each element
 
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

  integer                       :: ne                   !! number of elments
  integer 		        :: i, j, k, fid, ierr, tmp_idx
  integer 		        :: is, ie, dd
  real(kind=8)			:: x(n)                      !! uniform and  LGL input mesh points  
  real(kind=8)			:: v1D(n)                    !! input data values
  real(kind=8)			:: xout(m)                   !! output points to be approximated 
  real(kind=8)			:: v1Dout(m)                 !! approximated output values
  real(kind=8)			:: v1Dout_true(m)            !! True values at output points
  real(kind=8)			:: dxn, dxm
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
  do i=1,n-1
    x(i) = a + real(i-1, kind=8)*dxn
  enddo
  x(n)= b

  dd = d_el
  ne = (n-1) / dd                        !! calculates the number of elements
  
  !!** output mesh points **!
  dxm = (b-a) /real(m-1, kind=8)
  do i=1,m-1
    xout(i) = a + real(i-1, kind=8)*dxm
  enddo
  xout(m) = b

  !!** Data values associated to input meshes **!!
  do i=1,n
    call evalFun1D(fun, x(i), v1D(i))
  enddo

  !!** True solution **!!
  do i=1,m
    call evalFun1D(fun, xout(i), v1Dout_true(i))
  enddo
 
  !!** interpolation using PCHIP **!!
  if(d ==3 ) then
    call pchez(n, x, v1D, d_tmp, spline, wk, nwk, ierr)
    call pchev(n, x, v1D, d_tmp, m, xout, v1Dout, fdl, ierr)
 
    !!** open file and write to file **!!
    fid = 10
    fname =trim("mapping_data/data/")//trim(fun_name)//trim("PCHIP")//trim(fnumber)
    open(unit=fid,file=fname, status='unknown')
    do i=1, m
      write(fid,'(3(3x,E30.15))') xout(i), v1Dout_true(i), v1Dout(i)
    enddo
    close(fid)
  endif


  !!** Interpolation using DBI **!!
  v1Dout =0.0
  call adaptiveInterpolation1D(x, v1D, n, xout, v1Dout, m, d, 1, sten, eps0, eps1) 

  !!** open file and write to file **!!
  fid = 10
  fname = trim("mapping_data/data/")//trim(fun_name)//trim("DBI")//trim(fnumber)//trim(sst)
  open(unit=fid,file=fname, status='unknown')
  do i=1, m
    write(fid,'(3(3x,E30.15))') xout(i), v1Dout_true(i), v1Dout(i)
  enddo
  close(fid)

  !!** Interpolation using PPI **!!
  call adaptiveInterpolation1D(x, v1D, n, xout, v1Dout, m, d, 2, sten, eps0, eps1) 

  !!** open file and write to file **!!
  fid = 10
  fname = trim("mapping_data/data/")//trim(fun_name)//trim("PPI")//trim(fnumber)//trim(sst)
  open(unit=fid,file=fname, status='unknown')
  do i=1, m
    write(fid,'(3(3x,E30.15))') xout(i), v1Dout_true(i), v1Dout(i)
  enddo
  close(fid)

end subroutine

subroutine approximations2D()
!!
!! approximation2D is used to set up the 
!! diffferent configurations used to produces
!! the approximation results for the 2D functions
!! presented in the manuscript.
!!
  implicit none

  integer                       :: nx(5), nx2(5)
  integer                       :: ny(5), ny2(5)
  integer                       :: d(3)
  integer                       :: fun(3)
  integer                       :: i, ii, j, k, kk
  integer                       :: sten 
  integer, parameter            :: m = 1000   !!CHANGE m TO SMALLER VALUE FOR LESS RUNTIME
  real(kind=8)                  :: ax(3), bx(3)
  real(kind=8)                  :: ay(3), by(3)
  real(kind=8)                  :: eps0, eps1, eps_test(6)


  d = (/3, 4, 8/)                                                       !! array with interpolants degrees
  nx = (/17, 33, 65, 129, 257/)                                                !! array with number of inputpoints
  ny = (/17, 33, 65, 129, 257/)                                                !! array with number of inputpoints

  eps_test = (/ 1.0,  0.1,  0.01, 0.001, 0.0001, 0.00 /)

  !!** set up interval x \in [ax(i), bx(i)] and y \in [ay(i), by(i)]**!! 
  ax = (/-1.00, -0.20, 0.00 /)
  bx = (/ 1.00,  0.20, 2.00 /)
  ay = (/-1.00, -0.20, 0.00 /)
  by = (/ 1.00,  0.20, 1.00 /)
  
  !!** function type 1=runge funtion , 2= heaviside, 3=Gelb Tanner **!! 
  fun = (/1, 2, 3/)                                                     !! function type

  !!**
  sten = 3
  eps0 = 0.01
  eps1 = 1.0
  call testepsilon2D(2,eps_test, eps1, d(3), nx(1), ny(1), ax, bx, ay, by, m)
  !!** comparing against PCHIP **!!
  do k=1,3 
      print*, '*****  fun=', fun(k), '*****'
      do j=1, 3  
        print*, '*****  d=', d(j), '*****'
        do i=1, 5
           !!** Perform interpolation and calculate different L2-error
           !    norms different methods **!!
           print*, 'nx=', nx(i), 'ny=', ny(i)
           print*, 'ax=', ax(k), 'bx=', bx(k)
           print*, 'ay=', ay(k), 'by=', by(k)
           call test002(d(j), eps0, eps1, sten, fun(k), nx(i), ny(i), ax(k), bx(k), ay(k), by(k), m, 8)
        enddo
      enddo
  enddo 


end subroutine 



subroutine testepsilon2D(sten, eps0, eps1, d, nx, ny, ax, bx, ay, by, m)
!!
!! testepsilon2D aprroximates the modified Runge, smoothed Heaviside, and
!! Gelb and Tanner functions with different values of eps0 that
!! are used to bound the interpolant in the case of the PPI method. 
!! This function produces the results used to build the 1D figures 
!! In the manuscript.
!!
!! INPUT
!! sten: stencil selction procedure (sten=1, sten=2, sten=3) 
!! eps0: array of values of eps0 
!! d:  traget polynomial degree for each interpolant
!! n: number of points
!! ax: left boundaries for x 
!! ay: left boundaries for y 
!! bx: right boundaries for x
!! bx: right boundaries for y
!! m: number of output points 
!!

  use mod_adaptiveInterpolation


  implicit none


  integer, intent(in)           :: nx                    !! number of input points
  integer, intent(in)           :: ny                    !! number of input points
  integer, intent(in)           :: m                    !! number of output points
  integer, intent(in)           :: d                    !! target interpolant degree
  integer, intent(in)           :: sten
  real(kind=8), intent(in)      :: ax(3), bx(3)                 !! interval [a, b]
  real(kind=8), intent(in)      :: ay(3), by(3)                 !! interval [a, b]
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

  do k=1, 3
    !!** calculates intreval sizes **!!
    dxn = (bx(k)-ax(k)) /real(nx-1, kind=8)
    dxm = (bx(k)-ax(k)) /real(m-1, kind=8)
    dyn = (by(k)-ay(k)) /real(ny-1, kind=8)
    dym = (by(k)-ay(k)) /real(m-1, kind=8)
  

     !!** uniform mesh **!!
     do i=1,nx-1
       x(i) = ax(k) + real(i-1, kind=8)*dxn
     enddo      
     x(nx)= bx(k)
     do i=1,ny-1  
       y(i) = ay(k) + real(i-1, kind=8)*dyn
     enddo
     y(ny) = by(k)

    !!** output mesh points **!
    do i=1,m-1
      xout(i) = ax(k) + real(i-1, kind=8)*dxm
      yout(i) = ay(k) + real(i-1, kind=8)*dym
    enddo
    xout(m) = bx(k)
    yout(m) = by(k)



    !!** Data values associated to input meshes **!!
    do j=1, ny
      do i=1,nx
        call evalFun2D(k, x(i), y(j), v2D(i, j))
      enddo
    enddo

    !!** True solution **!!
    do j=1, m
      do i=1, m
        call evalFun2D(k, xout(i), yout(j), v2Dout_true(i, j))
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
      if(kk == 7)then
        call adaptiveInterpolation2D(x, y, nx, ny, v2D,  xout, yout, m, m, v2Dout, d, 1, sten, eps1, eps1)
      else
        if(k== 4) then
          call adaptiveInterpolation2D(x, y, nx, ny, v2D,  xout, yout, m, m, v2Dout, d, 2, sten, eps0(kk), eps0(kk))
        else
          call adaptiveInterpolation2D(x, y, nx, ny, v2D,  xout, yout, m, m, v2Dout, d, 2, sten, eps0(kk), eps1)
        endif
      endif

      ii = 1
      do j=1, m
        do i=1, m
          v2D_s(ii, kk+3) = v2Dout(i, j)
          ii = ii+1
        enddo
      enddo
    enddo

    !!** Open file **!! 
    fid = 10                                                      !! file ID
    if(k ==1 )then
      open(unit=fid, file='mapping_data/data/Runge2DEps', status='unknown')
    elseif(k ==2 )then
      open(unit=fid, file='mapping_data/data/Heaviside2DEps', status='unknown')
    elseif(k ==3 )then
      open(unit=fid, file='mapping_data/data/Surface1Eps', status='unknown')
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
!! test002 is used to approximate the Runge, smoothed Heaviside
!! and 2D Terrain functions using different interpolation 
!! methods. This function is used toproduce the 2D results presented
!! in the manuscript.
!!
!! INPUT
!! d: maximum polynomial degree for each interval
!! eps0: positive user-supplied value used to bound interpolant for 
!!       intervalswith no extrema.
!! eps1: positive user-supplied value used to bound interpolant for 
!!       intervals with extrema.
!! sten: user-supplied value used to indicate stencil selection process
!!       possible choices are sten=1, sten=2, sten=3.
!! fun: used to indicate function used
!! nx: number of point in x direction
!! ny: number of points in y direction
!! ax: global interval left boundary in the x direction
!! bx: global right interval boundary in the x direction
!! ay: global interval left boundary in the y direction
!! by: global right interval boundary in the y direction
!! d_el: number of lgl points in each element
!!

  use mod_adaptiveInterpolation


  implicit none


  integer, intent(in)           :: fun                   !! function type 
  integer, intent(in)           :: nx                    !! number of input points
  integer, intent(in)           :: ny                    !! number of input points
  integer, intent(in)           :: m                      !! number of output points
  integer, intent(in)           :: d                      !! target interpolant degree
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
  real(kind=8)			:: y(ny)                 !! input mesh points  
  integer 			:: degx(nx-1, ny)        !! input mesh points  
  real(kind=8)			:: v2D(nx, ny)           !! input data values
  real(kind=8)			:: xout(m)               !! output points to be approximated 
  real(kind=8)			:: yout(m)               !! output points to be approximated 
  real(kind=8)			:: v2Dout(m, m)          !! approximated output values
  real(kind=8)			:: v2Dout_true(m, m)       !! True values at output points
  real(kind=8)			:: v2D_tmp(m, ny)        !! True values at output points
  real(kind=8)			:: dxn, dxm, dyn, dym, err_L2, start_t, end_t


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
    fun_name = "Heaviside2D"
  elseif(fun ==3)then
    fun_name = "T1"
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
  do i=1,nx-1
    x(i) = ax + real(i-1, kind=8)*dxn
  enddo
  x(nx) = bx
  do i=1,ny-1
    y(i) = ay + real(i-1, kind=8)*dyn
  enddo
  y(ny) = by

  !!** number of elements **!!
  dd = d_el

  !** output mesh points **!
  do i=1,m-1
    xout(i) = ax + real(i-1, kind=8)*dxm
    yout(i) = ay + real(i-1, kind=8)*dym
  enddo
  xout(m) = bx
  yout(m) = by


  !!** Data values associated to input meshes **!!
  do j=1, ny
    do i=1,nx
      call evalFun2D(fun, x(i), y(j), v2D(i, j))
    enddo
  enddo

  !!** True solution **!!
  do j=1, m
    do i=1, m
      call evalFun2D(fun, xout(i), yout(j), v2Dout_true(i, j))
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

    !!** Open file **!! 
    fid = 10                                                      !! file ID
    fname = trim("mapping_data/data/")//trim(fun_name)//trim("PCHIP")//trim(fnumber)
    open(unit=fid,file=fname, status='unknown')
    !!** Write to open file **!!
    do j=1, m
      do i=1, m
        write(fid,'(4(3x,E30.15))') xout(i), yout(j), v2Dout_true(i, j), v2Dout(i, j)
      enddo
    enddo
    !!** close file **!!
    close(fid)
  endif
  !!**  Interpolation using Tensor product and DBI **!!
  v2Dout =0.0
  call adaptiveInterpolation2D(x, y, nx, ny, v2D,  xout, yout, m, m, v2Dout, d, 1, sten, eps0, eps1)

  !!** Open file **!! 
  fid = 10                                                      !! file ID
  fname = trim("mapping_data/data/")//trim(fun_name)//trim("DBI")//trim(fnumber)//trim(sst)
  open(unit=fid,file=fname, status='unknown')
  !!** Write to open file **!!
  do j=1, m
    do i=1, m
      write(fid,'(4(3x,E30.15))') xout(i), yout(j), v2Dout_true(i, j), v2Dout(i, j)
    enddo
  enddo
  !!** close file **!!
  close(fid)

  !!**  Interpolation using Tensor product and PPI **!!
  v2Dout = 0.0
  call adaptiveInterpolation2D(x, y, nx, ny, v2D,  xout, yout, m, m, v2Dout, d, 2, sten, eps0, eps1)

  !!** Open file **!! 
  fid = 10                                                      !! file ID
  fname = trim("mapping_data/data/")//trim(fun_name)//trim("PPI")//trim(fnumber)//trim(sst)
  open(unit=fid,file=fname, status='unknown')
  !!** Write to open file **!!
  do j=1, m
    do i=1, m
      write(fid,'(4(3x,E30.15))') xout(i), yout(j), v2Dout_true(i, j), v2Dout(i, j)
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

  qc2 = abs(minval(qc2))+ qc2
  qcp2 = abs(minval(qcp2))+ qcp2

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
  do i=1, nz
    call evalFun1D(1, zd_runge(i), runge2(i))
    call evalFun1D(1, zp_runge(i), rungep2(i))
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

  !!** save Initial data **!!
  do i=1, nz
    ud_out(i,1) = zd(i)
    up_out(i,1) = zp(i)
    ud_pchip_out(i,1) = zd(i)
    up_pchip_out(i,1) = zp(i)
    ud_dbi_out(i,1) = zd(i)
    up_dbi_out(i,1) = zp(i)
  enddo

  iter=1
  !!** save data on physics grid **!!
  do i=1, nz
   up_out(i,iter+1) = u2(i)
   up_pchip_out(i,iter+1) = u2(i)
   up_dbi_out(i,iter+1) = u2(i)
  enddo
  

  !!** Mapping data values from zd (dynamics mesh) to zp (physics mesh) using DBI  **!!
  call adaptiveInterpolation1D(zd, ud_dbi, nz, zp, up_dbi, nz, dd, 1, st, eps0, eps1) 
  
  !!** Mapping data values from zd (dynamics mesh) to zp (physics mesh) using PPI  **!!
  call adaptiveInterpolation1d(zd, ud, nz, zp, up, nz, dd, 2, st, eps0, eps1)

  !!** Mapping data values from zd (dynamics mesh) to zp (physics mesh) using PCHIP  **!!
  call pchez(nz, zd, ud_pchip, d_tmp, spline, wk, nwk, ierr)
  call pchev(nz, zd, ud_pchip, d_tmp, nz, zp, up_pchip, fdl, ierr)


  !!** save data that is on dynamics grid**!!
  do i=1, nz
   ud_out(i,iter+1) = ud(i) 
   ud_pchip_out(i,iter+1) = ud_pchip(i)
   ud_dbi_out(i,iter+1) = ud_dbi(i)
  enddo

  iter = 2
  !!** save data on physics grid **!!
  do i=1, nz
   up_out(i,iter+1) = up(i)
   up_pchip_out(i,iter+1) = up_pchip(i)
   up_dbi_out(i,iter+1) = up_dbi(i)
  enddo


  !!** Mapping data values from zp (physics mesh) to zd (dynamics mesh)  using DBI  **!!
  call adaptiveInterpolation1D(zp, up_dbi, nz, zd(2:nz-1), ud_dbi(2:nz-1), nz-2, dd, 1, st, eps0, eps1) 

  !!** Mapping data values from zp (physics mesh) to zd (dynamics mesh)  using PPI  **!!
  call adaptiveInterpolation1D(zp, up, nz, zd(2:nz-1), ud(2:nz-1), nz-2, dd, 2, st, eps0, eps1)

  !!** Mapping data values from zp (physics mesh) to zd (dynamics mesh)  using PCHIP  **!!
  call pchez(nz, zp, up_pchip, d_tmp, spline, wk, nwk, ierr)
  call pchev(nz, zp, up_pchip, d_tmp, nz-2, zd(2:nz-1), ud_pchip(2:nz-1), fdl, ierr)


  iter = 2
 
  !!** save data that is on dynamics grid **!!
  do i=1, nz
   ud_out(i,iter+1) = ud(i) 
   ud_pchip_out(i,iter+1) = ud_pchip(i)
   ud_dbi_out(i,iter+1) = ud_dbi(i)
  enddo

  write(tmp_str, '("dPPI", i5.5)')fnumber
  fname = trim("mapping_data/data/")//trim(profile_name)//trim(tmp_str)//trim(sst)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(3(1x,E30.15))') (ud_out(i, j), j=1, 3)
  enddo
  close(100)

  write(tmp_str, '("pPPI", i5.5)')fnumber
  fname = trim("mapping_data/data/")//trim(profile_name)//trim(tmp_str)//trim(sst)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.15))') (up_out(i, j), j=1, 3)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname

  write(tmp_str, '("dDBI", i5.5)')fnumber
  fname = trim("mapping_data/data/")//trim(profile_name)//trim(tmp_str)//trim(sst)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.15))') (ud_dbi_out(i, j), j=1, 3)
  enddo
  close(100)

  write(tmp_str, '("pDBI", i5.5)')fnumber
  fname = trim("mapping_data/data/")//trim(profile_name)//trim(tmp_str)//trim(sst)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.15))') (up_dbi_out(i, j), j=1, 3)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname

  write(tmp_str, '("dPCHIP", i5.5)') 3*1000+nz
  fname = trim("mapping_data/data/")//trim(profile_name)//trim(tmp_str)//trim(sst)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.15))') (ud_pchip_out(i, j), j=1, 3)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname


  write(tmp_str, '("pPCHIP", i5.5)') 3*1000+nz
  fname = trim("mapping_data/data/")//trim(profile_name)//trim(tmp_str)//trim(sst)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.15))') (up_pchip_out(i, j), j=1, 3)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname


end subroutine


subroutine evalFun1D(fun, x, v)
!!
!! To evaluate different function. fun determine
!! the function that will be evaluated at the points x
!!
!! INPUT
!! fun: function type
!! x: function input values
!!
!! OUTPUT
!! v: value of function evaluated at  
!!

  implicit none 

  integer, intent(in)                    :: fun                  !! function type
  real(kind=8), intent(in)               :: x                    !! point
  real(kind=8), intent(out)              :: v                    !! point
  real(kind=8)                           :: pi, k                !! temporary variables
  
  
  !!** intialize variables **!!
  k = 100
  pi = 4.0*atan(1.0) 
  
  !!** 1D runge function **!!
  if(fun .eq. 1) then
    v = 0.1 / (0.1 + 25.0 * x * x)

  !!** heaviside function **!!
  else if(fun .eq. 2)then
    v = 1.0/(1.0 + exp(-2*k*x))

  !!** Gelb and Tanner function **!!
  else if(fun .eq. 3)then
    if(x < -0.5) then
      v = 1.0 + (2.0* exp(2.0 * pi *(x+1.0)) - 1.0 -exp(pi)) /(exp(pi)-1.0)
    else
      v = 1.0 - sin(2.0*pi*x / 3.0 + pi/3.0)
    endif
  endif

end subroutine evalFun1D

subroutine evalFun2D(fun, x, y, v)
!!
!! To evaluate different fuunction. fun determine
!! the function that will be evaluated at the points x
!! fun: function type
!! x: function input values
!! y: function input values
!!
!! OUTPUT
!! v: value of function evaluated at  
!! 
  implicit none 

  integer, intent(in)           :: fun                  !! function type
  real(kind=8), intent(in)      :: x                    !! point
  real(kind=8), intent(in)      :: y                    !! point
  real(kind=8), intent(out)     :: v                    !! point
  real(kind=8)                  :: pi                  !! temporary variables
  
  !!** intialize variables **!!
  k = 100
  pi = 4.0*atan(1.0) 
  

  !!** 1D runge function **!!
  if(fun .eq. 1) then
    v = 0.1 / ( 0.1 + 25.0 * ( x*x + y*y) )
  !!** Smoothed Heaviside function
  else if(fun .eq. 2)then
    v = 1.0 / ( 1.0 + exp(-2*k*(x+y)*sqrt(2.0)/2.0) )
  !!** Surface fucntion **!!
  else if(fun .eq. 3)then
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

subroutine performanceEvaluation()
!!
!! Subroutine used to evaluate vectorization for 1D and 2D examples on Intel compiler 
!!
  implicit none 


  integer                       :: n(5)  			!! total number points used 
  integer                       :: d(3)				!! target degree for each interpolant
  integer                       :: i, j, k		
  integer                       :: sten			        !! stencil selection procedure
  real(kind=8)                  :: eps0, eps1                   !! parameters used to bound interpolants
  real(kind=8)                  :: run_time(3), run_times(5, 7)

  n = (/17, 33, 65, 129, 257/)                                              
  d = (/4, 8, 16/)                                                       

  !!** modify eps0 and eps1 to change the bounds on the interpolant **!!
  eps0 = 0.01
  eps1 = 1.0
  sten = 3
  write(*,*) '1D performance results'
  do j=1, 3
    if(j==2) then
      k = 1
    elseif(j==3) then
      k = 3
    elseif(j==4) then
      k = 5
    endif
    do i=1,5
      call performance1D(d(j), n(i), sten, eps0, eps1, n(i)+1, run_time)
      if(j==2) then
        run_times(i,1) = run_time(3)
      endif
      run_times(i,k+1) = run_time(1)
      run_times(i,k+2) = run_time(2)
    enddo
  enddo
  do i=1, 5
    write(*,'(I8, 7(4x, ES15.5))') n(i), run_times(i, 1), run_times(i, 2), run_times(i, 3), &
              run_times(i, 4), run_times(i, 5), run_times(i, 6), run_times(i, 7)
  enddo

  write(*,*) '2D performance results'
  do j=2, 4
    if(j==2) then
      k = 1
    elseif(j==3) then
      k = 3
    elseif(j==4) then
      k = 5
    endif
 
    do i=1,5
      call performance2D(d(j), n(i), sten, eps0, eps1, n(i)+1, run_time)
      if(j==2) then
        run_times(i,1) = run_time(3)
      endif
      run_times(i,k+1) = run_time(1)
      run_times(i,k+2) = run_time(2)
 
    enddo
  enddo
  do i=1, 5
    write(*,'(I8, 7(4x, ES15.5))') n(i), run_times(i, 1), run_times(i, 2), run_times(i, 3), &
              run_times(i, 4), run_times(i, 5), run_times(i, 6), run_times(i, 7)
  enddo


end subroutine

subroutine  performance2D(d, n, sten, eps0, eps1, m, time_data)
!!
!! subroutine used to test vectorization with Intel compiler
!!

  use omp_lib
  use mod_adaptiveInterpolation
  
  implicit none
 
  integer, intent(in) 			:: d
  integer, intent(in) 			:: n
  integer, intent(in) 			:: sten
  integer, intent(in) 			:: m
  real(kind=8), intent(out)		:: time_data(3)
  real(kind=8) 				:: eps0
  real(kind=8) 				:: eps1

  integer 				:: i, j, k
  integer 				:: deg(n-1)
  real(kind=8)				:: runtime
  real(kind=8)				:: dx
 
  real(kind=8)				:: x(n), y(n), v2D(n, n)
  real(kind=8)				:: xout(m), yout(m), v2Dout(m, m), v_tmp(m,n)

  !!** Local variables need for PCHIP **!!
  integer 		        :: nwk, ierr
  real(kind=8)			:: wk((n+1)*2), d_tmp(n+1)
  real(kind=8)			:: fdl(m)
  logical                       :: spline

  spline = .false.  !! needed for PCHIP
  nwk = (n+1)*2     !! needed for PCHIP


  !! 
  dx = 2.0 / real(n-1, kind=8)
  do i=1,n-1
   
    x(i) =-1.0 dx * real(i-1, kind=8)
  enddo
  x(n) = 1.0
  y = x
  do j=1, n
    do i=1, n
      v2D(i,j) = 0.1/(0.1 + 25.0* (x(i)*x(i) + y(j)*y(j)))
    enddo
  enddo
  !! Output mesh !!
  dx = 2.0 / real(m-1, kind=8)
  do i=1,m-1
    xout(i) = -1.0 +dx * real(i-1, kind=8)
  enddo
  xout(m) = 1.0
  yout = xout

  runtime = omp_get_wtime()
  do i=1, 100
    call adaptiveInterpolation2D(x, y, n, n, v2D,  xout, yout, m, m, v2Dout, d, 2, sten, eps0, eps1)
  enddo
  time_data(1) = omp_get_wtime() - runtime
  runtime = omp_get_wtime()
  do i=1, 100
    call adaptiveInterpolation2D_vec(x, y, n, n, v2D,  xout, yout, m, m, v2Dout, d, 2, sten, eps0, eps1)
  enddo
  time_data(2) = omp_get_wtime() - runtime
 
  runtime = omp_get_wtime()
  do i=1, 100
    do j=1, n
      call pchez(n, x, v2D(:,j), d_tmp, spline, wk, nwk, ierr)
      call pchev(n, x, v2D(:, j), d_tmp, m, xout, v_tmp(:, j), fdl, ierr)
    enddo
    do j=1,m
      call pchez(n, y, v_tmp(j,:), d_tmp, spline, wk, nwk, ierr)
      call pchev(n, y, v_tmp(j,:), d_tmp, m, yout, v2Dout(j, :), fdl, ierr)
    enddo
  enddo

  time_data(3) = omp_get_wtime() - runtime

end subroutine 

subroutine  performance1D(d, n, sten, eps0, eps1, m, time_data)
!!
!! subroutine used to test vectorization with Intel compiler
!!

  use omp_lib
  use mod_adaptiveInterpolation
  
  implicit none
 
  integer, intent(in) 			:: d
  integer, intent(in) 			:: n
  integer, intent(in) 			:: sten
  integer, intent(in) 			:: m
  real(kind=8), intent(out)		:: time_data(3)
  real(kind=8) 				:: eps0
  real(kind=8) 				:: eps1

  integer 				:: i, j, k
  integer 				:: deg(n-1)
  real(kind=8)				:: runtime
  real(kind=8)				:: dx
 
  real(kind=8)				:: x(n), v1D(n)
  real(kind=8)				:: xout(m), v1Dout(m)

  !!** Local variables need for PCHIP **!!
  integer 		        :: nwk, ierr
  real(kind=8)			:: wk((n+1)*2), d_tmp(n+1)
  real(kind=8)			:: fdl(m)
  logical                       :: spline

  spline = .false.  !! needed for PCHIP
  nwk = (n+1)*2     !! needed for PCHIP


  !! 
  dx = 2.0 / real(n-1, kind=8)
  !$OMP SIMD
  do i=1,n-1
    x(i) = -1.0 + dx * real(i-1, kind=8)
    v1D(i) = 0.1/(0.1 + 25.0*(x(i)*x(i))
  enddo
  x(n) = 1.0

  !! Output mesh !!
  dx = 2.0 / real(m-1, kind=8)
  do i=1,m-1
    xout(i) = dx * real(i-1, kind=8)
  enddo
  xout(m) = 1.0
  v1Dout = 0.0

  runtime = omp_get_wtime()
  do i=1, 1000
    call adaptiveInterpolation1D(x, v1D, n, xout, v1Dout, m, d, 2, sten, eps0, eps1, deg ) 
  enddo
  time_data(1) = omp_get_wtime() - runtime

  runtime = omp_get_wtime()
  do i=1, 1000
    call adaptiveInterpolation1D_vec(x, v1D, n, xout, v1Dout, m, d, 2, sten, eps0, eps1, deg  ) 
  enddo
  time_data(2) = omp_get_wtime() - runtime

  runtime = omp_get_wtime()
  do i=1, 1000
    call pchez(n, x, v1D, d_tmp, spline, wk, nwk, ierr)
    call pchev(n, x, v1D, d_tmp, m, xout, v1Dout, fdl, ierr)
  enddo
  time_data(3) = omp_get_wtime() - runtime

end subroutine 

subroutine scaleab(vin, vout, n, v_min, v_max, a, b)
!! Scale input data from [v_min v_max] to [a, b]
!!
!! INPUT:
!! n:      number of input elements 
!! vin(n): Input data of size n
!! v_min:  left boundary of input interval
!! v_max:  right boundary of input interval
!! a:      left boundary of output interval
!! b:      right boundary of output interval
!!
!! OUTPUT:
!! vout(n): output data of size n
!!
  integer, intent(in)                :: n                    !! number of elements in vin and vout
  real(kind = 8), intent(in)         :: a, b                 !! [a,b] inerval to scale to 
  real(kind = 8), intent(in)         :: vin(n)               !! input data scaled 
  real(kind = 8), intent(out)        :: vout(n)              !! output data that have been scale to interval [a, b]
  real(kind = 8), intent(in)         :: v_min, v_max
  
  !!** local variables **!!
  integer                               :: i

 do i=1, n
    !! map from [v_min, v_max] to [a, b]
    !! \forall x \in [v_min, v_max], 
    !! map(x) = a + (b-a)/(v_max -v_min)*(x-v_min) 
    vout(i) = a + (b-a)/(v_max-v_min)*(vin(i)-v_min)
  enddo

end subroutine scaleab



