program main
!!
!! Driver use to produce the approximation and mapping result based on the
!! different function and the TWP-ICE example.
!!
  implicit none
  integer           :: nz(3)
  integer           :: k


  call approximations1D()

  nz = (/64, 127, 253/)
  do k=1,3
    call mapping(nz(k))
  enddo
   
  call approximations2D()

  !! comparing vectorized and unvectorized code on KNL using AVX512 
  !! Intel compiler required 
  call performanceEvaluation()

end program 

subroutine approximations1D()
!!
!! approximation1D is used to set up the 
!! diffferent configuration used to produces
!! the approximation results for the 1D functions
!! presented in the manuscript.
!!
  use mod_adaptiveInterpolation, only: dp
  implicit none 


  integer                       :: n(5)                     !! total number points used 
  integer                       :: d(3)                     !! target degree for each interpolant
  integer                       :: fun(3)                   !! functions used
  integer                       :: i, ii, j, k, kk
  integer                       :: sten                     !! stencil selection procedure
  integer, parameter            :: m = 10000                !! number of output points
  real(dp)                  :: a(3)                     !! intervals left boundary
  real(dp)                  :: b(3)                     !! intervals right boundary
  real(dp)                  :: eps0, eps1, eps_test(6)  !! parameters used to bound interpolants


  !!** Initialization **!!
  n = (/17, 33, 65, 129, 257/)                                              
  d = (/3, 4, 8/)                                                       
  a = (/-1.0_dp, -0.2_dp, -1.0_dp/)
  b = (/ 1.0_dp,  0.2_dp,  1.0_dp/)
  eps_test = (/ 1.0_dp,  0.1_dp,  0.01_dp, 0.001_dp, 0.0001_dp, 0.00_dp /)
  sten = 3

  !!** Modify eps0 and eps1 to change the bounds on the interpolant **!!
  eps0 = 0.01_dp
  eps1 = 1.0_dp

  !!** functions 1=Runge , 2= heaviside, 3=Gelb Tanner **!! 
  fun = (/1, 2, 3/)           

  !!** Used to evaluate different choices of eps0 **!!
  call testepsilon1D(2, eps_test, eps1, d(3), n(1), a, b,  m)
  call testepsilon1D(2, eps_test, eps1, d(2), n(1), a, b,  m)

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
  real(dp), intent(in)      :: a(3), b(3)           !! interval [a, b]
  real(dp), intent(in)      :: eps0(6), eps1        !! test values used for eps0

  integer                       :: i, j, k, fid
  real(dp)                  :: x(n)                 !! uniform input mesh points  
  real(dp)                  :: v1D(n)               !! input data values
  real(dp)                  :: v1Dout(m, 11)         !! output values
  real(dp)                  :: dxn, dxm             !! interval sizes 

  character(len=16)             :: sd
  character(len=64)             :: fname

  !!** Local variables need for PCHIP **!!
  integer                  :: nwk, ierr
  real(dp)             :: wk((n+1)*2), d_tmp(n+1)
  real(dp)             :: fdl(m)
  logical                  :: spline

  !!** Local variables need for MQSI **!!
  real(dp)             :: bcoef(3*n)
  real(dp)             :: t(3*n)
  real(dp)             :: tmp(m)

  spline = .false.  !! needed for PCHIP
  nwk = (n+1)*2     !! needed for PCHIP



  !!** Initialize parameters **!!
  do k=1, 3
    
    !!** calculates intreval sizes **!!
    dxn = (b(k)-a(k)) /real(n-1, dp)
    dxm = (b(k)-a(k)) /real(m-1, dp)
  

    !!** uniform mesh **!!
    do i=1,n-1
      x(i) = a(k) + real(i-1, dp)*dxn
    enddo
    x(n) = b(k)

    !!** output mesh points **!
    do i=1,m-1
      v1Dout(i, 1) = a(k) + real(i-1, dp)*dxm
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

    do i=1, 9
      if(i==7)then
      call adaptiveInterpolation1D(x, v1D, n, v1Dout(:,1), v1Dout(:,2+i), m, d, 1, sten, eps1, eps1 ) 
      elseif(i==8) then
      call pchez(n, x, v1D, d_tmp, spline, wk, nwk, ierr)
      call pchev(n, x, v1D, d_tmp, m, v1Dout(:,1), v1Dout(:,2+i), fdl, ierr)
      elseif(i==9) then
      call mqsi_wrapper(x, v1D, n,  v1Dout(:,1), v1Dout(:,2+i), m)
      else
      call adaptiveInterpolation1D(x, v1D, n, v1Dout(:,1), v1Dout(:,2+i), m, d, 2, sten, eps0(i), eps1) 
      endif
    enddo

    !!** Initialize degree **!!
    if(d == 4) then
      sd = '_4'
    elseif(d ==8) then
      sd= '_8'
    else
      write(*,*) 'The parameter d must be set to 4 or 8 for the varying eps test.'
      stop
    endif

    !!** open file **!!
    fid = 10
    if( k == 1)then
       fname =trim("mapping_data/data/RungeEps")//trim(sd)
       open(unit=fid, file=fname, status='unknown')
    elseif( k == 2)then
       fname =trim("mapping_data/data/HeavisideEps")//trim(sd)
      open(unit=fid, file=fname, status='unknown')
    elseif( k == 3)then
       fname =trim("mapping_data/data/GelbTEps")//trim(sd)
      open(unit=fid, file=fname, status='unknown')
    endif
    !!** write to file **!!
    do i=1, m
      write(fid,'(11(3x,E30.16))') ( v1Dout(i, j), j=1, 11 )
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

  integer, intent(in)      :: fun                  !! function type 
  integer, intent(in)      :: n                    !! number of input points
  integer, intent(in)      :: m                    !! number of output points
  integer, intent(in)      :: d                    !! target interpolant degree
  integer, intent(in)      :: sten
  real(dp), intent(in) :: a                    !! left bounary
  real(dp), intent(in) :: b                    !! right boundary
  real(dp), intent(in) :: eps0 !! parameters used to bound intervals with no hidden extrema  
  real(dp), intent(in) :: eps1 !! parameters used to bound intervals with hidden extrema  
  integer, intent(in)      :: d_el

  integer                  :: ne                   !! number of elments
  integer                  :: i, j, k, fid, ierr, tmp_idx
  integer                  :: is, ie, dd
  real(dp)             :: x(n)                      !! uniform and  LGL input mesh points  
  real(dp)             :: v1D(n)                    !! input data values
  real(dp)             :: xout(m)                   !! output points to be approximated 
  real(dp)             :: v1Dout(m)                 !! approximated output values
  real(dp)             :: v1Dout_true(m)            !! True values at output points
  real(dp)             :: dxn, dxm
  character(len=16)        :: fun_name
  character(len=16)        :: fnumber
  character(len=16)        :: fnumber_mqsi
  character(len=16)        :: sst
  character(len=64)        :: fname

  !!** Local variables need for PCHIP **!!
  integer                  :: nwk
  real(dp)             :: wk((n+1)*2), d_tmp(n+1)
  real(dp)             :: fdl(m)
  logical                  :: spline

  spline = .false.  !! needed for PCHIP
  nwk = (n+1)*2     !! needed for PCHIP

  write(fnumber, '("", i5.5)') d*1000+n
  write(fnumber_mqsi, '("", i5.5)') 5*1000+n

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
    !!call exit(0)
    stop
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
    !!call exit(0)
    stop
  endif
  

  !!** uniform mesh **!!
  dxn = (b-a) /real(n-1, dp)
  do i=1,n-1
    x(i) = a + real(i-1, dp)*dxn
  enddo
  x(n)= b

  dd = d_el
  ne = (n-1) / dd                        !! calculates the number of elements
  
  !!** output mesh points **!
  dxm = (b-a) /real(m-1, dp)
  do i=1,m-1
    xout(i) = a + real(i-1, dp)*dxm
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
 
  if(d ==3 ) then
    !!** interpolation using PCHIP **!!
    call pchez(n, x, v1D, d_tmp, spline, wk, nwk, ierr)
    call pchev(n, x, v1D, d_tmp, m, xout, v1Dout, fdl, ierr)
    !!** open file and write to file **!!
    fid = 10
    fname =trim("mapping_data/data/")//trim(fun_name)//trim("PCHIP")//trim(fnumber)
    open(unit=fid,file=fname, status='unknown')
    do i=1, m
      write(fid,'(3(3x,E30.16))') xout(i), v1Dout_true(i), v1Dout(i)
    enddo
    close(fid)

    !!** interpolation using MQSI **!!
    v1Dout =0.0_dp
    call mqsi_wrapper(x, v1D, n, xout, v1Dout, m)
    !!** open and write to file **!!
    fid = 10
    fname =trim("mapping_data/data/")//trim(fun_name)//trim("MQSI")//trim(fnumber_mqsi)
    open(unit=fid,file=fname, status='unknown')
    do i=1, m
      write(fid,'(3(3x,E30.16))') xout(i), v1Dout_true(i), v1Dout(i)
    enddo
    close(fid)

  endif


  !!** Interpolation using DBI **!!
  v1Dout =0.0_dp
  call adaptiveInterpolation1D(x, v1D, n, xout, v1Dout, m, d, 1, sten, eps0, eps1) 

  !!** open and write to file **!!
  fid = 10
  fname = trim("mapping_data/data/")//trim(fun_name)//trim("DBI")//trim(fnumber)//trim(sst)
  open(unit=fid,file=fname, status='unknown')
  do i=1, m
    write(fid,'(3(3x,E30.16))') xout(i), v1Dout_true(i), v1Dout(i)
  enddo
  close(fid)

  !!** Interpolation using PPI **!!
  v1Dout =0.0_dp
  call adaptiveInterpolation1D(x, v1D, n, xout, v1Dout, m, d, 2, sten, eps0, eps1) 

  !!** open and write to file **!!
  fid = 10
  fname = trim("mapping_data/data/")//trim(fun_name)//trim("PPI")//trim(fnumber)//trim(sst)
  open(unit=fid,file=fname, status='unknown')
  do i=1, m
    write(fid,'(3(3x,E30.16))') xout(i), v1Dout_true(i), v1Dout(i)
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
  use mod_adaptiveInterpolation, only: dp
  implicit none

  integer                       :: nx(5)
  integer                       :: ny(5)
  integer                       :: d(3)
  integer                       :: fun(3)
  integer                       :: i, ii, j, k, kk
  integer                       :: sten 
  integer, parameter            :: m = 1000   !!CHANGE m TO SMALLER VALUE FOR LESS RUNTIME
  real(dp)                  :: ax(3), bx(3)
  real(dp)                  :: ay(3), by(3)
  real(dp)                  :: eps0, eps1, eps_test(6)


  d = (/3, 4, 8/)                  !! array with interpolants degrees
  nx = (/17, 33, 65, 129, 257/)    !! array with number of inputpoints
  ny = (/17, 33, 65, 129, 257/)    !! array with number of inputpoints

  eps_test = (/ 1.0_dp,  0.1_dp,  0.01_dp, 0.001_dp, 0.0001_dp, 0.00_dp /)

  !!** set up interval x \in [ax(i), bx(i)] and y \in [ay(i), by(i)]**!! 
  ax = (/-1.0_dp, -0.2_dp, 0.0_dp /)
  bx = (/ 1.0_dp,  0.2_dp, 2.0_dp /)
  ay = (/-1.0_dp, -0.2_dp, 0.0_dp /)
  by = (/ 1.0_dp,  0.2_dp, 1.0_dp /)
  
  !!** function type 1=runge funtion , 2= heaviside, 3=Gelb Tanner **!! 
  fun = (/1, 2, 3/)               

  !!**
  sten = 3
  eps0 = 0.01_dp
  eps1 = 1.0_dp
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
  real(dp), intent(in)      :: ax(3), bx(3)                 !! interval [a, b]
  real(dp), intent(in)      :: ay(3), by(3)                 !! interval [a, b]
  real(dp), intent(in)      :: eps0(6), eps1

  integer                       :: i, ii, kk, j, k, fid
  real(dp)                  :: x(nx)                 !! input mesh points  
  real(dp)                  :: y(ny)                 !! input mesh points  
  integer                       :: degx2(nx-1, ny)       !! input mesh points  
  integer                       :: degy2(ny-1, m)        !! input mesh points  
  real(dp)                  :: v2D(nx, ny)           !! input data values
  real(dp)                  :: xout(m)               !! output points to be approximated 
  real(dp)                  :: yout(m)               !! output points to be approximated 
  real(dp)                  :: v2Dout(m, m)          !! approximated output values
  real(dp)                  :: v2Dout_true(m, m)       !! True values at output points
  real(dp)                  :: v2D_tmp(m, ny)        !! True values at output points

  real(dp)                  :: v2D_s(m*m, 12)        !! True values at output points
  real(dp)                  :: dxn, dxm, dyn, dym
  real(dp)                  :: h                    !! element spacing

  
  character(len=36)             :: fname
  character(len=16)             :: sd

  if(d == 4) then
    sd = "_4"
  elseif(d == 8) then
    sd = "_8"
  else
      write(*,*) 'The parameter d must be set to 4 or 8 for the varying eps test.'
      stop
  endif

  do k=1, 3
    !!** calculates intreval sizes **!!
    dxn = (bx(k)-ax(k)) /real(nx-1, dp)
    dxm = (bx(k)-ax(k)) /real(m-1, dp)
    dyn = (by(k)-ay(k)) /real(ny-1, dp)
    dym = (by(k)-ay(k)) /real(m-1, dp)
  

     !!** uniform mesh **!!
     do i=1,nx-1
       x(i) = ax(k) + real(i-1, dp)*dxn
     enddo      
     x(nx)= bx(k)
     do i=1,ny-1  
       y(i) = ay(k) + real(i-1, dp)*dyn
     enddo
     y(ny) = by(k)

    !!** output mesh points **!
    do i=1,m-1
      xout(i) = ax(k) + real(i-1, dp)*dxm
      yout(i) = ay(k) + real(i-1, dp)*dym
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

    do kk=1, 9
      v2Dout = 0.0_dp
      v2D_tmp = 0.0_dp
      !!**  Interpolation using Tensor product and DBI **!!
      if(kk == 7)then
        call adaptiveInterpolation2D(x, y, nx, ny, v2D,  xout, yout, m, m, v2Dout, d, 1, sten, eps1, eps1)
      !!**  Interpolation using Tensor product and PCHIP **!!
      elseif(kk == 8)then
        call pchip_wrapper2D(x, y, v2D, nx, ny,  xout, yout, v2Dout, m, m)
      !!**  Interpolation using Tensor product and MQS **!!
      elseif(kk == 9)then
        call mqsi_wrapper2D(x, y, v2D, nx, ny,  xout, yout, v2Dout, m, m)
      !!**  Interpolation using Tensor product and PPI **!!
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
       fname =trim("mapping_data/data/Runge2DEps")//trim(sd)
       open(unit=fid, file=fname, status='unknown')
    elseif(k ==2 )then
       fname =trim("mapping_data/data/Heaviside2DEps")//trim(sd)
       open(unit=fid, file=fname, status='unknown')
    elseif(k ==3 )then
       fname =trim("mapping_data/data/Surface1Eps")//trim(sd)
       open(unit=fid, file=fname, status='unknown')
    endif
    !!** Write to open file **!!
    ii =1
    do j=1, m
      do i=1, m
        write(fid,'(12(3x,E30.16))') ( v2D_s(ii, kk), kk=1, 12 )
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


  integer, intent(in)       :: fun                   !! function type 
  integer, intent(in)       :: nx                    !! number of input points
  integer, intent(in)       :: ny                    !! number of input points
  integer, intent(in)       :: m                     !! number of output points
  integer, intent(in)       :: d                     !! target interpolant degree
  integer, intent(in)       :: d_el                  !! target interpolant degree
  integer, intent(in)       :: sten
  real(dp), intent(in)      :: ax, bx                !! interval [a, b]
  real(dp), intent(in)      :: ay, by                !! interval [a, b]
  real(dp), intent(in)      :: eps0, eps1

  integer                   :: limiter             !!
  integer                   :: i, j, k, fid, ierr, tmp_idx
  integer                   :: ii, jj
  integer                   ::  nwk, seed
  integer                   :: is, ie, nnx, nny, dd
  integer                   :: nex, ney              !! number of elements in x and y directions respectively
  real(dp)                  :: x(nx)                 !! input mesh points  
  real(dp)                  :: y(ny)                 !! input mesh points  
  integer                       :: degx(nx-1, ny)        !! input mesh points  
  real(dp)                  :: v2D(nx, ny)           !! input data values
  real(dp)                  :: xout(m)               !! output points to be approximated 
  real(dp)                  :: yout(m)               !! output points to be approximated 
  real(dp)                  :: v2Dout(m, m)          !! approximated output values
  real(dp)                  :: v2Dout_true(m, m)       !! True values at output points
  real(dp)                  :: v2D_tmp(m, ny)        !! True values at output points
  real(dp)                  :: dxn, dxm, dyn, dym, err_L2, start_t, end_t


  real(dp)                  :: wk((nx+1)*2), d_tmp(nx+1)
  real(dp)                  :: wk2((ny+1)*2), d_tmp2(ny+1)
  real(dp)                  :: tmpin(ny), tmpout(m)
  real(dp)                  :: fdl(m)
  logical                       :: spline

 
  character(len=16)         :: fnumber
  character(len=16)         :: fnumber_mqsi
  character(len=16)         :: sst
  character(len=16)         :: fun_name
  character(len=64)         :: fname

  write(fnumber, '("", i2.2, i3.3, "x", i3.3)')d, nx, ny
  write(fnumber_mqsi, '("", i2.2, i3.3, "x", i3.3)')5, nx, ny

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
    !!call exit(0)
    stop
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
    !!call exit(0)
    stop
  endif
 
  !!** Initialize variables **!!
  spline = .false.


  !!** calculates intreval sizes **!!
  dxn = (bx-ax) /real(nx-1, dp)
  dxm = (bx-ax) /real(m-1, dp)
  dyn = (by-ay) /real(ny-1, dp)
  dym = (by-ay) /real(m-1, dp)
  

  !!** unifnorm mesh **!!
  do i=1,nx-1
    x(i) = ax + real(i-1, dp)*dxn
  enddo
  x(nx) = bx
  do i=1,ny-1
    y(i) = ay + real(i-1, dp)*dyn
  enddo
  y(ny) = by

  !!** number of elements **!!
  dd = d_el

  !** output mesh points **!
  do i=1,m-1
    xout(i) = ax + real(i-1, dp)*dxm
    yout(i) = ay + real(i-1, dp)*dym
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

  if(d == 3) then 
    !!**  Interpolation using Tensor product and PCHIP **!!
    nwk = (nx+1)*2
    do j=1, ny
      call pchez(nx, x, v2D(:,j), d_tmp, spline, wk, nwk, ierr)
      call pchev(nx, x, v2D(:,j), d_tmp, m, xout, v2D_tmp(:, j), fdl, ierr)
    enddo
    nwk = (ny+1)*2
    do i=1, m
      do j=1,ny
        tmpin(j) = v2D_tmp(i,j)
      enddo
      !call pchez(ny, y, v2D_tmp(i,:), d_tmp2, spline, wk2, nwk, ierr)
      !call pchev(ny, y, v2D_tmp(i,:), d_tmp2, m, yout, v2Dout(i, :), fdl, ierr)
      call pchez(ny, y, tmpin, d_tmp2, spline, wk2, nwk, ierr)
      call pchev(ny, y, tmpin, d_tmp2, m, yout, tmpout, fdl, ierr)
      do j=1,m
        v2Dout(i,j) = tmpout(j)
      enddo
    enddo

    !!** Open file **!! 
    fid = 10                                                      !! file ID
    fname = trim("mapping_data/data/")//trim(fun_name)//trim("PCHIP")//trim(fnumber)
    open(unit=fid,file=fname, status='unknown')
    !!** Write to open file **!!
    do j=1, m
      do i=1, m
        write(fid,'(4(3x,E30.16))') xout(i), yout(j), v2Dout_true(i, j), v2Dout(i, j)
      enddo
    enddo
    !!** close file **!!
    close(fid)


    !!**  Interpolation using 2D version of MQSI **!!
    v2Dout =0.0_dp
    call mqsi_wrapper2D(x, y, v2D, nx, ny,  xout, yout, v2Dout, m, m)
    !!** Open file **!! 
    fid = 10                                                      !! file ID
    fname = trim("mapping_data/data/")//trim(fun_name)//trim("MQSI")//trim(fnumber_mqsi)
    open(unit=fid,file=fname, status='unknown')
    !!** Write to open file **!!
    do j=1, m
      do i=1, m
        write(fid,'(4(3x,E30.16))') xout(i), yout(j), v2Dout_true(i, j), v2Dout(i, j)
      enddo
    enddo
    !!** close file **!!
    close(fid)
  endif


  !!**  Interpolation using Tensor product and DBI **!!
  v2Dout =0.0_dp
  call adaptiveInterpolation2D(x, y, nx, ny, v2D,  xout, yout, m, m, v2Dout, d, 1, sten, eps0, eps1)

  !!** Open file **!! 
  fid = 10                                                      !! file ID
  fname = trim("mapping_data/data/")//trim(fun_name)//trim("DBI")//trim(fnumber)//trim(sst)
  open(unit=fid,file=fname, status='unknown')
  !write(*,*) 'fname =', fname
  !!** Write to open file **!!
  do j=1, m
    do i=1, m
      write(fid,'(4(3x,E30.16))') xout(i), yout(j), v2Dout_true(i, j), v2Dout(i, j)
    enddo
  enddo
  !!** close file **!!
  close(fid)

  !!**  Interpolation using Tensor product and PPI **!!
  v2Dout = 0.0_dp
  call adaptiveInterpolation2D(x, y, nx, ny, v2D,  xout, yout, m, m, v2Dout, d, 2, sten, eps0, eps1)

  !!** Open file **!! 
  fid = 10                                                      !! file ID
  fname = trim("mapping_data/data/")//trim(fun_name)//trim("PPI")//trim(fnumber)//trim(sst)
  open(unit=fid,file=fname, status='unknown')
  !!** Write to open file **!!
  do j=1, m
    do i=1, m
      write(fid,'(4(3x,E30.16))') xout(i), yout(j), v2Dout_true(i, j), v2Dout(i, j)
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

  integer               :: nz                   !! number of point used 64 127 253
  integer               :: d(3)                 !! polynomial degree used 
  integer               :: i, j , k             !! iteration ideces
  real(dp)          :: zd(nz),   zp(nz)     !! uniform and LGL mesh
  real(dp)          :: qcp(nz),   qcp2(nz),   qc(nz),   qc2(nz)

  character(len=32)     :: name_runge, name_qc
  character(len=32)     :: sst
  real(dp)          :: zd_runge(nz),   zp_runge(nz)                 !! uniform and LGL mesh
  real(dp)          :: rungep(nz), rungep2(nz), runge(nz), runge2(nz)                 !! data on uniform and LGL mesh
  real(dp)          :: dx, a_runge , b_runge 
    

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
  a_runge = -1.0_dp
  b_runge = 1.0_dp
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

  real(dp), intent(in)  :: zd(nz), zp(nz)                 !! uniform and LGL mesh
  real(dp), intent(in)  :: u(nz), u2(nz)
  character(len=12), intent(in)  :: profile_name
  character(len=12)              :: sst 

  integer               :: fun                  !! determine which fucntion isused
  integer               :: i, j , k             !! iteration ideces
  integer               :: is, ie
  integer               :: limiter
  integer               :: sten
  integer               :: iter
  integer               :: deg(nz-1), deg_dbi(nz-1)
  real(dp)          :: up(nz), ud(nz)     !! data on uniform and LGL mesh
  real(dp)          :: up_pchip(nz), ud_pchip(nz)     !! data on uniform and LGL mesh
  real(dp)          :: up_mqsi(nz), ud_mqsi(nz)     !! data on uniform and LGL mesh
  real(dp)          :: up_dbi(nz), ud_dbi(nz)     !! data on uniform and LGL mesh
  real(dp)          :: up_makima(nz), ud_makima(nz)     !! data on uniform and LGL mesh
  real(dp)          ::  eps0, eps1

  real(dp)          :: dz, xl, xr

  real(dp)          :: wk(nz*2), d_tmp(nz)
  real(dp)          :: fdl(nz)
  logical               :: spline
  integer               :: nwk, ierr
  integer               :: fnumber
  character(len=80)    :: fname, tmp_str


  !!** To save data **!!
  real(dp)          :: ud_out(nz, 3), up_out(nz, 3)
  real(dp)          :: ud_pchip_out(nz, 3), up_pchip_out(nz, 3)
  real(dp)          :: ud_mqsi_out(nz, 3), up_mqsi_out(nz, 3)
  real(dp)          :: ud_dbi_out(nz, 3), up_dbi_out(nz, 3)


  spline = .false.
  nwk = nz*2
  eps0 = 0.01_dp
  eps1 = 1.00_dp

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
    !!call exit(1)
    stop
  endif
 
  ud = u
  ud_pchip = u
  ud_mqsi = u
  ud_dbi = u

  !!** save Initial data **!!
  do i=1, nz
    ud_out(i,1) = zd(i)
    up_out(i,1) = zp(i)
    ud_pchip_out(i,1) = zd(i)
    up_pchip_out(i,1) = zp(i)
    ud_mqsi_out(i,1) = zd(i)
    up_mqsi_out(i,1) = zp(i)
    ud_dbi_out(i,1) = zd(i)
    up_dbi_out(i,1) = zp(i)
  enddo

  iter=1
  !!** save data on physics grid **!!
  do i=1, nz
   up_out(i,iter+1) = u2(i)
   up_pchip_out(i,iter+1) = u2(i)
   up_dbi_out(i,iter+1) = u2(i)
   up_mqsi_out(i,iter+1) = u2(i)
  enddo
  

  !!** Mapping data values from zd (dynamics mesh) to zp (physics mesh) using DBI  **!!
  call adaptiveInterpolation1D(zd, ud_dbi, nz, zp, up_dbi, nz, dd, 1, st, eps0, eps1) 
  
  !!** Mapping data values from zd (dynamics mesh) to zp (physics mesh) using PPI  **!!
  call adaptiveInterpolation1d(zd, ud, nz, zp, up, nz, dd, 2, st, eps0, eps1)

  !!** Mapping data values from zd (dynamics mesh) to zp (physics mesh) using PCHIP  **!!
  call pchez(nz, zd, ud_pchip, d_tmp, spline, wk, nwk, ierr)
  call pchev(nz, zd, ud_pchip, d_tmp, nz, zp, up_pchip, fdl, ierr)


  !!** Mapping data values from zd (dynamics mesh) to zp (physics mesh) using MQSI  **!!
  call mqsi_wrapper(zd, ud_mqsi, nz, zp,  up_mqsi, nz)
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


  !!** Mapping data values from zp (physics mesh) to zd (dynamics mesh)  using MQSI  **!!
  call mqsi_wrapper(zp, up_mqsi, nz, zd(2:nz-1), ud_mqsi(2:nz-1), nz-2)
  iter = 2
 
  !!** save data that is on dynamics grid **!!
  do i=1, nz
   ud_out(i,iter+1) = ud(i) 
   ud_pchip_out(i,iter+1) = ud_pchip(i)
   ud_dbi_out(i,iter+1) = ud_dbi(i)
   ud_mqsi_out(i,iter+1) = ud_mqsi(i)
  enddo

  write(tmp_str, '("dPPI", i5.5)')fnumber
  fname = trim("mapping_data/data/")//trim(profile_name)//trim(tmp_str)//trim(sst)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(3(1x,E30.16))') (ud_out(i, j), j=1, 3)
  enddo
  close(100)

  write(tmp_str, '("pPPI", i5.5)')fnumber
  fname = trim("mapping_data/data/")//trim(profile_name)//trim(tmp_str)//trim(sst)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.16))') (up_out(i, j), j=1, 3)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname

  write(tmp_str, '("dDBI", i5.5)')fnumber
  fname = trim("mapping_data/data/")//trim(profile_name)//trim(tmp_str)//trim(sst)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.16))') (ud_dbi_out(i, j), j=1, 3)
  enddo
  close(100)

  write(tmp_str, '("pDBI", i5.5)')fnumber
  fname = trim("mapping_data/data/")//trim(profile_name)//trim(tmp_str)//trim(sst)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.16))') (up_dbi_out(i, j), j=1, 3)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname

  write(tmp_str, '("dPCHIP", i5.5)') 3*1000+nz
  fname = trim("mapping_data/data/")//trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.16))') (ud_pchip_out(i, j), j=1, 3)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname


  write(tmp_str, '("pPCHIP", i5.5)') 3*1000+nz
  fname = trim("mapping_data/data/")//trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.16))') (up_pchip_out(i, j), j=1, 3)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname

  write(tmp_str, '("dMQSI", i5.5)') 5*1000+nz
  fname = trim("mapping_data/data/")//trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.16))') (ud_pchip_out(i, j), j=1, 3)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname


  write(tmp_str, '("pMQSI", i5.5)') 5*1000+nz
  fname = trim("mapping_data/data/")//trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.16))') (up_pchip_out(i, j), j=1, 3)
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

  use mod_adaptiveInterpolation, only: dp
  implicit none 

  integer, intent(in)                    :: fun                  !! function type
  real(dp), intent(in)               :: x                    !! point
  real(dp), intent(out)              :: v                    !! point
  real(dp)                           :: pi, k                !! temporary variables
  
  
  !!** intialize variables **!!
  k = 100_dp
  pi = 4.0_dp*atan(1.0_dp) 
  
  !!** 1D runge function **!!
  if(fun .eq. 1) then
    v = 0.1_dp / (0.1_dp + 25.0_dp * x * x)

  !!** heaviside function **!!
  else if(fun .eq. 2)then
    v = 1.0_dp/(1.0_dp + exp(-2_dp*k*x))

  !!** Gelb and Tanner function **!!
  else if(fun .eq. 3)then
    if(x < -0.5_dp) then
      v = 1.0_dp + (2.0_dp* exp(2.0_dp * pi *(x+1.0_dp)) - 1.0_dp -exp(pi)) /(exp(pi)-1.0_dp)
    else
      v = 1.0_dp - sin(2.0_dp*pi*x / 3.0_dp + pi/3.0_dp)
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
  use mod_adaptiveInterpolation, only: dp

  implicit none 

  integer, intent(in)           :: fun                  !! function type
  real(dp), intent(in)      :: x                    !! point
  real(dp), intent(in)      :: y                    !! point
  real(dp), intent(out)     :: v                    !! point
  real(dp)                  :: pi,k                  !! temporary variables
  
  !!** intialize variables **!!
  k = 100_dp
  pi = 4.0_dp*atan(1.0_dp) 
  

  !!** 1D runge function **!!
  if(fun .eq. 1) then
    v = 0.1_dp / ( 0.1 + 25.0_dp * ( x*x + y*y) )
  !!** Smoothed Heaviside function
  else if(fun .eq. 2)then
    v = 1.0_dp / ( 1.0_dp + exp(-2_dp*k*(x+y)*sqrt(2.0_dp)/2.0_dp) )
  !!** Surface fucntion **!!
  else if(fun .eq. 3)then
    if( (x-1.5_dp)*(x-1.5_dp) + (y-0.5_dp)*(y-0.5_dp) .le. 1.0_dp/16.0_dp )then
      v = 2.0_dp*0.5_dp*( cos(8.0_dp*atan(1.0_dp)*sqrt((x-1.5_dp)*(x-1.5_dp) &
          + (y-0.5_dp)*(y-0.5_dp))))!!+1)
    else if(y-x .ge. 0.5_dp)then
      v = 1.0_dp
    else if(0.0_dp .le. y-x .and. y-x .le. 0.5_dp)then
      v = 2.0_dp*(y-x)
    else
      v = 0.0_dp
    endif
  endif

end subroutine evalFun2D

subroutine performanceEvaluation()
!!
!! Subroutine used to evaluate vectorization for 1D and 2D examples on Intel compiler 
!!
  use mod_adaptiveInterpolation, only: dp
  implicit none 


  integer                       :: n(5)          !! total number points used 
  integer                       :: d(3)          !! target degree for each interpolant
  integer                       :: i, j, k   
  integer                       :: sten          !! stencil selection procedure
  real(dp)                  :: eps0, eps1    !! parameters used to bound interpolants
  real(dp)                  :: run_time(3), run_times(5, 7)

  n = (/ 17, 33, 65, 127, 257 /)                                              
  d = (/4, 8, 16/)                                                       

  !!** modify eps0 and eps1 to change the bounds on the interpolant **!!
  eps0 = 0.01_dp
  eps1 = 1.0_dp
  sten = 3
  write(*,*) '1D performance results'
  do j=1, 3
    if(j==1) then
      k = 1
    elseif(j==2) then
      k = 3
    elseif(j==3) then
      k = 5
    endif
    do i=1,5
      call performance1D(d(j), n(i), sten, eps0, eps1, n(i)+1, run_time)
      if(j==1) then
        run_times(i,1) = run_time(3)
      endif
      run_times(i,k+1) = run_time(1)
      run_times(i,k+2) = run_time(2)
    enddo
  enddo

  open(100,file='vectorization_results', status='unknown')
  write(100,*) '----  1D runtimes in ms ---- '
  do i=1, 5
    write(100,'(I8, 7(4x, E16.5))') n(i), run_times(i, 1), run_times(i, 2), run_times(i, 3), &
              run_times(i, 4), run_times(i, 5), run_times(i, 6), run_times(i, 7)
  enddo
  close(100)

  write(*,*) '2D performance results'
  do j=1, 3
    if(j==1) then
      k = 1
    elseif(j==2) then
      k = 3
    elseif(j==3) then
      k = 5
    endif
 
    do i=1,5
      call performance2D(d(j), n(i), sten, eps0, eps1, n(i)+1, run_time)
      if(j==1) then
        run_times(i,1) = run_time(3)
      endif
      run_times(i,k+1) = run_time(1)
      run_times(i,k+2) = run_time(2)
 
    enddo
  enddo

  open(100,file='vectorization_results', position='append',status='old',action='write')
  write(100,*) '----  2D runtimes in ms ---- '
  do i=1, 5
    write(100,'(I8, 7(4x, E16.5))') n(i), run_times(i, 1), run_times(i, 2), run_times(i, 3), &
              run_times(i, 4), run_times(i, 5), run_times(i, 6), run_times(i, 7)
  enddo
  close(100)


end subroutine

subroutine  performance2D(d, n, sten, eps0, eps1, m, time_data)
!!
!! subroutine used to test vectorization with Intel compiler
!!

  use omp_lib
  use mod_adaptiveInterpolation
  
  implicit none
 
  integer, intent(in)           :: d
  integer, intent(in)           :: n
  integer, intent(in)           :: sten
  integer, intent(in)           :: m
  real(dp), intent(out)     :: time_data(3)
  real(dp)                  :: eps0
  real(dp)                  :: eps1

  integer                       :: i, j, k
  integer                       :: deg(n-1)
  real(dp)                  :: runtime
  real(dp)                  :: dx
 
  real(dp)                  :: x(n), y(n), v2D(n, n)
  real(dp)                  :: xout(m), yout(m), v2Dout(m, m), v_tmp(m,n)

  !!** Local variables need for PCHIP **!!
  integer                       :: nwk, ierr
  real(dp)                  :: wk((n+1)*2), d_tmp(n+1)
  real(dp)                  :: tmpin(n), tmpout(m)
  real(dp)                  :: fdl(m)
  logical                       :: spline

  spline = .false.  !! needed for PCHIP
  nwk = (n+1)*2     !! needed for PCHIP


  !! 
  dx = 2.0_dp / real(n-1, dp)
  do i=1,n-1
   
    x(i) =-1.0_dp + dx * real(i-1, dp)
  enddo
  x(n) = 1.0_dp
  y = x
  do j=1, n
    do i=1, n
      v2D(i,j) = 0.1_dp/(0.1_dp + 25.0_dp* (x(i)*x(i) + y(j)*y(j)))
    enddo
  enddo
  !! Output mesh !!
  dx = 2.0_dp / real(m-1, dp)
  do i=1,m-1
    xout(i) = -1.0 + dx * real(i-1, dp)
  enddo
  xout(m) = 1.0_dp
  yout = xout

  runtime = omp_get_wtime()
  do i=1, 100
    call adaptiveInterpolation2D(x, y, n, n, v2D,  xout, yout, m, m, v2Dout, d, 2, sten, eps0, eps1)
  enddo
  time_data(1) = (omp_get_wtime() - runtime)*10.0_dp
  runtime = omp_get_wtime()
  do i=1, 100
    call adaptiveInterpolation2D_vec(x, y, n, n, v2D,  xout, yout, m, m, v2Dout, d, 2, sten, eps0, eps1)
  enddo
  time_data(2) = (omp_get_wtime() - runtime)*10.0_dp
 
  if(d ==4) then !! compute once
  runtime = omp_get_wtime()
  do i=1, 100
    do j=1, n
      call pchez(n, x, v2D(:,j), d_tmp, spline, wk, nwk, ierr)
      call pchev(n, x, v2D(:, j), d_tmp, m, xout, v_tmp(:, j), fdl, ierr)
    enddo
    do j=1,m
      do k=1, n
        tmpin(k) = v_tmp(j,k)
      enddo
      !call pchez(n, y, v_tmp(j,:), d_tmp, spline, wk, nwk, ierr)
      !call pchev(n, y, v_tmp(j,:), d_tmp, m, yout, v2Dout(j, :), fdl, ierr)
      call pchez(n, y, tmpin, d_tmp, spline, wk, nwk, ierr)
      call pchev(n, y, tmpin, d_tmp, m, yout, tmpout, fdl, ierr)
      do k=1, m
        v2Dout(j,k) = tmpout(k)
      enddo
    enddo
  enddo
  endif

  time_data(3) = (omp_get_wtime() - runtime)*10.0_dp

end subroutine 

subroutine  performance1D(d, n, sten, eps0, eps1, m, time_data)
!!
!! subroutine used to test vectorization with Intel compiler
!!

  use omp_lib
  use mod_adaptiveInterpolation
  
  implicit none
 
  integer, intent(in)            :: d
  integer, intent(in)            :: n
  integer, intent(in)            :: sten
  integer, intent(in)            :: m
  real(dp), intent(out)      :: time_data(3)
  real(dp)                   :: eps0
  real(dp)                   :: eps1

  integer                        :: i, j, k
  integer                        :: deg(n-1)
  real(dp)                   :: runtime
  real(dp)                   :: dx
 
  real(dp)                   :: x(n), v1D(n)
  real(dp)                   :: xout(m), v1Dout(m)

  !!** Local variables need for PCHIP **!!
  integer                        :: nwk, ierr
  real(dp)                   :: wk((n+1)*2), d_tmp(n+1)
  real(dp)                   :: fdl(m)
  logical                        :: spline

  spline = .false.  !! needed for PCHIP
  nwk = (n+1)*2     !! needed for PCHIP


  !! 
  dx = 2.0_dp / real(n-1, dp)
  !$OMP SIMD
  do i=1,n-1
    x(i) = -1.0_dp + dx * real(i-1, dp)
    v1D(i) = 0.1_dp/(0.1_dp + 25.0_dp*(x(i)*x(i)))
  enddo
  x(n) = 1.0_dp

  !! Output mesh !!
  dx = 2.0_dp / real(m-1, dp)
  do i=1,m-1
    xout(i) = -1.0_dp + dx * real(i-1, dp)
  enddo
  xout(m) = 1.0_dp
  v1Dout = 0.0_dp

  runtime = omp_get_wtime()
  do i=1, 100
    call adaptiveInterpolation1D(x, v1D, n, xout, v1Dout, m, d, 2, sten, eps0, eps1, deg ) 
  enddo
  time_data(1) = (omp_get_wtime() - runtime)*10.0_dp

  runtime = omp_get_wtime()
  do i=1, 100
    call adaptiveInterpolation1D_vec(x, v1D, n, xout, v1Dout, m, d, 2, sten, eps0, eps1, deg  ) 
  enddo
  time_data(2) = (omp_get_wtime() - runtime)*10.0_dp

  if(d==4) then ! compute only once
  runtime = omp_get_wtime()
  do i=1, 100
    call pchez(n, x, v1D, d_tmp, spline, wk, nwk, ierr)
    call pchev(n, x, v1D, d_tmp, m, xout, v1Dout, fdl, ierr)
  enddo
  time_data(3) = (omp_get_wtime() - runtime)*10.0_dp
  endif

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
  use mod_adaptiveInterpolation, only: dp
  integer, intent(in)                :: n                    !! number of elements in vin and vout
  real(dp), intent(in)         :: a, b                 !! [a,b] inerval to scale to 
  real(dp), intent(in)         :: vin(n)               !! input data scaled 
  real(dp), intent(out)        :: vout(n)              !! output data that have been scale to interval [a, b]
  real(dp), intent(in)         :: v_min, v_max
  
  !!** local variables **!!
  integer                               :: i

 do i=1, n
    !! map from [v_min, v_max] to [a, b]
    !! \forall x \in [v_min, v_max], 
    !! map(x) = a + (b-a)/(v_max -v_min)*(x-v_min) 
    vout(i) = a + (b-a)/(v_max-v_min)*(vin(i)-v_min)
  enddo

end subroutine scaleab


subroutine pchip_wrapper(x, v, n,  xout, vout, m)
!!
!!
!!
  use mod_adaptiveInterpolation, only: dp

  implicit none 

  integer                      :: n           !! number of input point
  integer                      :: m           !! number of ouput points
  
  real(dp), intent(in)     :: x(n)        !! input points     
  real(dp), intent(inout)  :: v(n)        !! values at input points     
  real(dp), intent(in)     :: xout(m)     !! output points     
  real(dp), intent(out)    :: vout(m)     !! values at output points     


  !!** Local variables need for PCHIP **!!
  integer                        :: nwk, ierr
  real(dp)                   :: wk((n+1)*2), d_tmp(n+1)
  real(dp)                   :: fdl(m)
  logical                        :: spline

  spline = .false.  !! needed for PCHIP
  nwk = (n+1)*2     !! needed for PCHIP

  call pchez(n, x, v, d_tmp, spline, wk, nwk, ierr)
  call pchev(n, x, v, d_tmp, m, xout, vout, fdl, ierr)

end subroutine

subroutine pchip_wrapper2D(x, y, v, nx, ny,  xout, yout, vout, mx, my)
!!
!!
!!

  use mod_adaptiveInterpolation, only: dp
  implicit none
  
  integer, intent(in)          :: nx, ny      !! number of input point in x and y dimension 
  integer, intent(in)          :: mx, my      !! number of ouput points in xout and yout dimension 
  
  real(dp), intent(in)     :: x(nx)       !! input points in x dimension     
  real(dp), intent(in)     :: y(nx)       !! input points in y dimension    
  real(dp), intent(inout)  :: v(nx, ny)   !! values at input points     
  real(dp), intent(in)     :: xout(mx)    !! output points in x dimension    
  real(dp), intent(in)     ::  yout(my)   !! output points in y dimension    
  real(dp), intent(out)    :: vout(mx, my)!! values at output points
  
  
  integer                      :: i,j
  real(dp)                 :: tmp(mx, ny) 
  real(dp)                 :: tmpin(ny) 
  real(dp)                 :: tmpout(my) 

  !!** interpolate along x **!!
  do j=1,ny
    call  pchip_wrapper(x, v(:,1), nx,  xout, tmp(:,j), mx)
  enddo

  !!** interpolate along y **!!
  do i=1,mx
    do j=1, ny
      tmpin(j) = tmp(i,j)
    enddo
    tmpout = 0.0_dp
    call  pchip_wrapper(y, tmpin, ny,  yout, tmpout, my)
    do j=1, my
     vout(i,j) = tmpout(j)
    enddo
  enddo


end subroutine

subroutine mqsi_wrapper(x, v, n,  xout, vout, m)
!!
!!
!!

  use mod_adaptiveInterpolation, only: dp
  implicit none
  
  integer                      :: n           !! number of input point
  integer                      :: m           !! number of ouput points
  
  real(dp), intent(in)     :: x(n)        !! input points     
  real(dp), intent(inout)  :: v(n)        !! values at input points     
  real(dp), intent(in)     :: xout(m)     !! output points     
  real(dp), intent(out)    :: vout(m)     !! values at output points     
  
  !!** local variables for MQSI algortihm **!!
  real(dp)                 :: bcoef(3*n)
  real(dp)                 :: t(3*n+6)
  real(dp)                 :: uv(n,2)
  integer                      :: info
  
  
  ! Define the interfaces for relevant MQSI package subroutines.
  INTERFACE
   SUBROUTINE MQSI(X, Y, T, BCOEF, INFO, UV)
     use mod_adaptiveInterpolation, only: dp
     REAL(dp), INTENT(IN),  DIMENSION(:) :: X
     REAL(dp), INTENT(INOUT),  DIMENSION(:) :: Y
     REAL(dp), INTENT(OUT), DIMENSION(:) :: T, BCOEF
     INTEGER, INTENT(OUT) :: INFO
     REAL(dp), INTENT(OUT), DIMENSION(:,:), OPTIONAL :: UV
   END SUBROUTINE MQSI
   SUBROUTINE EVAL_SPLINE(T, BCOEF, XY, INFO, D)
     use mod_adaptiveInterpolation, only: dp
     REAL(dp), INTENT(IN), DIMENSION(:) :: T, BCOEF
     REAL(dp), INTENT(INOUT), DIMENSION(:) :: XY
     INTEGER, INTENT(OUT) :: INFO
     INTEGER, INTENT(IN), OPTIONAL :: D
   END SUBROUTINE EVAL_SPLINE
  END INTERFACE
  
  CALL MQSI(x,v,t,bcoef,info,uv) ! Compute monotone quintic spline interpolant
    if(info .ne. 0) then
      write (*,"(/A/)") "MQSI: This test data should not produce an error!"
      write(*,*) "info =", info
      stop
    endif
    vout(1:m) = xout(1:m)
    CALL EVAL_SPLINE(t,bcoef,vout, info,0) ! Evaluate d^(I-1)Q(x)/dx at XY(.).
    if(info .ne. 0) then
      write (*,"(/A/)") "EVAL_SPLINE This test data should not produce an error!"
      write(*,*) "info =", info
      stop
     endif

end subroutine



subroutine mqsi_wrapper2D(x, y, v, nx, ny,  xout, yout, vout, mx, my)
!!
!! This subroutine is wrapper that is used to interface with mqsi_wrapper
!! and the MQSI algorithm.
!! 
!! INPUT:
!!
  use mod_adaptiveInterpolation, only: dp

  implicit none
  
  integer, intent(in)          :: nx, ny      !! number of input point in x and y dimension 
  integer, intent(in)          :: mx, my      !! number of ouput points in xout and yout dimension 
  
  real(dp), intent(in)     :: x(nx)       !! input points in x dimension     
  real(dp), intent(in)     :: y(nx)       !! input points in y dimension    
  real(dp), intent(inout)  :: v(nx, ny)   !! values at input points     
  real(dp), intent(in)     :: xout(mx)    !! output points in x dimension    
  real(dp), intent(in)     ::  yout(my)   !! output points in y dimension    
  real(dp), intent(out)    :: vout(mx, my)!! values at output points
  
  
  integer                  :: i,j
  real(dp)                 :: tmp(mx, ny) 
  real(dp)                 :: tmpin(ny) 
  real(dp)                 :: tmpout(my) 

  !!** interpolate along x **!!
  do j=1,ny
    call  mqsi_wrapper(x, v(:,1), nx,  xout, tmp(:,j), mx)
  enddo

  !!** interpolate along y **!!
  do i=1,mx
    do j=1, ny
      tmpin(j) = tmp(i,j)
    enddo
    tmpout = 0.0_dp
    call  mqsi_wrapper(y, tmpin, ny,  yout, tmpout, my)
    do j=1, my
     vout(i,j) = tmpout(j)
    enddo
  enddo

end subroutine
