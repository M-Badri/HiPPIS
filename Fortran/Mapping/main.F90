program main
!!
!! Driver use to produce the approximation and mapping result based on the
!! different functions and the TWP-ICE example.
!!
  implicit none
  integer           :: nz(3)
  integer           :: k

  !! 1D examples (f_1, f_2, f_3) function approximations 
  call approximations1D()

  !! Mapping examples with f_1 and TWP-ICE data
  nz = (/64, 127, 253/)
  do k=1,3
    call mapping(nz(k))
  enddo
   
  !! 2D examples (f_4, f_5 f_6) function approximations
  call approximations2D()

  !! comparing vectorized and unvectorized code on KNL using AVX512 
  !! Intel compiler and -xMic-AVX512 flag required 
  call performanceEvaluation()

end program 

subroutine approximations1D()
!!
!! approximation1D is used to set up the 
!! different configuration used to produce
!! the approximation results for the 1D functions
!! presented in the manuscript.
!!
  use mod_adaptiveInterpolation, only: dp
  implicit none 


  integer                       :: n(5)                     !! total number of points used 
  integer                       :: d(3)                     !! target degree for each interpolant
  integer                       :: fun(3)                   !! functions used
  integer                       :: i, j, k
  integer                       :: sten                     !! stencil selection procedure
  integer, parameter            :: m = 10000                !! number of output points
  real(kind=dp)                 :: a(3)                     !! intervals left boundary
  real(kind=dp)                 :: b(3)                     !! intervals right boundary
  real(kind=dp)                 :: eps0, eps1, eps_test(6)  !! parameters used to bound interpolants


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
!! The subroutine testepsilon1D(...) approximates the Runge, the smoothed 
!! Heaviside, and the Gelb and Tanner functions with different values of 
!! eps0 that are used to bound the interpolant in the case of the PPI method. 
!! This function produces the results used to build the 1D figures 
!! In the manuscript.
!!
!! INPUT
!! sten: stencil selection procedure (sten=1, sten=2, sten=3) 
!! eps0(6): array of values of eps0 
!! d:  target polynomial degree for each interpolant
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
  real(kind=dp), intent(in)     :: a(3), b(3)           !! interval [a, b]
  real(kind=dp), intent(in)     :: eps0(6), eps1        !! test values used for eps0

  integer                       :: i, j, k, fid
  real(kind=dp)                 :: x(n)                 !! uniform input mesh points  
  real(kind=dp)                 :: v1D(n)               !! input data values
  real(kind=dp)                 :: v1Dout(m, 11)         !! output values
  real(kind=dp)                 :: dxn, dxm             !! interval sizes 

  character(len=16)             :: sd
  character(len=64)             :: fname


  !!** Initialize parameters **!!
  do k=1, 3
    
    !!** calculates intreval sizes **!!
    dxn = (b(k)-a(k)) /real(n-1, kind=dp)
    dxm = (b(k)-a(k)) /real(m-1, kind=dp)
  

    !!** uniform mesh **!!
    do i=1,n-1
      x(i) = a(k) + real(i-1, kind=dp)*dxn
    enddo
    x(n) = b(k)

    !!** output mesh points **!
    do i=1,m-1
      v1Dout(i, 1) = a(k) + real(i-1, kind=dp)*dxm
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
        call pchip_wrapper(x, v1D, n,  v1Dout(:,1), v1Dout(:,2+i), m)
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
!! methods. This function is used to produce the 1D results presented
!! in the manuscript.
!!
!! INPUT
!! d: maximum polynomial degree for each interval
!! eps0: positive user-supplied value used to bound interpolant for 
!!       intervals with no extrema.
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

  integer, intent(in)       :: fun                  !! function type 
  integer, intent(in)       :: n                    !! number of input points
  integer, intent(in)       :: m                    !! number of output points
  integer, intent(in)       :: d                    !! target interpolant degree
  integer, intent(in)       :: sten
  real(kind=dp), intent(in) :: a                    !! left bounary
  real(kind=dp), intent(in) :: b                    !! right boundary
  real(kind=dp), intent(in) :: eps0 !! parameters used to bound intervals with no hidden extrema  
  real(kind=dp), intent(in) :: eps1 !! parameters used to bound intervals with hidden extrema  
  integer, intent(in)       :: d_el

  !integer                   :: ne                   !! number of elements
  integer                   :: i, fid
  !integer                   :: dd
  real(kind=dp)             :: x(n)                  !! uniform and  LGL input mesh points  
  real(kind=dp)             :: v1D(n)                !! input data values
  real(kind=dp)             :: xout(m)               !! output points to be approximated 
  real(kind=dp)             :: v1Dout(m)             !! approximated output values
  real(kind=dp)             :: v1Dout_true(m)        !! True values at output points
  real(kind=dp)             :: dxn, dxm
  character(len=16)         :: fun_name
  character(len=16)         :: fnumber
  character(len=16)         :: fnumber_mqsi
  character(len=16)         :: sst
  character(len=64)         :: fname

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
    stop
  endif

  !!** Get stencil selection procedure **!!
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
  

  i = d_el  ! initalize
  !!** uniform mesh **!!
  dxn = (b-a) /real(n-1, kind=dp)
  do i=1,n-1
    x(i) = a + real(i-1, kind=dp)*dxn
  enddo
  x(n)= b

  !dd = d_el
  !!ne = (n-1) / dd                        !! calculates the number of elements
  
  !!** output mesh points **!
  dxm = (b-a) /real(m-1, kind=dp)
  do i=1,m-1
    xout(i) = a + real(i-1, kind=dp)*dxm
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
    call pchip_wrapper(x, v1D, n, xout, v1Dout, m)
    !!** Ope file and write to file **!!
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
    !!** onOand write to file **!!
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

  !!** pe and write to file **!!
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

  !!** nOp and write to file **!!
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
!! different configurations used to produce
!! the approximation results for the 2D functions
!! presented in the manuscript.
!!
  use mod_adaptiveInterpolation, only: dp
  implicit none

  integer                       :: nx(5)
  integer                       :: ny(5)
  integer                       :: d(3)
  integer                       :: fun(3)
  integer                       :: i, j, k
  integer                       :: sten 
  integer, parameter            :: m = 1000   !!CHANGE m TO SMALLER VALUE FOR LESS RUNTIME
  real(kind=dp)                 :: ax(3), bx(3)
  real(kind=dp)                 :: ay(3), by(3)
  real(kind=dp)                 :: eps0, eps1, eps_test(6)


  d = (/3, 4, 8/)                  !! array with interpolants degrees
  nx = (/17, 33, 65, 129, 257/)    !! array with the number of input points
  ny = (/17, 33, 65, 129, 257/)    !! array with the number of input points

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
  call testepsilon2D(2,eps_test, eps1, d(2), nx(1), ny(1), ax, bx, ay, by, 100)
  call testepsilon2D(2,eps_test, eps1, d(3), nx(1), ny(1), ax, bx, ay, by, 100)
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
!! The subroutine testepsilon2D(...) approximates the modified Runge, the 
!! smoothed Heaviside, and the Gelb and Tanner functions with different 
!! values of eps0 that are used to bound the interpolant in the case of
!! the PPI method. This function produces the results used to build the 
!! 1D figures In the manuscript.
!!
!! INPUT
!! sten: stencil selection procedure (sten=1, sten=2, sten=3) 
!! eps0: array of values of eps0 
!! d:  target polynomial degree for each interpolant
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
  real(kind=dp), intent(in)     :: ax(3), bx(3)                 !! interval [a, b]
  real(kind=dp), intent(in)     :: ay(3), by(3)                 !! interval [a, b]
  real(kind=dp), intent(in)     :: eps0(6), eps1

  integer                       :: i, ii, kk, j, k, fid
  real(kind=dp)                 :: x(nx)                 !! input mesh points  
  real(kind=dp)                 :: y(ny)                 !! input mesh points  
  !integer                       :: degx2(nx-1, ny)       !! input mesh points  
  !integer                       :: degy2(ny-1, m)        !! input mesh points  
  real(kind=dp)                 :: v2D(nx, ny)           !! input data values
  real(kind=dp)                 :: xout(m)               !! output points to be approximated 
  real(kind=dp)                 :: yout(m)               !! output points to be approximated 
  real(kind=dp)                 :: v2Dout(m, m)          !! approximated output values
  real(kind=dp)                 :: v2Dout_true(m, m)       !! True values at output points
  !real(kind=dp)                 :: v2D_tmp(m, ny)        !! True values at output points

  real(kind=dp)                 :: v2D_s(m*m, 12)        !! True values at output points
  real(kind=dp)                 :: dxn, dxm, dyn, dym
  !real(kind=dp)                 :: h                    !! element spacing

  
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
    dxn = (bx(k)-ax(k)) /real(nx-1, kind=dp)
    dxm = (bx(k)-ax(k)) /real(m-1, kind=dp)
    dyn = (by(k)-ay(k)) /real(ny-1, kind=dp)
    dym = (by(k)-ay(k)) /real(m-1, kind=dp)
  

     !!** uniform mesh **!!
     do i=1,nx-1
       x(i) = ax(k) + real(i-1, kind=dp)*dxn
     enddo      
     x(nx)= bx(k)
     do i=1,ny-1  
       y(i) = ay(k) + real(i-1, kind=dp)*dyn
     enddo
     y(ny) = by(k)

    !!** output mesh points **!
    do i=1,m-1
      xout(i) = ax(k) + real(i-1, kind=dp)*dxm
      yout(i) = ay(k) + real(i-1, kind=dp)*dym
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
      !v2D_tmp = 0.0_dp
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
!! The subroutine test002(...) is used to approximate the Runge, the smoothed 
!! Heaviside and 2D Terrain functions using different interpolation 
!! methods. This function is used in the 2D results presented
!! in the manuscript.
!!
!! INPUT
!! d: maximum polynomial degree for each interval
!! eps0: positive user-supplied value used to bound interpolant for 
!!       pe no extrema.
!! eps1: positive user-supplied value used to bound interpolant for 
!!       intervals with extrema.
!! sten: user-supplied value used to indicate stencil selection process
!!       possible choices are sten=1, sten=2, sten=3.
!! fun: used to indicate function used
!! nx: number of points in the x direction
!! ny: number of points in the y direction
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
  real(kind=dp), intent(in) :: ax, bx                !! interval [a, b]
  real(kind=dp), intent(in) :: ay, by                !! interval [a, b]
  real(kind=dp), intent(in) :: eps0, eps1

  integer                   :: i, j, fid!, ierr, tmp_idx
  real(kind=dp)             :: x(nx)                 !! input mesh points  
  real(kind=dp)             :: y(ny)                 !! input mesh points  
  real(kind=dp)             :: v2D(nx, ny)           !! input data values
  real(kind=dp)             :: xout(m)               !! output points to be approximated 
  real(kind=dp)             :: yout(m)               !! output points to be approximated 
  real(kind=dp)             :: v2Dout(m, m)          !! approximated output values
  real(kind=dp)             :: v2Dout_true(m, m)       !! True values at output points
  real(kind=dp)             :: dxn, dxm, dyn, dym!, err_L2, start_t, end_t


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
    stop
  endif
 
  !!** calculate spacing between points **!!
  dxn = (bx-ax) /real(nx-1, kind=dp)
  dxm = (bx-ax) /real(m-1, kind=dp)
  dyn = (by-ay) /real(ny-1, kind=dp)
  dym = (by-ay) /real(m-1, kind=dp)
  

  !!** unifnorm mesh **!!
  do i=1,nx-1
    x(i) = ax + real(i-1, kind=dp)*dxn
  enddo
  x(nx) = bx
  do i=1,ny-1
    y(i) = ay + real(i-1, kind=dp)*dyn
  enddo
  y(ny) = by

  !!** number of elements **!!
  !!dd = d_el
  i = d_el ! initialize

  !** output mesh points **!
  do i=1,m-1
    xout(i) = ax + real(i-1, kind=dp)*dxm
    yout(i) = ay + real(i-1, kind=dp)*dym
  enddo
  xout(m) = bx
  yout(m) = by


  !!** Data values associated with input meshes **!!
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
    v2Dout =0.0_dp
    call pchip_wrapper2D(x, y, v2D, nx, ny,  xout, yout, v2Dout, m, m)

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
!! The subroutine mapping(...) is used to set up the mapping for the Runge and TWP-ICE examples.
!! The following files below are required for the experiment.
!! 'zd_qc_qv_pres_u_v_T_zp_127' and 'zd_qc_qv_pres_u_v_T_zp_253' are obtained by fitting
!! 'zd_qc_qv_pres_u_v_T_zp_64' using a radial basis function interpolation and the evaluating
!! the fitted function at the desired points.
!!
!! FILES
!! 'mapping_data/zd_qc_qv_pres_u_v_T_zp_64': obtained directly from TWP-ICE simulation at
!!   at t = XX s.
!! 'mapping_data/zd_qc_qv_pres_u_v_T_zp_127': obtained by adding a point at the center of each interval 
!! 'mapping_data/zd_qc_qv_pres_u_v_T_zp_253': obtained by adding 3 uniformly spaced points inside each 
!!   interval.
!!  
!! INPUT
!! nz: number of points to be used for the Runge and TWP-ICE examples 
!!

  use mod_adaptiveInterpolation

  implicit none

  integer                :: nz                   !! number of points used 64 127 253
  integer                :: d(3)                 !! polynomial degree used 
  integer                :: i                    !! iteration ideces
  real(kind=dp)          :: zd(nz),   zp(nz)     !! uniform and LGL mesh
  real(kind=dp)          :: qcp(nz),   qcp2(nz),   qc(nz),   qc2(nz)

  character(len=32)      :: name_runge, name_qc
  real(kind=dp)          :: zd_runge(nz),   zp_runge(nz)                 !! uniform and LGL mesh
  real(kind=dp)          :: rungep(nz), rungep2(nz), runge(nz), runge2(nz)                 !! data on uniform and LGL mesh
  real(kind=dp)          :: a_runge , b_runge 
    

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

  !do j=1,3
    do i=1, 3
      runge = runge2
      rungep = rungep2
      qc=qc2
      qcp=qcp2
      write(*, *) '********** d= ', d(i), '**********'
      call mapping2(nz, zd_runge, runge, zp_runge, rungep, d(i), 3, name_runge)
      call mapping2(nz, zd, qc, zp, qcp, d(i), 3, name_qc)
    enddo
  !enddo


end subroutine 

subroutine mapping2(nz, zd, u, zp, u2, dd, st, profile_name)
!!
!!
!! The subroutine mapping2(...) maps data from mesh points zd to zp and back to zp
!!
!! INPUT
!! nz: number of points
!! zd: first mesh points (dynamics mesh points)
!! u: data values associated with the first mesh
!! zp: second mesh points (physics mesh points)
!! u2: data values associated with the second mesh
!! dd: maximum degree used for each interpolant
!! profile_name: profile name to be used to save results
!!
!!

  use mod_adaptiveInterpolation

  implicit none

  integer, intent(in)            :: nz        !! number of points
  integer, intent(in)            :: dd        !! number of points
  integer, intent(in)            :: st        !! number of points

  real(kind=dp), intent(in)      :: zd(nz), zp(nz)                 !! uniform and LGL mesh
  real(kind=dp), intent(in)      :: u(nz), u2(nz)
  character(len=12), intent(in)  :: profile_name
  character(len=12)              :: sst 

  integer                        :: i, j              !! iteration ideces
  integer                        :: iter
  real(kind=dp)                  :: up(nz), ud(nz)     !! data on uniform and LGL mesh
  real(kind=dp)                  :: up_pchip(nz), ud_pchip(nz)     !! data on uniform and LGL mesh
  real(kind=dp)                  :: up_mqsi(nz), ud_mqsi(nz)     !! data on uniform and LGL mesh
  real(kind=dp)                  :: up_dbi(nz), ud_dbi(nz)     !! data on uniform and LGL mesh
  real(kind=dp)                  ::  eps0, eps1

  integer                        :: fnumber
  character(len=80)              :: fname, tmp_str


  !!** To save data **!!
  real(kind=dp)                  :: ud_out(nz, 3), up_out(nz, 3)
  real(kind=dp)                  :: ud_pchip_out(nz, 3), up_pchip_out(nz, 3)
  real(kind=dp)                  :: ud_mqsi_out(nz, 3), up_mqsi_out(nz, 3)
  real(kind=dp)                  :: ud_dbi_out(nz, 3), up_dbi_out(nz, 3)


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
  call pchip_wrapper(zd, ud_pchip, nz,  zp, up_pchip, nz)


  !!** Mapping data values from zd (dynamics mesh) to zp (physics mesh) using MQSI  **!!
  call mqsi_wrapper(zd, ud_mqsi, nz, zp,  up_mqsi, nz)
  !!** save data that is on dynamics grid**!!
  do i=1, nz
   ud_out(i,iter+1) = ud(i) 
   ud_pchip_out(i,iter+1) = ud_pchip(i)
   ud_mqsi_out(i,iter+1) = ud_mqsi(i)
   ud_dbi_out(i,iter+1) = ud_dbi(i)
  enddo

  iter = 2
  !!** save data on physics grid **!!
  do i=1, nz
   up_out(i,iter+1) = up(i)
   up_pchip_out(i,iter+1) = up_pchip(i)
   up_mqsi_out(i,iter+1) = up_mqsi(i)
   up_dbi_out(i,iter+1) = up_dbi(i)
  enddo


  !!** Mapping data values from zp (physics mesh) to zd (dynamics mesh)  using DBI  **!!
  call adaptiveInterpolation1D(zp, up_dbi, nz, zd(2:nz-1), ud_dbi(2:nz-1), nz-2, dd, 1, st, eps0, eps1) 

  !!** Mapping data values from zp (physics mesh) to zd (dynamics mesh)  using PPI  **!!
  call adaptiveInterpolation1D(zp, up, nz, zd(2:nz-1), ud(2:nz-1), nz-2, dd, 2, st, eps0, eps1)

  !!** Mapping data values from zp (physics mesh) to zd (dynamics mesh)  using PCHIP  **!!
  call pchip_wrapper(zp, up_pchip, nz, zd(2:nz-1), ud_pchip(2:nz-1), nz-2)


  !!** Mapping data values from zp (physics mesh) to zd (dynamics mesh)  using MQSI  **!!
  call mqsi_wrapper(zp, up_mqsi, nz, zd(2:nz-1), ud_mqsi(2:nz-1), nz-2)
  iter = 2
 
  !!** Save data that is on the dynamics grid **!!
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
    write(100, '(1002(1x,E30.16))') (ud_mqsi_out(i, j), j=1, 3)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname


  write(tmp_str, '("pMQSI", i5.5)') 5*1000+nz
  fname = trim("mapping_data/data/")//trim(profile_name)//trim(tmp_str)
  open(100,file=fname, status='unknown')
  do i=1, nz
    write(100, '(1002(1x,E30.16))') (up_mqsi_out(i, j), j=1, 3)
  enddo
  close(100)
  write(*,*) 'Saved in file name ', fname


end subroutine


subroutine evalFun1D(fun, x, v)
!!
!! the subroutine evalFun1D computes different 1D functions. 
!! fun determine the function that will be evaluated at the points x
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
  real(kind=dp), intent(in)              :: x                    !! point
  real(kind=dp), intent(out)             :: v                    !! point
  real(kind=dp)                          :: pi, k                !! temporary variables
  
  
  !!** intialize variables **!!
  k = 100_dp
  pi = 4.0_dp*atan(1.0_dp) 
  
  !!** 1D runge function **!!
  if(fun .eq. 1) then
    v = 0.1_dp / (0.1_dp + 25.0_dp * x * x)
    !v = 1.0_dp / (1.0_dp + 25.0_dp * x * x)
    !v =  cos( pi *x) 

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
!! The subroutine evalFun2D computes different 2D functions. 
!! fun determine the function that will be evaluated at the points x
!!
!! INPUT
!!fun: function type
!! x: function input values
!! y: function input values
!!
!! OUTPUT
!! v: value of function evaluated at  
!! 
  use mod_adaptiveInterpolation, only: dp

  implicit none 

  integer, intent(in)            :: fun                  !! function type
  real(kind=dp), intent(in)      :: x                    !! point
  real(kind=dp), intent(in)      :: y                    !! point
  real(kind=dp), intent(out)     :: v                    !! point
  real(kind=dp)                  :: k                  !! temporary variables
  
  !!** intialize variables **!!
  k = 100_dp
  !pi = 4.0_dp*atan(1.0_dp) 
  

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
!! The subroutine performanceEvaluation() evaluates runtime for the 
!! unvectorized and vectorization for 1D and 2D examples on Intel compiler 
!!
  use mod_adaptiveInterpolation, only: dp
  implicit none 


  integer                   :: n(5)          !! total number of points used 
  integer                   :: d(3)          !! target degree for each interpolant
  integer                   :: i, j, k, fun   
  integer                   :: sten          !! stencil selection procedure
  real(kind=dp)             :: eps0, eps1    !! parameters used to bound interpolants
  real(kind=dp)             :: run_time(4), run_times(5, 8)

  n = (/ 17, 33, 65, 127, 257 /)                                              
  n = (/ 65, 129, 257, 512, 1025 /)                                              
  d = (/4, 8, 16/)                                                       

  !!** modify eps0 and eps1 to change the bounds on the interpolant **!!
  eps0 = 0.01_dp
  eps1 = 1.0_dp
  sten = 3
  do fun=1,3
  write(*,*) '1D performance results fun =', fun
  do j=1, 3
    if(j==1) then
      k = 1
    elseif(j==2) then
      k = 3
    elseif(j==3) then
      k = 5
    endif
    do i=1,5
      call performance1D(fun, d(j), n(i), sten, eps0, eps1, n(i)+1, run_time)
      if(j==1) then
        run_times(i,1) = run_time(3)
        run_times(i,8) = run_time(4)
      endif
      run_times(i,k+1) = run_time(1)
      run_times(i,k+2) = run_time(2)
    enddo
  enddo

  if(fun ==1)then
    open(100,file='vectorization_results', status='unknown')
  else
    open(100,file='vectorization_results', position='append',status='old',action='write')
  endif
  write(100,*) '----  1D runtimes in ms ---- '
  do i=1, 5
    write(100,'(I8, 8(4x, E16.5))') n(i), run_times(i, 1), run_times(i, 8), run_times(i, 2),&
            run_times(i, 3), run_times(i, 4), run_times(i, 5), run_times(i, 6), run_times(i, 7)
  enddo
  close(100)
  enddo
  !!
  !!
  do fun=1,3
  write(*,*) '2D performance results fun =', fun
  do j=1, 3
    if(j==1) then
      k = 1
    elseif(j==2) then
      k = 3
    elseif(j==3) then
      k = 5
    endif
 
    do i=1,5
      call performance2D(fun, d(j), n(i), sten, eps0, eps1, n(i)+1, run_time)
      if(j==1) then
        run_times(i,1) = run_time(3)
        run_times(i,8) = run_time(4)
      endif
      run_times(i,k+1) = run_time(1)
      run_times(i,k+2) = run_time(2)
 
    enddo
  enddo

  open(100,file='vectorization_results', position='append',status='old',action='write')
  write(100,*) '----  2D runtimes in ms ---- '
  do i=1, 5
    write(100,'(I8, 8(4x, E16.5))') n(i), run_times(i, 1), run_times(i, 8), run_times(i, 2),&
            run_times(i, 3), run_times(i, 4), run_times(i, 5), run_times(i, 6), run_times(i, 7)
  enddo
  close(100)
  enddo


end subroutine

subroutine  performance2D(fun, d, n, sten, eps0, eps1, m, time_data)
!!
!! The subroutine performance2D(...) evaluates the runtimes for the vectorized and unvectorized 
!! 2D examples.
!!

  use omp_lib
  use mod_adaptiveInterpolation
  
  implicit none
 
  integer, intent(in)            :: fun
  integer, intent(in)            :: d
  integer, intent(in)            :: n
  integer, intent(in)            :: sten
  integer, intent(in)            :: m
  real(kind=dp), intent(out)     :: time_data(4)
  real(kind=dp)                  :: eps0
  real(kind=dp)                  :: eps1

  integer                        :: i, j, idx
  !integer                        :: deg(n-1)
  real(kind=dp)                  :: runtime
 
  real(kind=dp)                  :: x(n), y(n), v2D(n, n)
  real(kind=dp)                  :: xout(m), yout(m), v2Dout(m, m)!, v_tmp(m,n)
  real(kind=dp)                  :: ax(3), bx(3), ay(3), by(3)
  real(kind=dp)                  :: dxn, dxm, dyn, dym

  idx = fun
  !!** set up interval x \in [ax(i), bx(i)] and y \in [ay(i), by(i)]**!! 
  ax = (/-1.0_dp, -0.2_dp, 0.0_dp /)
  bx = (/ 1.0_dp,  0.2_dp, 2.0_dp /)
  ay = (/-1.0_dp, -0.2_dp, 0.0_dp /)
  by = (/ 1.0_dp,  0.2_dp, 1.0_dp /)
  
  !!** calculates intreval sizes **!!
  dxn = (bx(idx)-ax(idx)) /real(n-1, kind=dp)
  dxm = (bx(idx)-ax(idx)) /real(m-1, kind=dp)
  dyn = (by(idx)-ay(idx)) /real(n-1, kind=dp)
  dym = (by(idx)-ay(idx)) /real(m-1, kind=dp)
  

  !!** unifnorm mesh **!!
  do i=1,n-1
    x(i) = ax(idx) + real(i-1, kind=dp)*dxn
  enddo
  x(n) = bx(idx)
  do i=1,n-1
    y(i) = ay(idx) + real(i-1, kind=dp)*dyn
  enddo
  y(n) = by(idx)

  !** output mesh points **!
  do i=1,m-1
    xout(i) = ax(idx) + real(i-1, kind=dp)*dxm
    yout(i) = ay(idx) + real(i-1, kind=dp)*dym
  enddo
  xout(m) = bx(idx)
  yout(m) = by(idx)


  !!** Data values associated to input meshes **!!
  do j=1, n
    do i=1,n
      call evalFun2D(fun, x(i), y(j), v2D(i, j))
    enddo
  enddo

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
      call pchip_wrapper2D(x, y, v2D, n, n,  xout, yout, v2Dout, m, m)
    enddo
    time_data(3) = (omp_get_wtime() - runtime)*10.0_dp

    runtime = omp_get_wtime()
    do i=1, 100
      call mqsi_wrapper2D(x, y, v2D, n, n,  xout, yout, v2Dout, m, m)
    enddo
    time_data(4) = (omp_get_wtime() - runtime)*10.0_dp
  endif

end subroutine 

subroutine  performance1D(fun, d, n, sten, eps0, eps1, m, time_data)
!!
!! The subroutine performance1D(...) evaluates the runtimes for the vectorized and unvectorized 
!! 2D examples.
!!

  use omp_lib
  use mod_adaptiveInterpolation
  
  implicit none
 
  integer, intent(in)             :: d
  integer, intent(in)             :: fun
  integer, intent(in)             :: n
  integer, intent(in)             :: sten
  integer, intent(in)             :: m
  real(kind=dp), intent(out)      :: time_data(4)
  real(kind=dp)                   :: eps0
  real(kind=dp)                   :: eps1

  integer                         :: i,idx
  integer                         :: deg(n-1)
  real(kind=dp)                   :: runtime
  !real(kind=dp)                   :: dx
 
  real(kind=dp)                   :: x(n), v1D(n)
  real(kind=dp)                   :: xout(m), v1Dout(m)
  real(kind=dp)                   :: a(3), b(3), dxn, dxm

  idx =fun
  a = (/-1.0_dp, -0.2_dp, -1.0_dp/)
  b = (/ 1.0_dp,  0.2_dp,  1.0_dp/)

  !!** uniform mesh **!!
  dxn = (b(idx)-a(idx)) /real(n-1, kind=dp)
  do i=1,n-1
    x(i) = a(idx) + real(i-1, kind=dp)*dxn
  enddo
  x(n)= b(idx)

  !!** output mesh points **!
  dxm = (b(idx)-a(idx)) /real(m-1, kind=dp)
  do i=1,m-1
    xout(i) = a(idx) + real(i-1, kind=dp)*dxm
  enddo
  xout(m) = b(idx)

  !!** Data values associated to input meshes **!!
  do i=1,n
    call evalFun1D(fun, x(i), v1D(i))
  enddo

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
    call pchip_wrapper(x, v1D, n,  xout, v1Dout, m)
  enddo
  time_data(3) = (omp_get_wtime() - runtime)*10.0_dp
 
  runtime = omp_get_wtime()
  do i=1, 100
    call mqsi_wrapper(x, v1D, n,  xout, v1Dout, m)
  enddo
  time_data(4) = (omp_get_wtime() - runtime)*10.0_dp
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
!! a:      left boundary of the output interval
!! b:      right boundary of the output interval
!!
!! OUTPUT:
!! vout(n): output data of size n
!!
  use mod_adaptiveInterpolation, only: dp

  integer, intent(in)               :: n                    !! number of elements in vin and vout
  real(kind=dp), intent(in)         :: a, b                 !! [a,b] inerval to scale to 
  real(kind=dp), intent(in)         :: vin(n)               !! input data scaled 
  real(kind=dp), intent(out)        :: vout(n)              !! output data that have been scaled to interval [a, b]
  real(kind=dp), intent(in)         :: v_min, v_max
  
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
!! The subroutine pchip_wrapper(...) is used to interface with
!! piece wise cubic interpolation (PCHIP ) algorithm
!! to construct a polynomial for each interval and evaluate 
!! the constructed polynomial at the desired output points xout
!!
!! INPUT
!! x: 1D vector that holds input mesh points
!! v: 1D vector that holds data values associated to x
!! n: number of pints in x
!! xout: 1D vector that holds the output mesh points
!! m: number of points in xout
!!
!! OUTPUT
!! vout: 1D vector to  hold the data values associated with xout
!!

  use mod_adaptiveInterpolation, only: dp

  implicit none 

  integer                       :: n           !! number of input point
  integer                       :: m           !! number of ouput points
  
  real(kind=dp), intent(in)     :: x(n)        !! input points     
  real(kind=dp), intent(inout)  :: v(n)        !! values at input points     
  real(kind=dp), intent(in)     :: xout(m)     !! output points     
  real(kind=dp), intent(out)    :: vout(m)     !! values at output points     


  !!** Local variables need for PCHIP **!!
  integer                       :: nwk, ierr
  real(kind=dp)                 :: wk((n+1)*2), d_tmp(n+1)
  real(kind=dp)                 :: fdl(m)
  logical                       :: spline

  spline = .false.  !! needed for PCHIP
  nwk = (n+1)*2     !! needed for PCHIP

  call pchez(n, x, v, d_tmp, spline, wk, nwk, ierr)
  call pchev(n, x, v, d_tmp, m, xout, vout, fdl, ierr)

end subroutine

subroutine pchip_wrapper2D(x, y, v, nx, ny,  xout, yout, vout, mx, my)
!!
!! This subroutine is a wrapper that is used to interface with pchip_wrapper
!! and for 2D piecewise bi-cubic spline interpolation.
!! 
!! INPUT
!! nx: number of points in the 1D vector x
!! ny: number of points in the 1D vector y
!! x: 1D vector with discretization along x-axis 
!! y: 1D vector with discretization along y-axis 
!! v: 2D vector of size nx*nz with the datavalues associated with the 
!!    mesh obtained from the tensor product nx*ny
!! mx: number of points in the 1D vector xout
!! my: number of points in the 1D vector yout
!! xout: 1D vector with output points along the x-axis
!! yout: 1D vector with output points along the y-axis
!!
!!OUTPUT
!! vout: 2D vector of size mx*my to hold interpolated results
!!
 
  use mod_adaptiveInterpolation, only: dp
  implicit none
  
  integer, intent(in)           :: nx, ny      !! number of input point in x and y dimension 
  integer, intent(in)           :: mx, my      !! number of output points in xout and yout dimension 
  
  real(kind=dp), intent(in)     :: x(nx)       !! input points in x dimension     
  real(kind=dp), intent(in)     :: y(nx)       !! input points in y dimension    
  real(kind=dp), intent(inout)  :: v(nx, ny)   !! values at input points     
  real(kind=dp), intent(in)     :: xout(mx)    !! output points in x dimension    
  real(kind=dp), intent(in)     ::  yout(my)   !! output points in y dimension    
  real(kind=dp), intent(out)    :: vout(mx, my)!! values at output points
  
  
  integer                       :: i,j
  real(kind=dp)                 :: tmp(mx, ny) 
  real(kind=dp)                 :: tmpin(ny) 
  real(kind=dp)                 :: tmpout(my) 

  !!** interpolate along x **!!
  do j=1,ny
    call  pchip_wrapper(x, v(:,j), nx,  xout, tmp(:,j), mx)
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
!! The subroutine mqsi_wrapper(...) is used to interface with
!! the monotonic quintic spline interpolation (MQSI) algorithm
!! to construct a polynomial for each interval and evaluate 
!! the constructed polynomial at the desired output points xout
!!
!! INPUT
!! x: 1D vector that holds input mesh points
!! v: 1D vector that holds data values associated to x
!! n: number of pints in x
!! xout: 1D vector that holds the output mesh points
!! m: number of points in xout
!!
!! OUTPUT
!! vout: 1D vector to  hold the data values associated with xout
!!
  use mod_adaptiveInterpolation, only: dp
  implicit none
  
  integer                       :: n       !! number of input point
  integer                       :: m       !! number of ouput points
  
  real(kind=dp), intent(in)     :: x(n)        !! input points     
  real(kind=dp), intent(inout)  :: v(n)        !! values at input points     
  real(kind=dp), intent(in)     :: xout(m)     !! output points     
  real(kind=dp), intent(out)    :: vout(m)     !! values at output points     
  
  !!** local variables for MQSI algortihm **!!
  real(kind=dp)                 :: bcoef(3*n)
  real(kind=dp)                 :: t(3*n+6)
  !real(kind=dp)                 :: uv(n,2)
  integer                       :: info
  
  
  ! Define the interfaces for relevant MQSI package subroutines.
  interface
   subroutine MQSI(x, y, t, bcoef, info, uv)
     use mod_adaptiveInterpolation, only: dp
     real(kind=dp), intent(in),  dimension(:) :: x
     real(kind=dp), intent(inout),  dimension(:) :: y
     real(kind=dp), intent(out), dimension(:) :: t, bcoef
     integer, intent(out) :: info
     real(kind=dp), intent(out), dimension(:,:), optional :: uv
   end subroutine MQSI
   subroutine EVAL_SPLINE(t, bcoef, xy, info, d)
     use mod_adaptiveInterpolation, only: dp
     real(kind=dp), intent(in), dimension(:) :: t, bcoef
     real(kind=dp), intent(inout), dimension(:) :: xy
     integer, intent(out) :: info
     integer, intent(in), optional :: d
   end subroutine EVAL_SPLINE
  end interface
  
  CALL MQSI(x,v,t,bcoef,info) ! Compute monotone quintic spline interpolant
    if(info .ne. 0) then
      write (*,"(/A/)") "MQSI: This test data should not produce an error!"
      write(*,*) "info =", info
      stop
    endif
    vout(1:m) = xout(1:m)
    CALL EVAL_SPLINE(t,bcoef,vout, info) ! Evaluate d^(I-1)Q(x)/dx at XY(.).
    if(info .ne. 0) then
      write (*,"(/A/)") "EVAL_SPLINE This test data should not produce an error!"
      write(*,*) "info =", info
      stop
     endif

end subroutine



subroutine mqsi_wrapper2D(x, y, v, nx, ny,  xout, yout, vout, mx, my)
!!
!! This subroutine is a wrapper that is used to interface with mqsi_wrapper
!! and the MQSI algorithm.
!! 
!! INPUT
!! nx: number of points in the 1D vector x
!! ny: number of points in the 1D vector y
!! x: 1D vector with discretization along x-axis 
!! y: 1D vector with discretization along y-axis 
!! v: 2D vector of size nx*nz with the datavalues associated with the 
!!    mesh obtained from the tensor product nx*ny
!! mx: number of points in the 1D vector xout
!! my: number of points in the 1D vector yout
!! xout: 1D vector with output points along the x-axis
!! yout: 1D vector with output points along the y-axis
!!
!!OUTPUT
!! vout: 2D vector of size mx*my to hold interpolated results
!!
  use mod_adaptiveInterpolation, only: dp

  implicit none
  
  integer, intent(in)           :: nx, ny      !! number of input point in x and y dimension 
  integer, intent(in)           :: mx, my      !! number of output points in xout and yout dimension 
  
  real(kind=dp), intent(in)     :: x(nx)       !! input points in x dimension     
  real(kind=dp), intent(in)     :: y(nx)       !! input points in y dimension    
  real(kind=dp), intent(inout)  :: v(nx, ny)   !! values at input points     
  real(kind=dp), intent(in)     :: xout(mx)    !! output points in x dimension    
  real(kind=dp), intent(in)     ::  yout(my)   !! output points in y dimension    
  real(kind=dp), intent(out)    :: vout(mx, my)!! values at output points
  
  
  integer                       :: i,j
  real(kind=dp)                 :: tmp(mx, ny) 
  real(kind=dp)                 :: tmpin(ny) 
  real(kind=dp)                 :: tmpout(my) 

  !!** interpolate along x **!!
  do j=1,ny
    call  mqsi_wrapper(x, v(:,j), nx,  xout, tmp(:,j), mx)
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
