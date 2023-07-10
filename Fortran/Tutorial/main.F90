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
  integer, parameter                    :: m=33
  real(kind=dp)                          :: x(n), y(n), z(n)
  real(kind=dp)                          :: v(n), v2D(n,n), v3D(n,n,n)
  real(kind=dp)                          :: xout(m), yout(m), zout(m)
  real(kind=dp)                          :: vout_apprx(m), vout_apprx2D(m,m), vout_apprx3D(m,m,m)
  real(kind=dp)                          :: vt(m), vt2D(m,m), vt3D(m,m,m)
  real(kind=dp)                          :: eps0, eps1 
  real(kind=dp)                          :: dx 
  real(kind=dp)                          :: tmp 


  
  !-- 1D Tutorial -- !
  dx = 2.0_dp/real(n-1, kind=dp) 
  do i=1, n-1
    x(i) = -1.0_dp + real(i-1, kind=dp)*dx  ! input mesh points
    v(i) = 0.1_dp / (0.1_dp + 25_dp*x(i)**2)                   ! input data values
  enddo
  x(n) = 1.0_dp
  v(n) =  0.1_dp / (0.1_dp + 25_dp*x(n)**2); 
  dx = 2.0_dp/real(m-1, kind=dp) 
  do i=1, m-1
    xout(i) = -1.0_dp + real(i-1, kind=dp)*dx  ! output points
  enddo
  xout(m)= 1.0_dp
  
  do i=1, m
    vt(i) = 0.1_dp / (0.1_dp + 25_dp*xout(i)**2)                   ! input data values
  enddo
  d = 8                           ! target and maximum polynomial degree used for each interval
  interpolation_type = 2          ! 1 for DBI and 2 for PPI
  sten = 1                        ! optional parameter to guide stencil selection 1, 2, and 3
  eps0 = 0.01_dp                  ! optional positive parameter to bound interpolant in PPI
  eps1 = 1.0_dp                   ! optional positive parameter to bound interpolant in PPI

  call adaptiveInterpolation1D_vec(x, v, n, xout, vout_apprx, m, d, interpolation_type, sten, eps0, eps1 ) 
  !-- Display approximated results
  write(*,*) '-- 1D example -- '
  write(*,*) 'i         xout            v            v_apprx            error'
  write(*,*) '------------------------------------------------------------------'
  do i=1,m
      if(mod(i,2)==0 ) then
        write(*,fmt='(I3, 7x, F7.4, 9x, F7.4, 9x, F7.4,9x,  1ES8.2)') i,  &
                xout(i), vt(i), vout_apprx(i), abs(vt(i)-vout_apprx(i))
      endif
  enddo
  write(*,*) '------------------------------------------------------------------'
  tmp = 0.0_dp
  do i=1,m
      if(abs(vt(i)-vout_apprx(i)) > tmp )then
              tmp = abs(vt(i)-vout_apprx(i))
      endif
  enddo
  write(*,fmt='(''The maximum error is '', 1ES8.2)') tmp ;




  !-- 2D Tutorial -- !
  y = x                             
  do j=1,n
    do i=1,n
      v2D(i,j) = 0.1_dp / (0.1_dp + 25_dp*(x(i)**2+y(j)**2))                   ! input data values
    enddo
  enddo
  yout = xout
  
  do j=1,m
    do i=1,m
      vt2D(i,j) = 0.1_dp / (0.1_dp + 25_dp*(xout(i)**2+yout(j)**2)) ! solution data values
    enddo
  enddo
  d = 8                             ! target and maximum polynomial degree used for each interval
  interpolation_type = 2            ! 1 for DBI and 2 for PPI
  sten = 1                          ! optional parameter to guide stencil selection 1, 2, and 3
  eps0 = 0.01_dp                    ! optional positive parameter to bound interpolant in PPI
  eps1 = 1.0_dp                     ! optional positive parameter to bound interpolant in PPI

  call adaptiveInterpolation2D(x, y, n, n, v2D, xout, yout, m, m, vout_apprx2D, d, interpolation_type, sten, eps0, eps1 ) 

  !-- Display approximated results
  write(*,*) '-- 2D example -- '
  write(*,*) 'i        j         xout            yout            v            v_apprx            error'
  write(*,*) '-------------------------------------------------------------------------------------------'
  do i=1,m
    do j=1,m
      if(mod(i,2)==0 .and. mod(j,2) ==0 .and. i==j) then
        write(*,fmt='(I3, 5x, I3, 7x, F7.4, 9x, F7.4, 9x, F7.4, 9x, F7.4,9x,  1ES8.2)') i, j,  &
                xout(i), yout(j), vt2D(i,j), vout_apprx2D(i,j), abs(vt2D(i,j)-vout_apprx2D(i,j))
      endif
    enddo
  enddo
  write(*,*) '-------------------------------------------------------------------------------------------'
  tmp = 0.0_dp
  do i=1,m
    do j=1,m
      if(abs(vt2D(i,j)-vout_apprx2D(i,j)) > tmp )then
              tmp = abs(vt2D(i,j)-vout_apprx2D(i,j))
      endif
    enddo
  enddo
  write(*,fmt='(''The maximum error is '', 1ES8.2)') tmp ;

  

  !-- 3D Tutorial -- !
  z = x
  do k=1,n
    do j=1,n
      do i=1,n
        v3D(i,j,k) = 0.1_dp / (0.1_dp + 25_dp*(x(i)**2+y(j)**2+z(k)**2))                   ! input data values
      enddo
    enddo
  enddo
  zout =xout
  do k=1,m
    do j=1,m
      do i=1,m
        vt3D(i,j,k) = 0.1_dp/(0.1_dp + 25_dp*(xout(i)**2+yout(j)**2+zout(k)**2))   
      enddo
    enddo
  enddo
  
  d = 8                             ! target and maximum polynomial degree used for each interval
  interpolation_type = 2            ! 1 for DBI and 2 for PPI
  sten = 1                          ! optional parameter to guide stencil selection 1, 2, and 3
  eps0 = 0.01_dp                    ! optional positive parameter to bound interpolant in PPI
  eps1 = 1.0_dp                     ! optional positive parameter to bound interpolant in PPI

  call adaptiveInterpolation3D(x, y, z, n, n, n, v3D, xout, yout, zout, m, m, m, vout_apprx3D, &
                                         d, interpolation_type, sten, eps0, eps1 ); 
 
  !-- Display approximated results
  write(*,*) '-- 3D example -- '
  write(*,*) 'i        j        k         xout            yout            zout            ',&
             'v            v_apprx            error'
  write(*,*) '-----------------------------------------------------------------------------',&
             '-------------------------------------'
  do i=1,m
    do j=1,m
      do k=1,m
        if(mod(i,2)==0 .and. mod(j,2) ==0 .and. i==j .and. j==k) then
          write(*,fmt='(I3, 5x, I3, 5x, I3, 7x, F7.4, 9x, F7.4, 9x, F7.4, 9x, F7.4, 9x, F7.4,9x,  1ES8.2)') &
                  i, j, k, xout(i), yout(j), zout(k), vt3D(i,j,k), vout_apprx3D(i,j,k), &
                  abs(vt3D(i,j,k)-vout_apprx3D(i,j,k))
        endif
      enddo
    enddo
  enddo
  write(*,*) '-----------------------------------------------------------------------------',&
             '-------------------------------------'
  tmp = 0.0_dp
  do i=1,m
    do j=1,m
      do k=1,m
        if(abs(vt3D(i,j,k)-vout_apprx3D(i,j,k)) > tmp )then
                tmp = abs(vt3D(i,j,k)-vout_apprx3D(i,j,k))
        endif
      enddo
    enddo
  enddo
  write(*,fmt='(''The maximum error is '', 1ES8.2)') tmp ;


end program
