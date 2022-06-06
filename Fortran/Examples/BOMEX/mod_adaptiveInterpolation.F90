!!
!! Last modifided: 
!!


module mod_adaptiveInterpolation
!!
!!  Module for adaptive polynomial Interpolation  
!!

  implicit none


  contains

subroutine divdiff2(x, y, n, u)
!!
!! Computes divided difference using the expanded form
!!

  implicit none

  integer, intent(in)                   :: n            !! number points
  real(kind=8), intent(in)              :: x(n)         !! mesh points
  real(kind=8), intent(in)              :: y(n)         !! mesh points
  real(kind=8), intent(out)             :: u         !! mesh points
  

  integer                               :: j, k
  real(kind=8)                          :: tmp

  u = 0.0
  do j=1, n
     tmp = 1.0
     do k=1, n
       if(k .ne. j) then
         tmp = tmp * (x(j)-x(k))
       end if
     end do
     u = u + y(j)/tmp
  enddo

end subroutine divdiff2


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


subroutine divdiff(x, y, n, d, table)
!!! This subroutine computes the table of divided differences
!! 
!! INPUT:
!! x: 1D vector of x coorrdinates
!! y: 1D vector of y coorrdinates. y are the data values associated to the locations x.
!! d: maximum number of consecutive mesh points used to compute divided differences.
!! n: number of points in x.
!!
!! OUTPUT:
!!
!! table: array of dimension of n X (d+1)
!! table = u[1], u[1,2], u[1,3], ... u[1,d+1]
!!         u[2], u[2,3], u[2,4], ... u[2,d+2]
!!         u[3], u[3,4], u[3,5], ... u[2,d+3]
!!          .  ,    .  ,    .  , ...    .
!!          .  ,    .  ,    .  , ...    .
!!          .  ,    .  ,    .  , ...    .
!!         u[n],    0  ,    0  , ...    .

  integer, intent(in)		:: n,d
  real(kind=8), intent(in)		:: x(n), y(n)
  real(kind=8), intent(out)		:: table(n,d+1)
  integer			:: i, j
  real(kind=8)                  :: tmp


  do i=1,n
    table(i,1) = y(i);
  enddo
  do j=2,d+1
    do i=1,n-(j-1)
      call divdiff2(x(i:i+j-1), y(i:i+j-1), j, tmp)
      table(i,j) = tmp

      !table(i,j) = (table(i+1, j-1)-table(i, j-1)) / (x(i+j-1)-x(i))
      !!if(abs(table(i, j) - tmp) > 1e-10 ) then
      !!  write(*,*) 'TAJO tmp =', tmp, 'table(i,j) =', table(i,j)
      !!  write(*,*) 'i=',i, 'j=', j
      !!endif
    enddo
  enddo

end subroutine 
 


subroutine newtonPolyVal(x, u, d, xout, yout)
!!! This functioin builds up the newton interpolant and evaluates it at xout
!!
!! INPUT: 
!! x: mesh points to be used to build the interpolant.
!! u: divided difference need to build the interpolant.
!! d: interpolant degree.
!! xout: where we wish to evaluate  the interpolant.
!!
!! OUTPUT:
!! yout: result of evaluating the interpolant at xout.

  integer, intent(in)		:: d
  real(kind=8), intent(in)		:: xout, x(d+1), u(d+1)
  real(kind=8), intent(out)		:: yout 
  integer			:: i

  !!-TAJO: Original non-optimized
  yout = u(d+1)
  do i=d,1,-1
    yout = yout * (xout -x(i)) + u(i)
  enddo

end subroutine


subroutine adaptiveinterpolation1D(x, y, n, xout, yout, m, degree, interpolation_type, &
                                         st, eps0, eps1, deg, & 
                                         uumin, uumax, x_table, u_table, &
                                         lambda_table, sigma_table, prod_sigma_table)
!!!This routine adaptively build interpolants that are the evaluated at xout.
!!
!! INPUT: 
!! n: the number points in the 1D vector x.
!! m: the number of points in the 1D vector xout.
!! x: 1D mesh points of length n. For i=1, ..., n-1 x_{i} >  x_{i+1}
!! y: 1D vector that have the data values associated with the points x_{i} for i=1, ..., n
!! xout: 1D vector of length m that represent the locations where we which the interpolate to.
!! limiter: used to determine the type limiter to be used to build interpolant.
!!   - limiter=0: the interpolants are built by adaptively selecting the stencil points. This a essentially 
!!                non-oscillatory (ENO) approach.[https://www.sciencedirect.com/science/article/pii/S0021999196956326]).
!!   - limiter=1: a data-bounded limiter is used to ensure that the interpolant is bounded by the data values.
!!                The data bounded limiter is base on Berzins' work [https://epubs.siam.org/doi/pdf/10.1137/050625667].
!!   - limiter=2: a positivity-preserving limiter is used ensure that interpolant is not bounded by the data values but
!!                remains positive.
!! degree: degree of desired interpolant for each interval.
!!
!! OUTPUT:
!! yout: results of evaluating interpolants at the locations xout.
!! d(3): d(1) minmum order used, d(2) is the average order used, d(3) maximum order used.
!! deg: 1D vector that holds the degree of the interpolant used for each interval

  implicit none

  integer, intent(in)                   :: degree                       !! target polynomial degree for each subinterval 
  integer, intent(in)                   :: interpolation_type                      !! determines the type of interpolation to be used see begining of subroutine
  integer, intent(in)                   :: n                            !! number of input points
  integer, intent(in)                   :: m                            !! number of output points 

  real(kind=8), intent(in)              :: x(n)                       !! input points
  real(kind=8), intent(in)              :: y(n)                       !! data values associated with intput points
  real(kind=8), intent(in)              :: xout(m)                      !! output points
  real(kind=8), intent(out)             :: yout(m)                      !! output values associated with output points 

  integer, intent(out), optional        :: deg(n-1)                     !! to store degree used for each subinterval
  integer, intent(in), optional         :: st 
  real(kind=8), intent(in), optional    :: eps0				!! Used to constrain upper and lowrer bound of interpolant in interval with no hidden extrema        
  real(kind=8), intent(in), optional    :: eps1                         !! Used to constrain upper and lower bounds of interpolant in intervals with hidden extrema.
  real(kind=8), intent(inout), optional :: uumin(n-1), uumax(n-1)
  real(kind=8), intent(out), optional   :: u_table(n-1, degree+1)           !! table of devided diferences
  real(kind=8), intent(out), optional   :: lambda_table(n-1, degree+1)           !! table of devided diferences
  real(kind=8), intent(out), optional   :: sigma_table(n-1, degree+1)           !! table of devided diferences
  real(kind=8), intent(out), optional   :: prod_sigma_table(n-1, degree+1)           !! table of devided diferences
  real(kind=8), intent(out), optional   :: x_table(n-1, degree+1)           !! table of devided diferences


  !!** Local variables
  real(kind=8)                          :: u(degree+1)                  !! to store divided differences associated with selected points in xval
  real(kind=8)                          :: e(degree+1)        
  real(kind=8)                          :: error(degree+1)        
  real(kind=8)                          :: d(degree+1)        
  real(kind=8)                          :: up_b(degree+1)        
  real(kind=8)                          :: low_b(degree+1)        
  real(kind=8)                          :: lambda(degree+1)             !! values to tested againts DBI (limiter=1) and PPI (limiter=2)
  real(kind=8)                          :: sigma(degree+1)             !! values to tested againts DBI (limiter=1) and PPI (limiter=2)

  real(kind=8)                          :: sigma2(degree+1)             !! values to tested againts DBI (limiter=1) and PPI (limiter=2)
  real(kind=8)                          :: prod_sigma(degree+1)             !! values to tested againts DBI (limiter=1) and PPI (limiter=2)
  real(kind=8)                          :: prod_deltax(degree+1)
  real(kind=8)                          :: xval(degree+1)               !! to store slected points in order
  real(kind=8)                          :: table(n, degree+4)           !! table of devided diferences
  real(kind=8)                          :: ur, ul, ww, lambda_l, lambda_r, m_lambda, m_sigma
  real(kind=8)                          :: sigma_l, sigma_r
!!  real(kind=8)                          :: sigma_ll, sigma_lr, sigma_rl, sigma_rr
  real(kind=8)                          :: prod_sigma_l, prod_sigma_r
  real(kind=8)                          :: prod_sigma_l2, prod_sigma_r2
  real(kind=8)                          :: d_l, d_r, up_b_l, up_b_r
  real(kind=8)                          :: low_b_l, low_b_r, m_l, m_r
  real(kind=8)                          :: mm_l(n-1), mm_r(n-1), www(n-1)
  real(kind=8)                          :: prod_deltax_l, prod_deltax_r
  real(kind=8)                          :: a, b
  real(kind=8)                          :: slope(n+1), slope_i, slope_im1, slope_ip1 , tol
  real(kind=8)                          :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, delta
  real(kind=8)                          :: eps, inv_eps, eps2, eps3 
  real(kind=8)                          :: umax, umin                        !! parameter  used to for upper bound for each interval
  real(kind=8)                          :: xl , xr 
  integer                               :: i, j, k, kk!!, jj, ii
  integer                               :: si, ei!!, idxs, idxe
  integer                               :: tmp_si, tmp_ei
  integer                               :: tmp_idx, fid
!!!  integer                               :: umin_umax_type
  integer                               :: stencil_type
  integer, parameter                    :: Debug = 0
  logical                               :: bool_left, bool_right
  


  !!** Check input to make sure x_{i} < x_{i+1} **!!
  do i=1, n-1
    if( x(i) .ge. x(i+1) .or. abs(x(i+1)-x(i)) .le. eps ) then
     write(*,*)'ERROR: Incorrect input at i=', i, 'x(i)=', x(i),&
                 'x(i+1)=',x(i+1), 'x(i) must be less that x(i+1) and &
                |x(i+1)-x(i)| mus be greater than machine precision eps.'
     call exit(0)
    endif
  enddo

  !!** Check output to make sure x_{i} < x_{i+1} **!!
  do i=1, m-1
    if( xout(i) .ge. xout(i+1) ) then
      write(*,*)'ERROR: Incorrect output at k=', i, 'xout(k)=', xout(i),&
                 'xout(k+1)=',xout(i+1), 'xout(k) must be less that xout(k+1)'
      call exit(0)
    endif
  enddo

  !!** Set index for inspection  **!!!
  !! ii =18

  !!** Initialize variables **!!
  !!umin_umax_type = 2            !! used to select how to compute bounds u_{min} and u_{max} 
  !!stencil_type = 5              !! used to choose the stencil selection process
  !!tol = 1.0e-8
  table = 0.0  
  eps = 1e-30                   !! defined as epsilon 
  inv_eps = 1e+30               !! defined to be + infinity
  k = 1                         !! iteration idex used for output points
!!  m_lambda = 100.0              !! constant used to bound two consecutive \bar{lambda}_{j} and \bar{lambda}_{j+1}
  m_sigma = 100.0                !! constant used to bound two consecutive \bar{lambda}_{j} and \bar{lambda}_{j+1}

  !!** Using eps0 to set eps2. eps0 is a user defined parameter used in the PPI
  !!   method to relax the bounds the interpolant for the cases where hiden
  !!   extremum is dectected. eps0 and eps2 should be set to small values
  !!   otherwise this may lead to large oscillation for the PPI algorithm **!!
  if (present(eps0) )then     
    eps2 = eps0
  else
    eps2 = 0.01
  endif


  if(present(eps1)) then
    eps3 = eps1
  else
    eps3 = 1.0
  endif

  if(present(x_table) )then
    x_table = 0.0
  endif

  if(present(u_table) )then
    u_table = 0.0
  endif

  if(present(lambda_table) )then
    lambda_table = 0.0
  endif
  
  if(present(st) )then
    stencil_type  = st
  else
    stencil_type = 1
  endif 
  !!write(*,*) 'TAJO stencil_type =', stencil_type, 'eps0 =', eps0
  

  call divdiff(x, y, n, degree+3, table)    !!compute the table of divided differences 


  !!** Calculate slopes for each interval **!!
  slope(1) = (y(3)-y(2))/(x(3)-x(2))  !! left boundary
  do i=1, n-1
    slope(i+1) = (y(i+1)-y(i))/(x(i+1)-x(i))  !! right boundary
  enddo
  slope(n+1) = slope(n-1)
  

  !!!!** Mark intervals with extrema **!!
  !!extrema(1) = 0
  !!do i=1, n-1
  !!  if(slope(i+1)*slope(i) < 0 .or. slope(i)*slope(i+2) < 0 ) then
  !!    extrema(i+1) = 1   
  !!  else
  !!    extrema(i) = 0
  !!  endif
  !!enddo
  !!extrema(n+1) = 0

  !!** Calculate the polynomial bounds for each interval **!!
  if(degree > 1 .and. interpolation_type .eq. 2) then
  do i=1,n-1
    !!** the slople for interval i **!!
    slope_im1 = slope(i)
    slope_i   = slope(i+1)
    slope_ip1 = slope(i+2)

    umin = min(y(i), y(i+1))
    umax = max(y(i), y(i+1))
    u(2)= (y(i+1)-y(i))/(x(i+1)-x(i))    !! set second slected divided difference

!!!    !!if(degree > 1 ) then
!!!    !!** Option 1: umin_umax_type == 1. Note that if y(i) and y(i+1) are positive 
!!!    !!   u_min = 0.0 and u_max = 2* max(y(i), y(i+1)). **!!
!!!    if(umin_umax_type .eq. 1) then
!!!      tmp1 = min(y(i), y(i+1))
!!!      tmp2 = max(y(i), y(i+1))
!!!      umin = tmp1 - abs(tmp1)
!!!      umax = tmp2 + abs(tmp2)
!!!    endif
!!!  
!!!    !!** Option 2: umin_umax_type == 2. In intervals where a hidden local 
!!!    !!   extremum is detectected  umin and/or umax are Calculation of umin umax option 2 **!!
!!!    if(umin_umax_type .eq. 2) then

      tmp1 = min(y(i), y(i+1))
      tmp2 = max(y(i), y(i+1))

      !!** Calculcualtion of umin based of the existence of an extremum **!!
      if( (slope_im1*slope_ip1 < 0.0 .and. slope_im1 < 0.0) .or. &          !! Detects a minimum
          (slope_im1*slope_ip1 > 0.0 .and. slope_im1*slope_i < 0.0) ) then  !! Detects a maximum and/or minimum (ambiguous).
        umin = tmp1 - eps3*abs(tmp1)
      else                                                                  !! No extremum detected
        umin = tmp1 - eps2*abs(tmp1)
      endif 

      !!** Calculation of umax based on the existence of an extremum **!!
      if( (slope_im1*slope_ip1 < 0.0 .and. slope_im1 > 0.0) .or. &          !! Detects a maximum
          (slope_im1*slope_ip1 > 0.0 .and. slope_im1*slope_i < 0.0) ) then  !! Detects a minimum and/or a maxmimum
        umax = tmp2 + eps3*abs(tmp2)
      else                                                                  !! No extremum is detected
        umax = tmp2 + eps2*abs(tmp2)
      endif

!!!    endif 

    !!!!!** Calculation of umin umax option 3 **!!
    !!!if(umin_umax_type .eq. 3) then
    !!!  tmp1 = min(y(i), y(i+1))
    !!!  tmp2 = max(y(i), y(i+1))

    !!!  !!** Dectection of hidden local extremum **!!
    !!!  if(slope_im1*slope_ip1 < 0.0 .or. slope_i*slope_im1 < 0.0  ) then
    !!!    !!** Caluclate max and/or min of quadratic approximation **!!
    !!!    tmp_si = max(i-1, 1)
    !!!!!    ul = table(tmp_si, 3 )                    
    !!!!!    tmp_ei = min(i+2,n) 
    !!!!!    ur = table(i, min(i, tmp_ei-i+1) )
    !!!!!    xp = 0.5*(x(i) + x(i+1) - u(2) / (ul+eps))
    !!!!!    tmp3 = y(i) + u(2)* (xp-x(i+1)) + ul*(xp-x(i))*(xp-x(i+1))
    !!!!!    xp = 0.5*(x(i) + x(i+1) - u(2) / (ur+eps))
    !!!!!    tmp4 =y(i) + u(2)* (xp-x(i+1)) + ur*(xp-x(i))*(xp-x(i+1))
    !!!!!    
    !!!!!    !!** Evaluate linear interpolants at x_{i} and x_{i+1} **!!
    !!!!!    ul =table(tmp_si, 2)
    !!!    ur =table(tmp_ei-1, 2)
    !!!    tmp5 = y(tmp_si) + ul*(x(i+1) -x(tmp_si)) 
    !!!    tmp6 = y(tmp_ei-1) + ur * (x(i)-x(tmp_ei-1))

    !!!    !!** Update umin umax based on linear and quadratic resulst **!!
    !!!    ul = 2*abs( tmp1 - min(tmp3, tmp4, tmp5, tmp6, (1.0-eps2)*tmp1) )
    !!!    ur = 2*abs( tmp2 - max(tmp3, tmp4, tmp5, tmp6, (1.0-eps2)*tmp2) )
    !!!    umin = tmp1 - ul
    !!!    umax = tmp2 + ur
    !!!  else
    !!!    umin = tmp1 - eps2 * abs(tmp1)
    !!!    umax = tmp2 + eps2 * abs(tmp2)
    !!!  endif
    !!!endif

  
    !!!!! For debugging and investigations
    !!!if (present(uumin) ) then
    !!!  if( present(umin_umax_input) ) then
    !!!    !!** In case uumin and uumax are provided, used use provided values to
    !!!    !!   update umin and umax **!! 
    !!!    if( umin_umax_input) then
    !!!      umin = uumin(i)
    !!!    else
    !!!    !!** Incase uumin and uumax are not provide but the arrays uumin and uumax
    !!!    !!   are passed update the arrays using umin and umax **!! 
    !!!      uumin(i) = umin
    !!!    endif
    !!!  endif
    !!!endif
    !!!if (present(uumax) ) then
    !!!  if( present(umin_umax_input) ) then
    !!!    !!** In case uumin and uumax are provided, used use provided values to
    !!!    !!   update umin and umax **!! 
    !!!    if( umin_umax_input ) then
    !!!      umax = uumax(i)
    !!!    else
    !!!    !!** Incase uumin and uumax are not provide but the arrays uumin and uumax
    !!!    !!   are passed update the arrays using umin and umax **!! 
    !!!      uumax(i) = umax
    !!!    endif
    !!!  endif
    !!!endif

    !!** Compute the values of m_{\ell} and m_r for the positive-preserving 
    !!   method. This coresponds to the default setting **!!
    if(y(i) < y(i+1)) then
      ww = u(2)
      m_l = (umin-y(i)) / (y(i+1)-y(i))  
      m_l = min(0.0, m_l)
      m_r = (umax-y(i)) /  (y(i+1)-y(i)) 
      m_r = max(1.0, m_r)
    elseif(y(i) > y(i+1)) then
      ww = u(2)
      m_l = (umax-y(i)) /  (y(i+1)-y(i)) 
      m_l = min(0.0, m_l)
      m_r = (umin-y(i)) / (y(i+1)-y(i))  
      m_r = max(1.0, m_r)
    else !! This part deals with the special case where y(i)=y(i+1) 
      ww = u(2)
      tmp_si = max(i-1, 1)
      ul = table(tmp_si, 3 )    !! divided difference U[x_{i-1}, x_{i}, x_{i+1}] 
      tmp_ei = min(i+2,n) 
      ur = table(i, min(i, tmp_ei-i+1) )  !! divided difference U[x_{i}, x_{i+1}, x_{i+2}] 
      if( ul >  0.0 .and. ur .ne. 0.0)then
        ww = ul * (x(i+1)-x(i)) * (x(i+1)-x(tmp_si)) 
        !!u(2) = ww
        m_l = (umin-y(i)) / ww  
        m_l = min(0.0, m_l)
        m_r = (umax-y(i)) / ww 
        m_r = max(1.0, m_r)
      else if (ul < 0.0  .and. ur .ne. 0.0) then
        ww = ul * (x(i+1)-x(i)) * (x(i+1)-x(tmp_si)) 
        !!u(2) = ww
        m_l = (umax-y(i)) /  ww  
        m_l = min(0.0, m_l)
        m_r = (umin-y(i)) / ww  
        m_r = max(1.0, m_r)
      else if (ur > 0.0 .and. ul .ne. 0.0) then
        ww = ur * (x(i+1)-x(i)) * (x(tmp_ei)-x(i)) 
        !!u(2) = ww
        m_l = (umin-y(i)) / ww  
        m_l = min(0.0, m_l)
        m_r = (umax-y(i)) / ww 
        m_r = max(1.0, m_r)
      else if (ur < 0.0 .and. ul .ne. 0.0) then
        ww = ur * (x(i+1)-x(i)) * (x(tmp_ei)-x(i)) 
        !!u(2) = ww
        m_l = (umax-y(i)) /  ww  
        m_l = min(0.0, m_l)
        m_r = (umin-y(i)) / ww  
        m_r = max(1.0, m_r)
      else !! ur=0 or ul=0, the algorithm defaults to BDI.  
        ww = u(2)
        m_l = 0.0
        m_r = 1.0
      endif
    endif

    !!** Save the bunds for each interval and the www ***!!
    mm_l(i) = m_l
    mm_r(i) = m_r
    www(i) = ww
  enddo
  !!** Default case: DBI **!!
  else  
    do i=1, n-1
      !!** Compute the values of m_{\ell} and m_r for the data-bounded method 
      !!   if the limiter variable is set to 1 **!!
      mm_l(i) = 0.0
      mm_r(i) = 1.0
      www(i) = (y(i+1)-y(i)) / (x(i+1)-x(i)) 
    enddo
  endif

  !!write(*,*) ' TAJO ml=', mm_l(n-1), 'mr=', mm_r(n-1)
  !!** loop over each input intervals. For each  interval build an interpolant and 
  !!   evaluate the interpolant at the desired output points **!!
  do i=1,n-1
 
    !!** Initialize varibles for each interval**!!
    u = 0.0
    xval = 0.0
    up_b = 0.0
    low_b = 0.0
    prod_deltax = 1.0
    lambda = 0.0
    sigma = 1.0

    xval(1) = x(i)                       !! first point in stencil
    xval(2) = x(i+1)                     !! second point in stencil
    u(1) = y(i)                          !! set first selected divided difference
    u(2)= (y(i+1)-y(i))/(x(i+1)-x(i))    !! set second slected divided difference
    lambda(1) = 1.0                      !! set first ratio of divided difference  
    lambda(2) = 1.0                      !! set second ratio of divided difference  
    error(1) = u(2)*(x(i+1)-x(i))
    sigma(1) = 1.0
    sigma(2) = u(2)*(x(i+1)-x(i))
    sigma2(1) = 1.0
    sigma2(2) = 1.0
    prod_sigma(1) = 1.0
    prod_sigma(2) = (x(i+1)-x(i))
    up_b(1) = 1.0
    up_b(2) = 1.0
    low_b(1) = -1.0
    low_b(2) = -1.0
    m_l = mm_l(i)
    m_r = mm_r(i)
    ww = www(i)
    ei=i+1                               !! set e before entering the do loop
    si=i                                 !! set s before entering the do loop
   
    !!** Continue to build high degree interpolant only if the target degree is
    !!   greater than 1 and the both U{x_{i}, x_{i+1}, x_{i+2}} and U{x_{i-1},x_{i}, x_{i+1}}
    !!   are not zeros. **!!
    if(degree > 1 .and. ww .ne. 0.0)then
      do j=2, degree 

         !!** Initialize selection boolean variables **!!
         bool_left = .false.
         bool_right = .false.

         tmp_si = max(1, si-1)                          !! decrementing stencil left idex
         tmp_ei = min(n, ei+1)                          !! incrementing stencil right index
         
         !!** Calculate ul and ur (left and right divided differences
         !!   respectively)
         !!   Calculate lambda_l and lambda_r ( left and right ratio of 
         !!   of divided references respectively ) **!!
         prod_deltax_l = prod_deltax(j) * (x(ei)-x(tmp_si)) !! product of interval containing stencil
         if(xval(j) < x(i))then
           prod_sigma_l = prod_sigma(j)*(x(i+1)-xval(j))
         else
           prod_sigma_l = prod_sigma(j)*(x(i)-xval(j))
         endif
         prod_sigma_l2 = prod_sigma_l*(x(i+1)-x(tmp_si))
         prod_deltax_l = prod_deltax(j) * (x(ei)-x(tmp_si)) !! product of interval containing stencil
         if(si-1 > 0)then
           ul = table(tmp_si, ei-tmp_si+1)   !! get left divided difference   
           lambda_l = ul/ww * prod_deltax_l  !! calculate left lambda         
           sigma_l = ul * prod_sigma_l
           !!sigma_l =  prod_sigma_l
           xl = x(si-1)
         else
           ul = inv_eps  !! set left dividided difference to + infinity  
           lambda_l = inv_eps    !! set left lambda to infinity
           sigma_l = inv_eps
           xl = -inv_eps
         endif

         prod_deltax_r = prod_deltax(j) * (x(tmp_ei)-x(si))    !! product of inteval containing stencil
         if(xval(j) < x(i) )then
           prod_sigma_r = prod_sigma(j)*(x(i+1)-xval(j))
         else
           prod_sigma_r = prod_sigma(j)*(x(i)-xval(j))
         endif
     
         prod_sigma_r2 = prod_sigma_r*(x(tmp_ei)-x(i))
         if(ei+1 .le. n)then
           ur = table(si, tmp_ei-si+1)    !! get right divided difference 
           lambda_r =  ur/ ww * prod_deltax_r  !! calculate righ lambda
           sigma_r = ur * prod_sigma_r
           !!sigma_r = prod_sigma_r
           xr = x(ei+1)
         else
           ur = inv_eps    !! set righl divided difference to infinity
           lambda_r =  inv_eps   !! set righ lambda to infinity
           sigma_r = inv_eps
           xr = inv_eps
         endif

         e(j) = -(x(i)-xval(j))/(x(i+1)-x(i))  !! calculate e_j !!
         d_l = (x(ei)-x(tmp_si))/(x(i+1)-x(i)) !! calculate d_l 
         d_r = (x(tmp_ei)-x(si))/(x(i+1)-x(i)) !! calculate d_r

         !!** In the case where the points inserted to 
         !!   V_{j-1} form V_{j} is to the left. Calculate 
         !!   upper and lower bounds up_b_l and low_bl 
         !!   respectively  **!!
         if(j .eq. 2 )then 
           up_b_l = d_l*( -m_l*4.0 + 1.0 )
           low_b_l = d_l*( -(m_r-1.0)*4.0 - 1.0 ) 
         else
           if(e(j) <= 0.0) then
             up_b_l = (up_b(j) - lambda(j))* d_l / (1.0-e(j))
             low_b_l = (low_b(j) - lambda(j))*d_l / (1.0-e(j))
           elseif(e(j) > 0.0) then
             up_b_l = (low_b(j) - lambda(j))*d_l / (0.0-e(j))
             low_b_l = (up_b(j) - lambda(j))* d_l / (0.0-e(j))
           endif

         endif
     
         !!** In the case where the points inserted to 
         !!   V_{j-1} form V_{j} is to the right. Calculate 
         !!   upper and lower bounds up_b_r and low_b_r 
         !!   respectively  **!!
         if(j .eq. 2)then
            up_b_r = d_r*( -m_l*4.0 + 1.0 )
            low_b_r = d_r*( -(m_r-1.0)*4.0 - 1.0 )  
         else
           if(e(j) <= 0.0) then
             up_b_r = (up_b(j) - lambda(j))* d_r / (1.0-e(j))
             low_b_r = (low_b(j) - lambda(j))*d_r / (1.0-e(j))
           elseif(e(j) > 0.0) then
             up_b_r = (low_b(j) - lambda(j))*d_r / (0.0-e(j))
             low_b_r = (up_b(j) - lambda(j))* d_r / (0.0-e(j))
           endif

         endif
 	

         !!** Option 1: stencil_type = 1. In addition to positivity or 
         !!   data boundedness, the stencil selection is based on the ENO approach **!!
         if(stencil_type .eq. 1) then
           if( (low_b_l .le. lambda_l .and. lambda_l .le. up_b_l) .and. & !! Adding a point to left meets the requiremenst for DBI or PPI
               (low_b_r .le. lambda_r .and. lambda_r .le. up_b_r) )then   !! Adding a point to right meets the requiremenst for DBI or PPI
             !!** boolean variable is set to true based on the coresponding 
             !!   divided difference |ul| or |ur| is the smalest
             if(abs(ul) < abs(ur) )then
               bool_left = .true.
               bool_right = .false.
             else
               bool_left = .false.
               bool_right = .true.
             endif
           else if(low_b_r .le. lambda_r .and. lambda_r .le. up_b_r) then !! Adding a point to right meets the requiremenst for DBI or PPI
             bool_left = .false.
             bool_right = .true.
           else if(low_b_l .le. lambda_l .and. lambda_l .le. up_b_l) then !! Adding a point to left meets the requiremenst for DBI or PPI
             bool_left = .true.
             bool_right = .false.
           endif
         endif 

         !! Option 2: stencil_type = 2. In addition to DBI or PPI the 
         !! stencil selection prioritize a symetric stencil other others **!!
         if(stencil_type .eq. 2) then
           if( (low_b_l .le. lambda_l .and. lambda_l .le. up_b_l) .and. &
               (low_b_r .le. lambda_r .and. lambda_r .le. up_b_r) )then
             if( i-si < ei-i )then
               bool_left = .true.
               bool_right = .false.
             elseif( i-si > ei-i ) then
               bool_left = .false.
               bool_right = .true.
             else
               if(abs(lambda_l) < abs(lambda_r) )then
               !!if(abs(sigma_l) < abs(sigma_r) )then
                 bool_left = .true.
                 bool_right = .false.
               else
                 bool_left = .false.
                 bool_right = .true.
               endif
             endif
           else if(low_b_r .le. lambda_r .and. lambda_r .le. up_b_r ) then 
             bool_left = .false.
             bool_right = .true.
           else if(low_b_l .le. lambda_l .and. lambda_l .le. up_b_l) then 
             bool_left = .true.
             bool_right = .false.
           endif
         endif 

         !!** Stencil choice option 3 **!!
         if(stencil_type .eq. 3) then
           if( (low_b_l .le. lambda_l .and. lambda_l .le. up_b_l) .and. &
               (low_b_r .le. lambda_r .and. lambda_r .le. up_b_r) )then
             if( abs(x(i)-xl) < abs(xr-x(i+1)) )then
               bool_left = .true.
               bool_right = .false.
             elseif( abs(x(i)-xl) > abs(xr-x(i+1)) )then
               bool_left = .false.
               bool_right = .true.
             else
               if(abs(lambda_l) < abs(lambda_r) )then
                 bool_left = .true.
                 bool_right = .false.
               else
                 bool_left = .false.
                 bool_right = .true.
               endif
             endif
           else if(low_b_r .le. lambda_r .and. lambda_r .le. up_b_r ) then
             bool_left = .false.
             bool_right = .true.
           else if(low_b_l .le. lambda_l .and. lambda_l .le. up_b_l) then
             bool_left = .true.
             bool_right = .false.
           endif

         endif


         !!** Option 2: stencil_type = 2. In adddition to DBI or PPI the 
         !!   stencil selection is based on the smalest lambda **!!
         if(stencil_type .eq. 4) then
           if( (low_b_l .le. lambda_l .and. lambda_l .le. up_b_l) .and. &
               (low_b_r .le. lambda_r .and. lambda_r .le. up_b_r) )then
             if(abs(lambda_l) < abs(lambda_r) )then
               bool_left = .true.
               bool_right = .false.
             else
               bool_left = .false.
               bool_right = .true.
             endif
           else if(low_b_r .le. lambda_r .and. lambda_r .le. up_b_r) then 
             bool_left = .false.
             bool_right = .true.
           else if(low_b_l .le. lambda_l .and. lambda_l .le. up_b_l) then 
             bool_left = .true.
             bool_right = .false.
           endif
         endif 

     
         !!** Stencil choice option 4 **!!
         if(stencil_type .eq. 5) then
           if( (low_b_l .le. lambda_l .and. lambda_l .le. up_b_l) .and. &
               (low_b_r .le. lambda_r .and. lambda_r .le. up_b_r) )then
             if( abs(sigma_l) < abs(sigma_r) )then
               bool_left = .true.
               bool_right = .false.
             elseif( abs(sigma_l) >= abs(sigma_r) ) then
               bool_left = .false.
               bool_right = .true.
             else
               if(abs(lambda_l) < abs(lambda_r) )then
                 bool_left = .true.
                 bool_right = .false.
               else
                 bool_left = .false.
                 bool_right = .true.
               endif
             endif
           else if(low_b_r .le. lambda_r .and. lambda_r .le. up_b_r ) then 
             bool_left = .false.
             bool_right = .true.
           else if(low_b_l .le. lambda_l .and. lambda_l .le. up_b_l) then 
             bool_left = .true.
             bool_right = .false.
           endif

         endif


         if(stencil_type .eq. 6) then
           if( (low_b_l .le. lambda_l .and. lambda_l .le. up_b_l) .and. &
               (low_b_r .le. lambda_r .and. lambda_r .le. up_b_r) )then
             if( abs(prod_sigma_l2) < abs(prod_sigma_r2) )then
               if(abs(sigma_l) < 100*abs(sigma(j)) ) then
                 bool_left = .true.
                 bool_right = .false.
               elseif( abs(sigma_r) < 100*abs(sigma(j)) )then
                 bool_left = .false.
                 bool_right = .true.
               else
                 bool_left = .true.
                 bool_right = .false.
               endif
             elseif( abs(prod_sigma_l2) > abs(prod_sigma_r2) ) then
               if( abs(sigma_r) < 100*abs(sigma(j)) )then
                 bool_left = .false.
                 bool_right = .true.
               elseif(abs(sigma_l) < 100*abs(sigma(j)) ) then
                 bool_left = .true.
                 bool_right = .false.
               else
                 bool_left = .false.
                 bool_right = .true.
               endif
             else
               if(abs(sigma_l) < abs(sigma_r) )then
                 bool_left = .true.
                 bool_right = .false.
               else
                 bool_left = .false.
                 bool_right = .true.
               endif
             endif
           else if(low_b_r .le. lambda_r .and. lambda_r .le. up_b_r ) then 
             bool_left = .false.
             bool_right = .true.
           else if(low_b_l .le. lambda_l .and. lambda_l .le. up_b_l) then 
             bool_left = .true.
             bool_right = .false.
           endif

         endif


         !!** Add point to the left of current stencil and corresponding 
         !!   variables **!!
         if( (bool_left .eqv. .true.) .and. (bool_right .eqv. .false.)) then
           si = max(1, si-1)
           !!ei = ei
           lambda(j+1) = lambda_l
           u(j+1) = ul
           xval(j+1) = x(si)
           up_b(j+1) = up_b_l
           low_b(j+1) = low_b_l
           prod_deltax(j+1) = prod_deltax_l
           sigma(j+1) = sigma_l
           prod_sigma(j+1) = prod_sigma_l
           if(abs(u(j)) < eps )then
             sigma2(j+1) = inv_eps;
           else
             if(xval(j) > x(i)) then
               sigma2(j+1) = 1 + u(j+1)/u(j) * (xval(j) -x(i))
             else
               sigma2(j+1) = 1 + u(j+1)/u(j) * (xval(j) -x(i+1))
             endif
           endif
         !!** Add point to the right of current stencil and corresponding 
         !!   variables **!!
         elseif( (bool_left .eqv. .false.) .and. (bool_right .eqv. .true.) ) then
           !!si = si
           ei = min(ei+1, n)
           lambda(j+1) = lambda_r
           u(j+1) = ur
           xval(j+1) = x(ei)
           up_b(j+1) = up_b_r
           low_b(j+1) = low_b_r
           prod_deltax(j+1) = prod_deltax_r
           sigma(j+1) = sigma_r
           prod_sigma(j+1) = prod_sigma_r
           if(abs(u(j)) < eps )then
             sigma2(j+1) = inv_eps;
           else
             if(xval(j) > x(i)) then
               sigma2(j+1) = 1 + u(j+1)/u(j) * (xval(j) -x(i))
             else
               sigma2(j+1) = 1 + u(j+1)/u(j) * (xval(j) -x(i+1))
             endif
           endif
         !!*
         endif
         

         !!** write dat for investigations **!!
         if( i .eq. 24 .and.  m .eq. 62 ) then
           write(*,*) 'bool_left=', bool_left, 'bool_rigt =', bool_right
           write(*,*) 'i=', i, 'j =', j
           write(*,*) 'si =', si
           write(*,*) 'ei =', ei
           write(*,*) 'e =', e(j)
           write(*,*) 'low_b_r ', low_b_r
           write(*,*) 'lambda_r ', lambda_r
           write(*,*) 'up_b_r ', up_b_r
           write(*,*) 'ur ', ur
 
           write(*,*) 'low_b_l ', low_b_l
           write(*,*) 'lambda_l ', lambda_l
           write(*,*) 'up_b_l ', up_b_l
           write(*,*) 'ul ', ul
           write(*,*) 'u(j+1) ', u(j+1)
         endif
      enddo !! j loop
    endif !! end of if j > 1


    !!** save the interpolant degree used for the interval [x_{i}, x_{i+1}] **!!
    deg(i) = ei-si


    !!!!!!** **!!
    !!if(degree > 4)then
    !!do j=4, degree
    !!  if(abs(sigma(j-2)) < abs(sigma(j-1)) .and. abs(sigma(j-1)) < abs(sigma(j)) .and. abs(sigma(j)) < abs(sigma(j+1)) )then
    !!  !!if(abs(sigma(j-1)) < abs(sigma(j)) .and. abs(sigma(j)) < abs(sigma(j+1)) )then
    !!  !!if(abs(sigma(j)) > 10.0*abs(sigma(j+1)) ) then
    !!    u(j-2:degree+1) = 0.0
    !!    deg(i) = j-2
    !!    exit
    !!  endif
    !!enddo
    !!endif

    !!** Extrapolate to points that are to the left of the defined interval **!! 
    if(k <=m)then
      do while( xout(k) < x(1) )
        write(*,*) 'WARNING: Some of the output are obtained via &
                    extrapolation instead of interpolation. The desired &
                    proprety such as data-boundedness or positvity is not &
                    preserved in such case'
          write(*,*)  k, 1 
          write(*,*)  xout(k), x(1) 
        call newtonPolyVal(xval, u, degree, xout(k), yout(k))
        k = k+1
        if(k > m) exit
      enddo
    endif
 
    !!** Building and evaluating Interpolant at xout points **!!
    !! - do while( x(i) <= xout(k) .and. xout(k) <= x(i+1) .and. k <= m )
    if( k <=m)then
      do while( x(i) <= xout(k) .and. xout(k) <= x(i+1) )
        call newtonPolyVal(xval, u, degree, xout(k), yout(k))
        k = k+1
        if(k > m) exit
      enddo
    endif

    !!** Extrapolate to points that are to the right of the defined interval **!! 
    if(k <= m)then
      do while( xout(k) > x(n) )
          write(*,*) 'WARNING: Some of the output are obtained via &
                      extrapolation instead of interpolation. The desired &
                      proprety such as data-boundedness or positvity is not &
                      preserved in such case'
          write(*,*)  k, n 
          write(*,*)  xout(k), x(n) 
          call newtonPolyVal(xval, u, degree, xout(k), yout(k))
          k = k+1
          if(k > m) exit
      enddo
    endif

    if(n == 17 .and. degree ==4) then
       print*, si, '  ---- ',ei 
    endif
    !!** Save data into table for study **!!
    if (present(x_table) ) then     
    do kk=1, degree+1
      x_table(i, kk) = xval(kk)
    enddo
    endif
    if (present(u_table) ) then     
    do kk=1, degree+1
      u_table(i, kk) = u(kk)
    enddo
    endif
    if (present(lambda_table) ) then     
    do kk=1, degree+1
      !!lambda_table(i, kk) = lambda(kk)
      lambda_table(i, kk) = sigma2(kk)
    enddo
    endif
    if (present(sigma_table) ) then     
    do kk=1, degree+1
      sigma_table(i, kk) = sigma(kk)
    enddo
    endif
    if (present(prod_sigma_table) ) then     
    do kk=1, degree+1
      prod_sigma_table(i, kk) = prod_sigma(kk)
    enddo
    endif

  end do !! end of i loop



end subroutine adaptiveInterpolation1D


subroutine adaptiveInterpolation2D(x, y, nx, ny, v,  xout, yout, mx, my, vout, d, degree, limiter)
!!!This routine adaptively build interpolants that are the evaluated at xout, yout
!!
!! INPUT: 
!! x: 1D array of x coodinates
!! y: 1D array of y coodinates
!! v: 2D array. v is the data values associated with the locations (x(i), y(j))
!! xout: 2D array. Locations where we wihch the interpolate to.
!! yout: 2D array. Locations where we wihch the interpolate to.
!! limiter: used to determine the type limiter to be used to build interpolant.
!!   - limiter=0: the interpolants are built by adaptively selecting the stencil points. This a essentially 
!!                non-oscillatory (ENO) approach.[https://www.sciencedirect.com/science/article/pii/S0021999196956326]).
!!   - limiter=1: a data-bounded limiter is used to ensure that the interpolant is bounded by the data values.
!!                Berzins work [https://epubs.siam.org/doi/pdf/10.1137/050625667].
!!   - limiter=2: a positivity-preserving limiter is used ensure that interpolant is not bounded by the data values but
!!                remains positive.
!! nx: the number points in the x direction.
!! ny: the number points in the y direction.
!! mx: the number of points in the x direction where we want to evaluate the interpolant.
!! my: the number of points in the x direction where we want to evaluate the interpolant.
!! degree: degree of desired interpolant for each interval.
!!
!! OUTPUT:
!! vout: 2D array. results of evaluating interpolants at the locations (xout(i) yout(j)).
!! d(2,3): NOT FINISHED

  
  integer, intent(in)                   :: nx, ny, mx, my, degree, limiter
  real(kind=8), intent(in)              :: x(nx), y(ny), v(nx,ny), xout(mx), yout(my)
  real(kind=8), intent(out)             :: vout(mx, my), d(2,3)

  integer                               :: i, j
  real(kind=8)                          :: voutx(mx, ny), dx(3), dy(3)
  integer                               :: degx(nx), degy(ny)

  integer                               :: st
  real(kind=8)                          :: eps0, eps1 

  eps0 = 0.01
  eps1 = 1.0
  st = 1
  !!** initialize variable **!!
  voutx = 0.0
  dx = 0.0
  dy = 0.0
  degx = 0.0
  degy = 0.0

  !!** interpolate along x **!!
  do j=1,ny
    call adaptiveInterpolation1D(x, v(:,j), nx, xout, voutx(:,j), mx, degree, limiter, st, eps0, eps1, degx)
  enddo

  !!** interpolate along y **!!
  do i=1,mx
    call adaptiveInterpolation1D(y, voutx(i,:), ny, yout, vout(i,:), my, degree, limiter, st, eps0, eps1, degy)
  enddo

  !!** save order used for interpolations **!! 
  d(1,:) = dx(:)
  d(2,:) = dy(:)

end subroutine

!!--subroutine adaptiveInterpolation2DPCHIP(x, y, nx, ny, v,  xout, yout, mx, my, vout, d, degree, limiter)
!!--!!!This routine adaptively build interpolants that are the evaluated at xout, yout
!!--!!
!!--!! INPUT: 
!!--!! x: 1D array of x coodinates
!!--!! y: 1D array of y coodinates
!!--!! v: 2D array. v is the data values associated with the locations (x(i), y(j))
!!--!! xout: 2D array. Locations where we wihch the interpolate to.
!!--!! yout: 2D array. Locations where we wihch the interpolate to.
!!--!! limiter: used to determine the type limiter to be used to build interpolant.
!!--!!   - limiter=0: the interpolants are built by adaptively selecting the stencil points. This a essentially 
!!--!!                non-oscillatory (ENO) approach.[https://www.sciencedirect.com/science/article/pii/S0021999196956326]).
!!--!!   - limiter=1: a data-bounded limiter is used to ensure that the interpolant is bounded by the data values.
!!--!!                Berzins work [https://epubs.siam.org/doi/pdf/10.1137/050625667].
!!--!!   - limiter=2: a positivity-preserving limiter is used ensure that interpolant is not bounded by the data values but
!!--!!                remains positive.
!!--!! nx: the number points in the x direction.
!!--!! ny: the number points in the y direction.
!!--!! mx: the number of points in the x direction where we want to evaluate the interpolant.
!!--!! my: the number of points in the x direction where we want to evaluate the interpolant.
!!--!! degree: degree of desired interpolant for each interval.
!!--!!
!!--!! OUTPUT:
!!--!! vout: 2D array. results of evaluating interpolants at the locations (xout(i) yout(j)).
!!--!! d(2,3): NOT FINISHED
!!--
!!--  integer, intent(in)                   :: nx, ny, mx, my, degree, limiter
!!--  real(kind=8), intent(in)              :: x(nx), y(ny), v(nx,ny), xout(mx), yout(my)
!!--  real(kind=8), intent(out)             :: vout(mx, my), d(2,3)
!!--
!!--  integer                               :: i, j
!!--  real(kind=8)                          :: voutx(mx, ny), dx(3), dy(3)
!!--  real(kind=8)                          :: degx(nx), degy(ny)
!!--
!!--  integer                               :: nwk, ierr
!!--  real(kind=8)                          :: wk(nx*2), d_tmp(nx)
!!--  real(kind=8)                          :: fdl(mx)
!!--  character*12                          :: filename
!!--  logical                               :: spline
!!--
!!--  nwk = nx*2
!!--  spline = .false.
!!--
!!--  !!** interpolate along x **!!
!!--  do j=1,ny
!!--    !call adaptiveInterpolation1D(x, v(:,j), nx, xout, voutx(:,j), mx, dx, degree, limiter, degx)
!!--    call pchez(nx, x, v(:, j), d_tmp, spline, wk, nwk, ierr)
!!--    call pchev(nx, x, v(:,j), d_tmp, mx, xout, voutx(:,j), fdl, ierr)
!!--  enddo
!!--
!!--  !!** interpolate along y **!!
!!--  do i=1,mx
!!--    !call adaptiveInterpolation1D(y, voutx(i,:), ny, yout, vout(i,:), my, dy, degree, limiter, degy)
!!--    call pchez(ny, y, voutx(i,:), d_tmp, spline, wk, nwk, ierr)
!!--    call pchev(ny, y, voutx(i,:), d_tmp, my, yout, vout(i,:), fdl, ierr)
!!--  enddo
!!--
!!--  !!** save order used for interpolations **!! 
!!--  d(1,:) = dx(:)
!!--  d(2,:) = dy(:)
!!--
!!--end subroutine

subroutine adaptiveInterpolation3D(x, y, z, nx, ny, nz, v,  xout, yout, zout, mx, my, mz, vout, d, degree, limiter)
!!!This routine adaptively build interpolants that are the evaluated at xout, yout
!!
!! INPUT: 
!! x: 1D array of x coodinates
!! y: 1D array of y coodinates
!! z: 1D array of z coodinates
!! v: 3D array. v is the data values associated with the locations (x(i), y(j), z(k))
!! xout: 3D array. Locations where we wihch the interpolate to.
!! yout: 3D array. Locations where we wihch the interpolate to.
!! zout: 3D array. Locations where we wihch the interpolate to.
!! limiter: used to determine the type limiter to be used to build interpolant.
!!   - limiter=0: the interpolants are built by adaptively selecting the stencil points. This a essentially 
!!                non-oscillatory (ENO) approach.[https://www.sciencedirect.com/science/article/pii/S0021999196956326]).
!!   - limiter=1: a data-bounded limiter is used to ensure that the interpolant is bounded by the data values.
!!                Berzins work [https://epubs.siam.org/doi/pdf/10.1137/050625667].
!!   - limiter=2: a positivity-preserving limiter is used ensure that interpolant is not bounded by the data values but
!!                remains positive.
!! nx: the number points in the x direction.
!! ny: the number points in the y direction.
!! nz: the number points in the y direction.
!! mx: the number of points in the x direction where we want to evaluate the interpolant.
!! my: the number of points in the x direction where we want to evaluate the interpolant.
!! mz: the number of points in the x direction where we want to evaluate the interpolant.
!! degree: degree of desired interpolant for each interval.
!!
!! OUTPUT:
!! vout: 3D array. results of evaluating interpolants at the locations (xout(i) yout(j), zout(k)).
!! d(3,3)

  integer, intent(in)                   :: nx, ny, nz, mx, my, mz, degree, limiter
  real(kind=8), intent(in)              :: x(nx), y(ny), z(nz), v(nx,ny,nz), xout(mx), yout(my), zout(mz)
  real(kind=8), intent(out)             :: vout(mx,my,mz), d(3,3)

  integer                               :: i, j, k, ii, jj, kk
  real(kind=8)                          :: dx(3), dy(3), dz(3)
  integer                               :: degx(nx), degy(ny), degz(nz)
  real(kind=8)                          :: tmpin(max(nx, mx, ny, my, nz,mz))
  real(kind=8)                          :: tmpout(max(nx, mx, ny, my, nz,mz))
  real(kind=8)                          :: voutt(max(nx, mx), max(ny, my),max(nz,mz))
  real(kind=8)				:: eps0, eps1
  integer				:: st

  st=1
  eps0 = 0.01
  eps1 = 1.0

  !voutt = 0.0
  !!** interpolate along x **!!
  do k=1,nz
    do j=1,ny
      do ii=1, nx
        tmpin(ii) = v(ii,j,k)
      enddo
      call adaptiveInterpolation1D(x, tmpin(1:nx), nx, xout, tmpout(1:mx), mx, degree, limiter, st, eps0, eps1, degx)
      do ii=1, mx
        voutt(ii, j,k) = tmpout(ii)
      enddo
    enddo
  enddo

  !!** interpolate along y **!!
  do k=1,nz
    do i=1,mx
      do jj=1, ny
        tmpin(jj) = voutt(i, jj, k) 
      enddo
      call adaptiveInterpolation1D(y, tmpin(1:ny), ny, yout, tmpout(1:my), my, degree, limiter, st, eps0, eps1, degy)
      do jj=1, my
        voutt(i, jj, k) = tmpout(jj)
      enddo
    enddo
  enddo

  !!** interpolate along z **!!
  do j=1,my
    do i=1,mx
      do kk=1, nz
        tmpin(kk) = voutt(i,j, kk)
      enddo
      call adaptiveInterpolation1D(z, tmpin(1:nz), nz, zout, tmpout(1:mz), mz, degree, limiter, st, eps0, eps1, degz)
      do kk=1, mz
        vout(i, j, kk) = tmpout(kk)
      enddo
    enddo
  enddo

  !!** save order used for interpolations **!! 
  d(1,:) = dx(:)
  d(2,:) = dy(:)
  d(2,:) = dz(:)



end subroutine


subroutine adaptiveInterpolation3DPCHIP(x, y, z, nx, ny, nz, v,  xout, yout, zout, mx, my, mz, vout, d, degree, limiter)
!!!This routine adaptively build interpolants that are the evaluated at xout, yout
!!
!! INPUT: 
!! x: 1D array of x coodinates
!! y: 1D array of y coodinates
!! z: 1D array of z coodinates
!! v: 3D array. v is the data values associated with the locations (x(i), y(j), z(k))
!! xout: 3D array. Locations where we wihch the interpolate to.
!! yout: 3D array. Locations where we wihch the interpolate to.
!! zout: 3D array. Locations where we wihch the interpolate to.
!! limiter: used to determine the type limiter to be used to build interpolant.
!!   - limiter=0: the interpolants are built by adaptively selecting the stencil points. This a essentially 
!!                non-oscillatory (ENO) approach.[https://www.sciencedirect.com/science/article/pii/S0021999196956326]).
!!   - limiter=1: a data-bounded limiter is used to ensure that the interpolant is bounded by the data values.
!!                Berzins work [https://epubs.siam.org/doi/pdf/10.1137/050625667].
!!   - limiter=2: a positivity-preserving limiter is used ensure that interpolant is not bounded by the data values but
!!                remains positive.
!! nx: the number points in the x direction.
!! ny: the number points in the y direction.
!! nz: the number points in the y direction.
!! mx: the number of points in the x direction where we want to evaluate the interpolant.
!! my: the number of points in the x direction where we want to evaluate the interpolant.
!! mz: the number of points in the x direction where we want to evaluate the interpolant.
!! degree: degree of desired interpolant for each interval.
!!
!! OUTPUT:
!! vout: 3D array. results of evaluating interpolants at the locations (xout(i) yout(j), zout(k)).
!! d(3,3)

  integer, intent(in)                   :: nx, ny, nz, mx, my, mz, degree, limiter
  real(kind=8), intent(in)              :: x(nx), y(ny), z(nz), v(nx,ny,nz), xout(mx), yout(my), zout(mz)
  real(kind=8), intent(out)             :: vout(mx,my,mz), d(3,3)

  integer                               :: i, j, k
  real(kind=8)                          :: dx(3), dy(3), dz(3)
  real(kind=8)                          :: voutt(max(mx,nx), max(my,ny), max(mz,nz) )
  real(kind=8)                          :: degx(nx), degy(ny), degz(nz)
  real(kind=8)                          :: tmpin(max(nx, mx, ny, my, nz, mz))
  real(kind=8)                          :: tmpout(max(nx, mx, ny, my, nz, mz))

  integer                               :: nwk, ierr
  real(kind=8)                          :: wk(nx*2), d_tmp(nx)
  real(kind=8)                          :: fdl(mx)
  character*12                          :: filename
  logical                               :: spline

  nwk = nx*2
  spline = .false.


  !!** interpolate along x **!!
  do k=1,nz
    do j=1,ny
      tmpin(1:nx) = v(1:nx, j,k)
      call pchez(nx, x, tmpin(1:nx), d_tmp, spline, wk, nwk, ierr)
      call pchev(nx, x, tmpin(1:nx), d_tmp, mx, xout, tmpout(1:mx), fdl, ierr)
      voutt(1:mx, j, k) = tmpout(1:mx)
    enddo
  enddo

  !!** interpolate along y **!!
  do k=1,nz
    do i=1,mx
      tmpin(1:ny) = voutt(i, 1:ny, k)
      call pchez(ny, y, tmpin(1:ny), d_tmp, spline, wk, nwk, ierr)
      call pchev(ny, y, tmpin(1:ny), d_tmp, my, yout, tmpout(1:my), fdl, ierr)
      voutt(i, 1:my, k) = tmpout(1:my)
    enddo
  enddo

  !!** interpolate along z **!!
  do j=1,my
    do i=1,mx
      tmpin(1:nz) = voutt(i, j, 1:nz)
      call pchez(nz, z, tmpin(1:nz), d_tmp, spline, wk, nwk, ierr)
      call pchev(nz, z, tmpin(1:nz), d_tmp, mz, zout, tmpout(1:mz), fdl, ierr)
      vout(i, j, 1:mz) = tmpout(1:mz)
    enddo
  enddo

  !!** save order used for interpolations **!! 
  d(1,:) = dx(:)
  d(2,:) = dy(:)
  d(2,:) = dz(:)



end subroutine


subroutine lagrange_poly(x, y, xout, s, e, yout)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! contruct and evaluate lagrange polynomial at xout
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer, intent(in)                   :: s, e
  real(kind=8), intent(in)              :: x(e-s+1), y(e-s+1), xout
  real(kind=8), intent(out)             :: yout
  integer                               :: n, i, j
  real(kind=8)                          :: tmp1, tmp2

  n = e-s+1
  tmp2 = 0.0
  do i=1, n
    tmp1 = 1.0
    do j=1, n
      if(i /= j) then
       tmp1 = tmp1 * (xout - x(j)) / (x(i)-x(j))
      endif 
    enddo
    tmp2 = tmp2 + y(i) * tmp1
  enddo
  yout = tmp2
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  TAJO: SUBROUITINE AND FUNCTION ADDED FOR MAKIMA                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine diff(x, n, res)
!!
!! Calculates the difference between adjascent element 
!!

  implicit none 


  integer, intent(in)           :: n               !! number of elements in x
  real(kind=8), intent(in)      :: x(n)            !! array of values
  real(kind=8), intent(out)     :: res(n-1)        !! output array  

  integer                       :: i
  
  do i=1, n-1
    res(i) = x(i+1)-x(i)
  enddo

end subroutine


subroutine makimaSlopes(delta_in , m, s)
!!
!! Derivative values for modified Akima cubic Hermite interpolation
!!
!! Akima's derivative estimate at grid node x(i) requires the four finite
!! differences corresponding to the five grid nodes x(i-2:i+2).
!!
!! For boundary grid nodes x(1:2) and x(n-1:n), append finite differences
!! which would correspond to x(-1:0) and x(n+1:n+2) by using the following
!! uncentered difference formula correspondin to quadratic extrapolation
!! using the quadratic polynomial defined by data at x(1:3)
!! (section 2.3 in Akima's paper):
!!

  implicit none
  
  integer, intent(in)           :: m
  real(kind=8), intent(out)      :: s(m+1)
  real(kind=8), intent(in)      :: delta_in(m)

  integer                       :: n, i
  real(kind=8)                  :: delta_0
  real(kind=8)                  :: delta_m1
  real(kind=8)                  :: delta_n
  real(kind=8)                  :: delta_n1
  real(kind=8)                  :: delta(m+4)
  real(kind=8)                  :: weights(m+3)
  real(kind=8)                  :: weights1(m+1)
  real(kind=8)                  :: weights2(m+1)
  real(kind=8)                  :: weights12(m+1)
  real(kind=8)                  :: delta1(m+1)
  real(kind=8)                  :: delta2(m+1)
  real(kind=8)                  :: eps

  eps = 1.0e-20
  n = m+1
  delta_0  = 2.0*delta_in(1)   - delta_in(2)
  delta_m1 = 2.0*delta_0       - delta_in(1)
  delta_n  = 2.0*delta_in(n-1) - delta_in(n-2)
  delta_n1 = 2.0*delta_n      - delta_in(n-1)

  delta(1) = delta_m1
  delta(2) = delta_0
  delta(3:n+1) = delta_in
  delta(n+2) = delta_n
  delta(n+3) = delta_n1

  !! Akima's derivative estimate formula (equation (1) in the paper):
  !!
  !!       H. Akima, "A New Method of Interpolation and Smooth Curve Fitting
  !!       Based on Local Procedures", JACM, v. 17-4, p.589-602, 1970.
  !!
  !! s(i) = (|d(i+1)-d(i)| * d(i-1) + |d(i-1)-d(i-2)| * d(i))
  !!      / (|d(i+1)-d(i)|          + |d(i-1)-d(i-2)|)
  !!
  !! To eliminate overshoot and undershoot when the data is constant for more
  !! than two consecutive nodes, in MATLAB's 'makima' we modify Akima's
  !! formula by adding an additional averaging term in the weights.
  !! s(i) = ( (|d(i+1)-d(i)|   + |d(i+1)+d(i)|/2  ) * d(i-1) +
  !!          (|d(i-1)-d(i-2)| + |d(i-1)+d(i-2)|/2) * d(i)  )
  !!      / ( (|d(i+1)-d(i)|   + |d(i+1)+d(i)|/2  ) +
  !!          (|d(i-1)-d(i-2)| + |d(i-1)+d(i-2)|/2)
  call diff( delta, m+4, weights)
  weights = abs(weights) + abs((delta(1:m+3)+delta(2:m+4))*0.5) 

  
  weights1 = weights(1:n);   !! |d(i-1)-d(i-2)|
  weights2 = weights(3:n+2); !! |d(i+1)-d(i)|
  delta1 = delta(2:n+1);     !! d(i-1)
  delta2 = delta(3:n+2);     !! d(i)

  weights12 = weights1 + weights2;

  do i=1, n
    !!If the data is constant for more than four consecutive nodes, then the
    !!denominator is zero and the formula produces an unwanted NaN result.
    !!Replace this NaN with 0
    if(abs(weights12(i)) < eps) then
       s(i) = 0.00
    else
       s(i) = (weights2(i)/weights12(i)) * delta1(i) &
            + (weights1(i)/weights12(i)) * delta2(i);
    endif
  enddo

end subroutine

subroutine makima(x, v, n, xq, vq, m)
!!
!!
!!
  implicit none

  integer, intent(in)           :: n
  integer, intent(in)           :: m
  real(kind=8), intent(in)      :: x(n)
  real(kind=8), intent(in)      :: v(n)
  real(kind=8), intent(in)      :: xq(m)
  real(kind=8), intent(out)     :: vq(m)

  integer                       :: i,j
  real(kind=8)                  :: h(n-1)
  real(kind=8)                  :: delta(n-1)
  real(kind=8)                  :: slopes(n)
  real(kind=8)                  :: d(m)
  real(kind=8)                  :: s(m)
  real(kind=8)                  :: t(m)
  

  !! Calculate modified Akima slopes
  call diff(x, n, h)
  call diff(v, n, delta)

  do i=1, n-1
    delta(i) = delta(i) / h(i)
  enddo 

  call makimaSlopes(delta, n-1, slopes)

  !! Evaluate piece wise cubic Hermite interpolants

  do i=1, n-1
    do j=1, m
      if( xq(j) < x(1)) then
        call hermite_cubic_value( x(1), v(1), slopes(1), x(2), v(2), slopes(2), &
                                 1, xq(j), vq(j), d(j), s(j), t(j) )
      else if( xq(j) .eq. x(i) ) then
        vq(j) = v(i)
      else if( xq(j) .eq. x(i+1) ) then
        vq(j) = v(i+1)
      else if( x(i) < xq(j) .and. xq(j) .le. x(i+1) ) then
        call hermite_cubic_value( x(i), v(i), slopes(i), x(i+1), v(i+1), slopes(i+1), &
                                 1, xq(j), vq(j), d(j), s(j), t(j) )
      else if ( xq(j) > x(n) ) then
        call hermite_cubic_value( x(n-1), v(n-1), slopes(n-1), x(n), v(n), slopes(n), &
                                 1, xq(j), vq(j), d(j), s(j), t(j) )
      endif
   enddo
  enddo

end subroutine

subroutine hermite_cubic_value ( x1, f1, d1, x2, f2, d2, n, x, f, d, s, t )

!*****************************************************************************80
!
!! HERMITE_CUBIC_VALUE evaluates a Hermite cubic polynomial.
!
!  Discussion:
!
!    The input arguments can be vectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 2011
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, F1, D1, the left endpoint, function value
!    and derivative.
!
!    Input, real ( kind = 8 ) X2, F2, D2, the right endpoint, function value
!    and derivative.
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the points at which the Hermite cubic
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) F(N), D(N), S(N), T(N), the value and first
!    three derivatives of the Hermite cubic at X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c2
  real ( kind = 8 ) c3
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) df
  real ( kind = 8 ) dx(n)
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) h
  real ( kind = 8 ) s(n)
  real ( kind = 8 ) t(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2

  h =    x2 - x1
  df = ( f2 - f1 ) / h

  c2 = - ( 2.0 * d1 - 3.0 * df + d2 ) / h
  c3 =   (       d1 - 2.0 * df + d2 ) / h / h

  dx(1:n) = x(1:n) - x1
  f(1:n) = f1 + dx(1:n) &
             * ( d1 + dx(1:n) * (           c2 + dx(1:n) *           c3 ) )
  d(1:n) =       d1 + dx(1:n) * ( 2.0D+00 * c2 + dx(1:n) * 3.0D+00 * c3 )
  s(1:n) =                        2.0D+00 * c2 + dx(1:n) * 6.0D+00 * c3
  t(1:n) =                                                 6.0D+00 * c3

  return
end


end  !! end of module 


