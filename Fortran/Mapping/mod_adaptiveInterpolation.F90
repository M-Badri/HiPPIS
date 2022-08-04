
module mod_adaptiveInterpolation
!!
!!  Module for adaptive polynomial data-bounded and positivity-preserving interpolation.  
!!
!!  Below are are the subroutine in this module 
!!
!!    divdiff(...)
!!
!!    divdiff_vec(...)
!!
!!    newtonPolyVal(...)
!!
!!    adaptiveInterpolation1D(...)
!!
!!    adaptiveInterpolation1D_vec(...)
!!
!!    adaptiveInterpolation2D(...)
!!
!!    adaptiveInterpolation2D_vec(...)
!!
!!    adaptiveInterpolation3D(...)
!!
!!    adaptiveInterpolation3D_vec(...)
!!

  implicit none

  contains

subroutine divdiff_vec(x, y, n, d, table)
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


  !$OMP SIMD
  do i=1,n
    table(i,1) = y(i);
  enddo
  
  do j=2,d+1
    !$OMP SIMD
    do i=1,n-(j-1)
      table(i,j) = (table(i+1, j-1)-table(i, j-1)) / (x(i+j-1)-x(i))
    enddo
    
  enddo

endsubroutine

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
    do i=1,n
      if(i+j-1 <= n) then
        table(i,j) = (table(i+1, j-1)-table(i, j-1)) / (x(i+j-1)-x(i))
      endif
    enddo
  enddo


end subroutine 
 


subroutine newtonPolyVal(x, u, d, xout, yout)
!!! This function builds up the newton interpolant and evaluates it at xout
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
  !$OMP SIMD 
  do i=d,1,-1
    yout = yout * (xout -x(i)) + u(i)
  enddo

end subroutine

subroutine adaptiveinterpolation1D(x, y, n, xout, yout, m, degree, interpolation_type, st, eps0, eps1, deg ) 
!!
!! This function is a polynomial interpoaltion method that builds a piece-wise function based on the input (x,y).
!! The piece-wise function is then evaluate at the output points xout to give (xout, yout).
!! The interpolation method preserves positivity or data-boundedness. The data-bounded or
!! positivity-preserving interpolant is constructed for each interval
!! based on the theory in https://arxiv.org/abs/2204.06168
!! and the algorithm in the manuscript [REF]. 
!! 
!!
!! INPUT: 
!! n: the number points in the 1D vector x.
!! m: the number of points in the 1D vector xout.
!! x: 1D mesh points of length n. For i=1, ..., n-1 x_{i} >  x_{i+1}
!! y: 1D vector that have the data values associated with the points x_{i} for i=1, ..., n
!! xout: 1D vector of length m that represent the locations where we which the interpolate to.
!! Interpolation_type: used to determine the type of interpolation to be used to build interpolant.
!!   - interpolation_type=1: a data-bounded interpolant is built for each interpolant.
!!   - interpolation_type=2: a positivity-preserving interpolant is built for each interpolant.
!! degree: target polynomial degree and maximum polynomial degree used for each interval.
!! st (optional): used guide point selection process in cases when adding the next point to the 
!!   right or left both meet the requirements for positivity or datat-boundedness.
!!   - st=1 (default): the point with the smallest divided difference is added (ENO stencil).
!!   - st=2 the point to the left of current stencil is selected if the number of point to left
!!     of x_{i} is smaller than the number of points to right of x_{i} (i-si < ei-i). Similarly, 
!!     the point to the right is selected if the number of points to the right of x_{i} is smaller
!!     than the number of points to the left (i-si > ei-i). When both the number of points to right 
!!     and left are the same, the algorithm chooses the point with the smallest lambda.  
!!   - st=3 the point that is closest to the starting interval is chosen.
!! eps0 (optional): positive parameter use constrain the bound of the positive interpolant in intervals with no
!!   extremum detected.
!! eps1 (optional): positive parameter use constrain the bound of the positive interpolant in intervals with 
!!   extremum detected.
!!
!! OUTPUT:
!! yout: results of evaluating interpolants at the locations xout.
!! deg (optional): 1D vector that holds the degree of the interpolant used for each interval
!!

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


  !!** Local variables
  real(kind=8)                          :: u(degree+1)                  !! to store divided differences associated with selected points in xval
  real(kind=8)                          :: e
  real(kind=8)                          :: up_b
  real(kind=8)                          :: low_b
  real(kind=8)                          :: lambda
  real(kind=8)                          :: prod_deltax
  real(kind=8)                          :: xval(degree+1)               !! to store slected points in order
  real(kind=8)                          :: table(n, degree+4)           !! table of devided diferences
  real(kind=8)                          :: ur, ul, ww, lambda_l, lambda_r, m_lambda, m_sigma
  real(kind=8)                          :: sigma_l, sigma_r
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
  integer                               :: stencil_type
  integer, parameter                    :: Debug = 0
  logical                               :: bool_left, bool_right
  


  !!!** Check input to make sure x_{i} < x_{i+1} **!!
  !do i=1, n-1
  !  if( x(i) .ge. x(i+1) .or. abs(x(i+1)-x(i)) .le. eps ) then
  !   write(*,*)'ERROR: Incorrect input at i=', i, 'x(i)=', x(i),&
  !               'x(i+1)=',x(i+1), 'x(i) must be less that x(i+1) and &
  !              |x(i+1)-x(i)| mus be greater than machine precision eps.'
  !   call exit(0)
  !  endif
  !enddo

  !!!** Check output to make sure x_{i} < x_{i+1} **!!
  !do i=1, m-1
  !  if( xout(i) .ge. xout(i+1) ) then
  !    write(*,*)'ERROR: Incorrect output at k=', i, 'xout(k)=', xout(i),&
  !               'xout(k+1)=',xout(i+1), 'xout(k) must be less that xout(k+1)'
  !    call exit(0)
  !  endif
  !enddo


  !!** Initialize variables **!!
  eps = 1e-30                   !! defined as epsilon 
  inv_eps = 1e+30               !! defined to be + infinity
  k = 1                         !! iteration idex used for output points

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
 

  !call divdiff_vec(x, y, n, degree+3, table)    !!compute the table of divided differences 
  call divdiff(x, y, n, degree+1, table)    !!compute the table of divided differences 

  !!** Calculate slopes for each interval **!!
  slope(1) = (y(3)-y(2))/(x(3)-x(2))  !! left boundary
  do i=1, n-1
    slope(i+1) = (y(i+1)-y(i))/(x(i+1)-x(i))  !! right boundary
  enddo
  slope(n+1) = slope(n-1)
  

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

  !!!write(*,*) ' TAJO ml=', mm_l(n-1), 'mr=', mm_r(n-1)
  !!!** loop over each input intervals. For each  interval build an interpolant and 
  !!!   evaluate the interpolant at the desired output points **!!
  do i=1,n-1
 
    !!** Initialize varibles for each interval**!!
    u = 0.0
    xval = 0.0
    prod_deltax = 1.0
    !!lambda = 0.0
    !!sigma = 1.0

    xval(1) = x(i)                       !! first point in stencil
    xval(2) = x(i+1)                     !! second point in stencil
    u(1) = y(i)                          !! set first selected divided difference
    u(2)= (y(i+1)-y(i))/(x(i+1)-x(i))    !! set second slected divided difference
    lambda = 1.0                      !! set first ratio of divided difference  
    up_b = 1.0
    low_b = -1.0
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
         prod_deltax_l = prod_deltax * (x(ei)-x(tmp_si)) !! product of interval containing stencil
         if(si-1 > 0)then
           ul = table(tmp_si, ei-tmp_si+1)   !! get left divided difference   
           lambda_l = ul/ww * prod_deltax_l  !! calculate left lambda         
           xl = x(si-1)
         else
           ul = inv_eps  !! set left dividided difference to + infinity  
           lambda_l = inv_eps    !! set left lambda to infinity
           xl = -inv_eps
         endif

         prod_deltax_r = prod_deltax * (x(tmp_ei)-x(si))    !! product of inteval containing stencil
         if(ei+1 .le. n)then
           ur = table(si, tmp_ei-si+1)    !! get right divided difference 
           lambda_r =  ur/ ww * prod_deltax_r  !! calculate righ lambda
           xr = x(ei+1)
         else
           ur = inv_eps    !! set righl divided difference to infinity
           lambda_r =  inv_eps   !! set righ lambda to infinity
           xr = inv_eps
         endif

         e = -(x(i)-xval(j))/(x(i+1)-x(i))  !! calculate e_j !!
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
           if(e <= 0.0) then
             up_b_l = (up_b - lambda)* d_l / (1.0-e)
             low_b_l = (low_b - lambda)*d_l / (1.0-e)
           elseif(e > 0.0) then
             up_b_l = (low_b - lambda)*d_l / (0.0-e)
             low_b_l = (up_b - lambda)* d_l / (0.0-e)
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
           if(e <= 0.0) then
             up_b_r = (up_b - lambda)* d_r / (1.0-e)
             low_b_r = (low_b - lambda)*d_r / (1.0-e)
           elseif(e > 0.0) then
             up_b_r = (low_b - lambda)*d_r / (0.0-e)
             low_b_r = (up_b - lambda)* d_r / (0.0-e)
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

         !!** Add point to the left of current stencil and corresponding 
         !!   variables **!!
         if( (bool_left .eqv. .true.) .and. (bool_right .eqv. .false.)) then
           si = max(1, si-1)
           !!ei = ei
           lambda = lambda_l
           u(j+1) = ul
           xval(j+1) = x(si)
           up_b = up_b_l
           low_b = low_b_l
           prod_deltax = prod_deltax_l
         !!** Add point to the right of current stencil and corresponding 
         !!   variables **!!
         elseif( (bool_left .eqv. .false.) .and. (bool_right .eqv. .true.) ) then
           !!si = si
           ei = min(ei+1, n)
           lambda = lambda_r
           u(j+1) = ur
           xval(j+1) = x(ei)
           up_b = up_b_r
           low_b = low_b_r
           prod_deltax = prod_deltax_r
         endif

      enddo !! j loop
    endif !! end of if j > 1


    !!** save the interpolant degree used for the interval [x_{i}, x_{i+1}] **!!
    if(present(deg)) then
      deg(i) = ei-si
    endif

    !!!** Extrapolate to points that are to the left of the defined interval **!! 
    !if(k <=m)then
    !  do while( xout(k) < x(1) )
    !    write(*,*) 'WARNING: Some of the output are obtained via &
    !                extrapolation instead of interpolation. The desired &
    !                proprety such as data-boundedness or positvity is not &
    !                preserved in such case'
    !      write(*,*)  k, 1 
    !      write(*,*)  xout(k), x(1) 
    !    call newtonPolyVal(xval, u, degree, xout(k), yout(k))
    !    k = k+1
    !    if(k > m) exit
    !  enddo
    !endif
 
    !!** Building and evaluating Interpolant at xout points **!!
    !! - do while( x(i) <= xout(k) .and. xout(k) <= x(i+1) .and. k <= m )
    if( k <=m)then
      do while( x(i) <= xout(k) .and. xout(k) <= x(i+1) )
        call newtonPolyVal(xval, u, degree, xout(k), yout(k))
        k = k+1
        if(k > m) exit
      enddo
    endif

    !!!** Extrapolate to points that are to the right of the defined interval **!! 
    !if(k <= m)then
    !  do while( xout(k) > x(n) )
    !      write(*,*) 'WARNING: Some of the output are obtained via &
    !                  extrapolation instead of interpolation. The desired &
    !                  proprety such as data-boundedness or positvity is not &
    !                  preserved in such case'
    !      write(*,*)  k, n 
    !      write(*,*)  xout(k), x(n) 
    !      call newtonPolyVal(xval, u, degree, xout(k), yout(k))
    !      k = k+1
    !      if(k > m) exit
    !  enddo
    !endif

  end do !! end of i loop

end subroutine !!adaptiveInterpolation1D

subroutine adaptiveinterpolation1D_vec(x, y, n, xout, yout, m, degree, interpolation_type, st, eps0, eps1, deg) 
!!
!! This function is a polynomial interpoaltion method that builds a piece-wise function based on the input (x,y).
!! The piece-wise function is then evaluate at the output points xout to give (xout, yout).
!! The interpolation method preserves positivity or data-boundedness. The data-bounded or
!! positivity-preserving interpolant is constructed for each interval
!! based on the theory in https://arxiv.org/abs/2204.06168
!! and the algorithm in the manuscript [REF]. 
!! 
!!
!! INPUT: 
!! n: the number points in the 1D vector x.
!! m: the number of points in the 1D vector xout.
!! x: 1D mesh points of length n. For i=1, ..., n-1 x_{i} >  x_{i+1}
!! y: 1D vector that have the data values associated with the points x_{i} for i=1, ..., n
!! xout: 1D vector of length m that represent the locations where we which the interpolate to.
!! Interpolation_type: used to determine the type of interpolation to be used to build interpolant.
!!   - interpolation_type=1: a data-bounded interpolant is built for each interpolant.
!!   - interpolation_type=2: a positivity-preserving interpolant is built for each interpolant.
!! degree: target polynomial degree and maximum polynomial degree used for each interval.
!! st (optional): used guide point selection process in cases when adding the next point to the 
!!   right or left both meet the requirements for positivity or datat-boundedness.
!!   - st=1 (default): the point with the smallest divided difference is added (ENO stencil).
!!   - st=2 the point to the left of current stencil is selected if the number of point to left
!!     of x_{i} is smaller than the number of points to right of x_{i} (i-si < ei-i). Similarly, 
!!     the point to the right is selected if the number of points to the right of x_{i} is smaller
!!     than the number of points to the left (i-si > ei-i). When both the number of points to right 
!!     and left are the same, the algorithm chooses the point with the smallest lambda.  
!!   - st=3 the point that is closest to the starting interval is chosen.
!! eps0 (optional): positive parameter use constrain the bound of the positive interpolant in intervals with no
!!   extremum detected.
!! eps1 (optional): positive parameter use constrain the bound of the positive interpolant in intervals with 
!!   extremum detected.
!!
!! OUTPUT:
!! yout: results of evaluating interpolants at the locations xout.
!! deg (optional): 1D vector that holds the degree of the interpolant used for each interval
!!
!!
 
  use omp_lib

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

  !!** Local variables
  real(kind=8)                          :: u(degree+1)                  !! to store divided differences associated with selected points in xval
  real(kind=8)                          :: e
  real(kind=8)                          :: up_b
  real(kind=8)                          :: low_b
  real(kind=8)                          :: lambda
  real(kind=8)                          :: prod_deltax
  real(kind=8)                          :: xx(n-1)
  real(kind=8)                          :: xval(degree+1)               !! to store slected points in order
  real(kind=8)                          :: x_left(n-1), x_right(n-1)
  real(kind=8)                          :: table(n, degree+4)           !! table of devided diferences
  real(kind=8)                          :: ur, ul, ww, lambda_l, lambda_r, m_lambda, m_sigma
  real(kind=8)                          :: m_l, m_r
  real(kind=8)                          :: mm_l(n-1), mm_r(n-1), www(n-1)
  real(kind=8)                          :: wr1(n-1), wr2(n-1), wr3(n-1), wr4(n-1)
  real(kind=8)                          :: u_left(n-1), u_right(n-1)
  real(kind=8)                          :: lambda_new(n-1), u_new(n-1)
  real(kind=8)                          :: lambda_left(n-1), lambda_right(n-1)
  real(kind=8)                          :: prod_deltax_left(n-1)
  real(kind=8)                          :: prod_deltax_right(n-1)
  real(kind=8)                          :: prod_deltaxx(n-1)
  real(kind=8)                          :: B_minus_l(n-1)
  real(kind=8)                          :: B_minus_r(n-1)
  real(kind=8)                          :: B_plus_l(n-1)
  real(kind=8)                          :: B_plus_r(n-1)
  real(kind=8)                          :: B_plus(n-1)
  real(kind=8)                          :: B_minus(n-1)
  real(kind=8)                          :: prod_deltax_l, prod_deltax_r
  real(kind=8)                          :: a, b
  real(kind=8)                          :: slope(n+1), slope_i, slope_im1, slope_ip1 , tol
  real(kind=8)                          :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, delta
  real(kind=8)                          :: eps, inv_eps, eps2, eps3 
  real(kind=8)                          :: umax, umin                        !! parameter  used to for upper bound for each interval
  real(kind=8)                          :: xl , xr 
  real(kind=8)                          :: ml1, ml2, ml3, ml4, ml5, ml6 
  real(kind=8)                          :: mr1, mr2, mr3, mr4, mr5, mr6 
  real(kind=8)                          :: eps_l, eps_r, w1, w2 
  real(kind=8)                          :: left, right, neither, mask 
  integer                               :: bool(n-1)
  integer                               :: bool2(n-1)
  integer                               :: i, j, k, kk 
  integer                               :: si, ei 
  integer                               :: tmp_si, tmp_ei
  integer                               :: tmp_idx, fid
  integer                               :: stencil_type
  integer                               :: f_si(n-1), f_ei(n-1)
  logical                               :: b1(n-1), b2(n-1), b3(n-1), b4(n-1)
  logical                               :: b5(n-1), b6(n-1), b7(n-1)

  



  !!!** Check input to make sure x_{i} < x_{i+1} **!!
  !do i=1, n-1
  !  if( x(i) .ge. x(i+1) .or. abs(x(i+1)-x(i)) .le. eps ) then
  !   write(*,*)'ERROR: Incorrect input at i=', i, 'x(i)=', x(i),&
  !               'x(i+1)=',x(i+1), 'x(i) must be less that x(i+1) and &
  !              |x(i+1)-x(i)| mus be greater than machine precision eps.'
  !   call exit(0)
  !  endif
  !enddo

  !!!** Check output to make sure x_{i} < x_{i+1} **!!
  !do i=1, m-1
  !  if( xout(i) .ge. xout(i+1) ) then
  !    write(*,*)'ERROR: Incorrect output at k=', i, 'xout(k)=', xout(i),&
  !               'xout(k+1)=',xout(i+1), 'xout(k) must be less that xout(k+1)'
  !    call exit(0)
  !  endif
  !enddo


  !!** Initialize variables **!!
  eps = 1e-30                   !! defined as epsilon 
  inv_eps = 1e+30               !! defined to be + infinity
  k = 1                         !! iteration idex used for output points

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
 
  if(present(st) )then
    stencil_type  = st
  else
    stencil_type = 1
  endif 
  

  call divdiff_vec(x, y, n, degree+1, table)    !!compute the table of divided differences 

  !!** Calculate slopes for each interval **!!
  slope(1) = (y(3)-y(2))/(x(3)-x(2))  !! left boundary
  !$OMP SIMD PRIVATE(i)
  do i=1, n-1
    slope(i+1) = (y(i+1)-y(i))/(x(i+1)-x(i))  !! right boundary
  enddo
  
  slope(n+1) = slope(n-1)

  !!** Calculate the polynomial bounds for each interval **!!
  if(degree > 1 .and. interpolation_type .eq. 2) then
    !$OMP SIMD PRIVATE(i)
    do i=1,n-1
      bool(i) = ( (slope(i)*slope(i+2) < 0.0 .and. slope(i) < 0.0) .or. &          !! Detects a minimum
          (slope(i)*slope(i+2) > 0.0 .and. slope(i)*slope(i+1) < 0.0) )  !! Detects a maximum and/or minimum (ambiguous).
    enddo
    
    !$OMP SIMD PRIVATE(i)
    do i=1, n-1
      bool(i) = abs(bool(i))
    enddo
    
    !$OMP SIMD PRIVATE(i, tmp1, eps_l)
    do i=1,n-1
      eps_l = bool(i)*eps3 + (1-bool(i))*eps2
      tmp1 = min(y(i), y(i+1))
      mm_l(i) = tmp1 - eps_l*abs(tmp1) 
    enddo
    
    !$OMP SIMD PRIVATE(i)
    do i=1,n-1
      bool(i) = ( (slope(i)*slope(i+2) < 0.0 .and. slope(i) > 0.0) .or. &          !! Detects a minimum
          (slope(i)*slope(i+2) > 0.0 .and. slope(i)*slope(i+1) < 0.0) )  !! Detects a maximum and/or minimum (ambiguous).
    enddo
    
    !$OMP SIMD PRIVATE(i)
    do i=1,n-1
      bool(i) = abs(bool(i))
    enddo
    
    !$OMP SIMD PRIVATE(i, tmp2, eps_r)
    do i=1,n-1
      eps_r = bool(i)*eps3 + (1-bool(i))*eps2
      tmp2 = max(y(i), y(i+1))
      mm_r(i) = tmp2 + eps_r*abs(tmp2) 
    enddo
    


    !$OMP SIMD PRIVATE(i, umin, umax)
    do i=1, n-1
       umin = mm_l(i)
       umax = mm_r(i)
       wr1(i) = min( (umin-y(i)) / (y(i+1)-y(i)+eps),  (umax-y(i)) /  (y(i+1)-y(i)+eps) )
       wr2(i) = max( (umin-y(i)) / (y(i+1)-y(i)+eps),  (umax-y(i)) /  (y(i+1)-y(i)+eps) )
    enddo

    !$OMP SIMD PRIVATE(i)
    do i=1, n-1
      si = max(i-1,1)
      ei = min(i+2,n)
      bool(i) = (y(i) == y(i+1) .and. y(si) .ne. y(i) .and. y(i+1) .ne. y(ei) )
    enddo
    

    !$OMP SIMD PRIVATE(i, tmp_si)
    do i=1, n-1
      tmp_si = max(i-1, 1)
      u_new(i) = table(i,3)*(x(i+1)-x(i)) * (x(i+1)-x(tmp_si)) 
    enddo
    
    !$OMP SIMD PRIVATE(i, ww, umin, umax,tmp3)
    do i=1, n-2
      ww = u_new(i) ! ul*(x(i+1)-x(i)) * (x(i+1)-x(tmp_si)) 
      umin = mm_l(i)
      umax = mm_r(i)
      tmp3 = (y(i+1)-y(i)) / (x(i+1)-x(i))
      wr3(i) = (1-bool(i))*wr1(i) + bool(i)* min( (umin-y(i)) / (ww+eps),  (umax-y(i)) /  (ww+eps) )
      wr4(i) = (1-bool(i))*wr2(i) + bool(i)* max( (umin-y(i)) / (ww+eps),  (umax-y(i)) /  (ww+eps) )
      www(i) = (1-bool(i))*tmp3 + bool(i)*ww
    enddo
    
    ww = u_new(n-2) ! ul*(x(i+1)-x(i)) * (x(i+1)-x(tmp_si)) 
    umin = mm_l(n-1)
    umax = mm_r(n-1)
    tmp3 = (y(n)-y(n-1)) / (x(n)-x(n-1))
    wr3(n-1) = (1-bool(n-1))*wr1(n-1) + bool(n-1)* min( (umin-y(n-1)) / (ww+eps),  (umax-y(n-1)) /  (ww+eps) )
    wr4(n-1) = (1-bool(n-1))*wr2(n-1) + bool(n-1)* max( (umin-y(n-1)) / (ww+eps),  (umax-y(n-1)) /  (ww+eps) )
    www(n-1) = (1-bool(n-1))*tmp3 + bool(n-1)*ww
 

    !$OMP SIMD PRIVATE(i)
    do i=1, n-1
      mm_l(i) = min(0.0, wr3(i))
      mm_r(i) = max(1.0, wr4(i))
    enddo
    
  !!** Default case: DBI **!!
  else  
    !$OMP SIMD PRIVATE(i)
    do i=1, n-1
      !!** Compute the values of m_{\ell} and m_r for the data-bounded method 
      !!   if the limiter variable is set to 1 **!!
      mm_l(i) = 0.0
      mm_r(i) = 1.0
      www(i) = (y(i+1)-y(i)) / (x(i+1)-x(i)) 
    enddo
    
  endif

  !$OMP SIMD PRIVATE(i)
  do i=1, n-1
    f_si(i) = i
    f_ei(i) = i+1
    prod_deltaxx(i) = 1.0
    lambda_new(i) = 1.0
    B_plus(i) = 1.0
    B_minus(i) = -1.0
    xx(i) = x(i+1)
  enddo
  

  !!
  do j=2, degree

    !! Compute left
    !$OMP SIMD PRIVATE(i)
    do i=1, n-1
      bool(i) = (f_si(i)-1> 0)
      bool(i) = abs(bool(i))
    enddo
    
    !!$OMP SIMD 
    !do i=1, n-1
    !  bool(i) = abs(bool(i))
    !enddo
    
    !!
    !$OMP SIMD PRIVATE(i, tmp_si, ei)
    do i=1, n-1
      tmp_si = max(f_si(i)-1,1)
      ei = f_ei(i)
      u_left(i) = bool(i)*table(tmp_si, ei-tmp_si+1) + (1-bool(i))*inv_eps
      x_left(i) = bool(i)*x(tmp_si) + (1-bool(i))*inv_eps
      prod_deltax_left(i) = prod_deltaxx(i) * (x(ei)-x(tmp_si)) !! product of interval containing stencil
    enddo
    
    !!
    !$OMP SIMD PRIVATE(i)
    do i=1, n-1
      lambda_left(i) = bool(i)*u_left(i)/(www(i)+eps) * prod_deltax_left(i)+ &  !! calculate left lambda         
                       (1-bool(i))*inv_eps
    enddo 
    

    !! Compute Right
    !$OMP SIMD PRIVATE(i)
    do i=1, n-1
      bool(i) = (f_ei(i)+1<= n)
      bool(i) = abs(bool(i))
    enddo
    !!$OMP SIMD 
    !do i=1, n-1
    !  bool(i) = abs(bool(i))
    !enddo
    
    !!
    !$OMP SIMD PRIVATE(i, tmp_ei, si)
    do i=1, n-1
      si = f_si(i)
      tmp_ei = min(f_ei(i)+1,n)
      u_right(i) = bool(i)*table(si, tmp_ei-si+1) + (1-bool(i))*inv_eps
      prod_deltax_right(i) = prod_deltaxx(i) * (x(tmp_ei)-x(si)) !! product of interval containing stencil
      x_right(i) = bool(i)*x(tmp_ei) + (1-bool(i))*inv_eps
    enddo
    
    !!
    !$OMP SIMD PRIVATE(i)
    do i=1, n-1
      lambda_right(i) = bool(i)*u_right(i)/(www(i)+eps) * prod_deltax_right(i) + &
                        (1-bool(i))*inv_eps !! calculate left lambda         
    enddo
    

    if(j ==2) then
      !$OMP SIMD PRIVATE(i, tmp_si, tmp_ei, si, ei)
      do i=1, n-1
        tmp_si = max(f_si(i)-1,1)
        tmp_ei = min(f_ei(i)+1,n)
        si = f_si(i)
        ei = f_ei(i)
        wr1(i) = 1.0  !! calculate e_j !!
        wr2(i) = (x(ei)-x(tmp_si))/(x(i+1)-x(i)) !! calculate d_l 
        wr3(i) = (x(tmp_ei)-x(si))/(x(i+1)-x(i)) !! calculate d_r

      enddo
      

      !$OMP SIMD PRIVATE(i)
      do i=1,n-1
         B_plus_r(i)  = wr3(i)*( -mm_l(i)*4.0 + 1.0 )
         B_minus_r(i) = wr3(i)*( -(mm_r(i)-1.0)*4.0 - 1.0 )  
         B_plus_l(i)  = wr2(i)*( -mm_l(i)*4.0 + 1.0 )
         B_minus_l(i) = wr2(i)*( -(mm_r(i)-1.0)*4.0 - 1.0 )  
      enddo
      
    else ! j> 2
      !$OMP SIMD PRIVATE(i, tmp_si, tmp_ei, si, ei)
      do i=1, n-1
        tmp_si = max(f_si(i)-1,1)
        tmp_ei = min(f_ei(i)+1,n)
        si = f_si(i)
        ei = f_ei(i)
        wr1(i) = -(x(i)-xx(i))/(x(i+1)-x(i))  !! calculate e_j !!
        wr2(i) = (x(ei)-x(tmp_si))/(x(i+1)-x(i)) !! calculate d_l 
        wr3(i) = (x(tmp_ei)-x(si))/(x(i+1)-x(i)) !! calculate d_r
      enddo
      
      !$OMP SIMD PRIVATE(i)
      do i=1, n-1
        bool(i) = (wr1(i) <= 0.0)
        bool(i) = abs(bool(i))
      enddo
      
      !!$OMP SIMD PRIVATE(i)
      !do i=1, n-1
      !  bool(i) = abs(bool(i))
      !enddo
      

      !$OMP SIMD PRIVATE(i)
      do i=1, n-1
        B_minus_l(i) = bool(i)*(B_minus(i)-lambda_new(i))*wr2(i)/(1.0-wr1(i)+eps) + &
                     (1-bool(i))*(B_plus(i)-lambda_new(i))*wr2(i)/(0.0-wr1(i)+eps)
        !!
        B_plus_l(i) = bool(i)*(B_plus(i)-lambda_new(i))*wr2(i)/(1.0-wr1(i)+eps) + &
                     (1-bool(i))*(B_minus(i)-lambda_new(i))*wr2(i)/(0.0-wr1(i)+eps)
        !!
        B_minus_r(i) = bool(i)*(B_minus(i)-lambda_new(i))*wr3(i)/(1.0-wr1(i)+eps) + &
                       (1-bool(i))*(B_plus(i)-lambda_new(i))*wr3(i)/(0.0-wr1(i)+eps)
        !!
        B_plus_r(i)= bool(i)*(B_plus(i)-lambda_new(i))*wr3(i)/(1.0-wr1(i)+eps) + &
                     (1-bool(i))*(B_minus(i) - lambda_new(i))*wr3(i)/(0.0-wr1(i)+eps)
      enddo
      

    endif

    !!** Option 1: stencil_type = 1. In addition to positivity or 
    !!   data boundedness, the stencil selection is based on the ENO approach **!!
    if(stencil_type .eq. 1) then
      !$OMP SIMD PRIVATE(i)
      do i=1, n-1
        b1(i) = (B_minus_l(i) .le. lambda_left(i)) .and. (lambda_left(i) .le. B_plus_l(i)) .and. & !! Adding a point to left meets the requiremenst for DBI or PPI
               (B_minus_r(i) .le. lambda_right(i)) .and. (lambda_right(i) .le. B_plus_r(i)) .and. &   !! Adding a point to right meets the requiremenst for DBI or PPI
                (abs(u_left(i)) < abs(u_right(i)) )
        !!
        b2(i) =(B_minus_l(i) .le. lambda_left(i)) .and. (lambda_left(i) .le. B_plus_l(i)) .and. & !! Adding a point to left meets the requiremenst for DBI or PPI
               (B_minus_r(i) .le. lambda_right(i)) .and. (lambda_right(i) .le. B_plus_r(i)) .and. &   !! Adding a point to right meets the requiremenst for DBI or PPI
                (abs(u_left(i)) >= abs(u_right(i)))
      enddo
      

      !$OMP SIMD PRIVATE(i)
      do i=1, n-1
          b3(i) = (B_minus_r(i) .le. lambda_right(i)) .and. (lambda_right(i) .le. B_plus_r(i)) .and. &   !! Adding a point to right meets the requiremenst for DBI or PPI
                      ( b1(i) .eqv. .false.)  .and. (b2(i) .eqv. .false.)  
      enddo
      
      !$OMP SIMD PRIVATE(i)
      do i=1, n-1
        b4(i) = (B_minus_l(i) .le. lambda_left(i)) .and. (lambda_left(i) .le. B_plus_l(i)) .and. & !! Adding a point to left meets the requiremenst for DBI or PPI
                 (b1(i) .eqv. .false.) .and. (b2(i) .eqv. .false.) .and. (b3(i) .eqv. .false.)  
      enddo
      

      !$OMP SIMD PRIVATE(i)
      do i=1, n-1
        bool(i) = (b1(i) .or. b4(i))
        bool2(i) = (b2(i) .or. b3(i))
        bool(i) = abs(bool(i))
        bool2(i) = abs(bool2(i))
      enddo
      
      !!$OMP SIMD 
      !do i=1, n-1
      !  bool(i) = abs(bool(i))
      !  bool2(i) = abs(bool2(i))
      !enddo
      

    !! Option 2: stencil_type = 2. In addition to DBI or PPI the 
    !! stencil selection prioritize a symetric stencil other others **!!
    elseif(stencil_type == 2) then
      !$OMP SIMD PRIVATE(si, ei, i)
      do i=1, n-1
        si = f_si(i)
        ei = f_ei(i)
        b1(i) =( (B_minus_l(i) .le. lambda_left(i) .and. lambda_left(i) .le. B_plus_l(i)) .and. & !! Adding a point to left meets the requiremenst for DBI or PPI
               (B_minus_r(i) .le. lambda_right(i) .and. lambda_right(i) .le. B_plus_r(i)) .and. &   !! Adding a point to right meets the requiremenst for DBI or PPI
                i-si < ei-i )
        !!
        b2(i) =( (B_minus_l(i) .le. lambda_left(i) .and. lambda_left(i) .le. B_plus_l(i)) .and. & !! Adding a point to left meets the requiremenst for DBI or PPI
               (B_minus_r(i) .le. lambda_right(i) .and. lambda_right(i) .le. B_plus_r(i)) .and. &   !! Adding a point to right meets the requiremenst for DBI or PPI
                i-si > ei-i )
        !!
        b3(i) =( (B_minus_l(i) .le. lambda_left(i) .and. lambda_left(i) .le. B_plus_l(i)) .and. & !! Adding a point to left meets the requiremenst for DBI or PPI
               (B_minus_r(i) .le. lambda_right(i) .and. lambda_right(i) .le. B_plus_r(i)) .and. &   !! Adding a point to right meets the requiremenst for DBI or PPI
                i-si == ei-i )
      enddo
      
      !$OMP SIMD PRIVATE(i)
      do i=1, n-1
        b4(i) = (abs(lambda_left(i)) < abs(lambda_right(i))) .and. (b3(i) .eqv. .true.) .and. &
                   (b2(i) .eqv. .false.) .and. (b1(i) .eqv. .false.)
        b5(i) = (abs(lambda_left(i)) >= abs(lambda_right(i))) .and. (b3(i) .eqv. .true.) .and. &
                   (b2(i) .eqv. .false.) .and. (b1(i) .eqv. .false.)
      enddo
      
      !$OMP SIMD PRIVATE(i)
      do i=1, n-1
          b6(i) = (B_minus_r(i) .le. lambda_right(i)) .and. (lambda_right(i) .le. B_plus_r(i)) .and. &   !! Adding a point to right meets the requiremenst for DBI or PPI
                    (b1(i) .eqv. .false.)  .and. (b2(i) .eqv. .false.) .and. (b3(i) .eqv. .false.)  
      enddo
      
      !$OMP SIMD PRIVATE(i) 
      do i=1, n-1
        b7(i) = (B_minus_l(i) .le. lambda_left(i)) .and. (lambda_left(i) .le. B_plus_l(i)) .and. & !! Adding a point to left meets the requiremenst for DBI or PPI
                (b1(i) .eqv. .false.)  .and. (b2(i) .eqv. .false.) .and. (b3(i) .eqv. .false.) .and. (b6(i) .eqv. .false. )
      enddo
      

      !$OMP SIMD PRIVATE(i)
      do i=1, n-1
        bool(i) = (b1(i) .or. b4(i) .or. b7(i))
        bool2(i) = (b2(i) .or. b5(i) .or. b6(i))
        bool(i) = abs(bool(i))
        bool2(i) = abs(bool2(i))
      enddo
      

      !!$OMP SIMD 
      !do i=1, n-1
      !  bool(i) = abs(bool(i))
      !  bool2(i) = abs(bool2(i))
      !enddo
      

    !! Option 3: stencil_type = 3. In addition to DBI or PPI the 
    !! stencil selection prioritize a locality around the starting interval **!!
    elseif(stencil_type == 3) then
      !$OMP SIMD PRIVATE(i, si, ei)
      do i=1, n-1
        si = f_si(i)
        ei = f_ei(i)
        b1(i) =( (B_minus_l(i) .le. lambda_left(i) .and. lambda_left(i) .le. B_plus_l(i)) .and. & !! Adding a point to left meets the requiremenst for DBI or PPI
               (B_minus_r(i) .le. lambda_right(i) .and. lambda_right(i) .le. B_plus_r(i)) .and. &   !! Adding a point to right meets the requiremenst for DBI or PPI
                abs(x(i)-x_left(i)) < abs(x_right(i)-x(i+1)) )
        !!
        b2(i) =( (B_minus_l(i) .le. lambda_left(i) .and. lambda_left(i) .le. B_plus_l(i)) .and. & !! Adding a point to left meets the requiremenst for DBI or PPI
               (B_minus_r(i) .le. lambda_right(i) .and. lambda_right(i) .le. B_plus_r(i)) .and. &   !! Adding a point to right meets the requiremenst for DBI or PPI
                abs(x(i)-x_left(i)) > abs(x_right(i)-x(i+1)) )
        !!
        b3(i) =( (B_minus_l(i) .le. lambda_left(i) .and. lambda_left(i) .le. B_plus_l(i)) .and. & !! Adding a point to left meets the requiremenst for DBI or PPI
               (B_minus_r(i) .le. lambda_right(i) .and. lambda_right(i) .le. B_plus_r(i)) .and. &   !! Adding a point to right meets the requiremenst for DBI or PPI
                abs(x(i)-x_left(i)) == abs(x_right(i)-x(i+1)) )
      enddo
      

      !$OMP SIMD PRIVATE(i)
      do i=1, n-1
        b4(i) = (abs(lambda_left(i)) < abs(lambda_right(i))) .and. (b3(i) .eqv. .true.) .and. &
                   (b2(i) .eqv. .false.) .and. (b1(i) .eqv. .false.)
        b5(i) = (abs(lambda_left(i)) >= abs(lambda_right(i))) .and. (b3(i) .eqv. .true.) .and. &
                   (b2(i) .eqv. .false.) .and. (b1(i) .eqv. .false.)
      enddo
      
      !$OMP SIMD PRIVATE(i)
      do i=1, n-1
          b6(i) = (B_minus_r(i) .le. lambda_right(i)) .and. (lambda_right(i) .le. B_plus_r(i)) .and. &   !! Adding a point to right meets the requiremenst for DBI or PPI
                    (b1(i) .eqv. .false.)  .and. (b2(i) .eqv. .false.) .and. (b3(i) .eqv. .false.)
      enddo
      
      !$OMP SIMD PRIVATE(i)
      do i=1, n-1
        b7(i) = (B_minus_l(i) .le. lambda_left(i)) .and. (lambda_left(i) .le. B_plus_l(i)) .and. & !! Adding a point to left meets the requiremenst for DBI or PPI
                (b1(i) .eqv. .false.) .and. (b2(i) .eqv. .false.) .and. (b3(i) .eqv. .false.) .and. (b6(i) .eqv. .false. )
      enddo
      

      !$OMP SIMD PRIVATE(i)
      do i=1, n-1
        bool(i) = (b1(i) .or. b4(i) .or. b7(i))
        bool2(i) = (b2(i) .or. b5(i) .or. b6(i))
        bool(i) = abs(bool(i))
        bool2(i) = abs(bool2(i))
      enddo
      

      !!$OMP SIMD 
      !do i=1, n-1
      !  bool(i) = abs(bool(i))
      !  bool2(i) = abs(bool2(i))
      !enddo
      
    endif

    !$OMP SIMD PRIVATE(i, si, ei, tmp_si, tmp_ei)
    do i=1, n-1
      si = f_si(i)
      ei = f_ei(i)
      tmp_si = max(f_si(i)-1,1)
      tmp_ei = min(f_ei(i)+1, n)
      !!
      lambda_new(i) = bool(i)*lambda_left(i) + bool2(i)*lambda_right(i) + &
                     (1-bool(i))*(1-bool2(i))*lambda_new(i)
      !!
      B_minus(i) = bool(i)*B_minus_l(i) + bool2(i)*B_minus_r(i) + &
                     (1-bool(i))*(1-bool2(i))*B_minus(i)
      !!
      B_plus(i) = bool(i)*B_plus_l(i) + bool2(i)*B_plus_r(i) + &
                     (1-bool(i))*(1-bool2(i))*B_plus(i)
      !!
      prod_deltaxx(i) = bool(i)*prod_deltax_left(i) + bool2(i)*prod_deltax_right(i) + &
                     (1-bool(i))*(1-bool2(i))*prod_deltaxx(i)
      !!
      xx(i) = bool(i)*x(tmp_si) + bool2(i)*x(tmp_ei) + &
              (1-bool(i))*(1-bool2(i))*xx(i)
      f_si(i) = bool(i)*tmp_si + bool2(i)*f_si(i) + &
                     (1-bool(i))*(1-bool2(i))*f_si(i)
      !!
      f_ei(i) = bool(i)*f_ei(i) + bool2(i)*tmp_ei + &
                     (1-bool(i))*(1-bool2(i))*f_ei(i)
    enddo
   
  enddo ! of j loop
   
  k =1
  do i=1,n-1
      
      si = f_si(i)
      ei = f_ei(i)
      !$OMP SIMD PRIVATE(j)
      do j=1, degree+1
        xval(j) = 0.0
        u(j) = 0.0
      enddo
      
      !$OMP SIMD PRIVATE(j)
      do j=1, ei-si+1
        u(j) = table(si, j)
        xval(j) = x(si+j-1)
      enddo
      
      !!** Building and evaluating Interpolant at xout points **!!
      !! - do while( x(i) <= xout(k) .and. xout(k) <= x(i+1) .and. k <= m )
      if( k <=m)then
        do while( x(i) <= xout(k) .and. xout(k) <= x(i+1) )
          call newtonPolyVal(xval, u, degree, xout(k), yout(k))
          k = k+1
          if(k > m) exit
        enddo
      endif
    enddo

end subroutine !!adaptiveInterpolation1D



subroutine adaptiveInterpolation2D_vec(x, y, nx, ny, v,  xout, yout, mx, my, vout, degree, interpolation_type, st, eps0, eps1)
!!!This routine adaptively build ia 2D tensor product interpoaltion based adaptiveInterpolation1D(...)
!! INPUT: 
!! nx: the number points in the 1D vector x.
!! ny: the number points in the 1D vector y.
!! mx: the number of points in the 1D vector xout.
!! my: the number of points in the 1D vector xout.
!! x: 1D mesh points of length nx used to build tensor product mesh. For i=1, ..., n-1 x_{i} <  x_{i+1}
!! y: 1D mesh points of length ny used to build tensor product mesh. For i=1, ..., n-1 y_{i} <  y_{i+1}
!! v: 2D array that have the data values associated with the tensor product mesh obtained from x and y
!! xout: 1D vector of length mx used to construct the output tensor product mesh.
!! yout: 1D vector of length my used to construct the output tensor product mesh.
!! Interpolation_type: used to determine the type of interpolation to be used to build interpolant.
!!   - interpolation_type=1: a data-bounded interpolant is built for each interpolant.
!!   - interpolation_type=2: a positivity-preserving interpolant is built for each interpolant.
!! degree: target polynomial degree and maximum polynomial degree used for each interval.
!! st (optional): used guide point selection process in cases when adding the next point to the 
!!   right or left both meet the requirements for positivity or datat-boundedness.
!!   - st=1 (default): the point with the smallest divided difference is added (ENO stencil).
!!   - st=2 the point to the left of current stencil is selected if the number of point to left
!!     of x_{i} is smaller than the number of points to right of x_{i} (i-si < ei-i). Similarly, 
!!     the point to the right is selected if the number of points to the right of x_{i} is smaller
!!     than the number of points to the left (i-si > ei-i). When both the number of points to right 
!!     and left are the same, the algorithm chooses the point with the smallest lambda.  
!!   - st=3 the point that is closest to the starting interval is chosen.
!! eps0 (optional): positive parameter use constrain the bound of the positive interpolant in intervals with no
!!   extremum detected.
!! eps1 (optional): positive parameter use constrain the bound of the positive interpolant in intervals with 
!!   extremum detected.
!!
!! OUTPUT:
!! vout: results of evaluating interpolants on tensor product mesh obtained from xout and yout.
!!
!!

  
  integer, intent(in)                   :: nx, ny
  integer, intent(in)                   :: mx, my
  integer, intent(in)                   :: degree
  integer, intent(in)                   :: interpolation_type 
  integer, intent(in), optional         :: st 
  real(kind=8), intent(in), optional    :: eps0, eps1
  real(kind=8), intent(in)              :: x(nx), y(ny), v(nx,ny), xout(mx), yout(my)
  real(kind=8), intent(out)             :: vout(mx, my)

  integer                               :: i, j
  integer                               :: sten
  real(kind=8)                          :: voutx(mx, ny)
  real(kind=8)                          :: eps2, eps3
 
  ! Set optional parameters
  if(present(st)) then
    sten = st;
  else
    sten = 3;
  endif
  if(present(eps0)) then
    eps2 = eps0;
  else
    eps2 = 0.01;
  endif
  if(present(eps1)) then
    eps3 = eps1;
  else
    eps3 = 1.0;
  endif

 

  !!** interpolate along x **!!
  do j=1,ny
    call adaptiveInterpolation1D_vec(x, v(:,j), nx, xout, voutx(:,j), mx, degree, interpolation_type, sten, eps2, eps3)
  enddo

  !!** interpolate along y **!!
  do i=1,mx
    call adaptiveInterpolation1D_vec(y, voutx(i,:), ny, yout, vout(i,:), my, degree, interpolation_type, sten, eps2, eps3)
  enddo


end subroutine

subroutine adaptiveInterpolation2D(x, y, nx, ny, v,  xout, yout, mx, my, vout, degree, interpolation_type, st, eps0, eps1)
!!!This routine adaptively build ia 2D tensor product interpoaltion based adaptiveInterpolation1D(...)
!! INPUT: 
!! nx: the number points in the 1D vector x.
!! ny: the number points in the 1D vector y.
!! mx: the number of points in the 1D vector xout.
!! my: the number of points in the 1D vector xout.
!! x: 1D mesh points of length nx used to build tensor product mesh. For i=1, ..., n-1 x_{i} <  x_{i+1}
!! y: 1D mesh points of length ny used to build tensor product mesh. For i=1, ..., n-1 y_{i} <  y_{i+1}
!! v: 2D array that have the data values associated with the tensor product mesh obtained from x and y
!! xout: 1D vector of length mx used to construct the output tensor product mesh.
!! yout: 1D vector of length my used to construct the output tensor product mesh.
!! Interpolation_type: used to determine the type of interpolation to be used to build interpolant.
!!   - interpolation_type=1: a data-bounded interpolant is built for each interpolant.
!!   - interpolation_type=2: a positivity-preserving interpolant is built for each interpolant.
!! degree: target polynomial degree and maximum polynomial degree used for each interval.
!! st (optional): used guide point selection process in cases when adding the next point to the 
!!   right or left both meet the requirements for positivity or datat-boundedness.
!!   - st=1 (default): the point with the smallest divided difference is added (ENO stencil).
!!   - st=2 the point to the left of current stencil is selected if the number of point to left
!!     of x_{i} is smaller than the number of points to right of x_{i} (i-si < ei-i). Similarly, 
!!     the point to the right is selected if the number of points to the right of x_{i} is smaller
!!     than the number of points to the left (i-si > ei-i). When both the number of points to right 
!!     and left are the same, the algorithm chooses the point with the smallest lambda.  
!!   - st=3 the point that is closest to the starting interval is chosen.
!! eps0 (optional): positive parameter use constrain the bound of the positive interpolant in intervals with no
!!   extremum detected.
!! eps1 (optional): positive parameter use constrain the bound of the positive interpolant in intervals with 
!!   extremum detected.
!!
!! OUTPUT:
!! vout: results of evaluating interpolants on tensor product mesh obtained from xout and yout.
!!
!!

  integer, intent(in)                   :: nx, ny
  integer, intent(in)                   :: mx, my
  integer, intent(in)                   :: degree
  integer, intent(in)                   :: interpolation_type 
  integer, intent(in), optional         :: st 
  real(kind=8), intent(in), optional    :: eps0, eps1
  real(kind=8), intent(in)              :: x(nx), y(ny), v(nx,ny), xout(mx), yout(my)
  real(kind=8), intent(out)             :: vout(mx, my)

  integer                               :: i, j
  integer                               :: sten
  real(kind=8)                          :: voutx(mx, ny), dx(3), dy(3)
  real(kind=8)                          :: eps2, eps3

  ! Set optional parameters
  if(present(st)) then
    sten = st;
  else
    sten = 3;
  endif
  if(present(eps0)) then
    eps2 = eps0;
  else
    eps2 = 0.01;
  endif
  if(present(eps1)) then
    eps3 = eps1;
  else
    eps3 = 1.0;
  endif


  !!** interpolate along x **!!
  do j=1,ny
    call adaptiveInterpolation1D(x, v(:,j), nx, xout, voutx(:,j), mx, degree, interpolation_type, sten, eps2, eps3)
  enddo

  !!** interpolate along y **!!
  do i=1,mx
    call adaptiveInterpolation1D(y, voutx(i,:), ny, yout, vout(i,:), my, degree, interpolation_type, sten, eps2, eps3)
  enddo

end subroutine


subroutine adaptiveInterpolation3D_vec(x, y, z, nx, ny, nz, v,  xout, yout, zout, mx, my, mz, vout, degree, &
                                   interpolation_type, st, eps0, eps1)
!!!This routine adaptively build ia 3D tensor product interpoaltion based adaptiveInterpolation1D(...)
!! INPUT: 
!! nx: the number points in the 1D vector x.
!! ny: the number points in the 1D vector y.
!! nz: the number points in the 1D vector z.
!! mx: the number of points in the 1D vector xout.
!! my: the number of points in the 1D vector yout.
!! mz: the number of points in the 1D vector zout.
!! x: 1D mesh points of length nx used to build tensor product mesh. For i=1, ..., n-1 x_{i} <  x_{i+1}
!! y: 1D mesh points of length ny used to build tensor product mesh. For i=1, ..., n-1 y_{i} <  y_{i+1}
!! z: 1D mesh points of length ny used to build tensor product mesh. For i=1, ..., n-1 z_{i} <  z_{i+1}
!! v: 3D array that have the data values associated with the tensor product mesh obtained from x, y, and z
!! xout: 1D vector of length mx used to construct the output tensor product mesh.
!! yout: 1D vector of length my used to construct the output tensor product mesh.
!! zout: 1D vector of length my used to construct the output tensor product mesh.
!! Interpolation_type: used to determine the type of interpolation to be used to build interpolant.
!!   - interpolation_type=1: a data-bounded interpolant is built for each interpolant.
!!   - interpolation_type=2: a positivity-preserving interpolant is built for each interpolant.
!! degree: target polynomial degree and maximum polynomial degree used for each interval.
!! st (optional): used guide point selection process in cases when adding the next point to the 
!!   right or left both meet the requirements for positivity or datat-boundedness.
!!   - st=1 (default): the point with the smallest divided difference is added (ENO stencil).
!!   - st=2 the point to the left of current stencil is selected if the number of point to left
!!     of x_{i} is smaller than the number of points to right of x_{i} (i-si < ei-i). Similarly, 
!!     the point to the right is selected if the number of points to the right of x_{i} is smaller
!!     than the number of points to the left (i-si > ei-i). When both the number of points to right 
!!     and left are the same, the algorithm chooses the point with the smallest lambda.  
!!   - st=3 the point that is closest to the starting interval is chosen.
!! eps0 (optional): positive parameter use constrain the bound of the positive interpolant in intervals with no
!!   extremum detected.
!! eps1 (optional): positive parameter use constrain the bound of the positive interpolant in intervals with 
!!   extremum detected.
!!
!! OUTPUT:
!! vout: results of evaluating interpolants on tensor product mesh obtained from xout and yout.
!!
!!

  integer, intent(in)                   :: interpolation_type
  integer, intent(in), optional                   :: st
  integer, intent(in)                   :: nx, ny, nz, mx, my, mz, degree
  real(kind=8), intent(in)              :: x(nx), y(ny), z(nz), v(nx,ny,nz), xout(mx), yout(my), zout(mz)
  real(kind=8), intent(in), optional             :: eps0, eps1
  real(kind=8), intent(out)             :: vout(mx,my,mz)

  integer                               :: i, j, k, ii, jj, kk
  integer                               :: sten
  real(kind=8)                          :: dx(3), dy(3), dz(3)
  integer                               :: degx(nx), degy(ny), degz(nz)
  real(kind=8)                          :: tmpin(max(nx, mx, ny, my, nz,mz))
  real(kind=8)                          :: tmpout(max(nx, mx, ny, my, nz,mz))
  real(kind=8)                          :: voutt(max(nx, mx), max(ny, my),max(nz,mz))
  real(kind=8)                          :: eps2, eps3

  !! set optional parameters 
  if(present(st)) then
    sten = st
  else
    sten = 3
  endif
  if(present(eps0)) then
    eps2 = eps0
  else
    eps2 = 0.01
  endif
  if(present(eps1)) then
    eps3 = eps1
  else
    eps3 = 1.0
  endif
 
  !!** interpolate along x **!!
  do k=1,nz
    do j=1,ny
      do ii=1, nx
        tmpin(ii) = v(ii,j,k)
      enddo
      call adaptiveInterpolation1D_vec(x, tmpin(1:nx), nx, xout, tmpout(1:mx), mx, degree, interpolation_type, sten, eps2, eps3)
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
      call adaptiveInterpolation1D_vec(y, tmpin(1:ny), ny, yout, tmpout(1:my), my, degree, interpolation_type, sten, eps2, eps3)
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
      call adaptiveInterpolation1D_vec(z, tmpin(1:nz), nz, zout, tmpout(1:mz), mz, degree, interpolation_type, sten, eps2, eps3)
      do kk=1, mz
        vout(i, j, kk) = tmpout(kk)
      enddo
    enddo
  enddo


end subroutine



subroutine adaptiveInterpolation3D(x, y, z, nx, ny, nz, v,  xout, yout, zout, mx, my, mz, vout, degree, &
                                   interpolation_type, st, eps0, eps1)
!!!This routine adaptively build ia 3D tensor product interpoaltion based adaptiveInterpolation1D(...)
!! INPUT: 
!! nx: the number points in the 1D vector x.
!! ny: the number points in the 1D vector y.
!! nz: the number points in the 1D vector z.
!! mx: the number of points in the 1D vector xout.
!! my: the number of points in the 1D vector yout.
!! mz: the number of points in the 1D vector zout.
!! x: 1D mesh points of length nx used to build tensor product mesh. For i=1, ..., n-1 x_{i} <  x_{i+1}
!! y: 1D mesh points of length ny used to build tensor product mesh. For i=1, ..., n-1 y_{i} <  y_{i+1}
!! z: 1D mesh points of length ny used to build tensor product mesh. For i=1, ..., n-1 z_{i} <  z_{i+1}
!! v: 3D array that have the data values associated with the tensor product mesh obtained from x, y, and z
!! xout: 1D vector of length mx used to construct the output tensor product mesh.
!! yout: 1D vector of length my used to construct the output tensor product mesh.
!! zout: 1D vector of length my used to construct the output tensor product mesh.
!! Interpolation_type: used to determine the type of interpolation to be used to build interpolant.
!!   - interpolation_type=1: a data-bounded interpolant is built for each interpolant.
!!   - interpolation_type=2: a positivity-preserving interpolant is built for each interpolant.
!! degree: target polynomial degree and maximum polynomial degree used for each interval.
!! st (optional): used guide point selection process in cases when adding the next point to the 
!!   right or left both meet the requirements for positivity or datat-boundedness.
!!   - st=1 (default): the point with the smallest divided difference is added (ENO stencil).
!!   - st=2 the point to the left of current stencil is selected if the number of point to left
!!     of x_{i} is smaller than the number of points to right of x_{i} (i-si < ei-i). Similarly, 
!!     the point to the right is selected if the number of points to the right of x_{i} is smaller
!!     than the number of points to the left (i-si > ei-i). When both the number of points to right 
!!     and left are the same, the algorithm chooses the point with the smallest lambda.  
!!   - st=3 the point that is closest to the starting interval is chosen.
!! eps0 (optional): positive parameter use constrain the bound of the positive interpolant in intervals with no
!!   extremum detected.
!! eps1 (optional): positive parameter use constrain the bound of the positive interpolant in intervals with 
!!   extremum detected.
!!
!! OUTPUT:
!! vout: results of evaluating interpolants on tensor product mesh obtained from xout and yout.
!!
!!

  integer, intent(in)                   :: interpolation_type
  integer, intent(in), optional                   :: st
  integer, intent(in)                   :: nx, ny, nz, mx, my, mz, degree
  real(kind=8), intent(in)              :: x(nx), y(ny), z(nz), v(nx,ny,nz), xout(mx), yout(my), zout(mz)
  real(kind=8), intent(in), optional             :: eps0, eps1
  real(kind=8), intent(out)             :: vout(mx,my,mz)

  integer                               :: i, j, k, ii, jj, kk
  integer                               :: sten
  real(kind=8)                          :: dx(3), dy(3), dz(3)
  integer                               :: degx(nx), degy(ny), degz(nz)
  real(kind=8)                          :: tmpin(max(nx, mx, ny, my, nz,mz))
  real(kind=8)                          :: tmpout(max(nx, mx, ny, my, nz,mz))
  real(kind=8)                          :: voutt(max(nx, mx), max(ny, my),max(nz,mz))
  real(kind=8)                          :: eps2, eps3

  !! set optional parameters 
  if(present(st)) then
    sten = st
  else
    sten = 3
  endif
  if(present(eps0)) then
    eps2 = eps0
  else
    eps2 = 0.01
  endif
  if(present(eps1)) then
    eps3 = eps1
  else
    eps3 = 1.0
  endif
 
  !!** interpolate along x **!!
  do k=1,nz
    do j=1,ny
      do ii=1, nx
        tmpin(ii) = v(ii,j,k)
      enddo
      call adaptiveInterpolation1D(x, tmpin(1:nx), nx, xout, tmpout(1:mx), mx, degree, interpolation_type, st, eps1, eps3)
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
      call adaptiveInterpolation1D(y, tmpin(1:ny), ny, yout, tmpout(1:my), my, degree, interpolation_type, sten, eps2, eps3)
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
      call adaptiveInterpolation1D(z, tmpin(1:nz), nz, zout, tmpout(1:mz), mz, degree, interpolation_type, sten, eps2, eps3)
      do kk=1, mz
        vout(i, j, kk) = tmpout(kk)
      enddo
    enddo
  enddo


end subroutine

end  !! end of module 


