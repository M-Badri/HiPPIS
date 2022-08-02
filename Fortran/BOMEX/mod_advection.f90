module mod_advection

  implicit none
  public :: flux, limiter

  contains

  subroutine flux(f, n, ff, kappa)
  !!! This subroutine computes third order flux
  !!
  !! INPUT: 
  !! f: data values to be used to compute fluxes
  !! n: number of data points in f
  !! kappa: 1 first order central, -1 second oder upwind, 1/3 third order upwind 
  !!
  !! OUTPUT:
  !! ff: flux at cell centered 
    
    integer, intent(in)         :: n
    integer                     :: i
    real(kind=8), intent(in)    :: f(n), kappa
    real(kind=8), intent(out)   :: ff(n-1)
    real(kind=8)                :: lim(n-1), tmp, r(n-1), du0, du1


    !!!--  kappa = 1 second order central, kappa = -1 second order upwind
    !!!--  kappa = 1/3 third order upwind biased
    !!kappa = 1.0/3.0
    !!!-- compute the ratio
    !!!-- r_{i+1/2} = (f_{i+1}-f_{i}) / (f_{i}-f_{i-1})
    !do i=2, n-2
    !  r(i) = (f(i)-f(i+1))/ (f(i+1)-f(i+2) + 1e-30)
    !enddo
   
    !!!-- Compute the limter 
    !!! K(r) = (1-k)/2 + (1+k)/2*r 
    !!! phi_{i+1/2} = max( 0, min(2r, min(2, K(r_{i+1/2})) ) )
    !do i=2, n-2
    !   tmp = (1.0-kappa)/2.0 + (1.0+kappa)/2.0*r(i)
    !   lim(i) = max(0.0, min(2.0*r(i), min(2.0, tmp)))
    !   !tmp = 0.25 + 0.5 * r(i)
    !   !lim(i) = max(min(min(tmp, 4.0), min(tmp, 2*r(i))), 1e-8)
    !enddo

    do i=2, n-2
       du0 = f(i+1)-f(i+2)
       du1 = f(i)-f(i+1)
       call limiter(8, du0, du1, kappa, tmp)
       lim(i) = tmp
    enddo
    !!-- compute the 2 fluxes 
    !! ff_{i+1/2} = f_{i+1} + 0.5 phi_{i+1/2}(f_{i+1}-f_{i+2})
    do i=2, n-2
      ff(i) = f(i+1) + 0.5 * lim(i)*(f(i+1)-f(i+2))
      if(abs(ff(i)) < 1.0e-20 .and. ff(i) .ne. 0.0)then
        ff(i) = 0.0
      endif
    enddo

  end subroutine

  subroutine limiter( i, du0, du1, kappa, phi)
  !! limiter function
  !! value depend on value of i
  !! i = 1 first order  
  !! i = 2 second order central
  !! i = 3 van leer harnonic
  !! i = 4 third order unlimited
  !! i = 4 CCCT

  integer, intent(in)         :: i
  real(kind=8), intent(in)    :: du0, du1, kappa
  real(kind=8), intent(out)   :: phi
  real(kind=8)                :: tol, r, tmp



  tol = 1.0e-8
  if(i==1) then
     phi = 0.0
  elseif(i == 2) then
     phi = 1.0
  endif 
  if(i == 3) then
     r = du1 / ( du0 + 1.0e-26)
     phi = ( r + abs(r) ) / ( 1.0 + abs(r))
  elseif (i == 4) then
     r = du1 / ( du0 + 1.0e-26)
     phi = 0.25 + 0.75*r
  endif
  if (i == 5) then
     r = du1 / ( du0 + 1.0e-26)
     phi = 0.25 + 0.75*r
      if(phi > 4.0) then
         phi = 4.0
     endif
     if (phi > (2.0*r)) then
        phi = 2.0*r
     endif
     if(phi < tol) then 
        phi = 0.0
     endif
   endif
  if (i == 6) then
     r = du1 / ( du0 + 1.0e-26)
     phi = r;
      if (phi > 1.0) then
         phi = 1.0
     endif
     if (phi < tol) then 
        phi = 0.0
     endif
   endif
  if (i == 7) then
     r = du1 / ( du0 + 1.0e-26)
     phi = 0.5 + 0.5*r
      if (phi > 2.0) then
         phi = 2.0
     endif
     if (phi > (2*r)) then
        phi = 2*r;
     endif
      if (phi > 2.0) then
         phi = 2.0
     endif
     if (phi < tol) then
        phi = 0.0
     endif
   endif
   if(i == 8)then
     r = du1 / ( du0 + 1.0e-26)
     tmp = (1.0-kappa)/2.0 + (1.0+kappa)/2.0*r
     phi = max(0.0, min(2.0*r, min(2.0, tmp)))
   endif

  end subroutine  !! limiter
end module
