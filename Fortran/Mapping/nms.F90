subroutine pchev ( n, x, f, d, nval, xval, fval, dval, ierr )

!*****************************************************************************80
!
!! PCHEV evaluates a piecewise cubic Hermite or spline function.
!
!  Discussion:
!
!    PCHEV evaluates the function and first derivative of a piecewise
!    cubic Hermite or spline function at an array of points XVAL.
!
!    The evaluation will be most efficient if the elements of XVAL are
!    increasing relative to X; that is, for all J <= K,
!      X(I) <= XVAL(J)
!    implies
!      X(I) <= XVAL(K).
!
!    If any of the XVAL are outside the interval [X(1),X(N)],
!    values are extrapolated from the nearest extreme cubic,
!    and a warning error is returned.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    Fred Fritsch,
!    Piecewise Cubic Hermite Interpolation Package, Final Specifications,
!    Lawrence Livermore National Laboratory,
!    Computer Documentation UCID-30194, August 1982.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.  
!    N must be at least 2.
!
!    Input, real ( kind=dp ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind=dp ) F(N), the function values.  F(I) is the value
!    corresponding to X(I).
!
!    Input, real ( kind=dp ) D(N), the derivative values.  D(i) is the value
!    corresponding to X(I).
!
!    Input, integer ( kind = 4 ) NVAL, the number of points at which the 
!    functions are to be evaluated.
!
!    Input, real ( kind=dp ) XVAL(NVAL), the points at which the functions
!    are to be evaluated.
!
!    Output, real ( kind=dp ) FVAL(NVAL), the values of the cubic Hermite
!    function at XVAL.
!
!    Output, real ( kind=dp ) DVAL(NVAL), the derivatives of the cubic
!    Hermite function at XVAL.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    positive, means that extrapolation was performed at IERR points.
!    -1, if N < 2.
!    -3, if the X array is not strictly increasing.
!    -4, if NVAL < 1.
!    -5, if an error has occurred in CHFDV.
!
  implicit none

  integer, parameter  :: dp=kind(0.d0)                 

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nval

  real ( kind=dp ) d(n)
  real ( kind=dp ) dval(nval)
  real ( kind=dp ) f(n)
  real ( kind=dp ) fval(nval)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ), save :: incfd = 1
  logical, save :: skip = .true.
  real ( kind=dp ) x(n)
  real ( kind=dp ) xval(nval)

  call pchfd ( n, x, f, d, incfd, skip, nval, xval, fval, dval, ierr )

  return
end

subroutine pchez ( n, x, f, d, spline, wk, lwk, ierr )

!*****************************************************************************80
!
!! PCHEZ carries out easy to use spline or cubic Hermite interpolation.
!
!  Discussion:
!
!    This routine sets derivatives for spline (two continuous derivatives)
!    or Hermite cubic (one continuous derivative) interpolation.
!    Spline interpolation is smoother, but may not "look" right if the
!    data contains both "steep" and "flat" sections.  Hermite cubics
!    can produce a "visually pleasing" and monotone interpolant to
!    monotone data.
!
!    This routine is an easy to use driver for the PCHIP routines.
!    Various boundary conditions are set to default values by PCHEZ.
!    Many other choices are available in the subroutines PCHIC,
!    PCHIM and PCHSP.
!
!    Use PCHEV to evaluate the resulting function and its derivative.
!
!    If SPLINE is TRUE, the interpolating spline satisfies the default
!    "not-a-knot" boundary condition, with a continuous third derivative
!    at X(2) and X(N-1).
!
!    If SPLINE is FALSE, the interpolating Hermite cubic will be monotone
!    if the input data is monotone.  Boundary conditions are computed from
!    the derivative of a local quadratic unless this alters monotonicity.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    Fred Fritsch, Judy Butland,
!    A Method for Constructing Local Monotone Piecewise
!    Cubic Interpolants,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 5, Number 2, 1984, pages 300-304.
!
!    Carl deBoor,
!    A Practical Guide to Splines, Chapter IV,
!    Springer-Verlag,
!    New York, 1978.
!
!    Fred Fritsch,
!    Piecewise Cubic Hermite Interpolation Package, Final Specifications,
!    Lawrence Livermore National Laboratory,
!    Computer Documentation UCID-30194, August 1982.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.  
!    N must be at least 2.
!
!    Input, real ( kind=dp ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind=dp ) F(N), the function values.  F(I) is the value
!    corresponding to X(I).
!
!    Output, real ( kind=dp ) D(N), the derivative values at the data points.
!
!    Input, logical SPLINE, specifies if the interpolant is to be a spline
!    with two continuous derivatives (SPLINE is TRUE), or a Hermite cubic
!    interpolant with one continuous derivative (SPLINE is FALSE).
!
!    Workspace, real ( kind=dp ) WK(LWK), required only if SPLINE is TRUE.
!
!    Input, integer ( kind = 4 ) LWK, the length of the work array WK, which 
!    must be at least 2*N.  However, WK is not needed if SPLINE is FALSE,
!    and in this case LWK is arbitrary.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    positive, can only occur when SPLINE is FALSE,  means that there were
!      IERR switches in the direction of monotonicity.  When SPLINE is
!      FALSE, PCHEZ guarantees that if the input data is monotone, the
!      interpolant will be too.  This warning is to alert you to the fact
!      that the input data was not monotone.
!    -1, if N < 2.
!    -3, if the X array is not strictly increasing.
!    -7, if LWK is less than 2*N and SPLINE is TRUE.
!
  implicit none

  integer, parameter  :: dp=kind(0.d0)                 

  integer ( kind = 4 ) lwk
  integer ( kind = 4 ) n

  real ( kind=dp ) d(n)
  real ( kind=dp ) f(n)
  integer ( kind = 4 ), save, dimension ( 2 ) :: ic = (/ 0, 0 /)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ), parameter :: incfd = 1
  logical spline
  real ( kind=dp ) vc(2)
  real ( kind=dp ) wk(lwk)
  real ( kind=dp ) x(n)

  if ( spline ) then
    call pchsp ( ic, vc, n, x, f, d, incfd, wk, lwk, ierr )
  else
    call pchim ( n, x, f, d, incfd, ierr )
  end if

  return
end
subroutine pchfd ( n, x, f, d, incfd, skip, ne, xe, fe, de, ierr )

!*****************************************************************************80
!
!! PCHFD evaluates a piecewise cubic Hermite function.
!
!  Discsussion:
!
!    PCHFD evaluates a piecewise cubic Hermite function and its first
!    derivative at an array of points.  PCHFD may be used by itself
!    for Hermite interpolation, or as an evaluator for PCHIM
!    or PCHIC.
!
!    PCHFD evaluates the cubic Hermite function and its first derivative
!    at the points XE.
!
!    If only function values are required, use PCHFE instead.
!
!    To provide compatibility with PCHIM and PCHIC, includes an
!    increment between successive values of the F and D arrays.
!
!    Most of the coding between the call to CHFDV and the end of
!    the IR loop could be eliminated if it were permissible to
!    assume that XE is ordered relative to X.
!
!    CHFDV does not assume that X1 is less than X2.  Thus, it would
!    be possible to write a version of PCHFD that assumes a strictly
!    decreasing X array by simply running the IR loop backwards
!    and reversing the order of appropriate tests.
!
!    The present code has a minor bug, which I have decided is not
!    worth the effort that would be required to fix it.
!    If XE contains points in [X(N-1),X(N)], followed by points less than
!    X(N-1), followed by points greater than X(N), the extrapolation points
!    will be counted (at least) twice in the total returned in IERR.
!
!    The evaluation will be most efficient if the elements of XE are
!    increasing relative to X; that is, for all J <= K,
!      X(I) <= XE(J)
!    implies
!      X(I) <= XE(K).
!
!    If any of the XE are outside the interval [X(1),X(N)],
!    values are extrapolated from the nearest extreme cubic,
!    and a warning error is returned.
!
!  Modified:
!
!    13 August 2005
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
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
!    Input, integer ( kind = 4 ) N, the number of data points.  
!    N must be at least 2.
!
!    Input, real ( kind=dp ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind=dp ) F(INCFD,N), the function values.
!    F(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Input, real ( kind=dp ) D(INCFD,N), the derivative values.
!    D(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Input, integer ( kind = 4 ) INCFD, increment between successive values in 
!    F and D.
!
!    Input/output, logical SKIP, controls whether data validity checks
!    should be made.  Setting the input value to FALSE will skip the checks.
!    On output with 0 <= IERR, SKIP will be set to TRUE.
!
!    Input, integer ( kind = 4 ) NE, the number of evaluation points.
!
!    Input, real ( kind=dp ) XE(NE), points at which the function is
!    to be evaluated.
!
!    Output, real ( kind=dp ) FE(NE), the values of the cubic Hermite
!    function at XE.
!
!    Output, real ( kind=dp ) DE(NE), the derivative of the cubic
!    Hermite function at XE.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    positive, means that extrapolation was performed at IERR points.
!    -1, if N < 2.
!    -2, if INCFD < 1.
!    -3, if the X array is not strictly increasing.
!    -4, if NE < 1.
!    -5, if an error has occurred in the lower-level routine CHFDV.
!
  implicit none
  integer, parameter  :: dp=kind(0.d0)                 

  integer ( kind = 4 ) incfd
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ne

  real ( kind=dp ) d(incfd,n)
  real ( kind=dp ) de(ne)
  real ( kind=dp ) f(incfd,n)
  real ( kind=dp ) fe(ne)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierc
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j_first
  integer ( kind = 4 ) j_new
  integer ( kind = 4 ) j_save
  integer ( kind = 4 ) next(2)
  integer ( kind = 4 ) nj
  logical skip
  real ( kind=dp ) x(n)
  real ( kind=dp ) xe(ne)
!
!  Check arguments.
!
  if ( .not. skip ) then

    if ( n < 2 ) then
      ierr = -1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PCHFD - Fatal error!'
      write ( *, '(a)' ) '  Number of data points less than 2.'
      return
    end if

    if ( incfd < 1 ) then
      ierr = -2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PCHFD - Fatal error!'
      write ( *, '(a)' ) '  Increment less than 1.'
      return
    end if

    do i = 2, n
      if ( x(i) <= x(i-1) ) then
        ierr = -3
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PCHFD - Fatal error!'
        write ( *, '(a)' ) '  X array not strictly increasing.'
        return
      end if
    end do

  end if

  if ( ne < 1 ) then
    ierr = -4
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PCHFD - Fatal error!'
    write ( *, '(a)' ) '  Number of evaluation points less than 1.'
    return
  end if

  ierr = 0
  skip = .true.
!
!  Loop over intervals.
!  The interval index is IL = IR - 1.
!  The interval is X(IL) <= X < X(IR).
!
  j_first = 1
  ir = 2

  do
!
!  Skip out of loop if have processed all evaluation points.
!
    if ( ne < j_first ) then
      exit
    end if
!
!  Locate all points in interval.
!
    j_save = ne + 1

    do j = j_first, ne
      if ( x(ir) <= xe(j) ) then
        j_save = j
        if ( ir == n ) then
          j_save = ne + 1
        end if
        exit
      end if
    end do
!
!  Have located first point beyond interval.
!
    j = j_save

    nj = j - j_first
!
!  Skip evaluation if no points in interval.
!
    if ( nj /= 0 ) then
!
!  Evaluate cubic at XE(J_FIRST:J-1).
!
      call chfdv ( x(ir-1), x(ir), f(1,ir-1), f(1,ir), d(1,ir-1), d(1,ir), &
        nj, xe(j_first), fe(j_first), de(j_first), next, ierc )

      if ( ierc < 0 ) then
        ierr = -5
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PCHFD - Fatal error!'
        write ( *, '(a)' ) '  Error return from CHFDV.'
        return
      end if
!
!  In the current set of XE points, there are NEXT(2) to the right of X(IR).
!
      if ( next(2) /= 0 ) then

        if ( ir < n ) then
          ierr = -5
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'PCHFD - Fatal error!'
          write ( *, '(a)' ) '  IR < N.'
          return
        end if
!
!  These are actually extrapolation points.
!
        ierr = ierr + next(2)

      end if
!
!  In the current set of XE points, there are NEXT(1) to the left of X(IR-1).
!
      if ( next(1) /= 0 ) then
!
!  These are actually extrapolation points.
!
        if ( ir <= 2 ) then
          ierr = ierr + next(1)
!
!  XE is not ordered relative to X, so must adjust evaluation interval.
!
!  First, locate first point to left of X(IR-1).
!
        else

          j_new = -1

          do i = j_first, j-1
            if ( xe(i) < x(ir-1) ) then
              j_new = i
              exit
            end if
          end do

          if ( j_new == -1 ) then
            ierr = -5
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PCHFD - Fatal error!'
            write ( *, '(a)' ) '  Could not bracket the data point.'
            return
          end if
!
!  Reset J.  This will be the new J_FIRST.
!
          j = j_new
!
!  Now find out how far to back up in the X array.
!
          do i = 1, ir-1
            if ( xe(j) < x(i) ) then
              exit
            end if
          end do
!
!  At this point, either XE(J) < X(1) or X(i-1) <= XE(J) < X(I) .
!
!  Reset IR, recognizing that it will be incremented before cycling.
!
          ir = max ( 1, i-1 )

        end if

      end if

      j_first = j

    end if

    ir = ir + 1

    if ( n < ir ) then
      exit
    end if

  end do

  return
end

subroutine pchim ( n, x, f, d, incfd, ierr )

!*****************************************************************************80
!
!! PCHIM sets derivatives for a piecewise cubic Hermite interpolant.
!
!  Discussion:
!
!    The routine set derivatives needed to determine a monotone piecewise
!    cubic Hermite interpolant to given data.
!
!    The interpolant will have an extremum at each point where
!    monotonicity switches direction.  See PCHIC if user control is desired
!    over boundary or switch conditions.
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the F and D arrays.
!
!    The resulting piecewise cubic Hermite function may be evaluated
!    by PCHFE or PCHFD.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    Fred Fritsch, Judy Butland,
!    A Method for Constructing Local Monotone Piecewise
!    Cubic Interpolants,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 5, Number 2, 1984, pages 300-304.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.  
!    N must be at least 2.
!
!    Input, real ( kind=dp ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind=dp ) F(INCFD,N), dependent variable values to be
!    interpolated.  F(1+(I-1)*INCFD) is the value corresponding to X(I).
!    PCHIM is designed for monotonic data, but it will work for any F-array.
!    It will force extrema at points where monotonicity switches direction.
!    If some other treatment of switch points is desired, PCHIC should be
!    used instead.
!
!    Output, real ( kind=dp ) D(INCFD,N), the derivative values at the
!    data points.  If the data are monotonic, these values will determine
!    a monotone cubic Hermite function.  The value corresponding to X(I)
!    is stored in D(1+(I-1)*INCFD).
!
!    Input, integer ( kind = 4 ) INCFD, the increment between successive
!    values in F and D.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    positive, means IERR switches in the direction of monotonicity
!    were detected.
!    -1, if N < 2.
!    -2, if INCFD < 1.
!    -3, if the X array is not strictly increasing.
!
  implicit none

  integer, parameter  :: dp=kind(0.d0)                 
  integer ( kind = 4 ) incfd
  integer ( kind = 4 ) n

  real ( kind=dp ) d(incfd,n)
  real ( kind=dp ) del1
  real ( kind=dp ) del2
  real ( kind=dp ) dmax
  real ( kind=dp ) dmin
  real ( kind=dp ) drat1
  real ( kind=dp ) drat2
  real ( kind=dp ) dsave
  real ( kind=dp ) f(incfd,n)
  real ( kind=dp ) h1
  real ( kind=dp ) h2
  real ( kind=dp ) hsum
  real ( kind=dp ) hsumt3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) nless1
  real ( kind=dp ) pchst
  real ( kind=dp ) temp
  real ( kind=dp ) w1
  real ( kind=dp ) w2
  !!real ( kind=dp ) x(n)
  real ( kind=dp ) x(n)
!
!  Check the arguments.
!
  if ( n < 2 ) then
    ierr = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PCHIM - Fatal error!'
    write ( *, '(a)' ) '  Number of data points less than 2.'
    return
  end if

  if ( incfd < 1 ) then
    ierr = -2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PCHIM - Fatal error!'
    write ( *, '(a)' ) '  Increment less than 1.'
    return
  end if

  do i = 2, n
    if ( x(i) <= x(i-1) ) then
      ierr = -3
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PCHIM - Fatal error!'
      write ( *, '(a)' ) '  X array not strictly increasing.'
        print*, 'i=', i, 'x(i)=', x(i), 'x(i+1)=', x(i+1)
      return
    end if
  end do

  ierr = 0
  nless1 = n - 1
  h1 = x(2) - x(1)
  del1 = ( f(1,2) - f(1,1) ) / h1
  dsave = del1
!
!  Special case N=2, use linear interpolation.
!
  if ( n == 2 ) then
    d(1,1) = del1
    d(1,n) = del1
    return
  end if
!
!  Normal case, 3 <= N.
!
  h2 = x(3) - x(2)
  del2 = ( f(1,3) - f(1,2) ) / h2
!
!  Set D(1) via non-centered three point formula, adjusted to be
!  shape preserving.
!
  hsum = h1 + h2
  w1 = ( h1 + hsum ) / hsum
  w2 = -h1 / hsum
  d(1,1) = w1 * del1 + w2 * del2

  if ( pchst ( d(1,1), del1 ) <= 0.0_dp ) then

    d(1,1) = 0.0_dp
!
!  Need do this check only if monotonicity switches.
!
  else if ( pchst ( del1, del2 ) < 0.0_dp ) then

     dmax = 3.0_dp * del1

     if ( abs ( dmax ) < abs ( d(1,1) ) ) then
       d(1,1) = dmax
     end if

  end if
!
!  Loop through interior points.
!
  do i = 2, nless1

    if ( 2 < i ) then
      h1 = h2
      h2 = x(i+1) - x(i)
      hsum = h1 + h2
      del1 = del2
      del2 = ( f(1,i+1) - f(1,i) ) / h2
    end if
!
!  Set D(I)=0 unless data are strictly monotonic.
!
    d(1,i) = 0.0_dp

    temp = pchst ( del1, del2 )

    if ( temp < 0.0_dp ) then

      ierr = ierr + 1
      dsave = del2
!
!  Count number of changes in direction of monotonicity.
!
    else if ( temp == 0.0_dp ) then

      if ( del2 /= 0.0_dp ) then
        if ( pchst ( dsave, del2 ) < 0.0_dp ) then
          ierr = ierr + 1
        end if
        dsave = del2
      end if
!
!  Use Brodlie modification of Butland formula.
!
    else

      hsumt3 = 3.0_dp * hsum
      w1 = ( hsum + h1 ) / hsumt3
      w2 = ( hsum + h2 ) / hsumt3
      dmax = max ( abs ( del1 ), abs ( del2 ) )
      dmin = min ( abs ( del1 ), abs ( del2 ) )
      drat1 = del1 / dmax
      drat2 = del2 / dmax
      d(1,i) = dmin / ( w1 * drat1 + w2 * drat2 )

    end if

  end do
!
!  Set D(N) via non-centered three point formula, adjusted to be
!  shape preserving.
!
  w1 = -h2 / hsum
  w2 = ( h2 + hsum ) / hsum
  d(1,n) = w1 * del1 + w2 * del2

  if ( pchst ( d(1,n), del2 ) <= 0.0_dp ) then
    d(1,n) = 0.0_dp
  else if ( pchst ( del1, del2 ) < 0.0_dp ) then
!
!  Need do this check only if monotonicity switches.
!
    dmax = 3.0_dp * del2

    if ( abs ( dmax ) < abs ( d(1,n) ) ) then
      d(1,n) = dmax
    end if

  end if

  return
end

subroutine pchsp ( ic, vc, n, x, f, d, incfd, wk, nwk, ierr )

!*****************************************************************************80
!
!! PCHSP sets derivatives needed for Hermite cubic spline interpolant.
!
!  Description:
!
!    PCHSP sets derivatives needed to determine the Hermite representation
!    of the cubic spline interpolant to given data, with specified boundary
!    conditions.
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the F and D arrays.
!
!    The resulting piecewise cubic Hermite function may be evaluated
!    by PCHFE or PCHFD.
!
!    This is a modified version of Carl de Boor's cubic spline routine CUBSPL.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Carl de Boor,
!    A Practical Guide to Splines,
!    Springer-Verlag (new york, 1978).
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IC(2), specifies desired boundary conditions:
!    IC(1) = IBEG, desired condition at beginning of data.
!    0, to set D(1) so that the third derivative is continuous at X(2).
!      This is the "not a knot" condition provided by de Boor's cubic spline
!      routine CUBSPL, and is the default boundary condition here.
!    1, if first derivative at X(1) is given in VC(1).
!    2, if second derivative at X(1) is given in VC(1).
!    3, to use the 3-point difference formula for D(1).
!      Reverts to the default boundary condition if N < 3.
!    4, to use the 4-point difference formula for D(1).
!      Reverts to the default boundary condition if N < 4.
!    For the "natural" boundary condition, use ibeg=2 and vc(1)=0.
!    IC(2) = IEND, desired condition at end of data.
!    IEND may take on the same values as IBEG, but applied to derivative at
!    X(N).  In case IEND = 1 or 2, the value is given in VC(2).
!
!    Input, real ( kind=dp ) VC(2), specifies desired boundary values,
!    as indicated above.  VC(1) need be set only if IC(1) = 1 or 2.
!    VC(2) need be set only if IC(2) = 1 or 2.
!
!    Input, integer ( kind = 4 ) N, the number of data points.  N must be
!    at least 2.
!
!    Input, real ( kind=dp ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind=dp ) F(INCFD,N), the dependent values to be
!    interpolated.  F(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Output, real ( kind=dp ) D(INCFD,N), the derivative values at the
!    data points.  These values will determine the cubic spline interpolant
!    with the requested boundary conditions.  The value corresponding to
!    X(I) is stored in D(1+(I-1)*INCFD).
!
!    Input, integer ( kind = 4 ) INCFD, increment between successive values
!    in F and D.
!
!    Workspace, real ( kind=dp ) WK(NWK).
!
!    Input, integer ( kind = 4 ) NWK, the size of WK, which must be 
!    at least 2 * N.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    -1, if N < 2.
!    -2, if INCFD < 1.
!    -3, if the X array is not strictly increasing.
!    -4, if IBEG < 0 or 4 < IBEG.
!    -5, if IEND < 0 or 4 < IEND.
!    -6, if both of the above are true.
!    -7, if NWK is too small.
!    -8, in case of trouble solving the linear system
!        for the interior derivative values.
!
  implicit none

  integer, parameter  :: dp=kind(0.d0)                 

  integer ( kind = 4 ) incfd
  integer ( kind = 4 ) n

  real ( kind=dp ) d(incfd,n)
  real ( kind=dp ) f(incfd,n)
  real ( kind=dp ) g
  integer ( kind = 4 ) ibeg
  integer ( kind = 4 ) ic(2)
  integer ( kind = 4 ) iend
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) index
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nwk
  real ( kind=dp ) pchdf
  real ( kind=dp ) stemp(3)
  real ( kind=dp ) vc(2)
  real ( kind=dp ) wk(2,n)
  real ( kind=dp ) x(n)
  real ( kind=dp ) xtemp(4)

  if ( n < 2 ) then
    ierr = -1
    call xerror ('pchsp -- number of data points less than two', ierr, 1)
    return
  end if

  if ( incfd < 1 ) then
    ierr = -2
    call xerror ('pchsp -- increment less than one', ierr, 1)
    return
  end if

  do j = 2, n
    if ( x(j) <= x(j-1) ) then
      ierr = -3
      call xerror ('pchsp -- x-array not strictly increasing', ierr, 1)
      return
    end if
  end do

  ibeg = ic(1)
  iend = ic(2)
  ierr = 0

  if ( ibeg < 0 .or. 4 < ibeg ) then
    ierr = ierr - 1
  end if

  if ( iend < 0 .or. 4 < iend ) then
    ierr = ierr - 2
  end if

  if ( ierr < 0 ) then
    go to 5004
  end if
!
!  Function definition is ok -- go on.
!
  if ( nwk < 2 * n ) then
    go to 5007
  end if
!
!  Compute first differences of X sequence and store in wk(1,.). also,
!  compute first divided difference of data and store in wk(2,.).
!
  do j = 2, n
    wk(1,j) = x(j) - x(j-1)
    wk(2,j) = ( f(1,j) - f(1,j-1) ) / wk(1,j)
  end do
!
!  Set to default boundary conditions if N is too small.
!
  if ( n < ibeg ) then
    ibeg = 0
  end if

  if ( n < iend ) then
    iend = 0
  end if
!
!  Set up for boundary conditions.
!
  if ( ibeg == 1 .or. ibeg == 2 ) then
     d(1,1) = vc(1)
  else if ( 2 < ibeg ) then
!
!  Pick up first IBEG points, in reverse order.
!
     do j = 1, ibeg
       index = ibeg - j + 1
       xtemp(j) = x(index)
       if ( j < ibeg ) then
         stemp(j) = wk(2,index)
       end if
     end do

     d(1,1) = pchdf ( ibeg, xtemp, stemp, ierr )
     if ( ierr /= 0 ) then
       go to 5009
     end if

     ibeg = 1
  end if

  if ( iend == 1 .or. iend == 2 ) then
     d(1,n) = vc(2)
  else if ( 2 < iend ) then
!
!  Pick up last IEND points.
!
     do j = 1, iend
       index = n - iend + j
       xtemp(j) = x(index)
       if ( j < iend ) then
         stemp(j) = wk(2,index+1)
       end if
     end do

     d(1,n) = pchdf ( iend, xtemp, stemp, ierr )

     if ( ierr /= 0 ) then
       go to 5009
     end if

     iend = 1

  end if
!
!  Begin coding from cubspl
!
!  A tridiagonal linear system for the unknown slopes S(1:N) of
!  F at X(1:N) is generated and then solved by Gauss elimination,
!  with s(j) ending up in d(1,j), all j.
!  wk(1,.) and wk(2,.) are used for temporary storage.
!
!  Construct first equation from first boundary condition, of the form
!    wk(2,1) * s(1) + wk(1,1) * s(2) = D(1,1)
!
  if ( ibeg == 0 ) then
!
!  No condition at left end and N = 2.
!
     if ( n == 2 ) then
        wk(2,1) = 1.0_dp
        wk(1,1) = 1.0_dp
        d(1,1) = 2.0_dp * wk(2,2)
!
!  Not-a-knot condition at left end and 2 < N.
!
     else
        wk(2,1) = wk(1,3)
        wk(1,1) = wk(1,2) + wk(1,3)
        d(1,1) =(( wk(1,2) + 2.0_dp * wk(1,1) ) * wk(2,2) * wk(1,3) &
                             + wk(1,2)**2 * wk(2,3) ) / wk(1,1)
     end if
  else if ( ibeg == 1 ) then
!
!  Slope prescribed at left end.
!
     wk(2,1) = 1.0_dp
     wk(1,1) = 0.0_dp
  else
!
!  Second derivative prescribed at left end.
!
     wk(2,1) = 2.0_dp
     wk(1,1) = 1.0_dp
     d(1,1) = 3.0_dp * wk(2,2) - 0.5_dp * wk(1,2) * d(1,1)
  end if
!
!  If there are interior knots, generate the corresponding equations and
!  carry out the forward pass of Gauss elimination, after which the J-th
!  equation reads
!
!    wk(2,j) * s(j) + wk(1,j) * s(j+1) = d(1,j).
!
  if ( 1 < n-1 ) then
    do j = 2, n-1
        if ( wk(2,j-1) == 0.0_dp ) then
          go to 5008
        end if
        g = -wk(1,j+1) / wk(2,j-1)
        d(1,j) = g * d(1,j-1) + 3.0_dp &
          * ( wk(1,j) * wk(2,j+1) + wk(1,j+1) * wk(2,j) )
        wk(2,j) = g * wk(1,j-1) + 2.0_dp * ( wk(1,j) + wk(1,j+1) )
    end do
  end if
!
!  Construct last equation from second boundary condition, of the form
!
!    (-g * wk(2,n-1)) * s(n-1) + wk(2,n) * s(n) = d(1,n)
!
!  If slope is prescribed at right end, one can go directly to back-
!  substitution, since arrays happen to be set up just right for it
!  at this point.
!
  if ( iend == 1 ) then
    go to 30
  end if

  if ( iend == 0 ) then
     if ( n == 2 .and. ibeg == 0 ) then
!
!  Not-a-knot at right endpoint and at left endpoint and N = 2.
!
        d(1,2) = wk(2,2)
        go to 30
     else if ( n == 2 .or. ( n == 3 .and. ibeg == 0 ) ) then
!
!  Either ( N = 3 and not-a-knot also at left) or (N=2 and *not*
!  not-a-knot at left end point).
!
        d(1,n) = 2.0_dp * wk(2,n)
        wk(2,n) = 1.0_dp
        if ( wk(2,n-1) == 0.0_dp ) then
          go to 5008
        end if
        g = -1.0_dp / wk(2,n-1)
     else
!
!  Not-a-knot and 3 <= N, and either 3 < N or also not-a-
!  knot at left end point.
!
        g = wk(1,n-1) + wk(1,n)
!
!  Do not need to check following denominators (x-differences).
!
        d(1,n) = ( ( wk(1,n) + 2.0_dp * g ) * wk(2,n) * wk(1,n-1) &
          + wk(1,n)**2 * ( f(1,n-1) - f(1,n-2) ) / wk(1,n-1) ) / g
        if ( wk(2,n-1) == 0.0_dp ) then
          go to 5008
        end if
        g = -g / wk(2,n-1)
        wk(2,n) = wk(1,n-1)
     end if
  else
!
!  Second derivative prescribed at right endpoint.
!
     d(1,n) = 3.0_dp *wk(2,n) + 0.5_dp * wk(1,n) * d(1,n)
     wk(2,n) = 2.0_dp
     if ( wk(2,n-1) == 0.0_dp ) then
       go to 5008
     end if
     g = -1.0_dp / wk(2,n-1)
  end if
!
!  Complete forward pass of Gauss elimination.
!
  wk(2,n) = g * wk(1,n-1) + wk(2,n)

  if ( wk(2,n) == 0.0_dp ) then
    go to 5008
  end if

  d(1,n) = ( g * d(1,n-1) + d(1,n) ) / wk(2,n)
!
!  Carry out back substitution.
!
   30 continue

  do j = n-1, 1, -1
    if ( wk(2,j) == 0.0_dp ) then
      go to 5008
    end if
    d(1,j) = ( d(1,j) - wk(1,j) * d(1,j+1) ) / wk(2,j)
  end do

  return
!
!  error returns.
!
 5004 continue
!
!  ic out of range return.
!
  ierr = ierr - 3
  call xerror ('pchsp -- ic out of range', ierr, 1)
  return

 5007 continue
!
!  nwk too small return.
!
  ierr = -7
  call xerror ('pchsp -- work array too small', ierr, 1)
  return

 5008 continue
!
!  singular system.
!  theoretically, this can only occur if successive x-values
!  are equal, which should already have been caught (ierr=-3).
!
  ierr = -8
  call xerror ('pchsp -- singular linear system', ierr, 1)
  return
!
 5009 continue
!
!  error return from pchdf.
!  this case should never occur.
!
  ierr = -9
  call xerror ('pchsp -- error return from pchdf', ierr, 1)

  return
end
subroutine xerprt ( messg, nmessg )

!*****************************************************************************80
!
!! XERPRT prints a message on each file indicated by xgetua.
!
!  Modified:
!
!    05 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, character ( len = * ) MESSG, the message to be printed.
!
!    Input, integer ( kind = 4 ) NMESSG, the actual number of characters 
!    in MESSG.
!
  implicit none

  !- Added 
  integer ( kind = 4 ) ichar
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) kunit
  integer ( kind = 4 ) last
  integer ( kind = 4 ) lenmes
  integer ( kind = 4 ) lun(5)
  character ( len = * ) messg
  integer ( kind = 4 ) nmessg
  integer ( kind = 4 ) nunit


  ! added to avoid warning for unsused variable
  nmessg = nmessg
 
!
!  Obtain unit numbers and write line to each unit
!
  call xgetua ( lun, nunit )

  lenmes = len ( messg )

  do kunit = 1, nunit

     iunit = lun(kunit)

     do ichar = 1, lenmes, 72
        last = min ( ichar+71 , lenmes )
        if ( iunit == 0 ) then
          write (*,'(1x,a)') messg(ichar:last)
        else
          write (iunit,'(1x,a)') messg(ichar:last)
        end if
    end do

  end do

  return
end

subroutine xerror ( messg, nerr, level )

!*****************************************************************************80
!
!! XERROR processes a diagnostic error message.
!
!  Discussion:
!
!    XERROR processes a diagnostic message, in a manner
!    determined by the value of level and the current value
!    of the library error control flag, kontrl.
!    See subroutine xsetf for details.
!
!  Example:
!
!    call xerror('smooth -- num was zero.',1,2)
!
!    call xerror('integ  -- less than full accuracy achieved.',2,1)
!
!    call xerror( &
!      'rooter -- actual zero of f found before interval fully collapsed.',3,0)
!
!    call xerror('exp    -- underflows being set to zero.',1,-1)
!
!  Modified:
!
!    13 August 2005
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, character ( len = * ) MESSG, the message to be processed,
!    containing no more than 72 characters.
!
!    Input, integer ( kind = 4 ) NERR, the error number associated with this 
!    message.  NERR must not be zero.
!
!    Input, integer ( kind = 4 ) LEVEL, the error category.
!    2 means this is an unconditionally fatal error.
!    1 means this is a recoverable error.  (i.e., it is
!      non-fatal if XSETF has been appropriately called.)
!    0 means this is a warning message only.
!    -1 means this is a warning message which is to be printed at most once,
!      regardless of how many times this call is executed.
!
  implicit none
  integer, parameter  :: dp=kind(0.d0)                 

  integer ( kind = 4 ) level
  character ( len = * ) messg
  integer ( kind = 4 ) nerr
  integer ( kind = 4 ) nmessg

  nmessg = len ( messg )

  call xerrwv ( messg, nmessg, nerr, level, 0, 0, 0, 0, 0.0_dp, 0.0_dp )

  return
end
subroutine xerrwv ( messg, nmessg, nerr, level, ni, i1, i2, nr, r1, r2 )

!*****************************************************************************80
!
!! XERRWV processes an error message that includes numeric information.
!
!  Discussion:
!
!    XERRWV processes a diagnostic message, in a manner
!    determined by the value of level and the current value
!    of the library error control flag, kontrl.
!    (see subroutine xsetf for details.)
!    in addition, up to two integer ( kind = 4 ) values and two real
!    values may be printed along with the message.
!
!  Example:
!
!    call xerrwv ( 'smooth -- num (=i1) was zero.',29,1,2,1,num,0,0,0.,0.)
!
!    call xerrwv ( &
!      'quadxy -- requested error (r1) less than minimum(r2).', &
!      54,77,1,0,0,0,2,errreq,errmin)
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, character ( len = * ) MESSG, the message to be processed.
!
!    Input, integer ( kind = 4 ) NMESSG, the number of characters in MESSG.
!
!    Input, integer ( kind = 4 ) NERR, the error number associated with
!    this message.  NERR must not be zero.
!
!    Input, integer ( kind = 4 ) LEVEL, the error category.
!    2 means this is an unconditionally fatal error.
!    1 means this is a recoverable error.  (i.e., it is
!      non-fatal if xsetf has been appropriately called.)
!    0 means this is a warning message only.
!    -1 means this is a warning message which is to be printed at most
!      once, regardless of how many times this call is executed.
!
!    Input, integer ( kind = 4 ) NI, the number of integer values to be
!    printed. (0 to 2)
!
!    Input, integer ( kind = 4 ) I1, I2, the first and second integer values.
!
!    Input, integer ( kind = 4 ) NR, the number of real values to be
!    printed. (0 to 2)
!
!    Input, real ( kind=dp ) R1, R2, the first and second real values.
!
  implicit none

  integer, parameter  :: dp=kind(0.d0)                 

  character ( len = 37 ) form
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i1mach
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ifatal
  integer ( kind = 4 ) isizei
  integer ( kind = 4 ) isizef
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j4save
  integer ( kind = 4 ) junk
  integer ( kind = 4 ) kdummy
  integer ( kind = 4 ) kount
  integer ( kind = 4 ) kunit
  integer ( kind = 4 ) lerr
  integer ( kind = 4 ) level
  character ( len = 20 ) lfirst
  integer ( kind = 4 ) lkntrl
  integer ( kind = 4 ) llevel
  integer ( kind = 4 ) lmessg
  integer ( kind = 4 ) lun(5)
  integer ( kind = 4 ) maxmes
  character ( len = * ) messg
  integer ( kind = 4 ) mkntrl
  integer ( kind = 4 ) nerr
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nmessg
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) nunit
  real ( kind=dp ) r1
  real ( kind=dp ) r2
!
!  Get flags
!
  lkntrl = j4save ( 2, 0, .false. )
  maxmes = j4save ( 4, 0, .false. )
!
!  Check for valid input
!
  if ( 0 < nmessg .and. nerr /= 0 .and. -1 <= level .and. level <= 2 ) then
    go to 10
  end if

    if ( 0 < lkntrl ) then
      call xerprt('fatal error in...',17)
    end if

    call xerprt( 'XERROR -- invalid input', 23 )

    if ( 0 < lkntrl ) then
      call xerprt('job abort due to fatal error.',29)
    end if

    if ( 0 < lkntrl ) then
      call xersav ( ' ', 0, 0, 0, kdummy )
    end if

    call xerabt('XERROR -- invalid input',23)
    return

   10 continue
!
!  Record the message.
!
  junk = j4save(1,nerr,.true.)
  call xersav ( messg, nmessg, nerr, level, kount )
!
!  Let user override
!
  lfirst = messg
  lmessg = nmessg
  lerr = nerr
  llevel = level
  call xerctl(lfirst,lmessg,lerr,llevel,lkntrl)
!
!  Reset to original values.
!
  lmessg = nmessg
  lerr = nerr
  llevel = level
  lkntrl = max ( -2, min ( 2, lkntrl ) )
  mkntrl = abs ( lkntrl )
!
!  Decide whether to print message
!
  if ( llevel < 2 .and. lkntrl == 0 ) then
    go to 100
  end if

  if (((llevel == (-1)) .and. ( min ( 1, maxmes ) < kount ) ) &
    .or.((llevel == 0) .and. ( maxmes < kount )) &
    .or.((llevel == 1) .and. ( maxmes < kount ).and.(mkntrl==1) ) &
    .or.((llevel == 2) .and. ( max ( 1, maxmes ) < kount ) ) ) then
    go to 100
  end if

  if ( 0 < lkntrl ) then

    call xerprt(' ',1)

    if ( llevel == -1 ) then

      call xerprt &
      ( 'warning message...this message will only be printed once.',57)

    end if

    if ( llevel == 0 ) then
      call xerprt ( 'warning in...', 13 )
    else if ( llevel == 1 ) then
      call xerprt ( 'recoverable error in...', 23 )
    else if ( llevel == 2 ) then
      call xerprt ( 'fatal error in...', 17 )
    end if

  end if
!
!  Message
!
     call xerprt(messg,lmessg)
     call xgetua(lun,nunit)
     isizei = log10 ( real ( i1mach(9), dp ) ) + 1.0_dp
     isizef = log10 ( real ( i1mach(10), dp )**i1mach(14) ) + 1.0_dp

     do kunit = 1, nunit

        iunit = lun(kunit)

        do i = 1, min ( ni, 2 )
           write (form,21) i,isizei
   21          format ('(11x,21hin above message, i',i1,'=,i',i2,')   ')
           if ( iunit == 0 ) then
             if (i == 1) write (*,form) i1
             if (i == 2) write (*,form) i2
           else
             if (i == 1) write (iunit,form) i1
             if (i == 2) write (iunit,form) i2
           end if
        end do

        do i = 1, min ( nr, 2 )
           write (form,23) i,isizef+10,isizef
   23          format ('(11x,21hin above message, r',i1,'=,e',i2,'.',i2,')')
           if ( iunit == 0 ) then
             if ( i == 1 ) write (*,form) r1
             if ( i == 2 ) write (*,form) r2
           else
             if (i == 1) write (iunit,form) r1
             if (i == 2) write (iunit,form) r2
           end if
        end do

        if ( lkntrl <= 0 ) then
          go to 40
        end if
!
!  error number
!
           if ( iunit == 0 ) then
             write(*,30) lerr
           else
             write (iunit,30) lerr
           end if
!! old submission !!   30          format (15h error number =,i10)
   30          format ('error number =',i10)
   40       continue

     end do
!
!  Traceback
!
  100 continue
  ifatal = 0
  if ((llevel == 2).or.((llevel==1) .and. (mkntrl==2))) then
    ifatal = 1
  end if
!
!  quit here if message is not fatal
!
  if ( ifatal <= 0 ) then
    return
  end if

  if ( lkntrl <= 0 .or. max ( 1, maxmes ) < kount ) then
    go to 120
  end if
!
!  Print reason for abort
!
     if ( llevel == 1 ) then
       call xerprt ('job abort due to unrecovered error.',35)
     end if

     if ( llevel == 2 ) then
       call xerprt('job abort due to fatal error.',29)
     end if
!
!  Print error summary
!
     call xersav ( ' ', -1, 0, 0, kdummy )

  120 continue
!
!  Abort
!
  if ( llevel == 2 .and. max ( 1, maxmes ) < kount ) then
    lmessg = 0
  end if

  call xerabt ( messg, lmessg )

  return
end
subroutine xersav ( messg, nmessg, nerr, level, icount )

!*****************************************************************************80
!
!! XERSAV records that an error occurred.
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, character ( len = * ) MESSG, as in XERROR.
!
!    Input, integer ( kind = 4 ) NMESSG, as in XERROR, except that, when 
!    NMESSG = 0, the tables will be dumped and cleared; and when NMESSG < 0,
!    the tables will be dumped, but not cleared.
!
!    Input, integer ( kind = 4 ) NERR, as in XERROR.
!
!    Input, integer ( kind = 4 ) LEVEL, as in XERROR.
!
!    Output, integer ( kind = 4 ) ICOUNT, the number of times this message has
!    been seen, or zero if the table has overflowed and
!    does not contain this message specifically.
!    when nmessg=0, icount will not be altered.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1mach
  integer ( kind = 4 ) icount
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ), save, dimension ( 10 ) :: kount = (/ &
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
  integer ( kind = 4 ), save :: kountx = 0
  integer ( kind = 4 ) kunit
  integer ( kind = 4 ) level
  integer ( kind = 4 ), save, dimension ( 10 ) :: levtab
  integer ( kind = 4 ) lun(5)
  character ( len = 20 ) mes
  character ( len = * ) messg
  character ( len = 20 ), save, dimension ( 10 ) :: mestab
  integer ( kind = 4 ) nerr
  integer ( kind = 4 ), save, dimension ( 10 ) :: nertab
  integer ( kind = 4 ) nmessg
  integer ( kind = 4 ) nunit
!
!  Dump the table
!
  if ( nmessg <= 0 ) then

     if ( kount(1) == 0 ) then
       return
     end if
!
!  Print to each unit
!
     call xgetua ( lun, nunit )

     do kunit = 1, nunit

       iunit = lun(kunit)

       if ( iunit == 0 ) then
         iunit = i1mach(4)
       end if
!
!  Print table header
!
       write ( iunit, 10 )
   10  format ('0          error message summary'/ &
       ' message start             nerr     level     count')
!
!  print body of table
!
        do i = 1, 10
          if ( kount(i) == 0 ) then
            exit
          end if
          write (iunit,15) mestab(i),nertab(i),levtab(i),kount(i)
   15     format (1x,a20,3i10)
        end do
!
!  Print number of other errors
!
        if ( kountx /= 0 ) then
          write (iunit,40) kountx
        end if

!! old submission !!   40       format (41h0other errors not individually tabulated=,i10)
   40       format ('errors not individually tabulated=',i10)
        write ( iunit, '(a)' ) ' '
     end do
!
!  Clear the error tables
!
    if ( nmessg == 0 ) then
      kount(1:10) = 0
      kountx = 0
    end if

    return

  end if
!
!  process a message...
!  search for this message, or else an empty slot for this messg,
!  or else determine that the error table is full.
!
  mes = messg

  do i = 1, 10

    ii = i

    if ( kount(i) == 0 ) then
      mestab(ii) = mes
      nertab(ii) = nerr
      levtab(ii) = level
      kount(ii)  = 1
      icount = 1
      return
    end if

    if ( mes /= mestab(i) ) then
      go to 90
    end if

    if (nerr /= nertab(i) ) then
      go to 90
    end if

    if (level /= levtab(i) ) then
      go to 90
    end if

    go to 100

90  continue

  end do
!
!  The table is full.
!
  kountx = kountx + 1
  icount = 1
  return
!
!  Message found in table
!
  100    continue

     kount(ii) = kount(ii) + 1
     icount = kount(ii)

  return
end
function j4save ( iwhich, ivalue, iset )

!*****************************************************************************80
!
!! J4SAVE saves variables needed by the library error handling routines.
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IWHICH, the index of the item desired.
!    1, the current error number.
!    2, the current error control flag.
!    3, the current unit number to which error messages are sent.
!       (0 means use standard.)
!    4, the maximum times any message is printed (as set by xermax).
!    5, the number of units to which each error message is written.
!    6, the 2nd unit for error messages.
!    7, the 3rd unit for error messages.
!    8, the 4th unit for error messages.
!    9, the 5th unit for error messages.
!
!    Input, integer ( kind = 4 ) IVALUE, the value to be set for the IWHICH-th 
!    parameter, if ISET is TRUE.
!
!    Input, logical ISET.
!    TRUE: the IWHICH-th parameter will be given the value, IVALUE.
!
!    Output, integer ( kind = 4 ) J4SAVE, the old value of the IWHICH-th 
!    parameter.
!
  implicit none

  integer ( kind = 4 ), save, dimension ( 9 ) :: iparam = (/ &
    0, 2, 0, 10, 1, 0, 0, 0, 0 /)
  logical iset
  integer ( kind = 4 ) ivalue
  integer ( kind = 4 ) iwhich
  integer ( kind = 4 ) j4save

  j4save = iparam(iwhich)

  if ( iset ) then
    iparam(iwhich) = ivalue
  end if

  return
end
subroutine xerctl ( messg1, nmessg, nerr, level, kontrl )

!*****************************************************************************80
!
!! XERCTL allows user control over handling of individual errors.
!
!  Discussion:
!
!    Allows user control over handling of individual errors.
!    Just after each message is recorded, but before it is
!    processed any further (i.e., before it is printed or
!    a decision to abort is made), a call is made to XERCTL.
!    If the user has provided his own version of XERCTL, he
!    can then override the value of KONTRL used in processing
!    this message by redefining its value.
!
!    KONTRL may be set to any value from -2 to 2.
!    The meanings for KONTRL are the same as in XSETF, except
!    that the value of KONTRL changes only for this message.
!    If KONTRL is set to a value outside the range from -2 to 2,
!    it will be moved back into that range.
!
!  Modified:
!
!    05 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, character ( len = * ) MESSG1, the first word (only) of the error
!    message.
!
!    Input, integer ( kind = 4 ) NMESSG, same as in the call to XERROR 
!    or XERRWV.
!
!    Input, integer ( kind = 4 ) NERR, same as in the call to XERROR or XERRWV.
!
!    Input, integer ( kind = 4 ) LEVEL, same as in the call to XERROR or XERRWV.
!
!    Input/output, integer ( kind = 4 ) KONTRL.  On input, the current value of 
!    the control flag as set by a call to XSETF.  On output, the new value of
!    kontrl.  If KONTRL is not defined, it will remain at its original value.
!    This changed value affects only the current occurrence of the current
!    message.
!
  implicit none

  integer ( kind = 4 ) kontrl
  integer ( kind = 4 ) level
  character ( len = * ) messg1
  integer ( kind = 4 ) nerr
  integer ( kind = 4 ) nmessg

  kontrl = kontrl
  level = level 
  messg1 = messg1
  nerr = nerr
  nmessg = nmessg
  
  return
end
subroutine xgetua ( iunita, n )

!*****************************************************************************80
!
!! XGETUA returns the unit number(s) to which error messages are being sent.
!
!  Discussion:
!
!    XGETUA may be called to determine the unit number or numbers to which
!    error messages are being sent.  These unit numbers may have been set
!    by a call to XSETUN, or a call to XSETUA, or may be a default value.
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNITA(N),  an array unit numbers,
!    A value of zero refers to the default unit, as defined by the
!    I1MACH machine constant routine.  Only IUNITA(1),..., IUNITA(N) are
!    defined by XGETUA.  The values of IUNITA(N+1),..., IUNITA(5) are
!    not defined (for N < 5) or altered in any way by XGETUA.
!
!    Output, integer ( kind = 4 ) N, the number of units to which copies of the
!    error messages are being sent.  N will be in the range from 1 to 5.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) index
  integer ( kind = 4 ) iunita(5)
  integer ( kind = 4 ) j4save
  integer ( kind = 4 ) n

  n = j4save ( 5, 0, .false. )

  do i = 1, n

    index = i+4
    if ( i == 1 ) then
      index = 3
    end if

    iunita(i) = j4save ( index, 0, .false. )

  end do

  return
end
function i1mach ( i )

!*****************************************************************************80
!
!! I1MACH returns integer ( kind = 4 ) machine constants.
!
!  I/O unit numbers.
!
!    I1MACH(1) = the standard input unit.
!    I1MACH(2) = the standard output unit.
!    I1MACH(3) = the standard punch unit.
!    I1MACH(4) = the standard error message unit.
!
!  Words.
!
!    I1MACH(5) = the number of bits per integer ( kind = 4 ) storage unit.
!    I1MACH(6) = the number of characters per integer ( kind = 4 ) storage unit.
!
!  Integers.
!
!  Assume integer ( kind = 4 )s are represented in the S digit base A form:
!
!  Sign * (X(S-1)*A**(S-1) + ... + X(1)*A + X(0))
!  where 0<=X(I)<A for I=0 to S-1.
!
!    I1MACH(7) = A, the base.
!    I1MACH(8) = S, the number of base A digits.
!    I1MACH(9) = A**S-1, the largest integer ( kind = 4 ).
!
!  Floating point numbers
!
!  Assume floating point numbers are represented in the T digit base B form:
!
!    Sign * (B**E) * ((X(1)/B) + ... + (X(T)/B**T) )
!
!  where 0<=X(I)<B for I = 1 to T, 0<X(1) and EMIN<=E<=EMAX
!
!    I1MACH(10) = B, the base.
!
!  Single precision
!
!    I1MACH(11) = T, the number of base B digits.
!    I1MACH(12) = EMIN, the smallest exponent E.
!    I1MACH(13) = EMAX, the largest exponent E.
!
!  Double precision
!
!    I1MACH(14) = T, the number of base B digits.
!    I1MACH(15) = EMIN, the smallest exponent E.
!    I1MACH(16) = EMAX, the largest exponent E.
!
!  To alter this function for a particular environment, the desired set of DATA
!  statements should be activated by removing the C from column 1.  On rare
!  machines, a STATIC statement may need to be added, but probably more systems
!  prohibit than require it.
!
!  Also, the values of I1MACH(1) through I1MACH(4) should be checked for
!  consistency with the local operating system.  For FORTRAN 77, you may wish
!  to adjust the data statement so imach(6) is set to 1, and then to comment
!  out the executable test on I.EQ.6 below.
!
!  For IEEE-arithmetic machines (binary standard), the first set of constants
!  below should be appropriate, except perhaps for IMACH(1) - IMACH(4).
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1mach
  integer ( kind = 4 ) imach(16)
  integer ( kind = 4 ) output

  equivalence (imach(4),output)
!
!  IEEE arithmetic machines, such as the ATT 3B series, Motorola
!  68000 based machines such as the SUN 3 and ATT PC 7300, and
!  8087 based micros such asthe IBM PC and ATT 6300.
!
   data imach( 1) /    5 /
   data imach( 2) /    6 /
   data imach( 3) /    7 /
   data imach( 4) /    6 /
   data imach( 5) /   32 /
   data imach( 6) /    4 /
   data imach( 7) /    2 /
   data imach( 8) /   31 /
   data imach( 9) / 2147483647 /
   data imach(10) /    2 /
   data imach(11) /   24 /
   data imach(12) / -125 /
   data imach(13) /  128 /
   data imach(14) /   53 /
   data imach(15) / -1021 /
   data imach(16) /  1024 /
!
!  ALLIANT FX/8 UNIX FORTRAN compiler.
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     6 /
!      data imach( 4) /     0 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    32 /
!      data imach( 9) /2147483647/
!      data imach(10) /     2 /
!      data imach(11) /    24 /
!      data imach(12) /  -126 /
!      data imach(13) /   128 /
!      data imach(14) /    53 /
!      data imach(15) / -1022 /
!      data imach(16) /  1024 /
!
!  AMDAHL machines.
!
!      data imach( 1) /   5 /
!      data imach( 2) /   6 /
!      data imach( 3) /   7 /
!      data imach( 4) /   6 /
!      data imach( 5) /  32 /
!      data imach( 6) /   4 /
!      data imach( 7) /   2 /
!      data imach( 8) /  31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /  16 /
!      data imach(11) /   6 /
!      data imach(12) / -64 /
!      data imach(13) /  63 /
!      data imach(14) /  14 /
!      data imach(15) / -64 /
!      data imach(16) /  63 /
!
!  BURROUGHS 1700 system.
!
!      data imach( 1) /    7 /
!      data imach( 2) /    2 /
!      data imach( 3) /    2 /
!      data imach( 4) /    2 /
!      data imach( 5) /   36 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   33 /
!      data imach( 9) / Z1FFFFFFFF /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -256 /
!      data imach(13) /  255 /
!      data imach(14) /   60 /
!      data imach(15) / -256 /
!      data imach(16) /  255 /
!
!  BURROUGHS 5700 system.
!
!      data imach( 1) /   5 /
!      data imach( 2) /   6 /
!      data imach( 3) /   7 /
!      data imach( 4) /   6 /
!      data imach( 5) /  48 /
!      data imach( 6) /   6 /
!      data imach( 7) /   2 /
!      data imach( 8) /  39 /
!      data imach( 9) / O0007777777777777 /
!      data imach(10) /   8 /
!      data imach(11) /  13 /
!      data imach(12) / -50 /
!      data imach(13) /  76 /
!      data imach(14) /  26 /
!      data imach(15) / -50 /
!      data imach(16) /  76 /
!
!  BURROUGHS 6700/7700 systems.
!
!      data imach( 1) /   5 /
!      data imach( 2) /   6 /
!      data imach( 3) /   7 /
!      data imach( 4) /   6 /
!      data imach( 5) /  48 /
!      data imach( 6) /   6 /
!      data imach( 7) /   2 /
!      data imach( 8) /  39 /
!      data imach( 9) / O0007777777777777 /
!      data imach(10) /   8 /
!      data imach(11) /  13 /
!      data imach(12) / -50 /
!      data imach(13) /  76 /
!      data imach(14) /  26 /
!      data imach(15) / -32754 /
!      data imach(16) /  32780 /
!
!  CDC CYBER 170/180 series using NOS
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   60 /
!      data imach( 6) /   10 /
!      data imach( 7) /    2 /
!      data imach( 8) /   48 /
!      data imach( 9) / O"00007777777777777777" /
!      data imach(10) /    2 /
!      data imach(11) /   48 /
!      data imach(12) / -974 /
!      data imach(13) / 1070 /
!      data imach(14) /   96 /
!      data imach(15) / -927 /
!      data imach(16) / 1070 /
!
!  CDC CYBER 170/180 series using NOS/VE
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     7 /
!      data imach( 4) /     6 /
!      data imach( 5) /    64 /
!      data imach( 6) /     8 /
!      data imach( 7) /     2 /
!      data imach( 8) /    63 /
!      data imach( 9) / 9223372036854775807 /
!      data imach(10) /     2 /
!      data imach(11) /    47 /
!      data imach(12) / -4095 /
!      data imach(13) /  4094 /
!      data imach(14) /    94 /
!      data imach(15) / -4095 /
!      data imach(16) /  4094 /
!
!  CDC CYBER 200 series
!
!      data imach( 1) /      5 /
!      data imach( 2) /      6 /
!      data imach( 3) /      7 /
!      data imach( 4) /      6 /
!      data imach( 5) /     64 /
!      data imach( 6) /      8 /
!      data imach( 7) /      2 /
!      data imach( 8) /     47 /
!      data imach( 9) / X'00007FFFFFFFFFFF' /
!      data imach(10) /      2 /
!      data imach(11) /     47 /
!      data imach(12) / -28625 /
!      data imach(13) /  28718 /
!      data imach(14) /     94 /
!      data imach(15) / -28625 /
!      data imach(16) /  28718 /
!
!  CDC 6000/7000 series using FTN4.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   60 /
!      data imach( 6) /   10 /
!      data imach( 7) /    2 /
!      data imach( 8) /   48 /
!      data imach( 9) / 00007777777777777777B /
!      data imach(10) /    2 /
!      data imach(11) /   47 /
!      data imach(12) / -929 /
!      data imach(13) / 1070 /
!      data imach(14) /   94 /
!      data imach(15) / -929 /
!      data imach(16) / 1069 /
!
!  CDC 6000/7000 series using FTN5.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   60 /
!      data imach( 6) /   10 /
!      data imach( 7) /    2 /
!      data imach( 8) /   48 /
!      data imach( 9) / O"00007777777777777777" /
!      data imach(10) /    2 /
!      data imach(11) /   47 /
!      data imach(12) / -929 /
!      data imach(13) / 1070 /
!      data imach(14) /   94 /
!      data imach(15) / -929 /
!      data imach(16) / 1069 /
!
!  CONVEX C-1.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   32 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -128 /
!      data imach(13) /  127 /
!      data imach(14) /   53 /
!      data imach(15) /-1024 /
!      data imach(16) / 1023 /
!
!  CONVEX C-120 (native mode) without -R8 option
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    0 /
!      data imach( 4) /    6 /
!      data imach( 5) /   32 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -127 /
!      data imach(13) /  127 /
!      data imach(14) /   53 /
!      data imach(15) / -1023 /
!      data imach(16) /  1023 /
!
!  CONVEX C-120 (native mode) with -R8 option
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     0 /
!      data imach( 4) /     6 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /     2 /
!      data imach(11) /    53 /
!      data imach(12) / -1023 /
!      data imach(13) /  1023 /
!      data imach(14) /    53 /
!      data imach(15) / -1023 /
!      data imach(16) /  1023 /
!
!  CONVEX C-120 (IEEE mode) without -R8 option
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    0 /
!      data imach( 4) /    6 /
!      data imach( 5) /   32 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -125 /
!      data imach(13) /  128 /
!      data imach(14) /   53 /
!      data imach(15) / -1021 /
!      data imach(16) /  1024 /
!
!  CONVEX C-120 (IEEE mode) with -R8 option
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     0 /
!      data imach( 4) /     6 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /     2 /
!      data imach(11) /    53 /
!      data imach(12) / -1021 /
!      data imach(13) /  1024 /
!      data imach(14) /    53 /
!      data imach(15) / -1021 /
!      data imach(16) /  1024 /
!
!  CRAY 1, 2, XMP and YMP.
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /   102 /
!      data imach( 4) /     6 /
!      data imach( 5) /    64 /
!      data imach( 6) /     8 /
!      data imach( 7) /     2 /
!      data imach( 8) /    63 /
!      data imach( 9) /  777777777777777777777B /
!      data imach(10) /     2 /
!      data imach(11) /    47 /
!      data imach(12) / -8189 /
!      data imach(13) /  8190 /
!      data imach(14) /    94 /
!      data imach(15) / -8099 /
!      data imach(16) /  8190 /
!
!  DATA GENERAL ECLIPSE S/200.
!
!      data imach( 1) /   11 /
!      data imach( 2) /   12 /
!      data imach( 3) /    8 /
!      data imach( 4) /   10 /
!      data imach( 5) /   16 /
!      data imach( 6) /    2 /
!      data imach( 7) /    2 /
!      data imach( 8) /   15 /
!      data imach( 9) /32767 /
!      data imach(10) /   16 /
!      data imach(11) /    6 /
!      data imach(12) /  -64 /
!      data imach(13) /   63 /
!      data imach(14) /   14 /
!      data imach(15) /  -64 /
!      data imach(16) /   63 /
!
!  ELXSI 6400
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     6 /
!      data imach( 4) /     6 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    32 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /     2 /
!      data imach(11) /    24 /
!      data imach(12) /  -126 /
!      data imach(13) /   127 /
!      data imach(14) /    53 /
!      data imach(15) / -1022 /
!      data imach(16) /  1023 /
!
!  HARRIS 220
!
!      data imach( 1) /       5 /
!      data imach( 2) /       6 /
!      data imach( 3) /       0 /
!      data imach( 4) /       6 /
!      data imach( 5) /      24 /
!      data imach( 6) /       3 /
!      data imach( 7) /       2 /
!      data imach( 8) /      23 /
!      data imach( 9) / 8388607 /
!      data imach(10) /       2 /
!      data imach(11) /      23 /
!      data imach(12) /    -127 /
!      data imach(13) /     127 /
!      data imach(14) /      38 /
!      data imach(15) /    -127 /
!      data imach(16) /     127 /
!
!  HARRIS SLASH 6 and SLASH 7.
!
!      data imach( 1) /       5 /
!      data imach( 2) /       6 /
!      data imach( 3) /       0 /
!      data imach( 4) /       6 /
!      data imach( 5) /      24 /
!      data imach( 6) /       3 /
!      data imach( 7) /       2 /
!      data imach( 8) /      23 /
!      data imach( 9) / 8388607 /
!      data imach(10) /       2 /
!      data imach(11) /      23 /
!      data imach(12) /    -127 /
!      data imach(13) /     127 /
!      data imach(14) /      38 /
!      data imach(15) /    -127 /
!      data imach(16) /     127 /
!
!  HONEYWELL DPS 8/70 and 600/6000 series.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /   43 /
!      data imach( 4) /    6 /
!      data imach( 5) /   36 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   35 /
!      data imach( 9) / O377777777777 /
!      data imach(10) /    2 /
!      data imach(11) /   27 /
!      data imach(12) / -127 /
!      data imach(13) /  127 /
!      data imach(14) /   63 /
!      data imach(15) / -127 /
!      data imach(16) /  127 /
!
!  HP 2100, 3 word double precision option with FTN4
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    4 /
!      data imach( 4) /    1 /
!      data imach( 5) /   16 /
!      data imach( 6) /    2 /
!      data imach( 7) /    2 /
!      data imach( 8) /   15 /
!      data imach( 9) / 32767 /
!      data imach(10) /    2 /
!      data imach(11) /   23 /
!      data imach(12) / -128 /
!      data imach(13) /  127 /
!      data imach(14) /   39 /
!      data imach(15) / -128 /
!      data imach(16) /  127 /
!
!  HP 2100, 4 word double precision option with FTN4
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    4 /
!      data imach( 4) /    1 /
!      data imach( 5) /   16 /
!      data imach( 6) /    2 /
!      data imach( 7) /    2 /
!      data imach( 8) /   15 /
!      data imach( 9) / 32767 /
!      data imach(10) /    2 /
!      data imach(11) /   23 /
!      data imach(12) / -128 /
!      data imach(13) /  127 /
!      data imach(14) /   55 /
!      data imach(15) / -128 /
!      data imach(16) /  127 /
!
!  HP 9000
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     6 /
!      data imach( 4) /     7 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    32 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /     2 /
!      data imach(11) /    24 /
!      data imach(12) /  -126 /
!      data imach(13) /   127 /
!      data imach(14) /    53 /
!      data imach(15) / -1015 /
!      data imach(16) /  1017 /
!
!  IBM 360/370 series, XEROX SIGMA 5/7/9, SEL systems 85/86, PERKIN ELMER 3230,
!  and PERKIN ELMER (INTERDATA) 3230.
!
!      data imach( 1) /   5 /
!      data imach( 2) /   6 /
!      data imach( 3) /   7 /
!      data imach( 4) /   6 /
!      data imach( 5) /  32 /
!      data imach( 6) /   4 /
!      data imach( 7) /   2 /
!      data imach( 8) /  31 /
!      data imach( 9) / Z7FFFFFFF /
!      data imach(10) /  16 /
!      data imach(11) /   6 /
!      data imach(12) / -64 /
!      data imach(13) /  63 /
!      data imach(14) /  14 /
!      data imach(15) / -64 /
!      data imach(16) /  63 /
!
!  IBM PC - Microsoft FORTRAN
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     6 /
!      data imach( 4) /     0 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /     2 /
!      data imach(11) /    24 /
!      data imach(12) /  -126 /
!      data imach(13) /   127 /
!      data imach(14) /    53 /
!      data imach(15) / -1022 /
!      data imach(16) /  1023 /
!
!  IBM PC - Professional FORTRAN and Lahey FORTRAN
!
!      data imach( 1) /     4 /
!      data imach( 2) /     7 /
!      data imach( 3) /     7 /
!      data imach( 4) /     0 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /     2 /
!      data imach(11) /    24 /
!      data imach(12) /  -126 /
!      data imach(13) /   127 /
!      data imach(14) /    53 /
!      data imach(15) / -1022 /
!      data imach(16) /  1023 /
!
!  INTERDATA 8/32 with the UNIX system FORTRAN 77 compiler.
!  For the INTERDATA FORTRAN VII compiler, replace the Z's specifying hex
!  constants with Y's.
!
!      data imach( 1) /   5 /
!      data imach( 2) /   6 /
!      data imach( 3) /   6 /
!      data imach( 4) /   6 /
!      data imach( 5) /  32 /
!      data imach( 6) /   4 /
!      data imach( 7) /   2 /
!      data imach( 8) /  31 /
!      data imach( 9) / Z'7FFFFFFF' /
!      data imach(10) /  16 /
!      data imach(11) /   6 /
!      data imach(12) / -64 /
!      data imach(13) /  62 /
!      data imach(14) /  14 /
!      data imach(15) / -64 /
!      data imach(16) /  62 /
!
!  PDP-10 (KA processor).
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   36 /
!      data imach( 6) /    5 /
!      data imach( 7) /    2 /
!      data imach( 8) /   35 /
!      data imach( 9) / "377777777777 /
!      data imach(10) /    2 /
!      data imach(11) /   27 /
!      data imach(12) / -128 /
!      data imach(13) /  127 /
!      data imach(14) /   54 /
!      data imach(15) / -101 /
!      data imach(16) /  127 /
!
!  PDP-10 (KI processor).
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   36 /
!      data imach( 6) /    5 /
!      data imach( 7) /    2 /
!      data imach( 8) /   35 /
!      data imach( 9) / "377777777777 /
!      data imach(10) /    2 /
!      data imach(11) /   27 /
!      data imach(12) / -128 /
!      data imach(13) /  127 /
!      data imach(14) /   62 /
!      data imach(15) / -128 /
!      data imach(16) /  127 /
!
!  PDP-11 FORTRANS supporting 32-bit integer ( kind = 4 ) arithmetic.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   32 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -127 /
!      data imach(13) /  127 /
!      data imach(14) /   56 /
!      data imach(15) / -127 /
!      data imach(16) /  127 /
!
!  PDP-11 FORTRANS supporting 16-bit integer ( kind = 4 ) arithmetic.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   16 /
!      data imach( 6) /    2 /
!      data imach( 7) /    2 /
!      data imach( 8) /   15 /
!      data imach( 9) / 32767 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -127 /
!      data imach(13) /  127 /
!      data imach(14) /   56 /
!      data imach(15) / -127 /
!      data imach(16) /  127 /
!
!  PRIME 50 series systems with 32-bit integers and 64V MODE instructions,
!  supplied by Igor Bray.
!
!      data imach( 1) /            1 /
!      data imach( 2) /            1 /
!      data imach( 3) /            2 /
!      data imach( 4) /            1 /
!      data imach( 5) /           32 /
!      data imach( 6) /            4 /
!      data imach( 7) /            2 /
!      data imach( 8) /           31 /
!      data imach( 9) / :17777777777 /
!      data imach(10) /            2 /
!      data imach(11) /           23 /
!      data imach(12) /         -127 /
!      data imach(13) /         +127 /
!      data imach(14) /           47 /
!      data imach(15) /       -32895 /
!      data imach(16) /       +32637 /
!
!  SEQUENT BALANCE 8000.
!
!      data imach( 1) /     0 /
!      data imach( 2) /     0 /
!      data imach( 3) /     7 /
!      data imach( 4) /     0 /
!      data imach( 5) /    32 /
!      data imach( 6) /     1 /
!      data imach( 7) /     2 /
!      data imach( 8) /    31 /
!      data imach( 9) /  2147483647 /
!      data imach(10) /     2 /
!      data imach(11) /    24 /
!      data imach(12) /  -125 /
!      data imach(13) /   128 /
!      data imach(14) /    53 /
!      data imach(15) / -1021 /
!      data imach(16) /  1024 /
!
!  SUN Microsystems UNIX F77 compiler.
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     6 /
!      data imach( 4) /     0 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    32 /
!      data imach( 9) /2147483647/
!      data imach(10) /     2 /
!      data imach(11) /    24 /
!      data imach(12) /  -126 /
!      data imach(13) /   128 /
!      data imach(14) /    53 /
!      data imach(15) / -1022 /
!      data imach(16) /  1024 /
!
!  SUN 3 (68881 or FPA)
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    6 /
!      data imach( 4) /    0 /
!      data imach( 5) /   32 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -125 /
!      data imach(13) /  128 /
!      data imach(14) /   53 /
!      data imach(15) / -1021 /
!      data imach(16) /  1024 /
!
!  UNIVAC 1100 series.
!  Note that the punch unit, I1MACH(3), has been set to 7, which is appropriate
!  for the UNIVAC-FOR system.  If you have the UNIVAC-FTN system, set it to 1
!  instead.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   36 /
!      data imach( 6) /    6 /
!      data imach( 7) /    2 /
!      data imach( 8) /   35 /
!      data imach( 9) / O377777777777 /
!      data imach(10) /    2 /
!      data imach(11) /   27 /
!      data imach(12) / -128 /
!      data imach(13) /  127 /
!      data imach(14) /   60 /
!      data imach(15) /-1024 /
!      data imach(16) / 1023 /
!
!  VAX.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   32 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -127 /
!      data imach(13) /  127 /
!      data imach(14) /   56 /
!      data imach(15) / -127 /
!      data imach(16) /  127 /
!
!  Z80 microprocessor.
!
!      data imach( 1) /    1 /
!      data imach( 2) /    1 /
!      data imach( 3) /    0 /
!      data imach( 4) /    1 /
!      data imach( 5) /   16 /
!      data imach( 6) /    2 /
!      data imach( 7) /    2 /
!      data imach( 8) /   15 /
!      data imach( 9) / 32767 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -127 /
!      data imach(13) /  127 /
!      data imach(14) /   56 /
!      data imach(15) / -127 /
!      data imach(16) /  127 /
!
  if ( i < 1 .or. 16 < i )then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I1MACH - Fatal error!'
    write ( *, '(a,i8)' ) '  I is out of bounds:', i
    i1mach = 0
    stop
  else
    i1mach = imach(i)
  end if

  return
end
subroutine xerabt ( messg, nmessg )

!*****************************************************************************80
!
!! XERABT aborts program execution and prints an error message.
!
!  Discussion:
!
!    XERABT aborts the execution of the program.  The error message causing
!    the abort is given in the calling sequence, in case one needs it for
!    printing on a dayfile, for example.
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, character ( len = * ) MESSG, the message to be processed,
!    containing no more than 72 characters.
!
!    Input, integer ( kind = 4 ) NMESSG, the actual number of characters 
!    in MESSG.  If NMESSG is 0, no message is being supplied.
!
  implicit none

  character ( len = * ) messg
  integer ( kind = 4 ) nmessg

  if ( 0 < nmessg ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XERABT - Termination after fatal error!'
    write ( *, '(a)' ) trim ( messg )
  end if

  stop
end
subroutine xerclr

!*****************************************************************************80
!
!! XERCLR resets the current error number to zero.
!
!  Discussion:
!
!    This routine simply resets the current error number to zero.
!    This may be necessary to do in order to determine that
!    a certain error has occurred again since the last time
!    NUMXER was referenced.
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    None
!
  implicit none

  integer ( kind = 4 ) j4save
  integer ( kind = 4 ) junk

  junk = j4save ( 1, 0, .true. )

  return
end
subroutine xerdmp

!*****************************************************************************80
!
!! XERDMP prints the error tables and then clears them.
!
!  Modified:
!
!    05 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    None
!
  implicit none

  integer ( kind = 4 ) kount

  call xersav ( ' ', 0, 0, 0, kount )

  return
end
subroutine chfdv ( x1, x2, f1, f2, d1, d2, ne, xe, fe, de, next, ierr )

!*****************************************************************************80
!
!! CHFDV evaluates a cubic polynomial and its derivative given in Hermite form.
!
!  Discussion:
!
!    CHFDV evaluates a cubic polynomial and its first derivative.
!    The cubic polynomial is given in Hermite form.  The evaluation
!    is carried out at an array of points.
!
!    This routine was designed for use by PCHFD, but it may also be
!    useful directly as an evaluator for a piecewise cubic Hermite
!    function in applications, such as graphing, where the interval
!    is known in advance.
!
!    If only function values are required, use CHFEV instead.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( dp ) X1, X2, the endpoints of the interval of
!    definition of  the cubic.  X1 and X2 must be distinct.
!
!    Input, real ( dp ) F1, F2, the values of the function at X1 and
!    X2, respectively.
!
!    Input, real ( dp ) D1, D2, the derivative values at the ends
!     of the interval.
!
!    Input, integer ( kind = 4 ) NE, the number of evaluation points.
!
!    Input, real ( dp ) XE(NE), the points at which the functions are to
!    be evaluated.  If any of the XE are outside the interval
!    [X1,X2], a warning error is returned in next.
!
!    Output, real ( dp ) FE(NE), DE(NE), the values of the cubic
!    function and its derivative at the points XE(*).
!
!    Output, integer ( kind = 4 ) NEXT(2), indicates the number of 
!    extrapolation points:
!    NEXT(1) = number of evaluation points to left of interval.
!    NEXT(2) = number of evaluation points to right of interval.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    -1, NE < 1.
!    -2, X1 == X2.
!
  implicit none

  integer, parameter  :: dp=kind(0.d0)                 

  integer ( kind = 4 ) ne

  real ( kind=dp ) c2
  real ( kind=dp ) c2t2
  real ( kind=dp ) c3
  real ( kind=dp ) c3t3
  real ( kind=dp ) d1
  real ( kind=dp ) d2
  real ( kind=dp ) de(ne)
  real ( kind=dp ) del1
  real ( kind=dp ) del2
  real ( kind=dp ) delta
  real ( kind=dp ) f1
  real ( kind=dp ) f2
  real ( kind=dp ) fe(ne)
  real ( kind=dp ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) next(2)
  real ( kind=dp ) x
  real ( kind=dp ) x1
  real ( kind=dp ) x2
  real ( kind=dp ) xe(ne)
  real ( kind=dp ) xma
  real ( kind=dp ) xmi
!
!  Check arguments.
!
  if ( ne < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHFDV - Fatal error!'
    write ( *, '(a)' ) '  The number of evaluation points was less than 1.'
    stop
  end if

  h = x2 - x1

  if ( h == 0.0_dp ) then
    ierr = -2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHFDV - Fatal error!'
    write ( *, '(a)' ) '  The interval endpoints are equal.'
    return
  end if
!
!  Initialize.
!
  ierr = 0
  next(1) = 0
  next(2) = 0
  xmi = min ( 0.0_dp, h )
  xma = max ( 0.0_dp, h )
!
!  Compute cubic coefficients expanded about X1.
!
  delta = ( f2 - f1 ) / h
  del1 = ( d1 - delta ) / h
  del2 = ( d2 - delta ) / h

  c2 = -( del1 + del1 + del2 )
  c2t2 = c2 + c2
  c3 = ( del1 + del2 ) / h
  c3t3 = c3 + c3 + c3
!
!  Evaluation loop.
!
  do i = 1, ne

    x = xe(i) - x1
    fe(i) = f1 + x * ( d1 + x * ( c2 + x * c3 ) )
    de(i) = d1 + x * ( c2t2 + x * c3t3 )
!
!  Count extrapolation points.
!
    if ( x < xmi ) then
      next(1) = next(1) + 1
    end if

    if ( xma < x ) then
      next(2) = next(2) + 1
    end if

  end do

  return
end
function pchdf ( k, x, s, ierr )

!*****************************************************************************80
!
!! PCHDF approximates a derivative using divided differences of data.
!
!  Discussion:
!
!    The routine uses a divided difference formulation to compute a K-point
!    approximation to the derivative at X(K) based on the data in X and S.
!
!    It is called by PCHCE and PCHSP to compute 3 and 4 point boundary
!    derivative approximations.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Pages 10-16,
!    Springer-Verlag, 1978.
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, is the order of the desired derivative 
!    approximation.  K must be at least 3.
!
!    Input, real ( dp ) X(K), contains the K values of the independent
!    variable.  X need not be ordered, but the values must be distinct.
!
!    Input/output, real ( dp ) S(K-1).  On input, the associated slope
!    values:
!      S(I) = ( F(I+1)-F(I))/(X(I+1)-X(I))
!    On output, S is overwritten.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no error.
!    -1, if K < 2.
!
!    Output, real ( dp ) PCHDF, the desired derivative approximation if
!    IERR=0 or to zero if IERR=-1.
!
  implicit none

  integer, parameter  :: dp=kind(0.d0)                 

  integer ( kind = 4 ) k

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  real ( kind=dp ) pchdf
  real ( kind=dp ) s(k-1)
  real ( kind=dp ) value
  real ( kind=dp ) x(k)
!
!  Check for legal value of K.
!
  if ( k < 3 ) then
    ierr = -1
    call xerror ( 'pchdf -- k less than three', ierr, 1 )
    pchdf = 0.0_dp
    return
  end if
!
!  Compute coefficients of interpolating polynomial.
!
  do j = 2, k-1
    do i = 1, k-j
      s(i) = ( s(i+1) - s(i) ) / ( x(i+j) - x(i) )
    end do
  end do
!
!  Evaluate the derivative at X(K).
!
  value = s(1)

  do i = 2, k-1
    value = s(i) + value * ( x(k) - x(i) )
  end do

  ierr = 0
  pchdf = value

  return
end
function pchst ( arg1, arg2 )

!*****************************************************************************80
!
!! PCHST: PCHIP sign-testing routine.
!
!  Discussion:
!
!    This routine essentially computes the sign of ARG1 * ARG2.
!
!    The object is to do this without multiplying ARG1 * ARG2, to avoid
!    possible over/underflow problems.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
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
!    Input, real ( dp ) ARG1, ARG2, two values to check.
!
!    Output, real ( dp ) PCHST,
!    -1.0, if ARG1 and ARG2 are of opposite sign.
!     0.0, if either argument is zero.
!    +1.0, if ARG1 and ARG2 are of the same sign.
!
  implicit none

  integer, parameter  :: dp=kind(0.d0)                 

  real ( kind=dp ) arg1
  real ( kind=dp ) arg2
  real ( kind=dp ) pchst

  pchst = sign ( 1.0_dp, arg1 ) * sign ( 1.0_dp, arg2 )

  if ( arg1 == 0.0_dp .or. arg2 == 0.0_dp ) then
    pchst = 0.0_dp
  end if

  return
end

