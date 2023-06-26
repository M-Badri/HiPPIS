#include "fintrf.h"
!!
!!     Gateway routine for adaptiveInterpolation1D(...)
!!
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      use mod_adaptiveInterpolation

!!     Declarations
      implicit none

!!     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs

!!     Function declarations:

#if MX_HAS_INTERLEAVED_COMPLEX
      mwPointer mxGetDoubles
#else
      mwPointer mxGetPr
#endif

      mwPointer mxCreateDoubleMatrix
      integer mxIsNumeric
      mwPointer mxGetM, mxGetN

!!    Variable for function calls
      real(dp), dimension(:), pointer :: xin(:), vin(:)
      real(dp), dimension(:), pointer :: xout(:), vout(:)
      integer, dimension(:), pointer :: deg(:)
!!    Pointers to input/output mxArrays:
      mwPointer xin_ptr, xout_ptr, vin_ptr, vout_ptr, deg_ptr
      mwPointer degree_ptr, interpolation_ptr
      mwPointer sten_ptr, eps0_ptr, eps1_ptr

!!     Array information:
      mwPointer m, n

!!     Arguments for computational routine:
      integer  interpolation_type, d, sten
      real(dp) degree, interpolation, stencil, eps0, eps1

!!-----------------------------------------------------------------------
!!     Check for proper number of arguments. 
      if(nrhs < 4 .or. nrhs > 8) then
         call mexErrMsgIdAndTxt ('MATLAB:timestwo:nInput', &
                                'One input required.')
      elseif(nlhs < 1 .or. nlhs > 2) then
         call mexErrMsgIdAndTxt ('MATLAB:timestwo:nOutput', &
                                'Too many output arguments.')
      endif

!!     Validate inputs
!!     Check that the input is a number.
      if(mxIsNumeric(prhs(1)) .eq. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:timestwo:NonNumeric', &
                                'Input must be a number.')
      endif


!!     Get the size of the input array.
      m = mxGetN(prhs(3))
      n = mxGetN(prhs(1))

!!    Allocate space for arrays
      allocate(xin(n))      
      allocate(vin(n))      
      allocate(xout(m))      
      allocate(vout(m))      
      if(nlhs == 2) then
        allocate(deg(n-1))      
      endif

!!     Create Fortran array from the input argument.

#if MX_HAS_INTERLEAVED_COMPLEX
      xin_ptr = mxGetDoubles(prhs(1))
      vin_ptr = mxGetDoubles(prhs(2))
      xout_ptr = mxGetDoubles(prhs(3))
      degree_ptr = mxGetDoubles(prhs(4)
      interpolation_ptr = mxGetDoubles(prhs(5))
#else
      xin_ptr = mxGetPr(prhs(1))
      vin_ptr = mxGetPr(prhs(2))
      xout_ptr = mxGetPr(prhs(3))
      degree_ptr = mxGetPr(prhs(4))
      interpolation_ptr = mxGetPr(prhs(5))
      if(nrhs > 5) then
        sten_ptr = mxGetPr(prhs(6))
      endif
      if(nrhs > 6) then
        eps0_ptr = mxGetPr(prhs(7))
      endif
      if(nrhs > 7) then
        eps1_ptr = mxGetPr(prhs(8))
      endif
#endif
      !!** Obtain the input information **!!
      call mxCopyPtrToReal8(degree_ptr, degree, 1)
      call mxCopyPtrToReal8(interpolation_ptr, interpolation, 1)
      d = int(degree)
      interpolation_type = int(interpolation)
      call mxCopyPtrToReal8(xin_ptr, xin, n)
      call mxCopyPtrToReal8(vin_ptr, vin, n)
      call mxCopyPtrToReal8(xout_ptr, xout, m)
      if(nrhs > 5) then
        call mxCopyPtrToReal8(sten_ptr, stencil, 1)
        sten = int(stencil)
      endif
      if(nrhs > 6) then
        call mxCopyPtrToReal8(eps0_ptr, eps0, 1)
      endif
      if(nrhs > 7) then
        call mxCopyPtrToReal8(eps1_ptr, eps1, 1)
      endif


!!     Create matrix for the return argument.
      plhs(1) = mxCreateDoubleMatrix(1,m,0)
      if(nlhs ==2) then
        plhs(2) = mxCreateDoubleMatrix(1,n-1,0)
      endif
#if MX_HAS_INTERLEAVED_COMPLEX
      vout_ptr = mxGetDoubles(plhs(1))
      if(nlhs ==2) then
        deg_ptr = mxGetDoubles(plhs(2))
      endif
#else
      vout_ptr = mxGetPr(plhs(1))
      if(nlhs ==2) then
        deg_ptr = mxGetPr(plhs(2))
      endif
#endif

!!     Call the computational subroutine.
      if(nrhs == 5) then
       call adaptiveinterpolation1D(xin, vin, n, &
          xout, vout, m, d, interpolation_type )
      elseif(nrhs == 6) then
       call adaptiveinterpolation1D(xin, vin, n, &
          xout, vout, m, d, interpolation_type, sten )
      elseif(nrhs == 7) then
       call adaptiveinterpolation1D(xin, vin, n, &
          xout, vout, m, d, interpolation_type, sten, eps0)
      elseif(nrhs == 8 .and. nlhs == 1) then
       call adaptiveinterpolation1D(xin, vin, n, &
          xout, vout, m, d, interpolation_type, sten, eps0, eps1)
      elseif(nrhs == 8 .and. nlhs == 2) then
       call adaptiveinterpolation1D(xin, vin, n, &
          xout, vout, m, d, interpolation_type, sten, eps0, eps1, deg)
      endif


!!!   !!** Load the data into y_ptr, which is the output to MATLAB.
      call mxCopyReal8ToPtr(vout,vout_ptr,m)     
      if(nlhs ==2) then
        call mxCopyInteger1ToPtr(real(deg, kind=dp),deg_ptr,n-1)     
      endif
!
      return
      end


