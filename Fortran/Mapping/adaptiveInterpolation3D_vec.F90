#include "fintrf.h"
!!
!!     Gateway routine for adaptiveInterpolation2D(...)
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
      integer i, j, k
      mwPointer mxGetM, mxGetN

!!     Pointers to input/output mxArrays:
      real(dp), dimension(:), allocatable :: xin(:), yin(:), zin(:)
      real(dp), dimension(:), allocatable :: xout(:), yout(:), zout(:)
      real(dp), dimension(:), allocatable :: vin(:,:,:), vout(:,:,:)

      mwPointer xin_ptr, xout_ptr, yin_ptr, yout_ptr, deg_ptr
      mwPointer vin_ptr, vout_ptr, zin_ptr, zout_ptr
      mwPointer degree_ptr, interpolation_ptr
      mwPointer sten_ptr, eps0_ptr, eps1_ptr

      mwPointer mxCreateNumericArray
      integer*4 mxClassIDFromClassName

      !! Arguments for mxCreateNumericArray
      integer*4 classid
      integer*4 complexflag
      mwSize ndim
      mwSize dims(3)

!!     Array information:
      mwPointer mx, my, mz, nx, ny, nz
!-      mwSize size

!!     Arguments for computational routine:
      integer  interpolation_type, d, sten
      real(dp) degree, interpolation, stencil, eps0, eps1

!!-----------------------------------------------------------------------
!!     Check for proper number of arguments. 
      if(nrhs < 7 .or. nrhs > 12) then
         call mexErrMsgIdAndTxt ('MATLAB:timestwo:nInput', &
                                'at least 7 and at most 10', &
                                'input arguments are required.')
      elseif(nlhs .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:adaptiveInterpolation2D:', &
                                 'nOutput must have only one output',&
                                 'arguments.')
      endif

!!     Validate inputs
!!     Check that the input is a number.
      if(mxIsNumeric(prhs(1)) .eq. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:timestwo:NonNumeric', &
                                'Input must be a number.')
      endif


!!     Get the size of the input array.
      nx = mxGetN(prhs(1))
      ny = mxGetN(prhs(2))
      nz = mxGetN(prhs(3))
      mx = mxGetN(prhs(5))
      my = mxGetN(prhs(6))
      mz = mxGetN(prhs(7))

!!    Allocate space for arrays
      allocate(xin(nx))      
      allocate(yin(ny))      
      allocate(zin(nz))      
      allocate(xout(mx))      
      allocate(yout(my))      
      allocate(zout(mz))      
      allocate(vin(nx,ny,nz))      
      allocate(vout(mx,my,mz))      

!!     Create Fortran array from the input argument.

#if MX_HAS_INTERLEAVED_COMPLEX
      xin_ptr = mxGetDoubles(prhs(1))
      yin_ptr = mxGetDoubles(prhs(2))
      zin_ptr = mxGetDoubles(prhs(3))
      vin_ptr = mxGetDoubles(prhs(4))
      xout_ptr = mxGetDoubles(prhs(5))
      yout_ptr = mxGetDoubles(prhs(6))
      zout_ptr = mxGetDoubles(prhs(7))
      degree_ptr = mxGetDoubles(prhs(8)
      interpolation_ptr = mxGetDoubles(prhs(9))
#else
      xin_ptr = mxGetPr(prhs(1))
      yin_ptr = mxGetPr(prhs(2))
      zin_ptr = mxGetPr(prhs(3))
      vin_ptr = mxGetPr(prhs(4))
      xout_ptr = mxGetPr(prhs(5))
      yout_ptr = mxGetPr(prhs(6))
      zout_ptr = mxGetPr(prhs(7))
      degree_ptr = mxGetPr(prhs(8))
      interpolation_ptr = mxGetPr(prhs(9))
      if(nrhs > 9) then
        sten_ptr = mxGetPr(prhs(10))
      endif
      if(nrhs > 10) then
        eps0_ptr = mxGetPr(prhs(11))
      endif
      if(nrhs > 11) then
        eps1_ptr = mxGetPr(prhs(12))
      endif
#endif
      !!** Obtain the input information **!!
      call mxCopyPtrToReal8(degree_ptr, degree, 1)
      call mxCopyPtrToReal8(interpolation_ptr, interpolation, 1)
      d = int(degree)
      interpolation_type = int(interpolation)
      call mxCopyPtrToReal8(xin_ptr, xin, nx)
      call mxCopyPtrToReal8(yin_ptr, yin, ny)
      call mxCopyPtrToReal8(zin_ptr, zin, nz)
      call mxCopyPtrToReal8(vin_ptr, vin, ny*nx*nz)
      call mxCopyPtrToReal8(xout_ptr, xout, mx)
      call mxCopyPtrToReal8(yout_ptr, yout, my)
      call mxCopyPtrToReal8(zout_ptr, zout, mz)
      if(nrhs > 9) then
        call mxCopyPtrToReal8(sten_ptr, stencil, 1)
        sten = int(stencil)
      endif
      if(nrhs > 10) then
        call mxCopyPtrToReal8(eps0_ptr, eps0, 1)
      endif
      if(nrhs > 11) then
        call mxCopyPtrToReal8(eps1_ptr, eps1, 1)
      endif


      !!** Create matrix for the return argument.
      print *, nx, ny, nz
      print *, mx, my, mz
      classid = mxClassIDFromclassName('double')
      complexflag = 0
      ndim = 3
      dims(1) = mx
      dims(2) = my
      dims(3) = mz
      plhs(1) = mxCreateNumericArray(ndim, dims, classid,&
                complexflag)

#if MX_HAS_INTERLEAVED_COMPLEX
      vout_ptr = mxGetDoubles(plhs(1))
#else
      vout_ptr = mxGetPr(plhs(1))
#endif

!!     Call the computational subroutine.
      if(nrhs == 9 .and. nlhs == 1) then
        call adaptiveinterpolation3D_vec(xin,yin,zin,nx,ny,nz,vin, &
          xout, yout, zout, mx, my, mz, vout, d, interpolation_type)
      elseif(nrhs == 10 .and. nlhs ==1) then
        call adaptiveinterpolation3D_vec(xin,yin,zin,nx,ny,nz,vin, &
          xout, yout, zout, mx, my, mz, vout, d, interpolation_type, & 
          sten)
      elseif(nrhs == 11 .and. nlhs == 1) then
        call adaptiveinterpolation3D_vec(xin,yin,zin,nx,ny,nz,vin,&
          xout, yout, zout, mx, my, mz, vout, d, interpolation_type, & 
          sten, eps0)
      elseif(nrhs == 12 .and. nlhs == 1) then
        call adaptiveinterpolation3D_vec(xin,yin,zin,nx,ny,nz,vin, &
          xout, yout, zout, mx, my, mz, vout, d, interpolation_type, & 
          sten, eps0, eps1)
      endif


!!!   !!** Load the data into y_ptr, which is the output to MATLAB.
      call mxCopyReal8ToPtr(vout,vout_ptr, mx*my*mz)     


      return
      end


