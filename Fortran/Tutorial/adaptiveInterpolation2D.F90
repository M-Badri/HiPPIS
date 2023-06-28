#include "fintrf.h"
!!
!!     Gateway routine for adaptiveInterpolation2D(...)
!!
!!      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      use mod_adaptiveInterpolation
      use omp_lib

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
!!-      integer i, j
      mwPointer mxGetM, mxGetN

!!     Pointers to input/output mxArrays:
      real(dp), dimension(:), allocatable :: xin(:), yin(:), vin(:,:)
      real(dp), dimension(:), allocatable :: xout(:), yout(:), vout(:,:)
      real(dp):: time
!!-      integer, dimension(:),  allocatable :: deg(:)
      !!real*8 xin_ptr(:), yin_ptr(:), xout_ptr(:), yout_ptr(:), degree_ptr(:)
      !!integer xin_ptr;
      mwPointer xin_ptr, xout_ptr, yin_ptr, yout_ptr, deg_ptr
      mwPointer vin_ptr, vout_ptr
      mwPointer degree_ptr, interpolation_ptr
      mwPointer sten_ptr, eps0_ptr, eps1_ptr

!!     Array information:
      mwPointer mx, my, nx, ny
!-      mwSize size

!!     Arguments for computational routine:
      integer  interpolation_type, d, sten
      real(dp) degree, interpolation, stencil, eps0, eps1
      !real(dp)  xin, yin, xout(m), yout(m)

      time = omp_get_wtime()
!!-----------------------------------------------------------------------
!!     Check for proper number of arguments. 
      if(nrhs < 7 .or. nrhs > 10) then
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
      mx = mxGetN(prhs(4))
      my = mxGetN(prhs(5))

!!    Allocate space for arrays
      allocate(xin(nx))      
      allocate(yin(ny))      
      allocate(xout(mx))      
      allocate(yout(my))      
      allocate(vin(nx,ny))      
      allocate(vout(mx,my))      
      !!if(nlhs == 2) then
      !!  allocate(deg(n-1))      
      !!endif

!!     Create Fortran array from the input argument.

#if MX_HAS_INTERLEAVED_COMPLEX
      xin_ptr = mxGetDoubles(prhs(1))
      yin_ptr = mxGetDoubles(prhs(2))
      vin_ptr = mxGetDoubles(prhs(3))
      xout_ptr = mxGetDoubles(prhs(4))
      yout_ptr = mxGetDoubles(prhs(5))
      degree_ptr = mxGetDoubles(prhs(6)
      interpolation_ptr = mxGetDoubles(prhs(7))
#else
      xin_ptr = mxGetPr(prhs(1))
      yin_ptr = mxGetPr(prhs(2))
      vin_ptr = mxGetPr(prhs(3))
      xout_ptr = mxGetPr(prhs(4))
      yout_ptr = mxGetPr(prhs(5))
      degree_ptr = mxGetPr(prhs(6))
      interpolation_ptr = mxGetPr(prhs(7))
      if(nrhs > 7) then
        sten_ptr = mxGetPr(prhs(8))
      endif
      if(nrhs > 8) then
        eps0_ptr = mxGetPr(prhs(9))
      endif
      if(nrhs > 9) then
        eps1_ptr = mxGetPr(prhs(10))
      endif
#endif
      !!** Obtain the input information **!!
      !size = mxGetN(prhs(4))*mxGetM(prhs(4))
      call mxCopyPtrToReal8(degree_ptr, degree, 1)
      call mxCopyPtrToReal8(interpolation_ptr, interpolation, 1)
      d = int(degree)
      interpolation_type = int(interpolation)
      call mxCopyPtrToReal8(xin_ptr, xin, nx)
      call mxCopyPtrToReal8(yin_ptr, yin, ny)
      call mxCopyPtrToReal8(vin_ptr, vin, ny*nx)
      call mxCopyPtrToReal8(xout_ptr, xout, mx)
      call mxCopyPtrToReal8(yout_ptr, yout, my)
      !!call mxCopyPtrToReal8(vout_ptr, vout, my*mx)
      if(nrhs > 7) then
        call mxCopyPtrToReal8(sten_ptr, stencil, 1)
        sten = int(stencil)
      endif
      if(nrhs > 8) then
        call mxCopyPtrToReal8(eps0_ptr, eps0, 1)
      endif
      if(nrhs > 9) then
        call mxCopyPtrToReal8(eps1_ptr, eps1, 1)
      endif


      !!** Create matrix for the return argument.
!!-      print *, mx, my
      plhs(1) = mxCreateDoubleMatrix(mx,my,0)

      !!if(nlhs ==2) then
      !!  !plhs(2) = mxCreateNumericMatrix(1,n-1,0)
      !!  plhs(2) = mxCreateDoubleMatrix(1,n-1,0)
      !!endif
#if MX_HAS_INTERLEAVED_COMPLEX
      vout_ptr = mxGetDoubles(plhs(1))
      !if(nlhs ==2) then
      !  deg_ptr = mxGetDoubles(plhs(2))
      !endif
#else
      vout_ptr = mxGetPr(plhs(1))
      !if(nlhs ==2) then
      !  deg_ptr = mxGetPr(plhs(2))
      !endif
#endif

!!     Call the computational subroutine.
      if(nrhs == 7 .and. nlhs == 1) then
        call adaptiveinterpolation2D(xin, yin, nx, ny, vin, &
          xout, yout, mx, my, vout, d, interpolation_type)
      elseif(nrhs == 8 .and. nlhs ==1) then
        call adaptiveinterpolation2D(xin, yin, nx, ny, vin, &
          xout, yout, mx, my, vout, d, interpolation_type, & 
          sten)
      elseif(nrhs == 9 .and. nlhs == 1) then
        call adaptiveinterpolation2D(xin, yin, nx, ny, vin, &
          xout, yout, mx, my, vout, d, interpolation_type, & 
          sten, eps0)
      elseif(nrhs == 10 .and. nlhs == 1) then
        call adaptiveinterpolation2D(xin, yin, nx, ny, vin, &
          xout, yout, mx, my, vout, d, interpolation_type, & 
          sten, eps0, eps1)
      else
        write(*,*) 'ERROR: Incorrect number inputs.'
        stop
      endif

      !do i=1,mx
      !do j=1,my
      !vout(i,j) = 0.0_dp;
      !enddo
      !enddo


!!!   !!** Load the data into y_ptr, which is the output to MATLAB.
      call mxCopyReal8ToPtr(vout,vout_ptr, mx*my)     
      !if(nlhs ==2) then
      !  call mxCopyInteger1ToPtr(real(deg, kind=dp),deg_ptr,n-1)     
      !endif
!
!!    Allocate space for arrays
      !deallocate(xin)      
      !deallocate(yin)      
      !deallocate(xout)      
      !deallocate(yout)      
      !deallocate(vin)      
      write(*,*) 'Total 2D time ', omp_get_wtime()-time
 
      return
      end


