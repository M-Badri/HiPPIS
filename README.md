# HiPPIS
HiPPIS is a polynomial-based data-bounded and positivity-preserving interpolation software for function approximation and mapping data values between meshes.


## Getting Started

### Dependencies

* Matlab is required for the version of the software implemented in Matlab. No special toolboxes are required.
* The Frotran version requires the installation of Intel (icc, ifort), gnu (gcc, gfortran) compilers with OpenMP4.
  The vectorized version of the code required the use of Intel compilers. 
  The Make files can be modified for different compilers including Intel, gnu, HPE cray.

### Installing
* Downloading the desired version (Matlab and/or Fortran) is sufficient for installation. 

### Executing program
* Fortran
```
make 
./a.out
``` 
* Matlab
```
 main
```
or 
```
tutorial
```

### Help

The file *tutorial.m* is a good place to start. 
*tutorial.m* provides examples showing how HiPPIS for 1D, 2D, and 3D  with structured tensor product meshes.

### Examples
More examples can be found in *main.m*.
These examples produce the results presented in a manuscript submitted for publication that is entitled "HiPPIS:A High-Order Positivity-Preserving Mapping Software for Structured Meshes". 

