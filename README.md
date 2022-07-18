# HPPIS
HPPIS is a polynomail-based data-bounded and positivity-preserving interpolation software for function approximation and mapping data values between meshes.


## Getting Started

### Dependencies

* Matlab is required for the version of the software implemeted in Matlab. No special toolboxes are required.
* The Frotran version requires the installation of gnu compilers (gcc, g++, gfortran) with OpenMP for the vectorization. 
  The Make files can be modified for different compiler such intel and cray.

### Installing
* Dowloading the desired version (Matlab and/or Fortran) is sufficient for installation. 

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
*tutorial.m* provides examples showing how HPPIS for 1D, 2D, and 3D  with structures tensor product meshes.

### Examples
More examples can be found in *main.m*.
These examples are part of manuscript in preparation that will be submitted for publication. 

