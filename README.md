# HiPPIS
HiPPIS is a polynomial-based data-bounded and positivity-preserving interpolation software for function approximation and mapping data values between structure (uniform and nonuniform) meshes.


## Getting Started

### Dependencies

* Matlab is required for the version of the software implemented in Matlab. No special toolboxes are required.
  The core files used for the implementation of the data-bounded and positivity-preserving interpolation methods are in the folder *Src*.
  The folder structure for the Malab version is as follows:
```
  Matlab
  |  Src
  |  |  adaptiveInterpolation1D.m
  |  |  adaptiveInterpolation2D.m
  |  |  adaptiveInterpolation3D.m
  |  |  divdiff.m
  |  |  newtonPolyVal.m
  |  Drivers
  |  |  Tutorial
  |  |  |  main.m
  |  |  |  testing.m 
  |  |  Mapping 
  |  |  |  main.m
  |  Extras
  |  |  ...
```
  The remaining folders and files in the folder *Matlab* are drivers, examples, data, and scripts for using the data-bounded and positivity-preserving interpolation methods.
* The Fortran version requires the installation of Intel (ifort), or gnu (gfortran) compilers with OpenMP4.
  The vectorized version of the code required the use of Intel compilers. 
  The Make files *Fortran/Drivers/Mapping/Makefile* and *Fortran/Drivers/BOMEX/Makefile* can be modified for different compilers including Intel, gnu, HPE cray.
  The core files used for the implementation and demo of the data-bounded and positivity-preserving interpolation methods are in the folder *Src*.
  The folder structure for the Matlab version is as follows:
```
  Fortran
  |  Src
  |  |  mod_adaptiveInterpolation.F90
  |  |  MExFiles
  |  |  |  apdaptiveInterpolation1D.F90
  |  |  |  apdaptiveInterpolation2D.F90
  |  |  |  apdaptiveInterpolation3D.F90
  |  Drivers
  |  |  Tutorial
  |  |  |  main.F90
  |  |  |  testing.F90
  |  |  |  main.m
  |  |  Mapping
  |  |  |  main.F90
  |  |  BOMEX
  |  |  |  main.F90
  |  Extras
  |  |  ...
```
  The remaining folders and files in the folder *Fortran* are drivers, examples, data, and scripts for using the data-bounded and positivity-preserving interpolation methods.

### Installing
* Downloading HiPPIS (Matlab and/or Fortran) is sufficient for installation. 

### Executing program
* Matlab

Open Matlab and run the script *tutorial.m* in folder *Matlab*.
```
 cd Matlab/Drivers/Tutorial/
 main
```


* Fortran
```
cd Fortran/Drivers/Tutorial/
make 
./main
``` 
### Help

The file *main.m* (in *Matlab/Drivers/Tutorial*) or *main.m* or *main.F90* (in *Fortran/Drivers/Tutorial*) is a good place to start. 
Both *main.m* or *main.F90* provide simple examples showing how to use HiPPIS with 1D, 2D, and 3D structured tensor product meshes.

## Manuscript Examples
More examples can be found in */Matlab/Drivers/Mappin/main.m* or */Fortran/Drivers/Mapping/main.F90*.
These examples are used to produce the results presented in a manuscript submitted for publication that is entitled "Algorithm xxxx: HiPPIS: A High-Order Positivity-Preserving Mapping Software for Mtructured Meshes". Each example produces results that are saved in files. The saved data are then used to calculate the errors and produce the figures and tables in the manuscript. 

**Note:** Producing the results in the Manuscript takes a few hours because 1) the PDE problem in folder *BOMEX* is executed multiple times and uses a small time step, 2) the 2D approximation examples are executed multiple times and use 1000x1000 output mesh.

### Matlab Examples
```
cd Matlab/Mapping
main
```
The driver for the examples is */Matlab/Drivers/Mapping/main.m*.
### Examples in Fortran
```
cd Fortran/Extras/
sh run_manuscript_examples.sh
```
Open the Matlab software and run 
```
plot_manuscript_examples
```
The drivers for the examples are *Fortran/Drivers/Mapping/main.F90* and *Fortran/Drivers/BOMEX/main.F90*.

## Testing
Supplemental tests are provided in the files *Mapping/Drivers/Tutorial/testing.m* and *Fortran/Drivers/Tutorial/testing.F90*
