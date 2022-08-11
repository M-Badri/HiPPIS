#!/bin/bash
echo ''
echo 'Running the examples in the manuscript [REF when published]. This will take'
echo 'a long time because the runs include 2D examples using 1000x1000 output meshes'
echo 'a time dependent PDE problem with small time steps'
echo ''
ulimit -s unlimited
echo 'Increased stack space for 2D examples using "ulimit -s unlimited" '
echo ''
echo 'The Makefiles in "HiPPIS/Fortran/Mapping" and "HiPPIS/Fortran/BOMEX"'
echo ' are set up for Intel compiler. For other compilers modify '
echo '"HiPPIS/Fortran/Mapping" and "HiPPIS/Fortran/BOMEX" for the desired compiler'
echo 'Each Makefile has a template for gfortran compiler than can be used.'
echo ''
echo 'Moving into the "Mapping" folder ... ' 
cd Mapping
echo ''
echo 'Building and running executable ... ' 
make clean
make
echo ''
echo 'Leaving Mapping folder and moving into the BOMEX folder ... '
cd ../BOMEX/
echo ''
echo 'Building and running executable ... ' 
make clean
make
echo ''
echo 'The vectorization results are in "vectorization_results"'
echo ''
echo 'Start Matlab and run the script "plot_manuscript_examples" to produce the '
echo 'the figures and the tables in the manuscript. '

