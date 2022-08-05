#!/bin/bash
echo ''
echo 'Running the examples in the manuscript [REF when published]. This will take'
echo 'a long time because the runs include 2D exmaples using 1000x1000 output meshes'
echo 'a time dependent PDE problem with small time steps'
echo ''
ulimit -s unlimited
echo 'Increased stack space for 2D examples using "ulimit -s unlimited" '
echo ''
echo 'The Makefiles in "HiPPIS/Fortran/Mapping" and "HiPPIS/Fortran/BOMEX"'
echo ' are set up for Intel compiler. For other compilers modify '
echo '"HiPPIS/Fortran/Mapping" and "HiPPIS/Fortran/BOMEX" for the desired compiler'
