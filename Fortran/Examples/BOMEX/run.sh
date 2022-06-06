#!/bin/bash

echo '... creating folder data to hold simulation data '
mkdir data
#
echo '... build excutable  a.out '
make
#
echo '... move executable in the folder data/'
mv a.out data/	
#
echo '... cd into folder data '
cd data
#
echo ' ... run simulation error and printouts are saved in the file out'
./a.out &> out
cp out ../
#
echo ' ... run the Matlab script plot_bomex_paper.m to visualize the simulation data '

