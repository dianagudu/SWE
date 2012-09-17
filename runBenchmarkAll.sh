#!/bin/bash

# set program arguments
blockX=8 # number of blocks in x direction
blockY=1 # number of blocks in y direction
nX=50	 # number of cells in x directions on the coarsest level
nY=4	 # number of cells in y direction on the coarsest level

# compile for benchmarking
scons netCDFDir=/usr/local writeNetCDF=True adaptive=True benchmark=True

# create config file
echo "1 2 1 1 1 1 1 1" > config

# create output directory
mkdir bench_${blockX}x${blockY}_${nX}x${nY}

# run all cases
./build/SWE_gnu_release_none_augrie $blockX $blockY $nX $nY build/out config GTS 0
./build/SWE_gnu_release_none_augrie $blockX $blockY $nX $nY build/out config GTS 1
./build/SWE_gnu_release_none_augrie $blockX $blockY $nX $nY build/out config GTS 2
./build/SWE_gnu_release_none_augrie $blockX $blockY $nX $nY build/out config LTS 0
./build/SWE_gnu_release_none_augrie $blockX $blockY $nX $nY build/out config LTS 1
./build/SWE_gnu_release_none_augrie $blockX $blockY $nX $nY build/out config LTS 2

