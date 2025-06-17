#!/bin/bash


module load zlib/1.2.11_parallel_studio-2017.1
module load hdf5/1.8.18_parallel_studio-2017.1
module load openmpi/2.0.1_parallel_studio-2017.1   

../configure --prefix=/Users/nicolas/Tiles/dist CFLAGS="-O3 -march=core-avx2 -DNDEBUG"


make
make install
