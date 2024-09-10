#!/bin/bash


module load zlib/1.2.11_parallel_studio-2017.1
module load hdf5/1.8.18_parallel_studio-2017.1
module load openmpi/2.0.1_parallel_studio-2017.1  

# module load stack/2021.1 
# module load zlib/1.2.11_intel-2021.2.0
# module load hdf5/1.10.7_intel-2021.2.0-mpi
# module load openmpi/4.0.5_intel-2021.2.0


#../configure --prefix=/Users/nicolas/executables/asynch CFLAGS="-O3 -march=core-avx2 -DNDEBUG"
##../configure --prefix=/Users/nicolas/2022_EKI/exec/asynch CFLAGS="-O3 -march=core-avx2 -DNDEBUG"
../configure --prefix=/Users/nicolas/executables/tempmodel CFLAGS="-O3 -march=core-avx2 -DNDEBUG"


make
make install
