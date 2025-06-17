#!/bin/sh
#$ -N clearcreek_run
#$ -j y
#$ -cwd
#$ -o out.txt
####$ -m e
#$ -pe smp 8
#$ -q all.q


###mpirun -np 8 ../build/asynch_1.4.5 clearcreek.gbl
mpirun -np 2 ../build/src/asynch clearcreek.gbl
