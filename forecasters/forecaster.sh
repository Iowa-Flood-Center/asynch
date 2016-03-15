#!/bin/sh
#$ -N Test
#$ -j y
#$ -cwd
####$ -m e
####$ -M your-email-address@uiowa.edu
####$ -pe orte 8
#$ -pe 8cpn 8
####$ -l mf=16G
#$ -l ib=1
####$ -q UI
#$ -q IFC
####$ -q all.q
####$ -q COE

/bin/echo Running on host: `hostname`.
/bin/echo In directory: `pwd`
/bin/echo Starting on: `date`
/bin/echo "Got $NSLOTS processors."
/bin/echo 
/bin/echo 
/bin/echo 


#Start an individual forecaster
#mpirun -np 8 ./FORECASTER_MAPS examples/GlobalForecast262_ifc1c.gbl examples/fcast_file.fcst

#Start a forecaster group
python forecaster_group_example.py 8

