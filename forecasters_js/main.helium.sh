#!/bin/bash
source /etc/profile

PATH=$PATH:/Groups/IFC/.local/bin:/Groups/IFC/.helium/bin

# Load module OpenMPI
module load openmpi_intel_1.6.3
#module load hdf5/1.8.17

# Set current working dir to the directory of this script
cd "$(dirname "$0")"

# Run the generator (if it not already running)
export DEBUG=forecaster,sge
/usr/bin/flock -n lock/$1 node --max-old-space-size=16384 $1
