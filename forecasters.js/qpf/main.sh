#!/bin/bash
source /etc/profile

PATH=$PATH:$HOME/.local/bin:/Groups/IFC/.local/bin

# Load module OpenMPI
module load openmpi/intel-composer_xe_2015.3.187-1.8.8

# Set current working dir to the directory of this script
cd "$(dirname "$0")"

# Run the generator
node $1