# Quickstart Guide and Example Simulations

In this section, we will give the steps and data needed for running basic simulations This section is intended for getting a new user up and running quickly, and to provide links in this documentation for more information when needed.

First, be sure to follow the steps in Section 2 3 for installing the ASYNCH solvers For those using Iowa HPC resources, downloading and compiling the source code is not necessary.

The ASYNCH directory contains a folder called examples, which contains several data files for starting sample simulations, as well as sample outputs for comparision For Iowa HPC users, the ASYNCH directory is

```
/Groups/IFC/Asynch/
```

The `examples` directory should be copied to a location where the user has write access (for example, the home directory). On Helium or Neon, this can be done with

```
cp -r /Groups/IFC/Asynch/examples/ ~/
```

For the first example, we will produce output for a small basin with 11 links using a hydrological model with constant runoff. The global file to setup this simulation is `Global190.gbl`. This uses the model given in Section 8 3 1. If using your own machine, the simulation can be run with the command

```
mpirun -np 1 <bin path>/ASYNCH Global190 gbl
```

If using Iowa HPC resources, use the appropriate binary path:

```
/Groups/IFC/Asynch/bin helium/
/Groups/IFC/Asynch/bin neon/
```

As calculations are performed, you will see output produced to the terminal window If using Helium or Neon (or any system using the Sun Grid Engine), the submit script test.sh can be used to run the simulation Use the command

```
qsub test sh
```

while in the examples directory to launch the job Output from the program will be produced in a log file with a name like test run.o###### Try using the command

```
qstat -u <username>
```

to monitor the progress of your job.

Warning: A submit script is needed to run a job on multiple machines of Iowa HPC resources. If you attempt to run an ASYNCH simulation using just `mpirun` at a terminal window, you are probably running ASYNCH on a login node. Doing this limits the number of cores available to 12, slows down all other users's connections, and is an easy way to be reported to the HPC admins for misuse of resources!

When the program is complete, the output results are written to the folder `examples` The global file causes the production of three output fìles: `test.dat`, `test.pea`, and `test.rec`. These files should be identical to those found in `examples/results`. The dat file contains the output hydrograph for links with link ids 3 and 80 The pea file contains the peakfow information for every link The rec file contains the fnal value of every state of every link at the end of the simulation For this simulation, all output files are small enough to view in a text editor.

The simulations performed will use only 1 MPI process To increase this number, use, for example,

```
mpirun -np 2 <ASYNCH directory>/ASYNCH Global190 gbl
```

or modify `test.sh` to use more processes This can be done by modifying

```
#$ -pe orte 1
```

to use 2 processes instead of 1 Also be sure to modify the last line with mpirun so MPI looks for 2 processes.
When using more than 1 process, your results may difer slightly from those in `examples/results`. In fact, the results may vary slightly from simulation to simulation, even if nothing changed in the global file. This is a result from the asynchronous communication used by ASYNCH for MPI processes and is an expected behavior.

As a second example, try the same procedure as before using the global file `Global254.gbl`. If using an Iowa HPC resource, the submit script `clearcreek.sh` can be used. The model for this simulation is the toplayer hydrological model using the Clear Creek river basin See Section 8 3 2. Results for the output discharge and basefow are given in Figure 3 This basin is larger than in the previous simulation as it contains about 6,000 links. This is a good example to experiment with the number of processes used A time series of the channel discharge and basefow at the outlet are given in Figure 3
