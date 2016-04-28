# Non persistent Forecasters

A forecast is a not so special type of ASYNCH style simulation. The simulation begins with a rainfall forcing for a period of time (e.g. 15 hours for QPF). This forcing is typically some type of rainfall forecast or scenario.

Forecasters ASYNCHPERSIS and FORECASTER_MAPS are programs that runs continuously (persistently) that is it waits actively for new data to consume between solver runs. That can be a problem on shared infrastructure like HPC clusters where nodes/cores are allocated for a given time to perform a job (that is not supposed to run forever).

Forecasters.js are scripts that prepare and submit a regular ASYNCH simulation. Since it is not computationally intensive (sequential work) it typically runs on login nodes. It can also be ran at a given time interval using cronjobs. Basically the scripts performs the following tasks:

 - Check if new data are available
 - Check the current status of a similar job if there is one in the queue
 - Prepare a simulation
  - fetch the data from the DB (eg write stormfiles)
  - render a global file 
  - render a job file
 - Submit a simuation to the grid
 
This approach has advantages:
 - adapted to real time systems where latency matters (a large amount of processing is required for a short period of time)
 - less sensitive to the power loss or reboot: the cronjobs restart with the system, the state of the system is saved to files at regular interval
 
## Cronjobs

Here is the current crontab configuration:

```
cd /Users/sdebionne/asynch/forecasters.js
*  * * * * /main.sh qpe.js    >> /log/qpe.txt 2>&1
40 * * * * /main.sh qpf.js    >> /log/qpf.txt 2>&1
40 * * * * /main.sh whatif.js >> /log/whatif.txt 2>&1
```

In plain english:
 - QPE check if new QPE are available every minutes
 - QPF and WHATIF scenario are run at minute 40 every hours
 
Since crontabs are run in a nacked environment, the `main.sh` set the PATH and load the [Environment Modules](https://wiki.uiowa.edu/display/hpcdocs/Environment+Modules) required for the simulation.
 
## Scripts

### QPE.js

`qpe.js` is a script that runs at the time step of the observation and compute the most state of the domain according to the observation. This current state is used as initial conditions for the forecasts.

Each run of ASYNCH load intial conditions from the previous QPE run and save the final conditions to a .rec file. The files are timestamped so that the systems knows what is the current time of the simulation.

## QPF and WHATIF

`qpf.js` and `whatif.js` are forecasters generate forecast or outlooks respectively. They used the latest state available and run a simulation up to 10 days in the future with the current configuration.

## Interaction with SGE

Submitting a job to the grid adds to the queue of waiting jobs for some time. Since 'some time'is rather fuzzy, we have to deal with the different situation that may arise when preparing a job:

 - A similar job is already running
 - A similar job is already waiting
 - A job requires an other job to be finished before running
 
The stategy used is:
 - Job already running
  - Wait until it is finished (postpone the job preparation)
 - Job already waiting 
  - put the job on hold
  - prepare the new job configuration (that will be run by the job waiting)
  - release the job
 - Job depedencies
  - QPF and WHATIF are dependent on QPE (if a QPE sim is running, the jobs will wait until it is finished)
 
## Debugging the forecasters

To print debugging trace in the console:

```
export DEBUG=forecaster,sge
```
