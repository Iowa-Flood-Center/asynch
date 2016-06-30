# Installation

This section provides a quick guide for compiling and running ASYNCH on a Unix based system, including Iowa's local HPC clusters Helium and Neon ASYNCH will not work as-is in a Windows environment.

## Requirements

In addition to the source code, several programs and libraries are needed:

 * A C compiler
 * GNU Make
 * An MPI implementation
 * Libpq libraries

A brief description of each is provided below.

### C Compiler

ASYNCH has been successfully used and tested with the GNU C compiler gcc (version 4 1 2 and later) ASYNCH is compiled using the gnu99 standard Using another compiler (such as Intel's compilers) is an option, but some minor tweaks to the compile standard and/or included libraries may be required Information about gcc can be found at https://gcc.gnu.org/

### GNU Make

GNU Make is a utility for directing the generation of binaries from source code Make is available in most Linux repositories Details can be found at the project's website http://www.gnu.org/software/make/


### MPI Implementation

The Message Passing Interface (MPI) is a standard for transferring data on parallel computers ASYNCH uses MPI for communication between processes to perform calculations in parallel ASYNCH has been successfully used and tested with OpenMPI Details about OpenMPI can be found at http://www.open-mpi.org/

From the OpenMPI webpage:

>The Open MPI Project is an open source Message Passing Interface implementation that is developed and maintained by a consortium of academic, research, and industry partners. OpenMPI is available on Iowa's HPC clusters In theory, ASYNCH should work properly with any other implementation adhering to the MPI standard.

If installing an MPI implementation on a machine for ASYNCH, be sure to install the development packages of the MPI implementation In Linux repositories,these packages are usually denoted with a `-dev` or similar in the package name. If in doubt, try typing "mpirun" and "mpicc" in a terminal. Both of these should be present to run ASYNCH. If mpirun is not present, you have not installed the MPI binaries (meaning you probably haven't tried installing MPI at all). If mpicc is not present, then you are missing the development package.

### PostgreSQL Library

`libpq` is a library for communicating with a PostgreSQL database, created by the makers of PostgreSQL. The libpq website is http://www.postgresql.org/does/9.1/statie/libpq.html

From the `libpq` webpage:

>libpq is the C application programmer's interface to PostgreSQL libpq is a set of library functions that allow client programs to pass queries to the PostgreSQL back-end server and to receive the results of these queries
libpq MUST be installed to compile ASYNCH, even if the database features are not needed libpq is available on Iowa's HPC clusters

If installing libpq on a machine for ASYNCH, be sure to install the development packages In Linux repositories, these packages are usually denoted by a "-dev" or similar in the package name

### HDF5 Library



## Optional Software

Some additional software is available that may be useful, depending upon what the user wishes to do. These include

 * Git
 * dos2unix
 * NoMachine NX

### Git

Git is a distributed revision control system Although not needed for running ASYNCH, the source code repository does require Git for access Git is in most Linux repositories GitHub ofers information about usage

### dos2unix

This is a useful utility if editing input text fles for ASYNCH from a Windows machine Unix and Windows use a slightly diferent format for text documents. Although a fle may look the same under both a Linux and Windows text editor, subtle diferences can still exist In general, editing a text fle from Linux on a Windows machine will convert the fle to the Windows format To change the format to Unix, use the utility dos2unix If a fle is already under Unix format, this utility will not modify the fle Using a text fle in Windows format with ASYNCH will result in errors This process can be slow, depending upon the size of the text fles involved As such, ASYNCH does not automatically check the format of input text fles
dox2unix can be found in most Linux repositories

### NoMachine NX

NoMachine NX is a remote descktop that allows users to graphically login to a Linux system from a Windows, Linux, or Mac machine. This can be useful for users wishing to access an Iowa HPC resource, though this is not the only way Information about using and obtaining NoMachine NX for Iowa HPC resources can be found at [Helium Cluster Overview a Quick Start Guide]( https://www.iets.uiowa.edu/eonfluenee/display/ICTSit/Helium+Cluster+Overview+and+Quiek+Start+Guide).

## Source Code, Compiling, and Running ASYNCH

Note: Users of Iowa's HPC resources should NOT need to download and compile source code on Helium and Neon.

The ASYNCH source code is available in a repository hosted by GitHub Downloading the code from the repository requires the use of Git See Section 2 2 1 The source code can also be downloaded directly from GitHub as a zip file.

Once the source code is downloaded (and extracted from the zip fle, if needed), go to the directory where you downloaded the ASYNCH source code The ASYNCH library can be compiled with make ASYNCHLIB will create the library libs/libasynch.so in the ASYNCH directory, which contains the C interface routines. Using the command make ASYNCHLIB PY will compile the source code and create the library libs/libasynch py.so in the ASYNCH directory, which contains the Python interface routines See Section 9 2 for a list and description of these routines. Typing "make ASYNCH" will compile the basic simulation program If using your own computer, modifcations to the makefle may be necessary depending upon the location and names of the programs and libraries installed from Section 2 1

If the source code is ever updated, you may want to run `make clean` before recompiling. This removes all binaries and object fles of the old version. Once compiled, ASYNCH can be run with the command:

```
mpirun -np <number of processes> <path>/asynch < gbl flename>
```

## Iowa HPC Clusters

Currently, the University of Iowa has two HPC clusters available: Helium and Neon ASYNCH is already compiled on these clusters and is available to anyone with access to IFC shared resources All required software should be available The makefle included with the source code should work without modifcation on these clusters

However, these clusters do use OpenMPI through modules The module for OpenMPI must be loaded once per login session to run ASYNCH On Helium, this can be done with the command
module load openmpi gnu 1 4 3
On Neon, use the command
module load openmpi/gcc44-1 6 5
This loads OpenMPI version 1 4 3 or version 1 6 5 for use with the GNU compiler gcc (which was used to compile the existing version of ASYNCH) Instead of loading these modules manually, the commands can be added to the end of the fle .bashrc in the user's home directory Note that Helium and Neon each have a separate .bashrc fle In addition, if using the Python interface functions on Helium, the appropriate Python module must be loaded This can be done with a call to
module load python27
This can also be added to the .bashrc fle to automate the loading process
Binaries for ASYNCH are located in /Groups/IFC/Asynch/bin/ on Helium and Neon Libraries for using ASYNCH with custom designed software are located in the directory /Groups/IFC/Asynch/libs/

## Updating the package

```
autoreconf --install
make dist
```

## Installing the package

These are the generic instruction for an out of source build (prefered method):

```
mkdir build && cd build
../configure CFLAGS=-DNDEBUG
make
make check
make install
```

### Installing the package on NEON

First, `git clone` the repository or `tar xzf` a released packages. To install the software for the IFC group, load the following modules:

```
module load openmpi/intel-composer_xe_2015.3.187-1.8.8
module load hdf5/1.8.17
```

Then run the class GNU tool chain:

```
mkdir build && cd build
../configure --prefix=/Groups/IFC/.neon CFLAGS="-O2 -DNDEBUG" CHECK_CFLAGS=-I/Groups/IFC/.local/include CHECK_LIBS=/Groups/IFC/.local/lib/libcheck.a
make
make check
make install
```

## Updating the package

```
autoreconf --install
mkdir build && cd build
make dist
```

## Standard Makefile Targets


 - `make all` Build programs, libraries, documentation, etc. (Same as `make`.)
 - `make install` Install what needs to be installed.
 - `make install-strip` Same as `make install`, then strip debugging symbols.
 - `make uninstall` The opposite of `make install`.
 - `make clean` Erase what has been built (the opposite of `make all`).
 - `make distclean` Additionally erase anything `./configure` created.
 - `make check` Run the test suite, if any.
 - `make installcheck` Check the installed programs or libraries, if supported.
 - `make dist` Create PACKAGE-VERSION.tar.gz.
