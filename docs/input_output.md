# Input/Output Formats

Numerous inputs are required to construct and solve the model equations This includes link connectivity, link parameters, numerical error tolerances, precipitation forcings, a mathematical model, etc Most inputs are specified on a link-by-link basis, and can be given as either a file or a database table.

In addition, ASYNCH must know information about data outputs This includes where to store the outputs, what the outputs are, at which links outputs will be given, etc.

The format for each input and output is described below For database inputs/outputs, ASYNCH can access PostgreSQL through the libpq libraries Obviously, the user must have read or write permissions to a PostgreSQL database to utilize database features.

## Global File Structure

Global files (.gbl) are used to specify ALL inputs for the solvers. This includes river network topology files, database connections, initial model states, what information is printed to output files, model forcings, etc. Global files have a very rigid structure, unlike XML files, and information must be specified in a particular order. *The description of each input of the global files below are given in the order in which they are to be specified in the actual file*.

Global files are always ASCII files, assumed to be in UNIX format The percent sign (%) is used for single line comments As described below, certain inputs are expected to be given within a single line Other than this restriction, white space is ignored Arguments surrounded by { } below are mandatory, while those surrounded by the square brackets [ ] are dependent upon values specified by other parameters.

Although each setting in a global file modifies values for the ASYNCH solvers, their exact use can be modified by the user altering the underlying source code. This can be done with calls to routines described in Section

### Model Type and Maxtime

Format:

```
{model id} {total simulation time}
```

The first value specifies the id for the model to be used. This is a non-negative integer value which corresponds to a particular system of ordinary-differential equations (or possibly DAEs). Examples of built-in models is given in Section 8 3. For custom models, this model id is ignored See Section 9 4.
The second is the total simulated time used by the solvers This is a floating point number. The units are generally in minutes, however, some models can get away with using other units (seconds, hours, etc) if carefully crafted. All built in models use minutes for the time.

### Parameters on Filenames

Format:

```
{parameter on output filename flag}
```

This is a boolean value (0 or 1) that indicates whether all output filenames should include the uniform in space and time parameters 0 indicates no, 1 indicates yes. This feature can be useful for keeping track of output files from multiple simulations.

### Solver Outputs

Format:

```
{number of outputs}
[output1]
[output2]
[...]
```

This set of input parameters specifies the names of all outputs from the solvers. Several built in outputs exist, and the user is able to construct his own outputs. Built in outputs are given in Section 7 2. Output names are case sensitive. The first required value is the number of outputs (>= 0), followed by the names of each output, on separate lines.

### Peakflow Statistics Function Name

Format:

```
{function name}
```

This sets the function for outputting peakflow information. The built in peakflow function "Classic" is one option, and the user is free to construct his own. A function name must be specified here, even if peakflow information is not requested for any links.

### Global Parameters

Format:

```
{number of parameters} [parameter1] [parameter2]
```

This is where model parameters which are constant in space and time are specified The first value is a nonnegative integer specifying the number of global parameters to follow. Every model requires a certain number of global parameters If the number given in the global file is less than expected for a particular model, an error occurs If the number is greater than expected, a warning is given. These "extra" parameters are available to the model for use This can sometimes be useful for quick tests, but should be avoided normally.

The parameter meanings depend upon the model used. The units of these parameters is also model dependent.

### Buffer Sizes

Format:

```
{steps stored at each link} {max number of steps transferred} {discontinuity
buffer size}
```

These nonnegative integer values allow the user to specify sizes of internal buffers. In general, as these numbers are increased, the solvers run faster, but more memory is required. A good starting point that works in most cases is the set `30 10 30`. Typically, if these values need to be less than 10 to run the solvers, a deeper issue with memory constraints should be addressed.

### Topology

Format:

```
{topology flag} [output link id] { rvr filename or dbc filename}
```

This is where connectivity of the river network is specified. This can be done in one of two ways If the topology flag is 0, a river topology file (.rvr) is used. If the topology flag is 1, then topology is downloaded from the database specified with the database file (.dbc). The database connection allows for one additional feature: a subbasin can be specified If the output link id is taken to be 0, all link ids found in the database are used. Otherwise, the link with link id specified and all upstream links are used. Pulling subbasins from a topology file is not currently supported.

### Link Parameters

Format:

```
{parameter flag} { prm filename or dbc filename}
```

This specifies where parameters which vary by link and not time, are specified If the parameter flag is 0, the parameters are given in a parameter (.prm) file. If the flag is 1, then the parameters are downloaded from the database specified by the database connection file (.dbc). The number, order, meaning, and units of these parameters varies from model to model.

### Initial States

Format:

```
{initial state flag} { ini, uini, rec, or dbc filename} [unix time]
```

This section specifies the initial state of the model. The values for the initial state flag can be 0, 1, 2, 3 or 4 corresponding, respectively, to a ini,
 uini, rec, dbc, h5 file. The unix time argument is used for database connections only. This value is available in the query of the database connection file and can be used for selecting values from tables.

### Forcings

Format:

```
{number of forcings}
[forcing1 flag] [forcing1 information]
[forcing2 flag] [forcing2 information]
[...]
```

Information about time dependent forcings is specified here. Each model has an expected number of forcings. If the number of forcings specified here is less than expected, an error is thrown. If the number of forcings is greater than expected, a warning is given. This warning allows for tests to be performed and implemented quickly. In general, this feature should be avoided.

Forcing information varies considerably based upon the corresponding forcing flag. Several forcing types require unix times to determine what forcing data to use. If a model requires multiple forcings with unix times, the times do not need to be consistent, i.e., one forcing could start on July 1st 2014 at midnight, while another forcing starts at April 5th 2008.

#### No Forcing

Format:

```
0
```

A forcing flag of 0 specifies no forcing input This is the same as a forcing value of 0 0 for all links and all time.

#### Storm File

Format:

```
1 {.str filename}
```

A forcing flag of 1 indicates the forcing is specified by a .str file. The filename and path of a valid storm (.str) file is required.

#### Binary Files

Format:

```
2 {binary file identifier}
{chunk size} {time resolution} {beginning file index} {ending file index}
```

A forcing flag of 2 indicates the forcing is specified by a collection of binary forcing files. The identifier can be adorned with a path to the binary files. The chunk size is a positive integer that indicates the number of binary files kept in memory at once. The time resolution indicates the amount of time between successively indexed binary files. This value is a floating point number with units equal to those of the time variable of the model used The beginning and ending file indices indicate the range of the binary files to be used. The indices are integer valued.
The simulation will begin using the binary file with index given by the beginning file index. If the total simulation time would require binary files with index greater than the ending file index, the forcing values are taken to be 0.0 for all such binary files.

#### Forcings from Databases

Format:

```
3 {.dbc filename}
{chunk size} {time resolution} {beginning unix time} {ending unix time}
```

A forcing flag of 3 indicates the forcing data will be pulled from a PostgreSQL database. The database connection filename can include a path. The chunk size is a positive integer representing the number of forcing values pulled from the database at once from each link. A chunk size of 10 tends to work well. A larger chunk size requires more memory and larger datasets returned from the database, but a small number of queries. The time resolution is a floating point number with units in minutes. This represents the time resolution of the data in the accessed table. The integrity of the database table is not thoroughly checked by the solvers.

The simulation will begin using the data from the database with unix time given by the beginning unix time. If the total simulation time would require data from the database with unix time greater than the ending unix time, the forcing values are taken to be 0.0 for times greater than the ending unix time.

#### Uniform Forcings

Format:

```
4 {.ustr filename}
```

A forcing flag of 4 indicates a forcing that is uniform in space. The forcings are given by a uniform storm file (.ustr).

#### GZipped Binary Files

Format:

```
6 {gzipped binary file identifier}
{chunk size} {time resolution} {beginning file index} {ending file index}
```

A forcing flag of 6 indicates the forcing is specified by a collection of binary forcing files that have been gzipped (compressed as .gz files). All parameters for this flag are identical to that of using binary files with forcing flag 3.

#### Monthly Forcings

Format:

```
7 { mon filename}
{beginning unix time} {ending unix time}
```

A forcing flag of 7 indicates a uniform in space forcing that recurs monthly. When the end of the calendar year is reached, the monthly forcing file (.mon) is read again from the beginning The beginning unix time is used to determine the month the simulation begins (for this forcing). If the total simulation time takes the simulation past the ending unix time, the forcing is assumed to be `0.0` for all locations and times beyond the ending unix time

#### Grid Cell Files

Format:

```
8 {index filename}
{chunk size} {beginning file index} {ending file index}
```

A forcing flag of 8 indicates the forcing is specified by a collection of grid cell forcing files. The index filename can be adorned with a path to the index file. The chunk size is a positive integer that indicates the number of grid cell files kept in memory at once. The beginning and ending file indices indicate the range of the grid cell files to be used. The indices are integer valued.

The simulation will begin using the grid cell file with index given by the beginning file index. If the total simulation time would require grid cell files with index greater than the ending file index, the forcing values are taken to be 0 0 for all such grid cell files. In addition, if a grid cell file is missing, all values at each cell are assumed to be 0.0.

### Dams

Format:

```
{dam flag} [ dam or qvs filename]
```

This section specifies whether dams will be used A dam flag of 0 means no dams are used. A flag of 1 indicates a dam file ( dam) will be used, and a flag value of 2 indicates a discharge vs storage file ( qvs) will be used. Some models do not support dams. For these models, the dam flag must be set to 0 or an error occurs.

### State Forcing Feeds

Format:

```
{reservoir flag} [ rsv or dbc filename] [forcing index]
```

This section specifies whether a provided forcing (Section 6 1 10) is to be used as a forcing of the states of differential or algebraic equations at some links. A reservoir flag of 0 indicates no forcing will by applied to system states. A flag of 1 indicates state forcings will be applied to all link ids in the specified .rsv file. A reservoir flag of 2 indicates state forcing will be applied to all link ids pulled from the database the given .dbc file. If the reservoir flag is not 0, then the index of the forcing must be specified.

### Time Series Location

Format:

```
{time series flag} [time resolution] [ dat / csv / dbc filename] [table name]
```

This section specifies where the final output time series will be saved. A time series flag value of 0 indicates no time series data will be produced. Any flag with value greater than 0 requires a time resolution for the data. This value has units equal to the units of total simulation time (typically minutes). A value of -1 uses a resolution which varies from link to link based upon the expression:


where A is the upstream of the link, measured in km2.

A time series flag of 1 indicates the results of the simulation will be saved as a .dat file. The filename complete with a path must be specified. If a file with the name and path given already exists, it is overwritten. A time series flag of 2 indicates the results will be stored as a csv file. A time series flag of 3 indicates the results will be uploaded into the database described by the given .dbc file. In this case, a table name accessible by the queries in the .dbc file must be specified.

This section is independent of the section for Link IDs to Save described below (Section 6 1 15) For example, if link ids are specified in the Link IDs to Save section and the time series flag in the Time Series Locations set to 0, no output is generated. Similarly, if *the time series id flag* is set to 0 in the Link IDs to Save section and the time series flag is set to 1, a .dat file with 0 time series is produced.

**Notice: the time resolution is entirely independent of the time step used by the numerical integrators. Reducing this value does NOT produce more accurate results. To improve accuracy, reduce the error tolerances described in Section 6.1.19. There is no built-in way to produce results at every time step, as this is a very easy way to crash a compute node or file system.**

### Peakflow Data Location

Format:

```
{peakflow flag} [.pea / .dbc filename] [table name]
```

This section specifies where the final peakflow output will be saved. A peakflow flag of 0 indicates no peakflow data is produced. A peakflow flag of 1 indicates the peakflow results of the simulation will be saved as a .pea file. The filename complete with a path from the binary file must be specified. A peakflow flag of 2 indicates the results will be uploaded into the database described by the given .dbc file. In this case, a table name accessible by the queries in the dbc file must be specified.

This section is independent of the section for Link IDs to Save described below (Section 6 1 15). For example, if link ids are specified in the Link IDs to Save section and the peakflow flag in the peakflow Data Location is set to 0, no output is generated. Similarly, if the peakflow id flag is set to 0 in the Link IDs to Save section and the peakflow flag is set to 1, a .pea file with 0 peakflows is produced.

### Link IDs to Save

Format:

```
{time series id flag} [.sav / .dbc filename]
{peakflow id flag} [.sav / .dbc filename]
```

This section provides the list of link ids in which data is produced. The first line is for the time series outputs, while the second is for the peakflow outputs. The time series ID flag and the peakflow ID flag take the same list of possible values. A flag of 0 indicates no link IDs for which to produce data. A flag of 1 indicates the list of link IDs is provided by the corresponding save file (.sav). A flag of 2 indicates the list of link IDs is provided by the database specified in the given database connection file (.dbc). A flag of 3 indicates that all links will have data outputted.

**Warning: a time series ID flag of 3 can easily wreak havoc on a file system for simulations with a large number of links. At the very least, extremely large output files and database tables will occur. Be very careful with this!
Typically, using a flag value of 3 for peakflow link ids, or for the time series ID flag for a very small basin (< 500 links) will not create any problems.**

This section is independent of the sections for Time Series Location and peakflow Data Location above (Sections 6 1 13 and 6 1 14). For example, if link ids are specified in the Link IDs to Save section and the time series flag in the Time Series Location set to 0, no output is generated. Similarly, if the time series id flag is set to 0 in the Link IDs to Save section and the time series flag is set to 1, a .dat file with zero time series is produced.

### Snapshot Information

Format:

```
{snapshot flag} [.rec / .dbc / .h5 filename]
```

This section specifies where snapshot information is produced A snapshot is a record of every state at every link in the network Snapshots are produced at the end of simulations. This is useful for beginning a new simulation where an old one ended. A snapshot flag of 0 indicates no snapshot is produced. A snapshot flag of 1 indicates the snapshot will be produced as a recovery (.rec) file with path and filename specified. A snapshot flag of 2 indicates the snapshot will be uploaded to the database specified by the database connectivity (.dbc) file. A snapshot flag of 3 indicates the snapshot will be produced as a HDF5 (.h5) file with path and filename specified.

### Scratch Work Location

Format:

```
{filename}
```

This section specifies the location of temporary files. These files are used to store intermediate calculations. The filename can include a path name. If the file already exists, the contents are overwritten. If a simulation is aborted, these files may not be removed. Otherwise, they are deleted at the end of the simulation.

### Error Control Parameters

Format:

```
{facmin} {facmax} {fac}
```

This section specifies parameters related to the error control strategy of the numerical integrators. The value facmin represents the largest allowed decrease in the stepsize of the integrators as a percent of the current step Similarly, facmax represents the largest allowed increase. The value fac represents the safety factor of the integrators Any accepted stepsize is multiplied by this value Good values of facmin, facmax, and fac to use are 0 1, 10 0, and 0 9, respectively

### Numerical Error Tolerances

Format:

```
{solver flag} [ rkd filename]
[rk solver index]
[absolute error tolerance 1] [absolute error tolerance 2]
[relative error tolerance 1] [relative error tolerance 2]
[dense absolute error tolerance 1] [dense absolute error tolerance 2]
[dense relative error tolerance 1] [dense relative error tolerance 2]
```

This section specifies error tolerances for the numerical integrators. A solver flag of 0 indicates the same tolerances will be used for all links. A solver flag of 1 indicates the tolerance info will be specified in the given RK data (.rkd) file.
If solver flag is 0, than an rk solver index must be specified. A list of Runge-Kutta methods is given in Section 7 4. Each error tolerance must have a value for each state of the system. The order of the tolerances must match the order of the states in the state vectors. The absolute and relative error tolerances are those typically used for RK methods. The dense tolerances are for the numerical solution produced between time steps. A numerical solution is rejected if either the error tolerances or dense error tolerances for any state is believed to be violated.

## Database Connection Files

Database connection files are ASCII text files with a .dbc extension which specify how to connect to a database, and the queries to pull/push data from/to the database. Although the format of database connection files is the same, the specific requirements of the queries varies with how the file is used For instance, queries for pulling link connectivity information is very different from queries for uploading peakflows to a database table.

Format:

```
dbname={db} host={host} user={user} password={pass}
{number of queries}
[query 1]
[query 2]
...
```

The first line of every database connection file specifies the information needed to make a connection. A user must have permission to read or write from the database at the given host; otherwise, queries sent to the database will fail The number of queries will vary depending upon how the database connection file is used. The appropriate number of queries and what they should return is documented in the remainder of Section 6. The number of queries may be zero.

Any queries listed *MUST* be ended with a semicolon (;). For some queries, further information outside the database connection file may be available, depending upon how the query is to be used. This additional information is explained in the appropriate section below for input formats. Such information includes link ids and unix times. To denote in a query where this information should be placed, use the symbol "%u" for integers and "%s" for names.

## Link Connectivity Input

Link connectivity information is used to specify how links in a network are connected Connectivity can be provided through either a river network file (.rvr) file or through a database table. When a river network file is used, every link in the file is used (i.e. no subnetworks) then pulling connectivity data from a database, a subset of the network can be used.

Regardless of the format, all link ids must be given in the same order in the link connectivity, link parameter, and initial state inputs.

River network files are ASCII text files with the following format:

```
{number of links}
{link id 1}
{number of parents} [parent id 1] [parent id 2]
{link id 2}
{number of parents} [parent id 1] [parent id 2]
```

White space can be used freely throughout the file. The layout in the above specification is purely optional; the order of the information is what is important. The file begins with the total number of links in the file. Then each link id is specified, followed by the number of parents for the link and each of their ids. A link id can appear in a list of parent link ids at most once. If a link does not have parents, it must still appear in this file with a 0 for the number of parents.

If the connectivity is pulled from a database, a corresponding database connection file is used This file requires three queries:

 - Query to pull all link ids from a table
   - Inputs: none
   - Returned tuples: (link id)
 - Query to pull all link id, parent link id pairs
   - Inputs: none
   - Returned tuples: (link id, parent link id)
 - Query to pull all link id, parent link id pairs upstream from a given outlet link id
   - Inputs: outlet link id
   - Returned tuples: (link id, parent link id)

The last two queries must return a tuple for each link id for each parent link. So a link with two parents should appear twice in the returned tuples, once for each parent link. The returned tuples must be grouped by the link id so all parent information appears consecutively.

## Link Parameter Input

Link parameter input specifies the parameters for the model that vary link to link This information can be provided in a parameter file ( prm) or through a database table The number of parameters for each link, their meaning, and their order depends upon the model used In particular, the value of disk params determines the number of parameters expected at each link See Section 8 1 1.

Regardless of the format, all link ids must be given in the same order in the link connectivity, link parameter, and initial state inputs.

A parameter file is an ASCII text file with the following format:

```
{number of links}
{link id 1}
{parameter 1} {parameter 2} {parameter 3}
{link id 2}
{parameter 1} {parameter 2} {parameter 3}
```

White space can be used freely throughout the file. The layout in the above specification is purely optional; the order of the information is what is important. The file begins with the total number of links. Then each link id is specified, followed by the parameters for that link.

If the parameters are pulled from a database, a corresponding database
connection file is used  This file requires two queries:

 - Query to pull all parameters
   - Inputs: none
   - Returned tuples: (link id, parameter 1, parameter 2,  )
 - Query to pull all parameters above an outlet
   - Inputs: outlet link id
   - Returned tuples: (link id, parameter 1, parameter 2,  )

## Initial Values Input

The link initial values input specifies the initial values for the states of the differential and algebraic model equations This information can be provided in several different formats: an initial value file (.ini), a uniform initial value file (.uini), a recovery file (.rec), and through a database table.

An initial value file is an ASCII text file that lists the initial values for each link. The format is:

```
{model type}
{number of links}
{initial time}
{link id 1}
{initial value 1} {initial value 2}
{link id 2}
{initial value 1} {initial value 2}
```

The model type is the number of the model to be used. This determines how many initial values are expected for the model. Initial states must be provided only for those states determined by differential equations, and only for those which require an initial condition. These are the states with index between `diff_start` and `no_ini_start` in the state vectors See Section 8 1 1.

A uniform initial value file is similar to an initial value file, but the initial values, when required, are the same at every link The format is given by:

```
{model type}
{initial time}
{initial value 1} {initial value 2}
```

The model type is the number of the model to be used. This determines how many initial values are expected for the model. Initial values must be provided only for those states determined by differential equations, and only for those which require an initial condition. These are the states with index between `diff_start` and `no_ini_start` in the state vectors. See Section 8 1 1. Notice that unlike an initial value file, no link ids are given, and only one set of initial values are given.

A recovery file is an ASCII text file that lists the initial values for each link. The format is:

```
{model type}
{number of links}
{initial time}
{link id 1}
{initial value 1} {initial value 2}
{link id 2}
{initial value 1} {initial value 2}
```

The format is identical to that of an initial value file, with one important exception The initial value of EVERY state must be provided at each link For models with `diff_start` set to 0and `no_ini_start` set to dim, a recovery file is identical to an initial value file See Section 8 1 1 Warning: For the initial values of algebraic equations, no checks on the input data are performed to ensure the solution is consistent.

If the initial values are pulled from a database, a corresponding database connection file is used. This file requires one query:

 - Query to pull all initial states for every link:
   - Inputs: integer value
   - Returned tuples: (link id, initial value 1, initial value 2, )

The query allows for one input to be used to obtain the needed information. This value could be, for example, an outlet link id or a unix time. Similar to recovery files, initial values must be provided for every link.

## Forcing Inputs

Numerous and diverse formats are implemented for feeding forcing inputs into a model. These formats vary considerably, and can have different impacts on performance.

### Storm files

Storm files (.str) provide an ASCII text format for setting a forcing at each link. The format of these files is:

```
{number of links}
{link id 1} {number of changes}
{time 1} {value 1}
{time 2} {value 2}
{time 3} {value 3}
{link id 2} {number of changes}
{time 1} {value 1}
{time 2} {value 2}
{time 3} {value 3}
```

The format requires a time series to be provided for every link. The number of values can vary from link to link, and the time steps do not need to be uniformly spaced or even the same for each link. The first time at each link must be the same however, and correspond to the beginning of the simulation (typically 0). The forcings are assumed to be constant between time steps. After the final time step, the forcing value is held at the last value for the remainder of the simulation. The data provided by a storm file is entirely read into memory at the beginning of a run. As such, this format may not be suitable for large or long simulations.

### Uniform storm files

Uniform storm files (.ustr) provide an ASCII text format for setting a forcing uniformly at each link. The format of these files is:

```
{number of changes}
{time 1} {value 1}
{time 2} {value 2}
{time 3} {value 3}
```

The format is similar to that of a storm file, but only one time series is specified, and is applied at every link. The time steps do not need to be uniformly spaced. The first time must correspond to the beginning of the simulation (typically 0) The forcing is assumed to be constant between time steps. After the fnal time step, the forcing value is held at the last value for the remainder of the simulation. The data provided by a uniform storm file is entirely read into memory at the beginning of a run. As such, this format may not be suitable for extremely long simulations.

### Binary storm files

Binary storm files (no extension) may also be used. Instead of a single file, these are a collection of files providing forcing values at different times. The format of these files is:

```
{link id 1 value}
{link id 2 value}
{link id 3 value}
```

Each file is simply a list of forcing values for each link. Because link ids are not present, the values are assumed to be in the same order as the link ids from the link connectivity input. Further, a value must be specified for every link.
The filename of each binary file should end with an integer value. The index of each file should be numbered consecutively. The time step size between files, as well as the range of files to use, are specified in the global file (see 6 1 10). If the simulation goes beyond the last file, all further forcing values are assumed to be 0 for all links.
The values in binary storm files are stored as single precision floating point numbers, with endianness different from the native machine. These files are read into memory chunk by chunk. This allows for a ceiling on the memory used, independent of the number of files.

### Gzipped binary storm files

Gzipped binary storm files (.gz) are also supported. The format of these files is identical to that of the binary storm files, as they are simply gzipped binary storm files. The files are uncompressed before use, making this format slower than the regular binary storm files.

Grid cell files group link ids together, and specify forcing values for each group (referred to as a cell). Although similar to binary files, this format requires two accompanying text files: an index file and a lookup file.

The index file specifies meta data for the grid cell data. The format for this ASCII file is:

```
{time resolution (mins)}
{conversion factor}
{number of cells}
{grid cell data file name prefix}
{lookup filename}
```

The time resolution specifies the amount of time between grid cell files. The resolution is typically given in minutes. The conversion factor is a floating point number. Each forcing value in the grid cell files is multiplied by this factor. The number of cells is an integer value. Each grid cell filename has a prefix, followed by a number. The prefix is specified in the index file. The prefix may include a path. If the path is relative (i e , does not begin with a '/'), the path is taken relative to the location of the index file. Lastly, the index file includes the filename for the lookup file. A path is allowed, but is taken relative to the location of the index file, unless given as an absolute path.

The lookup file specifies which link ids are located in each cell. The format for this text file is:

```
{link id 1} {cell id 1}
{link id 2} {cell id 2}
```

The cell ids are indexed starting from 0, and the cell index cannot be larger than the number of cells specified in the accompanying index file.

### Grid cell files

The grid cell files are binary files. Each gives forcing values for each cell at a moment of time. If a cell is omitted from a file, then the forcing value is assumed to be 0. The filename for each grid cell file begins with the same prefix, followed by an integer. This integer indicated the order in which the grid cell files should be used. Although the number on these files should be consecutive, a missing file indicates all cells take a forcing value of 0. The format of the binary grid cell files is:

```
{1} {cell id 1} {forcing value 1}
{cell id 2} {forcing value 2}
```

The grid cell files are binary, so the spacing above is purely for readability. Each file begins with the integer value 1, stored as a 32-bit integer. This is used for checking the file format. Each cell id is given as a 32-bit integer and each forcing value is given as a 16-bit integer. Before the forcing values are loaded into memory, they are multiplied by the conversion factor given in the index file. Again, every cell does not need to be given in every grid cell file; only when the forcing value is nonzero does a value need to be given.

### Monthly recurring forcing

Monthly recurring forcing files (.mon) allow forcings to be set monthly. These files are given as ASCII text files in the format:

```
{value for January}
{value for February}
{value for March}
...
{value for December}
```

A value must be given for each month. Over each simulated month, the forcing value is held constant, and is uniform over all links.

### Database forcing

If the forcing data is pulled from a database, a corresponding database connection file is used. This file requires three queries:

 - Query to pull rainfall data for all link ids present in the database table
   - Inputs: lower unix time, upper unix time -Returned tuples: (unix time, forcing value, link id)
 - Query to pull rainfall data for all link ids upstream from an outlet link
   - Inputs: outlet link id, lower unix time, upper unix time -Returned tuples: (unix time, forcing value, link id)
 - Query to pull a single forcing intensity from the database table
   - Inputs: none
   - Returned tuple: (unix time)

The first and second queries are similar, except that the second pulls only a subset of links from the database table. Forcing values are read in blocks from the table, with forcing values coming between the lower (inclusive) and upper (exclusive) unix times available to the queries. If a value for a link is not pulled from the database, the value is assumed to be 0.

The last query is used to find an actual valid timestamp in the database table. It needs to return only one unix time.

## Dam Parameters Input

Two formats currently exist for setting parameters at links with dams: dam parameter files (.dam) and discharge vs storage files (.qvs). The format of dam parameter files is similar to that of parameter files:

```
{number of links with dams}
{link id 1}
{parameter 1} {parameter 2} {parameter 3} ...
{link id 2}
{parameter 1} {parameter 2} {parameter 3} ...
...
```

The number of parameters needed for each link is model dependent and determined by the value dam params size. See Section 8 1 1. For dam parameter files, only the links with dams must be listed here. Only links with id appearing in this file will have dams.

Discharge vs storage files take a series of discharge values and a corresponding series of storage values to decide the relationship between two states. The format of these files is similar to storm files (see Section 6 6):

```
{number of links with dams}
{link id 1} {number of points for link 1}
{storage value 1} {discharge value 1}
{storage value 2} {discharge value 2}
...
{link id 2} {number of points for link 2}
{storage value 1} {discharge value 1}
{storage value 2} {discharge value 2}
...
```

The number of points at each link can vary For dam parameter files, only links with dams are listed here. Only links with id appearing in this file will have dams. Internally, discharges and storages with values between two points are interpolated. This interpolation process is model dependent.

## Time Series Output

Three formats are supported for outputting time series calculations: data files (.dat), comma-separated values (.csv), and a database table. The particular time series calculated is set in the global file (see Section 6 1 13). The structure of each format is considerably different.

Data files are in ASCII text format. These files are designed to be generic and flexible so as to be easily read by whatever data analysis program the user prefers. Data files are created with the format:

```
{number of links}
{number of output values}
{link id 1} {number of points for link id 1}
{value 1 for series 1} {value 1 for series 2} {value 1 for series 3} ...
{value 2 for series 1} {value 2 for series 2} {value 2 for series 3} ...
...
{link id 2} {number of points for link id 2}
{value 1 for series 1} {value 1 for series 2} {value 1 for series 3}
{value 2 for series 1} {value 2 for series 2} {value 2 for series 3}
...
```

The series for the links appear in a column The number of points can vary from link to link, depending upon the user's selection in the global file The number of output values determines how many values appear in each line of the time series

A CSV file is a typical format to make data easy to read in spreadsheet software. The structure of CSV files is:

```
{link id 1},, ... , {link id 2},
Output 1, Output 2,  , Output 1, Output 2,   
{value 1,1,1},{value 1,2,1},  , {value 1,1,2},{value 1,2,2},   
{value 2,1,1},{value 2,2,1},  , {value 2,1,2},{value 2,2,2},   
{value 3,1,1},{value 3,2,1},  , {value 3,1,2},{value 3,2,2},   
```

The series for the links appear in a row. Under link id 1, each requested series appears, followed by the series for link id 2, and so on.

A database connection file can be used to upload results into a database table This file requires only one query:
 - Query to create a table for uploading data
   - Inputs: table name
   - Returned tuples: none

The query should create the table where the series information is to be stored ASYNCH does NOT remove any existing data from the table, or check if the table exists already.

## Peakflow Output

Peakflow outputs can be created in two formats: peakflow files (.pea) and database tables.

Peakflow files created with the "Classic" peakflow function take the structure:

```
{number of link}
{model type}
{link id 1} {peakflow value} {time to peak} {area}
{link id 2} {peakflow value} {time to peak} {area}
{link id 3} {peakflow value} {time to peak} {area}
```

The time to peak is measured since the beginning of the simulation. The peakflow value for each link is the maximum value achieved over the simulation for the state with index 0 in the state vector. The area given is the parameter from the link parameters with index area idx. See Section 8 1 1.

Peakfow output may be written to a database table if a database connection file is specified. One query is required, and one additional query is optional:
 - Query to create a table for uploading data
   - Inputs: table name
   - Returned tuples: none
 - Query to delete contents from a table
   - Inputs: table name
   - Returned tuples: none

The first query should create the table where the peakflow information is to be stored ASYNCH does NOT remove any existing data from the table, or check if the table exists already. The second query is optional, and will be used to delete any existing contents from the table before uploading data. The particular values uploaded to the database are determined through the peakflow function defined in Section 6 1 4.

## Link IDs for Time Series and Peakfows

Link ids must be specified for time series output and peakflow outputs. This can be done in one of two formats: save files (.sav) and database tables. Each of these formats is effectively just a list of link ids.

The structure of save files is:

```
{link id 1}
{link id 2}
{link id 3}
...
```

If a link id is specified in a save file, but is not present in the network, a warning will be issued, and the link id is ignored.

For pulling links from a database table, only one query is required:

 - Query to pull link ids
   - Inputs: none
   - Returned tuples: (link id)

## Snapshot Output

Snapshot outputs can take two formats: recovery files and database tables. The format for recovery files is covered in Section 6 5 as an input.

For using a database table, a database connection file is specified The database connection file has three optional queries:

 - Query to create a database table before the first upload
   - Inputs: table name
   - Returned tuples: none
 - Query to delete a table before the first upload
   - Inputs: table name
   - Returned tuples: none
 - Query to truncate a table before every upload
   - Inputs: table name
   - Returned tuples: none

In practice, snapshots are often applied periodically as a way to create check points for the program. The third query allows the user to limit the number of snapshots in a table to one.

## Runge-Kutta Data Input

Runge-Kutta Data files (.rkd) allow information about the underlying numerical methods to be specified link by link. These files are ASCII. The structure is given by:

```
{number of links}
{number of states}
{link id 1}
[absolute error tolerance 1] [absolute error tolerance 2]
[relative error tolerance 1] [relative error tolerance 2]
[dense absolute error tolerance 1] [dense absolute error tolerance 2]
[dense relative error tolerance 1] [dense relative error tolerance 2]
{RK Index for link id 1}
{link id 2}
[absolute error tolerance 1] [absolute error tolerance 2]
[relative error tolerance 1] [relative error tolerance 2]
[dense absolute error tolerance 1] [dense absolute error tolerance 2]
[dense relative error tolerance 1] [dense relative error tolerance 2]
{RK Index for link id 2}
...
```

An error tolerance is specified for every state at every link. The order of the links must match with the order given by the topology input, and number of states must agree with what the model expects

## Temporary Files

In general, sufficient memory is not available to hold a history of states while further computations take place. Thus, temporary files are created by ASYNCH to store time series outputs. These files are aggregated (typically after the simulation has completed) into the final time series output files (see Section 6 8).

Most users will not need to concern themselves with the underlying structure of these files. However, some routines exist for manipulating how these files are used, and an understanding of the temporary file structure is necessary for using those routines.

Temporary files are always in binary format. Every MPI process has its own file, which contains time series outputs for links assigned to the process. The format of the data looks like:

```
{link id 1} {expected number of steps}
{step 0, value 0} {step 0, value 1}
{step 1, value 0} {step 1, value 1}
{step 2, value 0} {step 2, value 1}
{link id 2} {expected number of steps}
{step 0, value 0} {step 0, value 1}
{step 1, value 0} {step 1, value 1}
{step 2, value 0} {step 2, value 1}
...
```

Because these files are a binary format, the presentation above is simply for readability. No new lines or spaces are present in the actual files. Only link ids for which time series data has been requested will appear in the files. Before any data is written, dummy values are placed in the files when they are created to insure the files are large enough. The number of these dummy files for each link is given by the expected number of steps value. This number is determined based upon the values of `maxtime` and the time resolution of the output time series when the temporary files are created. **Warning: Modifcations to these values after creation of the temporary files could create a situation where the files are not large enough to contain every point in the time series.** This generally happens when `maxtime` is increased. Therefore, when the temporary files are created, the largest expected value of `maxtime` should be set. If the temporary files are not large enough, a warning will be displayed during calculations, and further outputs will be lost.

While calculations are performed, processes will write output data to these files, replacing whatever dummy values are present. Modifying the behavior of these files is generally not needed, but can be performed through various routines. See Section 9 2.
The link ids and expected number of steps in temporary files are represented by unsigned 32 bit integers. The data types and sizes for the time series data will vary depending upon the desired outputs. If a built-in output time series is used for the states, these will appear in the temporary files as 64 bit floating point numbers.
