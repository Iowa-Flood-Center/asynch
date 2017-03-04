#ifndef ASYNCH_INTERFACE_H
#define ASYNCH_INTERFACE_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <stdio.h>
#include <stdbool.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif 

//#include "comm.h"
//#include "riversys.h"
//#include "processdata.h"
//#include "structs.h"
//#include "solvers.h"
//#include "io.h"
#include "data_types.h"
#include "mathmethods.h"


//#define ASYNCH_MAX_DB_CONNECTIONS 20

// Forward definitions
typedef struct ErrorData ErrorData;
typedef struct GlobalVars GlobalVars;
typedef struct Link Link;
typedef struct RKMethod RKMethod;
typedef struct TransData TransData;
typedef struct TempStorage TempStorage;
typedef struct ConnData ConnData;
typedef struct QVSData QVSData;
typedef struct Forcing Forcing;
typedef struct AsynchSolver AsynchSolver;

typedef struct VEC VEC;
typedef struct MAT MAT;


//Callback signatures

/// 
///
/// \param t The current time.
/// \param y The current state vector at time *t* for the link.
/// \param global_params The vector of parameters uniform amongst all links.
/// \param params The vector of parameters for the link.
/// \param state : The current state of the state vector.
/// \param user: User defined data.
/// \return Returns the data to be written as output.
typedef int (OutputIntCallback)(unsigned int id, double t, VEC y_i, VEC global_params, VEC params, int state, void* user);

/// 
///
/// \param t The current time.
/// \param y The current state vector at time *t* for the link.
/// \param global_params The vector of parameters uniform amongst all links.
/// \param params The vector of parameters for the link.
/// \param state : The current state of the state vector.
/// \param user: User defined data.
/// \return Returns the data to be written as output.
typedef double (OutputDoubleCallback)(unsigned int id, double t, VEC y_i, VEC global_params, VEC params, int state, void* user);


/// 
///
/// \param t The current time.
/// \param y The current state vector at time *t* for the link.
/// \param global_params The vector of parameters uniform amongst all links.
/// \param params The vector of parameters for the link.
/// \param state : The current state of the state vector.
/// \param user: User defined data.
/// \return Returns the data to be written as output.
typedef float (OutputFloatCallback)(unsigned int id, double t, VEC y_i, VEC global_params, VEC params, int state, void* user);


typedef union OutputCallback {
    OutputIntCallback *out_int;
    OutputDoubleCallback *out_double;
    OutputFloatCallback *out_float;
} OutputCallback;


typedef void (PeakflowOutputCallback)(unsigned int, double, VEC, VEC, VEC, double, unsigned int, void*, char*);

//Function signatures

/// These are the right-hand side functions for the ODEs.
///
/// \param t The current time
/// \param y_i The approximate value of the solution to the ODEs at the current link at time t
/// \param y_p y_p[j] has the approximate values for the immediately upstream link j to the current link at time t
/// \param num_parents The number of upstream links (parents) to link i
/// \param global_params The global parameters
/// \param forcings The rain fall values for link i
/// \param params The parameters for link i
/// \param state The current state of the system
/// \param user A pointer to user specified data
/// \param ans (set by method, assumed that space is allocated): The value returned by the right-hand side function
typedef void (DifferentialFunc)(double t, VEC y_i, VEC* y_p, unsigned short int num_parents, VEC global_params, double* forcings, QVSData*, VEC params, int state, void* user, VEC ans);                    //!< Right-hand side function for ODE                                            
typedef void (AlgebraicFunc)(VEC, VEC, VEC, QVSData*, int, void*, VEC);                                                                  //!< Function for algebraic variables
typedef int (CheckStateFunc)(VEC, VEC, VEC, QVSData*, unsigned int);                                                                     //!< Function to check what "state" the state variables are in (for discontinuities)
typedef void (JacobianFunc)(double, VEC, VEC*, unsigned short int, VEC, double*, VEC, MAT*);                                                    //!< Jacobian of right-hand side function
typedef int (RKSolverFunc)(Link*, GlobalVars*, int*, bool, FILE*, ConnData*, Forcing*, TempStorage*);   //!< RK solver to use
typedef void (CheckConsistencyFunc)(VEC y, VEC, VEC);                                                                                           //!< Function to check state variables

// Custom models signatures
typedef void (SetParamSizesFunc)(GlobalVars*, void*);
typedef void (ConvertFunc)(VEC, unsigned int, void*);
typedef void (RoutinesFunc)(Link*, unsigned int, unsigned int, unsigned short int, void*);
typedef void (PrecalculationsFunc)(Link*, VEC, VEC, unsigned int, unsigned int, unsigned short int, unsigned int, void*);
typedef int (InitializeEqsFunc)(VEC, VEC, QVSData*, unsigned short int, VEC, unsigned int, unsigned int, unsigned int, void*, void*);
typedef int* (PartitionFunc)(Link*, unsigned int, Link**, unsigned int, unsigned int**, unsigned int*, TransData*, short int*);

//Constructor / Destructor related routings

/// This routine initializes an instance of an AsynchSolver. Multiple solvers could be created if multiple
/// problems are to be solved. An instance of an AsynchSolver contains all the relevant data structures and
/// information to solve the equations for an underlying model. An AsynchSolver object should be destroyed
/// with a call to *Asynch_Free*.
///
/// \pre MPI should be properly initialized with MPI_Init before calling Asynch_Init.
/// \param comm  The MPI communicator to use with this solver object.
/// \return A pointer to an AsynchSolver object.
/// \see Asynch_Free
AsynchSolver* Asynch_Init(MPI_Comm comm);

/// This routine deallocates the memory occupied by an AsynchSolver object created with a call to Asynch_Init.
///
/// \param asynch  A pointer to a AsynchSolver object to free.
/// \see Asynch_Init
void Asynch_Free(AsynchSolver* asynch);

//Customization related routings

int Asynch_Custom_Model(
    AsynchSolver* asynch,
    SetParamSizesFunc *set_param_sizes,
    ConvertFunc *convert,
    RoutinesFunc *routines,
    PrecalculationsFunc *precalculations,
    InitializeEqsFunc *initialize_eqs);

int Asynch_Custom_Partitioning(AsynchSolver* asynch, PartitionFunc *partition);

//Routines to intialize network and model

/// This routine opens and processes a global file. It reads all specified database connection files,
/// but does not process any other input file. An error in this routine is considered fatal, and results
/// in a call to the routine *MPI_Abort* on the communicator used to create asynch.
///
/// \pre asynch should be initialized with Asynch_Init before calling Asynch_Parse_GBL.
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param filename  Path of a global file.
void Asynch_Parse_GBL(AsynchSolver* asynch, char* filename);

/// This routine processes topology inputs for the AsynchSolver object as set in the global file read
/// by Asynch_Parse_GBL. This initializes each Link object and sets their parent and child information.
/// Generally, this is the first initialization routine to call after parsing a global.
///
/// \pre Asynch_Parse_GBL should be called before Asynch_Load_Network.
/// \param asynch A pointer to a AsynchSolver object to use.
void Asynch_Load_Network(AsynchSolver* asynch);

/// This routine assigns the Links of the *asynchsolver* object to the MPI processes.
///
/// \pre This routine must be called after *Asynch_Load_Network*.
/// \param asynch A pointer to a AsynchSolver object to use.
void Asynch_Partition_Network(AsynchSolver* asynch);

/// This routine processes parameter inputs for the AsynchSolver object as set in the global file read by
/// *Asynch_Parse_GBL*. Setting *load_all* to true causes every MPI process to store the parameters at every Link.
/// 
/// \pre This routine can be called before *Asynch_Partition_Network* only if *load_all* is set to true.
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param load_all  Setting load_all to false causes the MPI processes to only store parameters for Links assigned to them.
void Asynch_Load_Network_Parameters(AsynchSolver* asynch, short int load_all);

/// This routine processes the dam inputs for the AsynchSolver object as set in the global file read by
// *Asynch_Parse_GBL*.
/// 
/// \pre This routine should be called after *Asynch_Partition_Network* and *Asynch_Load_Network_Parameters* have been called.
/// \param asynch A pointer to a AsynchSolver object to use.
void Asynch_Load_Dams(AsynchSolver* asynch);

/// This routine processes numerical solver inputs for the AsynchSolver object as set in the global file read by *Asynch_Parse_GBL*.
/// 
/// \pre This routine should be called after *Asynch_Partition_Network* has been called.
/// \param asynch A pointer to a AsynchSolver object to use.
void Asynch_Load_Numerical_Error_Data(AsynchSolver* asynch);

/// This routine sets the model specific routines for each link for the AsynchSolver object as set in the global file read by *Asynch_Parse_GBL*.
///
/// \pre This routine should be called after *Asynch_Partition_Network* and *Asynch_Load_Network_Parameters* have been called.
/// \param asynch A pointer to a AsynchSolver object to use.
void Asynch_Initialize_Model(AsynchSolver* asynch);

/// This routine processes the initial condition inputs for the AsynchSolver object as set in the global file read by *Asynch_Parse_GBL*.
///
/// \pre This routine should be called after *Asynch_Partition_Network* and *Asynch_Initialize_Model* have been called.
/// \param asynch A pointer to a AsynchSolver object to use.
void Asynch_Load_Initial_Conditions(AsynchSolver* asynch);

/// This routine processes the forcing inputs for the AsynchSolver object as set in the global file read by *Asynch_Parse_GBL*.
///
/// \pre This routine should be called after *Asynch_Partition_Parameters* has been called.
/// \param asynch A pointer to a AsynchSolver object to use.
void Asynch_Load_Forcings(AsynchSolver* asynch);

/// This routine processes save list inputs for the AsynchSolver object as set in the global file read by *Asynch_Parse_GBL*.
///
/// \pre This routine should be called after *Asynch_Partition_Network* has been called.
/// \param asynch A pointer to a AsynchSolver object to use.
void Asynch_Load_Save_Lists(AsynchSolver* asynch);

/// This routine checks that all inputs are loaded for the AsynchSolver object. Some small final initializations are also performed.
///
/// \pre This routine should be called as the last initialization routine.
/// \param asynch A pointer to a AsynchSolver object to use.
void Asynch_Finalize_Network(AsynchSolver* asynch);

// Integration related routines

/// This routine processes calculates appropriate step sizes for the integrators at each link in the AsynchSolver object.
///
/// \pre This routine must be called before a call to *Asynch_Advance*, and after all initializations are performed.
/// \param asynch A pointer to a AsynchSolver object to use.
void Asynch_Calculate_Step_Sizes(AsynchSolver* asynch);

/// This routine advances the numerical solver up to the time set in *maxtime*. See Section [sec: model type and maxtime].
/// Calculations to solve the model differential and algebraic equations are performed, using forcing data as needed.
///
/// \pre This routine should be called after *Asynch_Partition_Network* has been called.
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param print_flag  If false, no time series information is produced. Otherwise, time series information is produced.
void Asynch_Advance(AsynchSolver* asynch, bool print_flag);

//Snapshot

/// This routine processes calculates appropriate step sizes for the integrators at each link in the AsynchSolver object.
///
/// \pre This routine must be called before a call to *Asynch_Advance*, and after all initializations are performed.
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param prefix String to prepend to the snapshot output filename.
/// \return An error code. Returns 0 if a snapshot was made, 1 if an error was encountered, -1 if no snapshot is made.
int Asynch_Take_System_Snapshot(AsynchSolver* asynch, char* prefix);

//Forcing related routines

/// This routine activates a forcing for use.This means when a call to the routine *Asynch_Advance* is made,
/// the forcing with index *idx* will be applied for calculations.By default, all forcings set in a global file
/// are initially active.This routine has the opposite effect of *Asynch_Deactivate_Forcing*.
///
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param idx The index of the forcing to activate.
/// \return  Returns 0 if the forcing was activated, 1 if an error occurred
int Asynch_Activate_Forcing(AsynchSolver* asynch, unsigned int idx);

/// This routine deactivates a forcing for use. This means when a call to the routine *Asynch_Advance* is made,
/// all values for the forcing with index *idx* will be taken as 0. By default, all forcings set in a global file
/// are initially active. This routine has the opposite effect of *Asynch_Activate_Forcing*.
///
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param idx The index of the forcing to activate.
/// \return  Returns 0 if the forcing was activated, 1 if an error occurred
int Asynch_Deactivate_Forcing(AsynchSolver* asynch, unsigned int idx);

//Data file routines

/// This routine prepares the output sources for time series data. Preparation includes creating files or database tables.
///
/// \pre This routine must be called before any time series data can be produced.
/// \param asynch A pointer to a AsynchSolver object to use.
void Asynch_Prepare_Output(AsynchSolver* asynch);

/// This routine prepares the output sources for the peakflow data. Preparation includes creating files or database tables.
///
/// \pre This routine must be called before any peakflow data can be produced.
/// \param asynch A pointer to a AsynchSolver object to use.
void Asynch_Prepare_Peakflow_Output(AsynchSolver* asynch);

/// This routine takes all data written to temporary files and moves them to a final output destination for time series data.
/// If the output format is data file or CSV file, then the string *additional_out* is appended to the filename.
/// If *additional_out* is *NULL*, then no appending to the filename occurs.
///
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param file_suffix String appended to the end of any output filename.
///\ return Returns 0 if output is written, -1 means there is no data to output.
int Asynch_Create_Output(AsynchSolver* asynch, char* file_suffix);

/// This routine takes all calculated peakflow data and writes them to a final output destination.
///
/// \param asynch A pointer to a AsynchSolver object to use.
///\ return Returns 0 if output is written, -1 means there is no data to output.
int Asynch_Create_Peakflows_Output(AsynchSolver* asynch);

/// This routine writes the current state of each link to temporary files for every link where a time series output has
/// been specified. Normally, calling *Asynch_Advance* with the *print_flag* set is enough to write output time series.
/// However, advancing an AsynchSolver object without the *print_flag* set or calls that modify how data is outputted to
/// the temporary files will cause data to not be written for some times. This routine can be called to commit missing data.
///
/// \param asynch A pointer to a AsynchSolver object to use.
/// \return Returns 0 if the step is written, 1 if no temporary file is available.
int Asynch_Write_Current_Step(AsynchSolver* asynch);

/// This routine determines whether the time series with name *name* is being used for the current simulations.
/// If *name* is requested in the global file and has been set, the routine returns 1. If the *name* is requested in
/// the global file, but has not yet been set, this routine returns 0. If *name* was not requested in the global file,
/// the return value is -1. An output can be set with a call to *Asynch_Set_Output*.
/// See Section [sec: asynch\ :sub:`s`\ et\ :sub:`o`\ utput] for details about setting outputs.
///
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param name The name of the time series output to check.
/// \return Returns 1 if the output *name* is set, 0 if it is not set, and -1 if it is not present.
int Asynch_Check_Output(AsynchSolver* asynch, char* name);

/// This routine determines whether the time series with name *name* is being used for the current simulations.
/// If *name* is requested in the global file and has been set, the routine returns 1. If the *name* is requested in
/// the global file, but has not yet been set, this routine returns 0. If *name* was not requested in the global file,
/// the return value is -1. An output can be set with a call to *Asynch_Set_Output*.
/// See Section [sec: asynch\ :sub:`s`\ et\ :sub:`o`\ utput] for details about setting outputs.
///
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param name The name of the time series output to check.
/// \return Returns 1 if the output *name* is set, 0 if it is not set, and -1 if it is not present.
int Asynch_Check_Peakflow_Output(AsynchSolver* asynch, char* name);

//Routines for custom output

/// This routine sets a custom output time series. A function is required that will be called every time
/// output data is to be written for a link. The function set in this routine should have the
/// specification
/// Asynch_Set_Output_Int must also be given an array *used_states*. This array contains the indices of
/// the states in the state vectors which are needed by the user specified routine *callback*. All states
/// listed in this array are guaranteed to be available for the routine *callback*.
///
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param name Name of the custom time series.
/// \param callback Routine to call when writing data for the custom output.
/// \param used_states Array of the all indices in the state vector used by *callback*.
/// \param num_states Numver of indiciced in used_states.
/// \see Asynch_Set_Output_Double
int Asynch_Set_Output_Int(AsynchSolver* asynch, char* name, OutputIntCallback* callback, unsigned int* used_states, unsigned int num_states);

/// This routine sets a custom output time series. A function is required that will be called every time
/// output data is to be written for a link. The function set in this routine should have the
/// specification
/// Asynch_Set_Output_Double must also be given an array *used_states*. This array contains the indices of
/// the states in the state vectors which are needed by the user specified routine *callback*. All states
/// listed in this array are guaranteed to be available for the routine *callback*.
///
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param name Name of the custom time series.
/// \param callback Routine to call when writing data for the custom output.
/// \param  used_states Array of the all indices in the state vector used by *callback*.
/// \param num_states Numver of indiciced in used_states.
/// \see Asynch_Set_Output_Int
int Asynch_Set_Output_Double(AsynchSolver* asynch, char* name, OutputDoubleCallback* callback, unsigned int* used_states, unsigned int num_states);

/// This routine sets a custom output time series. A function is required that will be called every time
/// output data is to be written for a link. The function set in this routine should have the
/// specification
/// Asynch_Set_Output_Double must also be given an array *used_states*. This array contains the indices of
/// the states in the state vectors which are needed by the user specified routine *callback*. All states
/// listed in this array are guaranteed to be available for the routine *callback*.
///
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param name Name of the custom time series.
/// \param callback Routine to call when writing data for the custom output.
/// \param  used_states Array of the all indices in the state vector used by *callback*.
/// \param num_states Numver of indiciced in used_states.
/// \see Asynch_Set_Output_Int
int Asynch_Set_Output_Float(AsynchSolver* asynch, char* name, OutputFloatCallback* callback, unsigned int* used_states, unsigned int num_states);

/// This routine sets the filename of the output peakflow file. If a database connection is used, then
/// no changes are made.
///
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param name The name of the time series output to check.
/// \param callback Pointer to the routine to call when writing peakflow data for the custom output.
/// \return 1 if the output is set, 0 if an error occurred.
int Asynch_Set_Peakflow_Output(AsynchSolver* asynch, char* name, PeakflowOutputCallback * callback);


int Asynch_Create_OutputUser_Data(AsynchSolver* asynch, unsigned int data_size);
int Asynch_Free_OutputUser_Data(AsynchSolver* asynch);
void Asynch_Copy_Local_OutputUser_Data(AsynchSolver* asynch, unsigned int location, void* source, unsigned int size);
void Asynch_Set_Size_Local_OutputUser_Data(AsynchSolver* asynch, unsigned int location, unsigned int size);


// Temp files routines

/// This routine prepares the temporary files for time series data. Preparation includes creating files or database tables.
///
/// \pre This routine must be called before any time series data can be calculated. A call to *Asynch_Advance* with the
/// *print_flag* set before preparing temporary files will create an error.
/// \param asynch A pointer to a AsynchSolver object to use.
void Asynch_Prepare_Temp_Files(AsynchSolver* asynch);

/// This routine deletes the tempfiles for *asynch*. This is useful for cleaning up temporary files at the end of a simulation,
/// or for deleting the files if they must be reconstructed.
///
/// \param asynch A pointer to a AsynchSolver object to use.
/// \return Returns 0 if the tempfiles were deleted, 1 if no tempfiles exist, and 2 if an error occurred.
int Asynch_Delete_Temporary_Files(AsynchSolver* asynch);

/// This routine moves the temporary file pointer to a previous value. All future values are deleted.
/// The point where the pointer is moved is determined by *set_value*, which is in the output series determined by *output_idx*.
/// The local time where the next data will be written for each link is *set_time*.
///
/// \pre A warning occurs if *set_value* is not found for a link, and that link's tempfile pointer is not changed.
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param set_time The local time corresponding the step where the temp files will be set.
/// \param set_value The value to which the tempfiles will be set.
/// \param output_idx The index of the series output of which *set_value* corresponds.
/// \return Returns 0 if the temporary files are set, 1 if a warning occurred, 2 if there was an error setting the files.
int Asynch_Set_Temp_Files(AsynchSolver* asynch, double set_time, void* set_value, unsigned int output_idx);

/// This routine moves the tempfile pointer to the beginning of the file for each link. All future values are deleted.
/// The local time where the next data will be written for each link is *set_time*.
///
/// \pre A warning occurs if *set_value* is not found for a link, and that link's tempfile pointer is not changed.
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param set_time The time from the beginning of the simulation corresponding the step where the temp files will be set.
/// \return Returns 0 if the temporary files are set, 1 if a warning occurred, 2 if there was an error setting the files.
int Asynch_Reset_Temp_Files(AsynchSolver* asynch, double set_time);


//Snapshot

/// This routine creates snapshot output data to either a recovery file or a database table. The current
/// value of every state at every link is outputted. If a recovery file was the format specified in the
/// global file, then the string *name* is appended to the end of the recovery filename. This appending
/// does not occur if *name* is *NULL* or if a database table is the selected format.
///
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param prefix String to prepend to the snapshot output filename.
/// \return Returns 0 if a snapshot was made, 1 if an error was encountered, -1 if no snapshot is made.
int Asynch_Take_System_Snapshot(AsynchSolver* asynch, char* prefix);

/// This routine gets the filename of the output snapshot file. If a database connection is used, then
/// the contents of *snapshotname* is not modified. The Python interface routine returns the string with
/// the filename instead of taking it as an argument.
///
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param filename Name of the snapshot output filename returned.
/// \return 0 if the filename was set successfully. 1 otherwise.
int Asynch_Get_Snapshot_Output_Name(AsynchSolver* asynch, char* filename);

/// This routine sets the filename of the output snapshot file. If a database connection is used, then
/// no changes are made.
///
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param filename Name to set the snapshot output filename.
/// \return 0 if the filename was set successfully. 1 otherwise.
int Asynch_Set_Snapshot_Output_Name(AsynchSolver* asynch, char* filename);

//Set and get routines

/// This routine sets the model type.
///
/// \param asynch A pointer to a AsynchSolver object to use.
/// \return The model type
unsigned short Asynch_Get_Model_Type(AsynchSolver* asynch);

/// This routine returns the model type.
///
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param type    The model type.
void Asynch_Set_Model_Type(AsynchSolver* asynch, unsigned short type);

/// This routine returns the value of *duration*, as defined in Section[sec:model type and duration].
///
/// \param asynch A pointer to a AsynchSolver object to use.
/// \return The current duration of the simulation time of the AsynchSolver object.
double Asynch_Get_Total_Simulation_Duration(AsynchSolver* asynch);

/// This routine returns the value of *duration*, as defined in Section[sec:model type and duration].
///
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param duration The value to set the maximum simulation time.
void Asynch_Set_Total_Simulation_Duration(AsynchSolver* asynch, double duration);

/// This routine returns the total number of links in the network of *asynch*.
///
/// \param asynch A pointer to a AsynchSolver object to use.
/// \return Number of links in the network
unsigned short Asynch_Get_Num_Links(AsynchSolver* asynch);

Link* Asynch_Get_Links(AsynchSolver* asynch);

/// This routine returns the total number of links assigned to the current MPI process in the network of *asynch*.
/// This number represents the total number of links whose differential and algebraic equations are solved by the
/// current process.
///
/// \param asynch A pointer to a AsynchSolver object to use.
/// \return Number of links assigned to the current MPI process
unsigned short Asynch_Get_Num_Links_Proc(AsynchSolver* asynch);

Link* Asynch_Get_Links_Proc(AsynchSolver* asynch);

/// This routine sets a new database for an input or output. If information for a database has already been set,
/// it is released, and the new connection information is set.
/// Database information includes hostname, username, password, etc. This is the same information that is available
/// in the header of database connection files (see Section [sec: database connection files]).
/// Several database connections exist for every AsynchSolver object. *conn_idx* can take the values:
///  - ASYNCH_DB_LOC_TOPO
///  - ASYNCH_DB_LOC_PARAMS
///  - ASYNCH_DB_LOC_INIT
///  - ASYNCH_DB_LOC_RSV
///  - ASYNCH_DB_LOC_HYDROSAVE
///  - ASYNCH_DB_LOC_PEAKSAVE
///  - ASYNCH_DB_LOC_HYDRO_OUTPUT
///  - ASYNCH_DB_LOC_PEAK_OUTPUT
///  - ASYNCH_DB_LOC_SNAPSHOT_OUTPUT
///  - ASYNCH_DB_LOC_FORCING_START
///
/// The last value for *conn_idx* is the database connection for the first forcing specified in the global file.
/// Database connections for other forcings can be access sequentially. For example, to access the forcing with index 2 
/// (the third forcing in a global file), set *conn_idx* as `ASYNCH_DB_LOC_FORCING_START + 2`
///
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param conn_string A connection string
/// \param conn_idx Connection index
void Asynch_Set_Database_Connection(AsynchSolver* asynch, const char* conn_string, unsigned int conn_idx);


/// This routine returns the first timestamp for a forcing.
///
/// \pre This can only be used if the forcing with index *forcing_idx* is using a format of binary files, gz binary files,
/// or database table. See Section [sec: forcing inputs] for adescription of these formats.
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param forcing_idx An index of a forcing.
/// \return The timestamp of the last timestamp for the forcing with index *forcing_idx*.
unsigned int Asynch_Get_First_Forcing_Timestamp(AsynchSolver* asynch, unsigned int forcing_idx);

/// This routine sets the first timestamp for a forcing.
///
/// \pre This can only be used if the forcing with index *forcing_idx* is using a format of binary files, gz binary files,
/// or database table. See Section [sec: forcing inputs] for adescription of these formats.
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param unix_time The value to set the first forcing timestamp.
/// \param forcing_idx An index of a forcing.
/// \return The timestamp of the last timestamp for the forcing with index *forcing_idx*.
void Asynch_Set_First_Forcing_Timestamp(AsynchSolver* asynch, unsigned int unix_time, unsigned int forcing_idx);

/// This routine returns the last timestamp for a forcing.
///
/// \pre This can only be used if the forcing with index *forcing_idx* is using a format of binary files, gz binary files,
/// or database table. See Section [sec: forcing inputs] for adescription of these formats.
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param forcing_idx An index of a forcing.
/// \return The timestamp of the last timestamp for the forcing with index *forcing_idx*.
unsigned int Asynch_Get_Last_Forcing_Timestamp(AsynchSolver* asynch, unsigned int forcing_idx);

/// This routine sets the last timestamp for a forcing.
///
/// \pre This can only be used if the forcing with index *forcing_idx* is using a format of binary files, gz binary files,
/// or database table. See Section [sec: forcing inputs] for adescription of these formats.
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param unix_time The value to set the last forcing timestamp.
/// \param forcing_idx An index of a forcing.
/// \return The timestamp of the last timestamp for the forcing with index *forcing_idx*.
void Asynch_Set_Last_Forcing_Timestamp(AsynchSolver* asynch, unsigned int unix_time, unsigned int forcing_idx);

/// This routine sets the start time used for a forcing. This value is used for converting between timestamps in a database
/// table and the local time of the solvers. This can only be used if the forcing with index *forcing_idx* is using a format
/// of database table. See Section [sec: forcing inputs] for a description of this format.
///
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param unix_time The value to set the start timestamp.
/// \param forcing_idx An index of a forcing.
void Asynch_Set_Forcing_DB_Starttime(AsynchSolver* asynch, unsigned int unix_time, unsigned int forcing_idx);


// State variable

/// This routine sets a file for reading initial value data. The *init_flag* is set based upon the extension of *filename*.
/// The initial data is **NOT** read while executing this routine. A call to *Asynch_Load_System* is needed to set the filename.
/// This routine cannot be used to set the format to a database connection.
///
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param filename Filename to use for initial value data.
void Asynch_Set_Init_File(AsynchSolver* asynch, char* filename);

/// This routine sets the current time of each link to *t_0* and the current state to the vector in the array *values*.
/// The vectors are assumed to be in the same order as the links provided by the topology source specified in the global file
/// (see Section [sec: topology]).
///
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param unix_time The timestamp to set at each link.
/// \param states The vector of the current state of each link.
void Asynch_Set_System_State(AsynchSolver* asynch, double unix_time, VEC* states);


/// This routine clears all peakflow data for each link. The time to peak is set to the current local
/// time of *asynch*, and the value of the peakflow is set to the current state vector of the link.
///
/// \param asynch A pointer to a AsynchSolver object to use.
void Asynch_Reset_Peakflow_Data(AsynchSolver* asynch);

int Asynch_Set_Forcing_State(AsynchSolver* asynch, unsigned int idx, double t_0, unsigned int first_file, unsigned int last_file);

/// This routine gets the filename of the output peakflow file. If a database connection is used, then
/// the contents of *peakflowname* is not modified. The Python interface routine returns the string with
/// the filename instead of taking it as an argument.
///
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param peakflowname Name of the peakflow output filename returned.
/// \return 0 if the filename was set successfully. 1 otherwise.
int Asynch_Get_Peakflow_Output_Name(AsynchSolver* asynch, char* peakflowname);

/// This routine sets the filename of the output peakflow file. If a database connection is used, then
/// no changes are made.
///
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param peakflowname Name to set the peakflow output filename.
/// \return 0 if the filename was set successfully. 1 otherwise.
int Asynch_Set_Peakflow_Output_Name(AsynchSolver* asynch, char* peakflowname);

unsigned int Asynch_Get_Local_LinkID(AsynchSolver* asynch, unsigned int location);

int Asynch_Set_Init_Timestamp(AsynchSolver* asynch, unsigned int unix_time);
unsigned int Asynch_Get_Init_Timestamp(AsynchSolver* asynch);

/// This routine sets the filename of the output snapshot file. If a database connection is used, then
/// no changes are made.
///
/// \param asynch A pointer to a AsynchSolver object to use.
/// \return Returns 1 if the model has reservoir forcing. 0 otherwise.
int Asynch_Get_Reservoir_Forcing(AsynchSolver* asynch);

///This routine returns the number of parameters for the differential equations that are available to
///every link. This number is either specified in a global file, or with a call to the routine
///*Asynch_Set_Global_Parameters*.
/// 
/// \param asynch A pointer to a AsynchSolver object to use.
/// \return The number of global parameters, typically specified in a global file.
unsigned int Asynch_Get_Size_Global_Parameters(AsynchSolver* asynch);

/// This routine gets the values of the global parameters. For the C function, the array *params*
/// should contain enough memory to hold the global parameters. The number of global parameters can be
/// determined with a call to *Asynch_Get_Size_Global_Parameters*. The obtained values are a copy of
/// Asynch's internal values of the parameters, so modifying them will have no effect on the solver
/// results.
/// 
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param params Vector of global parameter values to retrieve.
/// \return Vector of the global parameters.
void Asynch_Get_Global_Parameters(AsynchSolver* asynch, VEC params);

/// This routine sets the values of the global parameters. The number of global parameters may be
/// different from the number previously stored by Asynch. However, if the new global parameters should
/// be consistent with the model equations at each link (for example, if the model at each hillslope
/// expects 4 global parameters, using this routine to change the number of global parameters to 2 WILL
/// cause an error). The values are copied from the array (or list) into Asynch's internal memory. So
/// changing their values after this function call will NOT change the solver results.
///
/// \param asynch A pointer to a AsynchSolver object to use.
/// \param params Vector of global parameter values to set.
/// \param num_params The number of global parameters in params.
/// \return Returns 1 if an error occurred, 0 otherwise.
int Asynch_Set_Global_Parameters(AsynchSolver* asynch, VEC params, unsigned int num_params);



#endif

