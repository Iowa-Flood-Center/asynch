#ifndef STRUCTS_H
#define STRUCTS_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <stdio.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#include "libpq-fwd.h"

#include "asynch_interface.h"

//Note: For parser, uncomment out formula and equation data here.
//Also go through riversys.c for equation (3 instances).
//Also make sure headers are included everywhere. Uncomment the functions in problems.c and .h.
//And system.c (1 block).
//Did I forget anything?

//Constants
#define ASYNCH_MAX_DB_CONNECTIONS 20

#define ASYNCH_DB_LOC_TOPO 0
#define ASYNCH_DB_LOC_PARAMS 1
#define ASYNCH_DB_LOC_INIT 2
#define ASYNCH_DB_LOC_QVS 3
#define ASYNCH_DB_LOC_RSV 4
#define ASYNCH_DB_LOC_HYDROSAVE 5
#define ASYNCH_DB_LOC_PEAKSAVE 6
#define ASYNCH_DB_LOC_HYDRO_OUTPUT 7
#define ASYNCH_DB_LOC_PEAK_OUTPUT 8
#define ASYNCH_DB_LOC_SNAPSHOT_OUTPUT 9
#define ASYNCH_DB_LOC_FORCING_START 10

#define ASYNCH_MAX_QUERY_LENGTH 2048
#define ASYNCH_MAX_CONNSTRING_LENGTH 1024
#define ASYNCH_MAX_QUERIES 5

#define ASYNCH_MAX_NUM_FORCINGS 12
#define ASYNCH_MAX_TIMESTAMP_LENGTH 12
#define ASYNCH_MAX_PATH_LENGTH 1024

#define ASYNCH_MAX_LINE_LENGTH 1024
#define ASYNCH_MAX_SYMBOL_LENGTH 64
#define ASYNCH_MAX_QUERY_LENGTH 2048

#define ASYNCH_MAX_DIM 256              //!< Maximum number of Degree of Freedom

#define ASYNCH_LINK_MAX_PARENTS 8


//struct Forcing;
//struct Link;
//struct GlobalVars;
//struct TransData;

// Forward definitions
typedef struct ErrorData ErrorData;
typedef struct GlobalVars GlobalVars;
typedef struct Link Link;
typedef struct RKMethod RKMethod;
typedef struct TransData TransData;
typedef struct TempStorage TempStorage;
typedef struct ConnData ConnData;
typedef struct Forcing Forcing;
typedef struct AsynchSolver AsynchSolver;

/// Structure to store temporary memory needed for RK solvers.
///
typedef struct TempStorage
{
    //Memory for all Solvers
    VEC sum, temp, temp2, temp3;       //!< Vectors for summations and temp workspace. size = dim of problem at each link.
    VEC** temp_parent_approx;       //!< List of vectors to hold temporary work from parent links. size of each is dim of problem.
    VEC* temp_k;                    //!< List of vectors to hold temporary internal stage values. size of each is num_dense.

    //Memory for Implicit Solvers
    int* ipiv;              //!< Array to hold pivots from LU decomps. length = s*dim.
    VEC rhs;                //!< Holds right hand side of linear systems. size = s*dim.
    //MAT* CoefMat;         //!< Holds coefficient matrix of linear systems. size = s*dim x s*dim.
    MAT JMatrix;            //!< Holds jacobian matrix of the right hand side function of the ode. size = dim x dim.
    VEC* Z_i;               //!< Space for s internal stages. Each had size = dim.
    VEC err;                //!< Space for error approximations. size = dim.
} TempStorage;

/// Holds all information for an RK method.
/// These are intended for dense output methods, but regular RK methods could be stored here as well.
typedef struct RKMethod
{
    MAT A;                  //!< A coefficients
    VEC b;                  //! <b coefficients
    VEC b_theta;            //!< b coefficients evaluated at a value theta in [0,1]
    VEC b_theta_deriv;
    VEC c;                  //!< c coefficients
    void(*dense_b)(double, VEC);        //!< Function to evaluate b at a value theta in [0,1]
    void(*dense_bderiv)(double, VEC);   //!< Derivative of b polynomials
    VEC e;                              //!< e vector for error coefficients
    VEC d;                              //!< d vector for dense error coefficients
    unsigned short int s;               //!< Number of stages
    unsigned short int unique_c;        //!< Number of unique values in c
    unsigned short int e_order;         //!< Error order + 1
    unsigned short int d_order;         //!< Dense error order + 1
    double e_order_ratio;               //!< e_order / Error order
    double d_order_ratio;               //!< d_order / Dense error order
    VEC w;                              //!< Weights for lagrange polynomial
    unsigned short int exp_imp;         //!< 0 if method is explicit, 1 if implicit
    unsigned short int localorder;      //!< Local order of the method
} RKMethod;

/// Holds the error estimation information for a link.
/// See Hairer, E. and Norsett, S.P. and Wanner, G., Solving Ordinary Differential Equations I, Nonstiff Problems.
typedef struct ErrorData
{
    double facmax;          //!< Parameter for error estimation
    double facmin;          //!< Parameter for error estimation
    double fac;             //!< Parameter for error estimation
    VEC abstol;             //!< Absolute tolerance
    VEC reltol;             //!< Relative tolerance
    VEC abstol_dense;       //!< Absolute tolerance for dense output
    VEC reltol_dense;       //!< Relative tolerance for dense output
} ErrorData;

/// Node for a linked list of the numerical solution for a link. Each node holds the numerical solution at time t.
///
typedef struct RKSolutionNode
{
    VEC* k;                 //!< Array of all k values at time t
    VEC y_approx;           //!< Approximate solution at time t
    double t;               //!< The time to which the data in this node corresponds
    struct RKSolutionNode* next;    //!< Next node in the linked list
    struct RKSolutionNode* prev;    //!< Previous node in the linked list
    int state;              //!< State of the solution
} RKSolutionNode;

/// Linked list for the numerical solution of a link.
///
typedef struct RKSolutionList
{
    RKSolutionNode* list_data;  //!< A pointer to the nodes in this list. Used for allocation/deallocation.
    RKSolutionNode* head;       //!< The beginning of the list. This node has the small t value.
    RKSolutionNode* tail;       //!< The end of the list. This node has the largest t value.
    //unsigned int dim;         //!< The dimension of the problem. This is the size of all y_approx vectors in each node.
    unsigned short int s;       //!< The number of stages in the RK method used to create these approximations.
} RKSolutionList;

/// Structure to contain the forcing data of a link.
///
typedef struct ForcingData
{
    double** data;          //!< 2D array with 2 columns. First column is time the rainfall changes to the rate in the second column
    unsigned int nrows;     //!< Number of rows in rainfall
} ForcingData;

/// Structure to contain the discharge vs storage data of a link.
///
typedef struct QVSData
{
    double** points;            //!< 2D array with 2 columns. First column is time the rainfall changes to the rate in the second column
    double* points_array;       //!< 1D equivalent of points
    unsigned int n_values;      //!< Number of rows in points (2*n_values is the number entries in points_array)
} QVSData;

/*
//Structure to contain the diff eqs for the standardized template models
typedef struct Formula
{
    muParserHandle_t parser;
    VEC* variable_values;
} Formula;
*/

/// Structure to hold information about an PostgreSQL database
///
typedef struct ConnData
{
    PGconn* conn;                                   //!< Connection to a database
    char connectinfo[ASYNCH_MAX_CONNSTRING_LENGTH]; //!< Connection information for a database
    char query[ASYNCH_MAX_QUERY_LENGTH];            //!< Buffer space for making queries
    //char* submission;         //!< Buffer space for making large submissions to the database
    char* queries[ASYNCH_MAX_QUERIES];
    unsigned int num_queries;
    //unsigned int submission_size;        //Size of submission buffer
    //unsigned int submission_content;    //Amount of characters (bytes) currently stored in submission
    unsigned int time_offset;       //!< Added to the integration time to get the actual unix time

} ConnData;

typedef struct OutputFunc
{
    //Temporary Calculations
    FILE* (*PrepareTempOutput)(struct Link*, unsigned int, int*, struct GlobalVars*, unsigned int*, unsigned int, unsigned int, char*, unsigned int**);

    //Prepare Final Output
    void (*PrepareOutput)(struct GlobalVars*, struct ConnData*);
    int (*PreparePeakflowOutput)(struct GlobalVars*, unsigned int);

    //Create Final Output
    int (*CreateOutput)(struct Link*, struct GlobalVars*, unsigned int, unsigned int*, unsigned int, unsigned int, unsigned int**, int*, char*, char*, struct ConnData*, FILE**);
    int (*CreatePeakflowOutput)(struct Link*, struct GlobalVars*, unsigned int, int*, unsigned int*, unsigned int, unsigned int**, struct ConnData*);

    //Create Snapshot
    int (*CreateSnapShot)(struct Link*, unsigned int, int*, struct GlobalVars*, char*, struct ConnData*);
} OutputFunc;


/// Structure to contain all data that is global to the river system.
///
typedef struct GlobalVars
{
    unsigned short int type;        //!< Index for the model used

    double maxtime;                 //!< Integrate up to this time (duration)
    double t_0;                     //!< Initial time to start integration

    unsigned int start_time;        //!< Unix start time
    unsigned int end_time;          //!< Unix end time

    unsigned short int method;      //!< RK method to use (if it is the same for all links)
    unsigned short int max_s;       //!< The largest number of internal stages of any RK method used    !!!! Is this needed? !!!!
    unsigned short int max_parents; //!< The largest number of parents any link has
    int iter_limit;                 //!< If a link has >= iter_limit of steps stored, no new computations occur
    int max_transfer_steps;         //!< Maximum number of steps to communicate at once between processes
    //unsigned int dim;             //!< The dimension of the ODE to solve at each link
    //unsigned int problem_dim;     //!< Same as dim when not using data assimilation. Otherwise, it's the model dimension
    //unsigned int increment;       //!< Number of rainfall files to store in memory at a time
    //unsigned int first_file;      //!< The index of the first rainfall file
    //unsigned int last_file;       //!< The index of the last rainfall file
    unsigned int discont_size;      //!< Size of discont, discont_send, discont_order_send at each link
    //double file_time;             //!< The time duration that a rainfall file lasts
    unsigned int max_localorder;    //!< Max local order of implemented numerical methods
    //unsigned int diff_start;      //!< Starting index of differential variables in solution vectors
    //unsigned int no_ini_start;    //!< Starting index of differential variables not read from disk
    unsigned short int uses_dam;    //!< 1 if this type can use dams, 0 else
    //unsigned short int rain_flag; //!< 0 for no rain, 1 for .str file, 2 for binary files, 3 for SQL database, 4 for uniform rain
    VEC global_params;              //!< List of global parameters
    unsigned int params_size;       //!< Size of params at each link without a dam
    //unsigned int iparams_size;    //!< Size of iparams at each link
    unsigned int dam_params_size;   //!< Size of params at each link with a dam
    unsigned int disk_params;       //!< Number of parameters to read from disk
    unsigned int area_idx;          //!< Index of upstream area (A_i) in params
    unsigned int areah_idx;         //!< Index of hillslope area (A_h) in params
    char* rain_filename;
    char* init_filename;
    char* rvr_filename;
    char* prm_filename;
    unsigned short int init_flag;   //!< 0 if reading .ini file, 1 if reading .uini file, 2 if reading .rec file
    unsigned short int rvr_flag;    //!< 0 if reading .rvr file, 1 if using database
    unsigned short int prm_flag;    //!< 0 if reading .prm file, 1 if using database
    //unsigned short int output_flag;   //0 for matlab (.dat), 1 for .csv
    //char* results_folder;
    //char* temp_folder;
    char* temp_filename;
    char* dam_filename;
    double print_time;              //!< Each link will write state every print_time minutes. -1 uses a formula.
    unsigned short int print_par_flag;  //!< 0 to use specified name for output files, 1 to add parameters
    unsigned short int dam_flag;        //!< 0 if not using .dam file, 1 if using
    unsigned short int hydrosave_flag;  //!< 0 if not saving hydrographs, 1 if saving
    unsigned short int peaksave_flag;   //!< 0 if not saving peak flows, 1 if saving
    char* hydrosave_filename;
    char* peaksave_filename;
    char* peakfilename;             //!< Filename for .pea file
    //char* identifier;
    unsigned int max_dim;
    unsigned int outletlink;        //!< For database: holds the link id of the outlet. Use 0 if reading entire database.
    //unsigned int num_dense;       //!< Number of states where dense output is calculated
    //unsigned int* dense_indices;  //!< List of indices in solution where dense output is needed
    unsigned int string_size;       //!< Size of filename buffers
    unsigned int query_size;        //!< Size of database query buffers
    //unsigned int raindb_start_time;       //!< Starting time for database if using rainfall data
    //unsigned short int assim_flag;        //!< 0 if not using data assimilation, 1 else
    short int rkd_flag;
    unsigned short int convertarea_flag;    //!< 1 if hillslope and upstream areas are converted from km^2 to m^2, 0 if not
    double discont_tol;                     //!< The error tolerance to use for locating discontinuities.
    //unsigned short int template_flag;

    unsigned int min_error_tolerances;      //!< The minimum number of error tolerances needed at every link. Used for uniform error tolerances.
    unsigned int num_forcings;

    short unsigned int hydros_loc_flag;
    short unsigned int peaks_loc_flag;
    short unsigned int dump_loc_flag;
    double dump_time;                   //!< Each link states will dump every dump_time minutes.
    short unsigned int res_flag;
    char* hydros_loc_filename;
    char* peaks_loc_filename;
    char* dump_loc_filename;
    char* rsv_filename;
    unsigned int init_timestamp;
    short int res_forcing_idx;


    //Outputs
    unsigned int num_states_for_printing;   //!< Number of states used for printing
    unsigned int num_print;                 //!< Number of outputs
    unsigned int* print_indices;            //!< List of indices in solution vectors where data is written to output (size is num_states_for_printing)
    OutputIntCallback **outputs_i;
    OutputDoubleCallback **outputs_d;
    char** output_names;
    char** output_specifiers;
    enum AsynchTypes* output_types;
    short int* output_sizes;

    OutputFunc output_func;

    //Peakflow stuff
    char* peakflow_function_name;
    PeakflowOutputCallback* peakflow_output;

    //unsigned int num_rainsteps;
    char* hydro_table;
    char* peak_table;
    char* dump_table;
    //char* dump_location;
    //char* halt_filename;
} GlobalVars;


/// This structure holds all the data for a link in the river system.
///
typedef struct Link
{
    RKMethod* method;                   //!< RK method to use for solving the ODEs for this link
    RKSolutionList* list;               //!< The list for the calculated numerical solution
    ErrorData* errorinfo;               //!< Error estimiation information for this link
    VEC params;                         //!< Parameters unique for the ODE for this link
    //IVEC iparams;                     //!< Integer (long) parameters for the ODE for this link

    DifferentialFunc *f;    //!< Right-hand side function for ODE                                            
    AlgebraicFunc *alg;     //!< Function for algebraic variables
    CheckStateFunc *state_check;    //!< Function to check what "state" the state variables are in (for discontinuities)
    JacobianFunc *Jacobian;         //!< Jacobian of right-hand side function
    RKSolverFunc *RKSolver;         //!< RK solver to use
    CheckConsistencyFunc *CheckConsistency; //!< Function to check state variables

    double h;                           //!< Current step size
    double last_t;                      //!< Last time in which a numerical solution was calculated
    double print_time;                  //!< Numerical solution is written to disk in increments of print_time
    double next_save;                   //!< Next time to write numerical solution to disk
    unsigned int ID;                    //!< ID for the link. This is how a link is referenced in data files
    unsigned int location;              //!< Index of this link in the system array
    short int ready;                    //!< Flag that is 1 if a step can be taken, 0 if not
    unsigned short int num_parents;     //!< Number of upstream links
    int disk_iterations;                //!< Number of iterations stored on disk
    double peak_time;                   //!< The time at which the largest discharge has occurred for this link
    VEC peak_value;                     //!< The value of the largest discharge for this link
    struct Link** parents;             //!< An array of all upstream links (parents)
    struct Link* child;                //!< The downstream link (child)
    int current_iterations;             //!< Number of stored iterations in list
    int steps_on_diff_proc;             //!< Number of steps for this link that are stored on another process
    int iters_removed;                  //!< Total number of iterations removed that has not been sent
    unsigned int distance;              //!< Maximum number of links upstream to get to an external link
    int rejected;                       //!< 0 if the previous step was accepted, 1 if rejected, 2 for discontinuity
    unsigned short int save_flag;       //!< 1 if saving data for this link, 0 if not
    unsigned short int peak_flag;       //!< 1 if saving peak flow data for this link, 0 if not
    //unsigned int** upstream;          //!< upstream[i] is a list of links (loc) upstream (inclusive) to parent i
    //unsigned int* numupstream;        //!< numupstream[i] is size of upstream[i]
    struct QVSData* qvs;               //!< Holds the discharge vs storage data
    //fpos_t pos;                       //!< Current location in temporary output file
    long int pos_offset;
    unsigned int expected_file_vals;    //!< Expected number of entries in temp output file
    unsigned short int dam;             //!< 0 if no dam at the link, 1 if dam present
    unsigned short int res;             //!< 0 if this link has no reservoir feed, 1 if it does
    unsigned int dim;                   //!< Dimension of the problem at this link
    unsigned int diff_start;            //!< Starting index of differential variables in solution vectors
    unsigned int no_ini_start;          //!< Starting index of differential variables not read from disk
    unsigned int num_dense;             //!< Number of states where dense output is calculated
    unsigned int* dense_indices;        //!< List of indices in solution where dense output is needed

    //Forcings data
    //unsigned int num_forcings;
    ForcingData** forcing_buff;         //!< Forcing data for this link
    double* forcing_change_times;       //!< Next time in which there is a change in rainfall, relative to last_t
    double* forcing_values;             //!< The current rainfall value for this link at time last_t
    unsigned int* forcing_indices;      //!< forcing_indices[i] has index of forcing_buff[i]->rainfall[*][0] that is currently used

    //For output data
    void* output_user;
    void* peakoutput_user;

    //For custom data
    void* user;

    //Parser data
    //struct Formula* equations;

    //Implicit Solver
    double last_eta;
    MAT JMatrix;
    MAT CoefMat;
    VEC* Z_i;
    VEC sol_diff;
    double h_old;
    double value_old;
    short int compute_J;
    short int compute_LU;

    //Discontinuity tracking
    int state;                          //!< The current state of the solution
    double* discont;                    //!< List of discontinuity times to step on
    unsigned int discont_count;         //!< Number of times in discont
    unsigned int discont_start;         //!< Starting index in discont
    unsigned int discont_end;           //!< Last index of a time in discont
    unsigned int discont_send_count;    //!< Number of times in discont_send and discont_order_send
    double* discont_send;               //!< List of discontinuity times to send downstream
    unsigned int* discont_order_send;   //!< List of discontinuity derivative orders to send downstream
} Link;

typedef struct Model
{
    SetParamSizesFunc *set_param_sizes;
    ConvertFunc *convert;
    RoutinesFunc *routines;
    PrecalculationsFunc *precalculations;
    InitializeEqsFunc *initialize_eqs;
    PartitionFunc *partition;
} Model;


/// This structure holds all the data for a forcing in the river system.
///
typedef struct Forcing
{
    unsigned int(*GetPasses)(struct Forcing*, double maxtime, struct ConnData* conninfo);
    double(*GetNextForcing)(struct Link*, unsigned int, unsigned int*, unsigned int, int*, struct GlobalVars*, struct Forcing*, struct ConnData*, unsigned int**, unsigned int);
    unsigned short int flag;
    char* filename;
    unsigned int increment;
    double file_time;
    unsigned int first_file;
    unsigned int last_file;
    unsigned int raindb_start_time; //!< This is the unix time corresponding to the local intergrator time 0
    //char* hydro_table;
    //char* peak_table;
    //char* dump_location;
    //char* halt_filename;
    //unsigned int num_rainsteps;
    ForcingData* GlobalForcing;
    unsigned int passes;
    int maxfileindex;
    double maxtime;
    unsigned int iteration;
    unsigned short int active;
    unsigned int good_timestamp;    //!< A timestamp in the db where a forcing actually exists
    double factor;
    char* lookup_filename;
    char* fileident;
    unsigned int** grid_to_linkid;
    unsigned int* num_links_in_grid;
    char* received;
    float* intensities;
    unsigned int num_cells;

    //For irregular timesteps
    unsigned int next_timestamp;        //!< Holds the next timestep to use for pulling data.
    unsigned int lastused_first_file;   //!< The value of first_file when the GetPasses routine was last called. 0 if never set.
    unsigned int lastused_last_file;    //!< The value of last_file when the GetPasses routine was last called. 0 if never set.
    unsigned int number_timesteps;      //!< The number of times which feature a forcing at some link.
} Forcing;


/// This structure holds information about how data is to be transfered between processes.
///
typedef struct TransData
{
    Link*** send_data;              //!< 2D array. send_data[i][j] points to a link about which data will be sent to process i.
    Link*** receive_data;           //!< 2D array. receive_data[i][j] points to a link about which data will be received from process i.
    unsigned int* send_size;        //!< send_size[i] has the number of entries in send_data[i].
    unsigned int* receive_size;     //!< receive_size[i] has the number of entries in receive_data[i].
    MPI_Request** send_requests;    //!< send_requests[i] has the request for data sent to process i.
    MPI_Request** receive_requests; //!< receive_requests[i] has request for data received from process i.
    char** send_buffer;             //!< 2D array. send_buffer[i][j] is a buffer for sending data for link send_data[i][j].
    char** receive_buffer;          //!< A buffer for receiving data through MPI
    short int* sent_flag;           //!< sent_flag[i] is 1 if a message has been sent to process i, 0 if not.
    short int* receiving_flag;      //!< receiving_flag[i] is 1 if a message is being received from process i.
    unsigned int* send_buffer_size;     //!< send_buffer_size[i] has the number of bytes stored in send_buffer[i]
    unsigned int* receive_buffer_size;  //!< receive_buffer_size[i] has the number of bytes stored in receive_buffer[i]
    unsigned int* num_sent;         //!< num_sent[i] is number of messages sent to process i
    unsigned int* num_recv;         //!< num_recv[i] is number of messages received from process i
    unsigned int* totals;           //!< workspace for flushing of size np
} TransData;

/// This is the main structure that holds the state of the server and associated data structures for a simulation.
///
typedef struct AsynchSolver
{
    //MPI Stuff
    MPI_Comm comm;		//!< COMM on which the solver works
    int np;			    //!< Number of procs in the comm
    int my_rank;		//!< This processes rank in the comm (varies by proc)

    //Routines for checking what is initialized
    short int setup_gbl;
    short int setup_topo;
    short int setup_params;
    short int setup_partition;
    short int setup_rkdata;
    short int setup_initmodel;
    short int setup_initconds;
    short int setup_forcings;
    short int setup_dams;
    short int setup_stepsizes;
    short int setup_savelists;
    short int setup_finalized;

    //Solver Stuff
    ErrorData* errors_tol;	    //!< Object for global error data
    GlobalVars* globals;		//!< Global information
    Link* sys;			        //!< Network of links
    RKMethod** AllMethods;		//!< List of RK methods
    TransData* my_data;		    //!< Data for communication between procs
    short int *getting;		    //!< List of data links to get information about
    int *assignments;		    //!< Link with sys location i is assigned to proc assignments[i]
    unsigned int* my_sys;		//!< Location in sys of links assigned to this proc
    unsigned int my_N;		    //!< Number of links in sys assigned to this proc
    unsigned int N;			    //!< Number of links in sys
    unsigned int nummethods;	//!< Number of methods in AllMethods
    unsigned int my_save_size;	//!< Number of links assigned to this proc in save_list
    unsigned int save_size;		//!< Number of links in save_list
    unsigned int *save_list;	//!< List of link ids to print data
    unsigned int *peaksave_list;
    unsigned int peaksave_size;	//Number of links to print peakflow data
    unsigned int my_peaksave_size;
    unsigned int *res_list;
    unsigned int res_size;
    unsigned int my_res_size;

    unsigned int** id_to_loc;	//!< Lookup table to convert from ids to sys locations
    TempStorage* workspace;		//!< Temporary workspace
    char rkdfilename[ASYNCH_MAX_PATH_LENGTH];	//!< Filename for .rkd file
    FILE* outputfile;		    //!< File handle for outputing temporary data
    FILE* peakfile;			    //!< File handle for the peakflow data
    ConnData db_connections[ASYNCH_MAX_DB_CONNECTIONS];	//!< Database connection information
    Forcing forcings[ASYNCH_MAX_DB_CONNECTIONS - ASYNCH_DB_LOC_FORCING_START];	//!< Forcing information
    Model* custom_model;
    void* ExternalInterface;
} AsynchSolver;


#endif //STRUCTS_H
