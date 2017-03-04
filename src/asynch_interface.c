#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#if !defined(_MSC_VER)
#define ASYNCH_SLEEP sleep
#else
#include <windows.h>
#define ASYNCH_SLEEP Sleep
#endif

#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(HAVE_UNISTD_H)
#include <unistd.h>
#endif

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#include "structs.h"
#include "comm.h"
#include "riversys.h"
#include "processdata.h"
#include "solvers.h"
#include "io.h"
#include "data_types.h"

#include "asynch_interface.h"


//Initializes the asynch solver object.
AsynchSolver* Asynch_Init(MPI_Comm comm)
{
    AsynchSolver* res = NULL;

    unsigned int i;
    int init_flag;

    res = malloc(sizeof(AsynchSolver));
    memset(res, 0, sizeof(AsynchSolver));

    res->comm = comm;
    if (comm != MPI_COMM_WORLD)	printf("Warning: asynchsolver object my not work fully with in a comm other than MPI_COMM_WORLD.\n");

    //Initialize MPI stuff
    MPI_Initialized(&init_flag);
    if (!init_flag)
    {
        printf("Error: MPI must be initialized prior to Asynch_Init.\n");
        MPI_Abort(MPI_COMM_WORLD, 0);
    }

    MPI_Comm_rank(comm, &(res->my_rank));
    MPI_Comm_size(comm, &(res->np));

    //Sets the global variables
    np = res->np;
    my_rank = res->my_rank;

    //Initialize asynchsolver members
    for (i = 0; i < ASYNCH_MAX_DB_CONNECTIONS - ASYNCH_DB_LOC_FORCING_START; i++)
         Forcing_Init(&res->forcings[i]);

    return res;
}

//Creates and loads a custom model
//Return 1 if there is an error, 0 if everything is ok
int Asynch_Custom_Model(
    AsynchSolver* asynch,
    SetParamSizesFunc *set_param_sizes,
    ConvertFunc *convert,
    RoutinesFunc *routines,
    PrecalculationsFunc *precalculations,
    InitializeEqsFunc *initialize_eqs)
{
    if (!(asynch->custom_model))
    {
        asynch->custom_model = (Model*)malloc(sizeof(Model));
        asynch->custom_model->partition = NULL;
    }

    asynch->custom_model->set_param_sizes = set_param_sizes;
    asynch->custom_model->convert = convert;
    asynch->custom_model->routines = routines;
    asynch->custom_model->precalculations = precalculations;
    asynch->custom_model->initialize_eqs = initialize_eqs;

    return 0;
}

//Sets a routine to use for partitioning.
int Asynch_Custom_Partitioning(AsynchSolver* asynch, PartitionFunc *partition)
{
    if (!(asynch->custom_model))
    {
        asynch->custom_model = (Model*)malloc(sizeof(Model));
        memset(asynch->custom_model, 0, sizeof(Model));
    }

    asynch->custom_model->partition = partition;

    return 0;
}


// Start network setup ******************************************************************************************************

//Reads a .gbl file.
void Asynch_Parse_GBL(AsynchSolver* asynch, char* filename)
{
    //Read in .gbl file
    asynch->globals = Read_Global_Data(filename, &(asynch->errors_tol), asynch->forcings, asynch->db_connections, asynch->rkdfilename, asynch->custom_model, asynch->ExternalInterface);
    if (!asynch->globals)
    {
        if (my_rank == 0)
            printf("[%i]: An error occurred reading the .gbl file. See above messages for details.\n", my_rank);
        MPI_Abort(asynch->comm, 1);
    }

    asynch->setup_gbl = 1;
}

void Asynch_Load_Network(AsynchSolver* asynch)
{
    if (!asynch->setup_gbl)
    {
        if (my_rank == 0)
            printf("Error: global file must be read before loading the network topology.\n");
        MPI_Abort(asynch->comm, 1);
    }

    asynch->sys = Create_River_Network(asynch->globals, &(asynch->N), &(asynch->id_to_loc), asynch->db_connections);
    if (!asynch->sys)	MPI_Abort(asynch->comm, 1);
    asynch->setup_topo = 1;
    MPI_Barrier(asynch->comm);
}

//If load_all == 1, then the parameters for every link are available on every proc.
//If load_all == 0, then the parameters are only available for links assigned to this proc.
void Asynch_Load_Network_Parameters(AsynchSolver* asynch, short int load_all)
{
    int i;
    if (!asynch->setup_topo)
    {
        if (my_rank == 0)
            printf("Error: Topology data must be read before loading the link parameters.\n");
        else
            ASYNCH_SLEEP(1);
        MPI_Abort(asynch->comm, 1);
    }
    if (!asynch->setup_partition && !load_all)
    {
        if (my_rank == 0)
            printf("Error: Paritioning must be done before reading link parameters.\n");
        else
            ASYNCH_SLEEP(1);
        MPI_Abort(asynch->comm, 1);
    }

    i = Load_Local_Parameters(asynch->sys, asynch->N, asynch->my_sys, asynch->my_N, asynch->assignments, asynch->getting, asynch->id_to_loc, asynch->globals, asynch->db_connections, load_all, asynch->custom_model, asynch->ExternalInterface);
    if (i)	MPI_Abort(asynch->comm, 1);
    asynch->setup_params = 1;
    MPI_Barrier(asynch->comm);
}

void Asynch_Partition_Network(AsynchSolver* asynch)
{
    int i;
    if (!asynch->setup_topo)
    {
        if (my_rank == 0)
            printf("Error: Topology data must be read before partitioning the network.\n");
        else
            ASYNCH_SLEEP(1);
        MPI_Abort(asynch->comm, 1);
    }

    i = Partition_Network(asynch->sys, asynch->N, asynch->globals, &(asynch->my_sys), &(asynch->my_N), &(asynch->assignments), &(asynch->my_data), &(asynch->getting), asynch->custom_model);
    if (i)	MPI_Abort(asynch->comm, 1);
    asynch->setup_partition = 1;
    MPI_Barrier(asynch->comm);
}

void Asynch_Load_Numerical_Error_Data(AsynchSolver* asynch)
{
    int i;
    if (!asynch->setup_partition)
    {
        if (my_rank == 0)
            printf("Error: Partitioning must be done before setting up numerical error data.\n");
        else
            ASYNCH_SLEEP(1);
        MPI_Abort(asynch->comm, 1);
    }

    i = Build_RKData(asynch->sys, asynch->rkdfilename, asynch->N, asynch->my_sys, asynch->my_N, asynch->assignments, asynch->getting, asynch->globals, asynch->errors_tol, &(asynch->AllMethods), &(asynch->nummethods));
    if (i)	MPI_Abort(asynch->comm, 1);
    asynch->setup_rkdata = 1;
    MPI_Barrier(asynch->comm);
}

void Asynch_Initialize_Model(AsynchSolver* asynch)
{
    int i;
    if (!asynch->setup_partition)
    {
        if (my_rank == 0)
            printf("Error: Partitioning must be done before initializing model.\n");
        else
            ASYNCH_SLEEP(1);
        MPI_Abort(asynch->comm, 1);
    }
    if (!asynch->setup_dams)
    {
        if (my_rank == 0)
            printf("Error: Dam parameters must be loaded before initializing model.\n");
        else
            ASYNCH_SLEEP(1);
        MPI_Abort(asynch->comm, 1);
    }
    if (!asynch->setup_dams)
    {
        if (my_rank == 0)
            printf("Error: Save lists must be loaded before initializing model.\n");
        else
            ASYNCH_SLEEP(1);
        MPI_Abort(asynch->comm, 1);
    }
    i = Initialize_Model(asynch->sys, asynch->N, asynch->my_sys, asynch->my_N, asynch->assignments, asynch->getting, asynch->globals, asynch->custom_model, asynch->ExternalInterface);
    if (i)	MPI_Abort(asynch->comm, 1);
    asynch->setup_initmodel = 1;
    MPI_Barrier(asynch->comm);
}

void Asynch_Load_Initial_Conditions(AsynchSolver* asynch)
{
    int i;
    if (!asynch->setup_partition)
    {
        if (my_rank == 0)
            printf("Error: Partitioning must be done before loading initial conditions.\n");
        else
            ASYNCH_SLEEP(1);
        MPI_Abort(asynch->comm, 1);
    }
    if (!asynch->setup_params)
    {
        if (my_rank == 0)
            printf("Error: Link parameters must be loaded before loading initial conditions.\n");
        else
            ASYNCH_SLEEP(1);
        MPI_Abort(asynch->comm, 1);
    }
    if (!asynch->setup_dams)
    {
        if (my_rank == 0)
            printf("Error: Dam parameters must be loaded before loading initial conditions.\n");
        else
            ASYNCH_SLEEP(1);
        MPI_Abort(asynch->comm, 1);
    }

    i = Load_Initial_Conditions(asynch->sys, asynch->N, asynch->assignments, asynch->getting, asynch->id_to_loc, asynch->globals, asynch->db_connections, asynch->custom_model, asynch->ExternalInterface);
    if (i)	MPI_Abort(asynch->comm, 1);
    asynch->setup_initconds = 1;
    MPI_Barrier(asynch->comm);
}

void Asynch_Load_Forcings(AsynchSolver* asynch)
{
    int i;
    if (!asynch->setup_partition)
    {
        if (my_rank == 0)
            printf("Error: Partitioning must be done before loading forcings.\n");
        else
            ASYNCH_SLEEP(1);
        MPI_Abort(asynch->comm, 1);
    }

    i = Load_Forcings(asynch->sys, asynch->N, asynch->my_sys, asynch->my_N, asynch->assignments, asynch->getting, asynch->res_list, asynch->res_size, asynch->id_to_loc, asynch->globals, asynch->forcings, asynch->db_connections);
    if (i)	MPI_Abort(asynch->comm, 1);
    asynch->setup_forcings = 1;
    MPI_Barrier(asynch->comm);
}

void Asynch_Load_Dams(AsynchSolver* asynch)
{
    int i;
    if (!asynch->setup_partition)
    {
        if (my_rank == 0)
            printf("Error: Partitioning must be done before loading dam information.\n");
        else
            ASYNCH_SLEEP(1);
        MPI_Abort(asynch->comm, 1);
    }

    i = Load_Dams(asynch->sys, asynch->N, asynch->my_sys, asynch->my_N, asynch->assignments, asynch->getting, asynch->id_to_loc, asynch->globals, asynch->errors_tol, asynch->db_connections, &(asynch->res_list), &(asynch->res_size), &(asynch->my_res_size));
    if (i)	MPI_Abort(asynch->comm, 1);
    asynch->setup_dams = 1;
    MPI_Barrier(asynch->comm);
}

void Asynch_Calculate_Step_Sizes(AsynchSolver* asynch)
{
    int i;
    if (!asynch->setup_finalized)
    {
        if (my_rank == 0)
            printf("Error: Network must be finalized before calculating step sizes.\n");
        else
            ASYNCH_SLEEP(1);
        MPI_Abort(asynch->comm, 1);
    }

    i = CalculateInitialStepSizes(asynch->sys, asynch->my_sys, asynch->my_N, asynch->globals, asynch->workspace, 1);
    if (i)	MPI_Abort(asynch->comm, 1);
    asynch->setup_stepsizes = 1;
    MPI_Barrier(asynch->comm);
}

void Asynch_Load_Save_Lists(AsynchSolver* asynch)
{
    int i;
    if (!asynch->setup_partition)
    {
        if (my_rank == 0)
            printf("Error: Partitioning must be done before loading output data information.\n");
        else
            ASYNCH_SLEEP(1);
        MPI_Abort(asynch->comm, 1);
    }

    i = BuildSaveLists(asynch->sys, asynch->N, asynch->my_sys, asynch->my_N, asynch->assignments, asynch->id_to_loc, asynch->globals, &(asynch->save_list), &(asynch->save_size), &(asynch->my_save_size), &(asynch->peaksave_list), &(asynch->peaksave_size), &(asynch->my_peaksave_size), asynch->db_connections);
    if (i)	MPI_Abort(asynch->comm, 1);
    asynch->setup_savelists = 1;
    MPI_Barrier(asynch->comm);
}

void Asynch_Finalize_Network(AsynchSolver* asynch)
{
    int i;
    if (my_rank == 0)
    {
        if (!asynch->setup_gbl)
        {
            printf("Error: Cannot finalize network: Global file not loaded.\n");
            MPI_Abort(asynch->comm, 1);
        }
        else if (!asynch->setup_topo)
        {
            printf("Error: Cannot finalize network: Topology data not loaded.\n");
            MPI_Abort(asynch->comm, 1);
        }
        else if (!asynch->setup_params)
        {
            printf("Error: Cannot finalize network: Link parameters not loaded.\n");
            MPI_Abort(asynch->comm, 1);
        }
        else if (!asynch->setup_partition)
        {
            printf("Error: Cannot finalize network: Network partitioning not performed.\n");
            MPI_Abort(asynch->comm, 1);
        }
        else if (!asynch->setup_rkdata)
        {
            printf("Error: Cannot finalize network: Numerical error data not loaded.\n");
            MPI_Abort(asynch->comm, 1);
        }
        else if (!asynch->setup_initmodel)
        {
            printf("Error: Cannot finalize network: Model initialization not performed.\n");
            MPI_Abort(asynch->comm, 1);
        }
        else if (!asynch->setup_initconds)
        {
            printf("Error: Cannot finalize network: Initial conditions not loaded.\n");
            MPI_Abort(asynch->comm, 1);
        }
        else if (!asynch->setup_forcings)
        {
            printf("Error: Cannot finalize network: Forcings not loaded.\n");
            MPI_Abort(asynch->comm, 1);
        }
        else if (!asynch->setup_dams)
        {
            printf("Error: Cannot finalize network: Dam data not loaded.\n");
            MPI_Abort(asynch->comm, 1);
        }
        else if (!asynch->setup_savelists)
        {
            printf("Error: Cannot finalize network: Output information not loaded.\n");
            MPI_Abort(asynch->comm, 1);
        }
    }

    i = FinalizeSystem(asynch->sys, asynch->N, asynch->my_sys, asynch->my_N, asynch->assignments, asynch->getting, asynch->id_to_loc, asynch->my_data, asynch->globals, asynch->db_connections, &(asynch->workspace));
    if (i)	MPI_Abort(asynch->comm, 1);
    asynch->setup_finalized = 1;
    MPI_Barrier(asynch->comm);
}


//Trash an asynchsolver object
void Asynch_Free(AsynchSolver* asynch)
{
    unsigned int i;

    TransData_Free(asynch->my_data);
    for (i = 0; i < ASYNCH_MAX_DB_CONNECTIONS; i++)
        ConnData_Free(&asynch->db_connections[i]);
    Destroy_ErrorData(asynch->errors_tol);
    Destroy_Workspace(asynch->workspace, asynch->globals->max_s, asynch->globals->max_parents);
    free(asynch->workspace);
    free(asynch->getting);
    
    if (asynch->outputfile)
        fclose(asynch->outputfile);
    
    for (i = 0; i < asynch->N; i++)
        Destroy_Link(&asynch->sys[i], asynch->globals->iter_limit, asynch->rkdfilename[0] != '\0', asynch->forcings, asynch->globals);

    for (i = 0; i < ASYNCH_MAX_DB_CONNECTIONS - ASYNCH_DB_LOC_FORCING_START; i++)
        Forcing_Free(&asynch->forcings[i]);

    free(asynch->sys);
    free(asynch->my_sys);
    free(asynch->assignments);
    for (i = 0; i < asynch->nummethods; i++)	Destroy_RKMethod(asynch->AllMethods[i]);
    free(asynch->AllMethods);
    if (asynch->save_list)
        free(asynch->save_list);
    if (asynch->peaksave_list)
        free(asynch->peaksave_list);
    if (asynch->res_list)
        free(asynch->res_list);
    for (i = 0; i < asynch->N; i++)
        free(asynch->id_to_loc[i]);
    free(asynch->id_to_loc);
    Destroy_UnivVars(asynch->globals);
    if (asynch->custom_model)
        free(asynch->custom_model);
    free(asynch);
}

void Asynch_Advance(AsynchSolver* asynch, bool print_flag)
{
    if (print_flag && !AreOutputsSet(asynch->globals))
    {
        printf("[%i]: Warning: Solver advance requested with data output enabled, but not all outputs are initialized. Continuing solver without outputing data.\n", my_rank);
        print_flag = 0;
    }

    Advance(asynch->sys, asynch->N, asynch->my_sys, asynch->my_N, asynch->globals, asynch->assignments, asynch->getting, asynch->res_list, asynch->res_size,
        asynch->id_to_loc, asynch->workspace, asynch->forcings, asynch->db_connections, asynch->my_data, print_flag, asynch->outputfile);
}


//Returns 0 if snapshot was made, 1 if an error was encountered, -1 if no snapshot was taken
int Asynch_Take_System_Snapshot(AsynchSolver* asynch, char* preface)
{
    if (!(asynch->globals->output_func.CreateSnapShot))
        return -1;
    return asynch->globals->output_func.CreateSnapShot(asynch->sys, asynch->N, asynch->assignments, asynch->globals, preface, &asynch->db_connections[ASYNCH_DB_LOC_SNAPSHOT_OUTPUT]);
}


unsigned short Asynch_Get_Model_Type(AsynchSolver* asynch)
{
    return asynch->globals->type;
}

void Asynch_Set_Model_Type(AsynchSolver* asynch, unsigned short type)
{
    asynch->globals->type = type;
}

unsigned short Asynch_Get_Num_Links(AsynchSolver* asynch)
{
    if (!asynch)
        return 0;
    return asynch->N;
}

Link* Asynch_Get_Links(AsynchSolver* asynch)
{
    if (!asynch)
        return NULL;
    return asynch->sys;
}

unsigned short Asynch_Get_Num_Links_Proc(AsynchSolver* asynch)
{
    if (!asynch)
        return 0;
    return asynch->my_N;
}

Link* Asynch_Get_Links_Proc(AsynchSolver* asynch)
{
    if (!asynch)
        return NULL;
    return asynch->sys;
}


void Asynch_Set_Database_Connection(AsynchSolver* asynch, const char* connstring, unsigned int conn_idx)
{
    ConnData_Free(&asynch->db_connections[conn_idx]);
    ConnData_Init(&asynch->db_connections[conn_idx], connstring);
}

double Asynch_Get_Total_Simulation_Duration(AsynchSolver* asynch)
{
    return asynch->globals->maxtime;
}

void Asynch_Set_Total_Simulation_Duration(AsynchSolver* asynch, double new_time)
{
    asynch->globals->maxtime = new_time;
}

unsigned int Asynch_Get_Last_Forcing_Timestamp(AsynchSolver* asynch, unsigned int forcing_idx)
{
    return asynch->forcings[forcing_idx].first_file + (unsigned int)(60.0 * asynch->forcings[forcing_idx].maxtime);
    //return asynch->globals->first_file + (unsigned int) (60.0 * asynch->globals->maxtime);
}

void Asynch_Set_First_Forcing_Timestamp(AsynchSolver* asynch, unsigned int epoch_timestamp, unsigned int forcing_idx)
{
    asynch->forcings[forcing_idx].first_file = epoch_timestamp;
    //asynch->globals->first_file = epoch_timestamp;
}

void Asynch_Set_Last_Forcing_Timestamp(AsynchSolver* asynch, unsigned int epoch_timestamp, unsigned int forcing_idx)
{
    asynch->forcings[forcing_idx].last_file = epoch_timestamp;
    //asynch->globals->last_file = epoch_timestamp;
}

unsigned int Asynch_Get_First_Forcing_Timestamp(AsynchSolver* asynch, unsigned int forcing_idx)
{
    return asynch->forcings[forcing_idx].first_file;
    //return asynch->globals->first_file;
}

void Asynch_Set_Forcing_DB_Starttime(AsynchSolver* asynch, unsigned int epoch_timestamp, unsigned int forcing_idx)
{
    asynch->forcings[forcing_idx].raindb_start_time = epoch_timestamp;
    //asynch->globals->raindb_start_time = epoch_timestamp;
}

//Returns 0 if everything is good. Returns 1 if init_flag does not support a timestamp.
int Asynch_Set_Init_Timestamp(AsynchSolver* asynch, unsigned int epoch_timestamp)
{
    if (asynch->globals->init_flag != 3)	return 1;
    asynch->globals->init_timestamp = epoch_timestamp;
    return 0;
}

unsigned int Asynch_Get_Init_Timestamp(AsynchSolver* asynch)
{
    return asynch->globals->init_timestamp;
}

void Asynch_Set_Init_File(AsynchSolver* asynch, char* filename)
{
    sprintf(asynch->globals->init_filename, "%s", filename);

    int length;
    for (length = 0; length < 256; length++)
        if (filename[length] == '\0')	break;

    if (length < 3 || length == 256)
    {
        if (my_rank == 0)	printf("Error: Bad init filename: %s.\n", filename);
        MPI_Abort(asynch->comm, 1);
    }

    //Check what type of file
    if (filename[length - 1] == 'c')
    {
        if (filename[length - 2] == 'e' && filename[length - 3] == 'r' && filename[length - 4] == '.')
        {
            asynch->globals->init_flag = 2;
            return;
        }
    }

    if (filename[length - 1] == 'i' && filename[length - 2] == 'n' && filename[length - 3] == 'i')
    {
        if (filename[length - 4] == '.')
        {
            asynch->globals->init_flag = 0;
            return;
        }

        if (filename[length - 4] == 'u' && filename[length - 5] == '.')
        {
            asynch->globals->init_flag = 1;
            return;
        }
    }

    if (my_rank == 0)	printf("Error: Bad init filename: %s.\n", filename);
    MPI_Abort(asynch->comm, 1);
}


void Asynch_Prepare_Output(AsynchSolver* asynch)
{
    if (asynch->globals->output_func.PrepareOutput)
        asynch->globals->output_func.PrepareOutput(asynch->globals, &asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);
}


void Asynch_Prepare_Temp_Files(AsynchSolver* asynch)
{
    unsigned int i;

    //Check that all outputs are set
    for (i = 0; i < asynch->globals->num_outputs; i++)
    {
        if (asynch->globals->output_types[i] == ASYNCH_BAD_TYPE)
        {
            if (my_rank == 0)
                printf("Error: Time series output %s is not defined.\n", asynch->globals->output_names[i]);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    if (asynch->globals->output_func.PrepareTempOutput)
        asynch->outputfile = asynch->globals->output_func.PrepareTempOutput(asynch->sys, asynch->N, asynch->assignments, asynch->globals, asynch->save_list, asynch->save_size, asynch->my_save_size, NULL, asynch->id_to_loc);
    else
        asynch->outputfile = NULL;
}

//Writes the current state to the temp files, if the current state is at a print time.
//Returns 0 if ok, 1 if no temp file is open.
int Asynch_Write_Current_Step(AsynchSolver* asynch)
{
    Link *sys = asynch->sys, *current;
    unsigned int i, N = asynch->N, loc, save_size = asynch->save_size;
    int *assignments = asynch->assignments;
    double time_diff;

    if (asynch->my_save_size && !(asynch->outputfile))
    {
        printf("[%i]: Error writting step. No temporary file is open.\n", my_rank);
        return 1;
    }

    for (i = 0; i < save_size; i++)
    {
        loc = find_link_by_idtoloc(asynch->save_list[i], asynch->id_to_loc, N);
        current = &sys[loc];
        if (assignments[loc] == my_rank)
        {
            time_diff = fabs(current->next_save - current->last_t);
            if (time_diff / current->next_save < 1e-12 || ((fabs(current->next_save) < 1e-12) ? (time_diff < 1e-12) : 0))
            {
                //WriteStep(current->last_t,current->list->head->y_approx,asynch->globals,current->params,current->state,asynch->outputfile,current->output_user,&(current->pos));
                WriteStep(asynch->outputfile, current->ID, current->last_t, current->list->head->y_approx, asynch->globals, current->params, current->state, current->output_user, &(current->pos_offset));	//!!!! Should be tail? !!!!
                current->next_save += current->print_time;
                (current->disk_iterations)++;
            }
        }
    }

    return 0;
}

void Asynch_Prepare_Peakflow_Output(AsynchSolver* asynch)
{
    if (asynch->globals->peakflow_output == NULL)
    {
        if (my_rank == 0)	printf("Error: Peakflow function %s not defined.\n", asynch->globals->peakflow_function_name);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (asynch->globals->output_func.PreparePeakflowOutput)
        asynch->globals->output_func.PreparePeakflowOutput(asynch->globals, asynch->peaksave_size);
}

//Return 0 means ok, -1 means no data to output
int Asynch_Create_Output(AsynchSolver* asynch, char* additional_out)
{
    //Flush the transfer buffers
    //!!!! I'm really not sure if this should be here. Seems like either the sends in processdata should use a different tag, or
    // the flush should unpack data instead of trashing it and be put in the Asynch_Advance function call. !!!!
    Flush_TransData(asynch->my_data);

    if (asynch->globals->output_func.CreateOutput && asynch->globals->hydrosave_flag)
        return asynch->globals->output_func.CreateOutput(asynch->sys, asynch->globals, asynch->N, asynch->save_list, asynch->save_size, asynch->my_save_size, asynch->id_to_loc, asynch->assignments, NULL, additional_out, &asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT], &(asynch->outputfile));
    return -1;
}

//Return 0 means ok, -1 means no peakflows to output
int Asynch_Create_Peakflows_Output(AsynchSolver* asynch)
{
    if (asynch->globals->output_func.CreatePeakflowOutput && asynch->globals->peaksave_flag)
        return asynch->globals->output_func.CreatePeakflowOutput(asynch->sys, asynch->globals, asynch->N, asynch->assignments, asynch->peaksave_list, asynch->peaksave_size, asynch->id_to_loc, &asynch->db_connections[ASYNCH_DB_LOC_PEAK_OUTPUT]);
    return -1;
}

//Return 0 means ok, 1 means no peakflows to output, 2 means the filename is a bad format
int Asynch_Set_Peakflow_Output_Name(AsynchSolver* asynch, char* peakflowname)
{
    size_t l;

    if (asynch->globals->peaks_loc_filename)
    {
        l = strlen(peakflowname);
        if (l > 3 && peakflowname[l - 4] == '.' && peakflowname[l - 3] == 'p' && peakflowname[l - 2] == 'e' && peakflowname[l - 1] == 'a')
        {
            peakflowname[l - 4] = '\0';
            sprintf(asynch->globals->peaks_loc_filename, peakflowname);
            peakflowname[l - 4] = '.';
            return 0;
        }
        else	return 2;
    }
    else	return 1;
}

//Return 0 means ok, 1 means no peakflows to output
int Asynch_Get_Peakflow_Output_Name(AsynchSolver* asynch, char* peakflowname)
{
    if (asynch->globals->peaks_loc_filename)
    {
        sprintf(peakflowname, asynch->globals->peaks_loc_filename);
        return 0;
    }
    else
        return 1;
}

/*
unsigned int Asynch_Get_Number_Links(AsynchSolver* asynch)
{
    if (!asynch)	return 0;
    return asynch->N;
}

unsigned int Asynch_Get_Local_Number_Links(AsynchSolver* asynch)
{
    if (!asynch)	return 0;
    return asynch->my_N;
}

int Asynch_Upload_Hydrographs_Database(asynchsolver* asynch)
{
    return UploadHydrosDB(asynch->sys,asynch->globals,asynch->N,asynch->save_list,asynch->save_size,asynch->my_save_size,asynch->id_to_loc,asynch->assignments,NULL,asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);
}
*/

int Asynch_Set_Temp_Files(AsynchSolver* asynch, double set_time, void* set_value, unsigned int output_idx)
{
    return SetTempFiles(set_time, set_value, asynch->globals->output_types[output_idx], output_idx, asynch->sys, asynch->N, asynch->outputfile, asynch->globals, asynch->my_save_size, asynch->id_to_loc);
    //return SetTempFiles(set_time,asynch->sys,asynch->N,asynch->outputfile,asynch->globals,asynch->my_save_size,asynch->id_to_loc);
}

int Asynch_Reset_Temp_Files(AsynchSolver* asynch, double set_time)
{
    return ResetTempFiles(set_time, asynch->sys, asynch->N, asynch->outputfile, asynch->globals, asynch->my_save_size, asynch->id_to_loc);
}

int Asynch_Set_Output_Int(AsynchSolver* asynch, char* name, OutputIntCallback* callback, unsigned int* used_states, unsigned int num_states)
{
    Link *sys = asynch->sys, *current;
    unsigned int *my_sys = asynch->my_sys, my_N = asynch->my_N;
    unsigned int i, j, idx, loc, num_to_add, *states_to_add = NULL;
    GlobalVars* GlobalVars = asynch->globals;

    //Find index
    for (i = 0; i < asynch->globals->num_outputs; i++)
    {
        if (strcmp(name, asynch->globals->output_names[i]) == 0)
        {
            idx = i;
            break;
        }
    }

    if (i == asynch->globals->num_outputs)
    {
        printf("[%i]: Output %s not set.\n", my_rank, name);
        return 0;
    }

    //Add new output
    asynch->globals->output_types[idx] = ASYNCH_INT;
    asynch->globals->outputs[idx].out_int = callback;

    asynch->globals->output_sizes[idx] = GetByteSize(ASYNCH_INT);
    asynch->globals->output_specifiers[idx] = GetSpecifier(ASYNCH_INT);

    //Check if anything should be added to the dense_indices from used_states
    for (loc = 0; loc < my_N; loc++)
    {
        current = &sys[my_sys[loc]];
        num_to_add = 0;
        states_to_add = (unsigned int*)realloc(states_to_add, num_states * sizeof(unsigned int));
        for (i = 0; i < num_states; i++)
        {
            if (used_states[i] > current->dim)
                continue;	//State is not present at this link
            for (j = 0; j < current->num_dense; j++)
            {
                if (used_states[i] == current->dense_indices[j])
                {
                    states_to_add[num_to_add++] = used_states[i];
                    break;
                }
            }
        }

        if (num_to_add)
        {
            current->dense_indices = (unsigned int*)realloc(current->dense_indices, (current->num_dense + num_to_add) * sizeof(unsigned int));
            for (i = 0; i < num_to_add; i++)
                current->dense_indices[i + current->num_dense] = states_to_add[i];
            current->num_dense += num_to_add;
            merge_sort_1D(current->dense_indices, current->num_dense);
        }
    }

    if (states_to_add)	free(states_to_add);

    return 1;
}


int Asynch_Set_Output_Double(AsynchSolver* asynch, char* name, OutputDoubleCallback* callback, unsigned int* used_states, unsigned int num_states)
{
    Link *sys = asynch->sys, *current;
    unsigned int *my_sys = asynch->my_sys, my_N = asynch->my_N;
    unsigned int i, j, idx, loc, num_to_add, *states_to_add = NULL;
    GlobalVars* GlobalVars = asynch->globals;

    //Find index
    for (i = 0; i < asynch->globals->num_outputs; i++)
    {
        if (strcmp(name, asynch->globals->output_names[i]) == 0)
        {
            idx = i;
            break;
        }
    }

    if (i == asynch->globals->num_outputs)
    {
        printf("[%i]: Output %s not set.\n", my_rank, name);
        return 0;
    }

    //Add new output
    asynch->globals->output_types[idx] = ASYNCH_DOUBLE;
    asynch->globals->outputs[idx].out_double = callback;

    asynch->globals->output_sizes[idx] = GetByteSize(ASYNCH_DOUBLE);
    asynch->globals->output_specifiers[idx] = GetSpecifier(ASYNCH_DOUBLE);

    //Check if anything should be added to the dense_indices from used_states
    for (loc = 0; loc < my_N; loc++)
    {
        current = &sys[my_sys[loc]];
        num_to_add = 0;
        states_to_add = (unsigned int*)realloc(states_to_add, num_states * sizeof(unsigned int));
        for (i = 0; i < num_states; i++)
        {
            if (used_states[i] > current->dim)
                continue;	//State is not present at this link
            for (j = 0; j < current->num_dense; j++)
            {
                if (used_states[i] == current->dense_indices[j])
                {
                    states_to_add[num_to_add++] = used_states[i];
                    break;
                }
            }
        }

        if (num_to_add)
        {
            current->dense_indices = (unsigned int*)realloc(current->dense_indices, (current->num_dense + num_to_add) * sizeof(unsigned int));
            for (i = 0; i < num_to_add; i++)
                current->dense_indices[i + current->num_dense] = states_to_add[i];
            current->num_dense += num_to_add;
            merge_sort_1D(current->dense_indices, current->num_dense);
        }
    }

    if (states_to_add)	free(states_to_add);

    return 1;
}


int Asynch_Set_Output_Float(AsynchSolver* asynch, char* name, OutputFloatCallback* callback, unsigned int* used_states, unsigned int num_states)
{
    Link *sys = asynch->sys, *current;
    unsigned int *my_sys = asynch->my_sys, my_N = asynch->my_N;
    unsigned int i, j, idx, loc, num_to_add, *states_to_add = NULL;
    GlobalVars* GlobalVars = asynch->globals;

    //Find index
    for (i = 0; i < asynch->globals->num_outputs; i++)
    {
        if (strcmp(name, asynch->globals->output_names[i]) == 0)
        {
            idx = i;
            break;
        }
    }

    if (i == asynch->globals->num_outputs)
    {
        printf("[%i]: Output %s not set.\n", my_rank, name);
        return 0;
    }

    //Add new output
    asynch->globals->output_types[idx] = ASYNCH_FLOAT;
    asynch->globals->outputs[idx].out_float = callback;

    asynch->globals->output_sizes[idx] = GetByteSize(ASYNCH_FLOAT);
    asynch->globals->output_specifiers[idx] = GetSpecifier(ASYNCH_FLOAT);

    //Check if anything should be added to the dense_indices from used_states
    for (loc = 0; loc < my_N; loc++)
    {
        current = &sys[my_sys[loc]];
        num_to_add = 0;
        states_to_add = (unsigned int*)realloc(states_to_add, num_states * sizeof(unsigned int));
        for (i = 0; i < num_states; i++)
        {
            if (used_states[i] > current->dim)
                continue;	//State is not present at this link
            for (j = 0; j < current->num_dense; j++)
            {
                if (used_states[i] == current->dense_indices[j])
                {
                    states_to_add[num_to_add++] = used_states[i];
                    break;
                }
            }
        }

        if (num_to_add)
        {
            current->dense_indices = (unsigned int*)realloc(current->dense_indices, (current->num_dense + num_to_add) * sizeof(unsigned int));
            for (i = 0; i < num_to_add; i++)
                current->dense_indices[i + current->num_dense] = states_to_add[i];
            current->num_dense += num_to_add;
            merge_sort_1D(current->dense_indices, current->num_dense);
        }
    }

    if (states_to_add)	free(states_to_add);

    return 1;
}


//Returns 1 if output function set successfully, 0 if there was a problem
int Asynch_Set_Peakflow_Output(AsynchSolver* asynch, char* name, PeakflowOutputCallback* callback)
{
    if (!asynch->globals->peakflow_function_name)
        asynch->globals->peakflow_function_name = (char*)malloc(asynch->globals->string_size * sizeof(char));

    if (strcmp(name, asynch->globals->peakflow_function_name))
    {
        printf("[%i]: Error setting peakflow output function to %s. Function %s was previously specified.\n", my_rank, asynch->globals->peakflow_function_name, name);
        return 0;
    }

    asynch->globals->peakflow_output = callback;

    return 1;
}

//Returns the link id of the link at my_sys[location]
unsigned int Asynch_Get_Local_LinkID(AsynchSolver* asynch, unsigned int location)
{
    return asynch->sys[asynch->my_sys[location]].ID;
}


//Copies data dump filename to filename
//Returns 0 if filename was copied. 1 if an error occured.
//!!!! Make sure filename has enough space !!!!
int Asynch_Get_Snapshot_Output_Name(AsynchSolver* asynch, char* filename)
{
    if (asynch->globals->dump_loc_filename == NULL)	return 1;
    strcpy(filename, asynch->globals->dump_loc_filename);
    return 0;
}

//Copies filename to data dump filename
//Returns 0 if filename was copied. 1 if an error occured.
int Asynch_Set_Snapshot_Output_Name(AsynchSolver* asynch, char* filename)
{
    if (strlen(filename) > asynch->globals->string_size)	return 1;
    if (!asynch->globals->dump_loc_filename)
        asynch->globals->dump_loc_filename = (char*)malloc(asynch->globals->string_size * sizeof(char));
    strcpy(asynch->globals->dump_loc_filename, filename);
    return 0;
}


//Returns 1 if output name is set,
//0 if output name is not set,
//-1 if output name is not present
int Asynch_Check_Output(AsynchSolver* asynch, char* name)
{
    unsigned int i;

    for (i = 0; i < asynch->globals->num_outputs; i++)
    {
        if (strcmp(name, asynch->globals->output_names[i]) == 0)
        {
            if (asynch->globals->output_types[i] == ASYNCH_BAD_TYPE)	return 0;
            else								return 1;
        }
    }

    return -1;
}

//Returns 1 if using the peakflow output name is set,
//0 if the peakflow output name is not specified,
//-1 if the peakflow output name is not present
int Asynch_Check_Peakflow_Output(AsynchSolver* asynch, char* name)
{
    if (!asynch->globals->peakflow_function_name || !name || strcmp(name, asynch->globals->peakflow_function_name))	return -1;
    return PeakflowOutputsSet(asynch->globals);
}


int Asynch_Delete_Temporary_Files(AsynchSolver* asynch)
{
    int ret_val = RemoveTemporaryFiles(asynch->globals, asynch->my_save_size, NULL);
    //if(ret_val == 1)	printf("[%i]: Error deleting temp file. File does not exist.\n");
    if (ret_val != 0)	printf("[%i]: Error [%i] while deleting temp file.\n", my_rank, ret_val);
    return ret_val;
}

//Return 0 if ok, 1 if error
//!!!! This should set the current forcing to something. But what if maxtime is beyond the ceiling term? !!!!
int Asynch_Activate_Forcing(AsynchSolver* asynch, unsigned int idx)
{
    //unsigned int i,j,l,my_N = asynch->my_N;
    unsigned int my_N = asynch->my_N;
    /*	Link* current;*/

    if (idx >= asynch->globals->num_forcings)
    {
        printf("[%i]: Cannot activate forcing %u. Not enough forcings.\n", my_rank, idx);
        return 1;
    }

    asynch->forcings[idx].active = 1;

    /*
    printf("Here\n");
    getchar();
        //Set the forcing value at each link
        for(i=0;i<my_N;i++)
        {
            current = asynch->sys[asynch->my_sys[i]];
    printf("rainfall buffer  t = %e\n",current->last_t);
    printf("*********\n");
    for(l=0;l<current->forcing_buff[idx]->n_times;l++)
    {
    printf("%e %e\n",current->forcing_buff[idx]->rainfall[l][0],current->forcing_buff[idx]->rainfall[l][1]);
    getchar();
    }
    printf("*********\n");

            //Find the right index in rainfall
            for(l=0;l<current->forcing_buff[idx]->n_times;l++)
    {
    printf("l = %u/%u\n",l,current->forcing_buff[idx]->n_times);
    printf("last t = %e\n",current->last_t);
    printf("buffer time = %e\n",current->forcing_buff[idx]->rainfall[l][0]);
    getchar();
                if( fabs(current->last_t - current->forcing_buff[idx]->rainfall[l][0]) < 1e-8 )	break;
    }

            double forcing_buffer = current->forcing_buff[idx]->rainfall[l][1];
            current->forcing_values[idx] = forcing_buffer;

            //Find and set the new change in rainfall
            for(j=l+1;j<current->forcing_buff[idx]->n_times;j++)
            {
                if( fabs(current->forcing_buff[idx]->rainfall[j][1] - forcing_buffer) > 1e-8 )
                {
                    current->forcing_change_times[idx] = current->forcing_buff[idx]->rainfall[j][0];
                    break;
                }
            }
            if(j == current->forcing_buff[idx]->n_times)
                current->forcing_change_times[idx] = current->forcing_buff[idx]->rainfall[j-1][0];
        }
    */

    return 0;
}

//Returns 0 if ok, 1 if error
int Asynch_Deactivate_Forcing(AsynchSolver* asynch, unsigned int idx)
{
    unsigned int i, my_N = asynch->my_N, *my_sys = asynch->my_sys;
    Link* sys = asynch->sys;

    if (idx >= asynch->globals->num_forcings)
    {
        printf("[%i]: Cannot deactivate forcing %u. Not enough forcings.\n", my_rank, idx);
        return 1;
    }

    //Deactivate forcing
    asynch->forcings[idx].active = 0;

    //Clear forcing values from links
    for (i = 0; i < my_N; i++)
    {
        sys[my_sys[i]].forcing_values[idx] = 0.0;
        sys[my_sys[i]].forcing_change_times[idx] = asynch->globals->maxtime + 1.0;
    }

    return 0;
}

//Sets information for a forcing.
//Returns 0 if everything is ok, 1 if there is an error.
int Asynch_Set_Forcing_State(AsynchSolver* asynch, unsigned int idx, double t_0, unsigned int first_file, unsigned int last_file)
{
    if (asynch->globals->num_forcings <= idx)
    {
        printf("[%i]: Error setting forcing state. Bad forcing index %u / %u.\n", my_rank, idx, asynch->globals->num_forcings);
        return 1;
    }

    asynch->forcings[idx].raindb_start_time = first_file;
    asynch->forcings[idx].maxtime = t_0;
    asynch->forcings[idx].iteration = 0;
    asynch->forcings[idx].first_file = first_file;
    asynch->forcings[idx].last_file = last_file;

    return 0;
}

//Resets the peakflow information at each link.
//The current time is used for the time to peak, and the last state is used for the values.
void Asynch_Reset_Peakflow_Data(AsynchSolver* asynch)
{
    unsigned int i, my_N = asynch->my_N, *my_sys = asynch->my_sys;
    Link *current, *sys = asynch->sys;
    double t_0 = sys[my_sys[0]].last_t;

    for (i = 0; i < my_N; i++)
    {
        current = &sys[my_sys[i]];
        current->peak_time = t_0;
        v_copy(current->list->tail->y_approx, current->peak_value);
    }
}

void Asynch_Set_System_State(AsynchSolver* asynch, double unix_time, VEC* states)
{
    unsigned i, j, k, l;
    Link* current;

    //Unpack asynch
    Link* sys = asynch->sys;
    unsigned int N = asynch->N, num_forcings = asynch->globals->num_forcings;
    //Forcing** forcings = asynch->forcings;
    GlobalVars* GlobalVars = asynch->globals;

    //Reset some things
    Flush_TransData(asynch->my_data);
    Asynch_Reset_Temp_Files(asynch, unix_time);

    //Reset links
    for (i = 0; i < N; i++)
    {
        current = &sys[i];
        if (current->list != NULL)
        {
            while (current->current_iterations > 1)
            {
                Remove_Head_Node(current->list);
                (current->current_iterations)--;
            }
            current->list->head->t = unix_time;
            current->last_t = unix_time;
            current->steps_on_diff_proc = 1;
            current->iters_removed = 0;
            current->rejected = 0;
            if (current->num_parents == 0)
                current->ready = 1;
            else
                current->ready = 0;


            for (j = 0; j < current->dim; j++)
                current->list->head->y_approx.ve[j] = states[i].ve[j];
            v_copy(states[i], current->list->head->y_approx);

            /*
                        //Reset the next_save time
                        if(current->save_flag)
                        {
                            current->next_save = t_0;		//!!!! This forces the print times to match up with the assimilation times !!!!
                            current->disk_iterations = 1;
                        }
            */

            //Reset peak flow information
            current->peak_time = unix_time;
            v_copy(current->list->head->y_approx, current->peak_value);

            //Reset current state
            if (current->state_check != NULL)
                current->state = current->state_check(current->list->head->y_approx, GlobalVars->global_params, current->params, current->qvs, current->dam);
            current->list->head->state = current->state;

            //Write initial state
            //if(current->save_flag)
            //	WriteStep(t_0,current->list->head->y_approx,asynch->globals,current->params,current->state,asynch->outputfile,current->output_user,&(current->pos));

            //Set forcings
            //!!!! This block was not here before I started toying with data assimilation stuff. Perhaps it causes problems with the forecasters... !!!!
            if (current->forcing_buff)
            {
                for (k = 0; k < num_forcings; k++)
                {
                    if (!(asynch->forcings[k].flag))	continue;

                    //Find the right index in rainfall
                    for (l = 0; l < current->forcing_buff[k]->nrows - 1; l++)
                        if (current->forcing_buff[k]->data[l][0] <= unix_time && unix_time < current->forcing_buff[k]->data[l + 1][0])	break;
                    double rainfall_buffer = current->forcing_buff[k]->data[l][1];
                    current->forcing_values[k] = rainfall_buffer;
                    current->forcing_indices[k] = l;

                    //Find and set the new change in rainfall
                    for (j = l + 1; j < current->forcing_buff[k]->nrows; j++)
                    {
                        if (fabs(current->forcing_buff[k]->data[j][1] - rainfall_buffer) > 1e-12)
                        {
                            current->forcing_change_times[k] = current->forcing_buff[k]->data[j][0];
                            break;
                        }
                    }
                    if (j == current->forcing_buff[k]->nrows)
                        current->forcing_change_times[k] = current->forcing_buff[k]->data[j - 1][0];

                    //Select new step size
                    //current->h = InitialStepSize(current->last_t,current,GlobalVars,workspace);
                }
            }
        }
    }
}

//Allocates space for output_user at each link.
//Returns 0 if everything is good, 1 if space is already allocated, 2 if an error occurred.
int Asynch_Create_OutputUser_Data(AsynchSolver* asynch, unsigned int data_size)
{
    if (!asynch)	return 2;

    unsigned int my_N = asynch->my_N, i, *my_sys = asynch->my_sys;
    Link* sys = asynch->sys;

    if (sys[my_sys[0]].output_user)	return 1;

    for (i = 0; i < my_N; i++)
        sys[my_sys[i]].output_user = malloc(data_size);

    return 0;
}

//Deallocates space for output_user at each link.
//Returns 0 if everything is good, 1 if space is already deallocated, 2 if an error occurred.
int Asynch_Free_OutputUser_Data(AsynchSolver* asynch)
{
    if (!asynch)	return 2;

    unsigned int my_N = asynch->my_N, i, *my_sys = asynch->my_sys;
    Link* sys = asynch->sys;

    if (sys[my_sys[0]].output_user == NULL)	return 1;

    for (i = 0; i < my_N; i++)
    {
        free(sys[my_sys[i]].output_user);
        sys[my_sys[i]].output_user = NULL;
    }

    return 0;
}

//Copies source into output_user for the link with location my_sys[location]
void Asynch_Copy_Local_OutputUser_Data(AsynchSolver* asynch, unsigned int location, void* source, unsigned int size)
{
    //printf("Copying size %u %u %p\n",size,location,asynch->sys[asynch->my_sys[location]]->output_user);
    memcpy(&asynch->sys[asynch->my_sys[location]].output_user, source, size);
}

//Allocates space for output_user for the link with location my_sys[location]
void Asynch_Set_Size_Local_OutputUser_Data(AsynchSolver* asynch, unsigned int location, unsigned int size)
{
    asynch->sys[asynch->my_sys[location]].output_user = malloc(size);
    //printf("Setting to size %u %p\n",size,asynch->sys[asynch->my_sys[location]]->output_user);
}

//Returns the forcing index used for reservoirs. Returns -1 if reservoir forcing is not set.
int Asynch_Get_Reservoir_Forcing(AsynchSolver* asynch)
{
    if (asynch->globals->res_flag == 0)	return -1;
    else	return asynch->globals->res_forcing_idx;
}

//Returns the number of global parameters in the system.
unsigned int Asynch_Get_Size_Global_Parameters(AsynchSolver* asynch)
{
    return asynch->globals->global_params.dim;
}

//Returns in gparams the global parameters of the system.
//This assumes gparams has enough memory allocated.
//Returns 1 if an error occurred. 0 otherwise.
void Asynch_Get_Global_Parameters(AsynchSolver* asynch, VEC params)
{
    if (!asynch || !asynch->globals)
        return;

    v_copy(asynch->globals->global_params, params);
}

//Sets the global parameters. n is the number of new parameters.
//Internally, memory is allocated/deallocated automatically.
//Returns 1 if an error occurred. 0 otherwise.
int Asynch_Set_Global_Parameters(AsynchSolver* asynch, VEC gparams, unsigned int n)
{
    if (!asynch || !asynch->globals)	return 1;
    v_resize(&asynch->globals->global_params, n);
    v_copy_n(gparams, asynch->globals->global_params, n);
    return 0;
}
