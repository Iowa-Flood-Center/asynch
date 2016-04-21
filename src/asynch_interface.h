#ifndef ASYNCH_INTERFACE_H
#define ASYNCH_INTERFACE_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "comm.h"
#include "riversys.h"
#include "processdata.h"
#include <time.h>
#include "structs.h"
#include <mpi.h>
#include "solvers.h"
#include <stdio.h>
#include <stdlib.h>
#include "io.h"
#include "data_types.h"

int my_rank;
int np;

typedef struct
{
	//MPI Stuff
	MPI_Comm comm;		//COMM on which the solver works
	int np;			//Number of procs in the comm
	int my_rank;		//This processes rank in the comm (varies by proc)

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
	short int mpi_initialized;	//1 if MPI was initialized by asynch, 0 if not

	//Solver Stuff
	ErrorData* GlobalErrors;	//Object for global error data
	UnivVars* GlobalVars;		//Global information
	Link** sys;			//Network of links
	RKMethod** AllMethods;		//List of RK methods
	TransData* my_data;		//Data for communication between procs
	short int *getting;		//List of data links to get information about
	int *assignments;		//Link with sys location i is assigned to proc assignments[i]
	unsigned int* my_sys;		//Location in sys of links assigned to this proc
	unsigned int my_N;		//Number of links in sys assigned to this proc
	unsigned int N;			//Number of links in sys
	unsigned int nummethods;	//Number of methods in AllMethods
	unsigned int my_save_size;	//Number of links assigned to this proc in save_list
	unsigned int save_size;		//Number of links in save_list
	unsigned int *save_list;	//List of link ids to print data
	unsigned int *peaksave_list;
	unsigned int peaksave_size;	//Number of links to print peakflow data
	unsigned int my_peaksave_size;
	unsigned int *res_list;
	unsigned int res_size;
	unsigned int my_res_size;

	unsigned int** id_to_loc;	//Table to convert from ids to sys locations
	TempStorage* workspace;		//Temporary workspace
	char rkdfilename[256];		//Filename for .rkd file
	FILE* outputfile;		//File for outputing temporary data
	FILE* peakfile;			//File for the peakflow data
	char* peakfilename;		//Filename for .pea file
	ConnData* db_connections[ASYNCH_MAX_DB_CONNECTIONS];	//Database connection information
	Forcing* forcings[ASYNCH_MAX_DB_CONNECTIONS - ASYNCH_DB_LOC_FORCING_START];	//Forcing information
	DataTypes dt_info;
	model* custom_model;
	void* ExternalInterface;
} asynchsolver;

//Constructor / Destructor related routings
asynchsolver* Asynch_Init(MPI_Comm comm,int* argc,char** argv[]);
int Asynch_Custom_Model(asynchsolver* asynch,
    void (*SetParamSizes)(UnivVars*,void*),
    void (*Convert)(VEC,unsigned int,void*),
    void (*Routines)(Link*,unsigned int,unsigned int,unsigned short int,void*),
	void (*Precalculations)(Link*,VEC,VEC,unsigned int,unsigned int,unsigned short int,unsigned int,void*),
    int (*InitializeEqs)(VEC,VEC,QVSData*,unsigned short int,VEC,unsigned int,unsigned int,unsigned int,void*,void*));
int Asynch_Custom_Partitioning(asynchsolver* asynch,int* (*Partition_Routine)(Link**,unsigned int,Link**,unsigned int,unsigned int**,unsigned int*,TransData*,short int*));
void Asynch_Free(asynchsolver* asynch);

//Routines to intialize network and model
void Asynch_Parse_GBL(asynchsolver* asynch,char* gbl_filename);
void Asynch_Load_Network(asynchsolver* asynch);
void Asynch_Partition_Network(asynchsolver* asynch);
void Asynch_Load_Network_Parameters(asynchsolver* asynch,short int load_all);
void Asynch_Load_Numerical_Error_Data(asynchsolver* asynch);
void Asynch_Initialize_Model(asynchsolver* asynch);
void Asynch_Load_Initial_Conditions(asynchsolver* asynch);
void Asynch_Load_Forcings(asynchsolver* asynch);
void Asynch_Load_Dams(asynchsolver* asynch);
void Asynch_Load_Save_Lists(asynchsolver* asynch);
void Asynch_Finalize_Network(asynchsolver* asynch);
void Asynch_Calculate_Step_Sizes(asynchsolver* asynch);

//Forcing routines
int Asynch_Activate_Forcing(asynchsolver* asynch,unsigned int idx);
int Asynch_Deactivate_Forcing(asynchsolver* asynch,unsigned int idx);

//Advance solver
void Asynch_Advance(asynchsolver* asynch,short int print_flag);

//Data file routines
void Asynch_Prepare_Output(asynchsolver* asynch);
void Asynch_Prepare_Temp_Files(asynchsolver* asynch);
void Asynch_Prepare_Peakflow_Output(asynchsolver* asynch);
int Asynch_Create_Output(asynchsolver* asynch,char* additional_out);
int Asynch_Create_Peakflows_Output(asynchsolver* asynch);
int Asynch_Delete_Temporary_Files(asynchsolver* asynch);
int Asynch_Write_Current_Step(asynchsolver* asynch);

//Snapshot
int Asynch_Take_System_Snapshot(asynchsolver* asynch,char* preface);

//Set and get routines
void Asynch_Set_Database_Connection(asynchsolver* asynch,char* database_info,unsigned int conn_idx);
double Asynch_Get_Total_Simulation_Time(asynchsolver* asynch);
void Asynch_Set_Total_Simulation_Time(asynchsolver* asynch,double new_time);
unsigned int Asynch_Get_Last_Rainfall_Timestamp(asynchsolver* asynch,unsigned int forcing_idx);
void Asynch_Set_Last_Rainfall_Timestamp(asynchsolver* asynch,unsigned int epoch_timestamp,unsigned int forcing_idx);
unsigned int Asynch_Get_First_Rainfall_Timestamp(asynchsolver* asynch,unsigned int forcing_idx);
void Asynch_Set_First_Rainfall_Timestamp(asynchsolver* asynch,unsigned int epoch_timestamp,unsigned int forcing_idx);
void Asynch_Set_RainDB_Starttime(asynchsolver* asynch,unsigned int epoch_timestamp,unsigned int forcing_idx);
void Asynch_Set_Init_File(asynchsolver* asynch,char* filename);
unsigned int Asynch_Get_Number_Links(asynchsolver* asynch);
unsigned int Asynch_Get_Local_Number_Links(asynchsolver* asynch);
void Asynch_Set_System_State(asynchsolver* asynch,double t_0,VEC* backup);
void Asynch_Reset_Peakflow_Data(asynchsolver* asynch);
int Asynch_Set_Forcing_State(asynchsolver* asynch,unsigned int idx,double t_0,unsigned int first_file,unsigned int last_file);
int Asynch_Set_Temp_Files(asynchsolver* asynch,double set_time,void* set_value,unsigned int output_idx);
int Asynch_Reset_Temp_Files(asynchsolver* asynch,double set_time);
int Asynch_Get_Peakflow_Output_Name(asynchsolver* asynch,char* peakflowname);
int Asynch_Set_Peakflow_Output_Name(asynchsolver* asynch,char* peakflowname);
unsigned int Asynch_Get_Local_LinkID(asynchsolver* asynch,unsigned int location);
int Asynch_Set_Init_Timestamp(asynchsolver* asynch,unsigned int epoch_timestamp);
unsigned int Asynch_Get_Init_Timestamp(asynchsolver* asynch);
int Asynch_Get_Snapshot_Output_Name(asynchsolver* asynch,char* filename);
int Asynch_Set_Snapshot_Output_Name(asynchsolver* asynch,char* filename);
int Asynch_Get_Reservoir_Forcing(asynchsolver* asynch);
unsigned int Asynch_Get_Size_Global_Parameters(asynchsolver* asynch);
int Asynch_Get_Global_Parameters(asynchsolver* asynch,VEC gparams);
int Asynch_Set_Global_Parameters(asynchsolver* asynch,VEC gparams,unsigned int n);

//Routines for output
int Asynch_Set_Output(asynchsolver* asynch,char* name,short int data_type,void (*func)(double,VEC,VEC,VEC,int,void*),unsigned int* used_states,unsigned int num_states);
int Asynch_Check_Output(asynchsolver* asynch,char* name);
int Asynch_Check_Peakflow_Output(asynchsolver* asynch,char* name);
int Asynch_Set_Peakflow_Output(asynchsolver* asynch,char* name,void (*func)(unsigned int,double,VEC,VEC,VEC,double,unsigned int,void*,char*));
int Asynch_Create_OutputUser_Data(asynchsolver* asynch,unsigned int data_size);
int Asynch_Free_OutputUser_Data(asynchsolver* asynch);
void Asynch_Copy_Local_OutputUser_Data(asynchsolver* asynch,unsigned int location,void* source,unsigned int size);
void Asynch_Set_Size_Local_OutputUser_Data(asynchsolver* asynch,unsigned int location,unsigned int size);

#endif

