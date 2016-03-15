#include <stdio.h>
#include <time.h>
#include <libpq-fe.h>
#include <string.h>
#if !defined(_MSC_VER)
#include <unistd.h>
#endif
#include "asynch_interface.h"
#include "forecaster_methods.h"

int my_rank;
int np;

typedef struct CustomParams
{
	unsigned int ID;
	int offset;
} CustomParams;

int Output_Linkid(double t,VEC* y_i,VEC* global_params,VEC* params,int state,void* user);
int Output_Timestamp(double t,VEC* y_i,VEC* global_params,VEC* params,int state,void* user);

void Init_Output_User_forecastparams(asynchsolver* asynch);
void Free_Output_User_forecastparams(asynchsolver* asynch);
void Set_Output_User_forecastparams(asynchsolver* asynch,unsigned int offset);

void Init_Output_PeakflowUser_Offset(asynchsolver* asynch);
void Free_Output_PeakflowUser_Offset(asynchsolver* asynch);
void Set_Output_PeakflowUser_Offset(asynchsolver* asynch,unsigned int offset);


int main(int argc,char* argv[])
{
	//Initialize MPI stuff
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&np);

	//Parse input
	if(argc < 3)
	{
		if(my_rank == 0)
			printf("Command line parameter required:\nA universal variable file (.gbl),\nA forecast file (.fcst)\n");		
		MPI_Finalize();
		return 1;
	}

	//Declare variables
	unsigned int i,j,k,current_offset;
	double holder,longest,total_time = 0.0;
	time_t start,start2,stop;
	asynchsolver* asynch;
	PGresult *res;
	MPI_Status status;
	char* query = (char*) malloc(1024*sizeof(char));
	Link* current;

	if(my_rank == 0)
		printf("\nBeginning initialization...\n*****************************\n");
	MPI_Barrier(MPI_COMM_WORLD);
	start = time(NULL);

	//Init asynch object and the river network
	asynch = Asynch_Init(MPI_COMM_WORLD,&argc,&argv);
	Asynch_Parse_GBL(asynch,argv[1]);

	//Load Forecast related data
	ForecastData* Forecaster = Init_ForecastData(argv[2],asynch->GlobalVars->string_size);
	if(!Forecaster)
		MPI_Abort(MPI_COMM_WORLD,1);
	double forecast_time = Forecaster->forecast_window;

	holder = Asynch_Get_Total_Simulation_Time(asynch);
	longest = (holder < forecast_time) ? forecast_time : holder;
	Asynch_Set_Total_Simulation_Time(asynch,longest);
	if(my_rank == 0)	printf("Loading network...\n");
	Asynch_Load_Network(asynch);
	if(my_rank == 0)	printf("Partitioning network...\n");
	Asynch_Partition_Network(asynch);
	if(my_rank == 0)	printf("Loading parameters...\n");
	Asynch_Load_Network_Parameters(asynch,0);
	if(my_rank == 0)	printf("Reading dam and reservoir data...\n");
	Asynch_Load_Dams(asynch);
	if(my_rank == 0)	printf("Setting up numerical error data...\n");
	Asynch_Load_Numerical_Error_Data(asynch);
	if(my_rank == 0)	printf("Initializing model...\n");
	Asynch_Initialize_Model(asynch);
	if(my_rank == 0)	printf("Loading initial conditions...\n");
	Asynch_Load_Initial_Conditions(asynch);
	if(my_rank == 0)	printf("Loading forcings...\n");
	Asynch_Load_Forcings(asynch);
	if(my_rank == 0)	printf("Loading output data information...\n");
	Asynch_Load_Save_Lists(asynch);
	if(my_rank == 0)	printf("Finalizing network...\n");
	Asynch_Finalize_Network(asynch);
	if(my_rank == 0)	printf("Calculating initial step sizes...\n");
	Asynch_Calculate_Step_Sizes(asynch);

	Asynch_Set_Total_Simulation_Time(asynch,holder);

	//Setup output for link id, if needed
	int setup_id = Asynch_Check_Output(asynch,"LinkID");
	int setup_timestamp = Asynch_Check_Output(asynch,"Timestamp");
	if( (setup_id || setup_timestamp) != 0)
	{
		if(my_rank == 0)	printf("[%i]: Forecaster needs LinkID (%i), Timestamp (%i).\n",my_rank,setup_id,setup_timestamp);
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	Init_Output_User_forecastparams(asynch);
	Asynch_Set_Output(asynch,"LinkID",ASYNCH_INT,(void (*)(double,VEC*,VEC*,VEC*,int,void*)) &Output_Linkid,NULL,0);
	Asynch_Set_Output(asynch,"Timestamp",ASYNCH_INT,(void (*)(double,VEC*,VEC*,VEC*,int,void*)) &Output_Timestamp,NULL,0);

	Init_Output_PeakflowUser_Offset(asynch);

	//Get some values about the river system
	unsigned int N = Asynch_Get_Number_Links(asynch);
	unsigned int my_N = Asynch_Get_Local_Number_Links(asynch);
	char dump_filename[ASYNCH_MAX_PATH_LENGTH],filename[ASYNCH_MAX_PATH_LENGTH];

	//Create halt file
	CreateHaltFile(Forecaster->halt_filename);

	//Find the index of the forcing to use for forecasting
	unsigned int forecast_idx = Forecaster->forecasting_forcing;
	if(forecast_idx >= asynch->GlobalVars->num_forcings)
	{
		if(my_rank == 0)	printf("[%i]: Error: No forecasting forcing set.\n",my_rank);
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	//Reserve space for backups
	VEC** backup = (VEC**) malloc(N*sizeof(VEC*));
	for(i=0;i<N;i++)
	{
		if(asynch->assignments[i] == my_rank || asynch->getting[i] == 1)
			backup[i] = v_get(asynch->sys[i]->dim);
		else
			backup[i] = NULL;
	}

	if(my_rank == 0)
	{
		printf("\nModel type is %u.\nGlobal parameters are:\n",asynch->GlobalVars->type);
		Print_Vector(asynch->GlobalVars->global_params);
		printf("\n");
	}

	MPI_Barrier(MPI_COMM_WORLD);
	stop = time(NULL);
	total_time = difftime(stop,start);
	if(my_rank == 0)	printf("Finished initializations. Total time: %f\n",total_time);
	MPI_Barrier(MPI_COMM_WORLD);
	ASYNCH_SLEEP(1);

	//Make sure everyone is good before getting down to it...
	printf("Process %i (%i total) is good to go with %i links.\n",my_rank,np,my_N);
	MPI_Barrier(MPI_COMM_WORLD);
	start = time(NULL);

	//Make the initial solve
	Asynch_Advance(asynch,0);

	//Stop the clock
	MPI_Barrier(MPI_COMM_WORLD);
	stop = time(NULL);
	total_time += difftime(stop,start);

	//Output some data
	if(my_rank == 0)
	{
		printf("%i: The answer at ID %i at time %.12f is\n",my_rank,asynch->sys[asynch->my_sys[0]]->ID,asynch->sys[asynch->my_sys[0]]->last_t);
		Print_Vector(asynch->sys[asynch->my_sys[0]]->list->tail->y_approx);
		printf("Total time for calculations: %f\n",difftime(stop,start));
	}

	//Begin persistent calculations
	if(my_rank == 0)
		printf("\n\n===================================\nBeginning persistent calculations\n===================================\n");
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);

	//Make some initializations and checks
	short unsigned int hr1 = 0;	//Hour of the day to perform maintainance on database
	short unsigned int hr2 = 12;	//Hour of the day to perform maintainance on database
	unsigned int wait_time = 120;	//Time to sleep if no rainfall data is available
	unsigned int num_tables = 10;
	unsigned int num_rainsteps = Forecaster->num_rainsteps;	//Number of rainfall intensities to use for the next forecast
	if(my_rank == 0 && asynch->forcings[forecast_idx]->increment < num_rainsteps + 3)
		printf("Warning: Increment for rain should probably be %u.\n",num_rainsteps + 3);
	asynch->forcings[forecast_idx]->increment = num_rainsteps;	//!!!! Not necessary, but makes me feel better. The solvers should really not do the last step where they download nothing. !!!!

	unsigned int nextraintime,nextforcingtime;
	short int halt = 0;
	int isnull,repeat_for_errors;
	short int vac = 0;	//0 if no vacuum has occured, 1 if vacuum has occured (during a specific hour)
	unsigned int last_file = asynch->forcings[forecast_idx]->last_file;
	unsigned int first_file = asynch->forcings[forecast_idx]->first_file;
	k = 0;
	for(i=0;i<N;i++)
		if(backup[i] != NULL)	v_copy(asynch->sys[i]->list->tail->y_approx,backup[i]);

	double simulation_time_with_data = 0.0;
	simulation_time_with_data = max(simulation_time_with_data,asynch->forcings[forecast_idx]->file_time * Forecaster->num_rainsteps);

	//Setup temp files
	Set_Output_User_forecastparams(asynch,first_file);
	Set_Output_PeakflowUser_Offset(asynch,first_file);
	Asynch_Set_Total_Simulation_Time(asynch,forecast_time);
	Asynch_Prepare_Temp_Files(asynch);

	//Prepare snapshots
	Asynch_Get_Snapshot_Output_Name(asynch,dump_filename);
	dump_filename[strlen(dump_filename)-4] = '\0';	//Removes .rec	!!!! Uh, is this ok? No chance for memory corruption? !!!!

	//Check if there is a schema used for the hydrograph archive
	int place,tablename_len = strlen(asynch->GlobalVars->hydro_table);
	char schema[128]; schema[0] = '\0';
	for(place=tablename_len-1;place>-1;place--)
	{
		if(asynch->GlobalVars->hydro_table[place] == '.')
		{
			asynch->GlobalVars->hydro_table[place] = '\0';
			strcpy(schema,asynch->GlobalVars->hydro_table);
			asynch->GlobalVars->hydro_table[place] = '.';
			schema[place] = '.';
			schema[place+1] = '\0';
			break;
		}
	}

	//Make some initializations to the database
	if(my_rank == 0)
	{
		printf("Making initializations to tables.\n");
		start = time(NULL);

		//Connect to hydrograph database
		ConnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);

		//Make sure the hydrographs table exists
		sprintf(query,"SELECT 1 FROM pg_class WHERE relname='%s';",asynch->GlobalVars->hydro_table);
		res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
		if(!PQntuples(res))
		{
			PQclear(res);
			sprintf(query,"CREATE TABLE %s(link_id int,time int,ratio real,discharge real); ALTER TABLE %s SET (autovacuum_enabled = false, toast.autovacuum_enabled = false);",asynch->GlobalVars->hydro_table,asynch->GlobalVars->hydro_table);
			res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
			CheckResError(res,"creating hydrographs table");
		}
		else
		{
			PQclear(res);
			sprintf(query,"TRUNCATE %s;",asynch->GlobalVars->hydro_table);
			res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
			CheckResError(res,"truncating hydrographs table");
		}
		PQclear(res);

		//Make sure the hydrograph tables are set correctly
		CheckPartitionedTable(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT],asynch->GlobalVars,Forecaster,num_tables,"archive_hydroforecast","forecast_time",schema);

		//Clear the future hydrographs in archive
		DeleteFutureValues(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT],num_tables,asynch->GlobalVars,"archive_hydroforecast",Forecaster->model_name,first_file,1,schema);

		//Disconnect from hydrograph database, connect to peakflow database
		DisconnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);
		ConnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_PEAK_OUTPUT]);

		//Clear all future peakflows
		sprintf(query,"TRUNCATE %s;",asynch->GlobalVars->peak_table);
		res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_PEAK_OUTPUT]->conn,query);
		CheckResError(res,"truncating peakforecast table");
		PQclear(res);

		//Disconnect
		DisconnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_PEAK_OUTPUT]);

		stop = time(NULL);
		printf("Total time to initialize tables: %.2f.\n",difftime(stop,start));
	}

	MPI_Barrier(MPI_COMM_WORLD);

	//Start the main loop
	while(!halt)
	{
		if(my_rank == 0)
		{
			time_t now = time(NULL);
			struct tm* now_info = localtime(&now);
			printf("\n\nPass %u\n",k);
			printf("Current time is %s",asctime(now_info));
		}

		//Clear buffers
		Flush_TransData(asynch->my_data);

		//Make some initializations
		first_file = last_file;
		last_file = last_file + (unsigned int) asynch->forcings[forecast_idx]->file_time * 60 * num_rainsteps;
		nextforcingtime = first_file + 60 * (unsigned int) rint(asynch->forcings[forecast_idx]->file_time) * (num_rainsteps-1);	//This is the actual timestamp of the last needed forcing data. This will be downloaded (unlike last_file)

		//Reset each link
		Asynch_Set_System_State(asynch,0.0,backup);
		Set_Output_User_forecastparams(asynch,first_file);
		Set_Output_PeakflowUser_Offset(asynch,first_file);
		Asynch_Write_Current_Step(asynch);
		Asynch_Set_Forcing_State(asynch,forecast_idx,0.0,first_file,last_file);

		for(i=0;i<asynch->GlobalVars->num_forcings;i++)	//Set any other database forcings to begin at first_file
		{
			if(asynch->forcings[i]->flag == 3)
				Asynch_Set_Forcing_State(asynch,i,0.0,first_file,asynch->forcings[i]->last_file);
		}

		//Check if a vacuum should be done
		if(my_rank == 0)	PerformTableMaintainance(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT],asynch->GlobalVars,Forecaster,&vac,hr1,num_tables,"archive_hydroforecast",schema);

		//Make sure all buffer flushing is done
		MPI_Barrier(MPI_COMM_WORLD);

		//Dump data for debugging and recovery
		if(k % 96 == 0)
		{
			sprintf(filename,"%s%u.rec",dump_filename,first_file);
			Asynch_Set_Snapshot_Output_Name(asynch,filename);
			Asynch_Take_System_Snapshot(asynch,NULL);
		}

		//Find the next time where rainfall occurs
		do
		{
			if(my_rank == 0)
			{
				ConnectPGDB(Forecaster->rainmaps_db);

				//Find the next rainfall time
				time(&start);
				sprintf(query,Forecaster->rainmaps_db->queries[0],nextforcingtime);;
				res = PQexec(Forecaster->rainmaps_db->conn,query);
				CheckResError(res,"checking for new rainfall data");
				time(&stop);
				printf("Total time to check for new rainfall data: %f.\n",difftime(stop,start));
				isnull = PQgetisnull(res,0,0);

				PQclear(res);
				DisconnectPGDB(Forecaster->rainmaps_db);
			}
			MPI_Bcast(&isnull,1,MPI_INT,0,MPI_COMM_WORLD);

			if(isnull)
			{
				if(my_rank == 0)
				{
					printf("No rainfall values returned from SQL database for forcing %u. %u %u\n",forecast_idx,last_file,isnull);
					PerformTableMaintainance(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT],asynch->GlobalVars,Forecaster,&vac,hr1,num_tables,"archive_hydroforecast",schema);
				}

				halt = CheckFinished(Forecaster->halt_filename);
				if(halt)
				{
					sprintf(filename,"%s%u.rec",dump_filename,first_file);
					Asynch_Set_Snapshot_Output_Name(asynch,filename);
					Asynch_Take_System_Snapshot(asynch,NULL);
				}
				else
				{
					fflush(stdout);
					ASYNCH_SLEEP(wait_time);
				}
			}
		} while(isnull && !halt);

		if(halt)	break;

		//Read in next set of rainfall data

		//Initialize some data for the first phase of calculations
		Asynch_Set_Total_Simulation_Time(asynch,simulation_time_with_data);		// !!!! This may not work for multiple forcings for forecasting. How do you handle different time resolutions? !!!!
		current_offset = first_file;
		Set_Output_User_forecastparams(asynch,current_offset);
		Set_Output_PeakflowUser_Offset(asynch,current_offset);

		MPI_Barrier(MPI_COMM_WORLD);
		time(&start);
if(my_rank == 0)
printf("first: %u last: %u\n",first_file,last_file);

		Asynch_Advance(asynch,1);

		MPI_Barrier(MPI_COMM_WORLD);
		if(my_rank == 0)
		{
			time(&stop);
			printf("Time for first phase calculations: %.2f\n",difftime(stop,start));
		}

		//Flush communication buffers
		Flush_TransData(asynch->my_data);

		//Reset the links (mostly) and make a backup for the second phase
		for(i=0;i<N;i++)	//Set time to 0.0
		{
			current = asynch->sys[i];
			if(current->list->head != NULL)
			{
				while(current->current_iterations > 1)
				{
					Remove_Head_Node(current->list);
					(current->current_iterations)--;
				}

				current->steps_on_diff_proc = 1;
				current->iters_removed = 0;
				current->rejected = 0;
				if(current->numparents == 0)	current->ready = 1;
				else				current->ready = 0;
				v_copy(current->list->head->y_approx,backup[i]);
			}
		}

		//Make second phase calculations
		MPI_Barrier(MPI_COMM_WORLD);
		time(&start);
		Asynch_Set_Total_Simulation_Time(asynch,forecast_time);
		Asynch_Deactivate_Forcing(asynch,forecast_idx);
		Asynch_Advance(asynch,1);
		Asynch_Activate_Forcing(asynch,forecast_idx);
		MPI_Barrier(MPI_COMM_WORLD);
		if(my_rank == 0)
		{
			time(&stop);
			printf("Time for second phase calculations: %.2f\n",difftime(stop,start));
		}

		//Output some data
		if(my_rank == 0)
		{
			printf("[%i]: The answer at ID %i at time %.12f is\n",my_rank,asynch->sys[asynch->my_sys[0]]->ID,asynch->sys[asynch->my_sys[0]]->last_t);
			Print_Vector(asynch->sys[asynch->my_sys[0]]->list->tail->y_approx);
		}

		//Upload the peak data to the database **********************************************************************************************

		//Adjust the table peak
		MPI_Barrier(MPI_COMM_WORLD);
		start = time(NULL);

		repeat_for_errors = Asynch_Create_Peakflows_Output(asynch);
		while(repeat_for_errors > 0)
		{
			if(my_rank == 0)	printf("[%i]: Attempting resend of peakflow data.\n",my_rank);
			ASYNCH_SLEEP(5);
			repeat_for_errors = Asynch_Create_Peakflows_Output(asynch);
		}

		MPI_Barrier(MPI_COMM_WORLD);
		if(my_rank == 0)
		{
			stop = time(NULL);
			printf("[%i]: Total time to transfer peak flow data: %.2f\n",my_rank,difftime(stop,start));
		}

		//Upload the hydrographs to the database ********************************************************************************************
		MPI_Barrier(MPI_COMM_WORLD);
		start = time(NULL);

		//Adjust the table hydrographs
		if(my_rank == 0)
		{
			//Make sure database connection is still good
			ConnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);

			sprintf(query,"TRUNCATE %s;",asynch->GlobalVars->hydro_table);
			res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
			CheckResError(res,"deleting hydrographs");
			PQclear(res);

			DisconnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		repeat_for_errors = Asynch_Create_Output(asynch,NULL);
		while(repeat_for_errors > 0)
		{
			if(my_rank == 0)	printf("[%i]: Attempting resend of hydrographs data.\n",my_rank);
			ASYNCH_SLEEP(5);
			repeat_for_errors = Asynch_Create_Output(asynch,NULL);
		}

		//Call functions
		if(my_rank == 0)
		{
			//Connect to database
			ConnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);

			//Functions for displaying data on IFIS
			if(Forecaster->ifis_display)
			{
				//Stage
				repeat_for_errors = 1;
				while(repeat_for_errors)
				{
					repeat_for_errors = 0;
					sprintf(query,"SELECT get_stages_ifc01();");
					res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
					repeat_for_errors = repeat_for_errors || CheckResError(res,"calling stage function");
					PQclear(res);
					if(repeat_for_errors)
					{
						printf("[%i]: Attempting to call stage function again...\n",my_rank);
						ASYNCH_SLEEP(5);
						CheckConnConnection(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);
					}
				}
/*
				//Warnings
				repeat_for_errors = 1;
				while(repeat_for_errors)
				{
					repeat_for_errors = 0;
					sprintf(query,"SELECT update_warnings_%s();",Forecaster->model_name);
					res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
					repeat_for_errors = repeat_for_errors || CheckResError(res,"calling warnings function");
					PQclear(res);
					if(repeat_for_errors)
					{
						printf("[%i]: Attempting to call warning function again...\n",my_rank);
						ASYNCH_SLEEP(5);
						CheckConnConnection(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);
					}
				}
*/

				//Run php script
				do
				{
					repeat_for_errors = system("wget -O /dev/null http://ifisfe.its.uiowa.edu/ifc/php/acquisition/_new_get_ifc_forecast.php");
					if(repeat_for_errors == -1)
					{
						printf("[%i]: Attempting to launch php script again...\n",my_rank);
						ASYNCH_SLEEP(5);
					}
				} while(repeat_for_errors == -1);
			}

			//Stage archive
			repeat_for_errors = 1;
			while(repeat_for_errors)
			{
				repeat_for_errors = 0;
				sprintf(query,"ALTER TABLE master_archive_hydroforecast_%s ALTER COLUMN forecast_time SET DEFAULT %u;",Forecaster->model_name,current_offset);
				res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
				repeat_for_errors = repeat_for_errors || CheckResError(res,"setting default value");
				PQclear(res);

				sprintf(query,"SELECT copy_to_archive_hydroforecast_%s();",Forecaster->model_name);
				res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
				repeat_for_errors = repeat_for_errors || CheckResError(res,"calling stage archive function");
				PQclear(res);

				sprintf(query,"ALTER TABLE master_archive_hydroforecast_%s ALTER COLUMN forecast_time DROP DEFAULT;",Forecaster->model_name);
				res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
				repeat_for_errors = repeat_for_errors || CheckResError(res,"dropping default value");
				PQclear(res);

				if(repeat_for_errors)
				{
					printf("[%i]: Attempting to call stage archive function again...\n",my_rank);
					ASYNCH_SLEEP(5);
					CheckConnConnection(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);
				}
			}

			//Disconnect
			DisconnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);
		}

		if(my_rank == 0)
		{
			time(&stop);
			printf("[%i]: Total time to transfer hydrograph data: %.2f\n",my_rank,difftime(stop,start));
		}
		fflush(stdout);
		MPI_Barrier(MPI_COMM_WORLD);

		//Check if program has received a terminate signal **********************************************************************************
		k++;
		halt = CheckFinished(Forecaster->halt_filename);

		//If stopping, make a .rec file
		if(halt)
		{
			for(i=0;i<N;i++)
			{
				current = asynch->sys[i];
				if(current->list->head != NULL)
					v_copy(backup[i],current->list->tail->y_approx);
			}

			sprintf(filename,"%s%u.rec",dump_filename,first_file);
			Asynch_Set_Snapshot_Output_Name(asynch,filename);
			Asynch_Take_System_Snapshot(asynch,NULL);
		}
	}

	//Clean up *************************************************************************
	free(query);
	for(i=0;i<N;i++)	v_free(backup[i]);
	free(backup);
	Free_ForecastData(&Forecaster);
	Asynch_Delete_Temporary_Files(asynch);
	Free_Output_PeakflowUser_Offset(asynch);
	Free_Output_User_forecastparams(asynch);
	Asynch_Free(asynch);
	MPI_Finalize();
	return 0;
}



//Output functions ****************************************************************************
int Output_Linkid(double t,VEC* y_i,VEC* global_params,VEC* params,int state,void* user)
{
	CustomParams* forecastparams = (CustomParams*) user;
	return forecastparams->ID;
}

int Output_Timestamp(double t,VEC* y_i,VEC* global_params,VEC* params,int state,void* user)
{
	CustomParams* forecastparams = (CustomParams*) user;
	return (int)(round(t * 60.0 + forecastparams->offset) + 0.1);
}


//Custom parameters for forecasting ***********************************************************
void Init_Output_User_forecastparams(asynchsolver* asynch)
{
	unsigned int i,my_N = asynch->my_N,*my_sys = asynch->my_sys;
	Link** sys = asynch->sys;

	for(i=0;i<my_N;i++)
		sys[my_sys[i]]->output_user = malloc(sizeof(CustomParams));
}

void Free_Output_User_forecastparams(asynchsolver* asynch)
{
	unsigned int i,my_N = asynch->my_N,*my_sys = asynch->my_sys;
	Link** sys = asynch->sys;

	for(i=0;i<my_N;i++)
	{
		free(sys[my_sys[i]]->output_user);
		sys[my_sys[i]]->output_user = NULL;
	}
}

void Set_Output_User_forecastparams(asynchsolver* asynch,unsigned int offset)
{
	unsigned int i,my_N = asynch->my_N,*my_sys = asynch->my_sys;
	Link** sys = asynch->sys;
	CustomParams* forecastparams;

	for(i=0;i<my_N;i++)
	{
		forecastparams = (CustomParams*) sys[my_sys[i]]->output_user;
		forecastparams->ID = sys[my_sys[i]]->ID;
		forecastparams->offset = offset;
	}
}

void Init_Output_PeakflowUser_Offset(asynchsolver* asynch)
{
	unsigned int i,my_N = asynch->my_N,*my_sys = asynch->my_sys;
	Link** sys = asynch->sys;

	for(i=0;i<my_N;i++)
		sys[my_sys[i]]->peakoutput_user = malloc(sizeof(unsigned int));
}

void Free_Output_PeakflowUser_Offset(asynchsolver* asynch)
{
	unsigned int i,my_N = asynch->my_N,*my_sys = asynch->my_sys;
	Link** sys = asynch->sys;

	for(i=0;i<my_N;i++)
	{
		free(sys[my_sys[i]]->peakoutput_user);
		sys[my_sys[i]]->peakoutput_user = NULL;
	}
}

void Set_Output_PeakflowUser_Offset(asynchsolver* asynch,unsigned int offset)
{
	unsigned int i,my_N = asynch->my_N,*my_sys = asynch->my_sys;
	Link** sys = asynch->sys;

	for(i=0;i<my_N;i++)
		*(unsigned int*)(sys[my_sys[i]]->peakoutput_user) = offset;
}


