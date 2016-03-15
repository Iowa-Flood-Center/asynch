#include "forecaster_methods.h"


//Deletes all future values from a set of partitioned tables.
int DeleteFutureValues(ConnData* conninfo,unsigned int num_tables,UnivVars* GlobalVars,char* table_name,char* model_name,unsigned int clear_after,unsigned int equal,char* schema)
{
	if(!conninfo)	return 1;

	PGresult* res;
	int i,del_table,error = 0,last_table;
	time_t current_time,clear_time = (time_t) clear_after;
	char* query = (char*) malloc(GlobalVars->query_size*sizeof(char));
	char operation[3];
	if(equal)	sprintf(operation,">=");
	else		sprintf(operation,">");
	ConnectPGDB(conninfo);

	//Find the last table with values to destroy
	time(&current_time);
	del_table = (int) (difftime(current_time,clear_time) / (60.0*60.0*24.0) + 1e-6);

	//Truncate any tables with index lower than del_table
	last_table = (num_tables < del_table) ? num_tables : del_table;
	for(i=0;i<last_table;i++)
	{
		sprintf(query,"TRUNCATE %s%s_%s_%i;",schema,table_name,model_name,i);
		res = PQexec(conninfo->conn,query);
		error = CheckResError(res,"truncating table");
		PQclear(res);
		if(error)	return error;
	}

	//Delete del_table
	if(del_table < num_tables)
	{
		sprintf(query,"DELETE FROM %s%s_%s_%i WHERE forecast_time %s %u;",schema,table_name,model_name,i,operation,clear_after);
		res = PQexec(conninfo->conn,query);
		error = CheckResError(res,"deleting from table");
		PQclear(res);
		if(error)	return error;
	}

	//Clean up
	DisconnectPGDB(conninfo);
	free(query);
	return error;
}


//Checks if the time is right to perform maintainance on the database.
void PerformTableMaintainance(ConnData* conninfo_hydros,UnivVars* GlobalVars,ForecastData* Forecaster,short int* vac,short unsigned int hr1,unsigned int num_tables,char* tablename,char* schema)
{
	if(!conninfo_hydros)	return;

	int j;
	time_t start,stop;
	PGresult* res;
	char query[ASYNCH_MAX_QUERY_LENGTH];
	struct tm* timeinfo;
	time(&start);
	timeinfo = localtime(&start);

	if(timeinfo->tm_hour == hr1 && *vac == 0)
	{
		printf("[%i]: Performing maintainance. Current time is %s",my_rank,asctime(timeinfo));

		//Adjust partitioned hydroforecast tables
		ConnectPGDB(conninfo_hydros);
		sprintf(query,"DROP TABLE %s%s_%s_%u;",schema,tablename,Forecaster->model_name,num_tables-1);
		res = PQexec(conninfo_hydros->conn,query);
		CheckResError(res,"dropping end archive table");
		PQclear(res);

		for(j=num_tables-2;j>=0;j--)
		{
			sprintf(query,"ALTER TABLE %s%s_%s_%i RENAME TO %s_%s_%i;",schema,tablename,Forecaster->model_name,j,tablename,Forecaster->model_name,j+1);
			res = PQexec(conninfo_hydros->conn,query);
			CheckResError(res,"renaming table");
			PQclear(res);

			sprintf(query,"ALTER INDEX idx_%s_%s_%i_forecast_time_link_id RENAME TO idx_%s_%s_%i_forecast_time_link_id;",tablename,Forecaster->model_name,j,tablename,Forecaster->model_name,j+1);
			res = PQexec(conninfo_hydros->conn,query);
			CheckResError(res,"renaming index on archive table");
			PQclear(res);
		}

		sprintf(query,"CREATE TABLE %s%s_%s_0 ( ) INHERITS (%smaster_%s_%s);",schema,tablename,Forecaster->model_name,schema,tablename,Forecaster->model_name);
		res = PQexec(conninfo_hydros->conn,query);
		CheckResError(res,"creating archive table 0");
		PQclear(res);

		sprintf(query,"CREATE INDEX idx_%s_%s_0_forecast_time_link_id ON %s%s_%s_0 USING btree (forecast_time,link_id);",tablename,Forecaster->model_name,schema,tablename,Forecaster->model_name);
		res = PQexec(conninfo_hydros->conn,query);
		CheckResError(res,"creating index on archive table 0");
		PQclear(res);

		DisconnectPGDB(conninfo_hydros);

		//Set flag and print the total time
		*vac = 1;
		time(&stop);
		printf("[%i]: Database cleanup complete. Total time %.2f.\n\n",my_rank,difftime(stop,start));
	}
	else	if(timeinfo->tm_hour != hr1)	*vac = 0;
}

//Checks that the timestamps of the archive tables match up correctly with the trigger.
//If not, the tables are adjusted.
//Note: The total number of tables is hard coded below.
//Checks that the timestamps of a partitioned table match up correctly with the trigger.
//If not, the child tables are adjusted.
//void CheckPeakforecastTable(ConnData* conninfo,UnivVars* GlobalVars,ForecastData* Forecaster,unsigned int num_tables)
void CheckPartitionedTable(ConnData* conninfo,UnivVars* GlobalVars,ForecastData* Forecaster,unsigned int num_tables,char* tablename,char* colname,char* schema)
{
	if(!conninfo)	return;

	unsigned int i,current_time,table_time,table_index,correct_table_index,diff_table_index;
	int diff_time,j,last_table_index;
	PGresult *res;
	char* query = conninfo->query;

	//Connect to db
	ConnectPGDB(conninfo);

	//Grab the current time from the database
	sprintf(query,"SELECT EXTRACT('epoch' FROM current_date AT time zone 'UTC');");
	res = PQexec(conninfo->conn,query);
	current_time = (unsigned int) rint(atof(PQgetvalue(res,0,0)));
	PQclear(res);

	//Find the first table with something in it
	table_index = num_tables;
	for(i=0;i<num_tables;i++)
	{
		sprintf(query,"SELECT %s FROM %s%s_%s_%u LIMIT 1;",colname,schema,tablename,Forecaster->model_name,i);
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"checking table contents");
		if(PQntuples(res))
		{
			table_time = (unsigned int) atoi(PQgetvalue(res,0,0));
			diff_time = table_time - current_time;
			table_index = i;
			PQclear(res);
			break;
		}
		PQclear(res);
	}

	//If no data exists in any table, then all is well
	if(table_index == num_tables)
		goto finish_up;

	//If table_time is in the correct table, return
	//if(-86400*(int)table_index < (int) diff_time && (int) diff_time <= -86400*((int)table_index-1))	return;
	if(-86400*(int)table_index <= (int) diff_time && (int) diff_time < -86400*((int)table_index-1))
		goto finish_up;

	//Otherwise, find the correct table
	correct_table_index = num_tables;
	for(i=table_index+1;i<num_tables;i++)
	{
		//if(-86400*(int)i < (int)diff_time && (int)diff_time <= -86400*((int)i-1))
		if(-86400*(int)i <= (int)diff_time && (int)diff_time < -86400*((int)i-1))
		{
			correct_table_index = i;
			break;
		}
	}
	diff_table_index = correct_table_index - table_index;

	//Trash the tables at the end
	last_table_index = num_tables-diff_table_index - 1;
	for(j=num_tables-1;j>last_table_index;j--)
	{
		sprintf(query,"DROP TABLE %s%s_%s_%u;",schema,tablename,Forecaster->model_name,j);
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"dropping table");
		PQclear(res);
	}

	//Move all the tables
	for(j=last_table_index;j>=0;j--)
	{
		sprintf(query,"ALTER TABLE %s%s_%s_%i RENAME TO %s_%s_%i;",schema,tablename,Forecaster->model_name,j,tablename,Forecaster->model_name,j+diff_table_index);
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"renaming table");
		PQclear(res);

		sprintf(query,"ALTER INDEX idx_%s_%s_%i_forecast_time_link_id RENAME TO idx_%s_%s_%i_forecast_time_link_id;",tablename,Forecaster->model_name,j,tablename,Forecaster->model_name,j+diff_table_index);
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"renaming index on archive table");
		PQclear(res);
	}

	//Create new tables
	for(i=0;i<diff_table_index;i++)
	{
		sprintf(query,"CREATE TABLE %s%s_%s_%u ( ) INHERITS (%smaster_%s_%s);",schema,tablename,Forecaster->model_name,i,schema,tablename,Forecaster->model_name);
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"creating table");
		PQclear(res);

		sprintf(query,"CREATE INDEX idx_%s_%s_%u_forecast_time_link_id ON %s%s_%s_%u USING btree (forecast_time,link_id);",tablename,Forecaster->model_name,i,schema,tablename,Forecaster->model_name,i);
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"creating index on archive table");
		PQclear(res);
	}

	finish_up:
	DisconnectPGDB(conninfo);
}

//Creates the halt file and sets the value to 0
void CreateHaltFile(char* filename)
{
	FILE* outputfile;

	if(my_rank == 0)
	{
		outputfile = fopen(filename,"w");
		if(!outputfile)	printf("Warning: Could not create halt file %s.\n",filename);
		fprintf(outputfile,"0");
		fclose(outputfile);
	}

	MPI_Barrier(MPI_COMM_WORLD);
}

//Returns 0 if the program should continue, 1 if the program should terminate.
short int CheckFinished(char* filename)
{
	FILE* inputfile;
	short int halt;
	time_t now;
	struct tm* timeinfo;

	if(my_rank == 0)
	{
		//Check terminate file
		inputfile = fopen(filename,"r");
		if(!inputfile)	halt = 0;
		else		fscanf(inputfile,"%hu",&halt);

		//Check if signal received
		if(halt)
		{
			time(&now);
			timeinfo = localtime(&now);
			printf("\nReceived halt signal on %s",asctime(timeinfo));
		}

		fclose(inputfile);
	}

	//Notify all processes of how to proceed
	MPI_Bcast(&halt,1,MPI_SHORT,0,MPI_COMM_WORLD);

	return halt;
}

//This function holds a proc until it's safe to upload data to a database.
//This is used to prevent multiple forecasters from choking the database.
//Returns 0 if the proc is safe to upload, 1 if an error occurred.
int WaitForDB(ConnData* conninfo,unsigned int naptime,int stale_time,unsigned int query_size)
{
	char tablename[] = "forecaster_upload";
	char query[ASYNCH_MAX_QUERY_LENGTH];
	PGresult *res;
	int error = 0;
	int flag,time,now,wait;
	short int status;

	if(my_rank == 0)
	{
		//Connect to db
		ConnectPGDB(conninfo);

		//Check that the table exists
		sprintf(query,"CREATE TABLE IF NOT EXISTS %s(flag integer, set_time integer);",tablename);
		res = PQexec(conninfo->conn,query);
		if(CheckResError(res,"creating forecaster upload table"))	error = 1;
		PQclear(res);

		do
		{
			wait = 0;

			//See what's in the table
			sprintf(query,"SELECT flag,set_time,EXTRACT('epoch' FROM current_timestamp)::integer FROM %s;",tablename);
			res = PQexec(conninfo->conn,query);
			if(CheckResError(res,"checking forecaster upload table"))	error = 1;

			if(PQntuples(res) == 0)
			{
				PQclear(res);
				sprintf(query,"INSERT INTO %s VALUES (0,EXTRACT('epoch' FROM current_timestamp)::integer); SELECT flag,set_time,EXTRACT('epoch' FROM current_timestamp)::integer FROM %s;",tablename,tablename);
				res = PQexec(conninfo->conn,query);
				if(CheckResError(res,"loading values into forecaster upload table"))	error = 1;
			}

			//Check if time is stale
			flag = atoi(PQgetvalue(res,0,0));
			time = atoi(PQgetvalue(res,0,1));
			now = atoi(PQgetvalue(res,0,2));
			PQclear(res);

			if(now-time > stale_time || flag == 0)	//Timestamp is really old or no one is uploading. So try to go.
			{
				sprintf(query,"LOCK TABLE %s IN ACCESS EXCLUSIVE MODE NOWAIT; TRUNCATE TABLE %s; INSERT INTO %s VALUES (1,EXTRACT('epoch' FROM current_timestamp)::integer);",tablename,tablename,tablename);
				res = PQexec(conninfo->conn,query);

				status = PQresultStatus(res);
				if( !(status == PGRES_COMMAND_OK || status == PGRES_TUPLES_OK) )	//Another forecaster beat you to the lock (or maybe an error occurred). Go back.
				{
					printf("[%i]: Someone else got the lock. Trying again...\n",my_rank);
					wait = 1;
					ASYNCH_SLEEP(naptime);
				}
				PQclear(res);
			}
			else	//Another forecaster is actually uploading. Go back.
			{
				wait = 1;
				ASYNCH_SLEEP(naptime);
			}
		} while(wait && !error);

		//Disconnect from db
		DisconnectPGDB(conninfo);
	}

	MPI_Bcast(&error,1,MPI_INT,0,MPI_COMM_WORLD);
	return error;
}

//Releases the lock set by the routine WaitForDB.
void FreeDBLock(ConnData* conninfo,unsigned int query_size)
{
	char tablename[] = "forecaster_upload";
	char query[ASYNCH_MAX_QUERY_LENGTH];
	PGresult *res;
	int error;

	if(my_rank == 0)
	{
		//Connect to db
		ConnectPGDB(conninfo);

		sprintf(query,"LOCK TABLE %s IN ACCESS EXCLUSIVE MODE NOWAIT; TRUNCATE TABLE %s; INSERT INTO %s VALUES (0,EXTRACT('epoch' FROM current_timestamp)::integer);",tablename,tablename,tablename);
		do
		{
			res = PQexec(conninfo->conn,query);
			error = CheckResError(res,"unlocking forecaster upload table");
			PQclear(res);
			if(error)
			{
				ASYNCH_SLEEP(2*60);
				CheckConnConnection(conninfo);
			}
		} while(error);

		//Disconnect from db
		DisconnectPGDB(conninfo);
	}

	MPI_Barrier(MPI_COMM_WORLD);
}


ForecastData* Init_ForecastData(char* fcst_filename,unsigned int string_size)
{
	FILE* inputfile = NULL;
	ForecastData* Forecaster;
	int errorcode,valsread;
	char end_char;
	unsigned int buff_size = string_size + 20;
	char* linebuffer = (char*) malloc(buff_size*sizeof(char));
	MPI_Barrier(MPI_COMM_WORLD);

	if(my_rank == 0)
	{
		//Open file
		inputfile = fopen(fcst_filename,"r");
		errorcode = 0;
		if(!inputfile)
		{
			printf("[%i]: Error opening file %s.\n",my_rank,fcst_filename);
			errorcode = 1;
		}
	}

	//Check if forecast file was openned
	MPI_Bcast(&errorcode,1,MPI_INT,0,MPI_COMM_WORLD);
	if(errorcode)	return NULL;

	//Reserve space
	Forecaster = (ForecastData*) malloc(sizeof(ForecastData));
	Forecaster->model_name = (char*) malloc(string_size*sizeof(char));

	//Read table name
	//if(my_rank == 0)
	{
		ReadLineFromTextFile(inputfile,linebuffer,buff_size,string_size);
		valsread = sscanf(linebuffer,"%s",Forecaster->model_name);
		if(ReadLineError(valsread,1,"forecaster model name"))	return NULL;
		//length = strlen(Forecaster->model_name);
	}
	//MPI_Bcast(&length,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	//MPI_Bcast(Forecaster->model_name,length+1,MPI_CHAR,0,MPI_COMM_WORLD);

	//Read if data is displayed on ifis
	//if(my_rank == 0)
	{
		ReadLineFromTextFile(inputfile,linebuffer,buff_size,string_size);
		valsread = sscanf(linebuffer,"%hi",&(Forecaster->ifis_display));
		if(ReadLineError(valsread,1,"flag if displaying on ifis"))	return NULL;
	}
	//MPI_Bcast(&(Forecaster->ifis_display),1,MPI_SHORT,0,MPI_COMM_WORLD);

	//Read which forcing index is used for forecasting
	//if(my_rank == 0)
	{
		ReadLineFromTextFile(inputfile,linebuffer,buff_size,string_size);
		valsread = sscanf(linebuffer,"%u",&(Forecaster->forecasting_forcing));
		if(ReadLineError(valsread,1,"index of forecastin forcing"))	return NULL;
	}
	//MPI_Bcast(&(Forecaster->forecasting_forcing),1,MPI_UNSIGNED,0,MPI_COMM_WORLD);

	//Read number of rainfall steps to use per forecast
	//if(my_rank == 0)
	{
		ReadLineFromTextFile(inputfile,linebuffer,buff_size,string_size);
		valsread = sscanf(linebuffer,"%u",&(Forecaster->num_rainsteps));
		if(ReadLineError(valsread,1,"number of precipitation values"))	return NULL;
	}
	//MPI_Bcast(&(Forecaster->num_rainsteps),1,MPI_UNSIGNED,0,MPI_COMM_WORLD);

	//Read forecast window
	ReadLineFromTextFile(inputfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%lf",&(Forecaster->forecast_window));
	if(ReadLineError(valsread,1,"forecast window"))	return NULL;

	//Read and create a database connection for the rain maps
	Forecaster->rainmaps_filename = NULL;
	Forecaster->rainmaps_db = NULL;
	//if(my_rank == 0)
	{
		Forecaster->rainmaps_filename = (char*) malloc(string_size*sizeof(char));
		ReadLineFromTextFile(inputfile,linebuffer,buff_size,string_size);
		valsread = sscanf(linebuffer,"%s",Forecaster->rainmaps_filename);
		if(ReadLineError(valsread,1,"rain map filename"))	return NULL;

		Forecaster->rainmaps_db = ReadDBC(Forecaster->rainmaps_filename,string_size);
		if(!Forecaster->rainmaps_db)	return NULL;
	}

	//Read halt filename
	Forecaster->halt_filename = (char*) malloc(string_size*sizeof(char));
	//if(my_rank == 0)
	{
		ReadLineFromTextFile(inputfile,linebuffer,buff_size,string_size);
		valsread = sscanf(linebuffer,"%s",Forecaster->halt_filename);
		if(ReadLineError(valsread,1,"halt filename"))	return NULL;
		//length = strlen(Forecaster->halt_filename);
	}
	//MPI_Bcast(&length,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	//MPI_Bcast(Forecaster->halt_filename,length+1,MPI_CHAR,0,MPI_COMM_WORLD);

	//Read ending mark
	//if(my_rank == 0)
	{
		ReadLineFromTextFile(inputfile,linebuffer,buff_size,string_size);
		valsread = sscanf(linebuffer,"%c",&end_char);
		if(ReadLineError(valsread,1,"ending mark"))	return NULL;
	}
	//MPI_Bcast(&end_char,1,MPI_CHAR,0,MPI_COMM_WORLD);

	//Clean up
	free(linebuffer);
	if(my_rank == 0)	fclose(inputfile);

	MPI_Barrier(MPI_COMM_WORLD);
	if(end_char != '#')
	{
		if(my_rank == 0)
			printf("[%i]: Error: Ending mark not seen in %s.\n",my_rank,fcst_filename);
		return NULL;
	}
	return Forecaster;
}

void Free_ForecastData(ForecastData** Forecaster)
{
	if((*Forecaster)->rainmaps_filename)
	{
		free((*Forecaster)->rainmaps_filename);
		ConnData_Free((*Forecaster)->rainmaps_db);
	}
	free((*Forecaster)->model_name);
	free((*Forecaster)->halt_filename);
	free(*Forecaster);
	*Forecaster = NULL;
}


//Uses ssh to transfer a data file.
//I took most of this code from the example file scp_write.c in the libssh2 examples.
//!!!! This could really be generalized... !!!!
int SendFilesTo51(char* loclfile,char* serverlocation)
{
    int error_code = 0;

#if defined(_MSV_VER) && defined(WITH_LIBSSH2)
	char filename[1024];	
	if(FindFilename(loclfile,filename))
	{
		printf("Error: Bad filename for hydrograph file. (%s)\n",loclfile);
		return 1;
	}

	char scppath[1024];
	sprintf(scppath,"%s/%s",serverlocation,filename);
	//sprintf(scppath,"/data/ifc_01_hydro/%s",filename);

	//Uh, yeah, probably not very secure...
	char *username = "my_user_name",*password = "my_password";

	//Check out the file info
	struct stat fileinfo;
	stat(loclfile,&fileinfo);
	FILE* local = fopen(loclfile, "rb");
	if(!local)
	{
		printf("Can't open local file %s\n", loclfile);
		return -1;
	}

	//Init ssh stuff
	long unsigned int hostaddr = inet_addr("128.255.26.166");
	if(hostaddr < 0)
	{
		printf("Bad host address.\n");
		return 1;
	}
	int rc = libssh2_init(0);
	if(rc)
	{
		printf("Problem initializing libssh2 (%i)\n",rc);
		return 1;
	}
	int sock = socket(AF_INET,SOCK_STREAM,0);
	if(-1 == sock)
	{
		printf("Failed to create socket\n");
		return 1;
	}

	struct sockaddr_in sin;
	sin.sin_family = AF_INET;
	sin.sin_port = htons(22);
	sin.sin_addr.s_addr = hostaddr;
	if(connect(sock,(struct sockaddr*)(&sin),sizeof(struct sockaddr_in)) != 0)
	{
		printf("Failed to connect!\n");
		return 1;
	}

	//Create session instance and start it
	LIBSSH2_SESSION *session;
	session = libssh2_session_init();
	if(!session)	return 1;
	rc = libssh2_session_handshake(session, sock);
	if(rc)
	{
		printf("Failure establishing SSH session %i\n",rc);
		return 1;
	}

	//Authenticate
	const char* fingerprint = libssh2_hostkey_hash(session,LIBSSH2_HOSTKEY_HASH_SHA1);
	if(libssh2_userauth_password(session, username, password))
	{
		printf("Authentication by password failed.\n");
		error_code = 1;
		goto shutdown;
        }

	//Copy file
	LIBSSH2_CHANNEL *channel = libssh2_scp_send(session,scppath,fileinfo.st_mode & 0777,(unsigned long)fileinfo.st_size);
	if(!channel)
	{
		char *errmsg;
		int errlen;
		int err = libssh2_session_last_error(session, &errmsg, &errlen, 0);

		printf("Unable to open a session: (%d) %s\n", err, errmsg);
		error_code = 1;
		goto shutdown;
	}

	char mem[1024],*ptr;
	size_t nread;
	do
	{
		nread = fread(mem, 1, sizeof(mem), local);
		if (nread <= 0)
		{
			// end of file
			break;
		}
		ptr = mem;

		do
		{
			// write the same data over and over, until error or completion 
			rc = libssh2_channel_write(channel, ptr, nread);

			if (rc < 0)
			{
				printf("ERROR %d\n", rc);
				break;
			}
			else
			{
				// rc indicates how many bytes were written this time
				ptr += rc;
				nread -= rc;
			}
		} while(nread);
	} while(1);

	//Tell server we are done
	libssh2_channel_send_eof(channel);
	libssh2_channel_wait_eof(channel);
	libssh2_channel_wait_closed(channel);

	//Clean up
	libssh2_channel_free(channel);
	channel = NULL;

	shutdown:
	if(session)
	{
	 	libssh2_session_disconnect(session, "Normal Shutdown, Thank you for playing");
		libssh2_session_free(session);
	}
	close(sock);
	if(local)	fclose(local);
	if(!error_code && remove(loclfile))
		printf("[%i]: Error deleting file %s.\n",my_rank,loclfile);
	libssh2_exit();

	if(!error_code)	printf("File %s uploaded!\n",loclfile);
#endif 
	return error_code;
}


