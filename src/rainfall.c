#include "rainfall.h"



//This reads in a set of binary files for the rainfall at each link.
//Assumes the file is full of floats. Assumes no IDs are in the file and that IDs are consecutive starting from 0
//Link** sys: An array of links.
//int N: The number of links in sys.
//int my_N: The number of links assigned to this process.
//UnivVars* GlobalVars: Contains all the information that is shared by every link in the system.
//int* my_sys: Array of links assigned to this process (value is location in sys array).
//int* assignments (set by this method): Will be an array with N entries. assignments[i] will the process link sys[i] is assigned to.
//char strfilename[]: String of the filename for the rain file to read (NOT .str files). Filenames should be indexed.
//unsigned int first: The index of the file to read first.
//unsigned int last: The index of the file to read last.
//double t_0: The time at which the first file starts.
//double increment: The amount of time between consecutively indexed files.
//int** id_to_loc (set by this method): Will be an array with N rows and 2 columns, sorted by first col. First col is a link id and second is
//				the location of the id in sys.
//unsigned int max_files: The maximum number of files to be read.
int Create_Rain_Data_Par(Link** sys,unsigned int N,unsigned int my_N,UnivVars* GlobalVars,unsigned int* my_sys,int* assignments,char strfilename[],unsigned int first,unsigned int last,double t_0,double increment,Forcing* forcing,unsigned int** id_to_loc,unsigned int max_files,unsigned int forcing_idx)
{
	unsigned int i,j,curr_idx;
	unsigned int k;
	Link* current;
	float forcing_buffer;
	unsigned int holder;
	char filename[128];
	FILE* stormdata = NULL;
	unsigned int numfiles = last - first + 1;

	//This is a time larger than any time in which the integrator is expected to get
	double ceil_time = 1e300;
	if(sys[my_sys[0]]->last_t > ceil_time*0.1)
		printf("[%i]: Warning: integrator time is extremely large (about %e). Loss of precision may occur.\n",my_rank,sys[my_sys[i]]->last_t);

	//Check that space for rain data has been allocated.
	if(sys[my_sys[0]]->forcing_buff[forcing_idx] == NULL)
	{
		for(i=0;i<my_N;i++)
		{
			sys[my_sys[i]]->forcing_buff[forcing_idx] = (ForcingData*) malloc(sizeof(ForcingData));
			sys[my_sys[i]]->forcing_buff[forcing_idx]->rainfall = (double**) malloc((max_files + 1)*sizeof(double*));
			sys[my_sys[i]]->forcing_buff[forcing_idx]->n_times = numfiles + 1;
			for(j=0;j<max_files+1;j++)	sys[my_sys[i]]->forcing_buff[forcing_idx]->rainfall[j] = (double*) malloc(2*sizeof(double));
		}
	}

	//Read through the files.
	for(k=0;k<numfiles;k++)
	{
		//sprintf(filename,"%s%i",strfilename,first+k);		//This is the usual one to use
		sprintf(filename,"%s%i",strfilename,first+k);
		//sprintf(filename,"%srain%i",strfilename,first+k);
		//sprintf(filename,"%sfile-%i",strfilename,first+k);
		stormdata = fopen(filename,"r");
		if(stormdata == NULL)	printf("[%i]: Error opening file %s\n",my_rank,filename);

		for(i=0;i<N;i++)
		{
			//Find the location of ID i (should be i by assumption)
			curr_idx = id_to_loc[i][1];

			//Read the index
			if(assignments[curr_idx] == my_rank)	//If information about the link is needed on this process
			{
				//Read in the storm data for this link
				fread(&forcing_buffer,sizeof(float),1,stormdata);

				//This assumes different endianness
				holder = *(unsigned int*) &forcing_buffer;	//Pointers are fun
				holder = (((holder & 0x0000ffff)<<16) | ((holder & 0xffff0000)>>16));
				holder = (((holder & 0x00ff00ff)<<8) | ((holder & 0xff00ff00)>>8));
				forcing_buffer = *(float*) &holder;

				//Store the data
				sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k][0] = t_0 + k*increment;
				sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k][1] = forcing_buffer;
			}
			else		//This link data is not needed on this process.
				fread(&forcing_buffer,sizeof(float),1,stormdata);
		}

		fclose(stormdata);
	}

	if(my_rank == 0)
		printf("Read %i binary files.\n",numfiles);

	//Add in terms for no rainfall if max_files > numfiles
	for(i=0;i<my_N;i++)
	{
		curr_idx = sys[my_sys[i]]->location;
		for(j=numfiles;j<max_files;j++)
		{
			sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[j][0] = sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[j-1][0] + .0001;
			sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[j][1] = 0.0;
		}
	}

	//Add a ceiling term
	for(i=0;i<my_N;i++)
	{
		curr_idx = sys[my_sys[i]]->location;
		//sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[max_files][0] = GlobalVars->maxtime + 1.0;
		sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[max_files][0] = ceil_time;
		sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[max_files][1] = -1.0;
	}

	//Calculate the first rain change time and set rain_value
	for(i=0;i<my_N;i++)
	{
		current = sys[my_sys[i]];
		forcing_buffer = current->forcing_buff[forcing_idx]->rainfall[0][1];
		current->forcing_values[forcing_idx] = forcing_buffer;
		current->forcing_indices[forcing_idx] = 0;

		for(j=1;j<current->forcing_buff[forcing_idx]->n_times;j++)
		{
			if( fabs(forcing_buffer - current->forcing_buff[forcing_idx]->rainfall[j][1]) > 1e-14 )
			{
				current->forcing_change_times[forcing_idx] = current->forcing_buff[forcing_idx]->rainfall[j][0];
				break;
			}
		}
		if(j == current->forcing_buff[forcing_idx]->n_times)
		{
			current->forcing_change_times[forcing_idx] = current->forcing_buff[forcing_idx]->rainfall[j-1][0];
		}
	}

	return 0;
}

//This reads in a set of gzip compressed binary files for the rainfall at each link.
//Assumes the file is full of floats. Assumes no IDs are in the file and that IDs are consecutive starting from 0
//Link** sys: An array of links.
//int N: The number of links in sys.
//int my_N: The number of links assigned to this process.
//UnivVars* GlobalVars: Contains all the information that is shared by every link in the system.
//int* my_sys: Array of links assigned to this process (value is location in sys array).
//int* assignments (set by this method): Will be an array with N entries. assignments[i] will the process link sys[i] is assigned to.
//char strfilename[]: String of the filename for the rain file to read (NOT .str files). Filenames should be indexed.
//unsigned int first: The index of the file to read first.
//unsigned int last: The index of the file to read last.
//double t_0: The time at which the first file starts.
//double increment: The amount of time between consecutively indexed files.
//int** id_to_loc (set by this method): Will be an array with N rows and 2 columns, sorted by first col. First col is a link id and second is
//				the location of the id in sys.
//unsigned int max_files: The maximum number of files to be read.
int Create_Rain_Data_GZ(Link** sys,unsigned int N,unsigned int my_N,UnivVars* GlobalVars,unsigned int* my_sys,int* assignments,char strfilename[],unsigned int first,unsigned int last,double t_0,double increment,Forcing* forcing,unsigned int** id_to_loc,unsigned int max_files,unsigned int forcing_idx)
{
	unsigned int i,j,curr_idx;
	unsigned int k;
	Link* current;
	float rainfall_buffer;
	unsigned int holder;
	char filename[128];
	FILE* stormdata = NULL;
	unsigned int numfiles = last - first + 1;
	FILE* compfile = NULL;
	unsigned int buffer_size = sizeof(unsigned int)+sizeof(float);
	char transferbuffer[buffer_size];

	//This is a time larger than any time in which the integrator is expected to get
	double ceil_time = 1e300;
	if(sys[my_sys[0]]->last_t > ceil_time*0.1)
		printf("[%i]: Warning: integrator time is extremely large (about %e). Loss of precision may occur.\n",my_rank,sys[my_sys[i]]->last_t);

	//Check that space for rain data has been allocated.
	if(sys[my_sys[0]]->forcing_buff[forcing_idx] == NULL)
	{
		for(i=0;i<my_N;i++)
		{
			sys[my_sys[i]]->forcing_buff[forcing_idx] = (ForcingData*) malloc(sizeof(ForcingData));
			sys[my_sys[i]]->forcing_buff[forcing_idx]->rainfall = (double**) malloc((max_files + 1)*sizeof(double*));
			sys[my_sys[i]]->forcing_buff[forcing_idx]->n_times = numfiles + 1;
			for(j=0;j<max_files+1;j++)	sys[my_sys[i]]->forcing_buff[forcing_idx]->rainfall[j] = (double*) malloc(2*sizeof(double));
		}
	}

	//Read through the files.
	MPI_Barrier(MPI_COMM_WORLD);
	for(k=0;k<numfiles;k++)
	{
		if(my_rank == 0)
		{
			//For compressed file
			sprintf(filename,"%s%i.gz",strfilename,first+k);
			compfile = fopen(filename,"r");
			if(!compfile)
			{
				printf("[%i]: Error opening file %s.\n",my_rank,filename);
				MPI_Abort(MPI_COMM_WORLD,1);
			}

			sprintf(filename,"%s%i_%i",strfilename,first+k,my_rank);
			stormdata = fopen(filename,"w");
			if(!stormdata)	printf("[%i]: Error opening file %s.\n",my_rank,filename);

			int ret = uncompress_gzfile(compfile,stormdata);
			if(ret != Z_OK)		zerr(ret);
			fclose(compfile);
			fclose(stormdata);

			sprintf(filename,"%s%i_%i",strfilename,first+k,my_rank);
			stormdata = fopen(filename,"r");
			if(!stormdata)	printf("[%i]: Error opening file %s.\n",my_rank,filename);

			for(i=0;i<N;i++)
			{
				//Find the location of ID i (should be i by assumption)
				curr_idx = id_to_loc[i][1];
				fread(&rainfall_buffer,sizeof(float),1,stormdata);

				if(assignments[curr_idx] == my_rank)	//Data needed on proc 0; don't send
				{
					//This assumes the files have a different endianness from the system
					holder = *(unsigned int*) &rainfall_buffer;	//Pointers are fun!
					holder = (((holder & 0x0000ffff)<<16) | ((holder & 0xffff0000)>>16));
					holder = (((holder & 0x00ff00ff)<<8) | ((holder & 0xff00ff00)>>8));
					rainfall_buffer = *(float*) &holder;

					//Store the data
					sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k][0] = t_0 + k*increment;
					sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k][1] = rainfall_buffer;					
				}
				else	//Send it to the correct proc
				{
					int pos = 0;
					MPI_Pack(&curr_idx,1,MPI_INT,transferbuffer,buffer_size,&pos,MPI_COMM_WORLD);
					MPI_Pack(&rainfall_buffer,1,MPI_FLOAT,transferbuffer,buffer_size,&pos,MPI_COMM_WORLD);
					MPI_Send(transferbuffer,buffer_size,MPI_PACKED,assignments[curr_idx],N,MPI_COMM_WORLD);		//N = tag, as typical communication should not send this many links.
				}
			}

			fclose(stormdata);
			remove(filename);
		}
		else
		{
			//MPI_Status status;
			for(i=0;i<my_N;i++)
			{
				int pos = 0;
				MPI_Recv(transferbuffer,buffer_size,MPI_PACKED,0,N,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				MPI_Unpack(transferbuffer,buffer_size,&pos,&curr_idx,1,MPI_INT,MPI_COMM_WORLD);
				MPI_Unpack(transferbuffer,buffer_size,&pos,&rainfall_buffer,1,MPI_FLOAT,MPI_COMM_WORLD);

				//This assumes the files have a different endianness from the system
				holder = *(unsigned int*) &rainfall_buffer;	//Pointers are fun!
				holder = (((holder & 0x0000ffff)<<16) | ((holder & 0xffff0000)>>16));
				holder = (((holder & 0x00ff00ff)<<8) | ((holder & 0xff00ff00)>>8));
				rainfall_buffer = *(float*) &holder;

				//Store the data
				sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k][0] = t_0 + k*increment;
				sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k][1] = rainfall_buffer;
			}
		}
	}

	if(my_rank == 0)
		printf("Read %i binary files.\n",numfiles);

	//Add in terms for no rainfall if max_files > numfiles
	for(i=0;i<my_N;i++)
	{
		curr_idx = sys[my_sys[i]]->location;
		for(j=numfiles;j<max_files;j++)
		{
			sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[j][0] = sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[j-1][0] + .0001;
			sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[j][1] = 0.0;
		}
	}

	//Add a ceiling term
	for(i=0;i<my_N;i++)
	{
		curr_idx = sys[my_sys[i]]->location;
		//sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[max_files][0] = GlobalVars->maxtime + 1.0;
		sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[max_files][0] = ceil_time;
		sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[max_files][1] = -1.0;
	}

	//Calculate the first rain change time and set rain_value
	for(i=0;i<my_N;i++)
	{
		current = sys[my_sys[i]];
		rainfall_buffer = current->forcing_buff[forcing_idx]->rainfall[0][1];
		current->forcing_values[forcing_idx] = rainfall_buffer;
		current->forcing_indices[forcing_idx] = 0;

		for(j=1;j<current->forcing_buff[forcing_idx]->n_times;j++)
		{
			if(rainfall_buffer != current->forcing_buff[forcing_idx]->rainfall[j][1])
			{
				current->forcing_change_times[forcing_idx] = current->forcing_buff[forcing_idx]->rainfall[j][0];
				break;
			}
		}
		if(j == current->forcing_buff[forcing_idx]->n_times)
			current->forcing_change_times[forcing_idx] = current->forcing_buff[forcing_idx]->rainfall[j-1][0];
	}

	return 0;
}


//This reads in a set of binary files for the rainfall at each link.
//The data is given as intensities per grid cell.
//Link** sys: An array of links.
//int N: The number of links in sys.
//int my_N: The number of links assigned to this process.
//UnivVars* GlobalVars: Contains all the information that is shared by every link in the system.
//int* my_sys: Array of links assigned to this process (value is location in sys array).
//int* assignments (set by this method): Will be an array with N entries. assignments[i] will the process link sys[i] is assigned to.
//char strfilename[]: String of the filename for the rain file to read (NOT .str files). Filenames should be indexed.
//unsigned int first: The index of the file to read first.
//unsigned int last: The index of the file to read last.
//double t_0: The time at which the first file starts.
//double increment: The amount of time between consecutively indexed files.
//int** id_to_loc (set by this method): Will be an array with N rows and 2 columns, sorted by first col. First col is a link id and second is
//				the location of the id in sys.
//unsigned int max_files: The maximum number of files to be read.
int Create_Rain_Data_Grid(Link** sys,unsigned int N,unsigned int my_N,UnivVars* GlobalVars,unsigned int* my_sys,int* assignments,char strfilename[],unsigned int first,unsigned int last,double t_0,double increment,Forcing* forcing,unsigned int** id_to_loc,unsigned int max_files,unsigned int forcing_idx)
{
	unsigned int i,j,curr_idx,k,holder,endianness,cell;
	short unsigned int intensity;
	Link* current;
	char filename[128];
	float forcing_buffer;
	FILE* stormdata = NULL;
	unsigned int numfiles = last - first + 1;
	size_t result;

	//This is a time larger than any time in which the integrator is expected to get
	double ceil_time = 1e300;
	if(sys[my_sys[0]]->last_t > ceil_time*0.1)
		printf("[%i]: Warning: integrator time is extremely large (about %e). Loss of precision may occur.\n",my_rank,sys[my_sys[i]]->last_t);

	//Check that space for rain data has been allocated.
	if(sys[my_sys[0]]->forcing_buff[forcing_idx] == NULL)
	{
		for(i=0;i<my_N;i++)
		{
			sys[my_sys[i]]->forcing_buff[forcing_idx] = (ForcingData*) malloc(sizeof(ForcingData));
			sys[my_sys[i]]->forcing_buff[forcing_idx]->rainfall = (double**) malloc((max_files + 1)*sizeof(double*));
			sys[my_sys[i]]->forcing_buff[forcing_idx]->n_times = numfiles + 1;
			for(j=0;j<max_files+1;j++)	sys[my_sys[i]]->forcing_buff[forcing_idx]->rainfall[j] = (double*) malloc(2*sizeof(double));
		}
	}

	//Read through the files.
	for(k=0;k<numfiles;k++)
	{
		if(my_rank == 0)
		{
			sprintf(filename,"%s%i",strfilename,first+k);
			stormdata = fopen(filename,"r");
			if(stormdata)
			{
				for(i=0;i<forcing->num_cells;i++)
					forcing->received[i] = 0;

				//Check endianness
				fread(&i,sizeof(unsigned int),1,stormdata);
				if(i == 0x1)			endianness = 0;
				else if(i == 0x80000000)	endianness = 1;
				else
				{
					printf("Error: Cannot read endianness flag in binary file %s.\n",filename);
					MPI_Abort(MPI_COMM_WORLD,1);
				}

				//Read file
				if(endianness)
				{
					while(!feof(stormdata))
					{
						//Read intensity
						result = fread(&cell,sizeof(unsigned int),1,stormdata);
						if(!result)	break;
						fread(&intensity,sizeof(short unsigned int),1,stormdata);

						//Swap byte order
						holder = (((cell & 0x0000ffff)<<16) | ((cell & 0xffff0000)>>16));
						cell = (((holder & 0x00ff00ff)<<8) | ((holder & 0xff00ff00)>>8));
						intensity = (((intensity & 0x00ff00ff)<<8) | ((intensity & 0xff00ff00)>>8));

						if(cell < forcing->num_cells)
						{
							if(forcing->received[cell])	printf("Warning: Received multiple intensities for cell %u in file %s.\n",cell,filename);
							forcing->received[cell] = 1;
							forcing->intensities[cell] = (short int) intensity * forcing->factor;
						}
						else
							printf("Warning: bad grid cell id in file %s.\n",filename);
					}
				}
				else
				{
					while(!feof(stormdata))
					{
						//Read intensity
						result = fread(&cell,sizeof(unsigned int),1,stormdata);
						if(!result)	break;
						fread(&intensity,sizeof(short unsigned int),1,stormdata);

						if(cell < forcing->num_cells)
						{
							if(forcing->received[cell])	printf("Warning: Received multiple intensities for cell %u in file %s.\n",cell,filename);
							forcing->received[cell] = 1;
							forcing->intensities[cell] = (short int) intensity * forcing->factor;
						}
						else
							printf("Warning: bad grid cell id in file %s.\n",filename);
					}
				}

				fclose(stormdata);

				//Store 0's for remaining cells
				for(i=0;i<forcing->num_cells;i++)
					if(!forcing->received[i])	forcing->intensities[i] = 0.0;
			}
			else	//No file, no rain
			{
				for(i=0;i<forcing->num_cells;i++)
					forcing->intensities[i] = 0.0;
			}
		}

		MPI_Bcast(forcing->intensities,forcing->num_cells,MPI_FLOAT,0,MPI_COMM_WORLD);

		//Load the data
		for(cell=0;cell<forcing->num_cells;cell++)
		{
			for(i=0;i<forcing->num_links_in_grid[cell];i++)	//!!!! Assuming only links on this proc !!!!
			{
				curr_idx = forcing->grid_to_linkid[cell][i];
				sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k][0] = t_0 + k*increment;
				sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k][1] = forcing->intensities[cell];
			}
		}
	}

	//Add in terms for no rainfall if max_files > numfiles
	for(i=0;i<my_N;i++)
	{
		curr_idx = sys[my_sys[i]]->location;
		for(j=numfiles;j<max_files;j++)
		{
			sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[j][0] = sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[j-1][0] + .0001;
			sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[j][1] = 0.0;
		}
	}

	//Add a ceiling term
	for(i=0;i<my_N;i++)
	{
		curr_idx = sys[my_sys[i]]->location;
		//sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[max_files][0] = GlobalVars->maxtime + 1.0;
		sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[max_files][0] = ceil_time;
		sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[max_files][1] = -1.0;
	}

	//Calculate the first rain change time and set rain_value
	for(i=0;i<my_N;i++)
	{
		current = sys[my_sys[i]];
		forcing_buffer = current->forcing_buff[forcing_idx]->rainfall[0][1];
		current->forcing_values[forcing_idx] = forcing_buffer;
		current->forcing_indices[forcing_idx] = 0;

		for(j=1;j<current->forcing_buff[forcing_idx]->n_times;j++)
		{
			if( fabs(forcing_buffer - current->forcing_buff[forcing_idx]->rainfall[j][1]) > 1e-14 )
			{
				current->forcing_change_times[forcing_idx] = current->forcing_buff[forcing_idx]->rainfall[j][0];
				break;
			}
		}
		if(j == current->forcing_buff[forcing_idx]->n_times)
		{
			current->forcing_change_times[forcing_idx] = current->forcing_buff[forcing_idx]->rainfall[j-1][0];
		}
	}

	return 0;
}

//This reads in rainfall data at each link from an SQL database.
//Link** sys: An array of links.
//int N: The number of links in sys.
//int my_N: The number of links assigned to this process.
//UnivVars* GlobalVars: Contains all the information that is shared by every link in the system.
//int* my_sys: Array of links assigned to this process (value is location in sys array).
//int* assignments (set by this method): Will be an array with N entries. assignments[i] will the process link sys[i] is assigned to.
//char strfilename[]: String of the filename for the rain file to read (NOT .str files). Filenames should be indexed.
//unsigned int first: The index of the file to read first.
//unsigned int last: The index of the file to read last.
//double t_0: The time at which the first file starts.
//double increment: The amount of time between consecutively indexed files.
//int** id_to_loc (set by this method): Will be an array with N rows and 2 columns, sorted by first col. First col is a link id and second is
//				the location of the id in sys.
//unsigned int max_files: The maximum number of files to be read.
int Create_Rain_Database(Link** sys,unsigned int N,unsigned int my_N,UnivVars* GlobalVars,unsigned int* my_sys,int* assignments,ConnData *conninfo,unsigned int first,unsigned int last,Forcing* forcing,unsigned int** id_to_loc,double maxtime,unsigned int forcing_idx)
{
	unsigned int i,j,k,curr_idx,tuple_count;
	Link* current;
	float forcing_buffer;
	int received_time;
	char* query = conninfo->query;
	PGresult *res;
	unsigned int *db_unix_time,*db_link_id;
	float *db_rain_intens;

	unsigned int *total_times = (unsigned int*) calloc(my_N,sizeof(unsigned int));
	for(i=0;i<my_N;i++)
	{
		if(sys[my_sys[i]]->forcing_buff[forcing_idx])
		{
			total_times[i] = sys[my_sys[i]]->forcing_buff[forcing_idx]->n_times;
			sys[my_sys[i]]->forcing_buff[forcing_idx]->n_times = 1;
		}
	}

	//This is a time larger than any time in which the integrator is expected to get
	double ceil_time = 1e300;
	static short int gave_warning = 0;
	if(sys[my_sys[0]]->last_t > ceil_time*0.1 && !gave_warning)
	{
		gave_warning = 1;
		printf("[%i]: Warning: integrator time is extremely large (about %e). Loss of precision may occur.\n",my_rank,sys[my_sys[i]]->last_t);
	}

/*
	//!!!! Fix this (or don't?) !!!!
	total_times = sys[my_sys[0]]->forcing_buff[forcing_idx]->n_times;
	for(i=0;i<my_N;i++)	sys[my_sys[i]]->forcing_buff[forcing_idx]->n_times = 1;
*/

//GlobalVars->outletlink = 318213;

	//Query the database
	if(my_rank == 0)
	{
		//Connect to the database
		ConnectPGDB(conninfo);

		if(GlobalVars->outletlink == 0)
			sprintf(query,conninfo->queries[0],first,last);
		else
			sprintf(query,conninfo->queries[1],GlobalVars->outletlink,first,last);
//printf("*************************\n");
//printf("First = %u Last = %u t = %e increment = %u\n",first,last,sys[my_sys[0]]->last_t,forcing->increment);
//printf("*************************\n");
//printf("Gmax = %e maxtime = %e\n",GlobalVars->maxtime,maxtime);
//printf("query: %s\n",query);
//printf("*************************\n");
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"downloading rainfall data");
		tuple_count = PQntuples(res);
printf("Received %u intensities.\n",tuple_count);
		//Disconnect
		DisconnectPGDB(conninfo);

		//Allocate space
		MPI_Bcast(&tuple_count,1,MPI_INT,0,MPI_COMM_WORLD);
		db_unix_time = malloc(tuple_count*sizeof(unsigned int));
		db_rain_intens = malloc(tuple_count*sizeof(float));
		db_link_id = malloc(tuple_count*sizeof(unsigned int));

		//Load up the buffers
		for(i=0;i<tuple_count;i++)
		{
			db_unix_time[i] = atoi(PQgetvalue(res,i,0));
			db_rain_intens[i] = atof(PQgetvalue(res,i,1));
			db_link_id[i] = atoi(PQgetvalue(res,i,2));
		}

		//Broadcast the data
		MPI_Bcast(db_unix_time,tuple_count,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(db_rain_intens,tuple_count,MPI_FLOAT,0,MPI_COMM_WORLD);
		MPI_Bcast(db_link_id,tuple_count,MPI_INT,0,MPI_COMM_WORLD);

		//Clean up
		PQclear(res);
	}
	else
	{
		//Allocate space
		MPI_Bcast(&tuple_count,1,MPI_INT,0,MPI_COMM_WORLD);
		db_unix_time = malloc(tuple_count*sizeof(unsigned int));
		db_rain_intens = malloc(tuple_count*sizeof(float));
		db_link_id = malloc(tuple_count*sizeof(unsigned int));

		//Receive the data
		MPI_Bcast(db_unix_time,tuple_count,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(db_rain_intens,tuple_count,MPI_FLOAT,0,MPI_COMM_WORLD);
		MPI_Bcast(db_link_id,tuple_count,MPI_INT,0,MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);

//double stop = time(NULL);
//if(my_rank == 0)	printf("Total time to get rain: %f\n%u %u\n",difftime(stop,start),first,last);

	//Setup initial time in rainfall data with 0 rain
	for(i=0;i<my_N;i++)
	{
		sys[my_sys[i]]->forcing_buff[forcing_idx]->rainfall[0][0] = (int)(first - forcing->raindb_start_time)/60.0;
		sys[my_sys[i]]->forcing_buff[forcing_idx]->rainfall[0][1] = 0.0;
	}
//printf("Started with first = %u, raindb_start = %u, %f\n",first,forcing->raindb_start_time,(int)(first - forcing->raindb_start_time)/60.0);
	//Setup the data received
	for(i=0;i<tuple_count;i++)
	{
		//Find the location of ID i (ids should be numbered from 2 to whatever)
		//!!!! This will need to be a search at some point !!!!
/*
		if(GlobalVars->outletlink == 0)
		{
			curr_idx = id_to_loc[db_link_id[i]-2][1];
			//curr_idx = id_to_loc[atoi(PQgetvalue(res,i,2))-2][1];
			if(sys[curr_idx]->ID != db_link_id[i])
				printf("Indices do not match %u %u\n",sys[curr_idx]->ID,db_link_id[i]);
			//if(sys[curr_idx]->ID != atoi(PQgetvalue(res,i,2)))
				//printf("Indices do not match %u %u\n",sys[curr_idx]->ID,atoi(PQgetvalue(res,i,2)));
		}
		else
*/
			curr_idx = find_link_by_idtoloc(db_link_id[i],id_to_loc,N);
			//curr_idx = find_link_by_idtoloc(atoi(PQgetvalue(res,i,2)),id_to_loc,N);

		if(curr_idx < N && assignments[curr_idx] == my_rank)
		{
			k = sys[curr_idx]->forcing_buff[forcing_idx]->n_times;
			received_time = db_unix_time[i] - forcing->raindb_start_time;	//In seconds
//printf("Got k = %i received_time = %i in secs = %i in mins = %f\n",k,received_time,(int)(sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k-1][0] * 60.0 + 0.01),sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k-1][0]);
			if(received_time > (int) (sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k-1][0] * 60.0 + 0.01))
			{
				if( received_time <= (int) (sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k-1][0]*60.0) + (unsigned int) (forcing->file_time*60.0) || sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k-1][1] == 0.0)
				{
//printf("stored ID = %u i = %i k = %i received = %i unix_time = %i raindb_start = %i\n",sys[curr_idx]->ID,i,k,received_time,db_unix_time[i],forcing->raindb_start_time);
					sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k][0] = received_time / 60.0;
					sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k][1] = db_rain_intens[i];
					(sys[curr_idx]->forcing_buff[forcing_idx]->n_times)++;
//printf("k = %i (%f %f)  Also: received = %i next = %u\n",k,sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k][0],sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k][1],received_time,(int) (sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k-1][0]*60.0) + (unsigned int) (forcing->file_time*60.0));
				}
				else	//Add a 0 rainfall data
				{
//printf("stored 0 i = %i k = %i received = %i unix_time = %i raindb_start = %i\n",i,k,received_time,db_unix_time[i],forcing->raindb_start_time);
					sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k][0] = sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k-1][0] + forcing->file_time;
					sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k][1] = 0.0;
					sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k+1][0] = received_time / 60.0;
					sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k+1][1] = db_rain_intens[i];
					sys[curr_idx]->forcing_buff[forcing_idx]->n_times += 2;
//printf("(%f %f)\n",sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k][0],sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k][1]);
//printf("(%f %f)\n",sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k+1][0],sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k+1][1]);
				}
			}
			else if(received_time <= (int) (sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k-1][0] * 60.0 + 0.01))	//If the initial rate needs to be reset
{
//printf("stored init i = %i k = %i received = %i unix_time = %i raindb_start = %i\n",i,k,received_time,db_unix_time[i],forcing->raindb_start_time);
				sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k-1][0] = received_time / 60.0;
				sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k-1][1] = db_rain_intens[i];
//printf("k = %i Init (%f %f)\n",k,sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k-1][0],sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k-1][1]);
}
else
{
printf("!!!! Uh oh... !!!!\n");
printf("!!!! i = %i k = %i received = %i unix_time = %i raindb_start = %i\n",i,k,received_time,db_unix_time[i],forcing->raindb_start_time);
}
		}
	}
//printf("Got %u\n",sys[0]->forcing_buff[forcing_idx]->n_times);
	//Add ceiling terms
	for(i=0;i<my_N;i++)
	{
		curr_idx = sys[my_sys[i]]->location;
		k = sys[curr_idx]->forcing_buff[forcing_idx]->n_times;
		if(sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k-1][1] == 0.0)	//No rain, add just a ceiling
		{
			//sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k][0] = maxtime * (1.1) + 1.0;
			sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k][0] = ceil_time;
			sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k][1] = -1.0;
		}
		else	//Add a 0.0, and a ceiling
		{
			sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k][0] = sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k-1][0] + forcing->file_time;
			sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k][1] = 0.0;
			//sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k+1][0] = maxtime * (1.1) + 1.0;
			sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k+1][0] = ceil_time;
			sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k+1][1] = -1.0;
		}
	}

	//Reset n_times  !!!! Fix (well, this might be ok to do) !!!!
	//for(i=0;i<my_N;i++)	sys[my_sys[i]]->forcing_buff[forcing_idx]->n_times = total_times;
	for(i=0;i<my_N;i++)	sys[my_sys[i]]->forcing_buff[forcing_idx]->n_times = total_times[i];

	//Calculate the first rain change time and set rain_value
	for(i=0;i<my_N;i++)
	{
		current = sys[my_sys[i]];

		//Get to the time block that corresponds to the current time, set forcing value
		for(j=1;j<current->forcing_buff[forcing_idx]->n_times;j++)
		{
			if(current->last_t < current->forcing_buff[forcing_idx]->rainfall[j][0] - 1e-12)
				break;
		}

		forcing_buffer = current->forcing_buff[forcing_idx]->rainfall[j-1][1];
		current->forcing_values[forcing_idx] = forcing_buffer;
		current->forcing_indices[forcing_idx] = j-1;

		for(;j<current->forcing_buff[forcing_idx]->n_times;j++)
		{
			if( fabs(forcing_buffer - current->forcing_buff[forcing_idx]->rainfall[j][1]) > 1e-12 )
			{
				current->forcing_change_times[forcing_idx] = current->forcing_buff[forcing_idx]->rainfall[j][0];
				break;
			}
		}
		if(j == current->forcing_buff[forcing_idx]->n_times)
			current->forcing_change_times[forcing_idx] = current->forcing_buff[forcing_idx]->rainfall[j-1][0];
	}

	//Clean up
	free(total_times);
	free(db_unix_time);
	free(db_link_id);
	free(db_rain_intens);

	return 0;
}


//This reads in rainfall data at each link from an SQL database. The timestamps are assumed to be irregularly spaced.
//Link** sys: An array of links.
//int N: The number of links in sys.
//int my_N: The number of links assigned to this process.
//UnivVars* GlobalVars: Contains all the information that is shared by every link in the system.
//int* my_sys: Array of links assigned to this process (value is location in sys array).
//int* assignments (set by this method): Will be an array with N entries. assignments[i] will the process link sys[i] is assigned to.
//char strfilename[]: String of the filename for the rain file to read (NOT .str files). Filenames should be indexed.
//unsigned int first: The index of the file to read first.
//unsigned int last: The index of the file to read last.
//double t_0: The time at which the first file starts.
//double increment: The amount of time between consecutively indexed files.
//int** id_to_loc (set by this method): Will be an array with N rows and 2 columns, sorted by first col. First col is a link id and second is
//				the location of the id in sys.
//unsigned int max_files: The maximum number of files to be read.
int Create_Rain_Database_Irregular(Link** sys,unsigned int N,unsigned int my_N,UnivVars* GlobalVars,unsigned int* my_sys,int* assignments,ConnData *conninfo,unsigned int first,unsigned int last,Forcing* forcing,unsigned int** id_to_loc,double maxtime,unsigned int forcing_idx)
{
	unsigned int i,j,k,curr_idx,tuple_count,current_timestamp;
	Link* current;
	float forcing_buffer;
	int received_time;
	char* query = conninfo->query;
	PGresult *res;
	unsigned int *db_unix_time,*db_link_id;
	float *db_rain_intens;
	unsigned int *actual_timestamps,num_actual_timestamps,max_timestamps = forcing->increment;	//The maximum number of intensities to get for each link

/*
	unsigned int *total_times = (unsigned int*) calloc(my_N,sizeof(unsigned int));
	for(i=0;i<my_N;i++)
	{
		if(sys[my_sys[i]]->forcing_buff[forcing_idx])
		{
			total_times[i] = sys[my_sys[i]]->forcing_buff[forcing_idx]->n_times;
			sys[my_sys[i]]->forcing_buff[forcing_idx]->n_times = 1;
		}
	}
*/

	//This is a time larger than any time in which the integrator is expected to get
	double ceil_time = 1e300;
	static short int gave_warning = 0;
	if(sys[my_sys[0]]->last_t > ceil_time*0.1 && !gave_warning)
	{
		gave_warning = 1;
		printf("[%i]: Warning: integrator time is extremely large (about %e). Loss of precision may occur.\n",my_rank,sys[my_sys[i]]->last_t);
	}

	//Query the database
	if(my_rank == 0)
	{
		//Connect to the database
		ConnectPGDB(conninfo);

		//Download timestamps
		sprintf(query,conninfo->queries[3],first,max_timestamps);
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"downloading rainfall timestamps");
		num_actual_timestamps = PQntuples(res);

		//Unpack and broadcast the timestamps
		actual_timestamps = (unsigned int*) malloc(num_actual_timestamps*sizeof(unsigned int));
		for(i=0;i<num_actual_timestamps;i++)
			actual_timestamps[i] = atoi(PQgetvalue(res,i,0));
		MPI_Bcast(&num_actual_timestamps,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
		MPI_Bcast(actual_timestamps,num_actual_timestamps,MPI_UNSIGNED,0,MPI_COMM_WORLD);
		last = actual_timestamps[num_actual_timestamps-1];
		PQclear(res);

		//Download intensities
		if(GlobalVars->outletlink == 0)
			sprintf(query,conninfo->queries[0],first,last);
		else
			sprintf(query,conninfo->queries[1],GlobalVars->outletlink,first,last);
//printf("*************************\n");
//printf("First = %u Last = %u t = %e increment = %u\n",first,last,sys[my_sys[0]]->last_t,forcing->increment);
//printf("*************************\n");
//printf("Gmax = %e maxtime = %e\n",GlobalVars->maxtime,maxtime);
//printf("query: %s\n",query);
//printf("*************************\n");
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"downloading rainfall data");
		tuple_count = PQntuples(res);
printf("Received %u intensities.\n",tuple_count);
		//Disconnect
		DisconnectPGDB(conninfo);

		//Allocate space
		MPI_Bcast(&tuple_count,1,MPI_INT,0,MPI_COMM_WORLD);
		db_unix_time = malloc(tuple_count*sizeof(unsigned int));
		db_rain_intens = malloc(tuple_count*sizeof(float));
		db_link_id = malloc(tuple_count*sizeof(unsigned int));

		//Load up the buffers
		for(i=0;i<tuple_count;i++)
		{
			db_unix_time[i] = atoi(PQgetvalue(res,i,0));
			db_rain_intens[i] = atof(PQgetvalue(res,i,1));
			db_link_id[i] = atoi(PQgetvalue(res,i,2));
		}

		//Broadcast the data
		MPI_Bcast(db_unix_time,tuple_count,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(db_rain_intens,tuple_count,MPI_FLOAT,0,MPI_COMM_WORLD);
		MPI_Bcast(db_link_id,tuple_count,MPI_INT,0,MPI_COMM_WORLD);

		//Clean up
		PQclear(res);
	}
	else
	{
		MPI_Bcast(&num_actual_timestamps,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
		actual_timestamps = (unsigned int*) malloc(num_actual_timestamps*sizeof(unsigned int));
		MPI_Bcast(actual_timestamps,num_actual_timestamps,MPI_UNSIGNED,0,MPI_COMM_WORLD);
		last = actual_timestamps[num_actual_timestamps-1];

		//Allocate space
		MPI_Bcast(&tuple_count,1,MPI_INT,0,MPI_COMM_WORLD);
		db_unix_time = malloc(tuple_count*sizeof(unsigned int));
		db_rain_intens = malloc(tuple_count*sizeof(float));
		db_link_id = malloc(tuple_count*sizeof(unsigned int));

		//Receive the data
		MPI_Bcast(db_unix_time,tuple_count,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(db_rain_intens,tuple_count,MPI_FLOAT,0,MPI_COMM_WORLD);
		MPI_Bcast(db_link_id,tuple_count,MPI_INT,0,MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);

//double stop = time(NULL);
//if(my_rank == 0)	printf("Total time to get rain: %f\n%u %u\n",difftime(stop,start),first,last);

/*
printf("+++++++++\n");
for(i=0;i<num_actual_timestamps;i++)
{
printf("%u\n",actual_timestamps[i]);
}
printf("+++++++++\n");
*/

	//Set the times and zero out the intensities (some of these might change later)
	for(i=0;i<my_N;i++)
	{
		for(j=0;j<num_actual_timestamps;j++)
		{
			sys[my_sys[i]]->forcing_buff[forcing_idx]->rainfall[j][0] = (double)(actual_timestamps[j] - forcing->raindb_start_time)/60.0;
			sys[my_sys[i]]->forcing_buff[forcing_idx]->rainfall[j][1] = 0.0;
		}
	}

	//Setup the forcing values
	i = 0;
	for(j=0;j<num_actual_timestamps;j++)
	{
		current_timestamp = actual_timestamps[j];
		for(;i < tuple_count && db_unix_time[i] <= current_timestamp;i++)
		{
			received_time = db_unix_time[i] - forcing->raindb_start_time;	//In seconds
			curr_idx = find_link_by_idtoloc(db_link_id[i],id_to_loc,N);
			if(curr_idx < N && assignments[curr_idx] == my_rank)
			{
				sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[j][0] = received_time / 60.0;
				sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[j][1] = db_rain_intens[i];
//if(sys[curr_idx]->ID == 456117)
//printf("j = %u received_time = %u intensity = %f\n",j,received_time,sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[j][1]);
			}
//else
//printf("!!!! Huh, wtf? !!!!\n");
		}
	}

/*
	current_timestamp = actual_timestamps[0];
	for(i=0;i < tuple_count && current_timestamp <= db_unix_time[i];i++)
	{
		received_time = db_unix_time[i] - forcing->raindb_start_time;	//In seconds
		//if(received_time > (int) (sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k-1][0] * 60.0 + 0.01))	break;
		curr_idx = find_link_by_idtoloc(db_link_id[i],id_to_loc,N);
		if(curr_idx < N && assignments[curr_idx] == my_rank)
		{
			//k = sys[curr_idx]->forcing_buff[forcing_idx]->n_times;
			//if(received_time <= (int) (sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k-1][0] * 60.0 + 0.01))	//If the initial rate needs to be reset
			{
				sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[0][0] = received_time / 60.0;
				sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[0][1] = db_rain_intens[i];
			}
		}
	}
	
	//!!!! Uh, any difference between this and j=0? !!!!
	for(j=1;j<num_actual_timestamps;j++)
	{
		current_timestamp = actual_timestamps[j];
		for(;i < tuple_count && current_timestamp <= db_unix_time[i];i++)
		{
			received_time = db_unix_time[i] - forcing->raindb_start_time;	//In seconds
			curr_idx = find_link_by_idtoloc(db_link_id[i],id_to_loc,N);

			if(curr_idx < N && assignments[curr_idx] == my_rank)
			{
				sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[j][0] = received_time / 60.0;
				sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[j][1] = db_rain_intens[i];
			}
		}
	}
*/

/*
	//Setup the data received
	for(i=0;i<tuple_count;i++)
	{
		//Find the location of ID i (ids should be numbered from 2 to whatever)
		curr_idx = find_link_by_idtoloc(db_link_id[i],id_to_loc,N);

		if(curr_idx < N && assignments[curr_idx] == my_rank)
		{
			k = sys[curr_idx]->forcing_buff[forcing_idx]->n_times;
			received_time = db_unix_time[i] - forcing->raindb_start_time;	//In seconds

			if(received_time > (int) (sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k-1][0] * 60.0 + 0.01))
			{
				if( received_time <= (int) (sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k-1][0]*60.0) + (unsigned int) (forcing->file_time*60.0) || sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k-1][1] == 0.0)
				{
					sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k][0] = received_time / 60.0;
					sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k][1] = db_rain_intens[i];
					(sys[curr_idx]->forcing_buff[forcing_idx]->n_times)++;
				}
				else	//Add a 0 rainfall data before adding this data
				{
					sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k][0] = sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k-1][0] + forcing->file_time;
					sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k][1] = 0.0;
					sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k+1][0] = received_time / 60.0;
					sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k+1][1] = db_rain_intens[i];
					sys[curr_idx]->forcing_buff[forcing_idx]->n_times += 2;
				}
			}
			else //if(received_time <= (int) (sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k-1][0] * 60.0 + 0.01))	//If the initial rate needs to be reset
			{
				sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k-1][0] = received_time / 60.0;
				sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k-1][1] = db_rain_intens[i];
			}
		}
	}
*/
	//Add ceiling terms
	for(i=0;i<my_N;i++)
	{
		curr_idx = sys[my_sys[i]]->location;
		//k = sys[curr_idx]->forcing_buff[forcing_idx]->n_times;
		k = num_actual_timestamps;
		//if(sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k-1][1] == 0.0)	//No rain, add just a ceiling
		{
			sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k][0] = ceil_time;
			sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k][1] = -1.0;
		}
/*
		else	//Add a 0.0, and a ceiling
		{
			sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k][0] = sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k-1][0] + forcing->file_time;
			sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k][1] = 0.0;
			//sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k+1][0] = maxtime * (1.1) + 1.0;
			sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k+1][0] = ceil_time;
			sys[curr_idx]->forcing_buff[forcing_idx]->rainfall[k+1][1] = -1.0;
		}
*/
	}

	//Reset n_times  !!!! Fix (well, this might be ok to do) !!!!
	//for(i=0;i<my_N;i++)	sys[my_sys[i]]->forcing_buff[forcing_idx]->n_times = total_times[i];

	//Calculate the first rain change time and set rain_value
	for(i=0;i<my_N;i++)
	{
		current = sys[my_sys[i]];

		//Get to the time block that corresponds to the current time, set forcing value
		for(j=1;j<current->forcing_buff[forcing_idx]->n_times;j++)
		{
			if(current->last_t < current->forcing_buff[forcing_idx]->rainfall[j][0] - 1e-12)
				break;
		}

		forcing_buffer = current->forcing_buff[forcing_idx]->rainfall[j-1][1];
		current->forcing_values[forcing_idx] = forcing_buffer;
		current->forcing_indices[forcing_idx] = j-1;

		for(;j<current->forcing_buff[forcing_idx]->n_times;j++)
		{
			if( fabs(forcing_buffer - current->forcing_buff[forcing_idx]->rainfall[j][1]) > 1e-12 )
			{
				current->forcing_change_times[forcing_idx] = current->forcing_buff[forcing_idx]->rainfall[j][0];
				break;
			}
		}
		if(j == current->forcing_buff[forcing_idx]->n_times)
			current->forcing_change_times[forcing_idx] = current->forcing_buff[forcing_idx]->rainfall[j-1][0];
	}

	//Clean up
	//free(total_times);
	free(actual_timestamps);
	free(db_unix_time);
	free(db_link_id);
	free(db_rain_intens);

	return last;
}


//Sets the rain data of link current to 0 for all times.
void SetRain0(Link** sys,unsigned int my_N,double maxtime,unsigned int* my_sys,UnivVars* GlobalVars,Forcing* forcing,unsigned int forcing_idx)
{
	unsigned int i,j,k;
	Link* current;

	//Check that space for rain data has been allocated.
	if(sys[my_sys[0]]->forcing_buff[forcing_idx] == NULL)
	{
		//k = 2;
		k = forcing->increment + 3;
		for(i=0;i<my_N;i++)
		{
			sys[my_sys[i]]->forcing_buff[forcing_idx] = (ForcingData*) malloc(sizeof(ForcingData));
			sys[my_sys[i]]->forcing_buff[forcing_idx]->rainfall = (double**) malloc(k*sizeof(double*));
			sys[my_sys[i]]->forcing_buff[forcing_idx]->n_times = k;
			for(j=0;j<k;j++)	sys[my_sys[i]]->forcing_buff[forcing_idx]->rainfall[j] = (double*) malloc(2*sizeof(double));
		}
	}	

	for(i=0;i<my_N;i++)
	{
		sys[my_sys[i]]->forcing_buff[forcing_idx]->rainfall[0][0] = 0.0;
		sys[my_sys[i]]->forcing_buff[forcing_idx]->rainfall[0][1] = 0.0;
		sys[my_sys[i]]->forcing_buff[forcing_idx]->rainfall[1][0] = maxtime * (1.1) + 1.0;
		sys[my_sys[i]]->forcing_buff[forcing_idx]->rainfall[1][1] = -1.0;
	}

	//Calculate the first rain change time and set rain_value
	for(i=0;i<my_N;i++)
	{
		current = sys[my_sys[i]];
		current->forcing_values[forcing_idx] = 0.0;
		current->forcing_indices[forcing_idx] = 0;
		current->forcing_change_times[forcing_idx] = current->forcing_buff[forcing_idx]->rainfall[1][0];
	}
}


//Sets a forcing based on monthly data.
//Assumes the rates are already set.
double CreateForcing_Monthly(Link** sys,unsigned int my_N,unsigned int* my_sys,UnivVars* GlobalVars,ForcingData* GlobalForcing,unsigned int forcing_idx,struct tm *current_time,time_t first_time,time_t last_time,double t_0)
{
	unsigned int j,k,num_months = 12;
	int i;
	int month_0,current_year;
	Link* current;
	char buffer[4];
	double t = t_0;

	//Find the current month
	strftime(buffer,4,"%m",current_time);
	month_0 = atoi(buffer) - 1;
	if(month_0 < 0 || month_0 > num_months)
	{
		printf("[%i]: Error: Bad month %i.\n",my_rank,month_0);
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	//Set the (local) times for the current month and previous months
	GlobalForcing->rainfall[month_0][0] = t_0;
	for(i=month_0-1;i>-1;i--)	GlobalForcing->rainfall[i][0] = GlobalForcing->rainfall[i+1][0] - 1.0;
	current_year = current_time->tm_year + 1900;

	//Find days until next month
	int days_this_month = days_in_month(month_0,current_year);
	int curr_day = current_time->tm_mday;
	t += (days_this_month - curr_day) * (24.0*60.0);	//Convert to mins

	//Add in hours and mins and seconds
	t += (double) ((23-current_time->tm_hour)*60 + (59-current_time->tm_min)) + (60.0-current_time->tm_sec)/60.0;

	//Save t
	GlobalForcing->rainfall[month_0+1][0] = t;

	//Set the (local) times for each future month
	for(i=month_0+2;i<=num_months;i++)	//This should set the ceiling term
	{
		t += days_in_month(i-1,current_year) * (24.0*60.0);
		GlobalForcing->rainfall[i][0] = t;
	}

	//Set ceiling term
	GlobalForcing->rainfall[num_months][1] = 0.0;

	//Check if this data goes past the last time
	double final_time = (last_time - first_time)/60.0;
	if(t > final_time)	//Need to set 0s past final_time
	{
		i = num_months-1;
		GlobalForcing->rainfall[i][0] = t+1.0;
		GlobalForcing->rainfall[i][1] = 0.0;

		for(i-=1;GlobalForcing->rainfall[i][0] > final_time;i--)
		{
			GlobalForcing->rainfall[i][0] = final_time;
			GlobalForcing->rainfall[i][1] = 0.0;
			GlobalForcing->rainfall[i+1][0] = t+1.0;
			GlobalForcing->rainfall[i+1][1] = 0.0;
		}
	}

	//Set the current forcing value at each link
	for(i=0;i<my_N;i++)
	{
		sys[my_sys[i]]->forcing_values[forcing_idx] = GlobalForcing->rainfall[month_0][1];
		sys[my_sys[i]]->forcing_indices[forcing_idx] = month_0;
		sys[my_sys[i]]->forcing_change_times[forcing_idx] = GlobalForcing->rainfall[month_0+1][0];
	}

	return t;
}

