#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <math.h>
#include <memory.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(HAVE_UNISTD_H)
#include <unistd.h>
#endif

#if defined(_MSC_VER)
#include <process.h>
#endif

#if defined(HAVE_POSTGRESQL)
#include <libpq-fe.h>
#endif

#if defined(HAVE_HDF5)
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

#include "riversys.h"

//Read topo data and build the network.
//Also creates id_to_loc.
Link** Create_River_Network(UnivVars* GlobalVars,unsigned int* N,unsigned int*** id_to_loc,ConnData** db_connections)
{
	Link **system;
	FILE* riverdata = NULL;
	PGresult *mainres,*res;
	unsigned int *link_ids,*dbres_link_id,*dbres_parent,sizeres = 0,i,j,k,**loc_to_children,*numparents,id,max_children = 10,*loc_to_children_array,curr_loc;
	link_ids = dbres_link_id = dbres_parent = NULL;
	int db;

	if(GlobalVars->rvr_flag == 0)	//Read topo data from file
	{
		if(my_rank == 0)
		{
			riverdata = fopen(GlobalVars->rvr_filename,"r");
			if(!riverdata)
			{
				if(my_rank == 0)	printf("Error: file %s not found for .rvr file.\n",GlobalVars->rvr_filename);
				*N = 0;
				return NULL;
			}
			if(CheckWinFormat(riverdata))
			{
				printf("Error: File %s appears to be in Windows format. Try converting to unix format using 'dos2unix' at the command line.\n",GlobalVars->rvr_filename);
				fclose(riverdata);
				return NULL;
			}

			fscanf(riverdata,"%u",N);
			link_ids = (unsigned int*) malloc(*N*sizeof(unsigned int));
			loc_to_children_array = (unsigned int*) calloc(*N*max_children,sizeof(unsigned int));
			loc_to_children = (unsigned int**) malloc(*N*sizeof(unsigned int*));	//This holds the ID of the children
			for(i=0;i<*N;i++)	loc_to_children[i] = &(loc_to_children_array[i*max_children]);
			numparents = (unsigned int*) malloc(*N*sizeof(unsigned int));

			for(i=0;i<*N;i++)
			{
				fscanf(riverdata,"%u %u",&(link_ids[i]),&(numparents[i]));
				for(j=0;j<numparents[i];j++)
				{
					fscanf(riverdata,"%u",&id);
					loc_to_children[i][j] = id;
				}

				if(numparents[i] > max_children)
				{
					printf("Error: assumed no link has more than %u parents, but link %u has %u.\n",max_children,link_ids[i],numparents[i]);
					printf("If this is not an error in the input data, modify max_children in riversys.c, function Create_River_Network.\n");
					return NULL;
				}
			}

			fclose(riverdata);

			//Broadcast data
			MPI_Bcast(N,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
			MPI_Bcast(link_ids,*N,MPI_UNSIGNED,0,MPI_COMM_WORLD);
			MPI_Bcast(loc_to_children_array,*N*max_children,MPI_UNSIGNED,0,MPI_COMM_WORLD);	//Yeah, this is probably not the most efficient...
			MPI_Bcast(numparents,*N,MPI_UNSIGNED,0,MPI_COMM_WORLD);
		}
		else
		{
			MPI_Bcast(N,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
			link_ids = (unsigned int*) malloc(*N*sizeof(unsigned int));
			loc_to_children_array = (unsigned int*) malloc(*N*max_children*sizeof(unsigned int));
			loc_to_children = (unsigned int**) malloc(*N*sizeof(unsigned int*));	//This holds the ID of the children
			for(i=0;i<*N;i++)	loc_to_children[i] = &(loc_to_children_array[i*max_children]);
			numparents = (unsigned int*) malloc(*N*sizeof(unsigned int));
			MPI_Bcast(link_ids,*N,MPI_UNSIGNED,0,MPI_COMM_WORLD);
			MPI_Bcast(loc_to_children_array,*N*max_children,MPI_UNSIGNED,0,MPI_COMM_WORLD);
			MPI_Bcast(numparents,*N,MPI_UNSIGNED,0,MPI_COMM_WORLD);
		}
	}
	else if(GlobalVars->rvr_flag == 1)	//Download topo data from database
	{
		//if(my_rank == 0)	printf("\nTransferring topology data from database...\n");
		//MPI_Barrier(MPI_COMM_WORLD);
		//start = time(NULL);

		if(my_rank == 0)
		{
			if(GlobalVars->outletlink == 0)	//Grab entire network
			{
				//Note: parent_id on the SQL database is like the child id here
				db = ConnectPGDB(db_connections[ASYNCH_DB_LOC_TOPO]);
				if(db)
				{
					printf("[%i]: Error connecting to the topology database.\n",my_rank);
					return NULL;
				}
				mainres = PQexec(db_connections[ASYNCH_DB_LOC_TOPO]->conn,db_connections[ASYNCH_DB_LOC_TOPO]->queries[0]);	//For all link ids
				res = PQexec(db_connections[ASYNCH_DB_LOC_TOPO]->conn,db_connections[ASYNCH_DB_LOC_TOPO]->queries[1]);		//For parent data
				if(CheckResError(res,"querying connectivity") || CheckResError(mainres,"querying DEM data"))
					return NULL;
				DisconnectPGDB(db_connections[ASYNCH_DB_LOC_TOPO]);

				*N = PQntuples(mainres);

				//Get the list of link ids
				link_ids = (unsigned int*) malloc(*N*sizeof(unsigned int));
				for(i=0;i<*N;i++)	link_ids[i] = atoi(PQgetvalue(mainres,i,0));
				PQclear(mainres);

				sizeres = PQntuples(res);
			}
			else	//Grab a sub basin (Note: all link ids (except the outlet) is in the second column of res)
			{
				db = ConnectPGDB(db_connections[ASYNCH_DB_LOC_TOPO]);
				if(db)
				{
					printf("[%i]: Error connecting to the topology database.\n",my_rank);
					return NULL;
				}

				//Make the queries
				//Be careful to not overload the database
				sprintf(db_connections[ASYNCH_DB_LOC_TOPO]->query,db_connections[ASYNCH_DB_LOC_TOPO]->queries[2],GlobalVars->outletlink);	//For parent data
				res = PQexec(db_connections[ASYNCH_DB_LOC_TOPO]->conn,db_connections[ASYNCH_DB_LOC_TOPO]->query);
				if(CheckResError(res,"querying connectivity"))	return NULL;
				DisconnectPGDB(db_connections[ASYNCH_DB_LOC_TOPO]);

				*N = PQntuples(res) + 1;

				//Get the list of link ids
				link_ids = (unsigned int*) malloc(*N*sizeof(unsigned int));
				for(i=0;i<*N-1;i++)	link_ids[i] = atoi(PQgetvalue(res,i,1));
				link_ids[i] = GlobalVars->outletlink;
				merge_sort_1D(link_ids,*N);

				sizeres = *N - 1;
			}

			dbres_link_id = (unsigned int*) malloc(sizeres*sizeof(unsigned int));
			dbres_parent = (unsigned int*) malloc(sizeres*sizeof(unsigned int));

			for(i=0;i<sizeres;i++)
			{
				dbres_link_id[i] = atoi(PQgetvalue(res,i,0));
				dbres_parent[i] = atoi(PQgetvalue(res,i,1));
			}

			PQclear(res);

			//Modify the data format
			loc_to_children_array = (unsigned int*) calloc(*N*max_children,sizeof(unsigned int));
			loc_to_children = (unsigned int**) malloc(*N*sizeof(unsigned int*));	//This holds the IDs of the parents. Really needs a better name...
			for(i=0;i<*N;i++)	loc_to_children[i] = &(loc_to_children_array[i*max_children]);
			numparents = (unsigned int*) malloc(*N*sizeof(unsigned int));

			curr_loc = 0;
			for(i=0;i<sizeres;i+=j)
			{
				//Select the next link
				id = dbres_link_id[i];

				if(id == link_ids[curr_loc])	//Extract parents
				{
					//Count the parents
					for(j=0;i+j<sizeres;j++)
						if(dbres_link_id[i+j] != id)	break;

					//Set the number of parents
					numparents[curr_loc] = j;

					//Set the information to find the child links
					for(k=0;k<j;k++)
						loc_to_children[curr_loc][k] = dbres_parent[i+k];
				}
				else	//No parents
				{
					numparents[curr_loc] = 0;
					j = 0;
				}

				//Set the next location
				curr_loc++;
			}
			for(;curr_loc<*N;curr_loc++)	//Finish any leaves at the end of link_ids
				numparents[curr_loc] = 0;

			//Send sizes and data to other processes
			MPI_Bcast(N,1,MPI_INT,0,MPI_COMM_WORLD);
			MPI_Bcast(link_ids,*N,MPI_UNSIGNED,0,MPI_COMM_WORLD);
			MPI_Bcast(loc_to_children_array,*N*max_children,MPI_UNSIGNED,0,MPI_COMM_WORLD);
			MPI_Bcast(numparents,*N,MPI_UNSIGNED,0,MPI_COMM_WORLD);

			//Clean up
			free(dbres_link_id);
			free(dbres_parent);
		}
		else
		{
			//Receive data
			MPI_Bcast(N,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
			link_ids = (unsigned int*) malloc(*N*sizeof(unsigned int));
			loc_to_children_array = (unsigned int*) calloc(*N*max_children,sizeof(unsigned int));
			loc_to_children = (unsigned int**) malloc(*N*sizeof(unsigned int*));	//This holds the ID of the children
			for(i=0;i<*N;i++)	loc_to_children[i] = &(loc_to_children_array[i*max_children]);
			numparents = (unsigned int*) calloc(*N,sizeof(unsigned int));
			MPI_Bcast(link_ids,*N,MPI_UNSIGNED,0,MPI_COMM_WORLD);
			MPI_Bcast(loc_to_children_array,*N*max_children,MPI_UNSIGNED,0,MPI_COMM_WORLD);
			MPI_Bcast(numparents,*N,MPI_UNSIGNED,0,MPI_COMM_WORLD);
		}


		//MPI_Barrier(MPI_COMM_WORLD);
		//stop = time(NULL);
		//if(my_rank == 0)	printf("Time to receive data: %f\n",difftime(stop,start));
	}
	else
	{
		if(my_rank == 0)	printf("Error: Bad topology flag %hi in .gbl file.\n",GlobalVars->rvr_flag);
		*N = 0;
		return NULL;
	}

	//Make a list of ids and locations, sorted by id
	*id_to_loc = malloc(*N*sizeof(int*));
	for(i=0;i<*N;i++)	(*id_to_loc)[i] = malloc(2*sizeof(int));	//!!!! This should really have an array... !!!!
	for(i=0;i<*N;i++)
	{
		(*id_to_loc)[i][0] = link_ids[i];
		(*id_to_loc)[i][1] = i;
	}
	merge_sort_ids(*id_to_loc,*N);

	//Check for crapiness
	if(my_rank == 0 && *N < (unsigned int)np)
		printf("\nWarning: using more processes (%i) than links (%u).\n",np,*N);

	//Allocate some space for the network
	system = (Link**) malloc(*N*sizeof(Link*));
	for(i=0;i<*N;i++)	system[i] = (Link*) malloc(sizeof(Link));

	//Build the network
	GlobalVars->max_parents = 0;
	for(i=0;i<*N;i++)
	{
		system[i]->location = i;
		system[i]->ID = link_ids[i];

		//Set the parents
		system[i]->numparents = numparents[i];
		system[i]->parents = (Link**) malloc(system[i]->numparents * sizeof(Link*));
		system[i]->c = NULL;
		GlobalVars->max_parents = (GlobalVars->max_parents > numparents[i]) ? GlobalVars->max_parents : numparents[i];

		//Set a few other data
		system[i]->output_user = NULL;
		system[i]->peakoutput_user = NULL;
		system[i]->f = NULL;
		system[i]->params = v_get(0);
		system[i]->dam = 0;
		system[i]->method = NULL;
		system[i]->errorinfo = NULL;
		system[i]->qvs = NULL;
		system[i]->discont = NULL;
		system[i]->discont_start = 0;
		system[i]->discont_end = GlobalVars->discont_size-1;
		system[i]->discont_count = 0;
		system[i]->discont_send = NULL;
		system[i]->discont_order_send = NULL;
		system[i]->discont_send_count = 0;
		system[i]->res = 0;
		system[i]->num_dense = 0;
		system[i]->dense_indices = NULL;
		system[i]->dim = 0;
		system[i]->diff_start = 0;
		system[i]->no_ini_start = 0;
		system[i]->disk_iterations = 0;
		system[i]->forcing_buff = NULL;
		system[i]->forcing_values = NULL;
		system[i]->forcing_change_times = NULL;
		system[i]->forcing_indices = NULL;
		system[i]->save_flag = 0;
		system[i]->peak_flag = 0;
		system[i]->user = NULL;
		//system[i]->equations = NULL;
	}

	//Setup the child and parent information
	for(i=0;i<*N;i++)
	{
		for(j=0;j<system[i]->numparents;j++)
		{
			curr_loc = find_link_by_idtoloc(loc_to_children[i][j],*id_to_loc,*N);
			if(curr_loc > *N)
			{
				if(my_rank == 0)	printf("Error: Invalid id in topology data (%u).\n",loc_to_children[i][j]);
				*N = 0;
				return NULL;
			}
			system[i]->parents[j] = system[curr_loc];
			system[curr_loc]->c = system[i];
		}
	}

	//Set an outletlink id. This only sets one outlet, and only if rvr_flag is not set.
	if(GlobalVars->prm_flag == 1 && GlobalVars->rvr_flag == 0)
	{
		for(i=0;i<*N;i++)
			if(system[i]->c == NULL)	break;
		GlobalVars->outletlink = system[i]->ID;
	}

	//Clean up
	free(loc_to_children);
	free(loc_to_children_array);
	free(link_ids);
	free(numparents);

	return system;
}



//Read in the local paramters for the network.
//Returns 1 if there is an error, 0 otherwise.
//If load_all == 1, then the parameters for every link are available on every proc.
//If load_all == 0, then the parameters are only available for links assigned to this proc.
int Load_Local_Parameters(Link** system,unsigned int N,unsigned int* my_sys,unsigned int my_N,int* assignments,short int* getting,unsigned int** id_to_loc,UnivVars* GlobalVars,ConnData** db_connections,short int load_all,model* custom_model,void* external)
{
	unsigned int i,j,*db_link_id,curr_loc;
	int db;
	double *db_params_array,**db_params;
	PGresult *res;
	FILE* paramdata;

	//Error checking
	if(!load_all && (!assignments || !getting))
	{
		if(my_rank == 0)
		{
			printf("Error loading link parameters: network partitioning must occur before loading parameters.\n");
			printf("(Hint: either partition the network before loading parameters, or load parameters at every link.)\n");
		}
		return 1;
	}

	//Allocate space
	db_link_id = (unsigned int*) malloc(N*sizeof(unsigned int));
	db_params_array = (double*) malloc(N*GlobalVars->disk_params*sizeof(double));
	db_params = (double**) malloc(N*sizeof(double*));
	for(i=0;i<N;i++)	db_params[i] = &(db_params_array[i*GlobalVars->disk_params]);

	//Read parameters
	if(my_rank == 0)
	{
		if(GlobalVars->prm_flag == 0)
		{
			paramdata = fopen(GlobalVars->prm_filename,"r");
			if(paramdata == NULL)
			{
				printf("Error: file %s not found for .prm file.\n",GlobalVars->prm_filename);
				return 1;
			}
			if(CheckWinFormat(paramdata))
			{
				printf("Error: File %s appears to be in Windows format. Try converting to unix format using 'dos2unix' at the command line.\n",GlobalVars->prm_filename);
				fclose(paramdata);
				return 1;
			}

			fscanf(paramdata,"%u",&i);
			if(i != N)
			{
				printf("Error: expected %u links in parameter file. Got %u.\n",N,i);
				return 1;
			}

			for(i=0;i<N;i++)
			{
				fscanf(paramdata,"%u",&(db_link_id[i]));
				for(j=0;j<GlobalVars->disk_params;j++)
				{
					if(fscanf(paramdata,"%lf",&(db_params[i][j])) == 0)
					{
						printf("Error reading from parameter file %s.\n",GlobalVars->prm_filename);
						return 1;
					}
				}
			}

			fclose(paramdata);
		}
		else if(GlobalVars->prm_flag == 1)
		{
			//printf("\nTransferring parameter data from database...\n");
			//start = time(NULL);

			if(GlobalVars->outletlink == 0)	//Grab entire network
			{
				db = ConnectPGDB(db_connections[ASYNCH_DB_LOC_PARAMS]);
				if(db)
				{
					printf("[%i]: Error connecting to the parameter database.\n",my_rank);
					return 1;
				}
				res = PQexec(db_connections[ASYNCH_DB_LOC_PARAMS]->conn,db_connections[ASYNCH_DB_LOC_PARAMS]->queries[0]);
				if(CheckResError(res,"querying DEM data"))	return 1;
				DisconnectPGDB(db_connections[ASYNCH_DB_LOC_PARAMS]);
			}
			else	//Grab a sub basin
			{
				db = ConnectPGDB(db_connections[ASYNCH_DB_LOC_PARAMS]);
				if(db)
				{
					printf("[%i]: Error connecting to the parameter database.\n",my_rank);
					return 1;
				}

				//Make the queries
				sprintf(db_connections[ASYNCH_DB_LOC_PARAMS]->query,db_connections[ASYNCH_DB_LOC_PARAMS]->queries[1],GlobalVars->outletlink);
				res = PQexec(db_connections[ASYNCH_DB_LOC_PARAMS]->conn,db_connections[ASYNCH_DB_LOC_PARAMS]->query);
				if(CheckResError(res,"querying DEM data"))	return 1;
				DisconnectPGDB(db_connections[ASYNCH_DB_LOC_PARAMS]);
			}

			i = PQntuples(res);
			if(i != N)
			{
				printf("Error processing link parameters: Got %u, expected %u.\n(Hint: make sure your topology and parameter sources have the same number of links.)\n",i,N);
				return 1;
			}

			//Load buffers
			for(i=0;i<N;i++)	db_link_id[i] = atoi(PQgetvalue(res,i,0));
			for(i=0;i<N;i++)
				for(j=0;j<GlobalVars->disk_params;j++)
					db_params[i][j] = atof(PQgetvalue(res,i,1+j));

			//Cleanup
			PQclear(res);
		}

		//stop = time(NULL);
		//printf("Time to receive data: %f\n",difftime(stop,start));
	}

	//Broadcast data
	MPI_Bcast(db_link_id,N,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Bcast(db_params_array,N*GlobalVars->disk_params,MPI_DOUBLE,0,MPI_COMM_WORLD);

	//Unpack the data
	if(load_all)
	{
		for(i=0;i<N;i++)
		{
			curr_loc = find_link_by_idtoloc(db_link_id[i],id_to_loc,N);
			if(curr_loc > N)
			{
				if(my_rank == 0)	printf("Error: link id %u appears in the link parameters, but not in the topology data.\n(Hint: make sure your topology and parameter sources are correct.)\n",db_link_id[i]);
				return 1;
			}
			system[curr_loc]->params = v_get(GlobalVars->params_size);
			for(j=0;j<GlobalVars->disk_params;j++)
				system[curr_loc]->params.ve[j] = db_params[i][j];

			if(custom_model)
                custom_model->Convert(system[curr_loc]->params,GlobalVars->type,external);
			else
                ConvertParams(system[curr_loc]->params,GlobalVars->type,external);
		}
	}
	else	//Just load the local parameters
	{
		for(i=0;i<N;i++)
		{
			curr_loc = find_link_by_idtoloc(db_link_id[i],id_to_loc,N);
			if(curr_loc > N)
			{
				if(my_rank == 0)	printf("Error: link id %u appears in the link parameters, but not in the topology data.\n(Hint: make sure your topology and parameter sources are correct.)\n",db_link_id[i]);
				return 1;
			}
			else
			{
				if(assignments[curr_loc] == my_rank || getting[curr_loc])
				{
					system[curr_loc]->params = v_get(GlobalVars->params_size);
					for(j=0;j<GlobalVars->disk_params;j++)
						system[curr_loc]->params.ve[j] = db_params[i][j];

					if(custom_model)	custom_model->Convert(system[curr_loc]->params,GlobalVars->type,external);
					else	ConvertParams(system[curr_loc]->params,GlobalVars->type,external);
				}
			}
		}
	}

	//Clean up
	free(db_link_id);
	free(db_params_array);
	free(db_params);

	return 0;
}



//Partitions the network amongst different MPI processes.
//!!!! Perhaps the leaves info could be moved deeper? How about errors here? !!!!
int Parition_Network(Link** system,unsigned int N,UnivVars* GlobalVars,unsigned int** my_sys,unsigned int* my_N,int** assignments,TransData** my_data,short int** getting,model* custom_model)
{
	unsigned int i,j;
	Link *current,*prev;//**upstream_order = (Link**) malloc(N*sizeof(Link*));

/*
	//Order the links by upstream area
	for(i=0;i<N;i++)
		upstream_order[i] = system[i];
	merge_sort(upstream_order,N,GlobalVars->area_idx);
*/

	//Perform a DFS to sort the leaves
	Link** stack = malloc(N*sizeof(Link*)); //Holds the index in system
	int stack_size = 0;

	Link** leaves = malloc(N*sizeof(Link*));
	unsigned int leaves_size = 0;
	unsigned short int numparents;

	for(j=0;j<N;j++)	//!!!! Iterate over a list of roots? No need to sort then... !!!!
	{
		if(system[j]->c == NULL)
		{
			stack[0] = system[j];
			stack_size = 1;
			while(stack_size > 0)
			{
				current = stack[stack_size-1];	//Top of stack
				numparents = current->numparents;

				if(numparents == 0)
				{
					stack_size--;
					leaves[leaves_size] = current;
					leaves_size++;
				}
				else
				{
					//If current is not a leaf, replace it with its parents
					for(i=0;i<numparents;i++)
					{
						stack[stack_size - 1 + i] = current->parents[numparents - 1 - i];
						//stack[stack_size - 1 + i]->c = current;
					}
					stack_size += numparents - 1;
				}
			}
		}
		//else
		//	break;
	}

	//!!!! Why is this needed? !!!!
	//Calculate the distance and number of upstream links for each link
	for(i=0;i<N;i++)	system[i]->distance = 0;

	for(i=0;i<leaves_size;i++)
	{
		prev = leaves[i];
		prev->distance = 1;
		for(current = prev->c; current != NULL; current = current->c)
		{
			if(current->distance > prev->distance+1)	break;
			else						current->distance = prev->distance + 1;
			prev = current;
		}
	}

	//Partition the system and assign the links
	*my_data = Initialize_TransData();
	*getting = (short int*) malloc(N*sizeof(short int));
	if(custom_model && custom_model->Partitioning)
		*assignments = custom_model->Partitioning(system,N,leaves,leaves_size,my_sys,my_N,*my_data,*getting);
	else
	{
		*assignments = Partition_System_By_Leaves(system,N,leaves,leaves_size,my_sys,my_N,*my_data,*getting);
		//*assignments = Partition_System_By_Leaves_2(system,N,leaves,leaves_size,my_sys,my_N,*my_data,*getting);
		//*assignments = Partition_METIS_Traditional(system,N,leaves,leaves_size,my_sys,my_N,*my_data,*getting,GlobalVars);
		//*assignments = Partition_METIS_RainChanges(system,N,leaves,leaves_size,my_sys,my_N,*my_data,*getting,GlobalVars);
		//*assignments = Partition_METIS_RainVolume(system,N,leaves,leaves_size,my_sys,my_N,*my_data,*getting,GlobalVars);	//!!!! Requires params for all links !!!!
	}

	//Clean up
	free(stack);
	free(leaves);
	//free(upstream_order);

	return 0;
}


//Reads numerical error tolerances. Builds RK methods.
//!!!! I'm not really sure how to handle specifying the dimension here. Should the rkd file allow a variable number of tols? !!!!
int Build_RKData(Link** system,char rk_filename[],unsigned int N,unsigned int* my_sys,unsigned int my_N,int* assignments,short int* getting,UnivVars* GlobalVars,ErrorData* GlobalErrors,RKMethod*** AllMethods,unsigned int* nummethods)
{
	unsigned int i,j,size,*link_ids,*methods;
	FILE* rkdata;
	double *filedata_abs,*filedata_rel,*filedata_abs_dense,*filedata_rel_dense;

	//Build all the RKMethods
	*nummethods = 4;
	*AllMethods = malloc(*nummethods * sizeof(RKMethod*));
	(*AllMethods)[0] = RKDense3_2();
	(*AllMethods)[1] = TheRKDense4_3();
	(*AllMethods)[2] = DOPRI5_dense();
	(*AllMethods)[3] = RadauIIA3_dense();
	GlobalVars->max_localorder = (*AllMethods)[0]->localorder;
	GlobalVars->max_s = (*AllMethods)[0]->s;
	for(i=1;i<*nummethods;i++)
	{
		GlobalVars->max_localorder = (GlobalVars->max_localorder < (*AllMethods)[i]->localorder) ? (*AllMethods)[i]->localorder : GlobalVars->max_localorder;
		GlobalVars->max_s = (GlobalVars->max_s > (*AllMethods)[i]->s) ? GlobalVars->max_s : (*AllMethods)[i]->s;
		//!!!! Note: Use a +1 for Radau solver? !!!!
	}

	if(rk_filename[0] != '\0')
	{
		link_ids = (unsigned int*) malloc(N*sizeof(unsigned int));
		methods = (unsigned int*) malloc(N*sizeof(unsigned int));

		if(my_rank == 0)
		{
			rkdata = fopen(rk_filename,"r");
			if(rkdata == NULL)
			{
				printf("Error: file %s not found for .rkd file\n",rk_filename);
				return 1;
			}
			if(CheckWinFormat(rkdata))
			{
				printf("Error: File %s appears to be in Windows format. Try converting to unix format using 'dos2unix' at the command line.\n",rk_filename);
				fclose(rkdata);
				return 1;
			}

			fscanf(rkdata,"%u %u",&i,&size);
			if(i != N)
			{
				printf("Error: the number of links in the rkd file differ from the number in the topology data (Got %u, expected %u).\n",i,N);
				return 1;
			}

			MPI_Bcast(&size,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
			filedata_abs = (double*) malloc(N*size*sizeof(double));
			filedata_rel = (double*) malloc(N*size*sizeof(double));
			filedata_abs_dense = (double*) malloc(N*size*sizeof(double));
			filedata_rel_dense = (double*) malloc(N*size*sizeof(double));
			methods = (unsigned int*) malloc(size*sizeof(unsigned int));

			//Read the file
			for(i=0;i<N;i++)
			{
				if(fscanf(rkdata,"%u",&(link_ids[i])) == 0)
				{
					printf("Error reading .rkd file: Not enough links in file (expected %u, got %u).\n",N,i);
					return 1;
				}
				for(j=0;j<size;i++)	fscanf(rkdata,"%lf",&(filedata_abs[i*size+j]));
				for(j=0;j<size;i++)	fscanf(rkdata,"%lf",&(filedata_rel[i*size+j]));
				for(j=0;j<size;i++)	fscanf(rkdata,"%lf",&(filedata_abs_dense[i*size+j]));
				for(j=0;j<size;i++)	fscanf(rkdata,"%lf",&(filedata_rel_dense[i*size+j]));
				fscanf(rkdata,"%u",&(methods[i]));
			}
		}
		else
		{
			MPI_Bcast(&size,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
			filedata_abs = (double*) malloc(N*size*sizeof(double));
			filedata_rel = (double*) malloc(N*size*sizeof(double));
			filedata_abs_dense = (double*) malloc(N*size*sizeof(double));
			filedata_rel_dense = (double*) malloc(N*size*sizeof(double));
		}

		//Broadcast data
		MPI_Bcast(link_ids,N,MPI_UNSIGNED,0,MPI_COMM_WORLD);
		MPI_Bcast(filedata_abs,N*size,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(filedata_rel,N*size,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(filedata_abs_dense,N*size,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(filedata_rel_dense,N*size,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(methods,N,MPI_UNSIGNED,0,MPI_COMM_WORLD);

		//Construct error data at each link
		for(i=0;i<N;i++)
		{
			if(assignments[i] == my_rank || getting[i])
			{
				system[i]->errorinfo = (ErrorData*) malloc(sizeof(ErrorData));
				system[i]->errorinfo->abstol = v_get(size);
				system[i]->errorinfo->reltol = v_get(size);
				system[i]->errorinfo->abstol_dense = v_get(size);
				system[i]->errorinfo->reltol_dense = v_get(size);
				system[i]->errorinfo->facmax = GlobalErrors->facmax;
				system[i]->errorinfo->facmin = GlobalErrors->facmin;
				system[i]->errorinfo->fac = GlobalErrors->fac;

				for(j=0;j<size;j++)
				{
					system[i]->errorinfo->abstol.ve[j] = filedata_abs[i*size+j];
					system[i]->errorinfo->reltol.ve[j] = filedata_rel[i*size+j];
					system[i]->errorinfo->abstol_dense.ve[j] = filedata_abs_dense[i*size+j];
					system[i]->errorinfo->reltol_dense.ve[j] = filedata_rel_dense[i*size+j];
				}
				system[i]->method = (*AllMethods)[methods[i]];
			}
		}
	}
	else
	{
		for(i=0;i<N;i++)
		{
			if(assignments[i] == my_rank || getting[i])
			{
				system[i]->errorinfo = GlobalErrors;
				system[i]->method = (*AllMethods)[GlobalVars->method];
			}
		}
	}

	return 0;
}

//Runs the init routine for the model. Also performs precalculations.
//Returns 0 if everything is ok, 1 if an error occurred.
int Initialize_Model(Link** system,unsigned int N,unsigned int* my_sys,unsigned int my_N,int* assignments,short int* getting,UnivVars* GlobalVars,model* custom_model,void* external)
{
	unsigned int i,j,max_dim = 0,smallest_dim;
	int my_error_code = 0,error_code;

	for(i=0;i<N;i++)
	{
		if(assignments[i] == my_rank || getting[i])
		{
			if(custom_model)
			{
				custom_model->Routines(system[i],GlobalVars->type,system[i]->method->exp_imp,system[i]->dam,external);
				custom_model->Precalculations(system[i],GlobalVars->global_params,system[i]->params,GlobalVars->disk_params,GlobalVars->params_size,system[i]->dam,GlobalVars->type,external);
			}
			else
			{
				InitRoutines(system[i],GlobalVars->type,system[i]->method->exp_imp,system[i]->dam,external);
				Precalculations(system[i],GlobalVars->global_params,system[i]->params,GlobalVars->disk_params,GlobalVars->params_size,system[i]->dam,GlobalVars->type,external);
			}

			max_dim = (max_dim < system[i]->dim) ? system[i]->dim : max_dim;

			//Be sure the problem dimension and number of error tolerances are compatible
			if(assignments[i] == my_rank)
			{
				smallest_dim = system[i]->errorinfo->abstol.dim;
				if(smallest_dim > system[i]->errorinfo->reltol.dim)
                    smallest_dim = system[i]->errorinfo->reltol.dim;
				if(smallest_dim > system[i]->errorinfo->abstol_dense.dim)
                    smallest_dim = system[i]->errorinfo->abstol_dense.dim;
				if(smallest_dim > system[i]->errorinfo->reltol_dense.dim)
                    smallest_dim = system[i]->errorinfo->reltol_dense.dim;
				if(GlobalVars->min_error_tolerances > smallest_dim)
				{
					printf("[%i] Error: link id %u does not have enough error tolerances (got %u, expected %u)\n",my_rank,system[i]->ID,smallest_dim, GlobalVars->min_error_tolerances);
					my_error_code = 1;
					break;
				}
			}
		}
	}

	//Check if an error occurred
	MPI_Allreduce(&my_error_code,&error_code,1,MPI_INT,MPI_LOR,MPI_COMM_WORLD);
	if(error_code)	return 1;

	//Make sure all procs know how large the problem is everywhere
	MPI_Allreduce(&max_dim,&(GlobalVars->max_dim),1,MPI_UNSIGNED,MPI_MAX,MPI_COMM_WORLD);
	for(i=0;i<N;i++)
		MPI_Bcast(&(system[i]->dim),1,MPI_UNSIGNED,assignments[i],MPI_COMM_WORLD);

	//Mix dense_indices with print_indices
	unsigned int loc,num_to_add;
	Link* current;
	unsigned int* states_to_add = (unsigned int*) malloc(GlobalVars->num_states_for_printing*sizeof(unsigned int));	//!!!! Only mix if save_flag set !!!!
	for(loc=0;loc<N;loc++)
	{
		if(assignments[loc] == my_rank || getting[loc])
		{
			current = system[loc];
			num_to_add = 0;

			for(i=0;i<GlobalVars->num_states_for_printing;i++)
			{
				if(GlobalVars->print_indices[i] > current->dim)	continue;	//State is not present at this link
				for(j=0;j<current->num_dense;j++)
					if(GlobalVars->print_indices[i] == current->dense_indices[j])	break;
				if(j == current->num_dense)
					states_to_add[num_to_add++] = GlobalVars->print_indices[i];
			}

			if(num_to_add)
			{
				current->dense_indices = (unsigned int*) realloc(current->dense_indices,(current->num_dense + num_to_add)*sizeof(unsigned int));
				for(i=0;i<num_to_add;i++)
					current->dense_indices[i+current->num_dense] = states_to_add[i];
				current->num_dense += num_to_add;
				merge_sort_1D(current->dense_indices,current->num_dense);
			}
		}
	}

	free(states_to_add);

	//Make sure all procs know the number of dense states at each link
	for(i=0;i<N;i++)
		MPI_Bcast(&(system[i]->num_dense),1,MPI_UNSIGNED,assignments[i],MPI_COMM_WORLD);

	return 0;
}


static int Load_Initial_Conditions_Ini(Link** system, unsigned int N, int* assignments, short int* getting, unsigned int** id_to_loc, UnivVars* GlobalVars, ConnData** db_connections, model* custom_model, void* external)
{
    unsigned int i, j, id, loc, no_ini_start, diff_start = 0, dim;
    FILE* initdata = NULL;
    short int *who_needs = NULL;
    short int my_need;
    VEC y_0 = v_get(0);

    //Proc 0 reads the file and sends the data to the other procs
    if (my_rank == 0)
    {
        initdata = fopen(GlobalVars->init_filename, "r");
        if (!initdata)
        {
            printf("Error: file %s not found for .ini file.\n", GlobalVars->init_filename);
            return 1;
        }
        if (CheckWinFormat(initdata))
        {
            printf("Error: File %s appears to be in Windows format. Try converting to unix format using 'dos2unix' at the command line.\n", GlobalVars->init_filename);
            fclose(initdata);
            return 1;
        }

        fscanf(initdata, "%*i %u %lf", &i, &(GlobalVars->t_0));	//Read model type, number of links, init time

        if (i != N)
        {
            printf("Error: the number of links in %s differs from the number in the topology data. (Got %u, expected %u)\n", GlobalVars->init_filename, i, N);
            return 1;
        }

        //Broadcast initial time
        MPI_Bcast(&(GlobalVars->t_0), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        //Read the .ini file
        who_needs = (short int*)malloc(np*sizeof(short int));
        for (i = 0; i < N; i++)
        {
            //Send current location
            fscanf(initdata, "%u", &id);
            loc = find_link_by_idtoloc(id, id_to_loc, N);
            if (loc > N)
            {
                printf("Error: link id %u in initial condition file, but not in network.\n", id);
                return 1;
            }
            MPI_Bcast(&loc, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

            //See who needs info about this link.
            //0 means the proc doesn't need it, 1 means link is assigned to proc, 2 means the link is a ghost to the proc.
            my_need = (assignments[loc] == my_rank) ? 1 : ((getting[loc]) ? 2 : 0);
            MPI_Gather(&my_need, 1, MPI_SHORT, who_needs, 1, MPI_SHORT, 0, MPI_COMM_WORLD);
            if (my_need == 1)
            {
                no_ini_start = system[loc]->no_ini_start;
                diff_start = system[loc]->diff_start;
                dim = system[loc]->dim;
            }
            else
            {
                MPI_Recv(&no_ini_start, 1, MPI_UNSIGNED, assignments[loc], 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&diff_start, 1, MPI_UNSIGNED, assignments[loc], 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&dim, 1, MPI_UNSIGNED, assignments[loc], 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            //Read init data
            v_resize(&y_0, dim);
            for (j = diff_start; j < no_ini_start; j++)
            {
                if (0 == fscanf(initdata, "%lf", &(y_0.ve[j])))
                {
                    printf("Error: not enough states in .ini file.\n");
                    return 1;
                }
            }

            //Send data to assigned proc and getting proc
            if (assignments[loc] == my_rank || getting[loc])
            {
                if (custom_model)
                    system[loc]->state = custom_model->InitializeEqs(GlobalVars->global_params, system[loc]->params, system[loc]->qvs, system[loc]->dam, y_0, GlobalVars->type, diff_start, no_ini_start, system[loc]->user, external);
                else
                    system[loc]->state = ReadInitData(GlobalVars->global_params, system[loc]->params, system[loc]->qvs, system[loc]->dam, y_0, GlobalVars->type, diff_start, no_ini_start, system[loc]->user, external);
                system[loc]->list = Create_List(y_0, GlobalVars->t_0, y_0.dim, system[loc]->num_dense, system[loc]->method->s, GlobalVars->iter_limit);
                system[loc]->list->head->state = system[loc]->state;
                system[loc]->last_t = GlobalVars->t_0;
            }

            if (assignments[loc] != my_rank)
                MPI_Send(&(y_0.ve[diff_start]), no_ini_start - diff_start, MPI_DOUBLE, assignments[loc], 2, MPI_COMM_WORLD);
            if (!(getting[loc]))
            {
                for (j = 0; j < np; j++)	if (who_needs[j] == 2)	break;
                if (j < np)
                    MPI_Send(&(y_0.ve[diff_start]), no_ini_start - diff_start, MPI_DOUBLE, (int)j, 2, MPI_COMM_WORLD);
            }
        }

        //Clean up
        fclose(initdata);
        free(who_needs);
    }
    else
    {
        //Get initial time
        MPI_Bcast(&(GlobalVars->t_0), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        who_needs = (short int*)malloc(np*sizeof(short int));
        for (i = 0; i < N; i++)
        {
            //Get link location
            MPI_Bcast(&loc, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

            //Is data needed for this link assigned at this proc?
            my_need = (assignments[loc] == my_rank) ? 1 : ((getting[loc]) ? 2 : 0);
            MPI_Gather(&my_need, 1, MPI_SHORT, who_needs, 1, MPI_SHORT, 0, MPI_COMM_WORLD);

            if (my_need)
            {
                no_ini_start = system[loc]->no_ini_start;
                diff_start = system[loc]->diff_start;
                dim = system[loc]->dim;

                if (assignments[loc] == my_rank)
                {
                    MPI_Send(&no_ini_start, 1, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD);
                    MPI_Send(&diff_start, 1, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD);
                    MPI_Send(&dim, 1, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD);	//!!!! Actually, this might be available everywhere now !!!!
                }

                v_resize(&y_0, dim);
                MPI_Recv(&(y_0.ve[diff_start]), no_ini_start - diff_start, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                if (custom_model)
                    system[loc]->state = custom_model->InitializeEqs(GlobalVars->global_params, system[loc]->params, system[loc]->qvs, system[loc]->dam, y_0, GlobalVars->type, diff_start, no_ini_start, system[loc]->user, external);
                else
                    system[loc]->state = ReadInitData(GlobalVars->global_params, system[loc]->params, system[loc]->qvs, system[loc]->dam, y_0, GlobalVars->type, diff_start, no_ini_start, system[loc]->user, external);
                system[loc]->list = Create_List(y_0, GlobalVars->t_0, y_0.dim, system[loc]->num_dense, system[loc]->method->s, GlobalVars->iter_limit);
                system[loc]->list->head->state = system[loc]->state;
                system[loc]->last_t = GlobalVars->t_0;
            }
        }

        free(who_needs);
    }

    //Clean up
    v_free(&y_0);

    return 0;
}

static int Load_Initial_Conditions_Uini(Link** system, unsigned int N, int* assignments, short int* getting, unsigned int** id_to_loc, UnivVars* GlobalVars, ConnData** db_connections, model* custom_model, void* external)
{
    unsigned int i, j, loc, no_ini_start, diff_start = 0;
    FILE* initdata = NULL;
    short int *who_needs = NULL;
    VEC y_0 = v_get(0);

    //Proc 0 reads the initial conds, and send them to the other procs
    if (my_rank == 0)
    {
        initdata = fopen(GlobalVars->init_filename, "r");
        if (!initdata)
        {
            printf("Error: file %s not found for .uini file.\n", GlobalVars->init_filename);
            return 1;
        }
        if (CheckWinFormat(initdata))
        {
            printf("Error: File %s appears to be in Windows format. Try converting to unix format using 'dos2unix' at the command line.\n", GlobalVars->init_filename);
            fclose(initdata);
            return 1;
        }

        fscanf(initdata, "%*i %lf", &(GlobalVars->t_0));	//Read model type, init time
    }

    //Broadcast the initial time
    MPI_Bcast(&(GlobalVars->t_0), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //Get number of values to read from disk (and error checking)
    for (i = 0; i<N; i++)
    {
        if (assignments[i] == my_rank)
        {
            loc = i;
            no_ini_start = system[loc]->no_ini_start;
            diff_start = system[loc]->diff_start;
            //y_0.dim = system[loc]->dim;
            break;
        }
    }
    for (; i<N; i++)
    {
        //if(assignments[i] == my_rank && y_0.dim != system[i]->dim)
        if (assignments[i] == my_rank && (no_ini_start != system[i]->no_ini_start || diff_start != system[i]->diff_start))
        {
            printf("[%i]: Error: model type %u does not support .uini files (because a variable number of states must be specified for the initial conditions).\n", my_rank, GlobalVars->type);
            return 1;
        }
    }

    //no_ini_start = system[loc]->no_ini_start;
    //diff_start = system[loc]->diff_start;
    VEC y_0_backup = v_get(no_ini_start - diff_start);
    //y_0_backup->dim = no_ini_start - diff_start;
    //y_0_backup.ve = (double*) calloc(y_0_backup->dim,sizeof(double));

    if (my_rank == 0)
    {
        //for(i=diff_start;i<no_ini_start;i++)
        for (i = 0; i<y_0_backup.dim; i++)
        {
            if (fscanf(initdata, "%lf", &(y_0_backup.ve[i])) == 0)
            {
                printf("Error reading .uini file: Not enough initial states.\n");
                return 1;
            }
        }

        //Done with file, so close it
        fclose(initdata);
    }

    MPI_Bcast(y_0_backup.ve, no_ini_start - diff_start, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //VEC* y_0_backup = v_get(y_0.dim);
    //v_copy(y_0,y_0_backup);

    //Store the data
    for (i = 0; i<N; i++)
    {
        if (assignments[i] == my_rank || getting[i])
        {
            v_resize(&y_0, system[i]->dim);
            for (j = diff_start; j<no_ini_start; j++)
                y_0.ve[j] = y_0_backup.ve[j - diff_start];

            if (custom_model)
                system[i]->state = custom_model->InitializeEqs(GlobalVars->global_params, system[i]->params, system[i]->qvs, system[i]->dam, y_0, GlobalVars->type, diff_start, no_ini_start, system[i]->user, external);
            else
                system[i]->state = ReadInitData(GlobalVars->global_params, system[i]->params, system[i]->qvs, system[i]->dam, y_0, GlobalVars->type, diff_start, no_ini_start, system[i]->user, external);
            system[i]->list = Create_List(y_0, GlobalVars->t_0, system[i]->dim, system[i]->num_dense, system[i]->method->s, GlobalVars->iter_limit);
            system[i]->list->head->state = system[i]->state;
            system[i]->last_t = GlobalVars->t_0;
            //v_copy(y_0_backup,y_0);
        }
    }

    //Clean up
    v_free(&y_0);
    v_free(&y_0_backup);

    return 0;
}

static int Load_Initial_Conditions_Rec(Link** system, unsigned int N, int* assignments, short int* getting, unsigned int** id_to_loc, UnivVars* GlobalVars, ConnData** db_connections, model* custom_model, void* external)
{
    unsigned int i, j, id, loc, dim;
    FILE* initdata = NULL;
    short int *who_needs = NULL;
    short int my_need;
    VEC y_0 = v_get(0);

    //Proc 0 reads the file and sends the data to the other procs
    if (my_rank == 0)
    {
        initdata = fopen(GlobalVars->init_filename, "r");
        if (!initdata)
        {
            printf("Error: file %s not found for .rec file.\n", GlobalVars->init_filename);
            return 1;
        }
        if (CheckWinFormat(initdata))
        {
            printf("Error: File %s appears to be in Windows format. Try converting to unix format using 'dos2unix' at the command line.\n", GlobalVars->init_filename);
            fclose(initdata);
            return 1;
        }

        fscanf(initdata, "%*i %u %lf", &i, &(GlobalVars->t_0));	//Read model type, number of links, init time

        if (i != N)
        {
            printf("Error: the number of links in %s differs from the number in the topology data. (Got %u, expected %u)\n", GlobalVars->init_filename, i, N);
            return 1;
        }

        //Broadcast the initial time
        MPI_Bcast(&(GlobalVars->t_0), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        //Read the .rec file
        who_needs = (short int*)malloc(np*sizeof(short int));
        for (i = 0; i<N; i++)
        {
            //Send current location
            fscanf(initdata, "%u", &id);
            loc = find_link_by_idtoloc(id, id_to_loc, N);
            if (loc > N)
            {
                printf("Error: link id %u in initial condition file, but not in network.\n", id);
                return 1;
            }
            MPI_Bcast(&loc, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

            //See who needs info about this link.
            //0 means the proc doesn't need it, 1 means link is assigned to proc, 2 means the link is a ghost to the proc.
            my_need = (assignments[loc] == my_rank) ? 1 : ((getting[loc]) ? 2 : 0);
            MPI_Gather(&my_need, 1, MPI_SHORT, who_needs, 1, MPI_SHORT, 0, MPI_COMM_WORLD);
            if (my_need == 1)
                dim = system[loc]->dim;
            else
                MPI_Recv(&dim, 1, MPI_UNSIGNED, assignments[loc], 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            //Read init data
            v_resize(&y_0, dim);
            for (j = 0; j<y_0.dim; j++)
            {
                if (0 == fscanf(initdata, "%lf", &(y_0.ve[j])))
                {
                    printf("Error: not enough states in .rec file.\n");
                    return 1;
                }
            }

            //Send data to assigned proc and getting proc
            if (assignments[loc] == my_rank || getting[loc])
            {
                if (system[loc]->state_check)
                    system[loc]->state = system[loc]->state_check(y_0, GlobalVars->global_params, system[loc]->params, system[loc]->qvs, system[loc]->dam);
                system[loc]->list = Create_List(y_0, GlobalVars->t_0, y_0.dim, system[loc]->num_dense, system[loc]->method->s, GlobalVars->iter_limit);
                system[loc]->list->head->state = system[loc]->state;
                system[loc]->last_t = GlobalVars->t_0;
            }

            if (assignments[loc] != my_rank)
                MPI_Send(y_0.ve, y_0.dim, MPI_DOUBLE, assignments[loc], 2, MPI_COMM_WORLD);
            if (!(getting[loc]))
            {
                for (j = 0; j<np; j++)	if (who_needs[j] == 2)	break;
                if (j < np)
                    MPI_Send(y_0.ve, y_0.dim, MPI_DOUBLE, (int)j, 2, MPI_COMM_WORLD);
            }
        }

        //Clean up
        fclose(initdata);
        free(who_needs);
    }
    else
    {
        //Get the initial time
        MPI_Bcast(&(GlobalVars->t_0), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        for (i = 0; i<N; i++)
        {
            //Get link location
            MPI_Bcast(&loc, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

            //Is data needed for this link assigned at this proc?
            my_need = (assignments[loc] == my_rank) ? 1 : ((getting[loc]) ? 2 : 0);
            MPI_Gather(&my_need, 1, MPI_SHORT, who_needs, 1, MPI_SHORT, 0, MPI_COMM_WORLD);

            if (my_need)
            {
                dim = system[loc]->dim;

                if (assignments[loc] == my_rank)
                    MPI_Send(&dim, 1, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD);

                v_resize(&y_0, dim);
                MPI_Recv(y_0.ve, y_0.dim, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                if (system[loc]->state_check)
                    system[loc]->state = system[loc]->state_check(y_0, GlobalVars->global_params, system[loc]->params, system[loc]->qvs, system[loc]->dam);
                system[loc]->list = Create_List(y_0, GlobalVars->t_0, y_0.dim, system[loc]->num_dense, system[loc]->method->s, GlobalVars->iter_limit);
                system[loc]->list->head->state = system[loc]->state;
                system[loc]->last_t = GlobalVars->t_0;
            }
        }
    }

    //Clean up
    v_free(&y_0);

    return 0;
}


static int Load_Initial_Conditions_Dbc(Link** system, unsigned int N, int* assignments, short int* getting, unsigned int** id_to_loc, UnivVars* GlobalVars, ConnData** db_connections, model* custom_model, void* external)
{
    unsigned int i, j, loc;
    short int *who_needs = NULL;
    short int my_need;
    VEC y_0 = v_get(0);
    PGresult *res;

    //!!!! Note: this assumes the database is like a .rec, with each state given. It also !!!!
    //!!!! assumes the same number of states at each link. !!!!

    //Set t_0 (I'm not really sure what else to do here...)
    GlobalVars->t_0 = 0.0;

    if (my_rank == 0)
    {
        //Download data
        if (ConnectPGDB(db_connections[ASYNCH_DB_LOC_INIT]))
        {
            printf("Error connecting to database for init conditions.\n");
            return 1;
        }
        sprintf(db_connections[ASYNCH_DB_LOC_INIT]->query, db_connections[ASYNCH_DB_LOC_INIT]->queries[0], GlobalVars->init_timestamp);
        res = PQexec(db_connections[ASYNCH_DB_LOC_INIT]->conn, db_connections[ASYNCH_DB_LOC_INIT]->query);
        if (CheckResError(res, "downloading init data"))	return 1;

        if (PQntuples(res) != N)
        {
            printf("Error downloading init data. Got %i conditions, expected %u.\n", PQntuples(res), N);
            return 1;
        }

        //Get dim
        y_0.dim = PQnfields(res) - 1;
        y_0.ve = (double*)realloc(y_0.ve, y_0.dim*sizeof(double));
        MPI_Bcast(&(y_0.dim), 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        who_needs = (short int*)malloc(np*sizeof(short int));

        //Read data
        for (i = 0; i<N; i++)
        {
            loc = find_link_by_idtoloc(atoi(PQgetvalue(res, i, 0)), id_to_loc, N);
            for (j = 0; j<y_0.dim; j++)	y_0.ve[j] = atof(PQgetvalue(res, i, j + 1));
            MPI_Bcast(&loc, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

            //See who needs info about this link.
            //0 means the proc doesn't need it, 1 means link is assigned to proc, 2 means the link is a ghost to the proc.
            my_need = (assignments[loc] == my_rank) ? 1 : ((getting[loc]) ? 2 : 0);
            MPI_Gather(&my_need, 1, MPI_SHORT, who_needs, 1, MPI_SHORT, 0, MPI_COMM_WORLD);

            //Send the data
            if (assignments[loc] == my_rank || getting[loc])
            {
                if (system[loc]->state_check != NULL)
                    system[loc]->state = system[loc]->state_check(y_0, GlobalVars->global_params, system[loc]->params, system[loc]->qvs, system[loc]->dam);
                system[loc]->list = Create_List(y_0, GlobalVars->t_0, y_0.dim, system[loc]->num_dense, system[loc]->method->s, GlobalVars->iter_limit);
                system[loc]->list->head->state = system[loc]->state;
                system[loc]->last_t = GlobalVars->t_0;
            }

            if (assignments[loc] != my_rank)
                MPI_Send(y_0.ve, y_0.dim, MPI_DOUBLE, assignments[loc], 2, MPI_COMM_WORLD);
            if (!(getting[loc]))
            {
                for (j = 0; j<np; j++)	if (who_needs[j] == 2)	break;
                if (j < np)	MPI_Send(y_0.ve, y_0.dim, MPI_DOUBLE, (int)j, 2, MPI_COMM_WORLD);
            }
        }

        //Clean up
        PQclear(res);
        DisconnectPGDB(db_connections[ASYNCH_DB_LOC_INIT]);
        v_free(&y_0);
        free(who_needs);
    }
    else
    {
        //Get dim
        MPI_Bcast(&(y_0.dim), 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        y_0.ve = (double*)realloc(y_0.ve, y_0.dim*sizeof(double));

        for (i = 0; i<N; i++)
        {
            //Get link location
            MPI_Bcast(&loc, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

            //Is data needed for this link assigned at this proc?
            my_need = (assignments[loc] == my_rank) ? 1 : ((getting[loc]) ? 2 : 0);
            MPI_Gather(&my_need, 1, MPI_SHORT, who_needs, 1, MPI_SHORT, 0, MPI_COMM_WORLD);

            if (my_need)
            {
                MPI_Recv(y_0.ve, y_0.dim, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                if (system[loc]->state_check != NULL)
                    system[loc]->state = system[loc]->state_check(y_0, GlobalVars->global_params, system[loc]->params, system[loc]->qvs, system[loc]->dam);
                system[loc]->list = Create_List(y_0, GlobalVars->t_0, y_0.dim, system[loc]->num_dense, system[loc]->method->s, GlobalVars->iter_limit);
                system[loc]->list->head->state = system[loc]->state;
                system[loc]->last_t = GlobalVars->t_0;
            }
        }

        //Clean up
        v_free(&y_0);
    }

    return 0;
}


static int Load_Initial_Conditions_H5(Link** system, unsigned int N, int* assignments, short int* getting, unsigned int** id_to_loc, UnivVars* GlobalVars, ConnData** db_connections, model* custom_model, void* external)
{
    unsigned char who_needs[ASYNCH_MAX_NUMBER_OF_PROCESS];
    unsigned char my_need;

    //Proc 0 reads the file and sends the data to the other procs
    if (my_rank == 0)
    {
        hid_t file_id = H5Fopen(GlobalVars->init_filename, H5F_ACC_RDONLY, H5P_DEFAULT);
        if (!file_id)
        {
            printf("Error: file %s not found for .h5 file.\n", GlobalVars->init_filename);
            return 1;
        }

        hsize_t dims[2];
        H5LTget_dataset_info(file_id, "/state", dims, NULL, NULL);

        if (dims[0] != N)
        {
            printf("Error: the number of links in %s differs from the number in the topology data. (Got %llu, expected %u)\n", GlobalVars->init_filename, dims[1], N);
            return 1;
        }

        //Assume that every links have the same dimension
        unsigned int dim = system[0]->dim;

        //Read model type, init time
        unsigned short type;
        H5LTget_attribute_ushort(file_id, "/", "model", &type);

        if (type != GlobalVars->type)
        {
            printf("Error: model type do no match. (Got %hu, expected %hu)\n", GlobalVars->type, type);
            return 1;
        }

        //Broadcast the initial time
        MPI_Bcast(&(GlobalVars->t_0), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        int *index = malloc(dims[0] * sizeof(unsigned int));
        double *data = malloc(dims[0] * dims[1] * sizeof(double));

        H5LTread_dataset_int(file_id, "/index", index);
        H5LTread_dataset_double(file_id, "/state", data);

        //Read the .h5 file
        for (unsigned int i = 0; i<N; i++)
        {
            //Send current location
            unsigned int id = index[i];
            unsigned int loc = find_link_by_idtoloc(id, id_to_loc, N);
            if (loc > N)
            {
                printf("Error: link id %u in initial condition file, but not in network.\n", id);
                return 1;
            }
            MPI_Bcast(&loc, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

            //See who needs info about this link.
            //0 means the proc doesn't need it, 1 means link is assigned to proc, 2 means the link is a ghost to the proc.
            my_need = (assignments[loc] == my_rank) ? 1 : ((getting[loc]) ? 2 : 0);
            MPI_Gather(&my_need, 1, MPI_UNSIGNED_CHAR, who_needs, 1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

            unsigned int dim;
            if (my_need == 1)
                dim = system[loc]->dim;
            else
                MPI_Recv(&dim, 1, MPI_UNSIGNED, assignments[loc], 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            //Read init data
            VEC y_0 = v_get(dim);
            memcpy(y_0.ve, &data[i * dim], dim * sizeof(double));

            //Send data to assigned proc and getting proc
            if (assignments[loc] == my_rank || getting[loc])
            {
                if (system[loc]->state_check)
                    system[loc]->state = system[loc]->state_check(y_0, GlobalVars->global_params, system[loc]->params, system[loc]->qvs, system[loc]->dam);
                system[loc]->list = Create_List(y_0, GlobalVars->t_0, y_0.dim, system[loc]->num_dense, system[loc]->method->s, GlobalVars->iter_limit);
                system[loc]->list->head->state = system[loc]->state;
                system[loc]->last_t = GlobalVars->t_0;
            }

            if (assignments[loc] != my_rank)
                MPI_Send(y_0.ve, y_0.dim, MPI_DOUBLE, assignments[loc], 2, MPI_COMM_WORLD);
            if (!(getting[loc]))
            {
                unsigned int j;
                for (j = 0; j<np; j++)	if (who_needs[j] == 2)	break;
                if (j < np)
                    MPI_Send(y_0.ve, y_0.dim, MPI_DOUBLE, (int)j, 2, MPI_COMM_WORLD);
            }

            //Clean up
            v_free(&y_0);
        }

        //Clean up
        H5Fclose(file_id);
    }
    else
    {
        //Get the initial time
        MPI_Bcast(&(GlobalVars->t_0), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        for (unsigned int i = 0; i<N; i++)
        {
            //Get link location
            unsigned int loc;
            MPI_Bcast(&loc, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

            //Is data needed for this link assigned at this proc?
            my_need = (assignments[loc] == my_rank) ? 1 : ((getting[loc]) ? 2 : 0);
            MPI_Gather(&my_need, 1, MPI_UNSIGNED_CHAR, who_needs, 1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

            if (my_need)
            {
                unsigned int dim = system[loc]->dim;

                if (assignments[loc] == my_rank)
                    MPI_Send(&dim, 1, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD);

                VEC y_0 = v_get(dim);
                MPI_Recv(y_0.ve, y_0.dim, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                if (system[loc]->state_check)
                    system[loc]->state = system[loc]->state_check(y_0, GlobalVars->global_params, system[loc]->params, system[loc]->qvs, system[loc]->dam);
                system[loc]->list = Create_List(y_0, GlobalVars->t_0, y_0.dim, system[loc]->num_dense, system[loc]->method->s, GlobalVars->iter_limit);
                system[loc]->list->head->state = system[loc]->state;
                system[loc]->last_t = GlobalVars->t_0;

                //Clean up
                v_free(&y_0);
            }
        }
    }

    return 0;
}


//Loads the initial conditions.
//Initial state (0 = .ini, 1 = .uini, 2 = .rec, 3 = .dbc)
int Load_Initial_Conditions(Link** system,unsigned int N,int* assignments,short int* getting,unsigned int** id_to_loc,UnivVars* GlobalVars,ConnData** db_connections,model* custom_model,void* external)
{
    int res = 0;

    switch (GlobalVars->init_flag)
    {
    //.ini
    case 0: 
        res = Load_Initial_Conditions_Ini(system, N, assignments, getting, id_to_loc, GlobalVars, db_connections, custom_model, external);
        break;
        
    //.uini
    case 1:
        res = Load_Initial_Conditions_Uini(system, N, assignments, getting, id_to_loc, GlobalVars, db_connections, custom_model, external);
        break;
    
    //.rec
    case 2:
        res = Load_Initial_Conditions_Rec(system, N, assignments, getting, id_to_loc, GlobalVars, db_connections, custom_model, external);
        break;
    
    //.dbc
    case 3:
        res = Load_Initial_Conditions_Dbc(system, N, assignments, getting, id_to_loc, GlobalVars, db_connections, custom_model, external);
        break;

    //.h5
    case 4:
        res = Load_Initial_Conditions_H5(system, N, assignments, getting, id_to_loc, GlobalVars, db_connections, custom_model, external);
        break;

    default:
        printf("Error: invalid intial condition file type.\n");
        res = -1;
    }

    return res;
}


//Loads the forcing data specified in the global file.
int Load_Forcings(Link** system,unsigned int N,unsigned int* my_sys,unsigned int my_N,int* assignments,short int* getting,unsigned int* res_list,unsigned int res_size,unsigned int** id_to_loc,UnivVars* GlobalVars,Forcing** forcings,ConnData** db_connections)
{
	unsigned int i,j,l,m,limit,id,loc;
	FILE* forcingfile = NULL;
	double *buffer = NULL,univ_forcing_change_time[ASYNCH_MAX_NUM_FORCINGS];
	PGresult *res;
	Link* current;

	//Reserve space for forcing data
	for(i=0;i<my_N;i++)
	{
		current = system[my_sys[i]];
		current->forcing_buff = (ForcingData**) malloc(GlobalVars->num_forcings*sizeof(ForcingData*));
		current->forcing_values = (double*) calloc(GlobalVars->num_forcings,sizeof(double));
		current->forcing_change_times = (double*) calloc(GlobalVars->num_forcings,sizeof(double));
		current->forcing_indices = (unsigned int*) malloc(GlobalVars->num_forcings*sizeof(double));
	}

	//Setup forcings. Read uniform forcing data and open .str files. Also initialize rainfall from database.
	for(l=0;l<GlobalVars->num_forcings;l++)
	{
		forcings[l]->maxtime = GlobalVars->t_0;
		forcings[l]->iteration = 0;
		forcings[l]->active = 1;

		//Go through each possible flag
		if(forcings[l]->flag == 0)	//0 forcing
		{
			//Set routines
			forcings[l]->GetPasses = &PassesOther;
			forcings[l]->GetNextForcing = &NextForcingOther;

			//Setup buffers at each link
			for(i=0;i<N;i++)
			{
				if(assignments[i] == my_rank)
				{
					system[i]->forcing_buff[l] = NULL;
					system[i]->forcing_values[l] = 0.0;
					system[i]->forcing_change_times[l] = GlobalVars->maxtime + 1.0;
				}
			}
		}
		else if(forcings[l]->flag == 1)	//Storm file
		{
			//Set routines
			forcings[l]->GetPasses = &PassesOther;
			forcings[l]->GetNextForcing = &NextForcingOther;

			if(my_rank == 0)
			{
				//Open .str file
				forcingfile = fopen(forcings[l]->filename,"r");
				if(!forcingfile)
				{
					printf("Error: cannot open forcing file %s.\n",forcings[l]->filename);
					return 1;
				}
				if(CheckWinFormat(forcingfile))
				{
					printf("Error: File %s appears to be in Windows format. Try converting to unix format using 'dos2unix' at the command line.\n",forcings[l]->filename);
					fclose(forcingfile);
					return 1;
				}

				fscanf(forcingfile,"%u",&limit);
				if(limit != N && (!(GlobalVars->res_flag) || l != GlobalVars->res_forcing_idx))
				{
					printf("Error: Number of links in .str file differs from number of links in network (%u vs %u).\n",limit,N);
					return 1;
				}
			}

			//Get total number of links with data
			MPI_Bcast(&limit,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);

			//Setup buffers at each link
			for(i=0;i<limit;i++)
			{
				if(my_rank == 0)
				{
					//Get location
					fscanf(forcingfile,"%i",&id);
					loc = find_link_by_idtoloc(id,id_to_loc,N);
					if(loc >= N)
					{
						printf("Error: forcing data provided for link id %u in forcing %u, but link id is not in the network.\n",m,l);
						return 1;
					}

					//Read values
					fscanf(forcingfile,"%i",&m);
					m++;		//Increase this by one to add a "ceiling" term
					buffer = (double*) realloc(buffer,2*m*sizeof(double));
					for(j=0;j<m-1;j++)
						fscanf(forcingfile,"%lf %lf",&(buffer[2*j]),&(buffer[2*j+1]));
					buffer[2*j] = GlobalVars->maxtime + 3.0;
					buffer[2*j+1] = -1.0;

					//Send data
					MPI_Bcast(&loc,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
					if(assignments[loc] != my_rank)
					{
						MPI_Send(&m,1,MPI_UNSIGNED,assignments[loc],1,MPI_COMM_WORLD);
						MPI_Send(buffer,2*m,MPI_DOUBLE,assignments[loc],1,MPI_COMM_WORLD);
					}
				}
				else
				{
					MPI_Bcast(&loc,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
					if(assignments[loc] == my_rank)
					{
						MPI_Recv(&m,1,MPI_UNSIGNED,0,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
						buffer = (double*) realloc(buffer,2*m*sizeof(double));
						MPI_Recv(buffer,2*m,MPI_DOUBLE,0,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					}
				}

				if(assignments[loc] == my_rank)
				{
					if(!(GlobalVars->res_flag) || !(l == GlobalVars->res_forcing_idx) || system[loc]->res)
					{
						system[loc]->forcing_buff[l] = (ForcingData*) malloc(sizeof(ForcingData));
						system[loc]->forcing_buff[l]->rainfall = (double**) malloc(m*sizeof(double*));
						system[loc]->forcing_buff[l]->n_times = m;
						for(j=0;j<m;j++)	system[loc]->forcing_buff[l]->rainfall[j] = (double*) malloc(2*sizeof(double));

						//Read in the storm data for this link
						for(j=0;j<m;j++)
						{
							system[loc]->forcing_buff[l]->rainfall[j][0] = buffer[2*j];
							system[loc]->forcing_buff[l]->rainfall[j][1] = buffer[2*j+1];
						}

						double rainfall_buffer = system[loc]->forcing_buff[l]->rainfall[0][1];
						system[loc]->forcing_values[l] = rainfall_buffer;
						system[loc]->forcing_indices[l] = 0;
						for(j=1;j<system[loc]->forcing_buff[l]->n_times;j++)
						{
							if(fabs(system[loc]->forcing_buff[l]->rainfall[j][1] - rainfall_buffer) > 1e-8)
							{
								system[loc]->forcing_change_times[l] = system[loc]->forcing_buff[l]->rainfall[j][0];
								break;
							}
						}
						if(j == system[loc]->forcing_buff[l]->n_times)
						{
							system[loc]->forcing_change_times[l] = system[loc]->forcing_buff[l]->rainfall[j-1][0];
							system[loc]->forcing_indices[l] = j-1;
						}
					}
					else	//No reservoir here
					{
						m = 2;	//Init value (assumed 0.0)
						
						system[loc]->forcing_buff[l] = (ForcingData*) malloc(sizeof(ForcingData));
						system[loc]->forcing_buff[l]->rainfall = (double**) malloc(m*sizeof(double*));
						system[loc]->forcing_buff[l]->n_times = m;
						for(j=0;j<m;j++)	system[loc]->forcing_buff[l]->rainfall[j] = (double*) malloc(2*sizeof(double));

						system[loc]->forcing_buff[l]->rainfall[0][0] = GlobalVars->t_0;
						system[loc]->forcing_buff[l]->rainfall[0][1] = 0.0;
						system[loc]->forcing_buff[l]->rainfall[1][0] = GlobalVars->maxtime + 3.0;
						system[loc]->forcing_buff[l]->rainfall[1][1] = -1.0;

						double rainfall_buffer = system[loc]->forcing_buff[l]->rainfall[0][1];
						system[loc]->forcing_values[l] = rainfall_buffer;
						system[loc]->forcing_indices[l] = 0;
						for(j=1;j<system[loc]->forcing_buff[l]->n_times;j++)
						{
							if(fabs(system[loc]->forcing_buff[l]->rainfall[j][1] - rainfall_buffer) > 1e-8)
							{
								system[loc]->forcing_change_times[l] = system[loc]->forcing_buff[l]->rainfall[j][0];
								break;
							}
						}
						if(j == system[loc]->forcing_buff[l]->n_times)
							system[loc]->forcing_change_times[l] = system[loc]->forcing_buff[l]->rainfall[j-1][0];
					}
				} //if(assigned to me)
			}

			//Allocate space for links without a reservoir
			if(GlobalVars->res_flag && l == GlobalVars->res_forcing_idx)
			{
				for(i=0;i<N;i++)
				{
					if(assignments[i] == my_rank && !(system[i]->forcing_buff[l]))
					{
						m = 2;	//Init value (assumed 0.0)

						system[i]->forcing_buff[l] = (ForcingData*) malloc(sizeof(ForcingData));
						system[i]->forcing_buff[l]->rainfall = (double**) malloc(m*sizeof(double*));
						system[i]->forcing_buff[l]->n_times = m;
						for(j=0;j<m;j++)	system[i]->forcing_buff[l]->rainfall[j] = (double*) malloc(2*sizeof(double));

						system[i]->forcing_buff[l]->rainfall[0][0] = GlobalVars->t_0;
						system[i]->forcing_buff[l]->rainfall[0][1] = 0.0;
						system[i]->forcing_buff[l]->rainfall[1][0] = GlobalVars->maxtime + 3.0;
						system[i]->forcing_buff[l]->rainfall[1][1] = -1.0;

						double rainfall_buffer = system[i]->forcing_buff[l]->rainfall[0][1];
						system[i]->forcing_values[l] = rainfall_buffer;
						system[i]->forcing_indices[l] = 0;
						for(j=1;j<system[i]->forcing_buff[l]->n_times;j++)
						{
							if(fabs(system[i]->forcing_buff[l]->rainfall[j][1] - rainfall_buffer) > 1e-8)
							{
								system[i]->forcing_change_times[l] = system[i]->forcing_buff[l]->rainfall[j][0];
								break;
							}
						}
						if(j == system[i]->forcing_buff[l]->n_times)
							system[i]->forcing_change_times[l] = system[i]->forcing_buff[l]->rainfall[j-1][0];
					}
				}
			}

			//Clean up
			if(forcingfile)	fclose(forcingfile);
		}
		else if(forcings[l]->flag == 2)	//Binary files
		{
			//Set routines
			forcings[l]->GetPasses = &PassesBinaryFiles;
			forcings[l]->GetNextForcing = &NextForcingBinaryFiles;

			//Setup buffers at each link
			//!!!! This should really be improved... !!!!
			for(i=0;i<N;i++)
			{
				if(assignments[i] == my_rank)
					system[i]->forcing_buff[l] = NULL;
			}
		}
		else if(forcings[l]->flag == 6)	//GZ binary files
		{
			//Set routines
			forcings[l]->GetPasses = &PassesBinaryFiles;
			forcings[l]->GetNextForcing = &NextForcingGZBinaryFiles;

			//Setup buffers at each link
			//!!!! This should really be improved... !!!!
			for(i=0;i<N;i++)
			{
				if(assignments[i] == my_rank)
					system[i]->forcing_buff[l] = NULL;
			}
		}
		else if(forcings[l]->flag == 8)	//Grid cell based files
		{
			//Set routines
			forcings[l]->GetPasses = &PassesBinaryFiles;
			forcings[l]->GetNextForcing = &NextForcingGridCell;

			//Setup buffers at each link
			//!!!! This should really be improved... !!!!
			for(i=0;i<N;i++)
			{
				if(assignments[i] == my_rank)
					system[i]->forcing_buff[l] = NULL;
			}

			//Read the index file
			if(my_rank == 0)
			{
				forcingfile = fopen(forcings[l]->filename,"r");
				if(!forcingfile)
				{
					printf("Error: forcing file %s not found.\n",forcings[l]->filename);
					return 1;
				}

				char* tempspace1 = (char*) malloc(1024*sizeof(char));
				char* tempspace2 = (char*) malloc(1024*sizeof(char));
				forcings[l]->lookup_filename = (char*) malloc(1024*sizeof(char));
				forcings[l]->fileident = (char*) malloc(1024*sizeof(char));
				FindPath(forcings[l]->filename,forcings[l]->fileident);
				i = strlen(forcings[l]->fileident);
				strcpy(forcings[l]->lookup_filename,forcings[l]->fileident);
				//fscanf(forcingfile,"%lf %lf %u %s %s",&(forcings[l]->file_time),&(forcings[l]->factor),&(forcings[l]->num_cells),&(forcings[l]->fileident[i]),&(forcings[l]->lookup_filename[i]));
				fscanf(forcingfile,"%lf %lf %u %s %s",&(forcings[l]->file_time),&(forcings[l]->factor),&(forcings[l]->num_cells),tempspace1,tempspace2);
				fclose(forcingfile);

				//Add path from index file?
				if(tempspace1[0] == '/')	//Use absolute path
					strcpy(forcings[l]->fileident,tempspace1);
				else		//Relative path
					strcpy(&(forcings[l]->fileident[i]),tempspace1);
				if(tempspace2[0] == '/')	//Use absolute path
					strcpy(forcings[l]->lookup_filename,tempspace2);
				else		//Relative path
					strcpy(&(forcings[l]->lookup_filename[i]),tempspace2);

				//Process the lookup file
				forcingfile = fopen(forcings[l]->lookup_filename,"r");
				if(!forcingfile)
				{
					printf("Error: lookup file %s not found.\n",forcings[l]->lookup_filename);
					return 1;
				}
				forcings[l]->grid_to_linkid = (unsigned int**) malloc(forcings[l]->num_cells*sizeof(unsigned int*));
				forcings[l]->num_links_in_grid = (unsigned int*) calloc(forcings[l]->num_cells,sizeof(unsigned int));

				while(!feof(forcingfile))	//Count the number of links in each grid cell
				{
					if(!fscanf(forcingfile,"%*u %u",&j))	break;
					(forcings[l]->num_links_in_grid[j])++;
				}

				rewind(forcingfile);
				for(i=0;i<forcings[l]->num_cells;i++)
					forcings[l]->grid_to_linkid[i] = (unsigned int*) malloc(forcings[l]->num_links_in_grid[i]*sizeof(unsigned int));
				unsigned int* counters = (unsigned int*) calloc(forcings[l]->num_cells,sizeof(unsigned int));

				while(!feof(forcingfile))	//Read in the link ids
				{
					if(!fscanf(forcingfile,"%u %u",&id,&j))	break;
					forcings[l]->grid_to_linkid[j][counters[j]++] = id;
				}

				fclose(forcingfile);
				free(counters);
				free(tempspace1);

				//Check if a grid cell file actually exists
				for(i=forcings[l]->first_file;i<forcings[l]->last_file;i++)
				{
					sprintf(tempspace2,"%s%u",forcings[l]->fileident,i);
					forcingfile = fopen(tempspace2,"rb");
					if(forcingfile)
					{
						fclose(forcingfile);
						break;
					}
				}
				if(i == forcings[l]->last_file)
					printf("Warning: No forcing files found at %s for forcing %u. Be sure this is the correct directory.\n",forcings[l]->fileident,l);

				free(tempspace2);
			}

			//Broadcast data
			MPI_Bcast(&(forcings[l]->file_time),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			MPI_Bcast(&(forcings[l]->factor),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			MPI_Bcast(&(forcings[l]->num_cells),1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
			if(my_rank != 0)
			{
				forcings[l]->grid_to_linkid = (unsigned int**) malloc(forcings[l]->num_cells*sizeof(unsigned int*));
				forcings[l]->num_links_in_grid = (unsigned int*) malloc(forcings[l]->num_cells*sizeof(unsigned int*));
			}
			MPI_Bcast(forcings[l]->num_links_in_grid,forcings[l]->num_cells,MPI_UNSIGNED,0,MPI_COMM_WORLD);
			if(my_rank != 0)
			{
				for(i=0;i<forcings[l]->num_cells;i++)
					forcings[l]->grid_to_linkid[i] = (unsigned int*) malloc(forcings[l]->num_links_in_grid[i]*sizeof(unsigned int));
			}
			for(i=0;i<forcings[l]->num_cells;i++)
				MPI_Bcast(forcings[l]->grid_to_linkid[i],forcings[l]->num_links_in_grid[i],MPI_UNSIGNED,0,MPI_COMM_WORLD);

			forcings[l]->received = (char*) malloc(forcings[l]->num_cells*sizeof(char));
			forcings[l]->intensities = (float*) malloc(forcings[l]->num_cells*sizeof(float));

			//Remove from grid_to_linkid all links not on this proc
			unsigned int drop;
			for(i=0;i<forcings[l]->num_cells;i++)
			{
				drop = 0;
				for(j=0;j<forcings[l]->num_links_in_grid[i];j++)
				{
					if(assignments[forcings[l]->grid_to_linkid[i][j]] != my_rank)
						drop++;
					else
						forcings[l]->grid_to_linkid[i][j-drop] = forcings[l]->grid_to_linkid[i][j];
				}
				forcings[l]->num_links_in_grid[i] -= drop;
				forcings[l]->grid_to_linkid[i] = (unsigned int*) realloc(forcings[l]->grid_to_linkid[i],forcings[l]->num_links_in_grid[i]*sizeof(unsigned int));
			}
		}
		else if(forcings[l]->flag == 3 || forcings[l]->flag == 9)	//Database (could be irregular timesteps)
		{
			//Set routines
			if(forcings[l]->flag == 3)
			{
				forcings[l]->GetPasses = &PassesDatabase;
				forcings[l]->GetNextForcing = &NextForcingDatabase;
			}
			else
			{
				forcings[l]->GetPasses = &PassesDatabase_Irregular;
				forcings[l]->GetNextForcing = &NextForcingDatabase_Irregular;
			}

			//Get good_timestamp
			unsigned int good_time,is_null;
			char* query = db_connections[ASYNCH_DB_LOC_FORCING_START+l]->query;
			if(my_rank == 0)
			{
				ConnectPGDB(db_connections[ASYNCH_DB_LOC_FORCING_START+l]);
				sprintf(query,db_connections[ASYNCH_DB_LOC_FORCING_START+l]->queries[2],forcings[l]->first_file);
				res = PQexec(db_connections[ASYNCH_DB_LOC_FORCING_START+l]->conn,query);
				if(CheckResError(res,"querying a valid forcing time"))	return 1;
				is_null = PQgetisnull(res,0,0);
				if(is_null)	printf("Warning: forcing %u may have no data...\n",i);
				else		good_time = atoi(PQgetvalue(res,0,0));
				PQclear(res);
				DisconnectPGDB(db_connections[ASYNCH_DB_LOC_FORCING_START+l]);
			}

			MPI_Bcast(&is_null,1,MPI_INT,0,MPI_COMM_WORLD);
			if(is_null)
				forcings[l]->good_timestamp = forcings[l]->raindb_start_time;
			else
			{
				MPI_Bcast(&good_time,1,MPI_INT,0,MPI_COMM_WORLD);
				forcings[l]->good_timestamp = good_time;
			}

			if(forcings[l]->flag == 9)	//Download some extra data if the first_file is not on a timestamp with a forcing
				forcings[l]->next_timestamp = forcings[l]->good_timestamp;	//!!!! Could probably do something similar for flag 3 !!!!

			//Allocate space
			for(i=0;i<N;i++)
			{
				if(assignments[i] == my_rank)
				{
					if(!(GlobalVars->res_flag) || !(l == GlobalVars->res_forcing_idx) || system[i]->res)
					{
						m = forcings[l]->increment + 4;	//+1 for init, +1 for ceiling, +2 for when init time doesn't line up with file_time
						system[i]->forcing_buff[l] = (ForcingData*) malloc(sizeof(ForcingData));
						system[i]->forcing_buff[l]->rainfall = (double**) malloc(m*sizeof(double*));
						system[i]->forcing_buff[l]->n_times = m;
						for(j=0;j<m;j++)	system[i]->forcing_buff[l]->rainfall[j] = (double*) calloc(2,sizeof(double));
							system[i]->forcing_buff[l]->rainfall[0][0] = GlobalVars->t_0;
						system[i]->forcing_values[l] = 0.0;
						system[i]->forcing_change_times[l] = fabs(GlobalVars->t_0 + GlobalVars->maxtime) + 10.0;	//Just pick something away from t_0, and positive
					}
					else	//Reservoir, so allocate only a little memory
					{
						m = 4;	//+1 for init, +1 for ceiling, +2 for when init time doesn't line up with file_time
						system[i]->forcing_buff[l] = (ForcingData*) malloc(sizeof(ForcingData));
						system[i]->forcing_buff[l]->rainfall = (double**) malloc(m*sizeof(double*));
						system[i]->forcing_buff[l]->n_times = m;
						for(j=0;j<m;j++)	system[i]->forcing_buff[l]->rainfall[j] = (double*) calloc(2,sizeof(double));
							system[i]->forcing_buff[l]->rainfall[0][0] = GlobalVars->t_0;
						system[i]->forcing_values[l] = 0.0;
						system[i]->forcing_change_times[l] = fabs(GlobalVars->t_0 + GlobalVars->maxtime) + 10.0;	//Just pick something away from t_0, and positive
					}
				}
			}
		}
		else if(forcings[l]->flag == 4)	//.ustr
		{
			//Set routines
			forcings[l]->GetPasses = &PassesOther;
			forcings[l]->GetNextForcing = &NextForcingOther;

			//Read uniform data
			if(my_rank == 0)
			{
				forcingfile = fopen(forcings[l]->filename,"r");
				if(!forcingfile)
				{
					printf("[%i]: Error: cannot open uniform forcing file %s.\n",my_rank,forcings[l]->filename);
					return 1;
				}
				if(CheckWinFormat(forcingfile))
				{
					printf("Error: File %s appears to be in Windows format. Try converting to unix format using 'dos2unix' at the command line.\n",forcings[l]->filename);
					fclose(forcingfile);
					return 1;
				}

				//Read the number of times for the rain data for this link
				fscanf(forcingfile,"%i",&m);
				m++;		//Increase this by one to add a "ceiling" term
				buffer = (double*) realloc(buffer,2*m*sizeof(double));

				//Read the data
				for(j=0;j<m-1;j++)
					fscanf(forcingfile,"%lf %lf",&(buffer[2*j]),&(buffer[2*j+1]));
				buffer[2*j] = GlobalVars->maxtime + 3.0;
				buffer[2*j+1] = -1.0;

				fclose(forcingfile);
			}

			MPI_Bcast(&m,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
			if(my_rank != 0)
				buffer = (double*) realloc(buffer,2*m*sizeof(double));
			MPI_Bcast(buffer,2*m,MPI_DOUBLE,0,MPI_COMM_WORLD);

			//Create a global forcing object
			forcings[l]->GlobalForcing = (ForcingData*) malloc(sizeof(ForcingData));
			forcings[l]->GlobalForcing->rainfall = (double**) malloc(m*sizeof(double*));
			forcings[l]->GlobalForcing->n_times = m;
			for(j=0;j<m;j++)	forcings[l]->GlobalForcing->rainfall[j] = (double*) malloc(2*sizeof(double));

			for(j=0;j<m;j++)
			{
				forcings[l]->GlobalForcing->rainfall[j][0] = buffer[2*j];
				forcings[l]->GlobalForcing->rainfall[j][1] = buffer[2*j+1];
			}

			double rainfall_buffer = forcings[l]->GlobalForcing->rainfall[0][1];
			for(j=1;j<forcings[l]->GlobalForcing->n_times;j++)
			{
				if(rainfall_buffer != forcings[l]->GlobalForcing->rainfall[j][1])
				{
					univ_forcing_change_time[l] = forcings[l]->GlobalForcing->rainfall[j][0];
					break;
				}
			}
			if(j == forcings[l]->GlobalForcing->n_times)
				univ_forcing_change_time[l] = forcings[l]->GlobalForcing->rainfall[j-1][0];

			//Load the links
			for(i=0;i<N;i++)
			{
				if(assignments[i] == my_rank)
				{
					system[i]->forcing_buff[l] = forcings[l]->GlobalForcing;
					system[i]->forcing_values[l] = system[i]->forcing_buff[l]->rainfall[0][1];
					system[i]->forcing_change_times[l] = univ_forcing_change_time[l];
					system[i]->forcing_indices[l] = 0;
				}
			}
		}
		else if(forcings[l]->flag == 7)	//monthly recurring files
		{
			//Set routines
			forcings[l]->GetPasses = &PassesRecurring;
			forcings[l]->GetNextForcing = &NextForcingRecurring;

			int num_months = 12;
			buffer = (double*) realloc(buffer,num_months*sizeof(double));

			//Read monthly file
			if(my_rank == 0)
			{
				forcingfile = fopen(forcings[l]->filename,"r");
				if(!forcingfile)
				{
					printf("Error: cannot open uniform forcing file %s.\n",forcings[l]->filename);
					return 1;
				}
				if(CheckWinFormat(forcingfile))
				{
					printf("Error: File %s appears to be in Windows format. Try converting to unix format using 'dos2unix' at the command line.\n",forcings[l]->filename);
					fclose(forcingfile);
					return 1;
				}

				//Read the data
				for(j=0;j<num_months;j++)	fscanf(forcingfile,"%lf",&(buffer[j]));

				fclose(forcingfile);
			}

			MPI_Bcast(buffer,num_months,MPI_DOUBLE,0,MPI_COMM_WORLD);

			//Create a global forcing object
			forcings[l]->GlobalForcing = (ForcingData*) malloc(sizeof(ForcingData));
			forcings[l]->GlobalForcing->rainfall = (double**) malloc((num_months+1)*sizeof(double*));
			forcings[l]->GlobalForcing->n_times = num_months + 1;
			for(j=0;j<num_months+1;j++)	forcings[l]->GlobalForcing->rainfall[j] = (double*) malloc(2*sizeof(double));
			for(j=0;j<num_months;j++)	forcings[l]->GlobalForcing->rainfall[j][1] = buffer[j];
			forcings[l]->GlobalForcing->rainfall[num_months][1] = -1.0;

			//Find the starting month, and use the next month for the change time
			time_t start_time_t = forcings[l]->first_file;
			struct tm *start_time = gmtime(&start_time_t);

			(start_time->tm_mon)++;
			start_time->tm_mday = 1;
			start_time->tm_hour = 0;
			start_time->tm_min = 0;
			start_time->tm_sec = 0;
			time_t next_month = mktime(start_time);
			univ_forcing_change_time[l] = difftime(next_month,forcings[l]->first_file)/60.0;

			//Load the links
			for(i=0;i<N;i++)
			{
				if(assignments[i] == my_rank)
				{
					system[i]->forcing_buff[l] = forcings[l]->GlobalForcing;
					system[i]->forcing_values[l] = system[i]->forcing_buff[l]->rainfall[0][1];
					system[i]->forcing_change_times[l] = univ_forcing_change_time[l];
					system[i]->forcing_indices[l] = 0;
				}
			}
		}
		else
		{
			if(my_rank == 0)	printf("Error: Bad forcing flag for forcing %u.\n",forcings[l]->flag);
			return 1;
		}

	}

	//Find links with state forcing
	if(GlobalVars->res_flag)
	{
		//Download forcing data	!!!! Not sure if this is the way to go. Maybe separate function? !!!!
		forcings[GlobalVars->res_forcing_idx]->passes = 1;
		unsigned int start_iteration = forcings[GlobalVars->res_forcing_idx]->iteration;
		forcings[GlobalVars->res_forcing_idx]->GetNextForcing(system,N,my_sys,my_N,assignments,GlobalVars,forcings[GlobalVars->res_forcing_idx],db_connections,id_to_loc,GlobalVars->res_forcing_idx);
		forcings[GlobalVars->res_forcing_idx]->iteration = start_iteration;	//Keep this the same

		//Setup init condition at links with forcing
		Exchange_InitState_At_Forced(system,N,assignments,getting,res_list,res_size,id_to_loc,GlobalVars);
	}

	//Clean up
	if(buffer)	free(buffer);
	return 0;
}


//Reads any dam sources. Sets up appropriate buffers. Also sets flag for reservoirs.
int Load_Dams(Link** system,unsigned int N,unsigned int* my_sys,unsigned int my_N,int* assignments,short int* getting,unsigned int** id_to_loc,UnivVars* GlobalVars,ErrorData* GlobalErrors,ConnData** db_connections,unsigned int** res_list,unsigned int* res_size,unsigned int* my_res_size)
{
	unsigned int i,j,k,m,num_dams,id,size,num_values;
	//int error;
	Link* current;
	FILE* damfile = NULL;
	double* buffer = NULL;
	PGresult *res;

	//Read dam file
	if(GlobalVars->uses_dam && GlobalVars->dam_flag == 1)	//.dam file
	{
		//VEC* buffer;
		size = GlobalVars->dam_params_size - GlobalVars->params_size;
		buffer = (double*) malloc(size*sizeof(double));

		//Setup needs array. This is the procs that have getting set to 1.
		int* needs = NULL;
		if(my_rank == 0) needs = (int*) malloc(N*sizeof(int));
		int not_needed = -1;
		for(i=0;i<N;i++)
		{
			if(getting[i])	MPI_Reduce(&my_rank,&(needs[i]),1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
			else		MPI_Reduce(&not_needed,&(needs[i]),1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
		}

		if(my_rank == 0)
		{
			damfile = fopen(GlobalVars->dam_filename,"r");
			if(damfile == NULL)
			{
				printf("Error: Dam file %s not found.\n",GlobalVars->dam_filename);
				return 1;
			}
			if(CheckWinFormat(damfile))
			{
				printf("Error: File %s appears to be in Windows format. Try converting to unix format using 'dos2unix' at the command line.\n",GlobalVars->dam_filename);
				fclose(damfile);
				return 1;
			}

			fscanf(damfile,"%u",&num_dams);
			MPI_Bcast(&num_dams,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);

			for(i=0;i<num_dams;i++)
			{
				fscanf(damfile,"%u",&id);
				m = find_link_by_idtoloc(id,id_to_loc,N);
				if(m > N)
				{
					printf("Error: Link id %u from .dam file %s not present in network.\n",id,GlobalVars->dam_filename);
					return 1;
				}

				for(j=0;j<size;j++)
					fscanf(damfile,"%lf",&(buffer[j]));

				MPI_Bcast(&m,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
				
				//Either store the parameters or send them
				if(my_rank == assignments[m] || my_rank == needs[m])
				{
					system[m]->params.ve = (double*) realloc(system[m]->params.ve,GlobalVars->dam_params_size*sizeof(double));
					system[m]->params.dim = GlobalVars->dam_params_size;

					for(j=GlobalVars->params_size;j<GlobalVars->dam_params_size;j++)
						system[m]->params.ve[j] = buffer[j-GlobalVars->params_size];
					system[m]->dam = 1;
				}

				if(my_rank != assignments[m])
					MPI_Send(buffer,size,MPI_DOUBLE,assignments[m],1,MPI_COMM_WORLD);

				if(needs[m] > 0)
					MPI_Send(buffer,size,MPI_DOUBLE,needs[m],1,MPI_COMM_WORLD);
			}

			fclose(damfile);
			free(needs);
		}
		else
		{
			MPI_Bcast(&num_dams,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);

			for(i=0;i<num_dams;i++)
			{
				MPI_Bcast(&m,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
				if(my_rank == assignments[m] || getting[m])
				{
					MPI_Recv(buffer,size,MPI_DOUBLE,0,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

					system[m]->params.ve = (double*) realloc(system[m]->params.ve,GlobalVars->dam_params_size*sizeof(double));
					system[m]->params.dim = GlobalVars->dam_params_size;

					for(j=GlobalVars->params_size;j<GlobalVars->dam_params_size;j++)
						system[m]->params.ve[j] = buffer[j-GlobalVars->params_size];
					system[m]->dam = 1;
				}
			}
		}

		//Clean up
		if(buffer)	free(buffer);

		//Set error tolerance
		if(GlobalVars->rkd_flag == 0)
		{
			GlobalVars->discont_tol = GlobalErrors->abstol.ve[0];
			if(GlobalVars->discont_tol < 1e-12 && my_rank == 0)
				printf("Warning: Discontinuity tolerance has been set to %e.\n",GlobalVars->discont_tol);
		}
		else
		{
			GlobalVars->discont_tol = 1e-8;
			if(my_rank == 0)	printf("Notice: Discontinuity tolerance has been set to %e.\n",GlobalVars->discont_tol);
		}
	}
	else if(GlobalVars->uses_dam && GlobalVars->dam_flag == 2)	//.qvs file
	{
		//Setup needs array. This is the procs that have getting set to 1.
		int* needs = NULL;
		if(my_rank == 0) needs = (int*) malloc(N*sizeof(int));
		int not_needed = -1;
		for(i=0;i<N;i++)
		{
			if(getting[i])	MPI_Reduce(&my_rank,&(needs[i]),1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
			else		MPI_Reduce(&not_needed,&(needs[i]),1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
		}

		if(my_rank == 0)
		{
			damfile = fopen(GlobalVars->dam_filename,"r");
			if(damfile == NULL)
			{
				printf("Error: Dam file %s not found.\n",GlobalVars->dam_filename);
				return 1;
			}
			if(CheckWinFormat(damfile))
			{
				printf("Error: File %s appears to be in Windows format. Try converting to unix format using 'dos2unix' at the command line.\n",GlobalVars->dam_filename);
				fclose(damfile);
				return 1;
			}
		
			fscanf(damfile,"%u",&num_dams);
			MPI_Bcast(&num_dams,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);

			for(i=0;i<num_dams;i++)
			{
				//Get link id and number of values
				fscanf(damfile,"%u %u",&id,&num_values);
				m = find_link_by_idtoloc(id,id_to_loc,N);
				if(m > N)
				{
					printf("Error: Link id %u from .qvs file %s not present in network.\n",id,GlobalVars->dam_filename);
					return 1;
				}

				//Read q vs s data
				buffer = realloc(buffer,2*num_values*sizeof(double));
				for(j=0;j<num_values;j++)
					fscanf(damfile,"%lf %lf",&(buffer[2*j]),&(buffer[2*j+1]));

				//Send location
				MPI_Bcast(&m,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
				
				//Either store the parameters or send them
				if(my_rank == assignments[m] || getting[m])
				{
					system[m]->dam = 1;
					system[m]->qvs = (QVSData*) malloc(sizeof(QVSData));
					system[m]->qvs->n_values = num_values;
					system[m]->qvs->points_array = (double*) malloc(2*num_values*sizeof(double));
					system[m]->qvs->points = (double**) malloc(num_values*sizeof(double*));
					for(j=0;j<num_values;j++)	system[m]->qvs->points[j] = &(system[m]->qvs->points_array[2*j]);

					for(j=0;j<num_values;j++)
					{
						system[m]->qvs->points[j][0] = buffer[2*j];
						system[m]->qvs->points[j][1] = buffer[2*j+1];
					}
				}

				if(my_rank != assignments[m])
				{
					MPI_Send(&num_values,1,MPI_UNSIGNED,assignments[m],1,MPI_COMM_WORLD);
					MPI_Send(buffer,2*num_values,MPI_DOUBLE,assignments[m],1,MPI_COMM_WORLD);
				}

				if(needs[m] > 0)
				{
					MPI_Send(&num_values,1,MPI_UNSIGNED,needs[m],1,MPI_COMM_WORLD);
					MPI_Send(buffer,2*num_values,MPI_DOUBLE,needs[m],1,MPI_COMM_WORLD);
				}
			}

			fclose(damfile);
			free(needs);
		}
		else
		{
			MPI_Bcast(&num_dams,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);

			for(i=0;i<num_dams;i++)
			{
				MPI_Bcast(&m,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
				if(my_rank == assignments[m] || getting[m])
				{
					MPI_Recv(&num_values,1,MPI_UNSIGNED,0,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					buffer = (double*) realloc(buffer,2*num_values*sizeof(double));
					MPI_Recv(buffer,2*num_values,MPI_DOUBLE,0,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

					system[m]->dam = 1;
					system[m]->qvs = (QVSData*) malloc(sizeof(QVSData));
					system[m]->qvs->n_values = num_values;
					system[m]->qvs->points_array = (double*) malloc(2*num_values*sizeof(double));
					system[m]->qvs->points = (double**) malloc(num_values*sizeof(double*));
					for(j=0;j<num_values;j++)	system[m]->qvs->points[j] = &(system[m]->qvs->points_array[2*j]);

					for(j=0;j<num_values;j++)
					{
						system[m]->qvs->points[j][0] = buffer[2*j];
						system[m]->qvs->points[j][1] = buffer[2*j+1];
					}
				}
			}
		}

		//Clean up
		if(buffer)	free(buffer);

		//Set error tolerance
		if(GlobalVars->rkd_flag == 0)
		{
			GlobalVars->discont_tol = GlobalErrors->abstol.ve[0];
			if(GlobalVars->discont_tol < 1e-12 && my_rank == 0)
				printf("Warning: Discontinuity tolerance has been set to %e.\n",GlobalVars->discont_tol);
		}
		else
		{
			GlobalVars->discont_tol = 1e-8;
			if(my_rank == 0)	printf("Notice: Discontinuity tolerance has been set to %e.\n",GlobalVars->discont_tol);
		}
	}
	else if(GlobalVars->uses_dam && GlobalVars->dam_flag == 3)	//database connection
	{
		Link* current;
		unsigned int num_pts,curr_loc;
		num_dams = 0;
		short int procs_sending_to[ASYNCH_MAX_NUMBER_OF_PROCESS],mine;
		double* array_holder;

		//Connect to the database and download data
		if(my_rank == 0)
		{
			ConnectPGDB(db_connections[ASYNCH_DB_LOC_QVS]);
			res = PQexec(db_connections[ASYNCH_DB_LOC_QVS]->conn,db_connections[ASYNCH_DB_LOC_QVS]->queries[0]);
			if(CheckResError(res,"querying qvs relations"))	return 1;
			num_pts = PQntuples(res);

			i = 0;
			while(i < num_pts)
			{
				//Get link id and location
				id = atoi(PQgetvalue(res,i,0));
				curr_loc = find_link_by_idtoloc(id,id_to_loc,N);

				//Count number of qvs values for current
				for(j=i;j<num_pts && id == atoi(PQgetvalue(res,j,0));j++);
				num_values = j-i;

				//Extract the data points
				array_holder = (double*) malloc(2*num_values*sizeof(double));
				for(j=0;j<num_values;j++)
				{
					array_holder[2*j] = atof(PQgetvalue(res,i+j,1));
					array_holder[2*j+1] = atof(PQgetvalue(res,i+j,2));
				}

				//Check for an error real quick
				for(j=1;j<num_values;j++)
				{
					if(array_holder[2*(j-1)] > array_holder[2*j] || array_holder[2*(j-1)+1] > array_holder[2*j+1])
					{
						printf("[%i]: Bad storage or discharge values found at link id %u. Check that the data is sorted correctly. (%u)\n",my_rank,id,j);
						//break;
						free(array_holder);
						return 1;
					}
				}

				if(curr_loc < N)
				{
					//Tell everyone what link has the current dam
					MPI_Bcast(&curr_loc,1,MPI_INT,0,MPI_COMM_WORLD);
					mine = (my_rank == assignments[curr_loc] || getting[curr_loc]);
					MPI_Gather(&mine,1,MPI_SHORT,procs_sending_to,1,MPI_SHORT,0,MPI_COMM_WORLD);

					//Send the data to whoever needs it
					for(j=1;j<np;j++)
					{
						if(procs_sending_to[j])
						{
							MPI_Send(&num_values,1,MPI_UNSIGNED,j,0,MPI_COMM_WORLD);
							MPI_Send(array_holder,2*num_values,MPI_DOUBLE,j,0,MPI_COMM_WORLD);
						}
					}

					//Check if proc 0 needs the data
					if(mine)
					{
						current = system[curr_loc];
						current->dam = 1;
						current->qvs = (QVSData*) malloc(sizeof(QVSData));
						current->qvs->n_values = num_values;
						current->qvs->points_array = array_holder;
						current->qvs->points = (double**) malloc(num_values*sizeof(double*));
						for(j=0;j<num_values;j++)	current->qvs->points[j] = &(current->qvs->points_array[2*j]);
					}
					else
						free(array_holder);
				}
				else
					free(array_holder);

				i += num_values;
			}

			//Finish up
			PQclear(res);
			DisconnectPGDB(db_connections[ASYNCH_DB_LOC_QVS]);
			curr_loc = -1;
			MPI_Bcast(&curr_loc,1,MPI_INT,0,MPI_COMM_WORLD);
		}
		else	//Receive dam data
		{
			curr_loc = 0;
			MPI_Bcast(&curr_loc,1,MPI_INT,0,MPI_COMM_WORLD);

			while((int)curr_loc != -1)
			{
				//Check if I need the data for this link
				mine = (my_rank == assignments[curr_loc] || getting[curr_loc]);
				MPI_Gather(&mine,1,MPI_SHORT,procs_sending_to,1,MPI_SHORT,0,MPI_COMM_WORLD);

				//Receive data
				if(mine)
				{
					MPI_Recv(&num_values,1,MPI_UNSIGNED,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					array_holder = (double*) malloc(2*num_values*sizeof(double));
					MPI_Recv(array_holder,2*num_values,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					current = system[curr_loc];
					current->dam = 1;
					current->qvs = (QVSData*) malloc(sizeof(QVSData));
					current->qvs->n_values = num_values;
					current->qvs->points_array = array_holder;
					current->qvs->points = (double**) malloc(num_values*sizeof(double*));
					for(j=0;j<num_values;j++)	current->qvs->points[j] = &(current->qvs->points_array[2*j]);
				}

				//Check next signal
				MPI_Bcast(&curr_loc,1,MPI_INT,0,MPI_COMM_WORLD);
			}
		}

		//Set error tolerance
		if(GlobalVars->rkd_flag == 0)
		{
			GlobalVars->discont_tol = GlobalErrors->abstol.ve[0];
			if(GlobalVars->discont_tol < 1e-12 && my_rank == 0)
				printf("Warning: Discontinuity tolerance has been set to %e.\n",GlobalVars->discont_tol);
		}
		else
		{
			GlobalVars->discont_tol = 1e-8;
			if(my_rank == 0)	printf("Notice: Discontinuity tolerance has been set to %e.\n",GlobalVars->discont_tol);
		}
	}
	else	//Some other type of discontinuity (or none at all)
	{
		//Set error tolerance
		if(GlobalVars->rkd_flag == 0)
		{
			GlobalVars->discont_tol = GlobalErrors->abstol.ve[0];
			if(GlobalVars->discont_tol < 1e-12 && my_rank == 0)
				printf("Warning: Discontinuity tolerance has been set to %e.\n",GlobalVars->discont_tol);
		}
		else
		{
			GlobalVars->discont_tol = 1e-8;
			//if(my_rank == 0)	printf("Notice: Discontinuity tolerance has been set to %e.\n",GlobalVars->discont_tol);
		}
	}

	//Check where state forcings should be set
	*my_res_size = 0;
	if(GlobalVars->res_flag)
	{
		if(Create_SAV_Data(GlobalVars->rsv_filename,system,N,res_list,res_size,db_connections[ASYNCH_DB_LOC_RSV],GlobalVars->res_flag))
			return 1;

		//Setup links with forcing
		for(j=0;j<*res_size;j++)
		{
			k = find_link_by_idtoloc((*res_list)[j],id_to_loc,N);

			if(k == N+1)
			{
				if(my_rank == 0)
					printf("Warning: Reservoir forcing requested for a link NOT in the network (link id = %u).\n",(*res_list)[j]);

				//Shift ids in peak save list to spot k, and decrement peaksave_size
				for(k=j+1;k<*res_size;k++)
					(*res_list)[k-1] = (*res_list)[k];
				j--;
				(*res_size)--;
			}
			else if(assignments[k] == my_rank || getting[k])
			{
				(*my_res_size)++;
				current = system[k];
				current->res = 1;
			}
		}
	}
	else
	{
		(*res_list) = NULL;
		*res_size = 0;
	}

	return 0;
}

//Calculates an estimate for the initial stepsizes at each link.
//Returns 0 if initial conditions set. 1 if a forcing is not ready.
int CalculateInitialStepSizes(Link** system,unsigned int* my_sys,unsigned int my_N,UnivVars* GlobalVars,TempStorage* workspace,short int watch_forcings)
{
	unsigned int i,j;

	for(i=0;i<my_N;i++)
		system[my_sys[i]]->h = InitialStepSize(GlobalVars->t_0,system[my_sys[i]],GlobalVars,workspace);
	if(watch_forcings)
	{
		for(i=0;i<my_N;i++)
		{
			for(j=0;j<GlobalVars->num_forcings;j++)
				if(system[my_sys[i]]->forcing_buff[j])
					system[my_sys[i]]->h = min(system[my_sys[i]]->h,system[my_sys[i]]->forcing_change_times[j] - system[my_sys[i]]->last_t);
		}
	}

	return 0;
}

//Reads the link ids for saving data.
//Sets in GlobalVars: save_list, save_size, peaksave_list, peaksave_size.
//Also sets: my_save_size.
int BuildSaveLists(Link** system,unsigned int N,unsigned int* my_sys,unsigned int my_N,unsigned int* assignments,unsigned int** id_to_loc,UnivVars* GlobalVars,unsigned int** save_list,unsigned int* save_size,unsigned int* my_save_size,unsigned int** peaksave_list,unsigned int* peaksave_size,unsigned int* my_peaksave_size,ConnData** db_connections)
{
	unsigned int j,k;
	Link* current;

	if(!save_list || !peaksave_list)
	{
		if(my_rank == 0)	printf("Error: id list for output data cannot be null.\n");
		return 1;
	}

	//For time series
	*my_save_size = 0;
	if(GlobalVars->hydros_loc_flag)
	{
		if(Create_SAV_Data(GlobalVars->hydrosave_filename,system,N,save_list,save_size,db_connections[ASYNCH_DB_LOC_HYDROSAVE],GlobalVars->hydrosave_flag))
			return 1;

		for(j=0;j<*save_size;j++)
		{
			k = find_link_by_idtoloc((*save_list)[j],id_to_loc,N);

			if(k == N+1)
			{
				if(my_rank == 0)
					printf("Warning: Time series output requested for a link NOT in the network (link id = %u).\n",(*save_list)[j]);

				//Shift ids in save list to spot k, and decrement save_size
				for(k=j+1;k<*save_size;k++)
					(*save_list)[k-1] = (*save_list)[k];
				j--;
				(*save_size)--;
			}
			else if(assignments[k] == my_rank)
			{
				(*my_save_size)++;
				current = system[k];
				current->save_flag = 1;
				current->next_save = GlobalVars->t_0;
				if(GlobalVars->print_time > 0.0)
					current->print_time = GlobalVars->print_time;
				else
					current->print_time = pow(current->params.ve[GlobalVars->area_idx]*0.1,0.5);
			}
		}
	}
	else
	{
		(*save_list) = NULL;
		*save_size = 0;
	}

	//For peakflows
	*my_peaksave_size = 0;
	if(GlobalVars->peaks_loc_flag)
	{
		if(Create_SAV_Data(GlobalVars->peaksave_filename,system,N,peaksave_list,peaksave_size,db_connections[ASYNCH_DB_LOC_PEAKSAVE],GlobalVars->peaksave_flag))
			return 1;

		for(j=0;j<*peaksave_size;j++)
		{
			k = find_link_by_idtoloc((*peaksave_list)[j],id_to_loc,N);

			if(k == N+1)
			{
				if(my_rank == 0)
					printf("Warning: Peakflow output requested for a link NOT in the network (link id = %u).\n",(*peaksave_list)[j]);

				//Shift ids in peak save list to spot k, and decrement peaksave_size
				for(k=j+1;k<*peaksave_size;k++)
					(*peaksave_list)[k-1] = (*peaksave_list)[k];
				j--;
				(*peaksave_size)--;
			}
			else if(assignments[k] == my_rank)
			{
				(*my_peaksave_size)++;
				current = system[k];
				current->peak_flag = 1;
			}
		}
	}
	else
	{
		(*peaksave_list) = NULL;
		*peaksave_size = 0;
	}

	return 0;
}

//Finalizes the system for calculations. Basically, lots of small miscellaneous initializations are made.
int FinalizeSystem(Link** system,unsigned int N,unsigned int* my_sys,unsigned int my_N,unsigned int* assignments,short int* getting,unsigned int** id_to_loc,TransData* my_data,UnivVars* GlobalVars,ConnData** db_connections,TempStorage** workspace)
{
	unsigned int i,j;
	int ii;
	Link* current;

	//Initialize workspace
	*workspace = Create_Workspace(GlobalVars->max_dim,GlobalVars->max_s,GlobalVars->max_parents);

	//Need space for nodes, number of iterations, discontinuities
	//Data: ( size(double)*(max_s*max_dim + max_dim + time)*# steps to transfer + size(int)(for stage)*# steps to transfer + size(int)*(location + # steps to transfer) ) * # of sending links
	//Upstream: + size(int) * (location + # of iterations) * # of receiving links
	//Discontinuities: + (size(int) + size(double)*discont_size + size(int)*discont_size) * # of sending links
	//unsigned int bytes1 = ( (sizeof(double)*(GlobalVars->max_s*2 + 2 + 1) + sizeof(int) )*GlobalVars->max_transfer_steps + sizeof(int)*2);
	unsigned int bytes2 = 2*sizeof(int);
	unsigned int bytes3 = sizeof(int) + (sizeof(int) + sizeof(double))*GlobalVars->discont_size;
	for(ii=0;ii<np;ii++)
	{
		my_data->send_buffer_size[ii] = bytes2 * my_data->receive_size[ii] + bytes3 * my_data->send_size[ii];
		for(j=0;j<my_data->send_size[ii];j++)
		{
			current = my_data->send_data[ii][j];
			my_data->send_buffer_size[ii] += (sizeof(double)*(current->method->s*current->num_dense + current->dim + 1) + sizeof(int) )*GlobalVars->max_transfer_steps + sizeof(int)*2;
		}

		my_data->receive_buffer_size[ii] = bytes2 * my_data->send_size[ii] + bytes3 * my_data->receive_size[ii];
		for(j=0;j<my_data->receive_size[ii];j++)
		{
			current = my_data->receive_data[ii][j];
			my_data->receive_buffer_size[ii] += ( sizeof(double)*(current->method->s*current->num_dense + current->dim + 1) + sizeof(int) )*GlobalVars->max_transfer_steps + sizeof(int)*2;
		}

		//(*my_data)->send_buffer_size[ii] = bytes1 * (*my_data)->send_size[ii] + bytes2 * (*my_data)->receive_size[ii] + bytes3 * (*my_data)->send_size[ii];
		//(*my_data)->receive_buffer_size[ii] = bytes1 * (*my_data)->receive_size[ii] + bytes2*(*my_data)->send_size[ii] + bytes3 * (*my_data)->receive_size[ii];

		if(my_data->send_buffer_size[ii])	my_data->send_buffer[ii] = (char*) malloc(my_data->send_buffer_size[ii]*sizeof(char));
		else					my_data->send_buffer[ii] = NULL;

		if(my_data->receive_buffer_size[ii])	my_data->receive_buffer[ii] = (char*) malloc(my_data->receive_buffer_size[ii]*sizeof(char));
		else					my_data->receive_buffer[ii] = NULL;
	}

	//Do initializations
	for(i=0;i<N;i++)
	{
		if(assignments[i] == my_rank || getting[i])
		{
			//Discontinuity information
			if(system[i]->numparents)
				system[i]->discont = (double*) malloc(GlobalVars->discont_size*sizeof(double));
			if(system[i]->c && my_rank != assignments[system[i]->c->location])
			{
				system[i]->discont_send = (double*) malloc(GlobalVars->discont_size*sizeof(double));
				system[i]->discont_order_send = (unsigned int*) malloc(GlobalVars->discont_size*sizeof(unsigned int));
				system[i]->discont_send_count = 0;
			}

			//Setup most of the remaining data
			//system[i]->last_t = GlobalVars->t_0;
			//system[i]->print_time = GlobalVars->t_0;
			system[i]->current_iterations = 1;
			if(!(system[i]->save_flag))	system[i]->next_save = GlobalVars->t_0 - 1.0;
			system[i]->peak_time = GlobalVars->t_0;
			//system[i]->save_flag = 0;
			//system[i]->peak_flag = 0;
			system[i]->peak_value = v_get(system[i]->dim);
			v_copy(system[i]->list->head->y_approx,system[i]->peak_value);
			if(system[i]->numparents)	system[i]->ready = 0;
			else				system[i]->ready = 1;
			system[i]->iters_removed = 0;
			system[i]->steps_on_diff_proc = 1;	//Note: This won't be used if link isn't stored on another proc
			system[i]->rejected = 0;
			system[i]->last_eta = 1e10;
			system[i]->compute_J = 1;
			system[i]->compute_LU = 1;

			if(system[i]->method->exp_imp == 1)
			{
				system[i]->JMatrix = m_get(system[i]->dim,system[i]->dim);
				system[i]->CoefMat = m_get(GlobalVars->max_s*system[i]->dim,GlobalVars->max_s*system[i]->dim);
				system[i]->Z_i = malloc(GlobalVars->max_s*sizeof(VEC*));
				for(j=0;j<GlobalVars->max_s;j++)
					system[i]->Z_i[j] = v_get(system[i]->dim);
				system[i]->sol_diff = v_get(system[i]->dim);
				system[i]->h_old = -1.0;
				system[i]->value_old = -1.0;
			}
			else
			{
				system[i]->JMatrix = m_get(0, 0);
				system[i]->CoefMat = m_get(0, 0);
				system[i]->Z_i = NULL;
				system[i]->sol_diff = v_get(0);
			}

		}
		else	//If link not assigned to this process, and not receiving any information from this link.
		{
			system[i]->method = NULL;
			system[i]->list = NULL;
			system[i]->peak_value = v_get(0);
			system[i]->params = v_get(0);

			//NULL out the unneeded data
			system[i]->f = NULL;
			system[i]->Jacobian = NULL;
			system[i]->JMatrix = m_get(0, 0);
			system[i]->CoefMat = m_get(0, 0);
			system[i]->Z_i = NULL;
			system[i]->sol_diff = v_get(0);
		}
	}

	return 0;
}


// *********************************************************************************************



//Reads in the data from a .gbl file. All data will be global data for the entire river system.
//char globalfilename[]: String with the filename of the .gbl file.
//ErrorData** GlobalErrors (set by this method): Will contain the error data for the entire river system, if the error data is global.
//PGconn** conn:	NULL pointer that will be set to an SQL database, if needed.
//char* rkdfilename (set by this method): Will be the filename of the .rkd file, if the error data is not global.
//Returns a UnivVars that contains all the global data read in from the file globalfilename.
UnivVars* Read_Global_Data(char globalfilename[],ErrorData** GlobalErrors,Forcing** forcings,ConnData** db_connections,char* rkdfilename,model* custom_model,void* external)
{
	unsigned int i,total,written;
	int flag,valsread;
	char endmark;
	unsigned int string_size = 256;	//This is the max size of all strings for filenames and location
	unsigned int buff_size = string_size + 20;
	char* linebuffer = (char*) malloc(buff_size*sizeof(char));
	UnivVars* GlobalVars = (UnivVars*) malloc(sizeof(UnivVars));
    memset(GlobalVars, 0, sizeof(UnivVars));

	GlobalVars->string_size = string_size;
	char* db_filename = (char*) malloc(string_size*sizeof(char));
	GlobalVars->query_size = 1024;

	*GlobalErrors = malloc(sizeof(ErrorData));
	FILE* globalfile = NULL;

	if(my_rank == 0)
	{
		globalfile = fopen(globalfilename,"r");
		if(globalfile == NULL)
		{
			printf("Error: Global file %s was not found.\n",globalfilename);
			return NULL;
		}

		if(CheckWinFormat(globalfile))
		{
			printf("Error: File %s appears to be in Windows format. Try converting to unix format using 'dos2unix' at the command line.\n",globalfilename);
			return NULL;
		}
	}

	GlobalVars->rain_filename = NULL;

	//Grab the type and maxtime
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%hu %lf",&(GlobalVars->type),&(GlobalVars->maxtime));
	if(ReadLineError(valsread,2,"type and maxtime"))	return NULL;

	//Grab the output filename info
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%hu",&(GlobalVars->print_par_flag));
	if(ReadLineError(valsread,1,"to print filename parameters"))	return NULL;

	//Grab components to print
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%u",&(GlobalVars->num_print));
	if(ReadLineError(valsread,1,"number of indices to print"))	return NULL;
	GlobalVars->output_names = (char**) malloc(GlobalVars->num_print*sizeof(char*));
	for(i=0;i<GlobalVars->num_print;i++)
	{
		GlobalVars->output_names[i] = (char*) malloc(string_size*sizeof(char));
		ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
		valsread = sscanf(linebuffer,"%s",GlobalVars->output_names[i]);
		if(ReadLineError(valsread,1,"a component to print"))	return NULL;
	}

	//Peakflow function
	GlobalVars->peakflow_function_name = (char*) malloc(string_size*sizeof(char));
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%s",GlobalVars->peakflow_function_name);
	if(ReadLineError(valsread,1,"peakflow function name"))	return NULL;
	SetPeakflowOutputFunctions(GlobalVars->peakflow_function_name,&(GlobalVars->peakflow_output));

	//Grab the parameters
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%u%n",&i,&total);
	if(ReadLineError(valsread,1,"number of global parameters"))	return NULL;
	GlobalVars->global_params = v_get(i);
	for(i=0;i<GlobalVars->global_params.dim;i++)
	{
		valsread = sscanf(&(linebuffer[total]),"%lf%n",&(GlobalVars->global_params.ve[i]),&written);
		if(ReadLineError(valsread,1,"a global parameter"))	return NULL;
		total += written;
	}

	//Set dim and other sizes
	if(custom_model)	custom_model->SetParamSizes(GlobalVars,external);
	else			SetParamSizes(GlobalVars,external);

	//Find the states needed for printing
	GlobalVars->num_states_for_printing = 0;
	GlobalVars->print_indices = (unsigned int*) calloc(GlobalVars->num_print,sizeof(unsigned int));
	GlobalVars->outputs_i = (int (**)(double,VEC,VEC,VEC,int,void*)) calloc( GlobalVars->num_print, sizeof(int (*)(double,VEC,VEC,VEC,int,void*)) );
	GlobalVars->outputs_d = (double (**)(double,VEC,VEC,VEC,int,void*)) calloc( GlobalVars->num_print, sizeof(double (*)(double,VEC,VEC,VEC,int,void*)) );
	GlobalVars->output_types = (short int*) malloc(GlobalVars->num_print*sizeof(short int));
	GlobalVars->output_sizes = (short int*) malloc(GlobalVars->num_print*sizeof(short int));
	GlobalVars->output_specifiers = (char**) malloc(GlobalVars->num_print*sizeof(char*));
	for(i=0;i<GlobalVars->num_print;i++)
	{
		GlobalVars->output_specifiers[i] = (char*) malloc(16*sizeof(char));
		SetOutputFunctions(GlobalVars->output_names[i],GlobalVars->output_specifiers[i],GlobalVars->print_indices,&(GlobalVars->num_states_for_printing),&(GlobalVars->output_sizes[i]),&(GlobalVars->output_types[i]),&(GlobalVars->outputs_i[i]),&(GlobalVars->outputs_d[i]));
	}
	GlobalVars->print_indices = (unsigned int*) realloc(GlobalVars->print_indices,GlobalVars->num_states_for_printing*sizeof(unsigned int));

	//Grab the stored steps limits
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%u %i %u",&(GlobalVars->iter_limit),&(GlobalVars->max_transfer_steps),&(GlobalVars->discont_size));
	if(ReadLineError(valsread,3,"steps stored, steps transfered, and discontinuity buffer size"))	return NULL;

	//Grab the topology data filename
	GlobalVars->outletlink = 0;
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%hu",&(GlobalVars->rvr_flag));
	if(ReadLineError(valsread,1,"topology data flag"))	return NULL;
	if(GlobalVars->rvr_flag == 0)
	{
		GlobalVars->rvr_filename = (char*) malloc(string_size*sizeof(char));
		valsread = sscanf(linebuffer,"%*u %s",GlobalVars->rvr_filename);
		if(ReadLineError(valsread,1,"filename for topology data"))	return NULL;
		if(!CheckFilenameExtension(GlobalVars->rvr_filename,".rvr"))	return NULL;
	}
	else	//Reading from database
	{
		valsread = sscanf(linebuffer,"%*u %u %s",&(GlobalVars->outletlink),db_filename);
		if(ReadLineError(valsread,2,"link id of downstream link for topology data or .dbc for topology"))	return NULL;
		if(!CheckFilenameExtension(db_filename,".dbc"))	return NULL;
		GlobalVars->rvr_filename = NULL;
		db_connections[ASYNCH_DB_LOC_TOPO] = ReadDBC(db_filename,string_size);
		if(db_connections[ASYNCH_DB_LOC_TOPO] == NULL)	return NULL;
	}

	//Grab the parameter data filename
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%hu",&(GlobalVars->prm_flag));
	if(ReadLineError(valsread,1,"parameter flag"))	return NULL;
	if(GlobalVars->prm_flag == 0)
	{
		GlobalVars->prm_filename = (char*) malloc(string_size*sizeof(char));
		valsread = sscanf(linebuffer,"%*u %s",GlobalVars->prm_filename);
		if(ReadLineError(valsread,1,".prm filename"))	return NULL;
		if(!CheckFilenameExtension(GlobalVars->prm_filename,".prm"))	return NULL;
	}
	else
	{
		valsread = sscanf(linebuffer,"%*u %s",db_filename);
		if(ReadLineError(valsread,1,".dbc for parameters"))	return NULL;
		if(!CheckFilenameExtension(db_filename,".dbc"))	return NULL;
		GlobalVars->prm_filename = NULL;
		db_connections[ASYNCH_DB_LOC_PARAMS] = ReadDBC(db_filename,string_size);
		if(!db_connections[ASYNCH_DB_LOC_PARAMS])	return NULL;
	}

	//Grab the initial data file
	GlobalVars->init_filename = NULL;
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%hu",&(GlobalVars->init_flag));
	if(ReadLineError(valsread,1,"initial data flag"))	return NULL;
	if(GlobalVars->init_flag == 0 || GlobalVars->init_flag == 1 || GlobalVars->init_flag == 2 || GlobalVars->init_flag == 4)
	{
		GlobalVars->init_filename = (char*) malloc(string_size*sizeof(char));
		valsread = sscanf(linebuffer,"%*u %s",GlobalVars->init_filename);
		if (ReadLineError(valsread,1,"initial data flag"))	return NULL;
		if (GlobalVars->init_flag == 0 && !CheckFilenameExtension(GlobalVars->init_filename,".ini")) return NULL;
		if (GlobalVars->init_flag == 1 && !CheckFilenameExtension(GlobalVars->init_filename,".uini")) return NULL;
		if (GlobalVars->init_flag == 2 && !CheckFilenameExtension(GlobalVars->init_filename,".rec")) return NULL;
        if (GlobalVars->init_flag == 4 && !CheckFilenameExtension(GlobalVars->init_filename, ".h5")) return NULL;
	}
	else if (GlobalVars->init_flag == 3)
	{
		valsread = sscanf(linebuffer,"%*u %s %u",db_filename,&(GlobalVars->init_timestamp));
		if(ReadLineError(valsread,1,".dbc for parameters"))	return NULL;
		if(!CheckFilenameExtension(db_filename,".dbc"))	return NULL;
		db_connections[ASYNCH_DB_LOC_INIT] = ReadDBC(db_filename,string_size);
		if(!db_connections[ASYNCH_DB_LOC_INIT])	return NULL;
	}

	//Grab number of forcings
	unsigned int got_forcings;
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%u",&got_forcings);
	if(ReadLineError(valsread,1,"rainfall flag"))	return NULL;
	if(got_forcings < GlobalVars->num_forcings && my_rank == 0)
	{
		printf("[%i]: Error: Got %u forcings in the .gbl file. Expected %u for model %u.\n",my_rank,got_forcings,GlobalVars->num_forcings,GlobalVars->type);
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	if(got_forcings > GlobalVars->num_forcings && my_rank == 0)
	{
		printf("[%i]: Warning: Got %u forcings in the .gbl file. Expected %u for model %u.\n",my_rank,got_forcings,GlobalVars->num_forcings,GlobalVars->type);
		GlobalVars->num_forcings = got_forcings;
	}

	//Grab the forcing parameters
	//0 for no rain, 1 for .str file, 2 for binary files, 3 for database, 4 for uniform rain (.ustr)
	GlobalVars->hydro_table = GlobalVars->peak_table = NULL;
	for(i=0;i<GlobalVars->num_forcings;i++)
	{
		ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
		valsread = sscanf(linebuffer,"%hi",&(forcings[i]->flag));
		if(ReadLineError(valsread,1,"forcings flag"))	return NULL;

		if(forcings[i]->flag == 1 || forcings[i]->flag == 2 || forcings[i]->flag == 4 || forcings[i]->flag == 6 || forcings[i]->flag == 8)
		{
			forcings[i]->filename = (char*) malloc(string_size*sizeof(char));
			valsread = sscanf(linebuffer,"%*i %s",forcings[i]->filename);
			if(ReadLineError(valsread,1,"forcing data filename"))	return NULL;
			if(forcings[i]->flag == 1 && !CheckFilenameExtension(forcings[i]->filename,".str"))	return NULL;
			if(forcings[i]->flag == 4 && !CheckFilenameExtension(forcings[i]->filename,".ustr"))	return NULL;

			if(forcings[i]->flag == 2 || forcings[i]->flag == 6)
			{
				ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
				valsread = sscanf(linebuffer,"%u %lf %u %u",&(forcings[i]->increment),&(forcings[i]->file_time),&(forcings[i]->first_file),&(forcings[i]->last_file));
				if(ReadLineError(valsread,4,"time increment, file time, first file, and last file"))	return NULL;
			}
			else if(forcings[i]->flag == 8)
			{
				ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
				valsread = sscanf(linebuffer,"%u %u %u",&(forcings[i]->increment),&(forcings[i]->first_file),&(forcings[i]->last_file));
				if(ReadLineError(valsread,3,"time increment, first file, and last file"))	return NULL;
			}
		}
		else if(forcings[i]->flag == 3)	//Database
		{
			valsread = sscanf(linebuffer,"%*i %s",db_filename);
			if(ReadLineError(valsread,1,".dbc for forcing"))	return NULL;
			if(!CheckFilenameExtension(db_filename,".dbc"))	return NULL;
			db_connections[ASYNCH_DB_LOC_FORCING_START+i] = ReadDBC(db_filename,string_size);

			ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
			valsread = sscanf(linebuffer,"%u %lf %u %u",&(forcings[i]->increment),&(forcings[i]->file_time),&(forcings[i]->first_file),&(forcings[i]->last_file));
			if(ReadLineError(valsread,4,"chunk size, time resolution, beginning unix time, and ending unix time"))	return NULL;
			forcings[i]->raindb_start_time = forcings[i]->first_file;
		}
		else if(forcings[i]->flag == 9)	//Database with irregular timesteps
		{
			valsread = sscanf(linebuffer,"%*i %s",db_filename);
			if(ReadLineError(valsread,1,".dbc for forcing"))	return NULL;
			if(!CheckFilenameExtension(db_filename,".dbc"))	return NULL;
			db_connections[ASYNCH_DB_LOC_FORCING_START+i] = ReadDBC(db_filename,string_size);

			ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
			valsread = sscanf(linebuffer,"%u %u %u",&(forcings[i]->increment),&(forcings[i]->first_file),&(forcings[i]->last_file));
			if(ReadLineError(valsread,3,"time increment, file time, first file, and last file"))	return NULL;
			forcings[i]->raindb_start_time = forcings[i]->first_file;

			forcings[i]->lastused_first_file = forcings[i]->lastused_last_file = 0;
			forcings[i]->next_timestamp = forcings[i]->first_file;
			forcings[i]->number_timesteps = 0;
		}
		else if(forcings[i]->flag == 7)	//Recurring
		{
			forcings[i]->filename = (char*) malloc(string_size*sizeof(char));
			valsread = sscanf(linebuffer,"%*i %s",forcings[i]->filename);
			if(ReadLineError(valsread,1,"recurring rainfall filename"))	return NULL;
			if(!CheckFilenameExtension(forcings[i]->filename,".mon"))	return NULL;

			ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
			valsread = sscanf(linebuffer,"%u %u",&(forcings[i]->first_file),&(forcings[i]->last_file));
			if(ReadLineError(valsread,2,"first time, and last time"))	return NULL;
			forcings[i]->raindb_start_time = forcings[i]->first_file;
		}
		else if(forcings[i]->flag == 0) //No forcing
		{
			forcings[i]->filename = NULL;
		}
		else
		{
			printf("[%i]: Error reading %s: Invalid forcing flag %i.\n",my_rank,globalfilename,forcings[i]->flag);
			return NULL;
		}
	}

	//Grab the .dam filename
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%hu",&(GlobalVars->dam_flag));
	if(ReadLineError(valsread,1,"dam flag"))	return NULL;

	if(GlobalVars->dam_flag)
	{
		GlobalVars->dam_filename = (char*) malloc(string_size*sizeof(char));
		valsread = sscanf(linebuffer,"%*u %s",GlobalVars->dam_filename);
		if(ReadLineError(valsread,1,"filename for dam info"))	return NULL;
		if(GlobalVars->dam_flag == 1 && !CheckFilenameExtension(GlobalVars->dam_filename,".dam"))	return NULL;
		if(GlobalVars->dam_flag == 2 && !CheckFilenameExtension(GlobalVars->dam_filename,".qvs"))	return NULL;

		if(GlobalVars->dam_flag == 3)
		{
			if(!CheckFilenameExtension(GlobalVars->dam_filename,".dbc"))	return NULL;
			db_connections[ASYNCH_DB_LOC_QVS] = ReadDBC(GlobalVars->dam_filename,string_size);
		}
	}
	else
		GlobalVars->dam_filename = NULL;

	//Get the link ids where reservoirs exist
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%hu",&(GlobalVars->res_flag));
	if(ReadLineError(valsread,1,"res flag"))	return NULL;

	if(GlobalVars->res_flag)
	{
		if(GlobalVars->res_flag == 1)
		{
			GlobalVars->rsv_filename = (char*) malloc(string_size*sizeof(char));
			valsread = sscanf(linebuffer,"%*u %s %hi",GlobalVars->rsv_filename,&(GlobalVars->res_forcing_idx));
			if(ReadLineError(valsread,2,".rsv filename"))	return NULL;
			if(!CheckFilenameExtension(GlobalVars->rsv_filename,".rsv"))	return NULL;
		}
		else	//Flag is 2
		{
			GlobalVars->rsv_filename = (char*) malloc(string_size*sizeof(char));
			valsread = sscanf(linebuffer,"%*u %s %hi",GlobalVars->rsv_filename,&(GlobalVars->res_forcing_idx));
			if(ReadLineError(valsread,2,".dbc for reservoirs"))	return NULL;
			if(!CheckFilenameExtension(GlobalVars->rsv_filename,".dbc"))	return NULL;
			db_connections[ASYNCH_DB_LOC_RSV] = ReadDBC(GlobalVars->rsv_filename,string_size);
		}

		if(GlobalVars->res_forcing_idx >= GlobalVars->num_forcings)
		{
			if(my_rank == 0)	printf("Bad forcing index for a reservoir feed (%hi). Only %i forcings available.\n",GlobalVars->res_forcing_idx,GlobalVars->num_forcings);
			return NULL;
		}
	}
	else
	{
		GlobalVars->rsv_filename = NULL;
		GlobalVars->res_forcing_idx = -1;
	}

	//Grab where to write the hydrographs
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%hu",&(GlobalVars->hydros_loc_flag));
	if(ReadLineError(valsread,1,"hydrographs location"))	return NULL;

	GlobalVars->hydros_loc_filename = NULL;
	GlobalVars->hydro_table = NULL;

	if(GlobalVars->hydros_loc_flag == 1 || GlobalVars->hydros_loc_flag == 2 || GlobalVars->hydros_loc_flag == 4)
	{
		GlobalVars->hydros_loc_filename = (char*) malloc(string_size*sizeof(char));
		valsread = sscanf(linebuffer,"%*u %lf %s",&(GlobalVars->print_time),GlobalVars->hydros_loc_filename);
		if(ReadLineError(valsread,2,"hydrographs location"))	return NULL;
		if(GlobalVars->hydros_loc_flag == 1 && !CheckFilenameExtension(GlobalVars->hydros_loc_filename,".dat"))	return NULL;
		if(GlobalVars->hydros_loc_flag == 2 && !CheckFilenameExtension(GlobalVars->hydros_loc_filename,".csv"))	return NULL;
		if(GlobalVars->hydros_loc_flag == 4 && !CheckFilenameExtension(GlobalVars->hydros_loc_filename,".rad"))	return NULL;
		//GlobalVars->output_flag = (GlobalVars->hydros_loc_flag == 1) ? 0 : 1;

		if(GlobalVars->hydros_loc_flag == 1)	RemoveSuffix(GlobalVars->hydros_loc_filename,".dat");
		else if(GlobalVars->hydros_loc_flag == 2)	RemoveSuffix(GlobalVars->hydros_loc_filename,".csv");
		else if(GlobalVars->hydros_loc_flag == 4)	RemoveSuffix(GlobalVars->hydros_loc_filename,".rad");
	}
	else if(GlobalVars->hydros_loc_flag == 3)
	{
		GlobalVars->hydros_loc_filename = (char*) malloc(string_size*sizeof(char));
		GlobalVars->hydro_table = (char*) malloc(string_size*sizeof(char));
		valsread = sscanf(linebuffer,"%*u %lf %s %s",&(GlobalVars->print_time),GlobalVars->hydros_loc_filename,GlobalVars->hydro_table);
		if(ReadLineError(valsread,3,"hydrographs location"))	return NULL;
		if(!CheckFilenameExtension(GlobalVars->hydros_loc_filename,".dbc"))	return NULL;
		db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT] = ReadDBC(GlobalVars->hydros_loc_filename,string_size);
	}

	//Grab where to write the peakflow data
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%hu",&(GlobalVars->peaks_loc_flag));
	if(ReadLineError(valsread,1,"peakflow location"))	return NULL;

	if(GlobalVars->peaks_loc_flag == 0)
	{
		GlobalVars->peaks_loc_filename = NULL;
	}
	else if(GlobalVars->peaks_loc_flag == 1)
	{
		GlobalVars->peaks_loc_filename = (char*) malloc(string_size*sizeof(char));
		valsread = sscanf(linebuffer,"%*u %s",GlobalVars->peaks_loc_filename);
		if(ReadLineError(valsread,1,"peakflow location"))	return NULL;
		if(!CheckFilenameExtension(GlobalVars->peaks_loc_filename,".pea"))	return NULL;

		RemoveSuffix(GlobalVars->peaks_loc_filename,".pea");
	}
	else if(GlobalVars->peaks_loc_flag == 2)
	{
		GlobalVars->peaks_loc_filename = (char*) malloc(string_size*sizeof(char));
		GlobalVars->peak_table = (char*) malloc(string_size*sizeof(char));
		valsread = sscanf(linebuffer,"%*u %s %s",GlobalVars->peaks_loc_filename,GlobalVars->peak_table);
		if(ReadLineError(valsread,2,"peakflow location"))	return NULL;
		if(!CheckFilenameExtension(GlobalVars->peaks_loc_filename,".dbc"))	return NULL;
		db_connections[ASYNCH_DB_LOC_PEAK_OUTPUT] = ReadDBC(GlobalVars->peaks_loc_filename,string_size);
	}

	//Grab the .sav files
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%hu",&(GlobalVars->hydrosave_flag));
	if(ReadLineError(valsread,1,"hydrographs save flag"))	return NULL;

	if(GlobalVars->hydrosave_flag == 1)
	{
		GlobalVars->hydrosave_filename = (char*) malloc(string_size*sizeof(char));
		valsread = sscanf(linebuffer,"%*u %s",GlobalVars->hydrosave_filename);
		if(ReadLineError(valsread,1,"hydrographs .sav filename"))	return NULL;
		if(!CheckFilenameExtension(GlobalVars->hydrosave_filename,".sav"))	return NULL;
	}
	else if(GlobalVars->hydrosave_flag == 2)
	{
		valsread = sscanf(linebuffer,"%*u %s",db_filename);
		if(ReadLineError(valsread,1,".dbc for hydrograph save ids"))	return NULL;
		if(!CheckFilenameExtension(db_filename,".dbc"))	return NULL;
		GlobalVars->hydrosave_filename = NULL;
		db_connections[ASYNCH_DB_LOC_HYDROSAVE] = ReadDBC(db_filename,string_size);
	}
	else
		GlobalVars->hydrosave_filename = NULL;

	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%hu",&(GlobalVars->peaksave_flag));
	if(ReadLineError(valsread,1,"peakflows save flag"))	return NULL;

	if(GlobalVars->peaksave_flag == 1)
	{
		GlobalVars->peaksave_filename = (char*) malloc(string_size*sizeof(char));
		valsread = sscanf(linebuffer,"%*u %s",GlobalVars->peaksave_filename);
		if(ReadLineError(valsread,1,"peakflows .sav filename"))	return NULL;
		if(!CheckFilenameExtension(GlobalVars->peaksave_filename,".sav"))	return NULL;
	}
	else if(GlobalVars->peaksave_flag == 2)
	{
		valsread = sscanf(linebuffer,"%*u %s",db_filename);
		if(ReadLineError(valsread,1,".dbc for peakflow save ids"))	return NULL;
		if(!CheckFilenameExtension(db_filename,".dbc"))	return NULL;
		GlobalVars->peaksave_filename = NULL;
		db_connections[ASYNCH_DB_LOC_PEAKSAVE] = ReadDBC(db_filename,string_size);
	}
	else
		GlobalVars->peaksave_filename = NULL;
	GlobalVars->peakfilename = NULL;

	//Grab data dump info
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%hu",&(GlobalVars->dump_loc_flag));
	if(ReadLineError(valsread,1,"snapshot save flag"))	return NULL;

	GlobalVars->dump_loc_filename = NULL;
	GlobalVars->dump_table = NULL;

	if(GlobalVars->dump_loc_flag == 1 || GlobalVars->dump_loc_flag == 3)
	{
		GlobalVars->dump_loc_filename = (char*) malloc(string_size*sizeof(char));
		valsread = sscanf(linebuffer,"%*u %s",GlobalVars->dump_loc_filename);
		if(ReadLineError(valsread,1,"snapshot filename"))	return NULL;
        if ((GlobalVars->dump_loc_flag == 1) && !CheckFilenameExtension(GlobalVars->dump_loc_filename, ".rec"))
		    return NULL;
        if ((GlobalVars->dump_loc_flag == 3) && !CheckFilenameExtension(GlobalVars->dump_loc_filename, ".h5"))
            return NULL;
	}
	else if(GlobalVars->dump_loc_flag == 2)
	{
		GlobalVars->dump_table = (char*) malloc(string_size*sizeof(char));
		valsread = sscanf(linebuffer,"%*u %s %s",db_filename,GlobalVars->dump_table);
		if(ReadLineError(valsread,2,".dbc for snapshots"))	return NULL;
		if(!CheckFilenameExtension(db_filename,".dbc"))	return NULL;
		//GlobalVars->dump_loc_filename = NULL;
		db_connections[ASYNCH_DB_LOC_SNAPSHOT_OUTPUT] = ReadDBC(db_filename,string_size);
	}
    else if (GlobalVars->dump_loc_flag == 4)
    {
        GlobalVars->dump_loc_filename = (char*)malloc(string_size * sizeof(char));
        valsread = sscanf(linebuffer, "%*u %lf %s", &GlobalVars->dump_time, GlobalVars->dump_loc_filename);
        if (ReadLineError(valsread, 2, "snapshot time and filename"))	return NULL;
        if (!CheckFilenameExtension(GlobalVars->dump_loc_filename, ".h5"))
            return NULL;
    }

	//Grab folder locations
	//GlobalVars->results_folder = (char*) malloc(string_size*sizeof(char));
	GlobalVars->temp_filename = (char*) malloc(string_size*sizeof(char));
	//ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	//valsread = sscanf(linebuffer,"%s",GlobalVars->results_folder);
	//if(ReadLineError(valsread,1,"results folder"))	return NULL;
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%s",GlobalVars->temp_filename);
	if(ReadLineError(valsread,1,"scratch work folder"))	return NULL;

	if(GlobalVars->print_par_flag)	//!!!! Is this needed? Why bother? !!!!
	{
		if(AttachParameters(GlobalVars->temp_filename,string_size,GlobalVars->global_params,string_size))
		{
			printf("[%i]: Error attaching global parameters to temporary filenames.\n",my_rank);
			return NULL;
		}
	}

	sprintf(db_filename,"_%i_%i",getpid(),my_rank);
	strcat(GlobalVars->temp_filename,db_filename);

	//Grab adapative data
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%lf %lf %lf",&((*GlobalErrors)->facmin),&((*GlobalErrors)->facmax),&((*GlobalErrors)->fac));
	if(ReadLineError(valsread,3,"facmin, facmax, fac"))	return NULL;

	//Read in the flag for the error tolerances
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%hi",&(GlobalVars->rkd_flag));
	if(ReadLineError(valsread,1,"error tolerance flag"))	return NULL;

	//Set some parameters
	GlobalVars->max_s = 0;
	GlobalVars->max_parents = 0;

	if(GlobalVars->rkd_flag == 0)	//Error data is found in the universal file
	{
		rkdfilename[0] = '\0';
		ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
		valsread = sscanf(linebuffer,"%u",&flag);
		if(ReadLineError(valsread,1,"RK method index"))	return NULL;
		GlobalVars->method = flag;

		//Count the number of states given
		ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
		total = 0;
		int error_dim = -1;
		double tempy;
		do
		{
			error_dim++;
			valsread = sscanf(&(linebuffer[total]),"%lf%n",&tempy,&written);	//Need to actually read to trigger valsread
			total += written;
		} while(valsread > 0);

		//Reserve memory
		(*GlobalErrors)->abstol = v_get(error_dim);
		(*GlobalErrors)->reltol = v_get(error_dim);
		(*GlobalErrors)->abstol_dense = v_get(error_dim);
		(*GlobalErrors)->reltol_dense = v_get(error_dim);

		//ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
		total = 0;
		for(i=0;i<error_dim;i++)	//Note: Don't read the line from disk again
		{
			valsread = sscanf(&(linebuffer[total]),"%lf%n",&((*GlobalErrors)->abstol.ve[i]),&written);
			//if(ReadLineError(valsread,1,"an abstol component"))	return NULL;
			total += written;
		}

		ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
		total = 0;
		for(i=0;i<error_dim;i++)
		{
			valsread = sscanf(&(linebuffer[total]),"%lf%n",&((*GlobalErrors)->reltol.ve[i]),&written);
			if(ReadLineError(valsread,1,"a reltol component"))	return NULL;
			total += written;
		}

		ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
		total = 0;
		for(i=0;i<error_dim;i++)
		{
			valsread = sscanf((&linebuffer[total]),"%lf%n",&((*GlobalErrors)->abstol_dense.ve[i]),&written);
			if(ReadLineError(valsread,1,"an abstol dense output component"))	return NULL;
			total += written;
		}

		ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
		total = 0;
		for(i=0;i<error_dim;i++)
		{
			valsread = sscanf(&(linebuffer[total]),"%lf%n",&((*GlobalErrors)->reltol_dense.ve[i]),&written);
			if(ReadLineError(valsread,1,"a reltol dense output component"))	return NULL;
			total += written;
		}
	}
	else if(GlobalVars->rkd_flag == 1)	//Error data is in an .rkd file
	{
		GlobalVars->method = -1;
		(*GlobalErrors)->abstol = v_get(0);
		(*GlobalErrors)->reltol = v_get(0);
		(*GlobalErrors)->abstol_dense = v_get(0);
		(*GlobalErrors)->reltol_dense = v_get(0);

		//ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
		valsread = sscanf(linebuffer,"%*i %s",rkdfilename);
		if(ReadLineError(valsread,1,".rkd filename"))
            return NULL;
		if(!CheckFilenameExtension(rkdfilename,".rkd"))
            return NULL;
	}
	else
	{
		printf("Error: bad flag for error tolerance. Got %hu.\n",GlobalVars->rkd_flag);
		return NULL;
	}

	//Check for end mark
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	sscanf(linebuffer,"%c",&endmark);
	if(endmark != '#')
	{
		printf("Error: an ending # not seen in %s file. Got %c.\n",globalfilename,endmark);
		return NULL;
	}

	//Setup an io object
	GlobalVars->output_data = BuildIO(GlobalVars);

	//Clean up
	free(db_filename);
	free(linebuffer);
	if(my_rank == 0)	fclose(globalfile);
	return GlobalVars;
}

//Reads in the contents of a .sav file.
//char filename[]: The filename of the .sav file.
//int N: The number of links in the river system.
//int** save_list (set by this method): Array of link ids where data will be written.
//int* size (set by this method): Will be the number of links for which data must be written to disk (number of links in the .sav file).
//Returns 1 if an error occured. 0 otherwise.
int Create_SAV_Data(char filename[],Link** sys,unsigned int N,unsigned int** save_list,unsigned int* size,ConnData *conninfo,unsigned short int flag)
{
	unsigned int i,id,error = 0,list_size = N;
	FILE* save_file = NULL;
	*size = 0;
	*save_list = NULL;

	//.sav file
	if(flag == 1)
	{
		if(my_rank == 0)
		{
			//Open save file
			save_file = fopen(filename,"r");
			if(!save_file)
			{
				printf("Error opening .sav file %s.\n",filename);
				error = 1;
			}
			else
			{
				if(CheckWinFormat(save_file))
				{
					printf("Error: File %s appears to be in Windows format. Try converting to unix format using 'dos2unix' at the command line.\n",filename);
					fclose(save_file);
					return 1;
				}

				*save_list = (unsigned int*) malloc(list_size*sizeof(unsigned int));

				//Read the save file
				while( fscanf(save_file,"%u",&id) != EOF )
				{
					(*save_list)[*size] = id;
					(*size)++;
					if((*size) == list_size)
					{
						list_size *= 2;
						*save_list = (unsigned int*) realloc(*save_list,list_size*sizeof(unsigned int));
					}
				}

				//Close save file
				fclose(save_file);
			}
		}

		//Send data
		MPI_Bcast(&error,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
		if(error)	return 1;
		MPI_Bcast(size,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
		*save_list = realloc(*save_list,*size * sizeof(unsigned int));
		MPI_Bcast(*save_list,*size,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	}
	else if(flag == 2)	//Grab from database
	{
		//char* query = conninfo->query;
		PGresult *res;

		if(my_rank == 0)
		{
			ConnectPGDB(conninfo);
			sprintf(conninfo->query,conninfo->queries[0]);
			res = PQexec(conninfo->conn,conninfo->query);
			error = CheckResError(res,"locating links with sensors");

			if(!error)
			{
				*size = PQntuples(res);
				*save_list = malloc(*size * sizeof(unsigned int));
				for(i=0;i<*size;i++)	(*save_list)[i] = atoi(PQgetvalue(res,i,0));
			}

			PQclear(res);
			DisconnectPGDB(conninfo);
		}

		MPI_Bcast(&error,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
		if(error)	return 1;
		MPI_Bcast(size,1,MPI_INT,0,MPI_COMM_WORLD);
		if(my_rank != 0)	*save_list = malloc(*size * sizeof(unsigned int));
		MPI_Bcast(*save_list,*size,MPI_INT,0,MPI_COMM_WORLD);
	}
	else if(flag == 3)	//All links
	{
		*size = N;
		*save_list = malloc(*size * sizeof(unsigned int));
		for(i=0;i<N;i++)	(*save_list)[i] = sys[i]->ID;
	}

	return 0;
}


void ReadLineFromTextFile(FILE* globalfile,char* linebuffer,unsigned int size,unsigned int string_size)
{
    size_t linebuffer_length;
	if(my_rank == 0)
	{
		linebuffer[0] = '%';
		while(!feof(globalfile) && (linebuffer[0] == '%' || linebuffer[0] == '\n'))	fgets(linebuffer,size,globalfile);
		linebuffer_length = strlen(linebuffer);
		if(linebuffer_length > string_size + 10)	printf("Warning: %zu %u Line in .gbl file may be too long. Read in the long line:\n\"%s\"\n",strlen(linebuffer),string_size,linebuffer);
	}

	MPI_Bcast(&linebuffer_length,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Bcast(linebuffer,linebuffer_length+1,MPI_CHAR,0,MPI_COMM_WORLD);
}

int ReadLineError(int valsread,int valswant,char message[])
{
	if(valsread < valswant)
	{
		if(my_rank == 0)	printf("Error: Did not get a value from .gbl file for %s. %i\n",message,valsread);
		return 1;
	}
	return 0;
}

//Removes a suffix from filename, if present.
//Returns 1 if suffix removed
//0 if not (not present)
int RemoveSuffix(char* filename, const char* suffix)
{
    size_t filename_length = strlen(filename);
    size_t suffix_length = strlen(suffix);

	if(suffix_length > filename_length)	return 0;

    char *dot = strrchr(filename, '.');
    if (!dot || dot == filename || strcmp(dot, suffix) != 0)
        return 0;

    *dot = '\0';
	return 1;
}

//Put a vector of global_params onto the end of a filename.
//Returns 1 if filename is not long enough to support this.
//Retruns 0 if the parameters are attached.
int AttachParameters(char* filename,unsigned int max_size,VEC v,unsigned int string_size)
{
	unsigned int i,count,total=0;
	char *buffer;
    size_t length = strlen(filename);

    buffer = (char *) malloc(string_size);

	for(i=0;i<v.dim;i++)
	{
		sprintf(buffer,"_%.4e%n",v.ve[i],&count);
		total += count;
		if(count+1 > string_size)	return 1;
		if(total+1 > max_size)		return 1;
		strcat(filename,buffer);
	}

    free(buffer);

	return 0;
}


//Checks that the filename ends with the approriate extension.
//Returns 1 if the extension is correct, 0 if not.
int CheckFilenameExtension(char* filename,char* extension)
{
    size_t l_filename,l_extension,i;
	
	if(!extension)	return 1;
	if(!filename)
	{
		if(my_rank == 0) printf("Error: no filename given.\n");
		return 0;
	}

	l_filename = strlen(filename);
	l_extension = strlen(extension);

	if(l_filename < l_extension)
	{
		if(my_rank == 0) printf("Error: filename %s is of wrong file type (expected %s).\n",filename,extension);
		return 0;
	}
	for(i=0;i<l_extension;i++)
	{
		if(extension[i] != filename[l_filename - l_extension + i])
		{
			if(my_rank == 0) printf("Error: filename %s is of wrong file type (expected %s).\n",filename,extension);
			return 0;
		}
	}

	return 1;
}

//Returns 1 if the ascii file is in Windows format. 0 otherwise.
int CheckWinFormat(FILE* file)
{
	char buffer[256];
	fpos_t pos;
	fgetpos(file,&pos);

	if(!file)	return 0;
	fgets(buffer,256,file);
	if(buffer[strlen(buffer)-2] == 13)	//CR
	{
		fsetpos(file,&pos);
		return 1;
	}
	fsetpos(file,&pos);
	return 0;
}

//Finds the path in filename. The path is copied into variable path.
//Returns 0 if the path is is extracted successfully.
//Returns 1 if filename does not contain a path (path variable set to empty).
//Returns 2 if filename is just a path (the path variable is still set).
int FindPath(char* filename,char* path)
{
    size_t i;
    size_t len = strlen(filename);
	char holder;

	if(len == 0)
	{
		path[0] = '\0';
		return 1;
	}

	if(filename[len-1] == '/')
	{
		strcpy(path,filename);
		return 2;
	}

	for(i=len-2;i>=0;i--)
	{
		if(filename[i] == '/')
		{
			holder = filename[i+1];
			filename[i+1] = '\0';
			strcpy(path,filename);
			filename[i+1] = holder;
			return 0;
		}
	}

	path[0] = '\0';
	return 1;
}

//Finds the filename from the full path. The filename is copied into variable filename.
//Returns 0 if the filename is is extracted successfully.
//Returns 1 if fullpath is just a path (the filename variable is set to empty).
int FindFilename(char* fullpath,char* filename)
{
    size_t i;
    size_t len = strlen(fullpath);

	if(len == 0 || fullpath[len-1] == '/')
	{
		filename[0] = '\0';
		return 1;
	}

	for(i=len-2;i>=0;i--)
	{
		if(fullpath[i] == '/')
			break;
	}

	i++;
	sprintf(filename,"%s",&(fullpath[i]));
	return 0;
}

