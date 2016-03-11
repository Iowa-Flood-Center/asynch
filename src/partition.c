#include "partition.h"

//Partitions a river system by first partitioning the leaves.
//Link** sys: The river system.
//int N: The number of links in the river system.
//Link** leaves: The list of leaves to the system.
//int numleaves: Number of leaves in the system.
//int** my_sys (set by this method): An array that will contain the location in the system of each link assigned to this process.
//int* my_N (set by this method): Will be the number of links assigned to this process.
//TransData* my_data (this is assumed to be allocated already, but contents will be set here): Information about the links that
//				each process will communicate information about will be stored in my_data.
//short int *getting (set by this method, but assumed space is allocated with N entries): getting[i] will equal 1 if this process needs
//				to receive information about link i, 0 otherwise.
//Returns an array of N integers that take values from 0 to np-1. The i-th entry contains which process the link in location i of
//				the system is assigned to.
int* Partition_System_By_Leaves(Link** sys,unsigned int N,Link** leaves,unsigned int numleaves,unsigned int** my_sys,unsigned int* my_N,TransData* my_data,short int *getting)
{
	//Assign leaves
	unsigned int i,start_index,end_index,extras;
	int j;
	unsigned int nodes_per_proc = numleaves/np;	//Number of leaves assigned to each process (except the last)

	start_index = nodes_per_proc * my_rank;
	if(my_rank == np-1)		end_index = numleaves;
	else				end_index = nodes_per_proc * (my_rank + 1);
	*my_N = end_index - start_index;
	unsigned int my_max_nodes = N - numleaves + nodes_per_proc;
	*my_sys = (unsigned int*) malloc(my_max_nodes * sizeof(unsigned int));	//The indices of this processes links (in sys)
	for(i=0;i<my_max_nodes;i++)	(*my_sys)[i] = -1;
	for(i=start_index;i<end_index;i++)	(*my_sys)[i-start_index] = leaves[i]->location;

	//Initialize assignments
	int* assignments = (int*) malloc(N*sizeof(int));
	for(i=0;i<N;i++)	assignments[i] = -1;

	//Calculate and store the assignments for the leaves
	extras = numleaves % np;	//This is how many extra leaves we have (they are assigned to proc np above). We assign them to np-1.
	end_index = numleaves - extras;
	for(i=0;i<end_index;i++)
		assignments[leaves[i]->location] = i/nodes_per_proc;
	for(i=end_index;i<numleaves;i++)
		assignments[leaves[i]->location] = np-1;

	//Initialize getting
	for(i=0;i<N;i++)	getting[i] = 0;

	//Assign the rest of the links
	//I do this by starting at each leaf, and assigning the leaf's child to its left parent's (0) process.
	Link *current,*q;
	unsigned int curr_idx;
	for(i=0;i<numleaves;i++)
	{
		current = leaves[i]->c;
		if(current != NULL)
		{
			curr_idx = current->location;
			q = leaves[i];
			while(assignments[curr_idx] == -1)
			{
				assignments[curr_idx] = assignments[current->parents[0]->location];
				if(assignments[curr_idx] == my_rank)	//If this node is assigned to this process
				{
					(*my_sys)[*my_N] = curr_idx;
					(*my_N)++;
				}
				q = current;
				current = current->c;
				if(current != NULL)	curr_idx = current->location;
				else			break;
			}

			if(current != NULL)
			{
				if(assignments[curr_idx] != assignments[q->location])
				{
					int q_proc = assignments[q->location];
					int curr_proc = assignments[curr_idx];

					//Check if this process will receive data
					if(my_rank == curr_proc)
					{
						//my_data->receive_data[q_proc][my_data->receive_size[q_proc]] = q;
						(my_data->receive_size[q_proc])++;
						getting[q->location] = 1;
					}

					//Check if this process will send data
					if(my_rank == q_proc)
					{
						//my_data->send_data[curr_proc][my_data->send_size[curr_proc]] = q;
						(my_data->send_size[curr_proc])++;
					}
				}
			}
		}
	}

	//Reorder my_sys so that the links with lower numbering are towards the beginning
	merge_sort_distance(sys,*my_sys,*my_N);

	//Allocate space in my_data for recieving and sending
	for(j=0;j<np;j++)
	{
		my_data->receive_data[j] = (Link**) malloc(my_data->receive_size[j] * sizeof(Link*));
		my_data->send_data[j] = (Link**) malloc(my_data->send_size[j] * sizeof(Link*));
	}
	int* current_receive_size = (int*) calloc(np,sizeof(int));
	int* current_send_size = (int*) calloc(np,sizeof(int));
	for(i=0;i<N;i++)	assignments[i] = -1;

	//Recalculate and store the assignments for the leaves
	extras = numleaves % np;
	end_index = numleaves - extras;
	for(i=0;i<end_index;i++)
		assignments[leaves[i]->location] = i/nodes_per_proc;
	for(i=end_index;i<numleaves;i++)
		assignments[leaves[i]->location] = np-1;

	//Assign links to receive_data and send_data
	for(i=0;i<numleaves;i++)
	{
		current = leaves[i]->c;
		if(current != NULL)
		{
			curr_idx = current->location;
			q = leaves[i];

			while(assignments[curr_idx] == -1)
			{
				assignments[curr_idx] = assignments[current->parents[0]->location];
				q = current;
				current = current->c;
				if(current != NULL)	curr_idx = current->location;
				else			break;
			}

			if(current != NULL)
			{
				if(assignments[curr_idx] != assignments[q->location])
				{
					int q_proc = assignments[q->location];
					int curr_proc = assignments[curr_idx];

					//Check if this process will receive data
					if(my_rank == curr_proc)
					{
						my_data->receive_data[q_proc][current_receive_size[q_proc]] = q;
						current_receive_size[q_proc]++;
					}

					//Check if this process will send data
					if(my_rank == q_proc)
					{
						my_data->send_data[curr_proc][current_send_size[curr_proc]] = q;
						current_send_size[curr_proc]++;
					}
				}
			}
		}
	}

	//Clean up
	free(current_receive_size);
	free(current_send_size);
	*my_sys = (unsigned int*) realloc(*my_sys,*my_N*sizeof(unsigned int));

	return assignments;
}


//Partitions a river system by first partitioning the leaves. Makes adjustments so all leaves are assigned to their child's process.
//Link** sys: The river system.
//int N: The number of links in the river system.
//Link** leaves: The list of leaves to the system.
//int numleaves: Number of leaves in the system.
//int** my_sys (set by this method): An array that will contain the location in the system of each link assigned to this process.
//int* my_N (set by this method): Will be the number of links assigned to this process.
//TransData* my_data (this is assumed to be allocated already, but contents will be set here): Information about the links that
//				each process will communicate information about will be stored in my_data.
//short int *getting (set by this method, but assumed space is allocated with N entries): getting[i] will equal 1 if this process needs
//				to receive information about link i, 0 otherwise.
//Returns an array of N integers that take values from 0 to np-1. The i-th entry contains which process the link in location i of
//				the system is assigned to.
int* Partition_System_By_Leaves_2(Link** sys,unsigned int N,Link** leaves,unsigned int numleaves,unsigned int** my_sys,unsigned int* my_N,TransData* my_data,short int *getting)
{
	//Assign leaves
	unsigned int i,j,k,l,start_index,end_index,extras;
	int ii;
	unsigned int nodes_per_proc = numleaves/np;	//Number of leaves assigned to each process (except the last)

	start_index = nodes_per_proc * my_rank;
	if(my_rank == np-1)		end_index = numleaves;
	else				end_index = nodes_per_proc * (my_rank + 1);
	*my_N = end_index - start_index;
	unsigned my_max_nodes = N - numleaves + nodes_per_proc;
	*my_sys = (unsigned int*) malloc(my_max_nodes * sizeof(int));	//The indices of this processes links (in sys)
	for(i=0;i<my_max_nodes;i++)	(*my_sys)[i] = -1;
	for(i=start_index;i<end_index;i++)	(*my_sys)[i-start_index] = leaves[i]->location;

	//Initialize assignments
	int* assignments = (int*) malloc(N*sizeof(int));
	for(i=0;i<N;i++)	assignments[i] = -1;

	//Calculate and store the assignments for the leaves
	extras = numleaves % np;	//This is how many extra leaves we have (they are assigned to proc np above). We assign them to np-1.
	end_index = numleaves - extras;
	for(i=0;i<end_index;i++)
		assignments[leaves[i]->location] = i/nodes_per_proc;
	for(i=end_index;i<numleaves;i++)
		assignments[leaves[i]->location] = np-1;

	//Initialize getting
	for(i=0;i<N;i++)	getting[i] = 0;

	//Assign the rest of the links
	//I do this by starting at each leaf, and assigning the leaf's child to its left parent's (0) process.
	Link *current,*q;
	int curr_idx;
	for(i=0;i<numleaves;i++)
	{
		current = leaves[i]->c;
		if(current != NULL)
		{
			curr_idx = current->location;
			q = leaves[i];

			while(assignments[curr_idx] == -1)
			{
				assignments[curr_idx] = assignments[current->parents[0]->location];
				if(assignments[curr_idx] == my_rank)	//If this node is assigned to this process
				{
					(*my_sys)[*my_N] = curr_idx;
					(*my_N)++;
				}

				//Check if any parents are leaves. If so, assign them to this process
				for(j=0;j<current->numparents;j++)
				{
					if(current->parents[j]->numparents == 0 && assignments[current->parents[j]->location] != assignments[curr_idx])
					{
						//Remove the parent from the my_sys it is currently in
						if(assignments[current->parents[j]->location] == my_rank)
						{
							for(k=0;k<*my_N;k++)
							{
								if((*my_sys)[k] == current->parents[j]->location)
								{
									for(l=k;l<*my_N-1;l++)
										(*my_sys)[l] = (*my_sys)[l+1];
									break;
								}
							}
							(*my_N)--;
						}

						//Assign the parent
						assignments[current->parents[j]->location] = assignments[current->location];
						if(assignments[current->parents[j]->location] == my_rank)
						{
							(*my_sys)[*my_N] = current->parents[j]->location;
							(*my_N)++;
						}
					}
				}

				q = current;
				current = current->c;
				if(current != NULL)	curr_idx = current->location;
				else			break;
			}

			if(current != NULL)
			{
				if(assignments[curr_idx] != assignments[q->location])
				{
					int q_proc = assignments[q->location];
					int curr_proc = assignments[curr_idx];

					//Check if this process will receive data
					if(my_rank == curr_proc)
					{
						(my_data->receive_size[q_proc])++;
						getting[q->location] = 1;
					}

					//Check if this process will send data
					if(my_rank == q_proc)	(my_data->send_size[curr_proc])++;
				}
			}
		}
	}

	//Reorder my_sys so that the links with lower numbering are towards the beginning
	merge_sort_distance(sys,*my_sys,*my_N);

	//Allocate space in my_data for recieving and sending
	for(ii=0;ii<np;ii++)
	{
		my_data->receive_data[ii] = (Link**) malloc(my_data->receive_size[ii] * sizeof(Link*));
		my_data->send_data[ii] = (Link**) malloc(my_data->send_size[ii] * sizeof(Link*));
	}
	int* current_receive_size = (int*) calloc(np,sizeof(int));
	int* current_send_size = (int*) calloc(np,sizeof(int));
	for(i=0;i<N;i++)	assignments[i] = -1;

	//Calculate and store the assignments for the leaves
	extras = numleaves % np;	//This is how many extra leaves we have (they are assigned to proc np above). We assign them to np-1.
	end_index = numleaves - extras;
	for(i=0;i<end_index;i++)
		assignments[leaves[i]->location] = i/nodes_per_proc;
	for(i=end_index;i<numleaves;i++)
		assignments[leaves[i]->location] = np-1;

	for(i=0;i<numleaves;i++)
	{
		current = leaves[i]->c;
		if(current != NULL)
		{
			curr_idx = current->location;
			q = leaves[i];

			while(assignments[curr_idx] == -1)
			{
				assignments[curr_idx] = assignments[current->parents[0]->location];

				//Check if any parents are leaves. If so, assign them to this process
				for(j=0;j<current->numparents;j++)
				{
					if(current->parents[j]->numparents == 0 && assignments[current->parents[j]->location] != assignments[curr_idx])
					{
						//Assign the parent
						assignments[current->parents[j]->location] = assignments[current->location];
					}
				}

				q = current;
				current = current->c;
				if(current != NULL)	curr_idx = current->location;
				else			break;
			}

			if(current != NULL)
			{
				if(assignments[curr_idx] != assignments[q->location])
				{
					int q_proc = assignments[q->location];
					int curr_proc = assignments[curr_idx];

					//Check if this process will receive data
					if(my_rank == curr_proc)
					{
						my_data->receive_data[q_proc][current_receive_size[q_proc]] = q;
						current_receive_size[q_proc]++;
					}

					//Check if this process will send data
					if(my_rank == q_proc)
					{
						my_data->send_data[curr_proc][current_send_size[curr_proc]] = q;
						current_send_size[curr_proc]++;
					}
				}
			}
		}
	}

	//Clean up
	free(current_receive_size);
	free(current_send_size);
	(*my_sys) = (unsigned int*) realloc(my_sys,*my_N*sizeof(unsigned int));

	return assignments;
}

/*
int* Partition_METIS_Traditional(Link** sys,unsigned int N,Link** leaves,unsigned int numleaves,unsigned int** my_sys,unsigned int* my_N,TransData* my_data,short int *getting,UnivVars* GlobalVars)
{
	unsigned int i,j,start_index,end_index,extras,partition,loc,retval;
	unsigned int nodes_per_proc = numleaves/np;	//Number of leaves assigned to each process (except the last)
	Link* current;

	start_index = nodes_per_proc * my_rank;
	if(my_rank == np-1)		end_index = numleaves;
	else				end_index = nodes_per_proc * (my_rank + 1);
//	*my_N = end_index - start_index;
	*my_N = 0;
	unsigned int my_max_nodes = N - numleaves + nodes_per_proc;
	*my_sys = (unsigned int*) malloc(my_max_nodes * sizeof(unsigned int));	//The indices of this processes links (in sys)
	for(i=0;i<my_max_nodes;i++)	(*my_sys)[i] = -1;
	for(i=start_index;i<end_index;i++)	(*my_sys)[i-start_index] = leaves[i]->location;
	for(i=0;i<N;i++)	getting[i] = 0;

	//Initialize assignments
	int* assignments = (int*) malloc(N*sizeof(int));
	for(i=0;i<N;i++)	assignments[i] = -1;

	//Form the graph to partition
	idx_t* xadj = malloc((N+1)*sizeof(idx_t));
	idx_t* adjncy = malloc(2*(N-1)*sizeof(idx_t));
	idx_t index = 0;

	for(i=0;i<N;i++)
	{
		xadj[i] = index;
		current = sys[i];
		if(current->c != NULL)
		{
			adjncy[index] = current->c->location;
			index++;
		}
		for(j=0;j<current->numparents;j++)
		{
			adjncy[index] = current->parents[j]->location;
			index++;
		}
	}
	xadj[N] = 2*(N-1);

	//Partition the system
	idx_t nverts = N;
	idx_t parts = np;
	idx_t ncon = 1;
	idx_t objval;
	idx_t* partitions = calloc(N,sizeof(idx_t));
	if(np != 1)
	{
		retval = METIS_PartGraphKway(&nverts,&ncon,xadj,adjncy,NULL,NULL,NULL,&parts,NULL,NULL,NULL,&objval,partitions);
		if(retval != METIS_OK)
		{
			printf("Error: METIS returned error code %i.\n",retval);
			return NULL;
		}
	}

	*my_N = 0;
	for(i=0;i<N;i++)
	{
		assignments[i] = partitions[i];	//!!!! Just use assignments? !!!!
		if(partitions[i] == my_rank)
		{
			(*my_sys)[*my_N] = i;
			(*my_N)++;
		}
	}

	//Set the getting array and determine number of sending and receiving links
	for(i=0;i<*my_N;i++)
	{
		//Receiving
		for(j=0;j<sys[(*my_sys)[i]]->numparents;j++)
		{
			loc = sys[(*my_sys)[i]]->parents[j]->location;
			if(assignments[loc] != my_rank)
			{
				getting[loc] = 1;
				my_data->receive_size[assignments[loc]]++;
			}
		}

		//Sending
		if(sys[(*my_sys)[i]]->c != NULL)
		{
			loc = sys[(*my_sys)[i]]->c->location;
			if(assignments[loc] != my_rank)
				my_data->send_size[assignments[loc]]++;
		}
	}

	//Reorder my_sys so that the links with lower numbering are towards the beginning
	merge_sort_distance(sys,*my_sys,*my_N);

	//Allocate space in my_data for recieving and sending
	for(j=0;j<np;j++)
	{
		my_data->receive_data[j] = (Link**) malloc(my_data->receive_size[j] * sizeof(Link*));
		my_data->send_data[j] = (Link**) malloc(my_data->send_size[j] * sizeof(Link*));
	}

	//Set the receive_data and send_data arrays
	int* current_receive_size = (int*) calloc(np,sizeof(int));
	int* current_send_size = (int*) calloc(np,sizeof(int));
	for(i=0;i<*my_N;i++)
	{
		//Receiving
		for(j=0;j<sys[(*my_sys)[i]]->numparents;j++)
		{
			loc = sys[(*my_sys)[i]]->parents[j]->location;
			if(assignments[loc] != my_rank)
			{
				my_data->receive_data[assignments[loc]][current_receive_size[assignments[loc]]] = sys[loc];
				current_receive_size[assignments[loc]]++;
			}
		}

		//Sending
		if(sys[(*my_sys)[i]]->c != NULL)
		{
			loc = sys[(*my_sys)[i]]->c->location;
			if(assignments[loc] != my_rank)
			{
				my_data->send_data[assignments[loc]][current_send_size[assignments[loc]]] = sys[(*my_sys)[i]];
				current_send_size[assignments[loc]]++;
			}
		}
	}

	//Clean up
	free(current_receive_size);
	free(current_send_size);
	free(xadj);
	free(adjncy);
	free(partitions);
	(*my_sys) = (unsigned int*) realloc(my_sys,*my_N*sizeof(unsigned int));

	return assignments;
}


int* Partition_METIS_RainChanges(Link** sys,unsigned int N,Link** leaves,unsigned int numleaves,unsigned int** my_sys,unsigned int* my_N,TransData* my_data,short int *getting,UnivVars* GlobalVars)
{
	unsigned int i,j,start_index,end_index,extras,partition,loc,retval;
	unsigned int nodes_per_proc = numleaves/np;	//Number of leaves assigned to each process (except the last)
	char filename[256];
	Link* current;

	start_index = nodes_per_proc * my_rank;
	if(my_rank == np-1)		end_index = numleaves;
	else				end_index = nodes_per_proc * (my_rank + 1);
//	*my_N = end_index - start_index;
	*my_N = 0;
	unsigned int my_max_nodes = N;
	*my_sys = (unsigned int*) malloc(my_max_nodes * sizeof(unsigned int));	//The indices of this processes links (in sys)
	for(i=0;i<my_max_nodes;i++)	(*my_sys)[i] = -1;
	for(i=start_index;i<end_index;i++)	(*my_sys)[i-start_index] = leaves[i]->location;
	for(i=0;i<N;i++)	getting[i] = 0;

	//Initialize assignments
	int* assignments = (int*) malloc(N*sizeof(int));
	for(i=0;i<N;i++)	assignments[i] = -1;

	//Form the graph to partition
	idx_t* xadj = malloc((N+1)*sizeof(idx_t));
	idx_t* adjncy = malloc(2*(N-1)*sizeof(idx_t));
	idx_t index = 0;

	for(i=0;i<N;i++)
	{
		xadj[i] = index;
		current = sys[i];
		if(current->c != NULL)
		{
			adjncy[index] = current->c->location;
			index++;
		}
		for(j=0;j<current->numparents;j++)
		{
			adjncy[index] = current->parents[j]->location;
			index++;
		}
	}
	xadj[N] = 2*(N-1);

	//Partition the system
	idx_t nverts = N;
	idx_t parts = np;
	idx_t ncon = 1;
	idx_t objval;
	idx_t* partitions = calloc(N,sizeof(idx_t));
	idx_t* vwgt = calloc(N,sizeof(idx_t));

	if(np != 1)
	{
		FILE* rainfile;

		//Check what type of rainfall data we are dealing with
		if(GlobalVars->rain_flag == 2)	//Binary
		{
			//Calculate the number of rainfall changes that occur
			unsigned int totalfiles = GlobalVars->last_file - GlobalVars->first_file + 1;
			float rainfall_buffer;
			unsigned int holder;
			float* last = malloc(N*sizeof(float));
			for(i=0;i<totalfiles;i++)
			{
				sprintf(filename,"%srain%i",GlobalVars->rain_filename,GlobalVars->first_file+i);
				rainfile = fopen(filename,"r");
				if(!rainfile)
				{
					printf("Error opening file %s.\n",filename);
					return;
				}

				//!!!! This loop should be done in parallel !!!!
				for(j=0;j<N;j++)
				{
					//Read in the storm data for this link
					fread(&rainfall_buffer,sizeof(float),1,rainfile);

					//This assumes different endianness
					holder = *(unsigned int*) &rainfall_buffer;
					holder = (((holder & 0x0000ffff)<<16) | ((holder & 0xffff0000)>>16));
					holder = (((holder & 0x00ff00ff)<<8) | ((holder & 0xff00ff00)>>8));
					rainfall_buffer = *(float*) &holder;

					if(i == 0 || !(.999 * rainfall_buffer <= last[j] && last[j] <= 1.001 * rainfall_buffer) )
					{
						last[j] = rainfall_buffer;
						vwgt[j]++;
					}
				}

				fclose(rainfile);
			}

			//Cleanup
			free(last);
		}
		else if(GlobalVars->rain_flag == 1)	//.str
		{
			unsigned int id,numtimes;
			rainfile = fopen(GlobalVars->rain_filename,"r");
			if(!rainfile)
			{
				printf("Error opening file %s.\n",GlobalVars->rain_filename);
				return;
			}

			fscanf(rainfile,"%*u");	//N
			for(i=0;i<N;i++)
			{
				fscanf(rainfile,"%u\n%u",&id,&numtimes);
				vwgt[i] = numtimes;
				for(j=0;j<numtimes;j++)
					fscanf(rainfile,"%*lf %*lf");
			}

			//Cleanup
			fclose(rainfile);
		}
		else	//Bad rain_flag
		{
			printf("Error: Invalid rain_flag in partitioning routine. Expected 1 or 2; got %hu.\n",GlobalVars->rain_flag);
			return NULL;
		}

		//Perform the partition
		retval = METIS_PartGraphKway(&nverts,&ncon,xadj,adjncy,vwgt,NULL,NULL,&parts,NULL,NULL, NULL,&objval,partitions);
		if(retval != METIS_OK)
		{
			printf("Error: METIS returned error code %i.\n",retval);
			return NULL;
		}
	}

	*my_N = 0;
	for(i=0;i<N;i++)
	{
		assignments[i] = partitions[i];	//!!!! Just use assignments? !!!!
		if(partitions[i] == my_rank)
		{
			(*my_sys)[*my_N] = i;
			(*my_N)++;
		}
	}

	//Set the getting array and determine number of sending and receiving links
	for(i=0;i<*my_N;i++)
	{
		//Receiving
		for(j=0;j<sys[(*my_sys)[i]]->numparents;j++)
		{
			loc = sys[(*my_sys)[i]]->parents[j]->location;
			if(assignments[loc] != my_rank)
			{
				getting[loc] = 1;
				my_data->receive_size[assignments[loc]]++;
			}
		}

		//Sending
		if(sys[(*my_sys)[i]]->c != NULL)
		{
			loc = sys[(*my_sys)[i]]->c->location;
			if(assignments[loc] != my_rank)
				my_data->send_size[assignments[loc]]++;
		}
	}

	//Reorder my_sys so that the links with lower numbering are towards the beginning
	merge_sort_distance(sys,*my_sys,*my_N);

	//Allocate space in my_data for recieving and sending
	for(j=0;j<np;j++)
	{
		my_data->receive_data[j] = (Link**) malloc(my_data->receive_size[j] * sizeof(Link*));
		my_data->send_data[j] = (Link**) malloc(my_data->send_size[j] * sizeof(Link*));
	}

	//Set the receive_data and send_data arrays
	int* current_receive_size = (int*) calloc(np,sizeof(int));
	int* current_send_size = (int*) calloc(np,sizeof(int));
	for(i=0;i<*my_N;i++)
	{
		//Receiving
		for(j=0;j<sys[(*my_sys)[i]]->numparents;j++)
		{
			loc = sys[(*my_sys)[i]]->parents[j]->location;
			if(assignments[loc] != my_rank)
			{
				my_data->receive_data[assignments[loc]][current_receive_size[assignments[loc]]] = sys[loc];
				current_receive_size[assignments[loc]]++;
			}
		}

		//Sending
		if(sys[(*my_sys)[i]]->c != NULL)
		{
			loc = sys[(*my_sys)[i]]->c->location;
			if(assignments[loc] != my_rank)
			{
				my_data->send_data[assignments[loc]][current_send_size[assignments[loc]]] = sys[(*my_sys)[i]];
				current_send_size[assignments[loc]]++;
			}
		}
	}

	//Clean up
	free(current_receive_size);
	free(current_send_size);
	free(xadj);
	free(adjncy);
	free(partitions);
	free(vwgt);
	(*my_sys) = (unsigned int*) realloc(my_sys,*my_N*sizeof(unsigned int));

	return assignments;
}


int* Partition_METIS_RainVolume(Link** sys,unsigned int N,Link** leaves,unsigned int numleaves,unsigned int** my_sys,unsigned int* my_N,TransData* my_data,short int *getting,UnivVars* GlobalVars)
{
	unsigned int i,j,start_index,end_index,extras,partition,loc,retval;
	unsigned int nodes_per_proc = numleaves/np;	//Number of leaves assigned to each process (except the last)
	char filename[256];
	Link* current;

	start_index = nodes_per_proc * my_rank;
	if(my_rank == np-1)		end_index = numleaves;
	else				end_index = nodes_per_proc * (my_rank + 1);
	*my_N = 0;
	unsigned int my_max_nodes = N;
	*my_sys = (unsigned int*) malloc(my_max_nodes * sizeof(unsigned int));	//The indices of this processes links (in sys)
	for(i=0;i<my_max_nodes;i++)	(*my_sys)[i] = -1;
	for(i=start_index;i<end_index;i++)	(*my_sys)[i-start_index] = leaves[i]->location;
	for(i=0;i<N;i++)	getting[i] = 0;

	//Initialize assignments
	int* assignments = (int*) malloc(N*sizeof(int));
	for(i=0;i<N;i++)	assignments[i] = -1;

	//Form the graph to partition
	idx_t* xadj = malloc((N+1)*sizeof(idx_t));
	idx_t* adjncy = malloc(2*(N-1)*sizeof(idx_t));
	idx_t index = 0;

	for(i=0;i<N;i++)
	{
		xadj[i] = index;
		current = sys[i];
		if(current->c != NULL)
		{
			adjncy[index] = current->c->location;
			index++;
		}
		for(j=0;j<current->numparents;j++)
		{
			adjncy[index] = current->parents[j]->location;
			index++;
		}
	}
	xadj[N] = 2*(N-1);

	//Partition the system
	idx_t nverts = N;
	idx_t parts = np;
	idx_t ncon = 1;
	idx_t objval;
	idx_t* partitions = calloc(N,sizeof(idx_t));
	idx_t* vwgt = calloc(N,sizeof(idx_t));
	float* sum = calloc(N,sizeof(float));

	if(np != 1)
	{
		FILE* rainfile;

		//Check what type of rainfall data we are dealing with
		if(GlobalVars->rain_flag == 2)	//Binary
		{
			//Calculate the number of rainfall changes that occur
			unsigned int totalfiles = GlobalVars->last_file - GlobalVars->first_file + 1;
			float rainfall_buffer;
			unsigned int holder;

			for(i=0;i<totalfiles;i++)
			{
				sprintf(filename,"%srain%i",GlobalVars->rain_filename,GlobalVars->first_file+i);
				rainfile = fopen(filename,"r");
				if(!rainfile)
				{
					printf("Error opening file %s.\n",filename);
					return;
				}

				//!!!! This loop should be done in parallel !!!!
				for(j=0;j<N;j++)
				{
					//Read in the storm data for this link
					fread(&rainfall_buffer,sizeof(float),1,rainfile);

					//This assumes different endianness
					holder = *(unsigned int*) &rainfall_buffer;
					holder = (((holder & 0x0000ffff)<<16) | ((holder & 0xffff0000)>>16));
					holder = (((holder & 0x00ff00ff)<<8) | ((holder & 0xff00ff00)>>8));
					rainfall_buffer = *(float*) &holder;

					sum[j] += rainfall_buffer;
				}

				fclose(rainfile);
			}

			//Compute the total volume (m^3)
			for(i=0;i<N;i++)
			{
				//Rainfall = mm/hr, A_h = km^2, file_time = mins
				sum[i] *= sys[i]->params->ve[GlobalVars->areah_idx] * GlobalVars->file_time * (1e3/60.0);
				vwgt[i] = (idx_t) ceil(sum[i] + 1.0);
			}
		}
		else if(GlobalVars->rain_flag == 1)	//.str
		{
			unsigned int id,numtimes;
			double prev_time,curr_time,prev_rain,curr_rain;
			rainfile = fopen(GlobalVars->rain_filename,"r");
			if(!rainfile)
			{
				printf("Error opening file %s.\n",GlobalVars->rain_filename);
				return;
			}

			fscanf(rainfile,"%*u");	//N

			for(i=0;i<N;i++)
			{
				fscanf(rainfile,"%u\n%u",&id,&numtimes);
				//vwgt[i] = numtimes;
				fscanf(rainfile,"%lf %lf",&prev_time,&prev_rain);
				for(j=1;j<numtimes;j++)
				{
					fscanf(rainfile,"%lf %lf",&curr_time,&curr_rain);
					sum[i] += prev_rain * (curr_time - prev_time) * sys[i]->params->ve[GlobalVars->areah_idx] * (1e3/60.0);
					prev_time = curr_time;
					prev_rain = curr_rain;
				}
				if(prev_time > GlobalVars->maxtime)
					sum[i] += prev_rain * (GlobalVars->maxtime - prev_time) * sys[i]->params->ve[GlobalVars->areah_idx] * (1e3/60.0) * 1e-6;
				vwgt[i] = (idx_t) ceil(sum[i] + 1.0);
			}

			//Cleanup
			fclose(rainfile);
		}
		else	//Bad rain_flag
		{
			printf("Error: Invalid rain_flag in partitioning routine. Expected 1 or 2; got %hu.\n",GlobalVars->rain_flag);
			return NULL;
		}
printf("%i: In\n",my_rank);
for(i=0;i<N;i++)
	printf("%u\n",vwgt[i]);
MPI_Barrier(MPI_COMM_WORLD);
		//Perform the partitioning
		retval = METIS_PartGraphKway(&nverts,&ncon,xadj,adjncy,vwgt,NULL,NULL,&parts,NULL,NULL, NULL,&objval,partitions);
printf("%i: Out\n",my_rank);
MPI_Barrier(MPI_COMM_WORLD);
		if(retval != METIS_OK)
		{
			printf("Error: METIS returned error code %i.\n",retval);
			return NULL;
		}
	}

	*my_N = 0;
	for(i=0;i<N;i++)
	{
		assignments[i] = partitions[i];	//!!!! Just use assignments? !!!!
		if(partitions[i] == my_rank)
		{
			(*my_sys)[*my_N] = i;
			(*my_N)++;
		}
	}
printf("%i: Here\n",my_rank);
MPI_Barrier(MPI_COMM_WORLD);
	//Set the getting array and determine number of sending and receiving links
	for(i=0;i<*my_N;i++)
	{
		//Receiving
		for(j=0;j<sys[(*my_sys)[i]]->numparents;j++)
		{
			loc = sys[(*my_sys)[i]]->parents[j]->location;
			if(assignments[loc] != my_rank)
			{
				getting[loc] = 1;
				my_data->receive_size[assignments[loc]]++;
			}
		}

		//Sending
		if(sys[(*my_sys)[i]]->c != NULL)
		{
			loc = sys[(*my_sys)[i]]->c->location;
			if(assignments[loc] != my_rank)
				my_data->send_size[assignments[loc]]++;
		}
	}
printf("%i: Here\n",my_rank);
MPI_Barrier(MPI_COMM_WORLD);
	//Reorder my_sys so that the links with lower numbering are towards the beginning
	merge_sort_distance(sys,*my_sys,*my_N);
printf("%i: Here\n",my_rank);
MPI_Barrier(MPI_COMM_WORLD);
	//Allocate space in my_data for recieving and sending
	for(j=0;j<np;j++)
	{
		my_data->receive_data[j] = (Link**) malloc(my_data->receive_size[j] * sizeof(Link*));
		my_data->send_data[j] = (Link**) malloc(my_data->send_size[j] * sizeof(Link*));
	}
printf("%i: Here\n",my_rank);
MPI_Barrier(MPI_COMM_WORLD);
	//Set the receive_data and send_data arrays
	int* current_receive_size = (int*) calloc(np,sizeof(int));
	int* current_send_size = (int*) calloc(np,sizeof(int));
	for(i=0;i<*my_N;i++)
	{
		//Receiving
		for(j=0;j<sys[(*my_sys)[i]]->numparents;j++)
		{
			loc = sys[(*my_sys)[i]]->parents[j]->location;
			if(assignments[loc] != my_rank)
			{
				my_data->receive_data[assignments[loc]][current_receive_size[assignments[loc]]] = sys[loc];
				current_receive_size[assignments[loc]]++;
			}
		}

		//Sending
		if(sys[(*my_sys)[i]]->c != NULL)
		{
			loc = sys[(*my_sys)[i]]->c->location;
			if(assignments[loc] != my_rank)
			{
				my_data->send_data[assignments[loc]][current_send_size[assignments[loc]]] = sys[(*my_sys)[i]];
				current_send_size[assignments[loc]]++;
			}
		}
	}
printf("%i: Here\n",my_rank);
MPI_Barrier(MPI_COMM_WORLD);

	//Clean up
	free(current_receive_size);
	free(current_send_size);
	free(xadj);
	free(adjncy);
	free(partitions);
	free(vwgt);
	free(sum);
	(*my_sys) = (unsigned int*) realloc(my_sys,*my_N*sizeof(unsigned int));

	return assignments;
}
*/
