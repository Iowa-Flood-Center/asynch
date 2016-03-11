#include "comm.h"

// **********  MPI related routines  **********


//Tranfers data amongst processes. Uses asynchronous communication scheme.
//TransData* my_data: Contains information about how the processes will communicate.
//Link** sys: The river system.
//UnivVars* GlobalVars: Contains all the information that is shared by every link in the system.
//Note: a few "silly" initializations take place before the loops. This prevents Valgrind from complaining on some systems about pointless errors.
void Transfer_Data(TransData* my_data,Link** sys,int* assignments,UnivVars* GlobalVars)
{
	int i,j,m,n,sender,steps_to_transfer = 0,curr_idx = 0,parval,s,dim,flag,position,total_links,count,removed = 0,order = 0,num_times = 0;
	unsigned int loc = 0,l,num_dense;
	double discont_time = 0.0;
	RKSolutionNode* node;
	Link *current,*next,*prev;
	MPI_Status status;

	//If sending
	for(i=0;i<np;i++)
	{
		if(my_data->send_size[i] != 0 || my_data->receive_size[i] != 0)
		{
			if(my_data->sent_flag[i])	MPI_Test(my_data->send_requests[i],&flag,MPI_STATUS_IGNORE);
			if(!my_data->sent_flag[i] || flag)
			{
				position = 0;

				//Pack data
				total_links = 0;
				for(l=0;l<my_data->send_size[i];l++)
				{
					current = my_data->send_data[i][l];
					curr_idx = current->location;

					//Figure out how many steps will be sent.
					steps_to_transfer = 0;
					node = current->list->head->next;

					while(steps_to_transfer < current->current_iterations-1 && steps_to_transfer < GlobalVars->max_transfer_steps && steps_to_transfer + current->steps_on_diff_proc < GlobalVars->iter_limit)
					{
						steps_to_transfer++;
						node = node->next;
					}

					//Pack all the data for each step and each discontinuity time
					if(steps_to_transfer > 0 || current->discont_send_count > 0)
					{
						total_links++;
						node = current->list->head->next;
						s = current->list->s;
						dim = current->dim;
						num_dense = current->num_dense;
						current->steps_on_diff_proc += steps_to_transfer;

						//Pack the steps
						MPI_Pack(&(current->location),1,MPI_INT,my_data->send_buffer[i],my_data->send_buffer_size[i],&position,MPI_COMM_WORLD);
						MPI_Pack(&steps_to_transfer,1,MPI_INT,my_data->send_buffer[i],my_data->send_buffer_size[i],&position,MPI_COMM_WORLD);
						for(m=0;m<steps_to_transfer;m++)
						{
							MPI_Pack(&(node->t),1,MPI_DOUBLE,my_data->send_buffer[i],my_data->send_buffer_size[i],&position,MPI_COMM_WORLD);
							MPI_Pack(node->y_approx->ve,dim,MPI_DOUBLE,my_data->send_buffer[i],my_data->send_buffer_size[i],&position,MPI_COMM_WORLD);
							for(n=0;n<s;n++)
								MPI_Pack(node->k[n]->ve,num_dense,MPI_DOUBLE,my_data->send_buffer[i],my_data->send_buffer_size[i],&position,MPI_COMM_WORLD);
							MPI_Pack(&(node->state),1,MPI_INT,my_data->send_buffer[i],my_data->send_buffer_size[i],&position,MPI_COMM_WORLD);

							Remove_Head_Node(current->list);
							node = node->next;
							(current->current_iterations)--;
						}

						//Pack the discontinuity times
						MPI_Pack(&(current->discont_send_count),1,MPI_INT,my_data->send_buffer[i],my_data->send_buffer_size[i],&position,MPI_COMM_WORLD);
						for(m=0;(unsigned int)m<current->discont_send_count;m++)
						{
							MPI_Pack(&(current->discont_send[m]),1,MPI_DOUBLE,my_data->send_buffer[i],my_data->send_buffer_size[i],&position,MPI_COMM_WORLD);
							MPI_Pack(&(current->discont_order_send[m]),1,MPI_INT,my_data->send_buffer[i],my_data->send_buffer_size[i],&position,MPI_COMM_WORLD);
						}

						current->discont_send_count = 0;
					}
				}

				//Pack iterations
				for(l=0;l<my_data->receive_size[i];l++)
				{
					current = my_data->receive_data[i][l];
					if(current->iters_removed > 0)
					{
						MPI_Pack(&(current->location),1,MPI_UNSIGNED,my_data->send_buffer[i],my_data->send_buffer_size[i],&position,MPI_COMM_WORLD);
						MPI_Pack(&(current->iters_removed),1,MPI_INT,my_data->send_buffer[i],my_data->send_buffer_size[i],&position,MPI_COMM_WORLD);
						current->iters_removed = 0;
					}
				}

				//If there's information to send, send it!
				if(position != 0)
				{
					my_data->sent_flag[i] = 1;
					(my_data->num_sent[i])++;
					MPI_Isend(my_data->send_buffer[i],position,MPI_PACKED,i,total_links,MPI_COMM_WORLD,my_data->send_requests[i]);
				}
			} //End if(flag)
		}
	} //End loop over processes (i)

	//Check if any incoming messages have been completely received
	for(i=0;i<np;i++)
	{
		if(my_data->receiving_flag[i])
		{
			MPI_Test(my_data->receive_requests[i],&flag,&status);
			if(flag)
			{
				sender = status.MPI_SOURCE;
				total_links = status.MPI_TAG;
				MPI_Get_count(&status,MPI_PACKED,&count);
				position = 0;

				//Unpack data
				for(j=0;j<total_links;j++)
				{
					MPI_Unpack(my_data->receive_buffer[i],count,&position,&curr_idx,1,MPI_INT,MPI_COMM_WORLD);

					//Unpack the steps
					MPI_Unpack(my_data->receive_buffer[i],count,&position,&steps_to_transfer,1,MPI_INT,MPI_COMM_WORLD);
					current = sys[curr_idx];
					s = current->list->s;
					dim = current->dim;
					num_dense = current->num_dense;

					for(m=0;m<steps_to_transfer;m++)
					{
						node = New_Step(current->list);
						node->t = 0.0; node->state = 0;
						MPI_Unpack(my_data->receive_buffer[i],count,&position,&(node->t),1,MPI_DOUBLE,MPI_COMM_WORLD);
						MPI_Unpack(my_data->receive_buffer[i],count,&position,node->y_approx->ve,dim,MPI_DOUBLE,MPI_COMM_WORLD);
						for(n=0;n<s;n++)
							MPI_Unpack(my_data->receive_buffer[i],count,&position,node->k[n]->ve,num_dense,MPI_DOUBLE,MPI_COMM_WORLD);
						MPI_Unpack(my_data->receive_buffer[i],count,&position,&(node->state),1,MPI_INT,MPI_COMM_WORLD);
					}

					if(steps_to_transfer > 0)
					{
						//Put the steps in place
						current->current_iterations += steps_to_transfer;
						current->last_t = node->t;
						parval = 0;
						for(n=0;n<current->c->numparents;n++)
							parval += (current->c->last_t < current->c->parents[n]->last_t);
						if(parval == current->c->numparents)
							current->c->ready = 1;

						//Make sure the child can take a step if current has reached limit
						if(current->current_iterations >= GlobalVars->iter_limit)
							current->c->h = min(current->c->h,current->last_t - current->c->last_t);
						if(current->c->h + current->c->last_t > current->last_t)	current->c->h *= .999;
					}

					//Unpack the discontinuity times
					MPI_Unpack(my_data->receive_buffer[i],count,&position,&num_times,1,MPI_INT,MPI_COMM_WORLD);
					for(m=0;m<num_times;m++)
					{
						MPI_Unpack(my_data->receive_buffer[i],count,&position,&discont_time,1,MPI_DOUBLE,MPI_COMM_WORLD);
						MPI_Unpack(my_data->receive_buffer[i],count,&position,&order,1,MPI_INT,MPI_COMM_WORLD);
						prev = current;
						next = current->c;
						for(n=order;(unsigned int)n<GlobalVars->max_localorder && next != NULL;n++)
						{
							if(my_rank == assignments[next->location] && n < next->method->localorder)
							{
								next->discont_end = Insert_Discontinuity(discont_time,next->discont_start,next->discont_end,&(next->discont_count),GlobalVars->discont_size,next->discont,next->ID);
							}
							else if(my_rank != assignments[next->location])
							{
								Insert_SendDiscontinuity(discont_time,n,&(prev->discont_send_count),GlobalVars->discont_size,prev->discont_send,prev->discont_order_send,prev->ID);
								break;
							}

							prev = next;
							next = next->c;
						}
					}
				}

				//Unpack iterations
				while(position < count)
				{
					MPI_Unpack(my_data->receive_buffer[i],count,&position,&loc,1,MPI_UNSIGNED,MPI_COMM_WORLD);
					MPI_Unpack(my_data->receive_buffer[i],count,&position,&removed,1,MPI_INT,MPI_COMM_WORLD);
					sys[loc]->steps_on_diff_proc -= removed;
				}

				my_data->receiving_flag[i] = 0;
				(my_data->num_recv[i])++;
			}
		}
	}

	//Begin receiving messages
	for(i=0;i<np;i++)
	{
		if(!(my_data->receiving_flag[i]))
		{
			MPI_Iprobe(i,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status);
			if(flag)
			{
				total_links = status.MPI_TAG;
				MPI_Get_count(&status,MPI_PACKED,&count);
				MPI_Irecv(my_data->receive_buffer[i],count,MPI_PACKED,i,total_links,MPI_COMM_WORLD,my_data->receive_requests[i]);
				my_data->receiving_flag[i] = 1;
			}
		}
	}
}

//Tranfers data amongst processes. Use for asynchronous communication scheme.
//Use for asynchronous communication and only after this process has finished all calculations. Sends all remaining data.
//TransData* my_data: Contains information about how the processes will communicate.
//Link** sys: The river system.
//UnivVars* GlobalVars: Contains all the information that is shared by every link in the system.
void Transfer_Data_Finish(TransData* my_data,Link** sys,int* assignments,UnivVars* GlobalVars)
{
	int i,j,m,n,steps_to_transfer = 0,curr_idx = 0,s,dim,flag,position,total_links,sender,count,removed = 0,order = 0,num_times = 0;
	unsigned int loc = 0,parval,l,num_dense;
	unsigned int data_to_send = 0;
	unsigned int data_sent = 0;
	double discont_time = 0.0;
	RKSolutionNode* node;
	Link *current,*next,*prev;
	MPI_Status status;

	//Check how much data still must be sent
	//Note: There should never be discontinuities to send AND no steps for a given link
	//Note2: I am pretty sure that there can never be a situation where a link on proc i has no steps
	//		to send but has iteration information to send to a proc j, unless i NEVER sends
	//		steps to j. The if(send_size[i] != 0 || receive_size[i] != 0) should take care of
	//		that situation. There could, however, be an issue if i had no steps to send to ANY
	//		process, but has iteration info to send to j (where j does not receive steps from i).
	//		So I have added iteration_removed to the data_to_send total.
	for(i=0;i<np;i++)
	{
		for(l=0;l<my_data->send_size[i];l++)
			data_to_send += my_data->send_data[i][l]->current_iterations - 1;
		for(l=0;l<my_data->receive_size[i];l++)		//See Note2 above
			if(my_data->receive_data[i][l]->iters_removed > 0)	data_to_send++;
	}

	//Send all remaining data
	while(data_sent < data_to_send)
	{
		for(i=0;i<np;i++)
		{
			if(my_data->send_size[i] != 0 || my_data->receive_size[i] != 0)
			{
				if(my_data->sent_flag[i])	MPI_Test(my_data->send_requests[i],&flag,MPI_STATUS_IGNORE);
				if(!my_data->sent_flag[i] || flag)
				{
					position = 0;

					//Pack data
					total_links = 0;
					for(l=0;l<my_data->send_size[i];l++)
					{
						current = my_data->send_data[i][l];
						curr_idx = current->location;

						//Figure out how many steps will be sent.
						steps_to_transfer = 0;
						node = current->list->head->next;

						while(steps_to_transfer < current->current_iterations-1 && steps_to_transfer < GlobalVars->max_transfer_steps && steps_to_transfer + current->steps_on_diff_proc < GlobalVars->iter_limit)
						{
							steps_to_transfer++;
							node = node->next;
						}

						//Pack all the data for each step.
						if(steps_to_transfer > 0 || current->discont_send_count > 0)
						{
							data_sent += steps_to_transfer;
							total_links++;
							node = current->list->head->next;
							s = current->list->s;
							dim = current->dim;
							num_dense = current->num_dense;
							current->steps_on_diff_proc += steps_to_transfer;

							//Pack the steps
							MPI_Pack(&(current->location),1,MPI_INT,my_data->send_buffer[i],my_data->send_buffer_size[i],&position,MPI_COMM_WORLD);
							MPI_Pack(&steps_to_transfer,1,MPI_INT,my_data->send_buffer[i],my_data->send_buffer_size[i],&position,MPI_COMM_WORLD);
							for(m=0;m<steps_to_transfer;m++)
							{
								MPI_Pack(&(node->t),1,MPI_DOUBLE,my_data->send_buffer[i],my_data->send_buffer_size[i], &position,MPI_COMM_WORLD);
								MPI_Pack(node->y_approx->ve,dim,MPI_DOUBLE,my_data->send_buffer[i],my_data->send_buffer_size[i],&position,MPI_COMM_WORLD);
								for(n=0;n<s;n++)
									MPI_Pack(node->k[n]->ve,num_dense,MPI_DOUBLE,my_data->send_buffer[i],my_data->send_buffer_size[i],&position,MPI_COMM_WORLD);
									//MPI_Pack(node->k[n]->ve,dim,MPI_DOUBLE,my_data->send_buffer[i],my_data->send_buffer_size[i],&position,MPI_COMM_WORLD);
								MPI_Pack(&(node->state),1,MPI_INT,my_data->send_buffer[i],my_data->send_buffer_size[i],&position,MPI_COMM_WORLD);

								Remove_Head_Node(current->list);
								node = node->next;
								(current->current_iterations)--;
							}

							//Pack the discontinuity times
							MPI_Pack(&(current->discont_send_count),1,MPI_INT,my_data->send_buffer[i],my_data->send_buffer_size[i],&position,MPI_COMM_WORLD);
							for(m=0;(unsigned int)m<current->discont_send_count;m++)
							{
								MPI_Pack(&(current->discont_send[m]),1,MPI_DOUBLE,my_data->send_buffer[i],my_data->send_buffer_size[i],&position,MPI_COMM_WORLD);
								MPI_Pack(&(current->discont_order_send[m]),1,MPI_INT,my_data->send_buffer[i],my_data->send_buffer_size[i],&position,MPI_COMM_WORLD);
							}

							current->discont_send_count = 0;
						}
					}

					//Pack iterations
					for(l=0;l<my_data->receive_size[i];l++)
					{
						current = my_data->receive_data[i][l];
						if(current->iters_removed > 0)
						{
							data_sent++;	//See Note2 above
							MPI_Pack(&(current->location),1,MPI_UNSIGNED,my_data->send_buffer[i],my_data->send_buffer_size[i], &position,MPI_COMM_WORLD);
							MPI_Pack(&(current->iters_removed),1,MPI_INT,my_data->send_buffer[i],my_data->send_buffer_size[i],&position,MPI_COMM_WORLD);
							current->iters_removed = 0;
						}
					}

					//If there's information to send, send it!
					if(position != 0)
					{
						my_data->sent_flag[i] = 1;
						(my_data->num_sent[i])++;
						MPI_Isend(my_data->send_buffer[i],position,MPI_PACKED,i,total_links,MPI_COMM_WORLD,my_data->send_requests[i]);
					}
				} //End if(flag)
			}
		} //End loop over processes (i)

		//Check if any incoming messages have been completely received
		for(i=0;i<np;i++)
		{
			if(my_data->receiving_flag[i])
			{
				MPI_Test(my_data->receive_requests[i],&flag,&status);
				if(flag)
				{
					sender = status.MPI_SOURCE;
					total_links = status.MPI_TAG;
					MPI_Get_count(&status,MPI_PACKED,&count);
					position = 0;

					//Unpack data
					for(j=0;j<total_links;j++)
					{
						MPI_Unpack(my_data->receive_buffer[i],count,&position,&curr_idx,1,MPI_INT,MPI_COMM_WORLD);
						MPI_Unpack(my_data->receive_buffer[i],count,&position,&steps_to_transfer,1,MPI_INT,MPI_COMM_WORLD);
						current = sys[curr_idx];
						s = current->list->s;
						dim = current->dim;
						num_dense = current->num_dense;

						for(m=0;m<steps_to_transfer;m++)
						{
							node = New_Step(current->list);
							node->t = 0.0; node->state = 0;
							MPI_Unpack(my_data->receive_buffer[i],count,&position,&(node->t),1,MPI_DOUBLE,MPI_COMM_WORLD);
							MPI_Unpack(my_data->receive_buffer[i],count,&position,node->y_approx->ve,dim,MPI_DOUBLE,MPI_COMM_WORLD);
							for(n=0;n<s;n++)
								MPI_Unpack(my_data->receive_buffer[i],count,&position,node->k[n]->ve,num_dense,MPI_DOUBLE,MPI_COMM_WORLD);
								//MPI_Unpack(my_data->receive_buffer[i],count,&position,node->k[n]->ve,dim,MPI_DOUBLE,MPI_COMM_WORLD);
							MPI_Unpack(my_data->receive_buffer[i],count,&position,&(node->state),1,MPI_INT,MPI_COMM_WORLD);
						}

						if(steps_to_transfer > 0)
						{
							//Put the steps in place
							current->current_iterations += steps_to_transfer;
							current->last_t = node->t;
							parval = 0;
							for(n=0;n<current->c->numparents;n++)
								parval += (current->c->last_t < current->c->parents[n]->last_t);
							if(parval == current->c->numparents)
								current->c->ready = 1;

							//Make sure the child can take a step if current has reached limit
							if(current->current_iterations >= GlobalVars->iter_limit)
								current->c->h = min(current->c->h,current->last_t - current->c->last_t);
							if(current->c->h + current->c->last_t > current->last_t)	current->c->h *= .999;
						}

						//Unpack the discontinuity times
						MPI_Unpack(my_data->receive_buffer[i],count,&position,&num_times,1,MPI_INT,MPI_COMM_WORLD);
						for(m=0;m<num_times;m++)
						{
							MPI_Unpack(my_data->receive_buffer[i],count,&position,&discont_time,1,MPI_DOUBLE,MPI_COMM_WORLD);
							MPI_Unpack(my_data->receive_buffer[i],count,&position,&order,1,MPI_INT,MPI_COMM_WORLD);

							prev = current;
							next = current->c;
							for(n=order;(unsigned int)n<GlobalVars->max_localorder && next != NULL;n++)
							{
								if(my_rank == assignments[next->location] && n < next->method->localorder)
									next->discont_end = Insert_Discontinuity(discont_time,next->discont_start,next->discont_end,&(next->discont_count),GlobalVars->discont_size,next->discont,next->ID);
								else if(my_rank != assignments[next->location])
								{
									Insert_SendDiscontinuity(discont_time,n,&(prev->discont_send_count),GlobalVars->discont_size,prev->discont_send,prev->discont_order_send,prev->ID);
									break;
								}

								prev = next;
								next = next->c;
							}
						}
					}

					//Unpack iterations (need this for rainfall so the last step will get sent)
					while(position < count)
					{
						MPI_Unpack(my_data->receive_buffer[i],count,&position,&loc,1,MPI_UNSIGNED,MPI_COMM_WORLD);
						MPI_Unpack(my_data->receive_buffer[i],count,&position,&removed,1,MPI_INT,MPI_COMM_WORLD);
						sys[loc]->steps_on_diff_proc -= removed;
					}

					my_data->receiving_flag[i] = 0;
					(my_data->num_recv[i])++;
				}
			}
		}

		//Begin receiving messages
		for(i=0;i<np;i++)
		{
			if(!my_data->receiving_flag[i])
			{
				MPI_Iprobe(i,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status);
				if(flag)
				{
					total_links = status.MPI_TAG;
					MPI_Get_count(&status,MPI_PACKED,&count);
					MPI_Irecv(my_data->receive_buffer[i],count,MPI_PACKED,i,total_links,MPI_COMM_WORLD,my_data->receive_requests[i]);
					my_data->receiving_flag[i] = 1;
				}
			}
		}
	} //End while
}


void Exchange_InitState_At_Forced(Link** system,unsigned int N,unsigned int* assignments,short int* getting,unsigned int* res_list,unsigned int res_size,unsigned int** id_to_loc,UnivVars* GlobalVars)
{
	unsigned int j,loc;

	//Find links with state forcing
	if(GlobalVars->res_flag)
	{
		//Setup links with forcing
		for(j=0;j<res_size;j++)
		{
			loc = find_link_by_idtoloc(res_list[j],id_to_loc,N);
			if(loc < N && assignments[loc] == my_rank)
			{
				//!!!! Not sure if this the way to go... !!!!
				system[loc]->f(system[loc]->last_t,system[loc]->list->tail->y_approx,NULL,0,GlobalVars->global_params,system[loc]->forcing_values,system[loc]->qvs,system[loc]->params,system[loc]->state,system[loc]->user,system[loc]->list->tail->y_approx);

				//Check if this initial state needs to be sent to other procs
				if(system[loc]->c && assignments[system[loc]->c->location] != my_rank)
					MPI_Send(system[loc]->list->tail->y_approx->ve,system[loc]->dim,MPI_DOUBLE,assignments[system[loc]->c->location],0,MPI_COMM_WORLD);
			}
			else if(getting[loc])
				MPI_Recv(system[loc]->list->tail->y_approx->ve,system[loc]->dim,MPI_DOUBLE,assignments[loc],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
	}
}


//Allocate space for a transmitting scheme
//Returns a pointer to a newly allocated TransData
TransData* Initialize_TransData()
{
	int i;
	TransData* data = (TransData*) malloc(sizeof(TransData));
	data->send_data = (Link***) malloc(np*sizeof(Link**));
	data->receive_data = (Link***) malloc(np*sizeof(Link**));
	data->send_size = (unsigned int*) calloc(np,sizeof(int));
	data->receive_size = (unsigned int*) calloc(np,sizeof(int));

	data->send_buffer = (char**) malloc(np*sizeof(char*));
	data->receive_buffer = (char**) malloc(np*sizeof(char*));
	data->send_requests = (MPI_Request**) malloc(np*sizeof(MPI_Request*));
	data->receive_requests = (MPI_Request**) malloc(np*sizeof(MPI_Request*));
	for(i=0;i<np;i++)
	{
		data->send_requests[i] = (MPI_Request*) malloc(sizeof(MPI_Request));
		data->receive_requests[i] = (MPI_Request*) malloc(sizeof(MPI_Request));
	}
	data->sent_flag = (short int*) calloc(np,sizeof(short int));
	data->receiving_flag = (short int*) calloc(np,sizeof(short int));
	data->send_buffer_size = (unsigned int*) malloc(np*sizeof(unsigned int));
	data->receive_buffer_size = (unsigned int*) malloc(np*sizeof(unsigned int));
	data->num_sent = (unsigned int*) calloc(np,sizeof(unsigned int));
	data->num_recv = (unsigned int*) calloc(np,sizeof(unsigned int));
	data->totals = (unsigned int*) malloc(np*sizeof(unsigned int));

	return data;
}


//Makes sure all sends have completed and receives any outstanding data that has been sent to
//this process. Resets the sent and receiving flags.
//TransData* data: The data to be cleaned
void Flush_TransData(TransData* data)
{
	int i,j;
	MPI_Status status;
	MPI_Barrier(MPI_COMM_WORLD);

	//Receive all outstanding messages where a recv has been posted
	for(i=0;i<np;i++)
	{
		if(data->receiving_flag[i])
		{
			MPI_Wait(data->receive_requests[i],MPI_STATUS_IGNORE);
			data->receiving_flag[i] = 0;
			(data->num_recv[i])++;
		}
	}

/*
	//Make sure all sends are complete
	for(i=0;i<np;i++)
	{
		if(data->sent_flag[i])
		{
			MPI_Wait(data->send_requests[i],MPI_STATUS_IGNORE);
			data->sent_flag[i] = 0;
		}
	}
*/

	//See how many messages remain to be received
	for(i=0;i<np;i++)
		MPI_Scatter(data->num_sent,1,MPI_INT,&(data->totals[i]),1,MPI_INT,i,MPI_COMM_WORLD);

	//Make sure all sends are complete
	for(i=0;i<np;i++)
	{
		if(data->sent_flag[i])
		{
			MPI_Wait(data->send_requests[i],MPI_STATUS_IGNORE);
			data->sent_flag[i] = 0;
		}
	}

	//Receive any remaining messages
	for(i=0;i<np;i++)
	{
		for(j=data->num_recv[i];(unsigned int)j<data->totals[i];j++)
			MPI_Recv(data->receive_buffer[i],data->receive_buffer_size[i],MPI_PACKED,i,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
	}

	for(i=0;i<np;i++)
	{
		data->num_recv[i] = 0;
		data->num_sent[i] = 0;
	}

	MPI_Barrier(MPI_COMM_WORLD);
}

//Free space for a transmitting scheme
//TransData* data: The data to be freed
void TransData_Free(TransData* data)
{
	int i;

	//Clean out any left over messages.
	MPI_Finalized(&i);
	if(!i)	Flush_TransData(data);

	//Free memory
	for(i=0;i<np;i++)	free(data->send_data[i]);
	free(data->send_data);
	for(i=0;i<np;i++)	free(data->receive_data[i]);
	free(data->receive_data);

	for(i=0;i<np;i++)
	{
		free(data->send_requests[i]);
		free(data->receive_requests[i]);
		if(data->send_buffer[i] != NULL)	free(data->send_buffer[i]);
		if(data->receive_buffer[i] != NULL)	free(data->receive_buffer[i]);
	}
	free(data->send_requests);
	free(data->receive_requests);
	free(data->send_buffer);
	free(data->receive_buffer);
	free(data->sent_flag);
	free(data->receiving_flag);
	free(data->send_buffer_size);
	free(data->receive_buffer_size);
	free(data->send_size);
	free(data->receive_size);
	free(data->num_sent);
	free(data->num_recv);
	free(data->totals);
	free(data);
}


// **********  Postgresql related routines  **********

//Create a ConnData object
ConnData* CreateConnData(char* connectstring)
{
	ConnData* conninfo = (ConnData*) malloc(sizeof(ConnData));
	conninfo->conn = NULL;
	conninfo->query = (char*) malloc(1024*sizeof(char));
	conninfo->connectinfo = (char*) malloc(1024*sizeof(char));
	strcpy(conninfo->connectinfo,connectstring);
	conninfo->time_offset = 0;
	conninfo->num_queries = 0;
	conninfo->queries = NULL;
	return conninfo;
}

//Destroy a ConnData object
void ConnData_Free(ConnData* conninfo)
{
	int i;
	if(conninfo)
	{
		//if(my_rank == 0 && conninfo->conn != NULL)	PQfinish(conninfo->conn);
		if(conninfo->conn && PQstatus(conninfo->conn) == CONNECTION_OK)
			PQfinish(conninfo->conn);
		free(conninfo->query);
		for(i=0;i<conninfo->num_queries;i++)	free(conninfo->queries[i]);
		if(conninfo->queries)	free(conninfo->queries);
		free(conninfo->connectinfo);
		free(conninfo);
	}
}

//Switch the database conninfo connects to
void SwitchDB(ConnData* conninfo,char connectinfo[])
{
	if(conninfo->conn && PQstatus(conninfo->conn) == CONNECTION_OK)
		PQfinish(conninfo->conn);
	conninfo->conn = NULL;
	sprintf(conninfo->connectinfo,"%s",connectinfo);
}


//Connect to the database with information stored in connectinfo
int ConnectPGDB(ConnData* conninfo)
{
	conninfo->conn = PQconnectdb(conninfo->connectinfo);
	if(PQstatus(conninfo->conn) == CONNECTION_BAD)
	{
		printf("[%i]: Error: Unable to connect to the database.\n",my_rank);
		return 1;
	}
	PQsetNoticeProcessor(conninfo->conn,ShutUp,NULL);	//Disable annoying notices
	return 0;
}

//Disconnect from the database
void DisconnectPGDB(ConnData* conninfo)
{
	if(PQstatus(conninfo->conn) == CONNECTION_OK)
	{
		PQfinish(conninfo->conn);
		conninfo->conn = NULL;
	}
}

//Check if an error related to an sql query occurred.
int CheckResError(PGresult* res,char* event)
{
	short int status = PQresultStatus(res);
	if( !(status == PGRES_COMMAND_OK || status == PGRES_TUPLES_OK) )
	{
		printf("[%i]: SQL error encountered while %s. %hi\n",my_rank,event,status);
		printf("[%i]: %s\n",my_rank,PQresultErrorMessage(res));

		return 1;
	}

	return 0;
}

//Check if a query returned a certain value
int CheckResState(PGresult* res,short int error_code)
{
	short int status = PQresultStatus(res);
	if(status == error_code)	return 0;
	else
	{
		printf("[%i]: Error: did not get error code %hi. Got %hi.\n",my_rank,error_code,status);
		return 1;
	}
}

//Check if connection to SQL database is still good
void CheckConnConnection(ConnData* conninfo)
{
	if(PQstatus(conninfo->conn) == CONNECTION_BAD)
	{
		printf("[%i]: Connection to database lost. Attempting to reconnect...\n",my_rank);
		PQreset(conninfo->conn);
		printf("[%i]: Connection reestablished.\n",my_rank);
	}
}

//Disables postgresql notices
void ShutUp(void *arg, const char *message)
{
	return;
}


