#include "solvers.h"

void AsynchSolver(Link** sys,unsigned int N,unsigned int* my_sys,unsigned int my_N,UnivVars* GlobalVars,int* assignments,short int* getting,unsigned int* res_list,unsigned int res_size,unsigned int** id_to_loc,TempStorage* workspace,Forcing** forcings,ConnData** db_connections,TransData* my_data,short int print_flag,FILE* outputfile)
{
	unsigned int i,k;
	
	//Initialize remaining data
	short int* done = (short int*) malloc(my_N*sizeof(short int));
	short int parentsval;
	unsigned int alldone;
	RKSolutionNode *roottail;
	Link* current;
	unsigned int last_idx,curr_idx,around;
	unsigned int two_my_N = 2*my_N;
	int error_code;

	//Initialize values for forcing data
	unsigned int passes = 1;
	double maxtime;
	for(i=0;i<GlobalVars->num_forcings;i++)
	{
		if(forcings[i]->active)
		{
			forcings[i]->passes = forcings[i]->GetPasses(forcings[i],GlobalVars->maxtime,db_connections[ASYNCH_DB_LOC_FORCING_START+i]);
			passes = max(passes,forcings[i]->passes);
//printf("Before: %u, %u, %u\n",i,forcings[i]->passes,passes);
		}
	}

	//Start the main loop
	for(k=0;k<passes;k++)
	{
		alldone = 0;
		around = 0;
		current = sys[my_sys[my_N-1]];
		curr_idx = 0;
		last_idx = my_N - 1;
		for(i=0;i<my_N;i++)	done[i] = 0;

		//Read in next set of forcing data
		maxtime = GlobalVars->maxtime;
		for(i=0;i<GlobalVars->num_forcings;i++)
		{
			if(forcings[i]->active)
			{
//printf("Forcing %u is active  %e %e\n",i,sys[my_sys[0]]->last_t,forcings[i]->maxtime);
				if( fabs(sys[my_sys[0]]->last_t - forcings[i]->maxtime) < 1e-14 )
				{
					forcings[i]->maxtime = forcings[i]->GetNextForcing(sys,N,my_sys,my_N,assignments,GlobalVars,forcings[i],db_connections,id_to_loc,i);
					//(forcings[i]->iteration)++;	if flag is 3 (dbc), this happens in GetNextForcing
//printf("setting forcing maxtime to %e, iteration = %u\n",forcings[i]->maxtime,forcings[i]->iteration);
				}
				maxtime = min(maxtime,forcings[i]->maxtime);
			}
		}

		//If a state forcing is used, previously outputted data may need to be rewritten
		if(GlobalVars->res_flag)
		{
			for(i=0;i<my_N;i++)	//!!!! Can we loop over just the links with reservoirs? Improve id_to_loc. !!!!
			{
				current = sys[my_sys[i]];
				//if(current->res && fabs( (current->last_t) - (current->next_save - current->print_time) ) < 1e-12)
				if(current->res)
				{
					current->f(current->last_t,current->list->tail->y_approx,NULL,current->numparents,GlobalVars->global_params,current->forcing_values,current->qvs,current->params,current->state,current->user,current->list->tail->y_approx);
					if(current->save_flag && fabs(current->last_t - (current->next_save - current->print_time))/(current->last_t+1e-12) < 1e-6)
					{
						error_code = overwrite_last_step(current,GlobalVars,outputfile);
						if(error_code)
							printf("[%i]: Error overwritting last step at link %u: Got error code %u.\n",my_rank,current->ID,error_code);
					}
				}
			}
		}

		//Set a new step size
		//CalculateInitialStepSizes(sys,N,my_sys,my_N,assignments,getting,res_list,res_size,id_to_loc,GlobalVars,forcings,workspace,db_connections);
		Exchange_InitState_At_Forced(sys,N,assignments,getting,res_list,res_size,id_to_loc,GlobalVars);
		for(i=0;i<my_N;i++)
			sys[my_sys[i]]->h = InitialStepSize(sys[my_sys[i]]->last_t,sys[my_sys[i]],GlobalVars,workspace);

		//This might be needed. Sometimes some procs get stuck in Finish for communication, but makes runs much slower.
		MPI_Barrier(MPI_COMM_WORLD);

		if(sys[my_sys[0]]->last_t < GlobalVars->maxtime)
		while(alldone < my_N)
		{
			//Find the next link to iterate
			do
			{
				curr_idx = (curr_idx > 0) ? curr_idx - 1 : last_idx;
				around++;
			}while( (around < two_my_N) && (sys[my_sys[curr_idx]]->ready == 0 || done[curr_idx] == 1));
			current = sys[my_sys[curr_idx]];
/*
//if(current->h < 1e-10)
{
	printf("[%i]: ID = %u time = %e step = %e\n",my_rank,current->ID,current->last_t,current->h);
	Print_Vector(current->list->tail->y_approx);
//	sleep(1);
}
*/
			if(around >= two_my_N)
			{
				//Communicate with other processes
				Transfer_Data(my_data,sys,assignments,GlobalVars);
				around = 0;
				curr_idx = 0;
			}
			else	//Compute an iteration
			{
				//If the current link is not too far ahead, it can compute some iterations
				if(current->current_iterations < GlobalVars->iter_limit)
				{
					//Solve a few steps of the current link
					if(current->numparents == 0)	//Leaf
					{
						while(current->last_t + current->h < maxtime && current->current_iterations<GlobalVars->iter_limit)
						{
							for(i=0;i<GlobalVars->num_forcings;i++)		//!!!! Put this in RKSolver !!!!
								if(forcings[i]->active && current->last_t < current->forcing_change_times[i])	current->h = min(current->h,current->forcing_change_times[i] - current->last_t);
							current->rejected = current->RKSolver(current,GlobalVars,assignments,print_flag,outputfile,db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT],forcings,workspace);
						}

						if(current->last_t + current->h >= maxtime  && current->current_iterations < GlobalVars->iter_limit && current->last_t < maxtime)	//If less than a full step is needed, just finish up
						{
							for(i=0;i<GlobalVars->num_forcings;i++)
								if(forcings[i]->active && current->last_t < current->forcing_change_times[i])	current->h = min(current->h,current->forcing_change_times[i] - current->last_t);
							current->h = min(current->h,maxtime - current->last_t);
							current->rejected = current->RKSolver(current,GlobalVars,assignments,print_flag,outputfile,db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT],forcings,workspace);

							while(current->rejected == 0)
							{
								for(i=0;i<GlobalVars->num_forcings;i++)
									if(forcings[i]->active && current->last_t < current->forcing_change_times[i])	current->h = min(current->h,current->forcing_change_times[i] - current->last_t);
								current->h = min(current->h,maxtime - current->last_t);
								current->rejected = current->RKSolver(current,GlobalVars,assignments,print_flag,outputfile,db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT],forcings,workspace);
							}
						}
					}
					else	//Has parents
					{
						parentsval = 0;
						for(i=0;i<current->numparents;i++)
							parentsval += (current->last_t + current->h <= current->parents[i]->last_t);

						while(parentsval == current->numparents && current->current_iterations<GlobalVars->iter_limit)
						{
							for(i=0;i<GlobalVars->num_forcings;i++)
								if(forcings[i]->active && current->last_t < current->forcing_change_times[i])	current->h = min(current->h,current->forcing_change_times[i] - current->last_t);

							if(current->discont_count > 0 && current->h > current->discont[current->discont_start] - current->last_t)
							{
								current->h_old = current->h;
								current->h = current->discont[current->discont_start] - current->last_t;
							}

							current->rejected = current->RKSolver(current,GlobalVars,assignments,print_flag,outputfile,db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT],forcings,workspace);

							parentsval = 0;
							for(i=0;i<current->numparents;i++)
								parentsval += (current->last_t + current->h <= current->parents[i]->last_t);
						}

						parentsval = 0;
						for(i=0;i<current->numparents;i++)
							parentsval += (current->parents[i]->last_t >= maxtime);
						if(parentsval == current->numparents && current->current_iterations < GlobalVars->iter_limit && current->last_t < maxtime)		//If all parents are done, then current should finish up too
						{
							current->h = min(current->h,maxtime - current->last_t);
							for(i=0;i<GlobalVars->num_forcings;i++)
								if(forcings[i]->active && current->last_t < current->forcing_change_times[i])	current->h = min(current->h,current->forcing_change_times[i] - current->last_t);
							if(current->discont_count > 0 && current->h > current->discont[current->discont_start] - current->last_t)
							{
								current->h_old = current->h;
								current->h = current->discont[current->discont_start] - current->last_t;
							}

							current->rejected = current->RKSolver(current,GlobalVars,assignments,print_flag,outputfile,db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT],forcings,workspace);

							while(current->last_t < maxtime && current->current_iterations < GlobalVars->iter_limit)
							{
								for(i=0;i<GlobalVars->num_forcings;i++)
									if(forcings[i]->active && current->last_t < current->forcing_change_times[i])	current->h = min(current->h,current->forcing_change_times[i] - current->last_t);

								if(current->discont_count > 0 && current->h > current->discont[current->discont_start] - current->last_t)
								{
									current->h_old = current->h;
									current->h = current->discont[current->discont_start] - current->last_t;
								}

								current->h = min(current->h,maxtime - current->last_t);
								current->rejected = current->RKSolver(current,GlobalVars,assignments,print_flag,outputfile,db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT],forcings,workspace);
							}
						}

						if(current->current_iterations < GlobalVars->iter_limit)
							current->ready = 0;
					}

					//Notify the child that a parent made progress
					Link* child = current->c;
					if(child != NULL && assignments[child->location] == my_rank)
					{
						//Make sure the child can take a step if current has reached limit
						if(current->current_iterations >= GlobalVars->iter_limit)
							child->h = min(child->h,current->last_t - child->last_t);

						if(child->h + child->last_t > current->last_t)	child->h *= .999;

						double next_t = child->last_t + child->h;
						parentsval = 0;
						for(i=0;i<child->numparents;i++)
							parentsval += (child->parents[i]->last_t >= next_t)||(child->parents[i]->last_t >= maxtime);

						if(parentsval == child->numparents)
							child->ready = 1;
						else
							child->ready = 0;
					}

					//See if current has parents that hit their limit
					for(i=0;i<current->numparents;i++)
					{
						if(current->parents[i]->current_iterations >= GlobalVars->iter_limit)
							current->h = min(current->h,current->parents[i]->last_t - current->last_t);

						if(current->h + current->last_t > current->parents[i]->last_t)	current->h *= .999;
					}

					parentsval = 0;
					for(i=0;i<current->numparents;i++)
						parentsval += (current->last_t + current->h <= current->parents[i]->last_t);
					if(parentsval == current->numparents && current->current_iterations < GlobalVars->iter_limit)
						current->ready = 1;

					//Check if current is done
					if(current->last_t >= maxtime)
					{
						alldone++;
						done[curr_idx] = 1;
						current->last_t = maxtime;	//In case of roundoff errors
					}

					//Reduce last_idx, if possible
					while(done[last_idx] == 1 && last_idx > 0)
						last_idx--;

					//If current is a root link, trash its data
					if(current->c == NULL)
					{
						roottail = current->list->tail;
						while(current->list->head != roottail)
						{
							Remove_Head_Node(current->list);
							(current->current_iterations)--;
						}
					}
				}
			}
		}//endwhile
		else
		{
//			if(my_rank == 0)	printf("%i: Should be done, k is %i/%i\n",my_rank,k,passes-1);
			break;
		}

		Transfer_Data_Finish(my_data,sys,assignments,GlobalVars);

		//Ensure all data is received !!!! This is sloppy. Transfer_Data_Finish should handle this. !!!!
		MPI_Barrier(MPI_COMM_WORLD);
		Transfer_Data_Finish(my_data,sys,assignments,GlobalVars);

		//if((rain_flag == 2 || rain_flag == 3) && my_rank == 0)
//		if(my_rank == 0)
//			printf("%i: Going to next set of forcing data, k is %i/%i\n",my_rank,k,passes-1);
	}

	//Cleanup
	free(done);
}




