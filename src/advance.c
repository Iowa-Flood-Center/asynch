#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <assert.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>

#include <minmax.h>
#include <processdata.h>
#include <rksteppers.h>
#include <structs.h>


void Advance(
    Link *sys, unsigned int N,
    Link **my_sys, unsigned int my_N,
    GlobalVars* globals,
    int* assignments, short int* getting, unsigned int* res_list, unsigned int res_size, const Lookup * const id_to_loc,
    Workspace* workspace,
    Forcing* forcings,
    ConnData* db_connections,
    TransData* my_data,
    int print_level,
    FILE* outputfile)
{
    //Initialize remaining data
    short int* done = (short int*)malloc(my_N * sizeof(short int));
    short int parentsval;
    Link* current;
    unsigned int last_idx, curr_idx, around;
    unsigned int two_my_N = 2 * my_N;
    int error_code;
	bool print_flag = false;

	if (print_level >= 1)
		print_flag = true;
	
    //Initialize values for forcing data
	if ((print_level >= 2) && (my_rank == 0))
		printf("[%i] Initializing forcings values...", my_rank);
    for (unsigned int i = 0; i < globals->num_forcings; i++)
    {
        if (forcings[i].active)
        {
            forcings[i].passes = forcings[i].GetPasses(&forcings[i], globals->maxtime, &db_connections[ASYNCH_DB_LOC_FORCING_START + i]);
            //passes = max(passes, forcings[i].passes);
            //printf("Before: %u, %u, %u\n",i,forcings[i].passes,passes);
        }
    }
	if ((print_level >= 2) && (my_rank == 0))
		printf("done.\n");

    ////Snapshot passes
    //if (globals->dump_loc_flag == 4)
    //    passes = max(passes, (unsigned int)ceil(globals->maxtime / globals->dump_time));

    //Start the main loop
    while (globals->t < globals->maxtime)
    {
        around = 0;
        current = my_sys[my_N - 1];
        curr_idx = 0;
        last_idx = my_N - 1;

        //Advance the current time
        globals->t = my_sys[0]->last_t;
		if ((print_level >= 2) && (my_rank == 0))
			printf("[%i] Solver at %0.2f of %0.2f minutes.\n", my_rank, globals->t, globals->maxtime);
        
        memset(done, 0, my_N * sizeof(short int));

        //Read in next set of forcing data
		if ((print_level >= 2) && (my_rank == 0))
		{
			printf("[%i] Reading next set of forcing data...", my_rank);
			fflush(stdout);
		}
        double maxtime = globals->maxtime;
        for (unsigned int i = 0; i < globals->num_forcings; i++)
        {
            if (forcings[i].active)
            {
                //printf("Forcing %u is active  %e %e\n",i,sys[my_sys[0]]->last_t,forcings[i].maxtime);
                if ((fabs(globals->t - forcings[i].maxtime) < 1e-14)  && (forcings[i].iteration < forcings[i].passes))
                {
                    forcings[i].maxtime = forcings[i].GetNextForcing(sys, N, my_sys, my_N, assignments, globals, &forcings[i], db_connections, id_to_loc, i);
                    //(forcings[i].iteration)++;	if flag is 3 (dbc), this happens in GetNextForcing
                    //printf("setting forcing maxtime to %e, iteration = %u\n",forcings[i].maxtime,forcings[i].iteration);
                }
                maxtime = min(maxtime, forcings[i].maxtime);
            }
        }
		if ((print_level >= 2) && (my_rank == 0))
		{
			printf("done.\n[%i] Check shapshot next time...", my_rank);
			fflush(stdout);
		}
		//Check shapshot next time
        if (globals->dump_loc_flag == 4)
        {
            double next_time = fmod(globals->t, globals->dump_time);
            if (next_time < 1e-14)
            {
                globals->output_func.CreateSnapShot(sys, N, assignments, globals, NULL, NULL);
                next_time = globals->t + globals->dump_time;
            }
            else
			{
                next_time = ceil(globals->t / globals->dump_time) * globals->dump_time;
			}

            maxtime = min(maxtime, next_time);
        }
		if ((print_level >= 2) && (my_rank == 0))
			printf("done.\n");

        //If a state forcing is used, previously outputted data may need to be rewritten
        if (globals->res_flag)
        {
            if ((print_level >= 2) && (my_rank == 0))
            {
                printf("[%i] Considering reservoir timeseries...", my_rank);
                fflush(stdout);
            }
            for (unsigned int i = 0; i < my_N; i++)	//!!!! Can we loop over just the links with reservoirs? Improve id_to_loc. !!!!
            {
                current = my_sys[i];
                assert(current->my != NULL);

                //if(current->res && fabs( (current->last_t) - (current->next_save - current->print_time) ) < 1e-12)                
                if (current->has_res)
                {
                    current->differential(current->last_t, current->my->list.tail->y_approx, current->dim, NULL, 0, 0, globals->global_params, current->params, current->my->forcing_values, current->qvs, current->state, current->user, current->my->list.tail->y_approx);
                    if (current->save_flag && fabs(current->last_t - (current->next_save - current->print_time)) / (current->last_t + 1e-12) < 1e-6)
                    {
                        error_code = overwrite_last_step(current, globals, outputfile);
                        if (error_code)
                            printf("[%i]: Error overwritting last step at link %u: Got error code %u.\n", my_rank, current->ID, error_code);
                    }
                }
            }
            if ((print_level >= 2) && (my_rank == 0))
                printf("done.\n");
        }
        
        // Update forcing
		if ((print_level >= 2) && (my_rank == 0))
		{
			printf("[%i] Updating forcings...", my_rank);
			fflush(stdout);
		}
        Exchange_InitState_At_Forced(sys, N, assignments, getting, res_list, res_size, id_to_loc, globals);
		
		if ((print_level >= 2) && (my_rank == 0))
		{
			printf("done.\n[%i] Solving...", my_rank);
			fflush(stdout);
		}
        
        for (unsigned int i = 0; i < my_N; i++)
        {
            my_sys[i]->h = InitialStepSize(my_sys[i]->last_t, my_sys[i], globals, workspace);
            assert(my_sys[i]->h > 1e-12);
        }

        //This might be needed. Sometimes some procs get stuck in Finish for communication, but makes runs much slower.
        MPI_Barrier(MPI_COMM_WORLD);
        if (globals->t < globals->maxtime)
        {
            unsigned int alldone = 0;
            while (alldone < my_N)
            {
                //Find the next link to iterate
                do
                {
                    curr_idx = (curr_idx > 0) ? curr_idx - 1 : last_idx;
                    around++;
                } while ((around < two_my_N) && (my_sys[curr_idx]->ready == 0 || done[curr_idx] == 1));
                current = my_sys[curr_idx];

                if (around >= two_my_N)
                {
                    //Communicate with other processes
                    Transfer_Data(my_data, sys, assignments, globals);
                    around = 0;
                    curr_idx = 0;
                }
                else	//Compute an iteration
                {
                    //If the current link is not too far ahead, it can compute some iterations
                    if (current->current_iterations < globals->iter_limit)
                    {
						// adlz
						// if ((print_level >= 2) && (my_rank == 0))
						//	printf("       ...computting iteration %i (max %i).\n", current->current_iterations, globals->iter_limit);
                        //Solve a few steps of the current link
                        if (current->num_parents == 0)	//Leaf
                        {
							// adlz
							// if ((print_level >= 2) && (my_rank == 0))
							//	printf("       ...solving a leaf...\n");
                            while (current->last_t + current->h < maxtime && current->current_iterations < globals->iter_limit)
                            {
								// adlz
								// if ((print_level >= 2) && (my_rank == 0))
								// {
								//	printf("       ...iterating 1 (%i of %i)...\n", current->current_iterations, globals->iter_limit);
								//	printf("         as %f + %f (%f) < %f.\n", current->last_t, current->h, current->last_t + current->h, maxtime);
								// }
                                for (unsigned int i = 0; i < globals->num_forcings; i++)		//!!!! Put this in solver !!!!
                                    if (forcings[i].active && current->last_t < current->my->forcing_change_times[i])
                                        current->h = min(current->h, current->my->forcing_change_times[i] - current->last_t);
                                current->rejected = current->solver(current, globals, assignments, print_flag, outputfile, &db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT], forcings, workspace);
                            }

							// adlz
							//if ((print_level >= 2) && (my_rank == 0))
							//	printf("       ...iterating 1 done.\n");

                            if (current->last_t + current->h >= maxtime  && current->current_iterations < globals->iter_limit && current->last_t < maxtime)	//If less than a full step is needed, just finish up
                            {
                                for (unsigned int i = 0; i < globals->num_forcings; i++)
                                    if (forcings[i].active && current->last_t < current->my->forcing_change_times[i])
                                        current->h = min(current->h, current->my->forcing_change_times[i] - current->last_t);
                                current->h = min(current->h, maxtime - current->last_t);
                                assert(current->h > 0);
                                current->rejected = current->solver(current, globals, assignments, print_flag, outputfile, &db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT], forcings, workspace);

                                while (current->rejected == 0)
                                {
									for (unsigned int i = 0; i < globals->num_forcings; i++)
                                        if (forcings[i].active && current->last_t < current->my->forcing_change_times[i])
                                            current->h = min(current->h, current->my->forcing_change_times[i] - current->last_t);
                                    current->h = min(current->h, maxtime - current->last_t);
                                    assert(current->h > 0);
                                    current->rejected = current->solver(current, globals, assignments, print_flag, outputfile, &db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT], forcings, workspace);
                                }
                            }
                        }
                        else	//Has parents
                        {
							parentsval = 0;
                            for (unsigned int i = 0; i < current->num_parents; i++)
                                parentsval += (current->last_t + current->h <= current->parents[i]->last_t);

                            while (parentsval == current->num_parents && current->current_iterations < globals->iter_limit)
                            {
                                for (unsigned int i = 0; i < globals->num_forcings; i++)
                                    if (forcings[i].active && current->last_t < current->my->forcing_change_times[i])
                                        current->h = min(current->h, current->my->forcing_change_times[i] - current->last_t);

                                if (current->discont_count > 0 && current->h > current->discont[current->discont_start] - current->last_t)
                                {
#if defined (ASYNCH_HAVE_IMPLICIT_SOLVER)
                                    current->h_old = current->h;
#endif
                                    current->h = current->discont[current->discont_start] - current->last_t;
                                    assert(current->h > 0);
                                }

                                current->rejected = current->solver(current, globals, assignments, print_flag, outputfile, &db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT], forcings, workspace);

                                parentsval = 0;
                                for (unsigned int i = 0; i < current->num_parents; i++)
                                    parentsval += (current->last_t + current->h <= current->parents[i]->last_t);
                            }

                            parentsval = 0;
                            for (unsigned int i = 0; i < current->num_parents; i++)
                                parentsval += (current->parents[i]->last_t >= maxtime);

                            if (parentsval == current->num_parents && current->current_iterations < globals->iter_limit && current->last_t < maxtime)		//If all parents are done, then current should finish up too
                            {
                                current->h = min(current->h, maxtime - current->last_t);
                                assert(current->h > 0);
                                for (unsigned int i = 0; i < globals->num_forcings; i++)
                                    if (forcings[i].active && current->last_t < current->my->forcing_change_times[i])
                                        current->h = min(current->h, current->my->forcing_change_times[i] - current->last_t);
                                if (current->discont_count > 0 && current->h > current->discont[current->discont_start] - current->last_t)
                                {
#if defined (ASYNCH_HAVE_IMPLICIT_SOLVER)
                                    current->h_old = current->h;
#endif
                                    current->h = current->discont[current->discont_start] - current->last_t;
                                    assert(current->h > 0);
                                }

                                current->rejected = current->solver(current, globals, assignments, print_flag, outputfile, &db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT], forcings, workspace);
                                assert(current->h > 0);

                                while (current->last_t < maxtime && current->current_iterations < globals->iter_limit)
                                {
                                    for (unsigned int i = 0; i < globals->num_forcings; i++)
                                        if (forcings[i].active && current->last_t < current->my->forcing_change_times[i])
                                            current->h = min(current->h, current->my->forcing_change_times[i] - current->last_t);

                                    if (current->discont_count > 0 && current->h > current->discont[current->discont_start] - current->last_t)
                                    {
#if defined (ASYNCH_HAVE_IMPLICIT_SOLVER)
                                        current->h_old = current->h;
#endif
                                        current->h = current->discont[current->discont_start] - current->last_t;

                                    }

                                    current->h = min(current->h, maxtime - current->last_t);
                                    current->rejected = current->solver(current, globals, assignments, print_flag, outputfile, &db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT], forcings, workspace);
                                }
                            }

                            if (current->current_iterations < globals->iter_limit)
                                current->ready = 0;
                        }

                        //Notify the child that a parent made progress
                        Link* child = current->child;
                        if (child != NULL && assignments[child->location] == my_rank)
                        {
                            assert(child->h > 0);

                            //Make sure the child can take a step if current has reached limit
                            if ((current->current_iterations >= globals->iter_limit) && (current->last_t > child->last_t))
                                child->h = min(child->h, current->last_t - child->last_t);

                            // TODO improve on this
                            if(child->h + child->last_t > current->last_t)
                                child->h *= .999;

                            assert(child->h > 0.0);

                            double next_t = child->last_t + child->h;
                            parentsval = 0;
                            for (unsigned int i = 0; i < child->num_parents; i++)
                                parentsval += (child->parents[i]->last_t >= next_t) || (child->parents[i]->last_t >= maxtime);

                            if (parentsval == child->num_parents)
                                child->ready = 1;
                            else
                                child->ready = 0;
                        }

                        //See if current has parents that hit their limit
                        for (unsigned int i = 0; i < current->num_parents; i++)
                        {
                            if (current->parents[i]->current_iterations >= globals->iter_limit)
                                current->h = min(current->h, current->parents[i]->last_t - current->last_t);

                            // TODO improve on this
                            if(current->h + current->last_t > current->parents[i]->last_t)
                                current->h *= .999;
                        }

                        parentsval = 0;
                        for (unsigned int i = 0; i < current->num_parents; i++)
                            parentsval += (current->last_t + current->h <= current->parents[i]->last_t);
                        if (parentsval == current->num_parents && current->current_iterations < globals->iter_limit)
                            current->ready = 1;

                        //Check if current is done
                        if (current->last_t >= maxtime)
                        {
                            alldone++;
                            done[curr_idx] = 1;
                            current->last_t = maxtime;	//In case of roundoff errors
                        }

                        //Reduce last_idx, if possible
                        while (done[last_idx] == 1 && last_idx > 0)
                            last_idx--;

						// adlz
						//if ((print_level >= 2) && (my_rank == 0))
						//	printf("       about to trash...\n");

                        //If current is a root link, trash its data
                        if (current->child == NULL)
                        {
                            RKSolutionNode *roottail = current->my->list.tail;
                            while (current->my->list.head != roottail)
                            {
                                Remove_Head_Node(&current->my->list);
                                (current->current_iterations)--;
                            }
                        }
                    }
                }

                assert(current->h > 0);
            }//endwhile
        }

		if ((print_level >= 2) && (my_rank == 0))
		{
			printf("done.\n[%i] Synchronizing...", my_rank);
			fflush(stdout);
		}

        Transfer_Data_Finish(my_data, sys, assignments, globals);

        //Ensure all data is received !!!! This is sloppy. Transfer_Data_Finish should handle this. !!!!
        MPI_Barrier(MPI_COMM_WORLD);
        Transfer_Data_Finish(my_data, sys, assignments, globals);

        //if((rain_flag == 2 || rain_flag == 3) && my_rank == 0)
//		if(my_rank == 0)
//			printf("%i: Going to next set of forcing data, k is %i/%i\n",my_rank,k,passes-1);

		if ((print_level >= 2) && (my_rank == 0))
		{
			printf("done.\n[%i] * * * * * * * * * * * * * * * *\n", my_rank);
			fflush(stdout);
		}
    }

    if (my_rank == 0)
        printf("\n");

    //Cleanup
    free(done);
}
