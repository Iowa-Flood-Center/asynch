#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <assert.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>

#include <blas.h>
#include <system.h>

//Frees link.
//Link* link: link to be freed.
//unsigned int list_length: the length of the list stored with link.
//int rkd_flag: should be 1 if an .rkd file was used, 0 if not.
void Destroy_Link(Link* link, int rkd_flag, Forcing* forcings, GlobalVars* global)
{
    unsigned int i;
    assert(link != NULL);

    free(link->params);

    if (link->my != NULL)
    {
        if (link->my->forcing_values)
            free(link->my->forcing_values);
        if (link->my->forcing_indices)
            free(link->my->forcing_indices);
        if (link->my->forcing_change_times)
            free(link->my->forcing_change_times);
        if (rkd_flag)
            Destroy_ErrorData(link->my->error_data);
        Destroy_List(&link->my->list);
        
        free(link->peak_value);
        if (link->discont != NULL)
            free(link->discont);
        if (link->discont_send != NULL)
        {
            free(link->discont_send);
            free(link->discont_order_send);
        }

        if (link->my->forcing_data)
        {
            for (i = 0; i < global->num_forcings; i++)
            {
                if (forcings[i].flag != 4 && forcings[i].flag != 7)
                    Destroy_ForcingData(&(link->my->forcing_data[i]));
            }
            free(link->my->forcing_data);
        }
        //if (link->qvs != NULL)
        //{
        //    /*
        //                for(i=0;i<link->qvs->n_values;i++)
        //                    free(link->qvs->points[i]);
        //    */
        //    free(link->qvs->points);
        //    free(link->qvs->points_array);
        //    free(link->qvs);
        //}
#if defined (ASYNCH_HAVE_IMPLICIT_SOLVER)
        m_free(&link->JMatrix);
        m_free(&link->CoefMat);
        v_free(&link->sol_diff);
        if (link->Z_i != NULL)
        {
            for (i = 0; i < global->max_s; i++)
                v_free(&link->Z_i[i]);
            free(link->Z_i);
        }
#endif
        /*
                if(global->template_flag && link->equations != NULL)
                {
                    //mupRelease(link->equations->parser);
                    printf("Warning: Not freeing parser in system.c, Destroy_Link\n");
                    v_free(&link->equations->variable_values);
                    free(link->equations);
                }
        */
    }

    if (link->dense_indices)
        free(link->dense_indices);

    free(link->parents);
}

//Frees rain
//ForcingData** forcing_buff: forcing data to be freed
void Destroy_ForcingData(TimeSerie* forcing_buff)
{
    if (forcing_buff)
    {
        if (forcing_buff->data)
            free(forcing_buff->data);
    }
}


/// Creates a list to hold the data for an ODE.
///
/// \param t0: the initial time.
/// \param y0: Vector of the initial data [num_dof]
/// \param num_dof: the number of degree of freedom of the ODE.
/// \param num_dense_dof
/// \param num_stages: the number of stages in the RKMethod.
/// \param list_length: the maximum number of steps to store in the list.
void Init_List(RKSolutionList* list, double t0, double *y0, unsigned int num_dof, unsigned int num_dense_dof, unsigned short int num_stages, unsigned int list_length)
{
    assert(list_length > 0);
    if (list_length < 2)
        printf("Warning in Create_List: list_length is %u.\n", list_length);

    memset(list, 0, sizeof(RKSolutionList));
    list->nodes = (RKSolutionNode*)calloc(list_length, sizeof(RKSolutionNode));

    //Set the next and prev ptrs for each node
    list->nodes[0].next = &list->nodes[1];
    list->nodes[0].prev = &list->nodes[list_length - 1];

    for (unsigned int i = 1; i < list_length - 1; i++)
    {
        list->nodes[i].next = &list->nodes[i + 1];
        list->nodes[i].prev = &list->nodes[i - 1];
    }

    list->nodes[list_length - 1].next = &list->nodes[0];
    list->nodes[list_length - 1].prev = &list->nodes[list_length - 2];

    //Allocate space for all the vectors
    list->y_storage = malloc(list_length * num_dof * sizeof(double));
    list->k_storage = malloc(list_length * num_stages * num_dense_dof * sizeof(double));

    for (unsigned int i = 0; i < list_length; i++)
    {
        list->nodes[i].y_approx = list->y_storage + i * num_dof;
        list->nodes[i].k = list->k_storage + i * num_stages * num_dense_dof;
    }

    //Set remaining fields
    list->head = &list->nodes[0];
    list->tail = &list->nodes[0];
    list->num_stages = num_stages;

    //Store the initial step
    list->head->t = t0;
    dcopy(y0, list->head->y_approx, 0, num_dof);
}

//Frees the data list.
void Destroy_List(RKSolutionList* list)
{
    free(list->nodes);
    free(list->y_storage);
    free(list->k_storage);
}

//Removes the first node in list.
void Remove_Head_Node(RKSolutionList* list)
{
    if (list->head == list->tail)	//Keep the tail at head if the list has no data
        list->tail = list->tail->next;

    list->head = list->head->next;
}

//Adds a new step to list, after tail. The new node is returned.
RKSolutionNode* New_Step(RKSolutionList* list)
{
    list->tail = list->tail->next;
    return list->tail;
}

//This method undoes the last step taken.
//Used for rejecting steps.
void Undo_Step(RKSolutionList* list)
{
    list->tail = list->tail->prev;
}

//Frees an RKMethod
void Destroy_RKMethod(RKMethod* method)
{
    assert(method != NULL);

    //if (method->b)
    //    free(method->b);
    //if (method->b_theta)
    //    free(method->b_theta);
    //if (method->b_theta_deriv)
    //    free(method->b_theta_deriv);
}

//Frees an ErrorData
void Destroy_ErrorData(ErrorData* error)
{
    assert(error != NULL);
    free(error->abstol);
    free(error->reltol);
    free(error->abstol_dense);
    free(error->reltol_dense);
}

//Allocates workspace for RK solvers
void Create_Workspace(Workspace *workspace, unsigned int max_dim, unsigned short num_stages, unsigned short max_parents)
{
    memset(workspace, 0, sizeof(Workspace));

    workspace->sum =malloc(max_dim * sizeof(double));
    workspace->temp = malloc(max_dim * sizeof(double));
    workspace->temp2 = malloc(max_dim * sizeof(double));
    workspace->temp3 = malloc(max_dim * sizeof(double));


    workspace->parents_approx = malloc(max_parents * max_dim * sizeof(double));
    workspace->stages_parents_approx = malloc(num_stages * max_parents * max_dim * sizeof(double));

    //workspace->temp_k = (VEC*)malloc(num_stages * sizeof(VEC));
    //for (unsigned int i = 0; i < num_stages; i++)
    //    workspace->temp_k[i] = v_init(dim);

    workspace->temp_k = malloc(num_stages * max_dim * sizeof(double));

    for (unsigned int i = 0; i < num_stages; i++)
        workspace->temp_k_slices[i] = workspace->temp_k + i * max_dim;

#if defined(ASYNCH_HAVE_IMPLICIT_SOLVER)
    workspace->ipiv = (int*)malloc(s*dim * sizeof(int));
    workspace->rhs = v_init(s*dim);
    //workspace->CoefMat = m_get(s*dim,s*dim);
    workspace->JMatrix = m_get(dim, dim);
    workspace->Z_i = (VEC*)malloc(s * sizeof(VEC));
    for (i = 0; i < s; i++)	workspace->Z_i[i] = v_init(dim);
    workspace->err = v_init(dim);
#endif // defined(ASYNCH_HAVE_IMPLICIT_SOLVER)
}

//Deallocates workspace for RK solvers
void Destroy_Workspace(Workspace* workspace, unsigned short int num_stages, unsigned short int max_parents)
{
    free(workspace->sum);
    free(workspace->temp);    
    free(workspace->temp2);
    free(workspace->temp3);

    //for (unsigned int i = 0; i < num_stages; i++)
    //{
    //    for (unsigned int j = 0; j < max_parents; j++)
    //        v_free(&workspace->temp_parent_approx[i][j]);
    //    free(workspace->temp_parent_approx[i]);
    //}
    free(workspace->parents_approx);
    free(workspace->stages_parents_approx);

    //for (unsigned int i = 0; i < num_stages; i++)
    //    v_free(&workspace->temp_k[i]);
    free(workspace->temp_k);
    
#if defined(ASYNCH_HAVE_IMPLICIT_SOLVER)
    free(workspace->ipiv);
    v_free(&workspace->rhs);
    //m_free(&workspace->CoefMat);
    m_free(&workspace->JMatrix);
    for (i = 0; i < s; i++)
        v_free(&workspace->Z_i[i]);
    free(workspace->Z_i);
    v_free(&workspace->err);
#endif // defined(ASYNCH_HAVE_IMPLICIT_SOLVER)
}


void Destroy_Outputs(Output* outputs, unsigned int num_outputs)
{
    for (unsigned int i = 0; i < num_outputs; i++)
        if (outputs[i].name)
            free(outputs[i].name);

    free(outputs);
}

//Deallocates UnivVars
void Destroy_UnivVars(GlobalVars* global)
{
    free(global->peakflow_function_name);

    Destroy_Outputs(global->outputs, global->num_outputs);

    if (global->rsv_filename)
        free(global->rsv_filename);
    if (global->rvr_filename)
        free(global->rvr_filename);
    if (global->prm_filename)
        free(global->prm_filename);
    if (global->init_filename)
        free(global->init_filename);
    if (global->rain_filename)
        free(global->rain_filename);
    if (global->dam_filename)
        free(global->dam_filename);
    if (global->hydrosave_filename)
        free(global->hydrosave_filename);
    if (global->peaksave_filename)
        free(global->peaksave_filename);
    if (global->temp_filename)
        free(global->temp_filename);
    if (global->peakfilename)
        free(global->peakfilename);
    if (global->hydros_loc_filename)
        free(global->hydros_loc_filename);
    if (global->peaks_loc_filename)
        free(global->peaks_loc_filename);
    if (global->hydro_table)
        free(global->hydro_table);
    if (global->peak_table)
        free(global->peak_table);
    if (global->dump_table)
        free(global->dump_table);
    if (global->dump_loc_filename)
        free(global->dump_loc_filename);
    //if (global->global_params)
    //    free(&global->global_params);
    if (global->print_indices)
        free(global->print_indices);
    free(global);
}

//Inserts a time into the list of discontinuities
unsigned int Insert_Discontinuity(double time, unsigned int start, unsigned int end, unsigned int* count, unsigned int size, double* array, unsigned int id)
{
    if (*count == 0)
    {
        array[start] = time;
        (*count)++;
        return (end + 1) % size;
    }
    else if ((end + 1) % size == start)
    {
        //printf("Warning: A discontinuity is being ignored. Increase the size of the discontinuity buffer. id = %u time = %f\n",id,time);
        return end;
    }

    unsigned int curr = end;
    unsigned int toofar = (start == 0) ? size - 1 : start - 1;

    while (curr != toofar && time < array[curr])
        curr = (curr == 0) ? size - 1 : curr - 1;
    if (curr != toofar && time == array[curr])	return end;

    curr = end;
    while (curr != toofar && time < array[curr])
    {
        array[(curr + 1) % size] = array[curr];
        curr = (curr == 0) ? size - 1 : curr - 1;
    }

    curr = (curr + 1) % size;
    end = (end + 1) % size;
    array[curr] = time;
    (*count)++;
    /*
        if(end == start)
        {
            printf("Warning: A discontinuity is being ignored. Increase the size of the discontinuity buffer.\n");
        }
    */
    return end;
}


//Inserts a time into the list of discontinuities to be send to another processes
void Insert_SendDiscontinuity(double time, unsigned int order, unsigned int* count, unsigned int size, double* array, unsigned int* order_array, unsigned int id)
{
    if (*count == 0)
    {
        array[0] = time;
        order_array[0] = order;
        *count = 1;
        return;
    }
    else if (*count + 1 == size)
    {
        //printf("Warning: A discontinuity is being ignored in send. Increase the size of the discontinuity buffer. id = %i time = %f\n",id,time);
        return;
    }

    unsigned int curr = *count - 1;

    while (curr < size && time < array[curr])	curr--;
    if (curr < size && time == array[curr])
    {
        if (order_array[curr] < order)	order_array[curr] = order;
        return;
    }

    curr = *count - 1;
    while (curr < size && time < array[curr])
    {
        array[curr + 1] = array[curr];
        order_array[curr + 1] = order_array[curr];
        curr--;
    }

    curr++;
    array[curr] = time;
    order_array[curr] = order;
    (*count)++;

    /*
        if(*count == size)
        {
            printf("Warning: A discontinuity is being ignored. Increase the size of the discontinuity buffer.\n");
        }
    */

}

