#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <assert.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>

#include "system.h"

//Frees link_i.
//Link* link_i: link to be freed.
//unsigned int list_length: the length of the list stored with link_i.
//int rkd_flag: should be 1 if an .rkd file was used, 0 if not.
void Destroy_Link(Link* link_i,unsigned int list_length,int rkd_flag,Forcing** forcings,UnivVars* GlobalVars)
{
	unsigned int i;
	if(link_i == NULL)	return;

	v_free(&link_i->params);

	if(link_i->method != NULL)
	{
		free(link_i->forcing_values);
		free(link_i->forcing_indices);
		free(link_i->forcing_change_times);
		if(rkd_flag)	Destroy_ErrorData(link_i->errorinfo);
		Destroy_List(link_i->list,list_length);
		//v_free(link_i->params);
		v_free(&link_i->peak_value);
		if(link_i->discont != NULL)	free(link_i->discont);
		if(link_i->discont_send != NULL)
		{
			free(link_i->discont_send);
			free(link_i->discont_order_send);
		}
		for(i=0;i<GlobalVars->num_forcings;i++)
		{
			if(link_i->forcing_buff && forcings[i]->flag != 4 && forcings[i]->flag != 7 && link_i->forcing_buff[i] != NULL)
				Destroy_ForcingData(&(link_i->forcing_buff[i]));
		}
		free(link_i->forcing_buff);
		if(link_i->qvs != NULL)
		{
/*
			for(i=0;i<link_i->qvs->n_values;i++)
				free(link_i->qvs->points[i]);
*/
			free(link_i->qvs->points);
			free(link_i->qvs->points_array);
			free(link_i->qvs);
		}
		m_free(&link_i->JMatrix);
		m_free(&link_i->CoefMat);
		v_free(&link_i->sol_diff);
        if (link_i->Z_i != NULL)
        {
		    for(i=0;i<GlobalVars->max_s;i++)
			    v_free(&link_i->Z_i[i]);
		    free(link_i->Z_i);
        }
/*
		if(GlobalVars->template_flag && link_i->equations != NULL)
		{
			//mupRelease(link_i->equations->parser);
			printf("Warning: Not freeing parser in system.c, Destroy_Link\n");
			v_free(&link_i->equations->variable_values);
			free(link_i->equations);
		}
*/
	}

	if(link_i->dense_indices)	free(link_i->dense_indices);

	free(link_i->parents);
	free(link_i);
}

//Frees rain
//ForcingData** forcing_buff: forcing data to be freed
void Destroy_ForcingData(ForcingData** forcing_buff)
{
	unsigned int i;
	if(forcing_buff && *forcing_buff)
	{
		for(i=0;i<(*forcing_buff)->n_times;i++)
			free((*forcing_buff)->rainfall[i]);
		free((*forcing_buff)->rainfall);
		free(*forcing_buff);
	}
}


//Creates a list to hold the data for an ODE.
//VEC* y0: the initial data.
//double t0: the initial time.
//int dim: the dimension of the ODE.
//int s: the number of steps in the RKMethod.
//unsigned int list_length: the maximum number of steps to store in the list.
//Returns the created list.
RKSolutionList* Create_List(VEC y0,double t0,int dim,unsigned int num_dense,unsigned short int s,unsigned int list_length)
{
	unsigned int i,j;

	if(list_length < 2)	printf("Warning in Create_List: list_length is %u.\n",list_length);

	//Allocate space
	RKSolutionList* list = (RKSolutionList*) malloc(sizeof(RKSolutionList));
    memset(list, 0, sizeof(RKSolutionList));
	list->list_data = (RKSolutionNode*) malloc(list_length * sizeof(RKSolutionNode));
    memset(list->list_data, 0, sizeof(list_length * sizeof(RKSolutionNode)));
	//for(i=0;i<list_length;i++)
    //  list->list_data[i] = (RKSolutionNode*) malloc(sizeof(RKSolutionNode));

	//Set the next and prev ptrs for each node
	list->list_data[0].next = &list->list_data[1];
	list->list_data[0].prev = &list->list_data[list_length-1];

	for(i=1;i<list_length-1;i++)
	{
		list->list_data[i].next = &list->list_data[i+1];
		list->list_data[i].prev = &list->list_data[i-1];
	}

	list->list_data[list_length-1].next = &list->list_data[0];
	list->list_data[list_length-1].prev = &list->list_data[list_length-2];

	//Allocate space for all the vectors
	for(i=0;i<list_length;i++)
	{
		list->list_data[i].y_approx = v_get(dim);
		list->list_data[i].k = (VEC*) malloc(s * sizeof(VEC));
		for(j=0;j<s;j++)
            list->list_data[i].k[j] = v_get(num_dense);
	}

	//Set remaining fields
	list->head = &list->list_data[0];
	list->tail = &list->list_data[0];
	list->s = s;

	//Store the initial step
	list->head->t = t0;
	v_copy(y0,list->head->y_approx);

	return list;
}

//Frees the data list.
void Destroy_List(RKSolutionList* list,unsigned int list_length)
{
	unsigned int i;
	for(i=0;i<list_length;i++)
	{
		Destroy_Node(&list->list_data[i],list->s);
	}
	free(list->list_data);
	free(list);
}

//Removes the first node in list.
void Remove_Head_Node(RKSolutionList* list)
{
	if(list->head == list->tail)	//Keep the tail at head if the list has no data
		list->tail = list->tail->next;

	list->head = list->head->next;
}

//Frees a node from a data list.
//RKSolutionNode* node: node to be freed.
//int s: size of node->k (number of steps).
void Destroy_Node(RKSolutionNode *node,unsigned short int s)
{
	unsigned int i;
	if(node->k != NULL)
	{
		for(i=0;i<s;i++)
			v_free(&node->k[i]);
		free(node->k);
	}

	v_free(&node->y_approx);
	node->next = NULL;
	node->prev = NULL;
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
	m_free(&method->A);
	v_free(&method->b);
	v_free(&method->b_theta);
	v_free(&method->b_theta_deriv);
	v_free(&method->c);
	v_free(&method->e);
	v_free(&method->d);
	v_free(&method->w);

	free(method);
}

//Frees an ErrorData
void Destroy_ErrorData(ErrorData* error)
{
    assert(error != NULL);
	v_free(&error->abstol);
	v_free(&error->reltol);
	v_free(&error->abstol_dense);
	v_free(&error->reltol_dense);
	free(error);
}

//Allocates workspace for RK solvers
TempStorage* Create_Workspace(unsigned int dim,unsigned short int s,unsigned short int max_parents)
{
	unsigned int i,j;

	TempStorage* workspace = (TempStorage*) malloc(sizeof(TempStorage));
    memset(workspace, 0, sizeof(TempStorage));

	workspace->temp = v_get(dim);
	workspace->sum = v_get(dim);
	workspace->temp2 = v_get(dim);
	workspace->temp3 = v_get(dim);

	workspace->temp_parent_approx = (VEC**) malloc(s*sizeof(VEC*));
	for(i=0;i<s;i++)
	{
		workspace->temp_parent_approx[i] = malloc(max_parents*sizeof(VEC));
		for(j=0;j<max_parents;j++)
			workspace->temp_parent_approx[i][j] = v_get(dim);
	}

	workspace->temp_k = (VEC*) malloc(s*sizeof(VEC));
	for(i=0;i<s;i++)
		workspace->temp_k[i] = v_get(dim);

	workspace->ipiv = (int*) malloc(s*dim*sizeof(int));
	workspace->RHS = v_get(s*dim);
	//workspace->CoefMat = m_get(s*dim,s*dim);
	workspace->JMatrix = m_get(dim,dim);
	workspace->Z_i = (VEC*) malloc(s*sizeof(VEC));
	for(i=0;i<s;i++)	workspace->Z_i[i] = v_get(dim);
	workspace->err = v_get(dim);

	return workspace;
}

//Deallocates workspace for RK solvers
void Destroy_Workspace(TempStorage* workspace,unsigned short int s,unsigned short int max_parents)
{
	unsigned int i,j;

	v_free(&workspace->temp);
	v_free(&workspace->sum);
	v_free(&workspace->temp2);
	v_free(&workspace->temp3);

	for(i=0;i<s;i++)
	{
		for(j=0;j<max_parents;j++)
            v_free(&workspace->temp_parent_approx[i][j]);
		free(workspace->temp_parent_approx[i]);
	}
	free(workspace->temp_parent_approx);

	for(i=0;i<s;i++)
		v_free(&workspace->temp_k[i]);
	free(workspace->temp_k);

	free(workspace->ipiv);
	v_free(&workspace->RHS);
	//m_free(&workspace->CoefMat);
	m_free(&workspace->JMatrix);
	for(i=0;i<s;i++)
        v_free(&workspace->Z_i[i]);
	free(workspace->Z_i);
	v_free(&workspace->err);
}

//Deallocates UnivVars
void Destroy_UnivVars(UnivVars* GlobalVars)
{
	unsigned int i;
	free(GlobalVars->output_data);
	free(GlobalVars->peakflow_function_name);
	free(GlobalVars->output_types);
	free(GlobalVars->output_sizes);
	for(i=0;i<GlobalVars->num_print;i++)
	{
		free(GlobalVars->output_names[i]);
		free(GlobalVars->output_specifiers[i]);
	}
	free(GlobalVars->output_names);
	free(GlobalVars->output_specifiers);
	free(GlobalVars->outputs_i);
	free(GlobalVars->outputs_d);
	if(GlobalVars->rsv_filename)	free(GlobalVars->rsv_filename);
	if(GlobalVars->rvr_filename)	free(GlobalVars->rvr_filename);
	if(GlobalVars->prm_filename)	free(GlobalVars->prm_filename);
	if(GlobalVars->init_filename)	free(GlobalVars->init_filename);
	if(GlobalVars->rain_filename)	free(GlobalVars->rain_filename);
	if(GlobalVars->dam_filename)	free(GlobalVars->dam_filename);
	if(GlobalVars->hydrosave_filename)	free(GlobalVars->hydrosave_filename);
	if(GlobalVars->peaksave_filename)	free(GlobalVars->peaksave_filename);
	if(GlobalVars->temp_filename)	free(GlobalVars->temp_filename);
	if(GlobalVars->peakfilename)	free(GlobalVars->peakfilename);
	if(GlobalVars->hydros_loc_filename)	free(GlobalVars->hydros_loc_filename);
	if(GlobalVars->peaks_loc_filename)	free(GlobalVars->peaks_loc_filename);
	if(GlobalVars->hydro_table)	free(GlobalVars->hydro_table);
	if(GlobalVars->peak_table)	free(GlobalVars->peak_table);
	if(GlobalVars->dump_table)	free(GlobalVars->dump_table);
	if(GlobalVars->dump_loc_filename)	free(GlobalVars->dump_loc_filename);
	v_free(&GlobalVars->global_params);
	free(GlobalVars->print_indices);
	free(GlobalVars);
}

//Inserts a time into the list of discontinuities
unsigned int Insert_Discontinuity(double time,unsigned int start,unsigned int end,unsigned int* count,unsigned int size,double* array,unsigned int id)
{
	if(*count == 0)
	{
		array[start] = time;
		(*count)++;
		return (end+1)%size;
	}
	else if((end+1)%size == start)
	{
		//printf("Warning: A discontinuity is being ignored. Increase the size of the discontinuity buffer. id = %u time = %f\n",id,time);
		return end;
	}

	unsigned int curr = end;
	unsigned int toofar = (start==0) ? size-1 : start-1;

	while(curr != toofar && time < array[curr])
		curr = (curr == 0) ? size-1 : curr-1;
	if(curr != toofar && time == array[curr])	return end;

	curr = end;
	while(curr != toofar && time < array[curr])
	{
		array[(curr+1)%size] = array[curr];
		curr = (curr == 0) ? size-1 : curr-1;
	}

	curr = (curr+1)%size;
	end = (end+1)%size;
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
void Insert_SendDiscontinuity(double time,unsigned int order,unsigned int* count,unsigned int size,double* array,unsigned int* order_array,unsigned int id)
{
	if(*count == 0)
	{
		array[0] = time;
		order_array[0] = order;
		*count = 1;
		return;
	}
	else if(*count+1 == size)
	{
		//printf("Warning: A discontinuity is being ignored in send. Increase the size of the discontinuity buffer. id = %i time = %f\n",id,time);
		return;
	}

	unsigned int curr = *count - 1;

	while(curr < size && time < array[curr])	curr--;
	if(curr < size && time == array[curr])
	{
		if(order_array[curr] < order)	order_array[curr] = order;
		return;
	}

	curr = *count - 1;
	while(curr < size && time < array[curr])
	{
		array[curr+1] = array[curr];
		order_array[curr+1] = order_array[curr];
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

