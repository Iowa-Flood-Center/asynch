#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <stdio.h>
#if defined(HAVE_UNISTD_H)
#  include <unistd.h>
#endif

#include "mpi.h"
#include "asynch_interface.h"

int my_rank;
int np;

int Output_Linkid(double t,VEC* y_i,VEC* global_params,VEC* params,int state,void* user);
void Set_Output_User_LinkID(asynchsolver* asynch);

typedef struct
{
	unsigned int* num_upstream;
	unsigned int** upstream;
} upstream_data;

//Custom model
void SetParamSizes_MyModel(UnivVars* GlobalVars,void* external);
void ConvertParams_MyModel(VEC* params,unsigned int type,void* external);
void InitRoutines_MyModel(Link* link,unsigned int type,unsigned int exp_imp,unsigned short int dam,void* external);
void Precalculations_MyModel(Link* link_i,VEC* global_params,VEC* params,unsigned int disk_params,unsigned int params_size,unsigned short int dam,unsigned int type,void* external);
int ReadInitData_MyModel(VEC* global_params,VEC* params,QVSData* qvs,unsigned short int dam,VEC* y_0,unsigned int type,unsigned int diff_start,unsigned int no_init_start,void* user,void* external);

//Assim model
void Setup_Errors(asynchsolver* asynch);
void Free_Upstream_Links(asynchsolver* asynch);
void Find_Upstream_Links(asynchsolver* asynch);
void SetParamSizes_Assim(UnivVars* GlobalVars,void* external);
void ConvertParams_Assim(VEC* params,unsigned int type,void* external);
void InitRoutines_Assim(Link* link,unsigned int type,unsigned int exp_imp,unsigned short int dam,void* external);
void Precalculations_Assim(Link* link_i,VEC* global_params,VEC* params,unsigned int disk_params,unsigned int params_size,unsigned short int dam,unsigned int type,void* external);
int ReadInitData_Assim(VEC* global_params,VEC* params,QVSData* qvs,unsigned short int dam,VEC* y_0,unsigned int type,unsigned int diff_start,unsigned int no_init_start,void* user,void* external);
void assim_river_rainfall_adjusted_custom(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,int state,void* user,VEC* ans);

int main(int argc,char* argv[])
{
	//Initialize MPI stuff
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&np);

	//Parse input
	if(argc < 2)
	{
		if(my_rank == 0)
		{
			printf("Command line parameter required:  A universal variable file (.gbl).\n");
			printf("\n");
		}
		MPI_Finalize();
		return 1;
	}

	//Declare variables
	time_t start,stop;
	double total_time;
	asynchsolver* asynch;

	start = time(NULL);

	if(my_rank == 0)
		printf("\nBeginning initialization...\n*****************************\n");
	MPI_Barrier(MPI_COMM_WORLD);

	//Init asynch object and the river network
	asynch = Asynch_Init(MPI_COMM_WORLD,&argc,&argv);

	//Model 191
	//Asynch_Custom_Model(asynch,&SetParamSizes_MyModel,&ConvertParams_MyModel,&InitRoutines_MyModel,&Precalculations_MyModel,&ReadInitData_MyModel);

	//Model 315
	Asynch_Custom_Model(asynch,&SetParamSizes_Assim,&ConvertParams_Assim,&InitRoutines_Assim,&Precalculations_Assim,&ReadInitData_Assim);

	if(my_rank == 0)	printf("Reading global file...\n");
	Asynch_Parse_GBL(asynch,argv[1]);
	if(my_rank == 0)	printf("Loading network...\n");
	Asynch_Load_Network(asynch);
	if(my_rank == 0)	printf("Partitioning network...\n");
	Asynch_Partition_Network(asynch);
	Find_Upstream_Links(asynch);
	if(my_rank == 0)	printf("Loading parameters...\n");
	Asynch_Load_Network_Parameters(asynch,0);
	if(my_rank == 0)	printf("Reading dam and reservoir data...\n");
	Asynch_Load_Dams(asynch);
	if(my_rank == 0)	printf("Setting up numerical error data...\n");
	Asynch_Load_Numerical_Error_Data(asynch);
	if(my_rank == 0)	printf("Initializing model...\n");
	Asynch_Initialize_Model(asynch);
	Setup_Errors(asynch);
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

	if(my_rank == 0)
	{
		printf("\nModel type is %u.\nGlobal parameters are:\n",asynch->GlobalVars->type);
		Print_Vector(asynch->GlobalVars->global_params);
		printf("\n");
	}

	//Setup output for link id, if needed
	int id_setup = Asynch_Check_Output(asynch,"LinkID");
	if(id_setup != -1)
	{
		Set_Output_User_LinkID(asynch);
		Asynch_Set_Output(asynch,"LinkID",ASYNCH_INT,(void (*)(double,VEC*,VEC*,VEC*,int,void*)) &Output_Linkid,NULL,0);
	}

	//Prepare output files
	Asynch_Prepare_Temp_Files(asynch);
	Asynch_Write_Current_Step(asynch);
	Asynch_Prepare_Peakflow_Output(asynch);
	Asynch_Prepare_Output(asynch);

	//Make sure everyone is good before getting down to it...
	printf("Process %i (%i total) is good to go with %i links.\n",my_rank,np,asynch->my_N);
	sleep(1);
	MPI_Barrier(MPI_COMM_WORLD);

	if(my_rank == 0)
	{
		stop = time(NULL);
		total_time = difftime(stop,start);
		printf("Finished initialization. Total time: %f\n\n\nComputing solution at each link...\n************************************\n",total_time);
	}
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);

	//Perform the calculations
	time(&start);
	Asynch_Advance(asynch,1);
	MPI_Barrier(MPI_COMM_WORLD);
	time(&stop);

	//Out information
	total_time += difftime(stop,start);
	if(my_rank == 0)	printf("\nComputations complete. Total time for calculations: %f\n",difftime(stop,start));
	if(asynch->sys[asynch->my_sys[0]]->c == NULL)
	{
		printf("[%i]: The solution at ID %i at time %.12f is\n",my_rank,asynch->sys[asynch->my_sys[0]]->ID,asynch->sys[asynch->my_sys[0]]->last_t);
		Print_Vector(asynch->sys[asynch->my_sys[0]]->list->tail->y_approx);
	}

	//Take a snapshot
	Asynch_Take_System_Snapshot(asynch,NULL);

	//Create output files
	Asynch_Create_Output(asynch,NULL);
	Asynch_Create_Peakflows_Output(asynch);

	//Cleanup
	Free_Upstream_Links(asynch);
	Asynch_Delete_Temporary_Files(asynch);
	Asynch_Free(asynch);
	return 0;
}



int Output_Linkid(double t,VEC* y_i,VEC* global_params,VEC* params,int state,void* user)
{
	return ((Link*)user)->ID;
}

//!!!! Gross, but not sure how else to handle this. Maybe with a lot of interface functions? !!!!
void Set_Output_User_LinkID(asynchsolver* asynch)
{
	unsigned int i,my_N = asynch->my_N,*my_sys = asynch->my_sys;
	Link** sys = asynch->sys;

	for(i=0;i<my_N;i++)
		sys[my_sys[i]]->output_user = (void*) sys[my_sys[i]];
}


//Custom model parameters

//Model 191 ************************************************************************************

void SetParamSizes_MyModel(UnivVars* GlobalVars,void* external)
{

	//num_global_params = 7;
	GlobalVars->uses_dam = 0;
	GlobalVars->params_size = 8;
	GlobalVars->dam_params_size = 0;
	GlobalVars->area_idx = 0;
	GlobalVars->areah_idx = 2;
	GlobalVars->disk_params = 3;
	GlobalVars->convertarea_flag = 0;
	GlobalVars->num_forcings = 3;
}


void ConvertParams_MyModel(VEC* params,unsigned int type,void* external)
{
	params->ve[1] *= 1000;	//L: km -> m
	params->ve[2] *= 1e6;	//A_h: km^2 -> m^2
}


void InitRoutines_MyModel(Link* link,unsigned int type,unsigned int exp_imp,unsigned short int dam,void* external)
{
	link->dim = 6;
	link->no_ini_start = 3;
	link->diff_start = 0;

	link->num_dense = 2;
	link->dense_indices = (unsigned int*) malloc(link->num_dense*sizeof(unsigned int));
	link->dense_indices[0] = 0;
	link->dense_indices[1] = 5;

	if(link->res)
	{
		link->f = &LinearHillslope_Reservoirs_extras;
		link->RKSolver = &ForcedSolutionSolver;
	}
	else
	{
		link->f = &LinearHillslope_MonthlyEvap_extras;
		link->RKSolver = &ExplicitRKSolver;
	}
	link->alg = NULL;
	link->state_check = NULL;
	link->CheckConsistency = &CheckConsistency_Nonzero_AllStates_q;
}


void Precalculations_MyModel(Link* link_i,VEC* global_params,VEC* params,unsigned int disk_params,unsigned int params_size,unsigned short int dam,unsigned int type,void* external)
{
	//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,c_1,c_2
	//The numbering is:	0   1   2   3  4    5    6   7
	//Order of global_params: v_r,lambda_1,lambda_2,RC,v_h,v_g (,v_B)
	//The numbering is:        0      1        2     3  4   5     6
	double* vals = params->ve;
	double A_i = params->ve[0];
	double L_i = params->ve[1];
	double A_h = params->ve[2];
	double v_r = global_params->ve[0];
	double lambda_1 = global_params->ve[1];
	double lambda_2 = global_params->ve[2];
	double RC = global_params->ve[3];
	double v_h = global_params->ve[4];
	double v_g = global_params->ve[5];

	vals[3] = v_h * L_i / A_h * 60.0;	//[1/min]  k2
	vals[4] = v_g * L_i / A_h * 60.0;	//[1/min]  k3
	vals[5] = 60.0*v_r*pow(A_i,lambda_2) / ((1.0-lambda_1)*L_i);	//[1/min]  invtau
	vals[6] = RC*(0.001/60.0);		//(mm/hr->m/min)  c_1
	vals[7] = (1.0-RC)*(0.001/60.0);	//(mm/hr->m/min)  c_2
}

int ReadInitData_MyModel(VEC* global_params,VEC* params,QVSData* qvs,unsigned short int dam,VEC* y_0,unsigned int type,unsigned int diff_start,unsigned int no_init_start,void* user,void* external)
{
	//For this type, the extra states need to be set (3,4,5)
	y_0->ve[3] = 0.0;
	y_0->ve[4] = 0.0;
	y_0->ve[5] = y_0->ve[0];	//I'm not really sure what to use here...
	return 0;
}


//Data assimilation model (Old Model 315) ************************************************************************************

void Setup_Errors(asynchsolver* asynch)
{
	UnivVars* GlobalVars = asynch->GlobalVars;
	ErrorData* GlobalErrors = asynch->GlobalErrors;
	unsigned int i,problem_dim = 2,max_dim = GlobalVars->max_dim;

	GlobalErrors->abstol->ve = realloc(GlobalErrors->abstol->ve,max_dim*sizeof(double));
	GlobalErrors->reltol->ve = realloc(GlobalErrors->reltol->ve,max_dim*sizeof(double));
	GlobalErrors->abstol_dense->ve = realloc(GlobalErrors->abstol_dense->ve,max_dim*sizeof(double));
	GlobalErrors->reltol_dense->ve = realloc(GlobalErrors->reltol_dense->ve,max_dim*sizeof(double));
	GlobalErrors->abstol->dim = GlobalErrors->reltol->dim = GlobalErrors->reltol_dense->dim = GlobalErrors->reltol_dense->dim = max_dim;

	//Setup error
	for(i=problem_dim+1;i<max_dim;i++)
	{
		GlobalErrors->abstol->ve[i] = GlobalErrors->abstol->ve[problem_dim];
		GlobalErrors->reltol->ve[i] = GlobalErrors->reltol->ve[problem_dim];
		GlobalErrors->abstol_dense->ve[i] = GlobalErrors->abstol_dense->ve[problem_dim];
		GlobalErrors->reltol_dense->ve[i] = GlobalErrors->reltol_dense->ve[problem_dim];
	}
}

void Free_Upstream_Links(asynchsolver* asynch)
{
	upstream_data *data;
	Link** sys = asynch->sys;
	unsigned int N = asynch->N,i,j;

	for(i=0;i<N;i++)
	{
		data = (upstream_data*) (sys[i]->user);
		if(data)
		{
			for(j=0;j<sys[i]->numparents;j++)
				if(data->upstream[j])	free(data->upstream[j]);
			free(data->upstream);
			free(data->num_upstream);
			free(data);
			sys[i]->user = NULL;
		}
	}
}

void Find_Upstream_Links(asynchsolver* asynch)
{
	Link **sys = asynch->sys,*current;
	unsigned int N = asynch->N,parentsval,leaves_size = 0,i,j,l,m;
	int *assignments = asynch->assignments;
	UnivVars *GlobalVars = asynch->GlobalVars;
	Link **leaves = (Link**) malloc(N*sizeof(Link*));
	Link **stack = (Link**) malloc(N*sizeof(Link*));
	unsigned short int* getting = asynch->getting;

	//Find leaves
	for(i=0;i<N;i++)
		if(sys[i]->numparents == 0)	leaves[leaves_size++] = sys[i];

	unsigned int* temp_numupstream = (unsigned int*) calloc(N,sizeof(unsigned int));
	for(i=0;i<leaves_size;i++)
		temp_numupstream[leaves[i]->location] = 1;

	//Count upstream links
	for(i=0;i<leaves_size;i++)
	{
		for(current = leaves[i]->c; current != NULL; current = current->c)
		{
			parentsval = 0;
			for(j=0;j<current->numparents;j++)	parentsval += (temp_numupstream[current->parents[j]->location] > 0);

			if(parentsval == current->numparents)	//All parents have temp_numupstream set
			{
				temp_numupstream[current->location] = 1;
				for(j=0;j<current->numparents;j++)
					temp_numupstream[current->location] += temp_numupstream[current->parents[j]->location];
			}
			else
				break;
		}
	}

	//Set the upstream links
	unsigned int** temp_upstream = (unsigned int**) malloc(N*sizeof(unsigned int*));	//temp_upstream[i] is list of all links upstream from link i
	for(i=0;i<N;i++)
		temp_upstream[i] = (unsigned int*) malloc(temp_numupstream[sys[i]->location] * sizeof(unsigned int));
	unsigned int* counter = (unsigned int*) calloc(N,sizeof(unsigned int));

	unsigned int stack_size = leaves_size;
	for(i=0;i<leaves_size;i++)	stack[i] = leaves[i];

	while(stack_size > 0)
	{
		current = stack[stack_size-1];
		l = current->location;

		//Add this link to its own upstream list
		temp_upstream[l][counter[l]] = l;
		counter[l]++;

		//Add each parents upstream list
		for(i=0;i<current->numparents;i++)
		{
			m = current->parents[i]->location;
			for(j=0;j<counter[m];j++)
				temp_upstream[l][counter[l]+j] = temp_upstream[m][j];
			counter[l] += counter[m];
		}

		stack_size--;

		//If every parent of current's child has an upstream list determined, add it to the stack
		if(current->c != NULL)
		{
			parentsval = 0;
			for(i=0;i<current->c->numparents;i++)
			{
				m = current->c->parents[i]->location;
				parentsval += (counter[m] > 0);
			}

			if(parentsval == current->c->numparents)
			{
				stack[stack_size] = current->c;
				stack_size++;
			}
		}
	}

	//Move the data from temp_upstream into the child upstream
	upstream_data* data;
	short int* used = (short int*) calloc(N,sizeof(short int));	//1 if temp_upstream[i] was used, 0 if not
	for(i=0;i<N;i++)
	{
		if(assignments[i] == my_rank || getting[i])
		{
			sys[i]->user = malloc(sizeof(upstream_data));
			data = (upstream_data*) (sys[i]->user);

			data->upstream = (unsigned int**) malloc(sys[i]->numparents * sizeof(unsigned int*));
			data->num_upstream = (unsigned int*) malloc(sys[i]->numparents * sizeof(unsigned int));
			for(j=0;j<sys[i]->numparents;j++)
			{
				data->upstream[j] = temp_upstream[sys[i]->parents[j]->location];
				data->num_upstream[j] = temp_numupstream[sys[i]->parents[j]->location];
				used[sys[i]->parents[j]->location] = 1;
			}
		}
	}

	//Cleanup
	for(i=0;i<N;i++)
	{
		if(!used[i])	free(temp_upstream[i]);
		//if(sys[i]->c == NULL)	free(temp_upstream[i]);
		//else	if(assignments[i] != my_rank && !getting[i])	free(temp_upstream[i]);
	}
	free(temp_upstream);
	free(temp_numupstream);
	free(counter);
	free(leaves);
	free(stack);
	free(used);
/*
	//Try removing low order links from the upstream lists
	printf("!!!! Removing low order links from upstream list...!!!!\n");
	unsigned int cut_off = 2,drop;
	unsigned int* order = (unsigned int*) malloc(N*sizeof(unsigned int));
	unsigned short int* complete = (unsigned short int*) malloc(N*sizeof(unsigned short int));
	CalcHortonOrder(sys,N,order,complete);

	for(i=0;i<N;i++)
	{
		for(j=0;j<sys[i]->numparents;j++)
		{
			drop = 0;
//printf("checking %i (%i)\n",sys[i]->ID,sys[i]->parents[j]->ID);
			for(k=0;k<sys[i]->numupstream[j];k++)
			{
				if(order[sys[i]->upstream[j][k]] <= cut_off)	//!!!! Once this happens, the rest of the links can be dropped !!!!
				{
					drop++;
				}
				else
				{
					sys[i]->upstream[j][k-drop] = sys[i]->upstream[j][k];
				}
			}
//printf("%i going to %i from %i\n",sys[i]->ID,sys[i]->numupstream[j]-drop,sys[i]->numupstream[j]);
			sys[i]->numupstream[j] -= drop;
		}
	}

	free(order);
	free(complete);
*/
}

void SetParamSizes_Assim(UnivVars* GlobalVars,void* external)
{
	GlobalVars->uses_dam = 0;
	GlobalVars->params_size = 20;
	GlobalVars->dam_params_size = 0;
	GlobalVars->area_idx = 2;
	GlobalVars->areah_idx = 1;
	GlobalVars->disk_params = 12;
	GlobalVars->convertarea_flag = 0;
	GlobalVars->num_forcings = 1;
}


void ConvertParams_Assim(VEC* params,unsigned int type,void* external)
{
	params->ve[0] *= 1000;	//L: km -> m
	params->ve[3] *= .001;	//h_b: mm -> m
	params->ve[4] *= .001;	//h_H: mm -> m
}

void InitRoutines_Assim(Link* link,unsigned int type,unsigned int exp_imp,unsigned short int dam,void* external)
{
	upstream_data* data = (upstream_data*) (link->user);
	unsigned int i,problem_dim = 2;	//Number of model eqs

	link->dim = problem_dim + problem_dim + (problem_dim-1)*(problem_dim-1);	//Model eqs + variational eqs from this link
	for(i=0;i<link->numparents;i++)
		link->dim += data->num_upstream[i] * problem_dim;	//Variational eqs from upstream !!!! Too high? !!!!
	link->no_ini_start = 2;
	link->diff_start = 0;

	link->num_dense = link->dim - 1;	//Take out s_p
	link->dense_indices = (unsigned int*) malloc(link->num_dense*sizeof(unsigned int));
	link->dense_indices[0] = 0;
	for(i=1;i<link->num_dense;i++)	link->dense_indices[i] = i+1;

	link->f = &assim_river_rainfall_adjusted_custom;
	link->alg = NULL;
	link->state_check = NULL;
	link->CheckConsistency = &CheckConsistency_Nonzero_2States;
	link->RKSolver = &ExplicitRKSolver;
}


void Precalculations_Assim(Link* link_i,VEC* global_params,VEC* params,unsigned int disk_params,unsigned int params_size,unsigned short int dam,unsigned int type,void* external)
{
	//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
	//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
	//Order of global_params: v_r,lambda_1,lambda_2,v_h,A_r,RC
	//The numbering is:        0      1        2     3   4   5
	//Need to set entries 12-19 of params.
	double* vals = params->ve;
	double K_T = 1.0;
	double s_r = 1.0;
	double rootS_h = pow(vals[7],.5);
	double L = params->ve[0];
	double A_h = params->ve[1] * 1e6;	//Put into m^2
	double eta = params->ve[8];
	double v_r = global_params->ve[0];
	double lambda_1 = global_params->ve[1];
	double lambda_2 = global_params->ve[2];
	double v_h = global_params->ve[3];
	double A_r = global_params->ve[4];
	double RC = global_params->ve[5];

	//!!!! Clean this model. You don't really need 20 parameters... !!!!
	vals[12] = 60.0*v_r*pow(vals[2]/A_r,lambda_2)/((1.0-lambda_1)*vals[0]);	//invtau [1/min]
	vals[13] = vals[3] / s_r; //epsilon
	vals[14] = v_h*L;	//c_1 [m^2/s]
	vals[15] = vals[6] * vals[0] * vals[3] / 3600.0; //c_2
	vals[16] = (1e-3/60.0) * RC;	//c_3
	vals[17] = 60.0*v_h*L/A_h;	//c_4 [1/min], A_h converted above
	vals[18] = K_T/60.0;
	vals[19] = vals[6]/(60.0*s_r);

	//iparams->ve[0] = link_i->location; //!!!! Is this even needed anywhere? !!!!
}

int ReadInitData_Assim(VEC* global_params,VEC* params,QVSData* qvs,unsigned short int dam,VEC* y_0,unsigned int type,unsigned int diff_start,unsigned int no_init_start,void* user,void* external)
{
	//For this type, all initial conditions for variational equation must be set here.
	//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
	//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
	//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
	//The numbering is:        0      1        2     3   4   5
	unsigned int i;
	unsigned int offset = 2;

	y_0->ve[offset] = 1.0;
	y_0->ve[offset + 1] = 1.0;
	y_0->ve[offset + 2] = 0.0;
	for(i=offset+3;i<y_0->dim;i++)	y_0->ve[i] = 0.0;

	return 0;
}

//Function for simple river system with data assimilation.
//Calculates the flow using simple parameters, using only the flow q.
//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
//Order of global_params: v_r,lambda_1,lambda_2,v_h,A_r,RC
//The numbering is:        0      1        2     3   4   5
//This uses the units and functions from September 18, 2011 document
//y_i[0] = q, y_i[1] = s, followed by N entries for the variational equation
void assim_river_rainfall_adjusted_custom(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,int state,void* user,VEC* ans)
{
	unsigned int i,j;
	unsigned int dim = ans->dim;
	unsigned int offset = 2;		//!!!! This needs to be num_dense, but without variational eqs !!!!
	unsigned int parent_offset;
	unsigned int problem_dim = 2;
	unsigned int all_states = (dim-offset)/problem_dim;
	double inflow = 0.0;
	upstream_data* data = (upstream_data*) user;

	double q = y_i->ve[0];
	double s_p = y_i->ve[1];

	double L = params->ve[0];
	double invtau = params->ve[12];
	double c_1 = params->ve[14];
	double c_3 = params->ve[16];
	double c_4 = params->ve[17];
	double lambda_1 = global_params->ve[1];

	double q_to_lambda_1 = pow(q,lambda_1);
	double q_to_lambda_1_m1 = (q > 1e-12) ? q_to_lambda_1 / q : pow(1e-12,lambda_1 - 1.0);
	double deriv_qpl = 1.0;

	double q_pl = s_p;

	//Flux equation (y_i[0])
	ans->ve[0] = -q + c_1 * q_pl;
	for(i=0;i<numparents;i++)
		inflow += y_p[i]->ve[0];
	ans->ve[0] = invtau * q_to_lambda_1 * (inflow + ans->ve[0]);

	//Ponded water equation (y_i[1])
	ans->ve[1] = c_3 * forcing_values[0] - c_4 * q_pl;
	//ans->ve[1] = c_3 * ( max(forcing_values[0] + 20.0*sin(t/5.0),0.0)) - c_4 * q_pl;

	//!!!! Pull if statements out of loops (should just need two cases total) !!!!
	//!!!! A lot of terms get repeated !!!!

	//Eqs for variational equations
	for(i=offset;i<dim;i++)	ans->ve[i] = 0.0;

	//s variable from this link
	ans->ve[offset] = -c_4*deriv_qpl*y_i->ve[offset];

	//q variables from this link
//	if(lambda_1 > 1e-12 && (inflow) > 1e-12)
		ans->ve[offset + 1] = (lambda_1 * invtau * q_to_lambda_1_m1 * (inflow + c_1*s_p) - (lambda_1 + 1) * invtau * q_to_lambda_1) * y_i->ve[offset + 1];
//	else
//		ans->ve[offset + 1] = -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i->ve[offset + 1];

//	if(lambda_1 > 1e-12 && (inflow) > 1e-12)
		ans->ve[offset + 2] = (lambda_1 * invtau * q_to_lambda_1_m1 * (inflow + c_1*s_p) - (lambda_1 + 1) * invtau * q_to_lambda_1) * y_i->ve[offset + 2] + invtau*c_1*q_to_lambda_1*deriv_qpl * y_i->ve[offset];
//	else
//		ans->ve[offset + 2] = -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i->ve[offset + 2] + invtau*c_1*deriv_qpl*y_i->ve[offset];

	//Adjust offset
	offset += 3;

	//Variables from parents
	for(i=0;i<numparents;i++)
	{
		parent_offset = 1 + problem_dim;

		for(j=0;j<data->num_upstream[i];j++)
		{
			ans->ve[offset] = invtau * q_to_lambda_1 * y_p[i]->ve[parent_offset];
//			if(lambda_1 > 1e-12 && (inflow) > 1e-12)
				ans->ve[offset] += (lambda_1 * invtau * q_to_lambda_1_m1 * (inflow + c_1*s_p) - (lambda_1 + 1) * invtau * q_to_lambda_1) * y_i->ve[offset];
//			else
//				ans->ve[offset] += -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i->ve[offset];

			ans->ve[offset + 1] = invtau * q_to_lambda_1 * y_p[i]->ve[parent_offset + 1];
//			if(lambda_1 > 1e-12 && (inflow) > 1e-12)
				ans->ve[offset + 1] += (lambda_1 * invtau * q_to_lambda_1_m1 * (inflow + c_1*s_p) - (lambda_1 + 1) * invtau * q_to_lambda_1) * y_i->ve[offset + 1];
//			else
//				ans->ve[offset + 1] += -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i->ve[offset + 1];

			offset += 2;
			parent_offset += 2;
		}
	}
}

