#include "asynch_interface_py.h"


//Utilities for mixing C and Python **************************************

void C_inc_ref(PyObject* obj)
{
	Py_INCREF(obj);
}

unsigned int* Allocate_CUINT_Array(unsigned int n)
{
	return (unsigned int*) malloc(n*sizeof(unsigned int));
}

void Free_PythonInterface(asynchsolver* asynch)
{
	free((PythonInterface*)(asynch->ExternalInterface));
}


//Routines for the Python interface **************************************

asynchsolver* Asynch_Init_py(int numprocs,int* ranks)
{
	//!!!! Assumes MPI_COMM_WORLD. This should create a new communicator and pass that along. !!!!
	asynchsolver* asynch = Asynch_Init(MPI_COMM_WORLD,NULL,NULL);

	asynch->ExternalInterface = malloc(sizeof(PythonInterface));

	return asynch;
}

void Asynch_Set_System_State_py(asynchsolver* asynch,double t_0,double* values)
{
	unsigned int i,j,N = asynch->N,*assignments = asynch->assignments;
	VEC** array = (VEC**) malloc(N*sizeof(VEC*));
	Link* current;

	//Convert input to vectors
	for(i=0;i<N;i++)
	{
		if(assignments[i] == my_rank)
		{
			current = asynch->sys[i];
			array[i] = v_get(current->dim);
			for(j=0;j<current->dim;j++)	array[i]->ve[j] = values[i*current->dim+j];
		}
		else
			array[i] = NULL;
	}

	Asynch_Set_System_State(asynch,t_0,array);

	//Trash vectors
	for(i=0;i<N;i++)	free(array[i]);
	free(array);
}

//Custom Model routines *******************************************


int Asynch_Custom_Model_py(asynchsolver* asynch,void (*SetParamSizes)(UnivVars*,PyObject*),void (*Convert)(VEC*,unsigned int,void*),void (*Routines)(Link*,unsigned int,unsigned int,unsigned short int,void*),
	void (*Precalculations)(Link*,VEC*,VEC*,unsigned int,unsigned int,unsigned short int,unsigned int,void*),int (*InitializeEqs)(VEC*,VEC*,QVSData*,unsigned short int,VEC*,unsigned int,unsigned int,unsigned int,void*,void*),PyObject* lib)
{
	//Setup the python interface data
	PythonInterface* external = (PythonInterface*) asynch->ExternalInterface;
	external->lib = lib;
	external->SetParamSizes_func = SetParamSizes;
	external->Routines_func = Routines;

	//Set the functions
	int val = Asynch_Custom_Model(asynch,&SetParamSizes_py,Convert,&InitRoutines_py,Precalculations,InitializeEqs);

	return val;
}

//Called by solvers
void SetParamSizes_py(UnivVars* GlobalVars,void* user)
{
	PythonInterface* external = (PythonInterface*) user;
	external->SetParamSizes_func(GlobalVars,external->lib);
}


void InitRoutines_py(Link* link,unsigned int type,unsigned int exp_imp,unsigned short int dam,void* user)
{
	PythonInterface* external = (PythonInterface*) user;
	external->Routines_func(link,type,exp_imp,dam,external->lib);
}

void Asynch_Copy_Local_OutputUser_Data_py(asynchsolver* asynch,unsigned int location,PyObject* source,unsigned int size)
{
	Asynch_Copy_Local_OutputUser_Data(asynch,location,(void*) source,size);
}


