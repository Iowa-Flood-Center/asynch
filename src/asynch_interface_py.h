#ifndef ASYNCH_INTERFACE_PY_H
#define ASYNCH_INTERFACE_PY_H

#include <Python.h>
#include "asynch_interface.h"

typedef struct PythonInterface
{
	PyObject* lib;
	void (*SetParamSizes_func)(UnivVars*,PyObject*);
	void (*Routines_func)(Link*,unsigned int,unsigned int,unsigned short int,void*);
} PythonInterface;

void C_inc_ref(PyObject* obj);
unsigned int* Allocate_CUINT_Array(unsigned int n);
void Free_PythonInterface(asynchsolver* asynch);

asynchsolver* Asynch_Init_py(int numprocs,int* ranks);
void Asynch_Set_System_State_py(asynchsolver* asynch,double t_0,double* values);

int Asynch_Custom_Model_py(asynchsolver* asynch,void (*SetParamSizes)(UnivVars*,PyObject*),void (*Convert)(VEC*,unsigned int,void*),void (*Routines)(Link*,unsigned int,unsigned int,unsigned short int,void*),
	void (*Precalculations)(Link*,VEC*,VEC*,unsigned int,unsigned int,unsigned short int,unsigned int,void*),int (*InitializeEqs)(VEC*,VEC*,QVSData*,unsigned short int,VEC*,unsigned int,unsigned int,unsigned int,void*,void*),PyObject* lib);
void SetParamSizes_py(UnivVars* GlobalVars,void* user);
void InitRoutines_py(Link* link,unsigned int type,unsigned int exp_imp,unsigned short int dam,void* user);
void Asynch_Copy_Local_OutputUser_Data_py(asynchsolver* asynch,unsigned int location,PyObject* source,unsigned int size);

#endif

