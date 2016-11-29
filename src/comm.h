#ifndef COMM_H
#define COMM_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#include "structs.h"
#include "system.h"
#include "sort.h"
#include "libpq-fwd.h"

extern int my_rank;
extern int np;

#define ASYNCH_MAX_NUMBER_OF_PROCESS 256

//MPI Related Methods
void Transfer_Data(TransData* my_data,Link* sys,int* assignments,GlobalVars* GlobalVars);
void Transfer_Data_Finish(TransData* my_data,Link* sys,int* assignments,GlobalVars* GlobalVars);
void Exchange_InitState_At_Forced(Link* system,unsigned int N,unsigned int* assignments,short int* getting,unsigned int* res_list,unsigned int res_size,unsigned int** id_to_loc,GlobalVars* GlobalVars);
TransData* Initialize_TransData();
void Flush_TransData(TransData* data);
void TransData_Free(TransData* data);


//PostgreSQL Database Methods
void ConnData_Init(ConnData* const conn, const char* connstring);
void ConnData_Free(ConnData* const conninfo);

int ConnectPGDB(ConnData* conninfo);
void DisconnectPGDB(ConnData* conninfo);

int CheckResError(PGresult* res,char* event);
int CheckResState(PGresult* res,short int error_code);
void CheckConnConnection(ConnData* conninfo);

#endif //COMM_H
