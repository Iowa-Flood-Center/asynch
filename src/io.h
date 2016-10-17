#ifndef IO_H
#define IO_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "structs.h"
#include "comm.h"
#include "modeloutputs.h"
#include "processdata.h"

#define ASYNCH_MAX_QUERY_SIZE 1024
#define ASYNCH_MAX_DB_CONNECTIONS 16
#define ASYNCH_DB_LOC_TOPO 0
#define ASYNCH_DB_LOC_PARAMS 1
#define ASYNCH_DB_LOC_INIT 2
#define ASYNCH_DB_LOC_QVS 3
#define ASYNCH_DB_LOC_RSV 4
#define ASYNCH_DB_LOC_HYDROSAVE 5
#define ASYNCH_DB_LOC_PEAKSAVE 6
#define ASYNCH_DB_LOC_HYDRO_OUTPUT 7
#define ASYNCH_DB_LOC_PEAK_OUTPUT 8
#define ASYNCH_DB_LOC_SNAPSHOT_OUTPUT 9
#define ASYNCH_DB_LOC_FORCING_START 10

extern int np;
extern int my_rank;

InOut* BuildIO(UnivVars* GlobalVars);

ConnData* ReadDBC(char* filename,unsigned int string_size);
void WriteValue(FILE* outputfile,char* specifier,char* data_storage,short int data_type,char* delim);
unsigned int WriteStep(double t,VEC y,UnivVars* GlobalVars,VEC params,unsigned int state,FILE* outputfile,void* user,long int* pos_offset);
unsigned int CatBinaryToString(char* submission,char* specifier,char* data_storage,short int data_type,char* delim);

#endif
