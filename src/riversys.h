#ifndef RIVERSYS_H
#define RIVERSYS_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <libpq-fe.h>
#include "structs.h"
#include "rkmethods.h"
#include "partition.h"
#include "comm.h"
#include "sort.h"
#include "definetype.h"
#include "misc.h"
#include "io.h"
#include "forcings.h"
#include "modeloutputs.h"

extern int np;
extern int my_rank;

//Routines to construct the system of equations to solve
//Link** Create_River_System_parallel(char rk_filename[],unsigned int* N,unsigned int** my_sys,unsigned int* my_N,unsigned int* my_max_nodes,TransData** my_data,int** assignments,short int** getting,RKMethod*** AllMethods,unsigned int* nummethods,UnivVars* GlobalVars,ErrorData* GlobalErrors,unsigned int** save_list,unsigned int* my_save_size,unsigned int* save_size,unsigned int* peaksave_size,unsigned int*** id_to_loc,Forcing** forcings,ConnData** db_connections,TempStorage** workspace,model* custom_model,void* external);
Link** Create_River_Network(UnivVars* GlobalVars,unsigned int* N,unsigned int*** id_to_loc,ConnData** db_connections);
int Load_Local_Parameters(Link** system,unsigned int N,unsigned int* my_sys,unsigned int my_N,int* assignments,short int* getting,unsigned int** id_to_loc,UnivVars* GlobalVars,ConnData** db_connections,short int load_all,model* custom_model,void* external);
int Parition_Network(Link** system,unsigned int N,UnivVars* GlobalVars,unsigned int** my_sys,unsigned int* my_N,int** assignments,TransData** my_data,short int** getting,model* custom_model);
int Build_RKData(Link** system,char rk_filename[],unsigned int N,unsigned int* my_sys,unsigned int my_N,int* assignments,short int* getting,UnivVars* GlobalVars,ErrorData* GlobalErrors,RKMethod*** AllMethods,unsigned int* nummethods);
int Initialize_Model(Link** system,unsigned int N,unsigned int* my_sys,unsigned int my_N,int* assignments,short int* getting,UnivVars* GlobalVars,model* custom_model,void* external);
int Load_Initial_Conditions(Link** system,unsigned int N,int* assignments,short int* getting,unsigned int** id_to_loc,UnivVars* GlobalVars,ConnData** db_connections,model* custom_model,void* external);
int Load_Forcings(Link** system,unsigned int N,unsigned int* my_sys,unsigned int my_N,int* assignments,short int* getting,unsigned int* res_list,unsigned int res_size,unsigned int** id_to_loc,UnivVars* GlobalVars,Forcing** forcings,ConnData** db_connections);
int Load_Dams(Link** system,unsigned int N,unsigned int* my_sys,unsigned int my_N,int* assignments,short int* getting,unsigned int** id_to_loc,UnivVars* GlobalVars,ErrorData* GlobalError,ConnData** db_connections,unsigned int** res_list,unsigned int* res_size,unsigned int* my_res_size);
int CalculateInitialStepSizes(Link** system,unsigned int* my_sys,unsigned int my_N,UnivVars* GlobalVars,TempStorage* workspace,short int watch_forcings);
int BuildSaveLists(Link** system,unsigned int N,unsigned int* my_sys,unsigned int my_N,unsigned int* assignments,unsigned int** id_to_loc,UnivVars* GlobalVars,unsigned int** save_list,unsigned int* save_size,unsigned int* my_save_size,unsigned int** peaksave_list,unsigned int* peaksave_size,unsigned int* my_peaksave_size,ConnData** db_connections);
int FinalizeSystem(Link** system,unsigned int N,unsigned int* my_sys,unsigned int my_N,unsigned int* assignments,short int* getting,unsigned int** id_to_loc,TransData* my_data,UnivVars* GlobalVars,ConnData** db_connections,TempStorage** workspace);

UnivVars* Read_Global_Data(char globalfilename[],ErrorData** GlobalErrors,Forcing** forcings,ConnData** db_connections,char* rkdfilename,model* custom_model,void* external);

int Create_SAV_Data(char filename[],Link** sys,unsigned int N,unsigned int** save_list,unsigned int* size,ConnData *conninfo,unsigned short int flag);

void ReadLineFromTextFile(FILE* globalfile,char* linebuffer,unsigned int size,unsigned int string_size);

int ReadLineError(int valsread,int valswant,char message[]);

int RemoveSuffix(char* filename,const char* suffix);

int AttachParameters(char* filename,unsigned int max_size,VEC v,unsigned int string_size);

int CheckFilenameExtension(char* filename,char* extension);
int CheckWinFormat(FILE* file);
int FindPath(char* filename,char* path);
int FindFilename(char* fullpath,char* filename);

#endif

