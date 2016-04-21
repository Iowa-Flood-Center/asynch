#ifndef RAINFALL_H
#define RAINFALL_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <stdio.h>
#include <stdlib.h>
#include <libpq-fe.h>
#include "structs.h"
#include "compression.h"
#include "comm.h"
#include "sort.h"
#include <time.h>
#include "date_manip.h"

extern int my_rank;
extern int np;


int Create_Rain_Data_Par(Link** sys,unsigned int N,unsigned int my_N,UnivVars* GlobalVars,unsigned int* my_sys,int* assignments,char strfilename[],unsigned int first,unsigned int last,double t_0,double increment,Forcing* forcing,unsigned int** id_to_loc,unsigned int max_files,unsigned int forcing_idx);

int Create_Rain_Data_GZ(Link** sys,unsigned int N,unsigned int my_N,UnivVars* GlobalVars,unsigned int* my_sys,int* assignments,char strfilename[],unsigned int first,unsigned int last,double t_0,double increment,Forcing* forcing,unsigned int** id_to_loc,unsigned int max_files,unsigned int forcing_idx);

int Create_Rain_Data_Grid(Link** sys,unsigned int N,unsigned int my_N,UnivVars* GlobalVars,unsigned int* my_sys,int* assignments,char strfilename[],unsigned int first,unsigned int last,double t_0,double increment,Forcing* forcing,unsigned int** id_to_loc,unsigned int max_files,unsigned int forcing_idx);

int Create_Rain_Database(Link** sys,unsigned int N,unsigned int my_N,UnivVars* GlobalVars,unsigned int* my_sys,int* assignments,ConnData *conninfo,unsigned int first,unsigned int last,Forcing* forcing,unsigned int** id_to_loc,double maxtime,unsigned int forcing_idx);

void SetRain0(Link** sys,unsigned int my_N,double maxtime,unsigned int* my_sys,UnivVars* GlobalVars,Forcing* forcing,unsigned int forcing_idx);

double CreateForcing_Monthly(Link** sys,unsigned int my_N,unsigned int* my_sys,UnivVars* GlobalVars,ForcingData* GlobalForcing,unsigned int forcing_idx,struct tm *current_time,time_t first_time,time_t last_time,double t_0);

int Create_Rain_Database_Irregular(Link** sys,unsigned int N,unsigned int my_N,UnivVars* GlobalVars,unsigned int* my_sys,int* assignments,ConnData *conninfo,unsigned int first,unsigned int last,Forcing* forcing,unsigned int** id_to_loc,double maxtime,unsigned int forcing_idx);

#endif

