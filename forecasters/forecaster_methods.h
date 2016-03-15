#ifndef FORECASTER_METHODS_H
#define FORECASTER_METHODS_H

#include "structs.h"
#include "comm.h"
#include "riversys.h"
#include <time.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#if !defined(_MSC_VER)
#include <unistd.h>
#include <libssh2.h>
#include <arpa/inet.h>
#else
#include <sys/stat.h>
#endif


typedef struct ForecastData
{
	char* model_name;
	char* halt_filename;
	unsigned int forecasting_forcing;
	unsigned int num_rainsteps;
	short int ifis_display;
	char* rainmaps_filename;
	ConnData* rainmaps_db;
	double forecast_window;
} ForecastData;

int DeleteFutureValues(ConnData* conninfo,unsigned int num_tables,UnivVars* GlobalVars,char* table_name,char* model_name,unsigned int clear_after,unsigned int equal,char* schema);
void PerformTableMaintainance(ConnData* conninfo_hydros,UnivVars* GlobalVars,ForecastData* Forecaster,short int* vac,short unsigned int hr1,unsigned int num_tables,char* tablename,char* schema);
void CheckPartitionedTable(ConnData* conninfo,UnivVars* GlobalVars,ForecastData* Forecaster,unsigned int num_tables,char* tablename,char* colname,char* schema);
void CreateHaltFile(char* filename);
short int CheckFinished(char* filename);
int WaitForDB(ConnData* conninfo,unsigned int naptime,int stale_time,unsigned int query_size);
void FreeDBLock(ConnData* conninfo,unsigned int query_size);
ForecastData* Init_ForecastData(char* fcst_filename,unsigned int string_size);
void Free_ForecastData(ForecastData** Forecaster);
int SendFilesTo51(char* loclfile,char* serverlocation);

#endif

