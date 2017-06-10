#ifndef PROCESSDATA_H
#define PROCESSDATA_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <structs.h>
#include <sort.h>
#include <comm.h>
#include <rkmethods.h>
#include <data_types.h>

#define DB_CONNS_AT_ONCE 10

extern int np;
extern int my_rank;

int DumpTimeSerieFile(Link* sys, GlobalVars* GlobalVars, unsigned int N, unsigned int* save_list, unsigned int save_size, unsigned int my_save_size, const Lookup * const id_to_loc, int* assignments, char* additional_temp, char* additional_out, ConnData* conninfo, FILE** my_tempfile);

int DumpTimeSerieDatFile(Link* sys, GlobalVars* globals, unsigned int N, unsigned int* save_list, unsigned int save_size, unsigned int my_save_size, const Lookup * const id_to_loc, int* assignments, char* additional_temp, char* additional_out);
int DumpTimeSerieCsvFile(Link* sys, GlobalVars* globals, unsigned int N, unsigned int* save_list, unsigned int save_size, unsigned int my_save_size, const Lookup * const id_to_loc, int* assignments, char* additional_temp, char* additional_out);
int DumpTimeSerieH5File(Link* sys, GlobalVars* globals, unsigned int N, unsigned int* save_list, unsigned int save_size, unsigned int my_save_size, const Lookup * const id_to_loc, int* assignments, char* additional_temp, char* additional_out);
int DumpTimeSerieNcFile(Link* sys, GlobalVars* globals, unsigned int N, unsigned int* save_list, unsigned int save_size, unsigned int my_save_size, const Lookup * const id_to_loc, int* assignments, char* additional_temp, char* additional_out);

#if defined(HAVE_POSTGRESQL)
void PrepareDatabaseTable(GlobalVars* GlobalVars, ConnData* conninfo);
int DumpTimeSerieDB(Link* sys, GlobalVars* GlobalVars, unsigned int N, unsigned int* save_list, unsigned int save_size, unsigned int my_save_size, const Lookup * const id_to_loc, int* assignments, char* additional_temp, char* additional_out, ConnData* conninfo, FILE** my_tempfile);
int DumpStateDB(Link* sys, unsigned int N, int* assignments, GlobalVars* GlobalVars, char* preface, ConnData* conninfo);
int DumpPeakFlowDB(Link* sys, GlobalVars* GlobalVars, unsigned int N, int* assignments, unsigned int* peaksave_list, unsigned int peaksave_size, const Lookup * const id_to_loc, ConnData* conninfo);
#endif //HAVE_POSTGRESQL

//void DataDump(Link* sys,unsigned int N,int* assignments,UnivVars* GlobalVars,unsigned int last_file);
int DumpStateText(Link* sys, unsigned int N, int* assignments, GlobalVars* GlobalVars, char* preface, ConnData* conninfo);
int DumpStateH5(Link* sys, unsigned int N, int* assignments, GlobalVars* GlobalVars, char* preface, ConnData* conninfo);

int PreparePeakFlowFiles(GlobalVars* GlobalVars, unsigned int peaksave_size);
int DumpPeakFlowText(Link* sys, GlobalVars* GlobalVars, unsigned int N, int* assignments, unsigned int* peaksave_list, unsigned int peaksave_size, const Lookup * const id_to_loc, ConnData* conninfo);

//Temporary files
FILE* PrepareTempFiles(Link* sys, unsigned int N, int* assignments, GlobalVars* GlobalVars, unsigned int* save_list, unsigned int save_size, unsigned int my_save_size, char* additional, const Lookup * const id_to_loc);
int RemoveTemporaryFiles(GlobalVars* GlobalVars, unsigned int my_save_size, char* additional_temp);
int ResetTempFiles(double set_time, Link* sys, unsigned int N, FILE* tempfile, GlobalVars* GlobalVars, unsigned int my_save_size, const Lookup * const id_to_loc);
int SetTempFiles(double set_time, void* set_value, enum AsynchTypes data_type, unsigned int component_idx, Link* sys, unsigned int N, FILE* tempfile, GlobalVars* GlobalVars, unsigned int my_save_size, const Lookup * const id_to_loc);

void LoadRecoveryFile(char* filename, Link* sys, unsigned int N, unsigned int my_N, unsigned int* assignments, GlobalVars* GlobalVars);

//Utility
//int ConvertBinaryToString(double* data_storage, char* submission, unsigned int blocks, unsigned int dimp1, unsigned int id);
int overwrite_last_step(Link* link_i, GlobalVars *GlobalVars, FILE* outputfile);

#endif

