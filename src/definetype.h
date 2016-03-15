#ifndef DEFINETYPE_H
#define DEFINETYPE_H

#include <stdio.h>
#include <stdlib.h>
#if !defined(_MSC_VER)
#include <unistd.h>
#endif
#include "structs.h"
#include "system.h"
#include "rkmethods.h"
#include "problems.h"
#include "math.h"

void SetParamSizes(UnivVars* GlobalVars,void* external);
void ConvertParams(VEC* params,unsigned int type,void* external);
void InitRoutines(Link* link,unsigned int type,unsigned int exp_imp,unsigned short int dam,void* external);
void Precalculations(Link* link_i,VEC* global_params,VEC* params,unsigned int disk_params,unsigned int params_size,unsigned short int dam,unsigned int type,void* external);
int ReadInitData(VEC* global_params,VEC* params,QVSData* qvs,unsigned short int dam,VEC* y_0,unsigned int type,unsigned int diff_start,unsigned int no_init_start,void* user,void* external);
//void AssimError(unsigned int N,UnivVars* GlobalVars,ErrorData* GlobalErrors);

#endif
