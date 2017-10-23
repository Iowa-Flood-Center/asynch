#if !defined(ASYNCH_MODELS_DEFINITIONS_H)
#define ASYNCH_MODELS_DEFINITIONS_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include <structs_fwd.h>


void SetParamSizes(
    GlobalVars* globals,
    void* external);

void SetOutputConstraints(
    GlobalVars* globals);

void ConvertParams(
    double *params,
    unsigned int type,
    void* external);

void InitRoutines(
    Link* link,
    unsigned int type,
    unsigned int exp_imp,
    unsigned short dam,
    void* external);

void Precalculations(
    Link* link_i,
    double *global_params, unsigned int num_global_params,
    double *params, unsigned int num_disk_params, unsigned int num_params,
    unsigned short dam,
    unsigned int type,
    void* external);

int ReadInitData(
    double *global_params, unsigned int num_global_params,
    double *params, unsigned int num_params,
    QVSData* qvs,
    unsigned short int dam,
    double *y_0, unsigned int dim,
    unsigned int type,
    unsigned int diff_start, unsigned int no_init_start,
    void* user,
    void* external);



//void AssimError(unsigned int N,UnivVars* GlobalVars,ErrorData* GlobalErrors);

#endif //ASYNCH_MODELS_DEFINITIONS_H
