#if !defined(ASYNCH_MODEL_CHECK_STATE_H)
#define ASYNCH_MODEL_CHECK_STATE_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include <stdbool.h>

#include <structs_fwd.h>


//Type 21
//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
//The numbering is:	0   1   2  3  4    5	       6      7       8     9	  10	    11       12  13  14
//Order of global_params: v_r,lambda_1,lambda_2,RC,S_0,v_h,v_g
//The numbering is:        0      1        2     3  4   5   6
int dam_check(
    double *y, unsigned int num_dof,
    const double * const params, unsigned int num_params,
    const double * const global_params, unsigned int num_global_params,
    QVSData *qvs,
    bool has_dam,
    void *user);

//Type 22
//Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
//The numbering is:	0   1   2  3   4   5   6  7   8          9	  10	  11   12     13      14        15  16   17
//Order of global_params: lambda_1,lambda_2,S_0,v_g
//The numbering is:         0        1       2   3
int dam_check2(
    double *y, unsigned int num_dof,
    const double * const params, unsigned int num_params,
    const double * const global_params, unsigned int num_global_params,
    QVSData *qvs,
    bool has_dam,
    void *user);

//Type 23
//Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
//The numbering is:	0   1   2  3   4   5   6  7   8          9	  10	  11   12     13      14        15  16   17
//Order of global_params: lambda_1,lambda_2,S_0,v_g
//The numbering is:         0        1       2   3
int dam_check3(
    double *y, unsigned int num_dof,
    const double * const params, unsigned int num_params,
    const double * const global_params, unsigned int num_global_params,
    QVSData *qvs,
    bool has_dam,
    void *user);

//Type 40 / 261 / 262
int dam_check_qvs(
    double *y, unsigned int num_dof,
    const double * const params, unsigned int num_params,
    const double * const global_params, unsigned int num_global_params,
    QVSData *qvs,
    bool has_dam,
    void *user);

//Type 402
int dam_check_qvs_402(
    double *y, unsigned int num_dof,
    const double * const params, unsigned int num_params,
    const double * const global_params, unsigned int num_global_params,
    QVSData *qvs,
    bool has_dam,
    void *user);

//Type 403
int dam_check_qvs_403(
    double *y, unsigned int num_dof,
    const double * const params, unsigned int num_params,
    const double * const global_params, unsigned int num_global_params,
    QVSData *qvs,
    bool has_dam,
    void *user);
    
#endif //!defined(ASYNCH_MODEL_CHECK_STATE_H)
