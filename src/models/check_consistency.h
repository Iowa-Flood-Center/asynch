#if !defined(ASYNCH_MODEL_CHECK_CONSISTENCY_H)
#define ASYNCH_MODEL_CHECK_CONSISTENCY_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


// Shared consistency checking functions
void CheckConsistency_Nonzero_1States(
    double *y, unsigned int num_dof,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,    
    void *user);

void CheckConsistency_Nonzero_2States(
    double *y, unsigned int num_dof,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    void *user);

void CheckConsistency_Nonzero_3States(
    double *y, unsigned int num_dof,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    void *user);

void CheckConsistency_Nonzero_4States(
    double *y, unsigned int num_dof,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    void *user);

void CheckConsistency_Nonzero_5States(
    double *y, unsigned int num_dof,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    void *user);

void CheckConsistency_Model5(
    double *y, unsigned int num_dof,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    void *user);

//TODO mode to specific model
void CheckConsistency_Model30(
    double *y, unsigned int num_dof,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    void *user);

void CheckConsistency_Nonzero_AllStates_q(
    double *y, unsigned int num_dof,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    void *user);

void CheckConsistency_608(
    double *y, unsigned int num_dof,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    void *user);

void CheckConsistency_Nonzero_AllStates_qs(
    double *y, unsigned int num_dof,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    void *user);

#endif //!defined(ASYNCH_MODEL_CHECK_CONSISTENCY_H)

