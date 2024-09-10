#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <assert.h>
#include <stddef.h>

#include <models/check_consistency.h>


void CheckConsistency_Nonzero_1States(
    double *y, unsigned int num_dof,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    void *user)
{
    assert(y != NULL);
    assert(num_dof >= 1);

    if (y[0] < 1e-14)
        y[0] = 1e-14;
}

void CheckConsistency_Nonzero_2States(
    double *y, unsigned int num_dof,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    void *user)
{
    assert(y != NULL);
    assert(num_dof >= 2);

    if (y[0] < 1e-14)
        y[0] = 1e-14;
    if (y[1] < 0.0)
        y[1] = 0.0;
}

void CheckConsistency_Nonzero_3States(
    double *y, unsigned int num_dof,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    void *user)
{
    assert(y != NULL);
    assert(num_dof >= 3);

    if (y[0] < 1e-14)
        y[0] = 1e-14;
    if (y[1] < 0.0)
        y[1] = 0.0;
    if (y[2] < 0.0)
        y[2] = 0.0;
}

void CheckConsistency_Nonzero_4States(
    double *y, unsigned int num_dof,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    void *user)
{
    assert(y != NULL);
    assert(num_dof >= 4);

    if (y[0] < 1e-14)
        y[0] = 1e-14;
    if (y[1] < 0.0)
        y[1] = 0.0;
    if (y[2] < 0.0)
        y[2] = 0.0;
    if (y[3] < 0.0)
        y[3] = 0.0;
}


void CheckConsistency_Nonzero_5States(
    double *y, unsigned int num_dof,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    void *user)
{
    assert(y != NULL);
    assert(num_dof >= 4);

    if (y[0] < 1e-14)
        y[0] = 1e-14;
    if (y[1] < 0.0)
        y[1] = 0.0;
    if (y[2] < 0.0)
        y[2] = 0.0;
    if (y[3] < 0.0)
        y[3] = 0.0;
    if (y[4] < 0.0)
        y[4] = 0.0;
}

void CheckConsistency_Model5(
    double *y, unsigned int num_dof,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    void *user)
{
    assert(y != NULL);
    assert(num_dof >= 4);

    if (y[0] < 1e-14)
        y[0] = 1e-14;
    if (y[1] < 0.0)
        y[1] = 0.0;
    if (y[2] < 0.0)
        y[2] = 0.0;
    if (y[3] < 0.0)
        y[3] = 0.0;
    if (y[3] > 1.0 - y[2])
        y[3] = 1.0 - y[2];
}

void CheckConsistency_Model30(
    double *y, unsigned int num_dof,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    void *user)
{
    assert(y != NULL);
    assert(num_dof >= 4);

    if (y[0] < 1e-14)
        y[0] = 1e-14;
    if (y[1] < 0.0)
        y[1] = 0.0;
    if (y[2] < 0.0)
        y[2] = 0.0;
    if (y[2] > params[4])
        y[2] = params[4];
    if (y[3] < 0.0)	y[3] = 0.0;
    else if (y[3] > 1.0)
        y[3] = 1.0;
}

void CheckConsistency_Nonzero_AllStates_q(
    double *y, unsigned int num_dof,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    void *user)
{
    assert(y != NULL);
    assert(num_dof >= 1);

    if (y[0] < 1e-14) 
        y[0] = 1e-14;
    for (unsigned int i = 1; i < num_dof; i++)
        if (y[i] < 0.0)
            y[i] = 0.0;
}

void CheckConsistency_608(
    double *y, unsigned int num_dof,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    void *user)
{
    assert(y != NULL);
    assert(num_dof >= 1);

    if (y[0] < 1e-14) 
        y[0] = 1e-14;
    for (unsigned int i = 1; i < num_dof; i++)
        if (y[i] < 0.0)
            y[i] = 0.0;
    if (y[3] <= 0) 
        y[3] = params[14];

}

void CheckConsistency_Nonzero_AllStates_qs(
    double *y, unsigned int num_dof,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    void *user)
{
    assert(y != NULL);
    assert(num_dof >= 2);

    if (y[0] < 1e-14)
        y[0] = 1e-14;
    if (y[1] < 1e-14)
        y[1] = 1e-14;
    for (unsigned int i = 2; i < num_dof; i++)
        if (y[i] < 0.0)
            y[i] = 0.0;
}


