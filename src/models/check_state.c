#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <assert.h>
#include <math.h>

#include <structs.h>
#include <models/check_state.h>


//Type 40 / 261 / 262 /
int dam_check_qvs(
    double *y, unsigned int num_dof,
    const double * const params, unsigned int num_params,
    const double * const global_params, unsigned int num_global_params,
    QVSData *qvs,
    bool has_dam,
    void *user)
{
    unsigned int i, iterations;
    double S = y[1];

    if (!has_dam)
        return -1;

    iterations = qvs->n_values - 1;
    for (i = 0; i<iterations; i++)
    {
        if (qvs->points[i][0] <= S && S < qvs->points[i + 1][0])
            return i;
    }

    return i;
}

int dam_check_qvs_402(double *y, unsigned int num_dof,
    const double * const params, unsigned int num_params,
    const double * const global_params, unsigned int num_global_params,
    QVSData *qvs,
    bool has_dam,
    void *user)
    {
    //unsigned int i, iterations;
    //double S = y[6]; //model 402 storage is state 6
    int debug =1;
    

    if (!has_dam)
        return 0;
    if(has_dam){
        if(debug) printf("found dam_check_qvs_402\n");
        return 1;
    }
        
    
    
}

int dam_check_qvs_403(double *y, unsigned int num_dof,
    const double * const params, unsigned int num_params,
    const double * const global_params, unsigned int num_global_params,
    QVSData *qvs,
    bool has_dam,
    void *user)
    {
    unsigned int i, iterations;
    double S = y[6]; //model 403 storage is state 6
    int debug =1;
    if(debug) printf("storage in dam_check_qvs_403 : %f\n", S);

    if (!has_dam)
        return -1;

    iterations = qvs->n_values - 1;
    for (i = 0; i<iterations; i++)
    {
        if (qvs->points[i][0] <= S && S < qvs->points[i + 1][0])
            return i;

    }
    if(debug) printf("dam_check_qvs_403 iterations found: %f\n", i);

    return i;
}