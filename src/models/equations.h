#if !defined(ASYNCH_MODELS_EQUATIONS_H)
#define ASYNCH_MODELS_EQUATIONS_H


#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <structs.h>
#include <rkmethods.h>

//RHS

//Tibebu's Models
void nodam_rain_hillslope(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void dam_rain_hillslope(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void nodam_rain_hillslope2(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void dam_rain_hillslope2(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void nodam_rain_hillslope3(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void dam_rain_hillslope3(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void nodam_rain_hillslope_qsv(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void dam_rain_hillslope_qsv(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void TopLayerHillslope_variable(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void dam_TopLayerHillslope_variable(const double * const y_i, unsigned int dim, const double * const global_params, const double * const params, const QVSData * const qvs, int state, void* user, double *ans);

void simple_river(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void river_rainfall(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void simple_hillslope(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void simple_soil(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void soil_rainfall(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void qsav_rainfall(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void river_rainfall_adjusted(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);

//Data Assimilation Models
void assim_simple_river(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void assim_river_rainfall(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void assim_river_rainfall_adjusted(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);

//Forecast Models
void LinearHillslope_Evap(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
int LinearHillslope_Evap_Check(const double * const y_i, const double * const params, const double * const global_params, const QVSData * const qvs, unsigned int dam);

void Hillslope_Toy(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void LinearHillslope_Evap_RC(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void LinearHillslope_MonthlyEvap(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void LinearHillslope_MonthlyEvap_OnlyRouts(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void LinearHillslope_MonthlyEvap_extras(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void LinearHillslope_Reservoirs_extras(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void NonLinearHillslope(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void TopLayerHillslope(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void TopLayerHillslope_Reservoirs(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void TopLayerNonlinearExp(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void TopLayerHillslope_extras(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void TopLayerHillslope_even_more_extras(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void TopLayerNonlinearExpSoilvel(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void TopLayerNonlinearExpSoilvel_Reservoirs(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void TopLayerNonlinearExpSoilvel_ConstEta(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void TopLayerNonlinearExpSoilvel_ConstEta_Reservoirs(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);


//Tiling
void Tiling(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);

//Misc Models
void lcuencas_soilrain(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void river_rainfall_summary(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void Robertson(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);


//Jacobians
//double NormJx_simple_river(double *y_i,double *global_params,double *param);
void Jsimple_river(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, double *ans);
void Jsimple(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, double *ans);
void Jsimple_soil(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, double *ans);
void JRobertson(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, double *ans);

//Algebraic equations
void dam_q(const double * const y_i, unsigned int num_dof, const double * const global_params, const double * const params, const QVSData * const qvs, int state, void* user, double *ans);
void dam_q2(const double * const y_i, unsigned int num_dof, const double * const global_params, const double * const params, const QVSData * const qvs, int state, void* user, double *ans);
void dam_q3(const double * const y_i, unsigned int num_dof, const double * const global_params, const double * const params, const QVSData * const qvs, int state, void* user, double *ans);
void dam_q_qvs(const double * const y_i, unsigned int num_dof, const double * const global_params, const double * const params, const QVSData * const qvs, int state, void* user, double *ans);
void dam_TopLayerNonlinearExpSoilvel(const double * const y_i, unsigned int num_dof, const double * const global_params, const double * const params, const QVSData * const qvs, int state, void* user, double *ans);
void dam_TopLayerNonlinearExpSoilvel_ConstEta(const double * const y_i, unsigned int num_dof, const double * const global_params, const double * const params, const QVSData * const qvs, int state, void* user, double *ans);

#endif //ASYNCH_MODELS_EQUATIONS_H
