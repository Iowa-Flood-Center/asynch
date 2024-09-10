#if !defined(ASYNCH_MODELS_EQUATIONS_H)
#define ASYNCH_MODELS_EQUATIONS_H


#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <structs.h>
#include <rkmethods.h>

//RHS

//Tibebu's Models
void nodam_rain_hillslope(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void dam_rain_hillslope(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void nodam_rain_hillslope2(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void dam_rain_hillslope2(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void nodam_rain_hillslope3(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void dam_rain_hillslope3(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void nodam_rain_hillslope_qsv(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void dam_rain_hillslope_qsv(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void model255(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void dam_model255(const double * const y_i, unsigned int dim, const double * const global_params, const double * const params, const QVSData * const qvs, int state, void* user, double *ans);
void dam_TopLayerHillslope_variable(const double * const y_i, unsigned int dim, const double * const global_params, const double * const params, const QVSData * const qvs, int state, void* user, double *ans);

void dam_Tetis03(const double * const y_i, unsigned int dim, const double * const global_params, const double * const params, const QVSData * const qvs, int state, void* user, double *ans);


void simple_river(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void river_rainfall(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void simple_hillslope(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void simple_soil(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void soil_rainfall(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void qsav_rainfall(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void river_rainfall_adjusted(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);

//Data Assimilation Models
void assim_simple_river(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void assim_river_rainfall(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void assim_river_rainfall_adjusted(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);

//Forecast Models
void LinearHillslope_Evap(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
int LinearHillslope_Evap_Check(const double * const y_i, const double * const params, const double * const global_params, const QVSData * const qvs, unsigned int dam);

void Hillslope_Toy(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void LinearHillslope_Evap_RC(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void LinearHillslope_MonthlyEvap(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void LinearHillslope_MonthlyEvap_OnlyRouts(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void LinearHillslope_MonthlyEvap_OnlyRouts_NotReservoir(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void LinearHillslope_MonthlyEvap_OnlyRouts_HasReservoir(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void LinearHillslope_MonthlyEvap_extras(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void LinearHillslope_MonthlyEvap_kI_extras(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void LinearHillslope_Reservoirs_extras(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void NonLinearHillslope(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void model1(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);

void model252(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void TopLayerHillslope_Reservoirs(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void TopLayerNonlinearExp(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void model249(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void model249_reservoirs(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void model251(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);

void model253(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void model254(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void model256(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void model256_reservoirs(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void model263(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void model264(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void TopLayerHillslope_spatial_velocity(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void OfflineTopLayerHillslope(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void OfflineTopLayerHillslope_Reservoirs(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void OfflineTopLayerInterflowHillslope(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void OfflineTopLayerInterflowHillslope_Reservoirs(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void TopLayerNonlinearExpSoilvel(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void TopLayerNonlinearExpSoilvel_Reservoirs(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void TopLayerNonlinearExpSoilvel_ConstEta(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void TopLayerNonlinearExpSoilvel_ConstEta_Reservoirs(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void model400(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void model401(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void Tetis03(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void Tetis03_Reservoirs(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void sac_sma01(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void routing_runoff1(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void routing_runoff2(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);

//Tiling
void Tiling(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void TilingHillslope(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);

//Nicolas Models 60X family
//601
void distributed_v1(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
//602
void distributed_v2(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
//603
void dist_v2_temperature(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
//608
void TilesModel(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void Tiles_Reservoirs(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
//609
void TilesModel_Base(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void Tiles_Reservoirs_Base(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
//610
void ActiveLayerSnow(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void ActiveLayerSnow_Reservoir(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
//612
void tetis_nicoV1(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
//654
void Variable_TopLayer(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
//Snow processes
double snow_rainfall_partition(double temp_air, double temp_thres, double temp_range);
double snow_melt_degree_day(double swe, double temp_air, double temp_thres, double melt_factor);
//Evapotranspiration equations
void evp_weighted(double *split, double s_p, double s_t, double s_s, double Tl, double Ts);
//Infiltration equations 
double vertical_qts(double s_t, double s_s, double Ts, double ks, int method);


//Stream Temperature Model
//700
void Temperature_Model(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
//temperature model functions
double solar_radiation(double hsi, double sf);
double longwave_radiation(double tair, double rh, double cloud, double twater);
double evaporative_heat_transfer(double rh, double a, double b, double twater, double w2);
double convective_heat_transfer(double pa, double tair, double twater, double a, double b, double w2);

//Misc Models
void lcuencas_soilrain(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void river_rainfall_summary(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void Robertson(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);

//Jacobians
//double NormJx_simple_river(double *y_i,double *global_params,double *param);
void Jsimple_river(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, double *ans);
void Jsimple(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, double *ans);
void Jsimple_soil(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, double *ans);
void JRobertson(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, double *ans);

//Algebraic equations
void dam_q(const double * const y_i, unsigned int num_dof, const double * const global_params, const double * const params, const QVSData * const qvs, int state, void* user, double *ans);
void dam_q2(const double * const y_i, unsigned int num_dof, const double * const global_params, const double * const params, const QVSData * const qvs, int state, void* user, double *ans);
void dam_q3(const double * const y_i, unsigned int num_dof, const double * const global_params, const double * const params, const QVSData * const qvs, int state, void* user, double *ans);
void dam_q_qvs(const double * const y_i, unsigned int num_dof, const double * const global_params, const double * const params, const QVSData * const qvs, int state, void* user, double *ans);
void dam_TopLayerNonlinearExpSoilvel(const double * const y_i, unsigned int num_dof, const double * const global_params, const double * const params, const QVSData * const qvs, int state, void* user, double *ans);
void dam_TopLayerNonlinearExpSoilvel_ConstEta(const double * const y_i, unsigned int num_dof, const double * const global_params, const double * const params, const QVSData * const qvs, int state, void* user, double *ans);

#endif //ASYNCH_MODELS_EQUATIONS_H