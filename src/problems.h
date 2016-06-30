#ifndef PROBLEMS_H
#define PROBLEMS_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "structs.h"
#include "rkmethods.h"
#include "mathmethods.h"
//#include "muParserDLL.h"

//RHS

void parser_test(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);

//Tibebu's Models
void nodam_rain_hillslope(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void dam_rain_hillslope(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void nodam_rain_hillslope2(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void dam_rain_hillslope2(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void nodam_rain_hillslope3(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void dam_rain_hillslope3(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void nodam_rain_hillslope_qsv(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void dam_rain_hillslope_qsv(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void TopLayerHillslope_variable(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void dam_TopLayerHillslope_variable(VEC y,VEC global_params,VEC params,QVSData* qvs,int state,void* user,VEC ans);

//Rodica's Models
void simple_river(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void river_rainfall(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void simple_hillslope(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void simple_soil(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void soil_rainfall(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void qsav_rainfall(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void river_rainfall_adjusted(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);

//Data Assimilation Models
void assim_simple_river(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void assim_river_rainfall(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void assim_river_rainfall_adjusted(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);

//Forecast Models
void LinearHillslope_Evap(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
int LinearHillslope_Evap_Check(VEC y,VEC params,VEC global_params,QVSData* qvs,unsigned int dam);
void Hillslope_Toy(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void LinearHillslope_Evap_RC(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void LinearHillslope_MonthlyEvap(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void LinearHillslope_MonthlyEvap_extras(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void LinearHillslope_Reservoirs_extras(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void NonLinearHillslope(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void TopLayerHillslope(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void TopLayerHillslope_Reservoirs(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void TopLayerNonlinearExp(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void TopLayerHillslope_extras(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void TopLayerNonlinearExpSoilvel(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void TopLayerNonlinearExpSoilvel_Reservoirs(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void TopLayerNonlinearExpSoilvel_ConstEta(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void TopLayerNonlinearExpSoilvel_ConstEta_Reservoirs(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);


//Tiling
void Tiling(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);

//Misc Models
void lcuencas_soilrain(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void river_rainfall_summary(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);
void Robertson(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,QVSData* qvs,VEC params,int state,void* user,VEC ans);


//Jacobians
//double NormJx_simple_river(VEC y_i,VEC global_params,VEC param);
void Jx_simple_river(VEC y_i,VEC global_params,VEC params,VEC ans);
void Jsimple(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,VEC param,MAT* ans);
void Jsimple_soil(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,VEC param,MAT* ans);
void JRobertson(double t,VEC y_i,VEC* y_p,unsigned short int numparents,VEC global_params,double* forcing_values,VEC param,MAT* ans);

//Checks
int dam_check(VEC y,VEC global_params,VEC param,QVSData* qvs,unsigned int dam);
int dam_check2(VEC y,VEC global_params,VEC param,QVSData* qvs,unsigned int dam);
int dam_check3(VEC y,VEC global_params,VEC param,QVSData* qvs,unsigned int dam);
int dam_check_qvs(VEC y,VEC global_params,VEC params,QVSData* qvs,unsigned int dam);

//Algebraic equations
void dam_q(VEC y,VEC global_params,VEC param,QVSData* qvs,int state,void* user,VEC ans);
void dam_q2(VEC y,VEC global_params,VEC param,QVSData* qvs,int state,void* user,VEC ans);
void dam_q3(VEC y,VEC global_params,VEC param,QVSData* qvs,int state,void* user,VEC ans);
void dam_q_qvs(VEC y,VEC global_params,VEC params,QVSData* qvs,int state,void* user,VEC ans);
void dam_TopLayerNonlinearExpSoilvel(VEC y,VEC global_params,VEC params,QVSData* qvs,int state,void* user,VEC ans);
void dam_TopLayerNonlinearExpSoilvel_ConstEta(VEC y,VEC global_params,VEC params,QVSData* qvs,int state,void* user,VEC ans);

//Consistency
void CheckConsistency_Nonzero_1States(VEC y,VEC params,VEC global_params);
void CheckConsistency_Nonzero_2States(VEC y,VEC params,VEC global_params);
void CheckConsistency_Nonzero_3States(VEC y,VEC params,VEC global_params);
void CheckConsistency_Nonzero_4States(VEC y,VEC params,VEC global_params);
void CheckConsistency_Model5(VEC y,VEC params,VEC global_params);
void CheckConsistency_Model30(VEC y,VEC params,VEC global_params);
void CheckConsistency_Nonzero_AllStates_q(VEC y,VEC params,VEC global_params);
void CheckConsistency_Nonzero_AllStates_qs(VEC y,VEC params,VEC global_params);

#endif

