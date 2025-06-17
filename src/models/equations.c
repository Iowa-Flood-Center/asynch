#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#include <minmax.h>
#include <models/equations.h>


//extern int flaggy;

//These are the right-hand side functions for the ODEs. For each of these functions:
//double t: The current time
//double *y_i: The approximate value of the solution to the ODEs at the current link at time t
//double *y_p: y_p[j] has the approximate values for the immediately upstream link j to the current link at time t
//unsigned short num_parents: The number of upstream links (parents) to link i
//double *global_params: The global parameters
//RainData* rain: The rain fall values for link i
//double *params: The parameters for link i
//int state: The current state of the system
//double *ans (set by method, assumed that space is allocated): The value returned by the right-hand side function

/*
//Type 2000
//Test for equations from a parser. Same as Type 1.
//Calculates the flow using simple parameters, using only the flow q.
//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
//The numbering is:        0      1        2     3   4   5
//This uses the units and functions from September 18, 2011 document
//y_i[0] = q, y_i[1] = s
void parser_test(double t,double *y_i,double *y_p,unsigned short num_parents, unsigned int max_dim,double *global_params,double *forcing_values,QVSData *qvs,double *params,int state,void* user,double *ans)
{
unsigned int i;
unsigned int dim = y_i.dim;
double *tempresult;
double inflow = 0.0;

Formula* equations = (Formula*) qvs;	//!!!! Pass in somehow !!!!
muParserHandle_t* parser = equations->parser;
double *variable_values = equations->variable_values;
//int numvars = mupGetVarNum(parser);


//Set variables
for(i=0;i<dim;i++)	variable_values.ve[i] = y_i[i];	//States
variable_values.ve[dim] = t;					//Time
for(i=0;i<num_parents;i++)	inflow += y_p[i * dim];
variable_values.ve[dim+1] = inflow;				//Inflow
variable_values.ve[dim+2] = forcing_values[0];			//Rainfall

//Evaluate equations
tempresult = mupEvalMulti(parser,&i);
for(i=0;i<dim;i++)	ans[i] = tempresult[i];
}
*/


double sq(double x) { return x * x; }

//Type 1 (UID 1)
void model1(double t, \
		const double * const y_i, \
		unsigned int dim, \
		const double * const y_p, \
		unsigned short num_parents, \
		unsigned int max_dim, \
		const double * const global_params, \
		const double * const params, \
		const double * const forcing_values, \
		const QVSData * const qvs, \
		int state, \
		void* user, \
		double *ans)
{
    unsigned short i;

    double lambda_1 = global_params[1];
    double k_3 = global_params[4];  //[1/min]
    double h_b = global_params[6];  //[m]
    double S_L = global_params[7];  //[m]
    double A = global_params[8];
    double B = global_params[9];
    double exponent = global_params[10];
    double v_B = global_params[11];
    double k_tl = global_params[12];
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));   //[mm/month] -> [m/min]

    double L = params[1];   //[m]
    double A_h = params[2]; //[m^2]
    double invtau = params[3];  //[1/min]
    double k_2 = params[4];     //[1/min]
    double k_i = params[5];     //[1/min]
    double c_1 = params[6];
    double c_2 = params[7];

    double q = y_i[0];		//[m^3/s]
    double S = y_i[1];		//[m^3]
    double s_p = y_i[2];	//[m]
    double s_t = y_i[3];	//[m]
    double s_s = y_i[4];	//[m]
//    double q_b = y_i[7];    //[m^3/s]


    //Evaporation
    double e_p, e_t, e_s;
    double Corr = s_p + s_t / S_L + s_s / (h_b - S_L);
    if (e_pot > 0.0 && Corr > 1e-12)
    {
        e_p = s_p * e_pot / Corr;
        e_t = s_t / S_L * e_pot / Corr;
        e_s = s_s / (h_b - S_L) * e_pot / Corr;
    }
    else
    {
        e_p = 0.0;
        e_t = 0.0;
        e_s = 0.0;
    }

    double pow_term = (1.0 - s_t / S_L > 0.0) ? pow(1.0 - s_t / S_L, exponent) : 0.0;
    double k_t = (A + B * pow_term) * k_2;

    //Fluxes
    double q_pl = k_2 * s_p;
    double q_pt = k_t * s_p;
    double q_ts = k_i * s_t;
    double q_tl = k_tl * s_t;
    double q_sl = k_3 * s_s;    //[m/min]

    //Discharge
    ans[0] = -q + (q_pl + q_tl + q_sl) * c_2;
    for (i = 0; i < num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];    // discharge[0]


    //Discharge
    dam_TopLayerHillslope_variable(y_i, dim, global_params, params, qvs, state, user, ans);	//ans is used for convenience !!!!
    // Is q available in y_i? !!!!
    double qm = ans[0] * 60.0;

    //Storage
    ans[1] = (q_pl + q_sl) * A_h - qm;
    for (i = 0; i<num_parents; i++)
        ans[1] += y_p[i * dim] * 60.0;

    //Hillslope
    ans[2] = forcing_values[0] * c_1 - q_pl - q_pt - e_p; // pond[1]
    ans[3] = q_pt - q_ts - q_tl - e_t;  // toplayer[2]
    ans[4] = q_ts - q_sl - e_s;  // subsurface[3]
}

//Type 601
//This is the first model from a set of models with all parameters taken from the disk not from the globals,
// I will start with a modifiction of the 225 Tiling model of morgan.
//The variables are river flow, ponded, top layer, subsurface layer.
//Order of parameters: A_i,L_i,A_h,slope,tileQ,invtau,aNo,actT,LinSl,aTile,LinIn,
//The numbering is:	0   1   2    3     4     5     6   7     8   9     10
//Order of global_params: v_r,lambda_1,lambda_2,tL,bL,kdry,ki,k2
//The numbering is:        0      1        2     3 4   5   6   7 
void VariableHillslope(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i; 
    //Distributed variables
    double A_i = params[0];
    double L_i = params[1];
    double A_h = params[2];
    double slope = params[3];
    double v_r = params[4];
    double a_r = params[5];
    double v_s = params[6];
    double a_s = params[7];
    double k1 = params[8];
    double k2 = params[9];
    double t_L = params[10];
    double b_L = params[11];
    double lambda_1 = params[12];
    double invtau = params[13];
    //Variables
    double q = y_i[0];		                                        // [m^3/s]
    double s_p = y_i[1];	                                        // [m]
    double s_l = y_i[2];	                                        // [m]
    double s_s = y_i[3];
    //Fluxes
    double q_in = forcing_values[0] * (0.001/60);	//[m/min]
    double q_pl = k1*pow((1-s_l/t_L),3)*s_p;                       // (m/min)
    double q_ls = k2*s_l;			      // (m/min)
    double q_pLink = v_r*pow(s_p,a_r);                                   // (m/min)
    double q_out = v_s*pow(s_s,a_s);
    //Evaporation
    double C_p = s_p;
    double C_l = s_l/t_L;
    double C_s = s_s/b_L;
    double Corr_evap = 1/(C_p + C_l + C_s);
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]
    double e_p = Corr_evap * C_p * e_pot;
    double e_l = Corr_evap * C_l * e_pot;
    double e_s = Corr_evap * C_s * e_pot;
    //Update variables
	double q_parent;
	int q_pidx;
    //Discharge
    ans[0] = -q + ((q_pLink + q_out) * A_h / 60.0);
	for (i = 0; i < num_parents; i++) {
		q_pidx = i * dim;
		q_parent = y_p[q_pidx];
		ans[0] += q_parent;
	}
    ans[0] = invtau * pow(q, lambda_1) * ans[0];
    //Ponded
    ans[1] = q_in - q_pl - q_pLink - e_p;
    //Top Soil Layer
    ans[2] = q_pl - q_ls - e_l;	
    //Subsurface (saturated) soil
    ans[3] = q_ls - q_out - e_s;
}


//Type 602
//model with interflow, in this case we consider that there is also an interflow
//from the toplayer
//This is the first model from a set of models with all parameters taken from the disk not from the globals, I will starrt with a modifiction of the 225 Tiling model of morgan. 
//The variables are river flow, ponded, top layer, subsurface layer.
//Order of parameters: A_i,L_i,A_h,slope,tileQ,invtau,aNo,actT,LinSl,aTile,LinIn,
//The numbering is:	0   1   2    3     4     5     6   7     8   9     10
//Order of global_params: v_r,lambda_1,lambda_2,tL,bL,kdry,ki,k2
//The numbering is:        0      1        2     3 4   5   6   7 
void VariableInterflow(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i; 
    //Distributed variables
    double A_i = params[0];
    double L_i = params[1];
    double A_h = params[2];
    double slope = params[3];
    double v_r = params[4];
    double a_r = params[5]; 
    double v_l = params[6];
    double a_l = params[7];
    double v_s = params[8];
    double a_s = params[9];
    double k1 = params[10];
    double k2 = params[11];
    double t_L = params[12];
    double b_L = params[13];
    double lambda_1 = params[14];
    double invtau = params[15];
    //Variables
    double q = y_i[0];		                                        // [m^3/s]
    double s_p = y_i[1];	                                        // [m]
    double s_l = y_i[2];	                                        // [m]
    double s_s = y_i[3];
    //Fluxes
    double q_in = forcing_values[0] * (0.001/60);	//[m/min]
    double q_pl = k1*pow((1-s_l/t_L),3)*s_p;                       // (m/min)
    double q_ls = k2*s_l;			      // (m/min)
    double q_pLink = v_r*pow(s_p,a_r);                                   // (m/min)
    double q_lLink = v_l*pow(s_l,a_l);                                   // (m/min) 
    //double C1 = v_s*pow(s_s,a_s);
    //double C2 = 2.04e-6*s_s; 
    double q_sLink = v_s*pow(s_s,a_s);
    //Evaporation
    double C_p = s_p;
    double C_l = s_l/t_L;
    double C_s = s_s/b_L;
    double Corr_evap = 1/(C_p + C_l + C_s);
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]
    double e_p = Corr_evap * C_p * e_pot;
    double e_l = Corr_evap * C_l * e_pot;
    double e_s = Corr_evap * C_s * e_pot;
    //Update variables
	double q_parent;
	int q_pidx;
    //Discharge
    ans[0] = -q + ((q_pLink + q_lLink + q_sLink) * A_h / 60.0);
	for (i = 0; i < num_parents; i++) {
		q_pidx = i * dim;
		q_parent = y_p[q_pidx];
		ans[0] += q_parent;
	}
    ans[0] = invtau * pow(q, lambda_1) * ans[0];
    //Ponded
    ans[1] = q_in - q_pl - q_pLink - e_p;
    //Top Soil Layer
    ans[2] = q_pl - q_ls - q_lLink  - e_l;	
    //Subsurface (saturated) soil
    ans[3] = q_ls - q_sLink - e_s;
}

//Type 603
//VariableTriLayer: model with variable parameters and three layers of soil, the
//model has a total of five states: 0: streamflow, 1: runoff, 2 to 4 soil layers
//The variables in the param file (.prm) are:
// Ai, Li, Ah, So, v1, a1, v2, a2, v3, a3, v4, a4, h1, h2, h3, k1, k2, k3, lam1, lam2, v0.
void VariableTriLayer(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i; 
    //Distributed variables
    double Ai = params[0];
    double Li = params[1];
    double Ah = params[2];
    double slope = params[3];
    double v1 = params[4];
    double a1 = params[5]; 
    double v2 = params[6];
    double a2 = params[7];
    double v3 = params[8];
    double a3 = params[9];
    double v4 = params[10];
    double a4 = params[11];
    double h1 = params[12];
    double h2 = params[13];
    double h3 = params[14];
    double k1 = params[15];
    double k2 = params[16];
    double k3 = params[17];
    double lambda1 = params[18];
    double invtau = params[19];
    //Variables
    double q = y_i[0];		                                        // [m^3/s]
    double sp = y_i[1];	                                        // [m]
    double s1 = y_i[2];	                                        // [m]
    double s2 = y_i[3];
    double s3 = y_i[4];
    //Vertical fluxes
    double qin = forcing_values[0] * (0.001/60);	//[m/min]
    double qv1 = k1*pow((1-s1/h1),3)*sp;                       // (m/min)
    double qv2 = k2*s2;			      // (m/min)
    double qv3 = k3*s3;
    //Horizontal fluxes
    double qh1 = v1*pow(sp,a1);
    double qh2 = v2*pow(s1,a2);
    double qh3 = v3*pow(s2,a3);
    double qh4 = v4*pow(s3,a4);
    //Evaporation
    double Cp = sp;
    double C1 = s1/h1;
    double C2 = s2/h2;
    double C3 = s3/h3;
    double Corr_evap = 1/(Cp + C1 + C2 + C3);
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]
    double ep = Corr_evap * Cp * e_pot;
    double e1 = Corr_evap * C1 * e_pot;
    double e2 = Corr_evap * C2 * e_pot;
    double e3 = Corr_evap * C3 * e_pot;
    //Update variables
	double q_parent;
	int q_pidx;
    //Discharge
    ans[0] = -q + ((qh1+qh2+qh3) * Ah / 60.0);
	for (i = 0; i < num_parents; i++) {
		q_pidx = i * dim;
		q_parent = y_p[q_pidx];
		ans[0] += q_parent;
	}
    ans[0] = invtau * pow(q, lambda1) * ans[0];
    //Ponded
    ans[1] = qin - qv1 - qh1 - ep;
    //Soil layers are almost the same
    ans[2] = qv1 - qv2 - qh2 - e1;
    ans[3] = qv2 - qv3 - qh3 - e2;
    ans[4] = qv3 - qh4 - e3;
}


//Type 604
//VariableTreshold: model similar to 602 but with a change of subsurface runoff
//in function of the soil water storage.
//model has a total of five states: 0: streamflow 1: runoff 2: top layer 3:
//subsurface
//The variables in the param file (.prm) are:
// Ai, Li, Ah, So, v1, a1, v2, a2, v3, a3,   h1, h2, h3, k1, k2, k3, lam1, lam2, v0.
void VariableThreshold(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i; 
    //Distributed variables
    double A_i = params[0];
    double L_i = params[1];
    double A_h = params[2];
    double v_r = params[3];
    double a_r = params[4];
    double v_s1 = params[5];
    double v_s2 = params[6];
    double k1 = params[7];
    double k2 = params[8];
    double t_L = params[9];
    double Si = params[10];
    double Sn = params[11];
    double lambda_1 = params[12];
    double invtau = params[13];
    double Is = params[14];
    //Variables
    double q = y_i[0];		                                        // [m^3/s]
    double s_p = y_i[1];	                                        // [m]
    double s_l = y_i[2];	                                        // [m]
    double s_s = y_i[3];
    double q_s = y_i[4];
    //Fluxes
    double q_in = forcing_values[0] * (0.001/60);	//[m/min]
    double q_pl = k1*pow((1-s_l/t_L),3)*s_p;                       // (m/min)
    double q_ls = k2*s_l;			      // (m/min)
    double q_pLink = v_r*pow(s_p,a_r);                                   // (m/min)
    double q_sLink = 0.0;
    if (s_s < Si){
        q_sLink = v_s1*(s_s - Sn);
        ans[4] = 0.0;
    }
    else{
        q_sLink = v_s2*(s_s - Si) + Is;
        ans[4] = q_sLink;
    }
    ans[5] = q_sLink;
    //ans[4] = q_sLink;
    //Evaporation
    double C_p = s_p;
    double C_l = s_l/t_L;
    double C_s = s_s/Si;
    double Corr_evap = 1/(C_p + C_l + C_s);
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]
    double e_p = Corr_evap * C_p * e_pot;
    double e_l = Corr_evap * C_l * e_pot;
    double e_s = Corr_evap * C_s * e_pot;
    //Update variables
	double q_parent;
	int q_pidx;
    //Discharge
    ans[0] = -q + ((q_pLink + q_sLink) * A_h / 60.0);
	for (i = 0; i < num_parents; i++) {
		q_pidx = i * dim;
		q_parent = y_p[q_pidx];
		ans[0] += q_parent;
	}
    ans[0] = invtau * pow(q, lambda_1) * ans[0];
    //Ponded
    ans[1] = q_in - q_pl - q_pLink - e_p;
    //Top Soil Layer
    ans[2] = q_pl - q_ls - e_l;	
    //Subsurface (saturated) soil
    ans[3] = q_ls - q_sLink - e_s;
    //Record the total rainfall in the run
    ans[6] = q_in;
}

//Type 605
//VariableTreshold3: similar to 604, in this case there are 3 different slopes
//in function of the soil water storage.
//model has a total of five states: 0: streamflow 1: runoff 2: top layer 3:
//subsurface
//The variables in the param file (.prm) are:
// Ai, Li, Ah, So, v1, a1, v2, a2, v3, a3,   h1, h2, h3, k1, k2, k3, lam1, lam2, v0.
void VariableThreshold3(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i; 
    //Distributed variables
    double A_i = params[0];
    double L_i = params[1];
    double A_h = params[2];
    double v_r = params[18];
    double a_r = params[4];
    double v_s1 = params[5];
    double v_s2 = params[6];
    double v_s3 = params[7];
    double v_s4 = params[8];
    double k1 = params[9];
    double k2 = params[10];
    double t_L = params[11];
    double S1 = params[12];
    double S2 = params[13];
    double S3 = params[14];
    double S4 = params[15];
    double lambda_1 = params[16];
    double invtau = params[17];
    double I1 = params[19];
    double I2 = params[20];
    double I3 = params[21];
    //Variables
    double q = y_i[0];		                                        // [m^3/s]
    double s_p = y_i[1];	                                        // [m]
    double s_l = y_i[2];	                                        // [m]
    double s_s = y_i[3];
    double q_s = y_i[4];
    //Fluxes
    double q_in = forcing_values[0] * (0.001/60);	//[m/min]
    double q_pl = k1*pow((1-s_l/t_L),3)*s_p;                       // (m/min)
    double q_ls = k2*s_l;			      // (m/min)
    double q_pLink = v_r*pow(s_p,a_r);                                   // (m/min)
    //subsurface runoff
    double q_sLink = 0.0;
    if (s_s < S2){
        q_sLink = v_s1*(s_s - S1);
        ans[4] = 0.0;
    }
    else if (s_s >= S2 && s_s < S3 ){
        q_sLink = v_s2*(s_s - S2) + I1;
        ans[4] = q_sLink;
    }
    else if (s_s >= S3 && s_s < S4 ){
        q_sLink = v_s3*(s_s - S3) + I2;
        ans[4] = q_sLink;
    }
    else{
        q_sLink = v_s4*(s_s - S4) + I3;
        ans[4] = q_sLink;
    }
    ans[5] = q_sLink;
    
    //Evaporation
    double C_p = s_p;
    double C_l = s_l/t_L;
    double C_s = s_s/S2;
    double Corr_evap = 1/(C_p + C_l + C_s);
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]
    double e_p = Corr_evap * C_p * e_pot;
    double e_l = Corr_evap * C_l * e_pot;
    double e_s = Corr_evap * C_s * e_pot;
    //Update variables
	double q_parent;
	int q_pidx;
    //Discharge
    ans[0] = -q + ((q_pLink + q_sLink) * A_h / 60.0);
	for (i = 0; i < num_parents; i++) {
		q_pidx = i * dim;
		q_parent = y_p[q_pidx];
		ans[0] += q_parent;
	}
    ans[0] = invtau * pow(q, lambda_1) * ans[0];
    //Ponded
    ans[1] = q_in - q_pl - q_pLink - e_p;
    //Top Soil Layer
    ans[2] = q_pl - q_ls - e_l;	
    //Subsurface (saturated) soil
    ans[3] = q_ls - q_sLink - e_s;
    //Record the total rainfall in the run
    ans[6] = q_in;
}


//Type 606
//VariableTreshold3: similar to 604, in this case there are 3 different slopes
//in function of the soil water storage.
//model has a total of five states: 0: streamflow 1: runoff 2: top layer 3:
//subsurface
//The variables in the param file (.prm) are:
// Ai, Li, Ah, So, v1, a1, v2, a2, v3, a3,   h1, h2, h3, k1, k2, k3, lam1, lam2, v0.
void ExponentialExp(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i; 
    //Distributed variables
    double A_i = params[0];
    double L_i = params[1];
    double A_h = params[2];
    double v_r = params[17];
    double a_r = params[4];
    double a = params[5];
    double b = params[6];
    double m = params[7];
    double mt = params[8];
    double k1 = params[9];
    double k2 = params[10];
    double t_L = params[11];
    double NoFlow = params[12];
    double Td = params[13];
    double Beta = params[14];
    double lambda_1 = params[15];
    double invtau = params[16];
    //Variables
    double q = y_i[0];		                                        // [m^3/s]
    double s_p = y_i[1];	                                        // [m]
    double s_l = y_i[2];	                                        // [m]
    double s_s = y_i[3];
    double s_t = y_i[4];
    //Fluxes
    double q_in = forcing_values[0] * (0.001/60);	//[m/min]
    double q_pl = k1*pow((1-s_l/t_L),3)*s_p;                       // (m/min)
    double q_ls = k2*s_l;			      // (m/min)
    double q_pLink = v_r*pow(s_p,a_r);                                   // (m/min)
    //subsurface runoff
    double q_sLtemp = 0.0;  // Variable to record the amount of tile and exponential flows.
    double q_sLink = 0.0;
    double q_inT = 0.0;
    double q_outT = 0.0;
    if (s_s > NoFlow){
        q_sLtemp = m * (s_s - NoFlow); // Base flow or linear portion
        q_sLink += q_sLtemp;
        //ans[4] = q_sLtemp; // Total base flow
    }
    if (s_s > Td){
        q_inT = mt * (s_s - Td); // Tile flow in function of the tile act depth and tile slopei
        ans[6] += q_inT;
    }    
    if (s_t > 0.0){
        q_outT = 2500*v_r * pow(s_t, a_r); // Tile bring water regardless of the level of the subsurface level
        //q_outT = 0.9*s_t;
        ans[5] = q_outT;
    }
    if (s_s>Beta){
        q_sLtemp = (s_s - Beta) * a * exp(b * (s_s - Beta)); // Active runoff explained by an exponential func
        q_sLink += q_sLtemp;
        //ans[6] = q_sLtemp;  // Total active outflow
    }  
    ans[6] = q_sLink;
    //Evaporation
    double C_p = s_p;
    double C_l = s_l/t_L;
    double C_s = s_s/Beta;
    double Corr_evap = 1/(C_p + C_l + C_s);
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]
    double e_p = Corr_evap * C_p * e_pot;
    double e_l = Corr_evap * C_l * e_pot;
    double e_s = Corr_evap * C_s * e_pot;
    //Update variables
	double q_parent;
	int q_pidx;
    //Discharge
    ans[0] = -q + ((q_pLink + q_sLink + q_outT) * A_h / 60.0);
	for (i = 0; i < num_parents; i++) {
		q_pidx = i * dim;
		q_parent = y_p[q_pidx];
		ans[0] += q_parent;
	}
    ans[0] = invtau * pow(q, lambda_1) * ans[0];
    //Ponded
    ans[1] = q_in - q_pl - q_pLink - e_p;
    //Top Soil Layer
    ans[2] = q_pl - q_ls - e_l;	
    //Subsurface (saturated) soil
    ans[3] = q_ls - q_sLink - e_s - q_inT;
    //Tile storage
    ans[4] = q_inT - q_outT;
    
    //Record the total rainfall in the run
    //ans[7] = q_in;
}


//Type 608
//VariableTreshold3: similar to 604, in this case there are 3 different slopes
//in function of the soil water storage.
//model has a total of five states: 0: streamflow 1: runoff 2: top layer 3:
//subsurface
//The variables in the param file (.prm) are:
// Ai, Li, Ah, So, v1, a1, v2, a2, v3, a3,   h1, h2, h3, k1, k2, k3, lam1, lam2, v0.
void Tiles_Reservoirs(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    if(forcing_values[2] > 0){
        ans[0] = forcing_values[2];
        //Copy Discharge open loop
        double invtau = params[16];
        double lambda_1 = params[15];
        double q_openloop = y_i[5];	
        //ans[5] = -q_openloop + ((q_pLink + q_sLink) * A_h / 60.0);
        unsigned short i;
        for (i =0; i<num_parents; i++)
            ans[5] += y_p[i*dim+5];
        
        ans[5] = invtau * pow(q_openloop, lambda_1) * ans[5];

    }
    if(forcing_values[2] <=0){
        unsigned short i;
        for (i =0; i<num_parents; i++)
            ans[0] += y_p[i*dim];
    }
    double Beta = params[14];
    ans[1] = 0.0;
    ans[2] = 0.0;
    ans[3] = Beta;
    //ans[4] = 0.0;
    //ans[5] = 0.0;
}

void TilesModel(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i; 
    //Distributed variables
    double A_i = params[0];
    double L_i = params[1];
    double A_h = params[2];
    double v_r = params[3];
    double a_r = params[4];
    double a = params[5];
    double b = params[6];
    double c = params[7];
    double d = params[8];
    double k3 = params[9];
    double ki_fac = params[10];
    double t_L = params[11];
    double NoFlow = params[12];
    double Td = params[13];
    double Beta = params[14];
    double lambda_1 = params[15];
    // Processed parameters
    double invtau = params[16];
    double k2 = params[17];
    //double expo = params[18];
    //Variables or sttates
    double q = y_i[0];	
    //double q_openloop = y_i[5];		                                        // [m^3/s]
    double s_p = y_i[1];	                                        // [m]
    double s_l = y_i[2];	                                        // [m]
    double s_s = y_i[3];

    //In forcings rainfall and snowmelt
    double q_in = forcing_values[0] * (0.001/60);	//[m/min]
    double snowmelt = forcing_values[3] * (0.001)/60; //mm/hour to m/min
    q_in += snowmelt;

    //Vertical fluxes
    double pow_t = (1.0 - s_l/t_L > 0.0)? pow(1.0 - s_l/t_L,3): 0.0;
    //double pow_t2 = (4.7 - 3.4*(s_s/1.67) > 0.0)? pow(4.7 - 3.4*(s_s/1.67),0.405): 0.0; // Exp kind of Green y Ampt 
    double q_pl = k2*99.0*pow_t*s_p;
    double q_ls = k2*ki_fac*s_l;
    //double q_ls = k2*ki_fac*pow_t2*s_l; //Exp Green y Ampt approach
    double q_pLink = k2*s_p;
    //subsurface runoff
    double q_sLink = 0.0;
    double q_inT = 0.0;
    //Base flow
    if (s_s > NoFlow){
        q_sLink += k3 * (s_s - NoFlow);                          // Base flow or linear portion
    } 
    //Active flow
    if (s_s>Beta){
        q_sLink += (s_s - Beta) * a * exp(b * (s_s - Beta));    // Active runoff explained by an exponential func
    }  
    //Tile flow
    if (s_s > Td){
        q_inT =  (s_s - Td) * c * exp(d *  (s_s - Td));
        q_sLink += q_inT;         // Tile flow in function of the tile act depth and tile slopei 
    }    
    //Evaporation
    double C_p = s_p;
    double C_l = s_l/t_L;
    double C_s = s_s/(Beta-NoFlow);
    double Corr_evap = 1/(C_p + C_l + C_s);
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]
    double e_p = Corr_evap * C_p * e_pot;
    double e_l = Corr_evap * C_l * e_pot;
    double e_s = Corr_evap * C_s * e_pot;
    //Update variables
	double q_parent;
	int q_pidx;
    //Discharge
    ans[0] = -q + ((q_pLink + q_sLink) * A_h / 60.0);	
    for (i = 0; i < num_parents; i++) {
		  q_pidx = i * dim;
		  q_parent = y_p[q_pidx];
		  ans[0] += q_parent;
	  }
    ans[0] = invtau * pow(q, lambda_1) * ans[0];
    //Ponded
    ans[1] = q_in - q_pl - q_pLink - e_p;
    //Top Soil Layer
    ans[2] = q_pl - q_ls - e_l;	
    //Subsurface (saturated) soil
    ans[3] = q_ls - q_sLink - e_s;
    //Tile storage
    //ans[6] = q_in;
    //ans[7] = q_pLink;
    //Record the total rainfall in the run
    //ans[7] = q_in;
}

//Type 609
//TilesModelBaseSep: Model 608 with baseflow separation
// Ai, Li, Ah, So, v1, a1, v2, a2, v3, a3,   h1, h2, h3, k1, k2, k3, lam1, lam2, v0.
void Tiles_Reservoirs_Base(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    if(forcing_values[2] > 0){
        ans[0] = forcing_values[2]; // Streamflow
    }
    if(forcing_values[2] <=0){
        unsigned short i;
        for (i =0; i<num_parents; i++)
            ans[0] += y_p[i*dim];
    }
    double Beta = params[14];
    ans[1] = 0.0; 
    ans[2] = 0.0;
    ans[3] = Beta; 
    ans[4] = 0.0; // Baseflow
    //ans[5] = 0.0; // OpenLoop streamflow

}

//Type 609
//TilesModelBase Model 608 with baseflow separation
void TilesModel_Base(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i; 
    //Distributed variables
    double A_i = params[0];
    double L_i = params[1];
    double A_h = params[2];
    double v_r = params[3];
    double a_r = params[4];
    double a = params[5];
    double b = params[6];
    double c = params[7];
    double d = params[8];
    double k3 = params[9];
    double ki_fac = params[10];
    double t_L = params[11]; //* forcing_values[3];
    double NoFlow = params[12];
    double Td = params[13];
    double Beta = params[14];
    double lambda_1 = params[15];
    // Processed parameters
    double invtau = params[16];
    double k2 = params[17];
    //Variables or sttates
    double q = y_i[0];		                                        // [m^3/s]
    double s_p = y_i[1];	                                        // [m]
    double s_l = y_i[2];	                                        // [m]
    double s_s = y_i[3];

    //double s_c = y_i[4];
    double q_b = (1.0e-7>y_i[4])? 1.0e-7: y_i[4];                                // for base flow separation
    //double q_b = y_i[4];


    //Fluxes
    double q_in = forcing_values[0] * (0.001/60);	//[m/min]
    //Crop (experimental)
    //double crop = forcing_values[3] * (1e-3 / (30.0*24.0*60.0)); //[mm/month] -> [m/min]
    //double q_cp = (s_c + q_in - crop > 0.0)? s_c + q_in - crop: 0.0; // Crop acting as a bucket
    //double q_cp = (crop > 0.0)? q_in*(1-s_c/crop): q_in

    double pow_t = (1.0 - s_l/t_L > 0.0)? pow(1.0 - s_l/t_L,3): 0.0;
    double q_pl = k2*99.0*pow_t*s_p;
    double q_ls = k2*ki_fac*s_l;
    double q_pLink = k2*s_p;
    //subsurface runoff
    double q_sLink = 0.0;
    double q_inT = 0.0;
    //Base flow
    if (s_s > NoFlow){
        q_sLink += k3 * (s_s - NoFlow);                          // Base flow or linear portion
    } 
    //Active flow
    if (s_s>Beta){
        q_sLink += (s_s - Beta) * a * exp(b * (s_s - Beta));    // Active runoff explained by an exponential func
    }  
    //Tile flow
    if (s_s > Td){
        q_inT =  (s_s - Td) * c * exp(d *  (s_s - Td));
        q_sLink += q_inT;         // Tile flow in function of the tile act depth and tile slopei 
    }    
    //Evaporation
    double C_p = s_p;
    double C_l = s_l/t_L;
    double C_s = s_s/(Beta-NoFlow);
    //double C_c = (crop > 0.0)? s_c: 0.0;
    double Corr_evap = 1/(C_p + C_l + C_s);
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]
    double e_p = Corr_evap * C_p * e_pot;
    double e_l = Corr_evap * C_l * e_pot;
    double e_s = Corr_evap * C_s * e_pot;
    //double e_c = Corr_evap * C_c * e_pot;
    //Update variables
	double q_parent;
	int q_pidx;
    //Discharge
    ans[0] = -q + ((q_pLink + q_sLink) * A_h / 60.0);
	ans[4] = -q_b + ((q_sLink) * A_h / 60.0);
    for (i = 0; i < num_parents; i++) {
		q_pidx = i * dim;
		q_parent = y_p[q_pidx];
		ans[0] += q_parent;
        q_parent = y_p[q_pidx+4];
		ans[4] += q_parent;
	}

    ans[0] = invtau * pow(q, lambda_1) * ans[0]; //
    ans[4] = invtau * pow(q, lambda_1) * ans[4];

    //ans[4] = (q_b/q)*ans[4];
    //Crops
    //ans[4] = q_in - q_cp - e_c;
    //Ponded
    ans[1] = q_in - q_pl - q_pLink - e_p;
    //Top Soil Layer
    ans[2] = q_pl - q_ls - e_l;	
    //Subsurface (saturated) soil
    ans[3] = q_ls - q_sLink - e_s;
}

//Type 654
//Variable_TopLayer: This is a version of the TopLayer model 254 with the
//difference that in this case all variables could be distributed 
//model has a total of four states: 0: streamflow 1: runoff 2: top layer 3:
//subsurface
//The variables in the param file (.prm) are:
// Ai, Li, Ah, So, v1, a1, v2, a2, v3, a3,   h1, h2, h3, k1, k2, k3, lam1, lam2, v0.
void Variable_TopLayer(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;
    //double A_i = params[0];
    double A_i = params[0];
    double L_i = params[1];
    double A_h = params[2];
    double v_h = params[3];
    double k_i_f = params[4];
    double k_3 = params[5];	//[1/min]
    double h_b = params[6];	//[m]
    double S_L = params[7];	//[m]
    double A = params[8];
    double B = params[9];
    double exponent = params[10];
    double lambda_1 = params[11];
    double lambda_2 = params[12];
    double v_0 = params[13];    
    //double invtau = params[12];

    //Precalculus    
    double invtau = 60.0*v_0*pow(A_i, lambda_2)/((1.0 - lambda_1)*L_i);
    double k_2 = v_h * L_i / A_h * 60.0;
    double k_i = k_2 * k_i_f;
    double c_2 = A_h / 60.0;
    double c_1 = 0.001/60.;

    //Initial conditions
    double q = y_i[0];		//[m^3/s]
    double s_p = y_i[1];	//[m]
    double s_t = y_i[2];	//[m]
    double s_s = y_i[3];	//[m]

    //Evaporation
    double e_p, e_t, e_s; 
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]
    double Corr = s_p + s_t / S_L + s_s / (h_b - S_L);
    if (e_pot > 0.0 && Corr > 1e-12)
    {
        e_p = s_p * e_pot / Corr;
        e_t = s_t / S_L * e_pot / Corr;
        e_s = s_s / (h_b - S_L) * e_pot / Corr;
    }
    else
    {
        e_p = 0.0;
        e_t = 0.0;
        e_s = 0.0;
    }

   // if (e_p > s_p)
    //{
     //   e_p = s_p * 0.99;
    //}

    //Fluxes
    double pow_term = (1.0 - s_t / S_L > 0.0) ? pow(1.0 - s_t / S_L, exponent) : 0.0;
    double k_t = (A + B * pow_term) * k_2;
    //double q_pl = k_2 * s_p;
    double q_pt = k_t * s_p;
    double q_ts = k_i * s_t;
    double q_sl = k_3 * s_s;	//[m/min]

    double q_pl = k_2*pow(s_p,1.67);
                                //Discharge
    ans[0] = -q + (q_pl + q_sl) * c_2;
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];

    //Hillslope
    ans[1] = forcing_values[0] * c_1 - q_pl - q_pt - e_p;
    ans[2] = q_pt - q_ts - e_t;
    ans[3] = q_ts - q_sl - e_s;

}




//Type 225
//This model is used for the tiling experiments. The rainfall will run through the hillslope and requires an if/then statement. The variables are river flow, ponded, top layer, subsurface layer
//Order of parameters: A_i,L_i,A_h,slope,tileQ,invtau,aNo,actT,LinSl,aTile,LinIn,
//The numbering is:	0   1   2    3     4     5     6   7     8   9     10
//Order of global_params: v_r,lambda_1,lambda_2,tL,bL,kdry,ki,k2
//The numbering is:        0      1        2     3 4   5   6   7 
void TilingHillslope(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;
    double lambda_1 = global_params[1];
    double tL = global_params[3];
    double bL = global_params[4];
    double kdry = global_params[5];
    double ki = global_params[6];
    double k2 = global_params[7];
    double A_h = params[2];
    double tileQ = params[4];
    double invtau = params[5];
    double aNo = params[6];
    double actTile = params[7];
    double linSl = params[8];
    double aTile = params[9];
    double linInt = params[10];
    //Variables
    double q = y_i[0];		                                        // [m^3/s]
    double s_p = y_i[1];	                                        // [m]
    double s_l = y_i[2];	                                        // [m]
    double s_s = y_i[3];
    //Fluxes
    double q_in = forcing_values[0] * (0.001/60);	//[m/min]
    double q_pl = kdry*pow((1-s_l/tL),3)*s_p;                       // (m/min)
    double q_ls = ki*s_l;			      // (m/min)
    double q_pLink = k2*pow(s_p,1.67);                                   // (m/min)
    double q_out = 0.0;
    if (s_s > actTile){
		if (tileQ > 0.5){
            q_out = aTile * exp(21 * s_s);
		} 
        else{
			q_out = aNo * s_s * exp(10 * s_s);
		}
	} 
    else{
        q_out = linSl * s_s + linInt;
        //if (q_out < 0){
         //   q_out = abs(linInt;
       //}
	}
   // q_out = k2*0.015*pow(s_s, 2);
    
    //Evaporation
    double C_p = s_p;
    double C_l = s_l/tL;
    double C_s = s_s/bL;
    double Corr_evap = 1/(C_p + C_l + C_s);
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]

    double e_p = Corr_evap * C_p * e_pot;
    double e_l = Corr_evap * C_l * e_pot;
    double e_s = Corr_evap * C_s * e_pot;
//    double e_s = 0.0

	double q_parent;
	int q_pidx;

    //Discharge
    ans[0] = -q + ((q_pLink + q_out) * A_h / 60.0);
	for (i = 0; i < num_parents; i++) {
		q_pidx = i * dim;
		q_parent = y_p[q_pidx];
		ans[0] += q_parent;
	}
    ans[0] = invtau * pow(q, lambda_1) * ans[0];
    //Ponded
    ans[1] = q_in - q_pl - q_pLink - e_p;
    //Top Soil Layer
    ans[2] = q_pl - q_ls - e_l;	
    //Subsurface (saturated) soil
    ans[3] = q_ls - q_out - e_s;
}

//Type 250
//Order of parameters: A_i,L_i,A_h,h_r,invtau,k_2,k_I,c_1,c_2
//The numbering is:	0   1   2   3    4     5   6   7   8
//E
//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,gamma,h_b,e_pot
//The numbering is:        0      1        2     3   4     5         6    7	8
void NonLinearHillslope(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    double lambda_1 = global_params[1];
    double k_3 = global_params[4];	//[1/min]
    double gamma = global_params[6];
    double h_b = global_params[7];	//[m]
                                        //double e_pot = global_params[8] * (1e-3/60.0);	//[mm/hr]->[m/min]
    double e_pot = global_params[8];	//[m/min]


    double L = params[1];	//[m]
    double A_h = params[2];	//[m^2]
    double h_r = params[3];	//[m]
    double invtau = params[4];	//[1/min]
    double k_2 = params[5];	//[1/min]
    double k_I = params[6];	//[1/min]
    double c_1 = params[7];
    double c_2 = params[8];

    double q = y_i[0];		//[m^3/s]
    double s_p = y_i[1];	//[m]
    double s_s = y_i[2];	//[m]

                            //Evaporation
    double C_p, C_s, C_T;
    if (e_pot > 0.0)
    {
        C_p = s_p / e_pot;
        C_s = s_s / e_pot;
        C_T = C_p + C_s;
    }
    else
    {
        C_p = 0.0;
        C_s = 0.0;
        C_T = 0.0;
    }

    double Corr_evap = (C_T > 1.0) ? 1.0 / C_T : 1.0;
    double e_p = Corr_evap * C_p * e_pot;
    double e_s = Corr_evap * C_s * e_pot;

    //Fluxes
    double q_pl = k_2 * s_p;
    double q_ps = k_I * s_p * pow(1.0 - s_s / h_b, gamma);
    double q_sl = k_3 * sq(s_s * L / A_h) * (h_r + h_b);

    //Discharge
    ans[0] = -q + (q_pl + q_sl) * c_2;
    for (unsigned short i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];

    //Hillslope
    ans[1] = forcing_values[0] * c_1 - q_pl - q_ps - e_p;
    ans[2] = q_ps - q_sl - e_s;
}


//Type 252
//TopLayer Hillslope
//Contains 3 layers on hillslope: ponded, top layer, soil
//Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
//The numbering is:	0   1   2     3    4   5   6   7
//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent
//The numbering is:        0      1        2     3   4     5        6   7  8 9  10
void model252(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;

    double lambda_1 = global_params[1];
    double k_3 = global_params[4];	//[1/min]
    double h_b = global_params[6];	//[m]
    double S_L = global_params[7];	//[m]
    double A = global_params[8];
    double B = global_params[9];
    double exponent = global_params[10];
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]
                                                                    //double e_pot = global_params[11];	//[m/min]
                                                                    //double e_pot = global_params[11] * (1e-3*60.0);	//[m/min]
                                                                    //double e_pot = 0.0;

    double L = params[1];	//[m]
    double A_h = params[2];	//[m^2]
                                //double h_r = params[3];	//[m]
    double invtau = params[3];	//[1/min]
    double k_2 = params[4];	//[1/min]
    double k_i = params[5];	//[1/min]
    double c_1 = params[6];
    double c_2 = params[7];

    double q = y_i[0];	//[m^3/s]
    double s_p = y_i[1];	//[m]
    double s_t = y_i[2];	//[m]
    double s_s = y_i[3];	//[m]

                            //Evaporation
    double e_p, e_t, e_s;
    double Corr = s_p + s_t / S_L + s_s / (h_b - S_L);
    if (e_pot > 0.0 && Corr > 1e-12)
    {
        e_p = s_p  * e_pot / Corr;
        e_t = s_t / S_L * e_pot / Corr;
        e_s = s_s / (h_b - S_L) * e_pot / Corr;
    }
    else
    {
        e_p = 0.0;
        e_t = 0.0;
        e_s = 0.0;
    }

    double pow_term = (1.0 - s_t / S_L > 0.0) ? pow(1.0 - s_t / S_L, exponent) : 0.0;
    double k_t = (A + B * pow_term) * k_2;

    //Fluxes
    double q_pl = k_2 * s_p;
    double q_pt = k_t * s_p;
    double q_ts = k_i * s_t;
    double q_sl = k_3 * s_s;

    //Discharge
    ans[0] = -q + (q_pl + q_sl) * c_2;
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];

    //Hillslope
    ans[1] = forcing_values[0] * c_1 - q_pl - q_pt - e_p;
    ans[2] = q_pt - q_ts - e_t;
    ans[3] = q_ts - q_sl - e_s;
}


//Type 253 , 255,256
//Contains 3 layers on hillslope: ponded, top layer, soil
//Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
//The numbering is:	0   1   2     3    4   5   6   7
//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent
//The numbering is:        0      1        2     3   4     5        6   7  8 9  10
void TopLayerHillslope_Reservoirs(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{	
	if(forcing_values[2] >0){
		ans[0] = forcing_values[2];
	}
    if(forcing_values[2] <=0){
		unsigned short i;
		for (i = 0; i<num_parents; i++)
			ans[0] += y_p[i * dim];
	}
    ans[1] = 0.0;
    ans[2] = 0.0;
    ans[3] = 0.0;
}
//Type 256
void model256_reservoirs(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{	
	if(forcing_values[2] >0){
		ans[0] = forcing_values[2];
        //Copy Discharge open loop
        double invtau = params[3];  //[1/min]
        double lambda_1 = global_params[1];
        double q_openloop = y_i[8];	
        //ans[5] = -q_openloop + ((q_pLink + q_sLink) * A_h / 60.0);
        unsigned short i;
        for (i =0; i<num_parents; i++)
            ans[8] += y_p[i*dim + 8];
        
        ans[8] = invtau * pow(q_openloop, lambda_1) * ans[8];
	}
    if(forcing_values[2] <=0){
		unsigned short i;
		for (i = 0; i<num_parents; i++)
			ans[0] += y_p[i * dim];
	}
    ans[1] = 0.0;
    ans[2] = 0.0;
    ans[3] = 0.0;
}

//type 249
void model249(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;

    double lambda_1 = global_params[1];
    double k_3 = global_params[4];	//[1/min]
    double h_b = global_params[6];	//[m]
    double S_L = global_params[7];	//[m]
    double A = global_params[8];
    double B = global_params[9];
    double exponent = global_params[10];
    double v_B = global_params[11];
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]

    double L = params[1];	//[m]
    double A_h = params[2];	//[m^2]
                                //double h_r = params[3];	//[m]
    double invtau = params[3];	//[1/min]
    double k_2 = params[4];	//[1/min]
    double k_i = params[5];	//[1/min]
    double c_1 = params[6];
    double c_2 = params[7];

    double q =   y_i[0];		//[m^3/s]
    double s_p = y_i[1];	//[m]
    double s_t = y_i[2];	//[m]
    double s_s = y_i[3];	//[m]
                            //double s_precip = y_i[4];	//[m]
                            //double V_r = y_i[5];	//[m^3]
    double q_b = (1.0e-7>y_i[4])? 1.0e-7: y_i[4];	//[m^3/s]
    double q_openloop = y_i[5];

                            //Evaporation
    double e_p, e_t, e_s;
    double Corr = s_p + s_t / S_L + s_s / (h_b - S_L);
    if (e_pot > 0.0 && Corr > 1e-12)
    {
        e_p = s_p * e_pot / Corr;
        e_t = s_t / S_L * e_pot / Corr;
        e_s = s_s / (h_b - S_L) * e_pot / Corr;
    }
    else
    {
        e_p = 0.0;
        e_t = 0.0;
        e_s = 0.0;
    }

    double pow_term = (1.0 - s_t / S_L > 0.0) ? pow(1.0 - s_t / S_L, exponent) : 0.0;
    double k_t = (A + B * pow_term) * k_2;

    //Fluxes
    double q_pl = k_2 * s_p;
    double q_pt = k_t * s_p;
    double q_ts = k_i * s_t;
    double q_sl = k_3 * s_s;	//[m/min]

                                //Discharge
    ans[0] = -q + (q_pl + q_sl) * c_2;
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];
    
    //ans[5] = ans[0];
    ans[5] = -q_openloop+ (q_pl + q_sl) * c_2;
    for (i = 0; i<num_parents; i++)
        ans[5] += y_p[i * dim+5];
    ans[5] = invtau * pow(q_openloop, lambda_1) * ans[5];

    //Hillslope
    ans[1] = forcing_values[0] * c_1 - q_pl - q_pt - e_p;
    ans[2] = q_pt - q_ts - e_t;
    ans[3] = q_ts - q_sl - e_s;
    ans[4] = q_sl * A_h - q_b*60.0;
    for (i = 0; i<num_parents; i++)
        ans[4] += y_p[i * dim + 4] ;
    //ans[6] += k_3*y_p[i].ve[3]*A_h;
    ans[4] = invtau * pow(q, lambda_1) * ans[4];
}
void model249_reservoirs(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;

    double lambda_1 = global_params[1];
    double k_3 = global_params[4];	//[1/min]
    double h_b = global_params[6];	//[m]
    double S_L = global_params[7];	//[m]
    double A = global_params[8];
    double B = global_params[9];
    double exponent = global_params[10];
    double v_B = global_params[11];
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]

    double L = params[1];	//[m]
    double A_h = params[2];	//[m^2]
                                //double h_r = params[3];	//[m]
    double invtau = params[3];	//[1/min]
    double k_2 = params[4];	//[1/min]
    double k_i = params[5];	//[1/min]
    double c_1 = params[6];
    double c_2 = params[7];

    //double q =   y_i[0];		//[m^3/s]
    double q =  y_i[0];
    double s_p = y_i[1];	//[m]
    double s_t = y_i[2];	//[m]
    double s_s = y_i[3];	//[m]
                            //double s_precip = y_i[4];	//[m]
                            //double V_r = y_i[5];	//[m^3]
    double q_b =  y_i[4];	//[m^3/s]
    //double q_openloop = (0.001>y_i[5])? 0.001 : y_i[5];
    //double q_openloop = y_i[5];
    double q_openloop = y_i[0];
    //double q_b = (1.0e-7>y_i[4])? 1.0e-7: y_i[4];   
                            //Evaporation
    double e_p, e_t, e_s;
    double Corr = s_p + s_t / S_L + s_s / (h_b - S_L);
    if (e_pot > 0.0 && Corr > 1e-12)
    {
        e_p = s_p * e_pot / Corr;
        e_t = s_t / S_L * e_pot / Corr;
        e_s = s_s / (h_b - S_L) * e_pot / Corr;
    }
    else
    {
        e_p = 0.0;
        e_t = 0.0;
        e_s = 0.0;
    }

    double pow_term = (1.0 - s_t / S_L > 0.0) ? pow(1.0 - s_t / S_L, exponent) : 0.0;
    double k_t = (A + B * pow_term) * k_2;

    //Fluxes
    double q_pl = k_2 * s_p;
    double q_pt = k_t * s_p;
    double q_ts = k_i * s_t;
    double q_sl = k_3 * s_s;	//[m/min]
    //Discharge data assim
    if(forcing_values[2] >0){
		ans[0] = forcing_values[2];
	}

    //Discharge open loop
    ans[5] = -q_openloop + (q_pl + q_sl) * c_2;
    for (i = 0; i<num_parents; i++)
        //ans[5] += y_p[i * dim+5];
        ans[5] += y_p[i * dim];
    ans[5] = invtau * pow(q_openloop, lambda_1) * ans[5];
    printf(" q_openloop: %f invtau: %f lambda1: %f\n", q_openloop, invtau, lambda_1);
    printf(" ans5: %f \n", ans[5]);
    // if(forcing_values[2] <=0){
    //     ans[0] =ans[5];
	// }
    //Hillslope
    ans[1] = forcing_values[0] * c_1 - q_pl - q_pt - e_p;
    ans[2] = q_pt - q_ts - e_t;
    ans[3] = q_ts - q_sl - e_s;
    ans[4] = q_sl * A_h - q_b*60.0;
    for (i = 0; i<num_parents; i++)
        ans[4] += y_p[i * dim + 4];
    ans[4] = invtau * pow(q_b, lambda_1) * ans[4];
}
//Type 251
// 4 states.
//soil layer thick is variable in time using a monthly forcing.
void model251(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;

    double lambda_1 = global_params[1];
    double k_3 = global_params[4];	//[1/min]
    double h_b = global_params[6];	//[m]
    //double S_L = global_params[7];	//[m]
    double S_L = forcing_values[4]; //[m]
    double A = global_params[8];
    double B = global_params[9];
    double exponent = global_params[10];
    double v_B = global_params[11];
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]
    double snowmelt = forcing_values[3] * (0.001)/60; //mm/hour to m/min
    double L = params[1];	//[m]
    double A_h = params[2];	//[m^2]
                                //double h_r = params[3];	//[m]
    double invtau = params[3];	//[1/min]
    double k_2 = params[4];	//[1/min]
    double k_i = params[5];	//[1/min]
    double c_1 = params[6];
    double c_2 = params[7];
	double rainfall = forcing_values[0] * c_1; //mm/hour to m/min
	rainfall += snowmelt;

    double q =   y_i[0];		//[m^3/s]
    double s_p = y_i[1];	//[m]
    double s_t = y_i[2];	//[m]
    double s_s = y_i[3];	//[m]
                            //double s_precip = y_i[4];	//[m]
                            //double V_r = y_i[5];	//[m^3]
    

                            //Evaporation
    double e_p, e_t, e_s;
    double Corr = s_p + s_t / S_L + s_s / (h_b - S_L);
    if (e_pot > 0.0 && Corr > 1e-12)
    {
        e_p = s_p * e_pot / Corr;
        e_t = s_t / S_L * e_pot / Corr;
        e_s = s_s / (h_b - S_L) * e_pot / Corr;
    }
    else
    {
        e_p = 0.0;
        e_t = 0.0;
        e_s = 0.0;
    }

    double pow_term = (1.0 - s_t / S_L > 0.0) ? pow(1.0 - s_t / S_L, exponent) : 0.0;
    double k_t = (A + B * pow_term) * k_2;

    //Fluxes
    double q_pl = k_2 * s_p;
    double q_pt = k_t * s_p;
    double q_ts = k_i * s_t;
    double q_sl = k_3 * s_s;	//[m/min]

                                //Discharge
    ans[0] = -q + (q_pl + q_sl) * c_2;
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];

    //Hillslope
    ans[1] = rainfall - q_pl - q_pt - e_p;
    ans[2] = q_pt - q_ts - e_t;
    ans[3] = q_ts - q_sl - e_s;

}
// void model249_reservoirs(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
// {
//     ans[1]=0;
//     ans[2]=0;
//     ans[3]=0;
//     double invtau = params[3];	//[1/min]
//     double lambda_1 = global_params[1];
//     double q_openloop = y_i[5];
//     //printf("time: %f\n", t);
//     //printf(" forc: %f %f\n", forcing_values[0], forcing_values[1]);
//     //printf(" flux: %f %f\n", ans[1], ans[2]);
//     double q_b = y_i[4];	//[m^3/s]

//     unsigned short i;
//     //Discharge data assim
//     if(forcing_values[2] >0){
// 		ans[0] = forcing_values[2];
// 	}

//     //Discharge open loop
//     ans[5] = -q_openloop ;
//     for (i = 0; i<num_parents; i++){
//         ans[5] += y_p[i * dim+5];
//     printf(" upstream state5: %f \n", y_p[i * dim+5]);
//     }
        
//     ans[5] = invtau * pow(q_openloop, lambda_1) * ans[5];
//     //printf(" q_openloop: %f invtau: %f lambda1: %f\n", q_openloop, invtau, lambda_1);
//     //printf(" ans5: %f \n", ans[5]);

//     // if(forcing_values[2] <=0){
//     //     ans[0] =ans[5];
//     //     //ans[0]=1;
// 	// }

//     ans[4] = - q_b;
//     for (i = 0; i<num_parents; i++)
//         ans[4] += y_p[i * dim + 4];
//     ans[4] = invtau * pow(q_b, lambda_1) * ans[4];

// }
//Type 253
// 4 states.
void model253(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;

    double lambda_1 = global_params[1];
    double k_3 = global_params[4];	//[1/min]
    double h_b = global_params[6];	//[m]
    double S_L = global_params[7];	//[m]
    double A = global_params[8];
    double B = global_params[9];
    double exponent = global_params[10];
    double v_B = global_params[11];
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]
    double snowmelt = forcing_values[3] * (0.001)/60; //mm/hour to m/min
    double L = params[1];	//[m]
    double A_h = params[2];	//[m^2]
                                //double h_r = params[3];	//[m]
    double invtau = params[3];	//[1/min]
    double k_2 = params[4];	//[1/min]
    double k_i = params[5];	//[1/min]
    double c_1 = params[6];
    double c_2 = params[7];
	double rainfall = forcing_values[0] * c_1; //mm/hour to m/min
	rainfall += snowmelt;

    double q =   y_i[0];		//[m^3/s]
    double s_p = y_i[1];	//[m]
    double s_t = y_i[2];	//[m]
    double s_s = y_i[3];	//[m]
                            //double s_precip = y_i[4];	//[m]
                            //double V_r = y_i[5];	//[m^3]
    

                            //Evaporation
    double e_p, e_t, e_s;
    double Corr = s_p + s_t / S_L + s_s / (h_b - S_L);
    if (e_pot > 0.0 && Corr > 1e-12)
    {
        e_p = s_p * e_pot / Corr;
        e_t = s_t / S_L * e_pot / Corr;
        e_s = s_s / (h_b - S_L) * e_pot / Corr;
    }
    else
    {
        e_p = 0.0;
        e_t = 0.0;
        e_s = 0.0;
    }

    double pow_term = (1.0 - s_t / S_L > 0.0) ? pow(1.0 - s_t / S_L, exponent) : 0.0;
    double k_t = (A + B * pow_term) * k_2;

    //Fluxes
    double q_pl = k_2 * s_p;
    double q_pt = k_t * s_p;
    double q_ts = k_i * s_t;
    double q_sl = k_3 * s_s;	//[m/min]

                                //Discharge
    ans[0] = -q + (q_pl + q_sl) * c_2;
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];

    //Hillslope
    ans[1] = rainfall - q_pl - q_pt - e_p;
    ans[2] = q_pt - q_ts - e_t;
    ans[3] = q_ts - q_sl - e_s;

}

//Type 254
//Contains 3 layers on hillslope: ponded, top layer, soil. Also has 3 extra states: total precip, total runoff, base flow
//Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
//The numbering is:	0   1   2     3    4   5   6   7
//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent,v_B
//The numbering is:        0      1        2     3   4     5        6   7  8 9  10       11
void model254(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;

    double lambda_1 = global_params[1];
    double k_3 = global_params[4];	//[1/min]
    double h_b = global_params[6];	//[m]
    double S_L = global_params[7];	//[m]
    double A = global_params[8];
    double B = global_params[9];
    double exponent = global_params[10];
    double v_B = global_params[11];
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]

    double L = params[1];	//[m]
    double A_h = params[2];	//[m^2]
                                //double h_r = params[3];	//[m]
    double invtau = params[3];	//[1/min]
    double k_2 = params[4];	//[1/min]
    double k_i = params[5];	//[1/min]
    double c_1 = params[6];
    double c_2 = params[7];

    double q =   y_i[0];		//[m^3/s]
    double s_p = y_i[1];	//[m]
    double s_t = y_i[2];	//[m]
    double s_s = y_i[3];	//[m]
                            //double s_precip = y_i[4];	//[m]
                            //double V_r = y_i[5];	//[m^3]
    double q_b = max(0.001,y_i[6]);	//[m^3/s]

                            //Evaporation
    double e_p, e_t, e_s;
    double Corr = s_p + s_t / S_L + s_s / (h_b - S_L);
    if (e_pot > 0.0 && Corr > 1e-12)
    {
        e_p = s_p * e_pot / Corr;
        e_t = s_t / S_L * e_pot / Corr;
        e_s = s_s / (h_b - S_L) * e_pot / Corr;
    }
    else
    {
        e_p = 0.0;
        e_t = 0.0;
        e_s = 0.0;
    }

    double pow_term = (1.0 - s_t / S_L > 0.0) ? pow(1.0 - s_t / S_L, exponent) : 0.0;
    double k_t = (A + B * pow_term) * k_2;

    //Fluxes
    double q_pl = k_2 * s_p;
    double q_pt = k_t * s_p;
    double q_ts = k_i * s_t;
    double q_sl = k_3 * s_s;	//[m/min]

                                //Discharge
    ans[0] = -q + (q_pl + q_sl) * c_2;
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];

    //Hillslope
    ans[1] = forcing_values[0] * c_1 - q_pl - q_pt - e_p;
    ans[2] = q_pt - q_ts - e_t;
    ans[3] = q_ts - q_sl - e_s;

    //Additional states

    ans[4] = forcing_values[0] * c_1;
    ans[5] = q_pl;


    ans[6] = q_sl * A_h - q_b*60.0;
    for (i = 0; i<num_parents; i++)
        ans[6] += y_p[i * dim + 6] * 60.0;
    //ans[6] += k_3*y_p[i].ve[3]*A_h;
    ans[6] *= v_B / L;
}


//Type 255, before called TopLayerHillslope_variable
//Contains 2 layers in the channel: discharge, storage. Contains 3 layers on hillslope: ponded, top layer, soil.
//Order of the states is:              0          1                                        2        3       4
//Order of parameters: A_i,L_i,A_h,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent | invtau,k_2,k_i,c_1,c_2
//The numbering is:	0   1   2   3   4      5       6   7  8 9   10        11    12  13  14  15
//Order of global_params: v_0,lambda_1,lambda_2
//The numbering is:        0      1        2
void model255(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;

    double lambda_1 = global_params[1];
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]

    double L = params[1];	//[m]
    double A_h = params[2];	//[m^2]
    double k_3 = params[4];	//[1/min]
    double h_b = params[6];	//[m]
    double S_L = params[7];	//[m]
    double A = params[8];
    double B = params[9];
    double exponent = params[10];
    double k_2 = params[12];	//[1/min]
    double k_i = params[13];	//[1/min]
    double c_1 = params[14];
    double c_2 = params[15];
    double invtau = params[11];	//[1/min]

    double rainfall = forcing_values[0] * c_1;
    double snowmelt = forcing_values[3] * (0.001)/60; //mm/hour to m/min
    rainfall+= snowmelt;
    
    double q = y_i[0];		//[m^3/s]
    double S = y_i[1];		//[m^3]
    double s_p = y_i[2];	//[m]
    double s_t = y_i[3];	//[m]
    double s_s = y_i[4];	//[m]

                            //Evaporation
    double e_p, e_t, e_s;
    double Corr = s_p + s_t / S_L + s_s / (h_b - S_L);
    if (e_pot > 0.0 && Corr > 1e-12)
    {
        e_p = s_p * e_pot / Corr;
        e_t = s_t / S_L * e_pot / Corr;
        e_s = s_s / (h_b - S_L) * e_pot / Corr;
    }
    else
    {
        e_p = 0.0;
        e_t = 0.0;
        e_s = 0.0;
    }

    double pow_term = (1.0 - s_t / S_L > 0.0) ? pow(1.0 - s_t / S_L, exponent) : 0.0;
    double k_t = (A + B * pow_term) * k_2;

    //Fluxes
    double q_pl = k_2 * s_p;
    double q_pt = k_t * s_p;
    double q_ts = k_i * s_t;
    double q_sl = k_3 * s_s;

    //Discharge
    dam_model255(y_i, dim, global_params, params, qvs, state, user, ans);	//ans is used for convenience !!!! Is q available in y_i? !!!!
   
    double qm = ans[0] * 60.0;

    //Storage
    ans[1] = (q_pl + q_sl) * A_h - qm;
    for (i = 0; i<num_parents; i++)
        ans[1] += y_p[i * dim] * 60.0;

    //Hillslope
    ans[2] = rainfall - q_pl - q_pt - e_p;
    ans[3] = q_pt - q_ts - e_t;
    ans[4] = q_ts - q_sl - e_s;
}

//Type 255 before called dam_TopLayerHillslope_variable
//Contains 2 layers in the channel: discharge, storage. Contains 3 layers on hillslope: ponded, top layer, soil.
//Order of the states is:              0          1                                        2        3       4
//Order of parameters: A_i,L_i,A_h,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent | invtau,k_2,k_i,c_1,c_2
//The numbering is:	0   1   2   3   4      5       6   7  8 9   10        11    12  13  14  15
//Order of global_params: v_0,lambda_1,lambda_2
//The numbering is:        0      1        2
void dam_model255(const double * const y_i, unsigned int num_dof, const double * const global_params, const double * const params, const QVSData * const qvs, int state, void* user, double *ans)
{
    double q1, q2, S1, S2, S_max, q_max, S;

    //Parameters
    double lambda_1 = global_params[1];
    double invtau = params[11];	//[1/min]

                                    //Find the discharge in [m^3/s]
    if (state == -1)
    {
        S = (y_i[1] < 0.0) ? 0.0 : y_i[1];
        ans[0] = invtau/60.0*pow(S,1.0/(1.0-lambda_1));
        //ans[0] = pow((1.0 - lambda_1)*invtau / 60.0 * S, 1.0 / (1.0 - lambda_1));
        
        
    }
    else if (state == (int)qvs->n_values - 1)
    {
        S_max = qvs->points[qvs->n_values - 1][0];
        q_max = qvs->points[qvs->n_values - 1][1];
        ans[0] = q_max;
    }
    else
    {
        S = (y_i[1] < 0.0) ? 0.0 : y_i[1];
        q2 = qvs->points[state + 1][1];
        q1 = qvs->points[state][1];
        S2 = qvs->points[state + 1][0];
        S1 = qvs->points[state][0];
        ans[0] = (q2 - q1) / (S2 - S1) * (S - S1) + q1;
    }
}

//Type 255 before called dam_TopLayerHillslope_variable
//Contains 2 layers in the channel: discharge, storage. Contains 3 layers on hillslope: ponded, top layer, soil.
//Order of the states is:              0          1                                        2        3       4
//Order of parameters: A_i,L_i,A_h,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent | invtau,k_2,k_i,c_1,c_2
//The numbering is:	0   1   2   3   4      5       6   7  8 9   10        11    12  13  14  15
//Order of global_params: v_0,lambda_1,lambda_2
//The numbering is:        0      1        2
void dam_TopLayerHillslope_variable(const double * const y_i, unsigned int num_dof, const double * const global_params, const double * const params, const QVSData * const qvs, int state, void* user, double *ans)
{
    double q1, q2, S1, S2, S_max, q_max, S;

    //Parameters
    double lambda_1 = global_params[1];
    double invtau = params[11];	//[1/min]

                                    //Find the discharge in [m^3/s]
    if (state == -1)
    {
        S = (y_i[1] < 0.0) ? 0.0 : y_i[1];
        //ans[0] = invtau/60.0*pow(S,1.0/(1.0-lambda_1));
        ans[0] = pow((1.0 - lambda_1)*invtau / 60.0 * S, 1.0 / (1.0 - lambda_1));
    }
    else if (state == (int)qvs->n_values - 1)
    {
        S_max = qvs->points[qvs->n_values - 1][0];
        q_max = qvs->points[qvs->n_values - 1][1];
        ans[0] = q_max;
    }
    else
    {
        S = (y_i[1] < 0.0) ? 0.0 : y_i[1];
        q2 = qvs->points[state + 1][1];
        q1 = qvs->points[state][1];
        S2 = qvs->points[state + 1][0];
        S1 = qvs->points[state][0];
        ans[0] = (q2 - q1) / (S2 - S1) * (S - S1) + q1;
    }
}
//Type 256
//TopLayerHillslope_even_more_extras.
//Contains 3 layers on hillslope: ponded, top layer, soil.
// Also has 3 extra states: total precip, total ET, total runoff, base flow
//Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
//The numbering is:	0   1   2     3    4   5   6   7
//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent,v_B,k_tl
//The numbering is:        0      1        2     3   4     5        6   7  8 9  10       11   12
void model256(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;

    double lambda_1 = global_params[1];
    double k_3 = global_params[4];  //[1/min]
    double h_b = global_params[6];  //[m]
    double S_L = global_params[7];  //[m]
    double A = global_params[8];
    double B = global_params[9];
    double exponent = global_params[10];
    double v_B = global_params[11];
    double k_tl = global_params[12];
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));   //[mm/month] -> [m/min]

    double L = params[1];   //[m]
    double A_h = params[2]; //[m^2]
                            //double h_r = params[3];	//[m]
    double invtau = params[3];  //[1/min]
    double k_2 = params[4];     //[1/min]
    double k_i = params[5];     //[1/min]
    double c_1 = params[6];
    double c_2 = params[7];

    double q = y_i[0];      //[m^3/s]
    double s_p = y_i[1];    //[m]
    double s_t = y_i[2];    //[m]
    double s_s = y_i[3];    //[m]
	double q_pl0=y_i[4];	//m3s-1 
	double q_tl0=y_i[5];	//m3s-1
	double q_sl0=y_i[6];	//m3s-1
                            //double s_precip = y_i[4];	//[m]
                            //double V_r = y_i[5];	//[m^3]
    double q_b = y_i[7];    //[m^3/s]
    double q_openloop = y_i[8];    //[m^3/s]

    //Evaporation
    double e_p, e_t, e_s;
    double Corr = s_p + s_t / S_L + s_s / (h_b - S_L);
    if (e_pot > 0.0 && Corr > 1e-12)
    {
        e_p = s_p * e_pot / Corr;
        e_t = s_t / S_L * e_pot / Corr;
        e_s = s_s / (h_b - S_L) * e_pot / Corr;
    }
    else
    {
        e_p = 0.0;
        e_t = 0.0;
        e_s = 0.0;
    }

    double pow_term = (1.0 - s_t / S_L > 0.0) ? pow(1.0 - s_t / S_L, exponent) : 0.0;
    double k_t = (A + B * pow_term) * k_2;

    //Fluxes
    double q_pl = k_2 * s_p;
    double q_pt = k_t * s_p;
    double q_ts = k_i * s_t;
    double q_tl = k_tl * s_t;
    double q_sl = k_3 * s_s;    //[m/min]

    //Discharge
    ans[0] = -q + (q_pl + q_tl + q_sl) * c_2;
    for (i = 0; i < num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];    // discharge[0]

    //Hillslope
    ans[1] = forcing_values[0] * c_1 - q_pl - q_pt - e_p;   // hillslope pond[1] depth m or m/min
    ans[2] = q_pt - q_ts - q_tl - e_t;                      // toplayer[2] depth m or m/min
    ans[3] = q_ts - q_sl - e_s;                             // subsurface[3] depth m or m/min

    //Additional states
    //ans[4] = forcing_values[0] * c_1;   // precip[4]
    //ans[5] = forcing_values[1] * c_1;   // et[5]
	//ans[6] = q_pl; //runoff
	
	ans[4] = q_pl *c_2 - q_pl0; //
	ans[5] = q_tl*c_2 - q_tl0; //
	ans[6] = q_sl*c_2 - q_sl0; //
	//ans[4] = q_pl *c_2 ; //hillslope surface flow m3s-1
	//ans[5] = q_tl*c_2 ; //hillslope interflow m3s-1
	//ans[6] = q_sl*c_2 ; //hillslope subsurface flow m3s-1
    ans[7] = q_sl * A_h - q_b*60.0;     // baseflow[7]
    for (i = 0; i < num_parents; i++)
        ans[7] += y_p[i * dim + 7] * 60.0;
    ans[7] *= v_B / L;
    //open loop
    ans[8] = -q_openloop + (q_pl + q_tl + q_sl) * c_2;
    for (i = 0; i < num_parents; i++)
        ans[8] += y_p[i * dim + 8];
    ans[8] = invtau * pow(q, lambda_1) * ans[8];
}



//Type 263: similar to model 256, with distributed parameters
//Contains 3 layers on hillslope: ponded, top layer, soil. Also has 3 extra states: total precip, total ET, total runoff, base flow
//Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
//The numbering is:	0   1   2     3    4   5   6   7
//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent,v_B,k_tl
//The numbering is:        0      1        2     3   4     5        6   7  8 9  10       11   12
void model263(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;

    //double lambda_1 = global_params[1];
	double A_i = params[0];
	double L = params[1];
	double A_h = params[2];
    double v_0=params[3];
    double lambda_1=params[4];
    double lambda_2=params[5];
    double invtau = 60.0 * v_0 * pow(A_i, lambda_2) / ((1.0 - lambda_1) * L);//[1/min]  invtau
    double v_h=params[6];
    double k_2=v_h * L / A_h * 60.0;	//[1/min] k_2
    //double k_i_factor1 = params[7];
    double k_i  = k_2*params[7];

    //double k_3 = global_params[4];  //[1/min]
    double k_3 = params[8];

    //double h_b = global_params[6];  //[m]
    double h_b = params[9];

    //double S_L = global_params[7];  //[m]
    double S_L = params[10];

    double A = params[11];//global_params[8];
    double B = params[12];//global_params[9];
    double exponent = params[13];//global_params[10];
    double v_B = params[14];//global_params[11];
    double k_tl = params[15];//]global_params[12];
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));   //[mm/month] -> [m/min]

    //double L = params[1];   //[m]
    //double A_h = params[2]; //[m^2]
                            //double h_r = params[3];	//[m]
    //double invtau = params[3];  //[1/min]
    //double k_2 = params[4];     //[1/min]
    //double k_i = params[5];     //[1/min]
    //double c_1 = params[6];
    //double c_2 = params[7];
	double c_1 = (0.001 / 60.0);		//(mm/hr->m/min)  c_1
	double c_2 = A_h / 60.0;	//  c_2

    double q = y_i[0];      //[m^3/s]
    double s_p = y_i[1];    //[m]
    double s_t = y_i[2];    //[m]
    double s_s = y_i[3];    //[m]
                            //double s_precip = y_i[4];	//[m]
                            //double V_r = y_i[5];	//[m^3]
    double q_b = y_i[4];    //[m^3/s]

    //Evaporation
    double e_p, e_t, e_s;
    double Corr = s_p + s_t / S_L + s_s / (h_b - S_L);
    if (e_pot > 0.0 && Corr > 1e-12)
    {
        e_p = s_p * e_pot / Corr;
        e_t = s_t / S_L * e_pot / Corr;
        e_s = s_s / (h_b - S_L) * e_pot / Corr;
    }
    else
    {
        e_p = 0.0;
        e_t = 0.0;
        e_s = 0.0;
    }

    double pow_term = (1.0 - s_t / S_L > 0.0) ? pow(1.0 - s_t / S_L, exponent) : 0.0;
    double k_t = (A + B * pow_term) * k_2;

    //Fluxes
    double q_pl = k_2 * s_p;
    double q_pt = k_t * s_p;
    double q_ts = k_i * s_t;
    double q_tl = k_tl * s_t;
    double q_sl = k_3 * s_s;    //[m/min]

    //Discharge
    ans[0] = -q + (q_pl + q_tl + q_sl) * c_2;
    for (i = 0; i < num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];    // discharge[0]

    //Hillslope
    ans[1] = forcing_values[0] * c_1 - q_pl - q_pt - e_p;   // pond[1]
    ans[2] = q_pt - q_ts - q_tl - e_t;                      // toplayer[2]
    ans[3] = q_ts - q_sl - e_s;                             // subsurface[3]

    //Additional states
    //ans[4] = forcing_values[0] * c_1;   // precip[4]
    //ans[5] = forcing_values[1] * c_1;   // et[5]
    //ans[6] = q_pl;                      // runoff[]6
    //ans[7] = q_sl * A_h - q_b*60.0;     // baseflow[7]
    for (i = 0; i < num_parents; i++)
        //ans[4] += y_p[i * dim + 4] * 60.0;
        ans[4] += y_p[i * dim + 4] ;
    //ans[4] *= v_B / L;
    //pow using q not q_b to move at same flow velocity that state0
    ans[4] = invtau * pow(q, lambda_1) * ans[4];    // baseflow[0]
}

//Type 264: similar to model 256, with a forcing for snowmelt
//Contains 3 layers on hillslope: ponded, top layer, soil. Also has 3 extra states: total precip, total ET, total runoff, base flow
//Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
//The numbering is:	0   1   2     3    4   5   6   7
//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent,v_B,k_tl
//The numbering is:        0      1        2     3   4     5        6   7  8 9  10       11   12
void model264(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
     unsigned short i;

    double lambda_1 = global_params[1];
    double k_3 = global_params[4];  //[1/min]
    double h_b = global_params[6];  //[m]
    double S_L = global_params[7];  //[m]
    double A = global_params[8];
    double B = global_params[9];
    double exponent = global_params[10];
    double v_B = global_params[11];
    double k_tl = global_params[12];
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));   //[mm/month] -> [m/min]
	double snowmelt = forcing_values[3] * (0.001)/60; //mm/hour to m/min
    double L = params[1];   //[m]
    double A_h = params[2]; //[m^2]
                            //double h_r = params[3];	//[m]
    double invtau = params[3];  //[1/min]
    double k_2 = params[4];     //[1/min]
    double k_i = params[5];     //[1/min]
    double c_1 = params[6];//converts mm/hour to m/min
    double c_2 = params[7];//converts m/min to m3/s
	double rainfall = forcing_values[0] * c_1; //mm/hour to m/min
	rainfall += snowmelt;

    double q = y_i[0];      //[m^3/s]
    double s_p = y_i[1];    //[m]
    double s_t = y_i[2];    //[m]
    double s_s = y_i[3];    //[m]
	double q_pl0=y_i[4];	//m3s-1 bad idea. makes everything 0
	double q_tl0=y_i[5];	//m3s-1
	double q_sl0=y_i[6];	//m3s-1
                            //double s_precip = y_i[4];	//[m]
                            //double V_r = y_i[5];	//[m^3]
    double q_b = y_i[7];    //[m^3/s]


    //Evaporation
    double e_p, e_t, e_s;
    double Corr = s_p + s_t / S_L + s_s / (h_b - S_L);
    if (e_pot > 0.0 && Corr > 1e-12)
    {
        e_p = s_p * e_pot / Corr;
        e_t = s_t / S_L * e_pot / Corr;
        e_s = s_s / (h_b - S_L) * e_pot / Corr;
    }
    else
    {
        e_p = 0.0;
        e_t = 0.0;
        e_s = 0.0;
    }

    double pow_term = (1.0 - s_t / S_L > 0.0) ? pow(1.0 - s_t / S_L, exponent) : 0.0;
    double k_t = (A + B * pow_term) * k_2;

    //Fluxes
    double q_pl = k_2 * s_p;
    double q_pt = k_t * s_p;
    double q_ts = k_i * s_t;
    double q_tl = k_tl * s_t;
    double q_sl = k_3 * s_s;    //[m/min]

    //Discharge
    ans[0] = -q + (q_pl + q_tl + q_sl) * c_2;
    for (i = 0; i < num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];    // discharge[0]

    //Hillslope
    ans[1] = rainfall - q_pl - q_pt - e_p;   // hillslope pond[1] depth m or m/min
    ans[2] = q_pt - q_ts - q_tl - e_t;                      // toplayer[2] depth m or m/min
    ans[3] = q_ts - q_sl - e_s;                             // subsurface[3] depth m or m/min

    //Additional states
    //ans[4] = forcing_values[0] * c_1;   // precip[4]
    //ans[5] = forcing_values[1] * c_1;   // et[5]
	//ans[6] = q_pl; //runoff
	
	ans[4] = q_pl *c_2 - q_pl0; //
	ans[5] = q_tl*c_2 - q_tl0; //
	ans[6] = q_sl*c_2 - q_sl0; //
	//ans[4] = q_pl *c_2 ; //hillslope surface flow m3s-1
	//ans[5] = q_tl*c_2 ; //hillslope interflow m3s-1
	//ans[6] = q_sl*c_2 ; //hillslope subsurface flow m3s-1
    ans[7] = q_sl * A_h - q_b*60.0;     // baseflow[7]
    for (i = 0; i < num_parents; i++)
        ans[7] += y_p[i * dim + 7] * 60.0;
    ans[7] *= v_B / L;
}
//Type 257
//Contains 3 layers on hillslope: ponded, top layer, soil. Also has 4 extra states: total precip, total ET, total runoff, base flow
//Order of parameters: A_i,L_i,A_h,horder,invtau,k_2,k_i,c_1,c_2
//The numbering is:    0   1   2   3      4      5   6   7...8
//Order of global_params: v_0_0,...,v_0_9,lambda_1_0,...,lambda_1_9,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A ,B, exponent,v_B,k_tl
//The numbering is:       0         9     10,            19         20       21  22  23         24  25  26 27 28       29  30
void TopLayerHillslope_spatial_velocity(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;

    double k_3 = global_params[22];  //[1/min]
    double h_b = global_params[24];  //[m]
    double S_L = global_params[25];  //[m]
    double A = global_params[26];
    double B = global_params[27];
    double exponent = global_params[28];
    double v_B = global_params[29];
    double k_tl = global_params[30];
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));   //[mm/month] -> [m/min]

    double L = params[1];   //[m]
    double A_h = params[2]; //[m^2]
                            //double h_r = params[3];	//[m]
    int h_order = (int)params[3];
    assert(0 < h_order && h_order <= 10);
    double invtau = params[4];  //[1/min]
    double k_2 = params[5];     //[1/min]
    double k_i = params[6];     //[1/min]
    double c_1 = params[7];
    double c_2 = params[8];

    double lambda_1 = global_params[10 + h_order - 1];

    double q = y_i[0];      //[m^3/s]
    double s_p = y_i[1];    //[m]
    double s_t = y_i[2];    //[m]
    double s_s = y_i[3];    //[m]
                            //double s_precip = y_i[4];	//[m]
                            //double V_r = y_i[5];	//[m^3]
    double q_b = y_i[7];    //[m^3/s]

                            //Evaporation
    double e_p, e_t, e_s;
    double Corr = s_p + s_t / S_L + s_s / (h_b - S_L);
    if (e_pot > 0.0 && Corr > 1e-12)
    {
        e_p = s_p * 1e3 * e_pot / Corr;
        e_t = s_t / S_L * e_pot / Corr;
        e_s = s_s / (h_b - S_L) * e_pot / Corr;
    }
    else
    {
        e_p = 0.0;
        e_t = 0.0;
        e_s = 0.0;
    }

    double pow_term = (1.0 - s_t / S_L > 0.0) ? pow(1.0 - s_t / S_L, exponent) : 0.0;
    double k_t = (A + B * pow_term) * k_2;

    //Fluxes
    double q_pl = k_2 * s_p;
    double q_pt = k_t * s_p;
    double q_ts = k_i * s_t;
    double q_tl = k_tl * s_t;
    double q_sl = k_3 * s_s;    //[m/min]

                                //Discharge
    ans[0] = -q + (q_pl + q_tl + q_sl) * c_2;
    for (i = 0; i < num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];    // discharge[0]

                                                    //Hillslope
    ans[1] = forcing_values[0] * c_1 - q_pl - q_pt - e_p;   // pond[1]
    ans[2] = q_pt - q_ts - q_tl - e_t;                      // toplayer[2]
    ans[3] = q_ts - q_sl - e_s;                             // subsurface[3]

                                                            //Additional states
    ans[4] = forcing_values[0] * c_1;   // precip[4]
    ans[5] = forcing_values[1] * c_1;   // et[5]
    ans[6] = q_pl;                      // runoff[]6
    ans[7] = q_sl * A_h - q_b*60.0;     // baseflow[7]
    for (i = 0; i < num_parents; i++)
        ans[7] += y_p[i * dim + 7] * 60.0;
    ans[7] *= v_B / L;
}

//Type 258
//Contains 3 layers on hillslope: ponded, top layer, soil. Also has 3 extra states: total precip, total runoff, base flow
//Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
//The numbering is:	    0   1   2     3    4   5   6  7
//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,v_B
//The numbering is:        0      1        2     3   4      5       6   7   8
void OfflineTopLayerHillslope(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;

    double lambda_1 = global_params[1];
    double k_3 = global_params[4];	                               // [1/min]
    double h_b = global_params[6];	                               // [m]
    double S_L = global_params[7];	                               // [m]
    double v_B = global_params[8];
    double e_pot = forcing_values[2] * (1e-3 / (30.0*24.0*60.0));  // [mm/month] -> [m/min]

    double L = params[1];	                                       // [m]
    double A_h = params[2];	                                       // [m^2]
    double invtau = params[3];	                                   // [1/min]
    double k_2 = params[4];	                                       // [1/min]
    double k_i = params[5];	                                       // [1/min]
    double c_1 = params[6];
    double c_2 = params[7];

    double q = y_i[0];		//[m^3/s]
    double s_p = y_i[1];	//[m]
    double s_t = y_i[2];	//[m]
    double s_s = y_i[3];	//[m]
                            //double s_precip = y_i[4];	//[m]
                            //double V_r = y_i[5];	//[m^3]
    double q_b = y_i[6];	//[m^3/s]

                            //Evaporation
    double e_p, e_t, e_s;
    double Corr = s_p + s_t / S_L + s_s / (h_b - S_L);
    if (e_pot > 0.0 && Corr > 1e-12)
    {
        e_p = s_p * 1e3 * e_pot / Corr;
        e_t = s_t / S_L * e_pot / Corr;
        e_s = s_s / (h_b - S_L) * e_pot / Corr;
    }
    else
    {
        e_p = 0.0;
        e_t = 0.0;
        e_s = 0.0;
    }

    // fluxes precip->ponds and ponds-> links
    double q_ra_po = forcing_values[0] * c_1;             // (mm/hr -> m/min)
    double q_po_li = k_2 * s_p;                           // (1/min * m)

    // fluxes precip->top and top->subsurface
    double q_ra_to = forcing_values[1] * c_1;             // (mm/hr -> m/min)
    double q_to_ss = k_i * s_t;
    
    // flux subsurface -> link
    double q_ss_li = k_3 * s_s;                           // (1/min * m)

    //Discharge
    ans[0] = -q + (q_po_li + q_ss_li) * c_2;              // main streamflow
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];

    //Hillslope
    ans[1] = q_ra_po - q_po_li - e_p;                     // ponds
    ans[2] = q_ra_to - q_to_ss - e_t;                     // top layer
    ans[3] = q_to_ss - q_ss_li - e_s;                     // sub surface

    //Additional states
    ans[4] = q_ra_po + q_ra_to;                           // acc precip
    ans[5] = q_ra_po;                                     // acc runoff
    ans[6] = e_pot;                                       // acc evaporation
    ans[7] = q_ss_li * A_h - q_b*60.0;                    // base flow
    for (i = 0; i<num_parents; i++)
        ans[7] += y_p[i * dim + 6] * 60.0;
    ans[7] *= v_B / L;
}

//Type 258
//Contains 3 layers on hillslope: ponded, top layer, soil
//Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
//The numbering is:	0   1   2     3    4   5   6   7
//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent
//The numbering is:        0      1        2     3   4     5        6   7  8 9  10
void OfflineTopLayerHillslope_Reservoirs(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    ans[0] = forcing_values[3];
    ans[1] = 0.0;
    ans[2] = 0.0;
    ans[3] = 0.0;
}

//Type 259
//Contains 3 layers on hillslope: ponded, top layer, soil. Also has 3 extra states: total precip, total runoff, base flow
//Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
//The numbering is:	    0   1   2     3    4   5   6  7
//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,v_B,k_tl
//The numbering is:        0      1        2     3   4      5       6   7   8   9
void OfflineTopLayerInterflowHillslope(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;

    double lambda_1 = global_params[1];
    double k_3 = global_params[4];	                               // [1/min]
    double h_b = global_params[6];	                               // [m]
    double S_L = global_params[7];	                               // [m]
    double v_B = global_params[8];
    double k_tl = global_params[9];

    double e_pot = forcing_values[2] * (1e-3 / (30.0*24.0*60.0));  // [mm/month] -> [m/min]

    double L = params[1];	                                       // [m]
    double A_h = params[2];	                                       // [m^2]
    double invtau = params[3];	                                   // [1/min]
    double k_2 = params[4];	                                       // [1/min]
    double k_i = params[5];	                                       // [1/min]
    double c_1 = params[6];
    double c_2 = params[7];

    double q = y_i[0];		//[m^3/s]
    double s_p = y_i[1];	//[m]
    double s_t = y_i[2];	//[m]
    double s_s = y_i[3];	//[m]
                            //double s_precip = y_i[4];	//[m]
                            //double V_r = y_i[5];	//[m^3]
    double q_b = y_i[6];	//[m^3/s]

                            //Evaporation
    double e_p, e_t, e_s;
    double Corr = s_p + s_t / S_L + s_s / (h_b - S_L);
    if (e_pot > 0.0 && Corr > 1e-12)
    {
        e_p = s_p * 1e3 * e_pot / Corr;
        e_t = s_t / S_L * e_pot / Corr;
        e_s = s_s / (h_b - S_L) * e_pot / Corr;
    }
    else
    {
        e_p = 0.0;
        e_t = 0.0;
        e_s = 0.0;
    }

    // fluxes precip->ponds and ponds-> links
    double q_ra_po = forcing_values[0] * c_1;             // (mm/hr -> m/min)
    double q_po_li = k_2 * s_p;                           // (1/min * m)
    double q_to_li = k_tl * s_t;                          // (1/min * m)

                                                          // fluxes precip->top and top->subsurface
    double q_ra_to = forcing_values[1] * c_1;             // (mm/hr -> m/min)
    double q_to_ss = k_i * s_t;

    // flux subsurface -> link
    double q_ss_li = k_3 * s_s;                           // (1/min * m)

                                                          //Discharge
    ans[0] = -q + (q_po_li + q_to_li + q_ss_li) * c_2;    // main streamflow
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];

    //Hillslope
    ans[1] = q_ra_po - q_po_li - e_p;                     // ponds
    ans[2] = q_ra_to - q_to_li - q_to_ss - e_t;           // top layer
    ans[3] = q_to_ss - q_ss_li - e_s;                     // sub surface

                                                          //Additional states
    ans[4] = q_ra_po + q_ra_to;                           // acc precip
    ans[5] = q_ra_po;                                     // acc runoff
    ans[6] = e_pot;                                       // acc evaporation
    ans[7] = q_ss_li * A_h - q_b*60.0;                    // base flow
    for (i = 0; i<num_parents; i++)
        ans[7] += y_p[i * dim + 6] * 60.0;
    ans[7] *= v_B / L;
}

//Type 259
//Contains 3 layers on hillslope: ponded, top layer, soil
//Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
//The numbering is:	0   1   2     3    4   5   6   7
//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent
//The numbering is:        0      1        2     3   4     5        6   7  8 9  10
void OfflineTopLayerInterflowHillslope_Reservoirs(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    ans[0] = forcing_values[3];
    ans[1] = 0.0;
    ans[2] = 0.0;
    ans[3] = 0.0;
}

//Type 400
//Tetis model structure for runoff generation + normal routing (NO stream order based velocity)
//Four layers
//Global parameters:
//The numbering is:    0   1   2   3      4      5   6   7...8

//y_i: vector with model states to be resolved by the solver
//dim:scalar with number of dimensions (states?) of the model
//y_p
//num_parents: number of tributary links
//max_dim
//global_params: global parameters applied to all hillslopes. See Precalculations in definitions.c
//params: distributed parameters per hillslope. see Precalculations in definitions.c
//forcing_values:
//qvs: unused
//state: unused
//user: unused
//ans: unused
void model400(double t, \
		const double * const y_i, \
		unsigned int dim, \
		const double * const y_p, \
		unsigned short num_parents, \
		unsigned int max_dim, \
		const double * const global_params, \
		const double * const params, \
		const double * const forcing_values, \
		const QVSData * const qvs, \
		int state, \
		void* user, \
		double *ans)
{

	 	unsigned short i; //auxiliary variable for loops
	    double L = params[1];   // Length of the channel [m]
	    double A_h = params[2]; //Area of the hillslopes [m^2]
	    double c_1 = params[4]; //factor .converts [mm/hr] to [m/min]
	    double rainfall = forcing_values[0] * c_1; //rainfall. from [mm/hr] to [m/min]
	    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));//potential et[mm/month] -> [m/min]
        double temperature = forcing_values[2]; //daily temperature in Celsius
        double temp_thres=global_params[10]; // celsius degrees
        double melt_factor = global_params[9] *(1/(24*60.0)) *(1/1000.0); // mm/day/degree to m/min/degree
        double frozen_ground = forcing_values[3]; // 1 if ground is frozen, 0 if not frozen 
        double x1 =0;

        //states
        unsigned int STATE_DISCHARGE=0;
        unsigned int STATE_STATIC= 1;
        unsigned int STATE_SURFACE=2;
        unsigned int STATE_SUBSURF=3;
        unsigned int STATE_GW = 4;
        unsigned int STATE_SNOW = 5;

        //INITIAL VALUES
        double h5 = y_i[STATE_SNOW];//snow storage [m]
		double h1 = y_i[STATE_STATIC]; //static storage [m]
		double h2 = y_i[STATE_SURFACE];//water in the hillslope surface [m]
		double h3 = y_i[STATE_SUBSURF]; //water in the gravitational storage in the upper part of soil [m]
		double h4 = y_i[STATE_GW]; //water in the aquifer storage [m]
	    double q = y_i[STATE_DISCHARGE];      //[m^3/s]

        //snow storage
        //temperature =0 is the flag for no forcing the variable. no snow process
        if(temperature==0){
            x1 = rainfall;
            ans[STATE_SNOW]=0;
        }
        else{
            if(temperature>=temp_thres){
                double snowmelt = min(h5,temperature * melt_factor); // in [m]
                ans[STATE_SNOW]=-snowmelt; //melting outs of snow storage
                x1 = rainfall + snowmelt; // in [m]
               // printf("temp > th: %f\n", temperature);
               // printf("snowmelt : %f\n", snowmelt);
            }
            if(temperature != 0 & temperature <temp_thres){
                ans[STATE_SNOW]=rainfall; //all precipitation is stored in the snow storage
                x1=0;
                //printf("temp < th: %f\n", temperature);
            }
        }
        

		//static storage
		double Hu = global_params[3]/1000; //max available storage in static tank [mm] to [m]
		double x2 = max(0,x1 + h1 - Hu ); //excedance flow to the second storage [m] [m/min] check units
        //if ground is frozen, x1 goes directly to the surface
        //therefore nothing is diverted to static tank
        if(frozen_ground == 1){
            x2 = x1;
        }
            
		double d1 = x1 - x2; // the input to static tank [m/min]
		double out1 = min(e_pot, h1); //evaporation from the static tank. it cannot evaporate more than h1 [m]
		//double out1 = (e_pot > h1) ? e_pot : 0.0;
		ans[STATE_STATIC] = d1 - out1; //differential equation of static storage


		//surface storage tank
		double infiltration = global_params[4]*c_1; //infiltration rate [m/min]
         if(frozen_ground == 1){
            infiltration = 0;
        }
		double x3 = min(x2, infiltration); //water that infiltrates to gravitational storage [m/min]
		double d2 = x2 - x3; // the input to surface storage [m] check units
        double alfa2 =global_params[6]; //velocity in m/s
        double w = alfa2 * L / A_h  * 60; // [1/min]
        w = min(1,w); //water can take less than 1 min (dt) to leave surface
        double out2 =0;
        out2  = h2 * w; //direct runoff [m/min]
		ans[STATE_SURFACE] = d2 - out2; //differential equation of surface storage


		// SUBSURFACE storage
		double percolation = global_params[5]*c_1; // percolation rate to aquifer [m/min]
		double x4 = min(x3,percolation); //water that percolates to aquifer storage [m/min]
		double d3 = x3 - x4; // input to gravitational storage [m/min]
		double alfa3 = global_params[7]* 24*60; //residence time [days] to [min].
        double out3=0;
        if(alfa3>=1)
		    out3 = h3/alfa3; //interflow [m/min]
		ans[STATE_SUBSURF] = d3 - out3; //differential equation for gravitational storage

		//aquifer storage
		double x5 = 0;//water loss to deeper aquifer [m]
		double d4 = x4 - x5;
		double alfa4 = global_params[8]* 24*60; //residence time [days] to [min].
        double out4=0;
        if(alfa4>=1)
		    out4 = h4/alfa4 ; //base flow [m/min]
		ans[STATE_GW] = d4 - out4; //differential equation for aquifer storage

		//channel storage

		double lambda_1 = global_params[1];
	    double invtau = params[3];// 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
	   	double c_2 = params[5];// = A_h / 60.0;	//  c_2

	    ans[STATE_DISCHARGE] = -q + (out2 + out3 + out4) * c_2; //[m/min] to [m3/s]
	    for (i = 0; i < num_parents; i++)
	        ans[STATE_DISCHARGE] += y_p[i * dim + STATE_DISCHARGE];
	    ans[STATE_DISCHARGE] = invtau * pow(q, lambda_1) * ans[STATE_DISCHARGE];    // discharge[0]

        // if (forcing_values[0]>1 && ratio<1) {
        //     printf("time: %f\n", t);
        //     printf(" rain in mm/hour: %f\n", forcing_values[0]);
        //     printf(" area hill, area basin, area ratio: %f %f %f\n", A_h,A_i,ratio);
        //     MPI_Abort(MPI_COMM_WORLD, 1);
        // }

}
//Type 401
//model 400 + tracking of total runoff, surface runoff, subsurface runoff and groundwater runoff
void model401(double t, \
		const double * const y_i, \
		unsigned int dim, \
		const double * const y_p, \
		unsigned short num_parents, \
		unsigned int max_dim, \
		const double * const global_params, \
		const double * const params, \
		const double * const forcing_values, \
		const QVSData * const qvs, \
		int state, \
		void* user, \
		double *ans)
{

	 	unsigned short i; //auxiliary variable for loops
	    double L = params[1];   // Length of the channel [m]
	    double A_h = params[2]; //Area of the hillslopes [m^2]
        double A_i = params[0]; //Area of the basin [km^2]
         A_i*=1e6; //[m^2]
	    double c_1 = params[4]; //factor .converts [mm/hr] to [m/min]
	    double rainfall = forcing_values[0] * c_1; //rainfall. from [mm/hr] to [m/min]
	    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));//potential et[mm/month] -> [m/min]
		double temperature = forcing_values[2]; //daily temperature in Celsius
        double temp_thres=global_params[10]; // celsius degrees
        double melt_factor = global_params[9] *(1/(24*60.0)) *(1/1000.0); // mm/day/degree to m/min/degree
        double frozen_ground = forcing_values[3]; // 1 if ground is frozen, 0 if not frozen 
        double x1 =0;

        // 9 states
        // i need to put the fluxes on top because cant print more than state7. bug
        //y0=q discharge[m3/s]
        //y1 = basin rainfall and snowmelt [m3/hour]
        //y2 = basin surface runoff [m3/hour]
        //y3= basin subsurface runoff [m3/hour]
        //y4 = basin gw rounoff [m3/hour]
        //y5= h1 static storage[m]
        //y6= h2 water hill surface[m]
        //y7 = h3 water upper soil [m]
        //y8 = h4  water lower soil [m]
        //y9 = h5 snow storage [m]
        unsigned int STATE_DISCHARGE=0;
        unsigned int STATE_CUMRAINFALL=1;
        unsigned int STATE_CUMSURF = 2;
        unsigned int STATE_CUMSUB=3;
        unsigned int STATE_CUMGW = 4;
        unsigned int STATE_STATIC=5;
        unsigned int STATE_SURFACE=6;
        unsigned int STATE_SUBSURF =7;
        unsigned int STATE_GW = 8;
        unsigned int STATE_SNOW = 9;

        #ifndef minf
        #define minf(a,b) ((a) < (b) ? (a) : (b))
        #endif

        #ifndef maxf
        #define maxf(a,b) ((a) > (b) ? (a) : (b))
        #endif

        //INITIAL VALUES
        double h5 = y_i[STATE_SNOW];//snow storage [m]
        double basin_rainfall = y_i[STATE_CUMRAINFALL]; //[m3/hour] see last lines
		double h1 = y_i[STATE_STATIC]; //static storage [m]
		double h2 = y_i[STATE_SURFACE];//water in the hillslope surface [m]
		double h3 = y_i[STATE_SUBSURF]; //water in the gravitational storage in the upper part of soil [m]
		double h4 = y_i[STATE_GW]; //water in the aquifer storage [m]
	    double q = y_i[STATE_DISCHARGE];      //[m^3/s]

        
        //snow storage
        //temperature =0 is the flag for no forcing the variable. no snow process
        if(temperature==0){
            x1 = rainfall;
            ans[STATE_SNOW]=0;
        }
        else{
            if(temperature>=temp_thres){
                double snowmelt = min(h5,temperature * melt_factor); // in [m/min]
                ans[STATE_SNOW]=-snowmelt; //melting outs of snow storage
                x1 = rainfall + snowmelt; // in [m/min]
                //printf("temp > th: %f\n", temperature);
                //printf("snowmelt : %f\n", snowmelt);
            }
            if(temperature != 0 & temperature <temp_thres){
                ans[STATE_SNOW]=rainfall; //all precipitation is stored in the snow storage
                x1=0;
                //printf("temp < th: %f\n", temperature);
            }
        }

		//static storage
		double Hu = global_params[3]/1000; //max available storage in static tank [mm] to [m]
        double aux1 = x1 + h1 - Hu;
		double x2 = max(0,aux1); //excedance flow to the second storage [m] [m/min] check units
        //double x2 = (aux1>0)? aux1: 0.0;
        //if ground is frozen, x1 goes directly to the surface
        //therefore nothing is diverted to static tank
        if(frozen_ground == 1){
            x2 = x1;
        }
		double d1 = x1 - x2; // the input to static tank [m/min]
		double out1 = min(e_pot, h1); //evaporation from the static tank. it cannot evaporate more than h1 [m]
		//double out1 = (e_pot > h1) ? e_pot : 0.0;
		ans[STATE_STATIC] = d1 - out1; //differential equation of static storage
        //printf("t %f\n",t);
        //printf("d1 %f\n",d1);
        //printf("x1 and x2 %f %f\n",x1,x2);
        //printf("out1 %f\n", out1);
        //printf(" rain in mm/hour: %f\n", forcing_values[0]);
        //printf("h1 %f\n ",h1);
        //MPI_Abort(MPI_COMM_WORLD, 1);

		//surface storage tank
		double infiltration = global_params[4]*c_1; //infiltration rate [m/min]
        if(frozen_ground == 1){
            infiltration = 0;
        }
		double x3 = min(x2, infiltration); //water that infiltrates to gravitational storage [m/min]
		double d2 = x2 - x3; // the input to surface storage [m] check units
		//double alfa2 = global_params[6]* 24*60; //residence time [days] to [min].
		double alfa2 =global_params[6]; //velocity in m/s
        double w = alfa2 * L / A_h  * 60; // [1/min]
        w = min(w,1); // water can take less than 1 min to
        double out2 =0;
        out2  = h2 * w; //direct runoff [m/min]
		ans[STATE_SURFACE] = d2 - out2; //differential equation of surface storage
        double surface_runoff = y_i[STATE_CUMSURF]; //[m3/hour]

		// SUBSURFACE storage
		double percolation = global_params[5]*c_1; // percolation rate to aquifer [m/min]
		double x4 = min(x3,percolation); //water that percolates to aquifer storage [m/min]
		double d3 = x3 - x4; // input to gravitational storage [m/min]
		double alfa3 = global_params[7]* 24*60; //residence time [days] to [min].
		double out3=0;
        if(alfa3>=1)
		    out3 = h3/alfa3; //interflow [m/min]
		ans[STATE_SUBSURF] = d3 - out3; //differential equation for gravitational storage
        double subsurface_runoff = y_i[STATE_CUMSUB];  //[m3/hour]

		//aquifer storage
		double x5 = 0;//water loss to deeper aquifer [m]
		double d4 = x4 - x5;
		double alfa4 = global_params[8]* 24*60; //residence time [days] to [min].
		double out4=0;
        if(alfa4>=1)
		    out4 = h4/alfa4 ; //base flow [m/min]
		ans[STATE_GW] = d4 - out4; //differential equation for aquifer storage
        double groundwater_runoff = y_i[STATE_CUMGW]; //[m3/hour]

		//channel storage
		double lambda_1 = global_params[1];
	    double invtau = params[3];// 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
	   	double c_2 = params[5];// = A_h / 60.0;	//  c_2

	    ans[STATE_DISCHARGE] = -q + (out2 + out3 + out4) * c_2; //[m/min] to [m3/s]

        double aux = forcing_values[0] *(1/1000.0) * A_h;//[mm/h] to [m3/h] 
        ans[STATE_CUMRAINFALL] = -basin_rainfall + aux; //[m3/hour] 
        ans[STATE_CUMSURF] = -surface_runoff + out2 * 60.0 * A_h; //[m/min] to [m3/hour]
        ans[STATE_CUMSUB] = -subsurface_runoff + out3 *60.0 * A_h ; //[m/min] to [m3/hour]
        ans[STATE_CUMGW] = -groundwater_runoff + out4 *60.0 *A_h ; //[m/min] to [m3/hour]

	    for (i = 0; i < num_parents; i++){
            ans[STATE_DISCHARGE] += y_p[i * dim + STATE_DISCHARGE];
            ans[STATE_CUMRAINFALL] += y_p[i * dim +STATE_CUMRAINFALL];
            ans[STATE_CUMSURF] += y_p[i * dim +STATE_CUMSURF];
            ans[STATE_CUMSUB] += y_p[i * dim +STATE_CUMSUB];
            ans[STATE_CUMGW] += y_p[i * dim +STATE_CUMGW];
        }
	        
	    ans[STATE_DISCHARGE] = invtau * pow(q, lambda_1) * ans[STATE_DISCHARGE];    // discharge[0]
}
void model401reservoir(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{	
    unsigned int FORCING_RESERVOIR = 4;
    unsigned int STATE_DISCHARGE=0;
	if(forcing_values[FORCING_RESERVOIR] >=0){
		ans[STATE_DISCHARGE] = forcing_values[FORCING_RESERVOIR];
	}
}

//Type 402
void model402(double t, \
		const double * const y_i, \
		unsigned int dim, \
		const double * const y_p, \
		unsigned short num_parents, \
		unsigned int max_dim, \
		const double * const global_params, \
		const double * const params, \
		const double * const forcing_values, \
		const QVSData * const qvs, \
		int state, \
		void* user, \
		double *ans)
{

	 	unsigned short i; //auxiliary variable for loops
	    double L = params[1];   // Length of the channel [m]
	    double A_h = params[2]; //Area of the hillslopes [m^2]
	    double c_1 = params[4]; //factor .converts [mm/hr] to [m/min]
	    double rainfall = forcing_values[0] * c_1; //rainfall. from [mm/hr] to [m/min]
	    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));//potential et[mm/month] -> [m/min]
        double temperature = forcing_values[2]; //daily temperature in Celsius
        double temp_thres=global_params[10]; // celsius degrees
        double melt_factor = global_params[9] *(1/(24*60.0)) *(1/1000.0); // mm/day/degree to m/min/degree
        double frozen_ground = forcing_values[3]; // 1 if ground is frozen, 0 if not frozen 
        double x1 =0;

        //states
        unsigned int STATE_DISCHARGE=0;
        unsigned int STATE_STATIC= 1;
        unsigned int STATE_SURFACE=2;
        unsigned int STATE_SUBSURF=3;
        unsigned int STATE_GW = 4;
        unsigned int STATE_SNOW = 5;
        unsigned int STATE_DAM_STORAGE=6;

        //INITIAL VALUES
        double h5 = y_i[STATE_SNOW];//snow storage [m]
		double h1 = y_i[STATE_STATIC]; //static storage [m]
		double h2 = y_i[STATE_SURFACE];//water in the hillslope surface [m]
		double h3 = y_i[STATE_SUBSURF]; //water in the gravitational storage in the upper part of soil [m]
		double h4 = y_i[STATE_GW]; //water in the aquifer storage [m]
	    double q = y_i[STATE_DISCHARGE];      //[m^3/s]
        double dam_storage =y_i[STATE_DAM_STORAGE] ;//m^3

        //snow storage
        //temperature =0 is the flag for no forcing the variable. no snow process
        if(temperature==0){
            x1 = rainfall;
            ans[STATE_SNOW]=0;
        }
        else{
            if(temperature>=temp_thres){
                double snowmelt = min(h5,temperature * melt_factor); // in [m]
                ans[STATE_SNOW]=-snowmelt; //melting outs of snow storage
                x1 = rainfall + snowmelt; // in [m]
               // printf("temp > th: %f\n", temperature);
               // printf("snowmelt : %f\n", snowmelt);
            }
            if(temperature != 0 & temperature <temp_thres){
                ans[STATE_SNOW]=rainfall; //all precipitation is stored in the snow storage
                x1=0;
                //printf("temp < th: %f\n", temperature);
            }
        }
        

		//static storage
		double Hu = global_params[3]/1000; //max available storage in static tank [mm] to [m]
		double x2 = max(0,x1 + h1 - Hu ); //excedance flow to the second storage [m] [m/min] check units
        //if ground is frozen, x1 goes directly to the surface
        //therefore nothing is diverted to static tank
        if(frozen_ground == 1){
            x2 = x1;
        }
            
		double d1 = x1 - x2; // the input to static tank [m/min]
		double out1 = min(e_pot, h1); //evaporation from the static tank. it cannot evaporate more than h1 [m]
		//double out1 = (e_pot > h1) ? e_pot : 0.0;
		ans[STATE_STATIC] = d1 - out1; //differential equation of static storage


		//surface storage tank
		double infiltration = global_params[4]*c_1; //infiltration rate [m/min]
         if(frozen_ground == 1){
            infiltration = 0;
        }
		double x3 = min(x2, infiltration); //water that infiltrates to gravitational storage [m/min]
		double d2 = x2 - x3; // the input to surface storage [m] check units
        double alfa2 =global_params[6]; //velocity in m/s
        double w = alfa2 * L / A_h  * 60; // [1/min]
        w = min(1,w); //water can take less than 1 min (dt) to leave surface
        double out2 =0;
        out2  = h2 * w; //direct runoff [m/min]
        //control surface storage to not become too small. if not, calculation takes longer
        if(h2<1E-3)
            h2 = 1E-2; // this wont make the model run faster or slower
		ans[STATE_SURFACE] = d2 - out2; //differential equation of surface storage


		// SUBSURFACE storage
		double percolation = global_params[5]*c_1; // percolation rate to aquifer [m/min]
		double x4 = min(x3,percolation); //water that percolates to aquifer storage [m/min]
		double d3 = x3 - x4; // input to gravitational storage [m/min]
		double alfa3 = global_params[7]* 24*60; //residence time [days] to [min].
        double out3=0;
        if(alfa3>=1)
		    out3 = h3/alfa3; //interflow [m/min]
		ans[STATE_SUBSURF] = d3 - out3; //differential equation for gravitational storage

		//aquifer storage
		double x5 = 0;//water loss to deeper aquifer [m]
		double d4 = x4 - x5;
		double alfa4 = global_params[8]* 24*60; //residence time [days] to [min].
        double out4=0;
        if(alfa4>=1)
		    out4 = h4/alfa4 ; //base flow [m/min]
		ans[STATE_GW] = d4 - out4; //differential equation for aquifer storage

		//channel storage

		double lambda_1 = global_params[1];
	    double invtau = params[3];// 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
	   	double c_2 = params[5];// = A_h / 60.0;	//  c_2

        if(state==0){ //no dams
            ans[STATE_DISCHARGE] = -q + (out2 + out3 + out4) * c_2; //[m/min] to [m3/s]
	        for (i = 0; i < num_parents; i++)
	            ans[STATE_DISCHARGE] += y_p[i * dim + STATE_DISCHARGE];
	        ans[STATE_DISCHARGE] = invtau * pow(q, lambda_1) * ans[STATE_DISCHARGE];    // discharge[0]
            //this lines were intented to keep storage at zero on no dam links, but is making the
           dam_storage=1e-1; // this is the model making run slower
           ans[STATE_DAM_STORAGE] = -dam_storage; // this is the model making run slower

        }
        if(state ==1){// dams
            double dam_input = 0;//m3
            double dam_output=0;//m3
            double dam_outflow=0; //m3s-1
            int debug = 0;
            if(debug) printf("time: %f\n", t);
            //inflow from upstream links
            for (i = 0; i < num_parents; i++)
	            dam_input += y_p[i * dim + STATE_DISCHARGE] * 60.0; //m3/s to m3
            if(debug) printf("inflow from upstream in m3: %f\n", dam_input);
            // storage inputs
            if(debug) printf("initial dam storage in m3: %f\n", dam_storage);
            //dam_storage= dam_storage + dam_input;
            if(debug) printf("step1 dam storage in m3: %f\n", dam_storage + dam_input);

            //calculate dam outflow
            dam_outflow = get_discharge_from_storage(dam_storage + dam_input); //m3s-1
            if(debug) printf("dam outflow in m3/s: %f\n", dam_outflow);

            ans[STATE_DISCHARGE] = dam_outflow - q;
            
            //storage output
            //dam_output =get_storage_from_discharge(dam_outflow);
            dam_output = dam_outflow * 60.0; //m3s-1 to m3
            if(debug) printf("dam output storage in m3: %f\n", dam_output);

            //dam_storage =dam_storage - dam_output;
            if(debug) printf("step2 dam storage in m3: %f\n", dam_storage + dam_input - dam_output);

            //dam_storage = max(dam_storage,0);
            //if(debug) printf("step3 dam storage in m3: %f\n", dam_storage);
            //ans[STATE_DAM_STORAGE] = dam_storage;
            ans[STATE_DAM_STORAGE] = dam_input - dam_output;

            if(debug) printf("\n");

            
        }
}

double get_discharge_from_storage(double storage){
    double discharge=0;
    if(storage<125000) //discharge = 1m3s-1
        discharge = (1-0)/(125000-0.0)*storage;
     if(storage>=125000)
        discharge = ((100.0-1.0)/(225000.0-125000.0))*(storage - 125000.0)+1.0;
    return discharge;
}

double get_storage_from_discharge(double discharge){
    double storage=0;
    if(discharge<1)
        storage = (125000-0)/(1-0)*discharge ;
    else
        storage = (225000-125000)/(100.0-1.0)*(discharge-1.0) + 125000; 
    return storage;
}

//Type 403
//model tetis with dams using specific individual ratings 
void model403(double t, \
		const double * const y_i, \
		unsigned int dim, \
		const double * const y_p, \
		unsigned short num_parents, \
		unsigned int max_dim, \
		const double * const global_params, \
		const double * const params, \
		const double * const forcing_values, \
		const QVSData * const qvs, \
		int state, \
		void* user, \
		double *ans)
{

	 	unsigned short i; //auxiliary variable for loops
	    double L = params[1];   // Length of the channel [m]
	    double A_h = params[2]; //Area of the hillslopes [m^2]
	    double c_1 = params[4]; //factor .converts [mm/hr] to [m/min]
	    double rainfall = forcing_values[0] * c_1; //rainfall. from [mm/hr] to [m/min]
	    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));//potential et[mm/month] -> [m/min]
        double temperature = forcing_values[2]; //daily temperature in Celsius
        double temp_thres=global_params[10]; // celsius degrees
        double melt_factor = global_params[9] *(1/(24*60.0)) *(1/1000.0); // mm/day/degree to m/min/degree
        double frozen_ground = forcing_values[3]; // 1 if ground is frozen, 0 if not frozen 
        double x1 =0;

        //states
        unsigned int STATE_DISCHARGE=0;
        unsigned int STATE_STATIC= 1;
        unsigned int STATE_SURFACE=2;
        unsigned int STATE_SUBSURF=3;
        unsigned int STATE_GW = 4;
        unsigned int STATE_SNOW = 5;
        unsigned int STATE_DAM_STORAGE=6;

        //INITIAL VALUES
        double h5 = y_i[STATE_SNOW];//snow storage [m]
		double h1 = y_i[STATE_STATIC]; //static storage [m]
		double h2 = y_i[STATE_SURFACE];//water in the hillslope surface [m]
		double h3 = y_i[STATE_SUBSURF]; //water in the gravitational storage in the upper part of soil [m]
		double h4 = y_i[STATE_GW]; //water in the aquifer storage [m]
	    double q = y_i[STATE_DISCHARGE];      //[m^3/s]
        double dam_storage =y_i[STATE_DAM_STORAGE] ;//m^3

        //snow storage
        //temperature =0 is the flag for no forcing the variable. no snow process
        if(temperature==0){
            x1 = rainfall;
            ans[STATE_SNOW]=0;
        }
        else{
            if(temperature>=temp_thres){
                double snowmelt = min(h5,temperature * melt_factor); // in [m]
                ans[STATE_SNOW]=-snowmelt; //melting outs of snow storage
                x1 = rainfall + snowmelt; // in [m]
               // printf("temp > th: %f\n", temperature);
               // printf("snowmelt : %f\n", snowmelt);
            }
            if(temperature != 0 & temperature <temp_thres){
                ans[STATE_SNOW]=rainfall; //all precipitation is stored in the snow storage
                x1=0;
                //printf("temp < th: %f\n", temperature);
            }
        }
        

		//static storage
		double Hu = global_params[3]/1000; //max available storage in static tank [mm] to [m]
		double x2 = max(0,x1 + h1 - Hu ); //excedance flow to the second storage [m] [m/min] check units
        //if ground is frozen, x1 goes directly to the surface
        //therefore nothing is diverted to static tank
        if(frozen_ground == 1){
            x2 = x1;
        }
            
		double d1 = x1 - x2; // the input to static tank [m/min]
		double out1 = min(e_pot, h1); //evaporation from the static tank. it cannot evaporate more than h1 [m]
		//double out1 = (e_pot > h1) ? e_pot : 0.0;
		ans[STATE_STATIC] = d1 - out1; //differential equation of static storage


		//surface storage tank
		double infiltration = global_params[4]*c_1; //infiltration rate [m/min]
         if(frozen_ground == 1){
            infiltration = 0;
        }
		double x3 = min(x2, infiltration); //water that infiltrates to gravitational storage [m/min]
		double d2 = x2 - x3; // the input to surface storage [m] check units
        double alfa2 =global_params[6]; //velocity in m/s
        double w = alfa2 * L / A_h  * 60; // [1/min]
        w = min(1,w); //water can take less than 1 min (dt) to leave surface
        double out2 =0;
        out2  = h2 * w; //direct runoff [m/min]
		ans[STATE_SURFACE] = d2 - out2; //differential equation of surface storage


		// SUBSURFACE storage
		double percolation = global_params[5]*c_1; // percolation rate to aquifer [m/min]
		double x4 = min(x3,percolation); //water that percolates to aquifer storage [m/min]
		double d3 = x3 - x4; // input to gravitational storage [m/min]
		double alfa3 = global_params[7]* 24*60; //residence time [days] to [min].
        double out3=0;
        if(alfa3>=1)
		    out3 = h3/alfa3; //interflow [m/min]
		ans[STATE_SUBSURF] = d3 - out3; //differential equation for gravitational storage

		//aquifer storage
		double x5 = 0;//water loss to deeper aquifer [m]
		double d4 = x4 - x5;
		double alfa4 = global_params[8]* 24*60; //residence time [days] to [min].
        double out4=0;
        if(alfa4>=1)
		    out4 = h4/alfa4 ; //base flow [m/min]
		ans[STATE_GW] = d4 - out4; //differential equation for aquifer storage

		//channel storage

		double lambda_1 = global_params[1];
	    double invtau = params[3];// 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
	   	double c_2 = params[5];// = A_h / 60.0;	//  c_2

        //printf("state: %d\n", state);
        if(state == 0){ //no dams
            ans[STATE_DISCHARGE] = -q + (out2 + out3 + out4) * c_2; //[m/min] to [m3/s]
	        for (i = 0; i < num_parents; i++)
	            ans[STATE_DISCHARGE] += y_p[i * dim + STATE_DISCHARGE];
	        ans[STATE_DISCHARGE] = invtau * pow(q, lambda_1) * ans[STATE_DISCHARGE];    // discharge[0]
            dam_storage=1e-1;
            ans[STATE_DAM_STORAGE] = -dam_storage;
        }
        //state is the array index corresponding to the current storage in the qvs 
        if(state == 1 ){// dams
            int debug = 0;
            
            if(debug) printf("state: %d\n", state);

            double dam_input = 0;//m3
            double dam_output=0;//m3
            double dam_outflow=0; //m3s-1
            
            
            if(debug) printf("time: %f\n", t);
            //inflow from upstream links
            for (i = 0; i < num_parents; i++)
	            dam_input += y_p[i * dim + STATE_DISCHARGE] * 60.0; //m3/s to m3
            if(debug) printf("inflow from upstream in m3: %f\n", dam_input);
            // storage inputs
            if(debug) printf("initial dam storage in m3: %f\n", dam_storage);
            //dam_storage= dam_storage + dam_input;
            if(debug) printf("step1 dam storage in m3: %f\n", dam_storage + dam_input);

            //if(debug) printf("state: %u\n", state);
            unsigned int max_storage_pond_index = (int)qvs->n_values - 1;
            double max_storage_pond = qvs->points[max_storage_pond_index][0];
            if (dam_storage + dam_input >= max_storage_pond)
            {
                if(debug) printf("storage >= max.sto.pond\n");
                //S_max = qvs->points[qvs->n_values - 1][0];
                double q_max = qvs->points[max_storage_pond_index][1];
                if(debug) printf("qmax : %f\n", q_max);
                dam_outflow = q_max;//m3s-1
            }
            else if (dam_storage + dam_input < max_storage_pond)
            {
                if(debug) printf("storage within rating range\n");

                //S = (y_i[1] < 0.0) ? 0.0 : y_i[1];
                double S = (dam_storage + dam_input < 0.0) ? 0.0 : dam_storage + dam_input;
                if(debug) printf("storage  : %f\n", S);

                //find discharge value in rating
                unsigned int ii;
                for (ii = 0; ii < qvs->n_values - 1; ii++){
				   // if(debug)printf("model403 storage qvs points: %f\n", qvs->points[ii][0]);
                    if (qvs->points[ii][0] <= S
						&& S < qvs->points[ii + 1][0])    
					break;
                }
                double q2 = qvs->points[ii + 1][1];
                double q1 = qvs->points[ii][1];
                double S2 = qvs->points[ii + 1][0];
                double S1 = qvs->points[ii][0];
                dam_outflow = (q2 - q1) / (S2 - S1) * (S - S1) + q1;//m3s-1
                
            }
            
            //calculate dam outflow
            //dam_outflow = get_discharge_from_storage(dam_storage + dam_input); //m3s-1
            if(debug) printf("dam outflow in m3/s: %f\n", dam_outflow);

            ans[STATE_DISCHARGE] = dam_outflow - q;
            
            //storage output
            //dam_output =get_storage_from_discharge(dam_outflow);
            dam_output = dam_outflow * 60.0; //m3s-1 to m3
            if(debug) printf("dam output storage in m3: %f\n", dam_output);

        
            if(debug) printf("step2 dam storage in m3: %f\n", dam_storage + dam_input - dam_output);

            //dam_storage = max(dam_storage,0);
            //if(debug) printf("step3 dam storage in m3: %f\n", dam_storage);
            //ans[STATE_DAM_STORAGE] = dam_storage;
            ans[STATE_DAM_STORAGE] = dam_input - dam_output;

            if(debug) printf("\n");

            
        }

        // if (forcing_values[0]>1 && ratio<1) {
        //     printf("time: %f\n", t);
        //     printf(" rain in mm/hour: %f\n", forcing_values[0]);
        //     printf(" area hill, area basin, area ratio: %f %f %f\n", A_h,A_i,ratio);
        //     MPI_Abort(MPI_COMM_WORLD, 1);
        // }
}

//Type 402
void Tetis03_Reservoirs(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    ans[0] = forcing_values[3];
    ans[1] = 0.0;
    ans[2] = 0.0;
    ans[3] = 0.0;
    ans[4] = 0.0;
    ans[5] = 0.0;
    ans[6] = 0.0;
    ans[7] = 0.0;
    ans[8] = 0.0;
}

//Type 404
//Tetis model structure for runoff generation + normal routing (NO stream order based velocity)
//Four layers
//Global parameters:
//The numbering is:    0   1   2   3      4      5   6   7...8

//y_i: vector with model states to be resolved by the solver
//dim:scalar with number of dimensions (states?) of the model
//y_p
//num_parents: number of tributary links
//max_dim
//global_params: global parameters applied to all hillslopes. See Precalculations in definitions.c
//params: distributed parameters per hillslope. see Precalculations in definitions.c
//forcing_values:
//qvs: unused
//state: unused
//user: unused
//ans: unused
void model404(double t, \
		const double * const y_i, \
		unsigned int dim, \
		const double * const y_p, \
		unsigned short num_parents, \
		unsigned int max_dim, \
		const double * const global_params, \
		const double * const params, \
		const double * const forcing_values, \
		const QVSData * const qvs, \
		int state, \
		void* user, \
		double *ans)
{

	 	unsigned short i; //auxiliary variable for loops

        double c_1 = (0.001 / 60.0);		//(mm/hr->m/min)  c_1
        double A_i = params[0]; //drainage area in km2
	    double L = params[1];   // Length of the channel [m]
	    double A_h = params[2]; //Area of the hillslopes [m^2]
        double c_2 = A_h / 60.0;	//  c_2
        double v_0 = params[3];
        double lambda_1 = params[4];
        double lambda_2 = params[5];
        double Hu = params[6]/1000; //[m]
        double infiltration = params[7]*c_1; //infiltration rate [m/min]
		double percolation = params[8]*c_1; // percolation rate to aquifer [m/min]
        double alfa2 =params[9]; //velocity in m/s
		double alfa3 = params[10]* 24*60; //residence time [days] to [min].
		double alfa4 = params[11]* 24*60; //residence time [days] to [min].
		double melt_factor = params[12] *(1/(24*60.0)) *(1/1000.0); // mm/day/degree to m/min/degree
        double temp_thres= params[13]; // celsius degrees
        double invtau = 60.0 * v_0 * pow(A_i, lambda_2) / ((1.0 - lambda_1) * L);//[1/min]  invtau params[3]


        double x1 =0;


        //forcings
	    double rainfall = forcing_values[0] * c_1; //rainfall. from [mm/hr] to [m/min]
	    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));//potential et[mm/month] -> [m/min]
        double temperature = forcing_values[2]; //daily temperature in Celsius
        double frozen_ground = forcing_values[3]; // 1 if ground is frozen, 0 if not frozen 
        
        //states
        unsigned int STATE_DISCHARGE=0;
        unsigned int STATE_STATIC= 1;
        unsigned int STATE_SURFACE=2;
        unsigned int STATE_SUBSURF=3;
        unsigned int STATE_GW = 4;
        unsigned int STATE_SNOW = 5;

        //INITIAL VALUES
        double h5 = y_i[STATE_SNOW];//snow storage [m]
		double h1 = y_i[STATE_STATIC]; //static storage [m]
		double h2 = y_i[STATE_SURFACE];//water in the hillslope surface [m]
		double h3 = y_i[STATE_SUBSURF]; //water in the gravitational storage in the upper part of soil [m]
		double h4 = y_i[STATE_GW]; //water in the aquifer storage [m]
	    double q = y_i[STATE_DISCHARGE];      //[m^3/s]

        //snow storage
        //temperature =0 is the flag for no forcing the variable. no snow process
        if(temperature==0){
            x1 = rainfall;
            ans[STATE_SNOW]=0;
        }
        else{
            if(temperature>=temp_thres){
                double snowmelt = min(h5,temperature * melt_factor); // in [m]
                ans[STATE_SNOW]=-snowmelt; //melting outs of snow storage
                x1 = rainfall + snowmelt; // in [m]
               // printf("temp > th: %f\n", temperature);
               // printf("snowmelt : %f\n", snowmelt);
            }
            if(temperature != 0 & temperature <temp_thres){
                ans[STATE_SNOW]=rainfall; //all precipitation is stored in the snow storage
                x1=0;
                //printf("temp < th: %f\n", temperature);
            }
        }
        

		//static storage
		 //max available storage in static tank [mm] to [m]
		double x2 = max(0,x1 + h1 - Hu ); //excedance flow to the second storage [m] [m/min] check units
        //if ground is frozen, x1 goes directly to the surface
        //therefore nothing is diverted to static tank
        if(frozen_ground == 1){
            x2 = x1;
        }
            
		double d1 = x1 - x2; // the input to static tank [m/min]
		double out1 = min(e_pot, h1); //evaporation from the static tank. it cannot evaporate more than h1 [m]
		//double out1 = (e_pot > h1) ? e_pot : 0.0;
		ans[STATE_STATIC] = d1 - out1; //differential equation of static storage


		//surface storage tank
         if(frozen_ground == 1){
            infiltration = 0;
        }
		double x3 = min(x2, infiltration); //water that infiltrates to gravitational storage [m/min]
		double d2 = x2 - x3; // the input to surface storage [m] check units
        //double alfa2 =global_params[6]; //velocity in m/s
        
        
        //double alfa2 = (pow(h2,2./3.) * pow(slope,0.5))/manning; //m/s

        double w = alfa2 * L / A_h  * 60; // [1/min]
        w = min(1,w); //water can take less than 1 min (dt) to leave surface
        double out2 =0;
        out2  = h2 * w; //direct runoff [m/min]
		ans[STATE_SURFACE] = d2 - out2; //differential equation of surface storage


		// SUBSURFACE storage
		double x4 = min(x3,percolation); //water that percolates to aquifer storage [m/min]
		double d3 = x3 - x4; // input to gravitational storage [m/min]
        double out3=0;
        if(alfa3>=1)
		    out3 = h3/alfa3; //interflow [m/min]
		ans[STATE_SUBSURF] = d3 - out3; //differential equation for gravitational storage

		//aquifer storage
		double x5 = 0;//water loss to deeper aquifer [m]
		double d4 = x4 - x5;
        double out4=0;
        if(alfa4>=1)
		    out4 = h4/alfa4 ; //base flow [m/min]
		ans[STATE_GW] = d4 - out4; //differential equation for aquifer storage

		//channel storage

		//double lambda_1 = params[1];
	    //double invtau = params[3];// 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
	   	//double c_2 = params[5];// = A_h / 60.0;	//  c_2

	    ans[STATE_DISCHARGE] = -q + (out2 + out3 + out4) * c_2; //[m/min] to [m3/s]
	    for (i = 0; i < num_parents; i++)
	        ans[STATE_DISCHARGE] += y_p[i * dim + STATE_DISCHARGE];
	    ans[STATE_DISCHARGE] = invtau * pow(q, lambda_1) * ans[STATE_DISCHARGE];    // discharge[0]

        // if (forcing_values[0]>1 && ratio<1) {
        //     printf("time: %f\n", t);
        //     printf(" rain in mm/hour: %f\n", forcing_values[0]);
        //     printf(" area hill, area basin, area ratio: %f %f %f\n", A_h,A_i,ratio);
        //     MPI_Abort(MPI_COMM_WORLD, 1);
        // }

}

//Type 405
//model tetis with dams using specific individual ratings, but using the flowchart algorithm . 
void model405(double t, \
		const double * const y_i, \
		unsigned int dim, \
		const double * const y_p, \
		unsigned short num_parents, \
		unsigned int max_dim, \
		const double * const global_params, \
		const double * const params, \
		const double * const forcing_values, \
		const QVSData * const qvs, \
		int state, \
		void* user, \
		double *ans)
{

	 	unsigned short i; //auxiliary variable for loops
	    double L = params[1];   // Length of the channel [m]
	    double A_h = params[2]; //Area of the hillslopes [m^2]
	    double c_1 = params[4]; //factor .converts [mm/hr] to [m/min]
	    double rainfall = forcing_values[0] * c_1; //rainfall. from [mm/hr] to [m/min]
	    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));//potential et[mm/month] -> [m/min]
        double temperature = forcing_values[2]; //daily temperature in Celsius
        double temp_thres=global_params[10]; // celsius degrees
        double melt_factor = global_params[9] *(1/(24*60.0)) *(1/1000.0); // mm/day/degree to m/min/degree
        double frozen_ground = forcing_values[3]; // 1 if ground is frozen, 0 if not frozen 
        double x1 =0;
        double factor_low_threshold = global_params[11];
        double factor_high_threshold = global_params[12];
        double factor_operating_window = global_params[13];

        //states
        unsigned int STATE_DISCHARGE=0;
        unsigned int STATE_STATIC= 1;
        unsigned int STATE_SURFACE=2;
        unsigned int STATE_SUBSURF=3;
        unsigned int STATE_GW = 4;
        unsigned int STATE_SNOW = 5;
        unsigned int STATE_DAM_STORAGE=6;

        //INITIAL VALUES
        double h5 = y_i[STATE_SNOW];//snow storage [m]
		double h1 = y_i[STATE_STATIC]; //static storage [m]
		double h2 = y_i[STATE_SURFACE];//water in the hillslope surface [m]
		double h3 = y_i[STATE_SUBSURF]; //water in the gravitational storage in the upper part of soil [m]
		double h4 = y_i[STATE_GW]; //water in the aquifer storage [m]
	    double q = y_i[STATE_DISCHARGE];      //[m^3/s]
        double dam_storage =y_i[STATE_DAM_STORAGE] ;//m^3

        //snow storage
        //temperature =0 is the flag for no forcing the variable. no snow process
        if(temperature==0){
            x1 = rainfall;
            ans[STATE_SNOW]=0;
        }
        else{
            if(temperature>=temp_thres){
                double snowmelt = min(h5,temperature * melt_factor); // in [m]
                ans[STATE_SNOW]=-snowmelt; //melting outs of snow storage
                x1 = rainfall + snowmelt; // in [m]
               // printf("temp > th: %f\n", temperature);
               // printf("snowmelt : %f\n", snowmelt);
            }
            if(temperature != 0 & temperature <temp_thres){
                ans[STATE_SNOW]=rainfall; //all precipitation is stored in the snow storage
                x1=0;
                //printf("temp < th: %f\n", temperature);
            }
        }
        

		//static storage
		double Hu = global_params[3]/1000; //max available storage in static tank [mm] to [m]
		double x2 = max(0,x1 + h1 - Hu ); //excedance flow to the second storage [m] [m/min] check units
        //if ground is frozen, x1 goes directly to the surface
        //therefore nothing is diverted to static tank
        if(frozen_ground == 1){
            x2 = x1;
        }
            
		double d1 = x1 - x2; // the input to static tank [m/min]
		double out1 = min(e_pot, h1); //evaporation from the static tank. it cannot evaporate more than h1 [m]
		//double out1 = (e_pot > h1) ? e_pot : 0.0;
		ans[STATE_STATIC] = d1 - out1; //differential equation of static storage


		//surface storage tank
		double infiltration = global_params[4]*c_1; //infiltration rate [m/min]
         if(frozen_ground == 1){
            infiltration = 0;
        }
		double x3 = min(x2, infiltration); //water that infiltrates to gravitational storage [m/min]
		double d2 = x2 - x3; // the input to surface storage [m] check units
        double alfa2 =global_params[6]; //velocity in m/s
        double w = alfa2 * L / A_h  * 60; // [1/min]
        w = min(1,w); //water can take less than 1 min (dt) to leave surface
        double out2 =0;
        out2  = h2 * w; //direct runoff [m/min]
		ans[STATE_SURFACE] = d2 - out2; //differential equation of surface storage


		// SUBSURFACE storage
		double percolation = global_params[5]*c_1; // percolation rate to aquifer [m/min]
		double x4 = min(x3,percolation); //water that percolates to aquifer storage [m/min]
		double d3 = x3 - x4; // input to gravitational storage [m/min]
		double alfa3 = global_params[7]* 24*60; //residence time [days] to [min].
        double out3=0;
        if(alfa3>=1)
		    out3 = h3/alfa3; //interflow [m/min]
		ans[STATE_SUBSURF] = d3 - out3; //differential equation for gravitational storage

		//aquifer storage
		double x5 = 0;//water loss to deeper aquifer [m]
		double d4 = x4 - x5;
		double alfa4 = global_params[8]* 24*60; //residence time [days] to [min].
        double out4=0;
        if(alfa4>=1)
		    out4 = h4/alfa4 ; //base flow [m/min]
		ans[STATE_GW] = d4 - out4; //differential equation for aquifer storage

		//channel storage

		double lambda_1 = global_params[1];
	    double invtau = params[3];// 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
	   	double c_2 = params[5];// = A_h / 60.0;	//  c_2

        //printf("state: %d\n", state);
        if(state == 0){ //no dams
            ans[STATE_DISCHARGE] = -q + (out2 + out3 + out4) * c_2; //[m/min] to [m3/s]
	        for (i = 0; i < num_parents; i++)
	            ans[STATE_DISCHARGE] += y_p[i * dim + STATE_DISCHARGE];
	        ans[STATE_DISCHARGE] = invtau * pow(q, lambda_1) * ans[STATE_DISCHARGE];    // discharge[0]
            dam_storage=1e-1;
            ans[STATE_DAM_STORAGE] = -dam_storage;
        }

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // DAM MODEL
        //state is the array index corresponding to the current storage in the qvs 
        if(state == 1 ){// dams
            int debug = 0;
            
            if(debug) printf("state: %d\n", state);

            double dam_input = 0;//m3
            double dam_output=0;//m3
            double dam_outflow=0; //m3s-1
            double rating_flow=0; //m3s-1
            double output_mode=2; //defines reservoir operation scheme. 0 = passive. 1 = conservation pools. 2=incremental outflow increases with storage            
            

            if(debug) printf("time: %f\n", t);
            //inflow from upstream links
            for (i = 0; i < num_parents; i++)
	            dam_input += y_p[i * dim + STATE_DISCHARGE] * 60.0; //m3/s to m3
            dam_input += ((out2 + out3 + out4) * c_2) * 60; //adds hillslope runoff to the channel input
            if(debug) printf("inflow from upstream in m3: %f\n", dam_input);
            // storage inputs
            if(debug) printf("initial dam storage in m3: %f\n", dam_storage);
            //dam_storage= dam_storage + dam_input;
            if(debug) printf("step1 dam storage in m3: %f\n", dam_storage + dam_input);

            //if(debug) printf("state: %u\n", state);
            unsigned int max_storage_pond_index = (int)qvs->n_values - 1;
            double max_storage_pond = qvs->points[max_storage_pond_index][0];
            // find the upper and lower storage threshold for this link
            double low_storage_threshold = max_storage_pond * factor_low_threshold;
            double high_storage_threshold = max_storage_pond * factor_high_threshold;
            //if(dam_storage + dam_input >0) debug=1; //only print results if there is water in the pond
            if(debug) printf("max.sto.pond %f\n",max_storage_pond);
            if(debug) printf("low_storage_factor and high_storage_factor: %f and %f\n", factor_low_threshold,factor_high_threshold);

            if(debug) printf("low_storage_threshold and high_storage_threshold: %f and %f\n", low_storage_threshold,high_storage_threshold);
            if (dam_storage + dam_input >= max_storage_pond)
            {
                if(debug) printf("storage >= max.sto.pond\n");
                //S_max = qvs->points[qvs->n_values - 1][0];
                double q_max = qvs->points[max_storage_pond_index][1];
                if(debug) printf("qmax : %f\n", q_max);
                dam_outflow = q_max;//m3s-1
            }
            else if (dam_storage + dam_input < max_storage_pond)
            {
                if(debug) printf("storage within rating range\n");
                double S = (dam_storage + dam_input < 0.0) ? 0.0 : dam_storage + dam_input;
                if(debug) printf("storage  : %f\n", S);
                //find discharge value in rating
                unsigned int ii;
                for (ii = 0; ii < qvs->n_values - 1; ii++){
				   // if(debug)printf("model403 storage qvs points: %f\n", qvs->points[ii][0]);
                    if (qvs->points[ii][0] <= S
						&& S < qvs->points[ii + 1][0])    
					break;
                }
                double q2 = qvs->points[ii + 1][1];
                double q1 = qvs->points[ii][1];
                double S2 = qvs->points[ii + 1][0];
                double S1 = qvs->points[ii][0];
                //dam_outflow = (q2 - q1) / (S2 - S1) * (S - S1) + q1;//m3s-1
                //rating_flow is Qrating in the flowchart
                rating_flow = (q2 - q1) / (S2 - S1) * (S - S1) + q1;//m3s-1

                if(output_mode == 0){ //passive
                    dam_outflow = rating_flow;
                }

                if(output_mode==1){
                    if(dam_storage + dam_input <= low_storage_threshold){
                        if(debug) printf("storage <= low_storage_threshold\n");
                    //is thhe pond storage too low, is it an emergency and we
                    //need to hurry to build storage
                        if(dam_storage + dam_input <= low_storage_threshold*0.8){
                            if(debug) printf("storage <= low_storage_threshold*.8\n");
                            dam_outflow = 1e-6;
                        }
                        if(dam_storage + dam_input > low_storage_threshold*0.8){
                            if(debug) printf("storage > low_storage_threshold*0.8\n");

                            //find ideal flow
                            double deltaS = low_storage_threshold-(dam_storage + dam_input);
                            double deltaQ = deltaS / 60.0;
                            double inflow = dam_input / 60.0;
                            double Qideal = inflow - deltaQ;
                            if(Qideal <=rating_flow){
                                dam_outflow = Qideal;
                            }
                            if(Qideal >rating_flow){
                                dam_outflow = rating_flow;
                            }   
                        }
                    } 
                    if(dam_storage + dam_input > low_storage_threshold){
                    //if the pond storage is greater than threslow but less than threshigh
                        if(dam_storage + dam_input < high_storage_threshold){
                            if(debug) printf("storage < high_storage_threshold\n");
                            double deltaS = (dam_storage + dam_input) - low_storage_threshold;
                            if(debug) printf("delS  : %f\n", deltaS);
                            double deltaQ = deltaS / 60.0;
                            if(debug) printf("delQ  : %f\n", deltaQ);
                            double inflow = dam_input / 60.0;
                            if(debug) printf("inflow  : %f\n", inflow);
                            double Qideal = factor_operating_window * (inflow + deltaQ);
                            if(debug) printf("Qideal  : %f\n", Qideal);
                            if(Qideal <=rating_flow){
                                dam_outflow = Qideal;
                            }
                            if(Qideal >rating_flow){
                                dam_outflow = rating_flow;
                            }
                        }
                        if(dam_storage + dam_input >= high_storage_threshold){
                            if(debug) printf("storage >= high_storage_threshold\n");
                            double deltaS = (dam_storage + dam_input) - low_storage_threshold;
                            double deltaQ = deltaS / 60.0;
                            double inflow = dam_input / 60.0;
                            double Qideal = inflow + deltaQ;
                            if(Qideal <=rating_flow){
                                dam_outflow = Qideal;
                            }
                            if(Qideal >rating_flow){
                                dam_outflow = rating_flow;
                            }
                        }
                    }
                }
                if(output_mode==2){
                    if(dam_storage + dam_input <= (max_storage_pond * 0.05)){
                        dam_outflow = rating_flow * 0.05;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.05) && (dam_storage + dam_input <= (max_storage_pond * 0.1))){
                        dam_outflow = rating_flow * 0.05;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.1) && (dam_storage + dam_input <= (max_storage_pond * 0.15))){
                        dam_outflow = rating_flow * 0.05;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.15) && (dam_storage + dam_input <= (max_storage_pond * 0.2))){
                        dam_outflow = rating_flow * 0.05;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.2) && (dam_storage + dam_input <= (max_storage_pond * 0.25))){
                        dam_outflow = rating_flow * 0.05;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.25) && (dam_storage + dam_input <= (max_storage_pond * 0.3))){
                        dam_outflow = rating_flow * 0.05;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.3) && (dam_storage + dam_input <= (max_storage_pond * 0.35))){
                        dam_outflow = rating_flow * 0.05;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.35) && (dam_storage + dam_input <= (max_storage_pond * 0.4))){
                        dam_outflow = rating_flow * 0.05;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.4) && (dam_storage + dam_input <= (max_storage_pond * 0.45))){
                        dam_outflow = rating_flow * 0.05;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.45) && (dam_storage + dam_input <= (max_storage_pond * 0.5))){
                        dam_outflow = rating_flow * 0.05;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.5) && (dam_storage + dam_input <= (max_storage_pond * 0.55))){
                        dam_outflow = rating_flow * 0.1;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.55) && (dam_storage + dam_input <= (max_storage_pond * 0.6))){
                        dam_outflow = rating_flow * 0.3;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.6) && (dam_storage + dam_input <= (max_storage_pond * 0.65))){
                        dam_outflow = rating_flow * 0.3;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.65) && (dam_storage + dam_input <= (max_storage_pond * 0.7))){
                        dam_outflow = rating_flow * 0.5;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.7) && (dam_storage + dam_input <= (max_storage_pond * 0.75))){
                        dam_outflow = rating_flow * 0.5;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.75) && (dam_storage + dam_input <= (max_storage_pond * 0.8))){
                        dam_outflow = rating_flow * 0.7;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.8) && (dam_storage + dam_input <= (max_storage_pond * 0.85))){
                        dam_outflow = rating_flow * 0.7;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.85) && (dam_storage + dam_input <= (max_storage_pond * 0.9))){
                        dam_outflow = rating_flow * 0.9;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.9) && (dam_storage + dam_input <= (max_storage_pond * 0.95))){
                        dam_outflow = rating_flow * 0.9;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.95)){
                        dam_outflow = rating_flow;
                    }                
                }
            
            //calculate dam outflow
            //dam_outflow = get_discharge_from_storage(dam_storage + dam_input); //m3s-1
            if(debug) printf("dam outflow in m3/s: %f\n", dam_outflow);

            ans[STATE_DISCHARGE] = dam_outflow - q;
            
            //storage output
            //dam_output =get_storage_from_discharge(dam_outflow);
            dam_output = dam_outflow * 60.0; //m3s-1 to m3
            if(debug) printf("dam output storage in m3: %f\n", dam_output);
            if(debug) printf("step2 dam storage in m3: %f\n", dam_storage + dam_input - dam_output);

            //dam_storage = max(dam_storage,0);
            //if(debug) printf("step3 dam storage in m3: %f\n", dam_storage);
            //ans[STATE_DAM_STORAGE] = dam_storage;
            ans[STATE_DAM_STORAGE] = dam_input - dam_output;

            if(debug) printf("\n");

            
        }

        // if (forcing_values[0]>1 && ratio<1) {
        //     printf("time: %f\n", t);
        //     printf(" rain in mm/hour: %f\n", forcing_values[0]);
        //     printf(" area hill, area basin, area ratio: %f %f %f\n", A_h,A_i,ratio);
        //     MPI_Abort(MPI_COMM_WORLD, 1);
        // }
    }
}



//Type 261
//Contains 2 states in the channel: dischage, storage
//Contains 3 layers on hillslope: ponded, top layer, soil
//Contains 3 extra states: 
//Order of parameters: A_i,L_i,A_h,S_h,T_L,h_b,k_D,k_dry,k_i | invtau,c_1,c_2,c_3
//The numbering is:	0   1   2   3   4   5   6   7     8  |   9     10  11  12
//Order of global_params: v_0,lambda_1,lambda_2,N,phi,v_B
//The numbering is:        0      1        2    3  4   5 
void TopLayerNonlinearExpSoilvel(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;
    double pers_to_permin = 60.0;

    //Global params
    double lambda_1 = global_params[1];
    double N = global_params[3];	//[]
    double phi = global_params[4];	//[]
    double v_B = global_params[5];	//[m/s]

                                        //Local params
    double L_i = params[1];	//[m]
    double A_h = params[2];	//[m^2]
    double S_h = params[3];	//[]
    double T_L = params[4];	//[m]
    double h_b = params[5];	//[m]
    double k_D = params[6];	//[1/s]
    double k_dry = params[7];	//[1/s]
    double k_i = params[8];	//[1/s]

                                //Precalculations
    double invtau = params[9];	//[1/min]
    double c_1 = params[10];
    double c_2 = params[11];
    double c_3 = params[12];

    //Forcings
    double data = forcing_values[0] * c_1;			//[mm/hr] -> [m/min]
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]
    double eta = forcing_values[2];					//[]

                                                    //System states
    double q = y_i[0];		//[m^3/s]
    double S = y_i[1];		//[m^3]
    double s_p = y_i[2];	//[m]
    double s_t = y_i[3];	//[m]
    double s_s = y_i[4];	//[m]
                            //double s_precip = y_i[5];	//[m]
                            //double V_r = y_i[6];	//[m^3]
    double q_b = y_i[7];	//[m^3/s]

                            //Evaporation
    double e_p, e_t, e_s;
    double Corr = s_p + s_t / T_L + s_s / (h_b - T_L);
    if (e_pot > 0.0 && Corr > 1e-12)
    {
        e_p = s_p * 1e3 * e_pot / Corr;
        e_t = s_t / T_L * e_pot / Corr;
        e_s = s_s / (h_b - T_L) * e_pot / Corr;
    }
    else
    {
        e_p = 0.0;
        e_t = 0.0;
        e_s = 0.0;
    }

    //Fluxes
    double k_2 = c_3 / eta;
    double q_pl = k_2 * pers_to_permin * pow(s_p, phi);	//[m/min]
    double q_pt = k_dry * pers_to_permin * ((1.0 - s_t / T_L > 0.0) ? pow(1.0 - s_t / T_L, N) : 0.0) * s_p;	//[m/min]
    double q_ts = k_i * pers_to_permin * s_t;	//[m/min]
    double q_sl = k_D * pers_to_permin * s_s;	//[m/min]

                                                //Discharge
    dam_TopLayerNonlinearExpSoilvel(y_i, dim, global_params, params, qvs, state, user, ans);	//ans is used for convenience
    double qm = ans[0] * 60.0;

    //Storage
    ans[1] = (q_pl + q_sl) * A_h - qm;
    for (i = 0; i<num_parents; i++)
        ans[1] += y_p[i * dim] * 60.0;

    //Hillslope
    ans[2] = data - q_pl - q_pt - e_p;
    ans[3] = q_pt - q_ts - e_t;
    ans[4] = q_ts - q_sl - e_s;

    //Additional states
    ans[5] = data;
    ans[6] = q_pl;
    ans[7] = q_sl * A_h - q_b*60.0;
    for (i = 0; i<num_parents; i++)
        ans[7] += y_p[i * dim + 7] * 60.0;
    ans[7] *= v_B / L_i;
}

//Type 261
//Contains 2 states in the channel: dischage, storage
//Contains 3 layers on hillslope: ponded, top layer, soil
//Order of parameters: A_i,L_i,A_h,S_h,T_L,h_b,k_D,k_dry,k_i | invtau,c_1,c_2,c_3
//The numbering is:	0   1   2   3   4   5   6   7     8  |   9     10  11  12
//Order of global_params: v_0,lambda_1,lambda_2,N,phi,v_B
//The numbering is:        0      1        2    3  4   5 
void dam_TopLayerNonlinearExpSoilvel(const double * const y_i, unsigned int num_dof, const double * const global_params, const double * const params, const QVSData * const qvs, int state, void* user, double *ans)
{
    double q1, q2, S1, S2, S_max, q_max, S;

    //Parameters
    double lambda_1 = global_params[1];
    double invtau = params[9];

    //Find the discharge in [m^3/s]
    if (state == -1)
    {
        S = (y_i[1] < 0.0) ? 0.0 : y_i[1];
        //ans[0] = invtau/60.0*pow(S,1.0/(1.0-lambda_1));
        ans[0] = pow((1.0 - lambda_1)*invtau / 60.0 * S, 1.0 / (1.0 - lambda_1));
    }
    else if (state == (int)qvs->n_values - 1)
    {
        S_max = qvs->points[qvs->n_values - 1][0];
        q_max = qvs->points[qvs->n_values - 1][1];
        ans[0] = q_max;
    }
    else
    {
        S = (y_i[1] < 0.0) ? 0.0 : y_i[1];
        q2 = qvs->points[state + 1][1];
        q1 = qvs->points[state][1];
        S2 = qvs->points[state + 1][0];
        S1 = qvs->points[state][0];
        ans[0] = (q2 - q1) / (S2 - S1) * (S - S1) + q1;
    }
}

//Type 261
//Contains 2 states in the channel: dischage, storage
//Contains 3 layers on hillslope: ponded, top layer, soil
//Order of parameters: A_i,L_i,A_h,S_h,T_L,h_b,k_D,k_dry,k_i | invtau,c_1,c_2,c_3
//The numbering is:	0   1   2   3   4   5   6   7     8  |   9     10  11  12
//Order of global_params: v_0,lambda_1,lambda_2,N,phi,v_B
//The numbering is:        0      1        2    3  4   5 
void TopLayerNonlinearExpSoilvel_Reservoirs(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    ans[0] = forcing_values[3];
    ans[1] = 0.0;
    ans[2] = 0.0;
    ans[3] = 0.0;
    ans[4] = 0.0;
    ans[5] = 0.0;
    ans[6] = 0.0;
    ans[7] = 0.0;
}

//Type 262
//Contains 2 states in the channel: dischage, storage
//Contains 3 layers on hillslope: ponded, top layer, soil
//Order of parameters: A_i,L_i,A_h,S_h,T_L,eta,h_b,k_D,k_dry,k_i | invtau,c_1,c_2,k_2
//The numbering is:	0   1   2   3   4   5   6   7   8     9  |   10    11  12  13
//Order of global_params: v_0,lambda_1,lambda_2,N,phi,v_B
//The numbering is:        0      1        2    3  4   5
void TopLayerNonlinearExpSoilvel_ConstEta(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;
    double pers_to_permin = 60.0;

    //Global params
    double lambda_1 = global_params[1];
    double N = global_params[3];	//[]
    double phi = global_params[4];	//[]
    double v_B = global_params[5];	//[m/s]

                                        //Local params
    double L_i = params[1];	//[m]
    double A_h = params[2];	//[m^2]
    double S_h = params[3];	//[]
    double T_L = params[4];	//[m]
    double eta = params[5];	//[]
    double h_b = params[6];	//[m]
    double k_D = params[7];	//[1/s]
    double k_dry = params[8];	//[1/s]
    double k_i = params[9];	//[1/s]

                                //Precalculations
    double invtau = params[10];	//[1/min]
    double c_1 = params[11];
    double c_2 = params[12];
    double k_2 = params[13];

    //Forcings
    double data = forcing_values[0] * c_1;			//[mm/hr] -> [m/min]
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]

                                                                    //System states
    double q = y_i[0];		//[m^3/s]
    double S = y_i[1];		//[m^3]
    double s_p = y_i[2];	//[m]
    double s_t = y_i[3];	//[m]
    double s_s = y_i[4];	//[m]
                            //double s_precip = y_i[5];	//[m]
                            //double V_r = y_i[6];	//[m^3]
    double q_b = y_i[7];	//[m^3/s]

                            //Evaporation
    double e_p, e_t, e_s;
    double Corr = s_p + s_t / T_L + s_s / (h_b - T_L);
    if (e_pot > 0.0 && Corr > 1e-12)
    {
        e_p = s_p * 1e3 * e_pot / Corr;
        e_t = s_t / T_L * e_pot / Corr;
        e_s = s_s / (h_b - T_L) * e_pot / Corr;
    }
    else
    {
        e_p = 0.0;
        e_t = 0.0;
        e_s = 0.0;
    }

    //Fluxes
    //double k_2 = c_3 / eta;
    double q_pl = k_2 * pers_to_permin * pow(s_p, phi);	//[m/min]
    double q_pt = k_dry * pers_to_permin * ((1.0 - s_t / T_L > 0.0) ? pow(1.0 - s_t / T_L, N) : 0.0) * s_p;	//[m/min]
    double q_ts = k_i * pers_to_permin * s_t;	//[m/min]
    double q_sl = k_D * pers_to_permin * s_s;	//[m/min]

                                                //Discharge
    dam_TopLayerNonlinearExpSoilvel_ConstEta(y_i, dim, global_params, params, qvs, state, user, ans);	//ans is used for convenience
    double qm = ans[0] * 60.0;

    //Storage
    ans[1] = (q_pl + q_sl) * A_h - qm;
    for (i = 0; i<num_parents; i++)
        ans[1] += y_p[i * dim] * 60.0;

    //Hillslope
    ans[2] = data - q_pl - q_pt - e_p;
    ans[3] = q_pt - q_ts - e_t;
    ans[4] = q_ts - q_sl - e_s;

    //Additional states
    ans[5] = data;
    ans[6] = q_pl;
    ans[7] = q_sl * A_h - q_b*60.0;
    for (i = 0; i<num_parents; i++)
        ans[7] += y_p[i * dim + 7] * 60.0;
    ans[7] *= v_B / L_i;
}


//Type 262
//Contains 2 states in the channel: dischage, storage
//Contains 3 layers on hillslope: ponded, top layer, soil
//Order of parameters: A_i,L_i,A_h,S_h,T_L,eta,h_b,k_D,k_dry,k_i | invtau,c_1,c_2,k_2
//The numbering is:	0   1   2   3   4   5   6   7   8     9  |   10    11  12  13
//Order of global_params: v_0,lambda_1,lambda_2,N,phi,v_B
//The numbering is:        0      1        2    3  4   5
void dam_TopLayerNonlinearExpSoilvel_ConstEta(const double * const y_i, unsigned int num_dof, const double * const global_params, const double * const params, const QVSData * const qvs, int state, void* user, double *ans)
{
    double q1, q2, S1, S2, S_max, q_max, S;

    //Parameters
    double lambda_1 = global_params[1];
    double invtau = params[10];

    //Find the discharge in [m^3/s]
    if (state == -1)
    {
        S = (y_i[1] < 0.0) ? 0.0 : y_i[1];
        //ans[0] = invtau/60.0*pow(S,1.0/(1.0-lambda_1));
        ans[0] = pow((1.0 - lambda_1)*invtau / 60.0 * S, 1.0 / (1.0 - lambda_1));
    }
    else if (state == (int)qvs->n_values - 1)
    {
        S_max = qvs->points[qvs->n_values - 1][0];
        q_max = qvs->points[qvs->n_values - 1][1];
        ans[0] = q_max;
    }
    else
    {
        S = (y_i[1] < 0.0) ? 0.0 : y_i[1];
        q2 = qvs->points[state + 1][1];
        q1 = qvs->points[state][1];
        S2 = qvs->points[state + 1][0];
        S1 = qvs->points[state][0];
        ans[0] = (q2 - q1) / (S2 - S1) * (S - S1) + q1;
    }
}

//Type 262
//Contains 2 states in the channel: dischage, storage
//Contains 3 layers on hillslope: ponded, top layer, soil
//Order of parameters: A_i,L_i,A_h,S_h,T_L,eta,h_b,k_D,k_dry,k_i | invtau,c_1,c_2,k_2
//The numbering is:	0   1   2   3   4   5   6   7   8     9  |   10    11  12  13
//Order of global_params: v_0,lambda_1,lambda_2,N,phi,v_B
//The numbering is:        0      1        2    3  4   5
void TopLayerNonlinearExpSoilvel_ConstEta_Reservoirs(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    ans[0] = forcing_values[2];
    ans[1] = 0.0;
    ans[2] = 0.0;
    ans[3] = 0.0;
    ans[4] = 0.0;
    ans[5] = 0.0;
    ans[6] = 0.0;
    ans[7] = 0.0;
}

//Type 19
//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,c_1,c_2
//The numbering is:	0   1   2   3  4    5    6   7
//Order of global_params: v_r,lambda_1,lambda_2,RC,v_h,v_g,e_pot
//The numbering is:        0      1        2     3  4   5    6
void LinearHillslope_Evap(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;

    double lambda_1 = global_params[1];
    double e_pot = global_params[6] * (1e-3 / 60.0);	//[mm/hr]->[m/min]

    double A_h = params[2];
    double k2 = params[3];
    double k3 = params[4];
    double invtau = params[5];
    double c_1 = params[6];
    double c_2 = params[7];

    double q = y_i[0];		//[m^3/s]
    double s_p = y_i[1];	//[m]
    double s_a = y_i[2];	//[m]

    double q_pl = k2 * s_p;
    double q_al = k3 * s_a;

    //Evaporation
    double C_p, C_a, C_T, Corr_evap;
    if (e_pot > 0.0)
    {
        C_p = s_p / e_pot;
        C_a = s_a / e_pot;
        C_T = C_p + C_a;
    }
    else
    {
        C_p = 0.0;
        C_a = 0.0;
        C_T = 0.0;
    }

    //Corr_evap = (!state && C_T > 0.0) ? 1.0/C_T : 1.0;
    Corr_evap = (C_T > 1.0) ? 1.0 / C_T : 1.0;

    double e_p = Corr_evap * C_p * e_pot;
    double e_a = Corr_evap * C_a * e_pot;

    //Discharge
    ans[0] = -q + (q_pl + q_al) * A_h / 60.0;
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];

    //Hillslope
    ans[1] = forcing_values[0] * c_1 - q_pl - e_p;
    ans[2] = forcing_values[0] * c_2 - q_al - e_a;
}


//Type 19
//State = 0: C_T >= 1.0
//State = 1: C_T < 1.0
int LinearHillslope_Evap_Check(const double * const y_i, const double * const params, const double * const global_params, const QVSData * const qvs, unsigned int dam)
{
    double s_p = y_i[1];	//[m]
    double s_a = y_i[2];	//[m]

    double e_pot = global_params[6] * (1e-3 / 60.0);	//[mm/hr]->[m/min]

    double C_p, C_a, C_T;
    if (e_pot > 0.0)
    {
        C_p = s_p / e_pot;
        C_a = s_a / e_pot;
        C_T = C_p + C_a;
    }
    else
        return 0;

    if (C_T >= 1.0)	return 0;
    else 		return 1;
}



//Type 190
//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,c_1,c_2
//The numbering is:	0   1   2   3  4    5    6   7
//Order of global_params: v_r,lambda_1,lambda_2,RC,v_h,v_g
//The numbering is:        0      1        2     3  4   5 
void LinearHillslope_MonthlyEvap(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
        unsigned short i;

    double lambda_1 = global_params[1];
    //double e_pot = global_params[6] * (1e-3/60.0);	//[mm/hr]->[m/min]

    double A_h = params[2];
    double k2 = params[3];
    double k3 = params[4];
    double invtau = params[5];
    double c_1 = params[6];
    double c_2 = params[7];

    double q = y_i[0];		//[m^3/s]
    double s_p = y_i[1];	//[m]
    double s_a = y_i[2];	//[m]

    double q_pl = k2 * s_p;
    double q_al = k3 * s_a;

    //Evaporation
    double C_p, C_a, C_T, Corr_evap;
    //double e_pot = forcing_values[1] * (1e-3/60.0);
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]
	
	double infiltration_eff = forcing_values[2] + 1;

    if (e_pot > 0.0)
    {
        C_p = s_p / e_pot;
        C_a = s_a / e_pot;
        C_T = C_p + C_a;
    }
    else
    {
        C_p = 0.0;
        C_a = 0.0;
        C_T = 0.0;
    }

    //Corr_evap = (!state && C_T > 0.0) ? 1.0/C_T : 1.0;
    Corr_evap = (C_T > 1.0) ? 1.0 / C_T : 1.0;

    double e_p = Corr_evap * C_p * e_pot;
    double e_a = Corr_evap * C_a * e_pot;

    //Discharge
    ans[0] = -q + (q_pl + q_al) * A_h / 60.0;
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];

    //Hillslope
    ans[1] = forcing_values[0] * c_1 - q_pl - e_p;
    ans[2] = forcing_values[0] * c_2 - q_al - e_a;

    /*
    if ((s_p > 0.000010) && (s_p < 0.000012) && (t < 200)) {
        printf("time: %f\n", t);
        printf(" forc: %f %f\n", forcing_values[0], forcing_values[1]);
        printf(" flux: %f %f\n", ans[1], ans[2]);
        printf(" state: %f %f %f\n", q, s_p, s_a);
    }
    */
}


//Type 191
//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,c_1,c_2
//The numbering is:	0   1   2   3  4    5    6   7
//Order of global_params: v_r,lambda_1,lambda_2,RC,v_h,v_g,v_B
//The numbering is:        0      1        2     3  4   5   6
void LinearHillslope_MonthlyEvap_extras(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    double lambda_1 = global_params[1];
    double v_B = global_params[6];

    double L = params[1];
    double A_h = params[2];
    double k2 = params[3];
    double k3 = params[4];
    double invtau = params[5];
    double c_1 = params[6];
    double c_2 = params[7];

    double q = y_i[0];		//[m^3/s]
    double s_p = y_i[1];	//[m]
    double s_a = y_i[2];	//[m]
    double q_b = y_i[5];	//[m^3/s]

    double q_pl = k2 * s_p;
    double q_al = k3 * s_a;

    //Evaporation
    double C_p, C_a, C_T, Corr_evap;
    //double e_pot = forcing_values[1] * (1e-3/60.0);
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]

    if (e_pot > 0.0)
    {
        C_p = s_p / e_pot;
        C_a = s_a / e_pot;
        C_T = C_p + C_a;
    }
    else
    {
        C_p = 0.0;
        C_a = 0.0;
        C_T = 0.0;
    }

    //Corr_evap = (!state && C_T > 0.0) ? 1.0/C_T : 1.0;
    Corr_evap = (C_T > 1.0) ? 1.0 / C_T : 1.0;

    double e_p = Corr_evap * C_p * e_pot;
    double e_a = Corr_evap * C_a * e_pot;

    //Discharge
    ans[0] = -q + (q_pl + q_al) * A_h / 60.0;
    for (unsigned short i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];

    //Hillslope
    ans[1] = forcing_values[0] * c_1 - q_pl - e_p;
    ans[2] = forcing_values[0] * c_2 - q_al - e_a;

    //Additional states
    ans[3] = forcing_values[0] * c_1;
    ans[4] = q_pl;
    ans[5] = q_al * A_h - q_b * 60.0;
    for (unsigned short i = 0; i<num_parents; i++)
        ans[5] += y_p[i * dim + 5] * 60.0;
    //ans[5] += k3*y_p[i].ve[2]*A_h;
    ans[5] *= v_B / L;
}


//Type 191 and 192
//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,c_1,c_2
//The numbering is:	0   1   2   3  4    5    6   7
//Order of global_params: v_r,lambda_1,lambda_2,RC,v_h,v_g,v_B
//The numbering is:        0      1        2     3  4   5   6
void LinearHillslope_Reservoirs_extras(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    ans[0] = forcing_values[2];
    ans[1] = 0.0;
    ans[2] = 0.0;
    ans[3] = 0.0;
    ans[4] = 0.0;
    ans[5] = 0.0;
    ans[6] = 0.0;
}

//Type 192
//Order of parameters: A_i,L_i,A_h,k_2,k_i,invtau,c_1
//The numbering is:     0   1   2   3   4    5     6
//Order of global_params: v_r,lambda_1,lambda_2,k_I_factor,v_h,k_3,v_B
//The numbering is:        0      1        2        3       4   5   6
void LinearHillslope_MonthlyEvap_kI_extras(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    double lambda_1 = global_params[1];
    double v_B = global_params[6];

    double L = params[1];
    double A_h = params[2];
    double k_2 = params[3];
    double k_i = params[4];
    double k_3 = global_params[5];
    double invtau = params[5];

    double c_1 = params[6];

    double q = y_i[0];	    // [m^3/s]
    double s_p = y_i[1];    // [m]
    double s_a = y_i[2];    // [m]
    double q_b = y_i[5];    // [m^3/s]

    //Fluxes
    double q_pl = k_2 * s_p;
    double q_pa = k_i * s_p;
    double q_al = k_3 * s_a;

    //Evaporation
    double C_p, C_a, C_T, Corr_evap;
    //double e_pot = forcing_values[1] * (1e-3/60.0);
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]

    if (e_pot > 0.0)
    {
        C_p = s_p / e_pot;
        C_a = s_a / e_pot;
        C_T = C_p + C_a;
    }
    else
    {
        C_p = 0.0;
        C_a = 0.0;
        C_T = 0.0;
    }

    //Corr_evap = (!state && C_T > 0.0) ? 1.0/C_T : 1.0;
    Corr_evap = (C_T > 1.0) ? 1.0 / C_T : 1.0;

    double e_p = Corr_evap * C_p * e_pot;
    double e_a = Corr_evap * C_a * e_pot;

    //Discharge
    ans[0] = -q + (q_pl + q_al) * A_h / 60.0;
    for (unsigned short i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];

    //Hillslope
    ans[1] = (forcing_values[0] * c_1) - q_pa - q_pl - e_p;
    ans[2] = q_pa - q_al - e_a;

    //Additional states
    ans[3] = forcing_values[0] * c_1;
    ans[4] = q_pl;
    ans[5] = q_al * A_h - q_b * 60.0;
    for (unsigned short i = 0; i<num_parents; i++)
        ans[5] += y_p[i * dim + 5] * 60.0;
    //ans[5] += k3*y_p[i].ve[2]*A_h;
    ans[5] *= v_B / L;
}

//Type 193
//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,c_1,c_2
//The numbering is:	0   1   2   3  4    5    6   7
//Order of global_params: v_r,lambda_1,lambda_2,RC,v_h,v_g
//The numbering is:        0      1        2     3  4   5
void routing_runoff2(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;

    double lambda_1 = global_params[1];

	double A_i = params[0];
	double L_i = params[1];
	double A_h = params[2];
	double invtau = params[3];
	//printf("invtau: %f\n", invtau);
    double q = y_i[0];		                                        // [m^3/s]

    double runoff = forcing_values[0]*(0.001 / 60.0);               //(mm/hr ->m/min)


    //in this model, all the water is routed downstream
    ans[0] = -q + (runoff * A_h / 60.0); //m3/min to m3/s
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    
	
    //Hillslope
    //ans[1] = q_rp - q_pl;

    //Sub-surface
   // ans[2] = q_ra - q_al - e_a;

    //Accumulated precip
    //ans[3] = q_rp + q_ra;


}

//Type 194
//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,c_1,c_2
//The numbering is:	0   1   2   3  4    5    6   7
//Order of global_params: v_r,lambda_1,lambda_2,RC,v_h,v_g
//The numbering is:        0      1        2     3  4   5
void routing_runoff1(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;

    double lambda_1 = global_params[1];

	double A_i = params[0];
	double L_i = params[1];
	double A_h = params[2];
	double invtau = params[3];
	//printf("invtau: %f\n", invtau);
    double q = y_i[0];		                                        // [m^3/s]

    double runoff = forcing_values[0]*(0.001 / 60.0);               //(mm/hr ->m/min)



    ans[0] = -q + (runoff * A_h / 60.0); //m3/min to m3/s
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];
	
    //Hillslope
    //ans[1] = q_rp - q_pl;

    //Sub-surface
   // ans[2] = q_ra - q_al - e_a;

    //Accumulated precip
    //ans[3] = q_rp + q_ra;


}

//Type 195
//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,c_1,c_2
//The numbering is:	0   1   2   3  4    5    6   7
//Order of global_params: v_r,lambda_1,lambda_2,RC,v_h,v_g
//The numbering is:        0      1        2     3  4   5 
void LinearHillslope_MonthlyEvap_OnlyRouts(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;

    double lambda_1 = global_params[1];

    double A_h = params[2];
    double k2 = params[3];
    double k3 = params[4];
    double invtau = params[5];

    double q = y_i[0];		                                        // [m^3/s]
	double s_p = y_i[1];	                                        // [m]
    double s_a = y_i[2];	                                        // [m]
	double acc = y_i[3];	                                        // [m]

	double q_rp = forcing_values[0] * (0.001 / 60.0);		        // (mm/hr -> m/min)
	double q_pl = k2 * s_p;                                         // (1/min * m)
	
	double q_ra = forcing_values[1] * (0.001 / 60.0);               // (mm/hr -> m/min)
    double q_al = k3 * s_a;                                         // (1/min * m)

    //Evaporation
    double C_a, C_T, Corr_evap;
    double e_pot = forcing_values[2] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]

    if (e_pot > 0.0){
        C_a = s_a / e_pot;
        C_T = C_a;
    } else {
        C_a = 0.0;
        C_T = 0.0;
    }

    Corr_evap = (C_T > 1.0) ? 1.0 / C_T : 1.0;

    double e_a = Corr_evap * C_a * e_pot;
	double q_parent;
	int q_pidx;

    //Discharge
    ans[0] = -q + ((q_al + q_pl) * A_h / 60.0);
	for (i = 0; i < num_parents; i++) {
		q_pidx = i * dim;
		q_parent = y_p[q_pidx];
		ans[0] += q_parent;
	}
    ans[0] = invtau * pow(q, lambda_1) * ans[0];

    //Hillslope
	ans[1] = q_rp - q_pl;
	
	//Sub-surface
    ans[2] = q_ra - q_al - e_a;	

	//Accumulated precip
	ans[3] = q_rp + q_ra;

    /*
    if (t < 0.030) {
        printf("time: %f\n", t);
    } else if ((t > 0.030) && (t < 0.035)) {
        printf("time: %f\n", t);
        printf(" forc: %f %f %f\n", forcing_values[0], forcing_values[1], forcing_values[2]);
        printf(" semi: %f %f\n", q_rp, q_ra);
        printf(" flux: %f %f\n", ans[1], ans[2]);
        printf(" state: %f %f %f %f\n", q, s_p, s_a, acc);
    }
    */
}


//Type 196
//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,c_1,c_2
//The numbering is:	0   1   2   3  4    5    6   7
//Order of global_params: v_r,lambda_1,lambda_2,RC,v_h,v_g
//The numbering is:        0      1        2     3  4   5 
void LinearHillslope_MonthlyEvap_OnlyRouts_NotReservoir(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
	unsigned short i;

	double lambda_1 = global_params[1];

	double A_h = params[2];
	double k2 = params[3];
	double k3 = params[4];
	double invtau = params[5];

	double q = y_i[0];		                                        // [m^3/s]
	double s_p = y_i[1];	                                        // [m]
	double s_a = y_i[2];	                                        // [m]
	double acc = y_i[3];	                                        // [m]

	double q_rp = forcing_values[0] * (0.001 / 60.0);		        // (mm/hr -> m/min)
	double q_pl = k2 * s_p;                                         // (1/min * m)

	double q_ra = forcing_values[1] * (0.001 / 60.0);               // (mm/hr -> m/min)
	double q_al = k3 * s_a;                                         // (1/min * m)

																	//Evaporation
	double C_a, C_T, Corr_evap;
	double e_pot = forcing_values[2] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]

	if (e_pot > 0.0) {
		C_a = s_a / e_pot;
		C_T = C_a;
	}
	else {
		C_a = 0.0;
		C_T = 0.0;
	}

	Corr_evap = (C_T > 1.0) ? 1.0 / C_T : 1.0;

	double e_a = Corr_evap * C_a * e_pot;
	double q_parent;
	int q_pidx;

	//Discharge
	ans[0] = -q + ((q_al + q_pl) * A_h / 60.0);
	for (i = 0; i < num_parents; i++) {
		q_pidx = i * dim;
		q_parent = y_p[q_pidx];
		ans[0] += q_parent;
	}
	ans[0] = invtau * pow(q, lambda_1) * ans[0];

	//Hillslope
	ans[1] = q_rp - q_pl;

	//Sub-surface
	ans[2] = q_ra - q_al - e_a;

	//Accumulated precip
	ans[3] = q_rp + q_ra;

	//Accumulated runoff
	ans[4] = q_rp;
}


//Type 196
//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,c_1,c_2
//The numbering is:	0   1   2   3  4    5    6   7
//Order of global_params: v_r,lambda_1,lambda_2,RC,v_h,v_g,v_B
//The numbering is:        0      1        2     3  4   5   6
void LinearHillslope_MonthlyEvap_OnlyRouts_HasReservoir(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
	ans[0] = forcing_values[3];
	ans[1] = 0.0;
	ans[2] = 0.0;
	ans[3] = 0.0;
	ans[4] = 0.0;
}


//Type 20
//Order of parameters: A_i,L_i,A_h,invtau,c_1,c_2
//The numbering is:	0   1   2    3     4   5
//Order of global_params: v_r,lambda_1,lambda_2,beta,k_p,k_a,theta_p,theta_a,scale_p,scale_a
//The numbering is:        0      1        2     3    4   5     6       7      8        9

//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,c_1,c_2
//The numbering is:	0   1   2   3  4    5    6   7
//Order of global_params: v_r,lambda_1,lambda_2,RC,v_h,v_g
//The numbering is:        0      1        2     3  4   5 
void Hillslope_Toy(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{

    unsigned short i;
    double lambda_1 = global_params[1];


    //For complex
    double omega = global_params[3];
    double alpha = global_params[4];
    double beta = global_params[5];
    double A_p = global_params[6];
    double A_a = global_params[7];

    /*
    //For gamma
    double beta = global_params[3];
    double k_p = global_params[4];
    double k_a = global_params[5];
    double theta_p = global_params[6];
    double theta_a = global_params[7];
    double scale_p = global_params[8];
    double scale_a = global_params[9];
    */
    /*
    //For lognormal
    double mu = global_params[3];
    double sigma2 = global_params[4];
    */
    /*
    //For exponential
    double lambda = global_params[3];
    double scale = global_params[4];
    */
    /*
    //For line and drop
    double t_p = global_params[3];
    double p = global_params[4];
    double alpha = global_params[5];
    */
    /*
    //For GEV
    double mu = global_params[3];
    double sigma = global_params[4];
    double ksi = global_params[5];
    double shift = global_params[6];
    */
    /*
    //For linear
    double k2 = global_params[3];
    double k3 = global_params[4];
    double scale = global_params[5];
    */
    //General use
    double A_h = params[2];
    double invtau = params[3];
    double c_1 = params[4];
    double c_2 = params[5];

    //States
    double q = y_i[0];
    double s_p = y_i[1];
    double s_a = y_i[2];

    //For Complex Model
    double Q_pl = A_p * pow(sin(omega*s_p), 2.0);
    double Q_al = A_a * (exp(alpha*s_a) - 1.0);
    double Q_pa = beta * s_p;

    /*
    //For Gamma Forcing Model
    //double Q_pl = A_p * pow(sin(omega*s_p),2.0);
    //double Q_al = A_a * (exp(alpha*s_a) - 1.0);
    double Q_pl = scale_p * c_1 * pow(t,k_p-1.0) * exp(-t/theta_p); //[m/min]
    //double Q_al = 0.0;
    double Q_al = scale_a * c_2 * pow(t,k_a-1.0) * exp(-t/theta_a); //[m/min]
    double Q_pa = beta * s_p; //[m/min]
    */
    /*
    //For Log Normal Forcing Model
    double Q_pl = (t <= 0.0) ? 0.0 : c_1/t * exp( -sq(log(t)-mu)/(2.0*sigma2) );
    double Q_al = 0.0;
    */
    /*
    //For Exponential Forcing Model
    double Q_pl = scale * lambda * exp(-lambda*t);
    double Q_al = 0.0;
    */
    /*
    //For Line and Drop Model
    double Q_pl = (t <= t_p) ? p/t_p * t : p * exp(-alpha*(t-t_p));
    double Q_al = 0.0;
    */
    /*
    //For GEV
    double holder = pow(1.0 + (t-shift-mu)/sigma*ksi,-1.0/ksi);
    //double holder = (t < 25.0) ? pow(1.0 + (t-shift-mu)/sigma*ksi,-1.0/ksi) : 0.0;
    //double holder = (t < 6.0 + shift) ? pow(1.0 + (t-shift-mu)/sigma*ksi,-1.0/ksi) : 0.0;
    //double holder = (t < 1000.0 + shift) ? pow(1.0 + (t-shift-mu)/sigma*ksi,-1.0/ksi) : 0.0;
    double Q_pl = 1.0/sigma * pow(holder,ksi+1) * exp(-holder);
    double Q_al = 0.0;
    */
    /*
    //For Linear Model
    double Q_pl = scale * (k2 * 0.015 * exp(-k2*t) + k3 * 0.015 * exp(-k3*t));
    double Q_al = 0.0;
    */
    //Discharge
    ans[0] = -q + (Q_pl + Q_al) * A_h / 60.0;
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];

    //Hillslope
    //ans[1] = 0.0;
    //ans[2] = 0.0;
    ans[1] = forcing_values[0] * (0.001 / 60.0) - Q_pl - Q_pa;
    ans[2] = -Q_al + Q_pa;

}

//Type 60
//Order of parameters: A_i,L_i,A_h,h_b,k2,k3,invtau,c_1
//The numbering is:	0   1   2   3   4  5    6    7 
//Order of global_params: v_r,lambda_1,lambda_2,v_h,v_g,e_pot
//The numbering is:        0      1        2     3   4   5  
void LinearHillslope_Evap_RC(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;

    double lambda_1 = global_params[1];
    double e_pot = global_params[5] * (1e-3 / 60.0);	//[mm/hr]->[m/min]

    double A_h = params[2];	//[m^2]
    double h_b = params[3];	//[m]
    double k2 = params[4];
    double k3 = params[5];
    double invtau = params[6];
    double c_1 = params[7];

    double q = y_i[0];		//[m^3/s]
    double s_p = y_i[1];	//[m]
    double s_a = y_i[2];	//[m]

    double q_pl = k2 * s_p;
    double q_al = k3 * s_a;

    //Run off
    double d_soil = (s_a > h_b) ? 0.0 : h_b - s_a;
    double RC = (s_p > 0.0) ? (s_p*(s_p + 2.0*d_soil)) / sq(s_p + d_soil) : 0.0;

    //Evaporation
    double C_p, C_a, C_T, Corr_evap;
    if (e_pot > 0.0)
    {
        C_p = s_p / e_pot;
        C_a = s_a / e_pot;
        C_T = C_p + C_a;
    }
    else
    {
        C_p = 0.0;
        C_a = 0.0;
        C_T = 0.0;
    }
    Corr_evap = (C_T > 1.0) ? 1.0 / C_T : 1.0;

    double e_p = Corr_evap * C_p * e_pot;
    double e_a = Corr_evap * C_a * e_pot;

    //Discharge
    ans[0] = -q + (q_pl + q_al) * A_h / 60.0;
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];

    //Hillslope
    ans[1] = forcing_values[0] * c_1*RC - q_pl - e_p;
    ans[2] = forcing_values[0] * c_1*(1.0 - RC) - q_al - e_a;
}


// Type 21 ************************************************

//Type 21
//Order of parameters: A_i,L_i,A_h,k2,k3,invtau
//The numbering is:	0   1   2  3  4    5
//Order of global_params: v_r,lambda_1,lambda_2,RC,S_0,v_h,v_g
//The numbering is:        0      1        2     3  4   5   6
void nodam_rain_hillslope(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;
    double lambda_1 = global_params[1];
    double RC = global_params[3];
    double A_h = params[2];
    double k2 = params[3];
    double k3 = params[4];
    double invtau = params[5];
    //double F_et = 1.0;
    double F_et = 0.05;

    double S = y_i[1];
    double Ss = y_i[2];
    double Sg = y_i[3];

    //double q = invtau*S;
    double qm = invtau*pow(S, 1.0 / (1.0 - lambda_1));

    ans[1] = k2*Ss + k3*Sg - qm;
    for (i = 0; i<num_parents; i++)	ans[1] += y_p[i * dim] * 60.0;
    //ans[1] = ans[1];

    ans[2] = RC*forcing_values[0] * A_h*(.001 / 60.0) - k2*Ss;
    ans[3] = (1.0 - RC)*F_et*forcing_values[0] * A_h*(.001 / 60.0) - k3*Sg;
    //ans[3] = 0.0;
}


//Type 21
//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
//The numbering is:	0   1   2  3  4    5	       6      7       8     9	  10	    11       12  13  14
//Order of global_params: v_r,lambda_1,lambda_2,RC,S_0,v_h,v_g
//The numbering is:        0      1        2     3  4   5   6
void dam_rain_hillslope(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;
    double lambda_1 = global_params[1];
    double RC = global_params[3];
    double A_h = params[2];
    double k2 = params[3];
    double k3 = params[4];
    double invtau = params[5];
    double orifice_area = params[6];
    double H_spill = params[7];
    double H_max = params[8];
    double S_max = params[9];
    double alpha = params[10];
    double diam = params[11];
    double c_1 = params[12];
    double c_2 = params[13];
    double L_spill = params[14];
    double g = 9.81;
    double F_et = 1.0;

    double qm;
    double S = y_i[1];
    double Ss = y_i[2];
    double Sg = y_i[3];
    double h = H_max * pow(S / S_max, alpha);
    double diff = (h - H_spill >= 0) ? h - H_spill : 0.0;

    if (state == 0)
        //q = invtau*S;
        qm = invtau*pow(S, 1.0 / (1.0 - lambda_1));
    else if (state == 1)
        qm = c_1*orifice_area*pow(2 * g*h, .5) * 60.0;
    //qm = 0.0;
    else if (state == 2)
        qm = c_1*orifice_area*pow(2 * g*h, .5) * 60.0 + c_2*L_spill*pow(diff, 1.5) * 60.0;
    else if (state == 3)
        qm = c_1*orifice_area*pow(2 * g*h, .5) * 60.0 + c_2*L_spill*pow(diff, 1.5) * 60.0 + invtau*pow(S - S_max, 1.0 / (1.0 - lambda_1));
    else if (state == 4)
    {
        double r = diam / 2.0;
        double frac = (h < 2 * r) ? (h - r) / r : 1.0;
        double A = -r*r*(acos(frac) - pow(1 - frac*frac, .5)*frac - 3.1415926535897932);
        qm = c_1*A*pow(2 * g*h, .5) * 60.0;
    }
    else printf("Error in function evaluation. Invalid state %u.\n", state);

    ans[1] = k2*Ss + k3*Sg - qm;
    for (i = 0; i<num_parents; i++)	ans[1] += y_p[i * dim] * 60.0;
    //ans[1] = ans[1];

    ans[2] = RC*forcing_values[0] * A_h*(.001 / 60.0) - k2*Ss;
    ans[3] = (1.0 - RC)*F_et*forcing_values[0] * A_h*(.001 / 60.0) - k3*Sg;
}

//Type 21
//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
//The numbering is:	0   1   2  3  4    5	       6      7       8     9	  10	    11       12  13  14
//Order of global_params: v_r,lambda_1,lambda_2,RC,S_0,v_h,v_g
//The numbering is:        0      1        2     3  4   5   6
void dam_q(const double * const y_i, unsigned int  num_dof, const double * const global_params, const double * const params, const QVSData * const qvs, int state, void* user, double *ans)
{
    double lambda_1 = global_params[1];
    double invtau = params[5];
    double S = (y_i[1] < 0.0) ? 0.0 : y_i[1];

    if (state == 0)
        //ans[0] = invtau*S;
        ans[0] = invtau / 60.0*pow(S, 1.0 / (1.0 - lambda_1));
    else
    {
        double orifice_area = params[6];
        double H_spill = params[7];
        double H_max = params[8];
        double S_max = params[9];
        double alpha = params[10];
        double diam = params[11];
        double c_1 = params[12];
        double c_2 = params[13];
        double L_spill = params[14];
        double g = 9.81;

        double h = H_max * pow(S / S_max, alpha);
        double diff = (h - H_spill >= 0) ? h - H_spill : 0.0;

        if (state == 1)
            ans[0] = c_1*orifice_area*pow(2 * g*h, .5);
        //ans[0] = 0.0;
        else if (state == 2)
            ans[0] = c_1*orifice_area*pow(2 * g*h, .5) + c_2*L_spill*pow(diff, 1.5);
        else if (state == 3)
            ans[0] = c_1*orifice_area*pow(2 * g*h, .5) + c_2*L_spill*pow(diff, 1.5) + invtau / 60.0*pow(S - S_max, 1.0 / (1.0 - lambda_1));
        else if (state == 4)
        {
            double r = diam / 2.0;
            double frac = (h < 2 * r) ? (h - r) / r : 1.0;
            double A = -r*r*(acos(frac) - pow(1.0 - frac*frac, .5)*frac - 3.141592653589);
            ans[0] = c_1*A*pow(2 * g*h, .5);
        }
        else printf("Error in dam evaluation. Invalid state %u.\n", state);
    }
}


//Type 21
//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
//The numbering is:	0   1   2  3  4    5	       6      7       8     9	  10	    11       12  13  14
//Order of global_params: v_r,lambda_1,lambda_2,RC,S_0,v_h,v_g
//The numbering is:        0      1        2     3  4   5   6
int dam_check(
    double *y, unsigned int num_dof,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    QVSData *qvs,
    bool has_dam,
    void *user)
{
    if (!has_dam)
        return 0;

    double H_spill = params[7];
    double H_max = params[8];
    double S_max = params[9];
    double alpha = params[10];
    double diam = params[11];
    double S = y[1];
    double h = H_max * pow(S / S_max, alpha);

    if (h < diam)		return 4;
    if (h <= H_spill)	return 1;
    if (h <= H_max)		return 2;
    return 3;
}

// Type 22 ************************************************

//Type 22
//Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau
//The numbering is:	0   1   2   3  4   5   6  7    8
//Order of global_params: lambda_1,lambda_2,S_0,v_g
//The numbering is:         0        1       2   3
void nodam_rain_hillslope2(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;
    double lambda_1 = global_params[0];
    double A_h = params[2];
    double RC = params[3];
    double k2 = params[6];
    double k3 = params[7];
    double invtau = params[8];

    double S = y_i[1];
    double Ss = y_i[2];
    double Sg = y_i[3];

    //double q = invtau*S;
    double qm = invtau*pow(S, 1.0 / (1.0 - lambda_1));

    ans[1] = k2*Ss + k3*Sg - qm;
    for (i = 0; i<num_parents; i++)	ans[1] += y_p[i * dim] * 60.0;
    //ans[1] = ans[1];

    ans[2] = RC*forcing_values[0] * A_h*(.001 / 60.0) - k2*Ss;
    ans[3] = (1.0 - RC)*forcing_values[0] * A_h*(.001 / 60.0) - k3*Sg;
    //ans[3] = 0.0;
}


//Type 22
//Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
//The numbering is:	0   1   2  3   4   5   6  7   8          9	  10	  11   12     13      14        15  16   17
//Order of global_params: lambda_1,lambda_2,S_0,v_g
//The numbering is:         0        1       2   3
void dam_rain_hillslope2(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;
    double lambda_1 = global_params[0];
    double A_h = params[2];
    double RC = params[3];
    double k2 = params[6];
    double k3 = params[7];
    double invtau = params[8];
    double orifice_area = params[9];
    double H_spill = params[10];
    double H_max = params[11];
    double S_max = params[12];
    double alpha = params[13];
    double diam = params[14];
    double c_1 = params[15];
    double c_2 = params[16];
    double L_spill = params[17];
    double g = 9.81;

    double qm;
    double S = y_i[1];
    double Ss = y_i[2];
    double Sg = y_i[3];
    double h = H_max * pow(S / S_max, alpha);
    double diff = (h - H_spill >= 0) ? h - H_spill : 0.0;

    if (state == 0)
        //q = invtau*S;
        qm = invtau*pow(S, 1.0 / (1.0 - lambda_1));
    else if (state == 1)
        qm = c_1*orifice_area*pow(2 * g*h, .5) * 60.0;
    //qm = 0.0;
    else if (state == 2)
        qm = c_1*orifice_area*pow(2 * g*h, .5) * 60.0 + c_2*L_spill*pow(diff, 1.5) * 60.0;
    else if (state == 3)
        qm = c_1*orifice_area*pow(2 * g*h, .5) * 60.0 + c_2*L_spill*pow(diff, 1.5) * 60.0 + invtau*pow(S - S_max, 1.0 / (1.0 - lambda_1));
    else if (state == 4)
    {
        double r = diam / 2.0;
        double frac = (h < 2 * r) ? (h - r) / r : 1.0;
        double A = -r*r*(acos(frac) - pow(1 - frac*frac, .5)*frac - 3.1415926535897932);
        qm = c_1*A*pow(2 * g*h, .5) * 60.0;
    }
    else printf("Error in function evaluation. Invalid state %u.\n", state);

    ans[1] = k2*Ss + k3*Sg - qm;
    for (i = 0; i<num_parents; i++)	ans[1] += y_p[i * dim] * 60.0;
    //ans[1] = ans[1];

    ans[2] = RC*forcing_values[0] * A_h*(.001 / 60.0) - k2*Ss;
    ans[3] = (1.0 - RC)*forcing_values[0] * A_h*(.001 / 60.0) - k3*Sg;
}


//Type 22
//Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
//The numbering is:	0   1   2  3   4   5   6  7   8          9	  10	  11   12     13      14        15  16   17
//Order of global_params: lambda_1,lambda_2,S_0,v_g
//The numbering is:         0        1       2   3
void dam_q2(const double * const y_i, unsigned int num_dof, const double * const global_params, const double * const params, const QVSData * const qvs, int state, void* user, double *ans)
{
    double lambda_1 = global_params[0];
    double invtau = params[8];
    double S = (y_i[1] < 0.0) ? 0.0 : y_i[1];

    if (state == 0)
        //ans[0] = invtau*S;
        ans[0] = invtau / 60.0*pow(S, 1.0 / (1.0 - lambda_1));
    else
    {
        double orifice_area = params[9];
        double H_spill = params[10];
        double H_max = params[11];
        double S_max = params[12];
        double alpha = params[13];
        double diam = params[14];
        double c_1 = params[15];
        double c_2 = params[16];
        double L_spill = params[17];
        double g = 9.81;

        double h = H_max * pow(S / S_max, alpha);
        double diff = (h - H_spill >= 0) ? h - H_spill : 0.0;

        if (state == 1)
            ans[0] = c_1*orifice_area*pow(2 * g*h, .5);
        //ans[0] = 0.0;
        else if (state == 2)
            ans[0] = c_1*orifice_area*pow(2 * g*h, .5) + c_2*L_spill*pow(diff, 1.5);
        else if (state == 3)
            ans[0] = c_1*orifice_area*pow(2 * g*h, .5) + c_2*L_spill*pow(diff, 1.5) + invtau / 60.0*pow(S - S_max, 1.0 / (1.0 - lambda_1));
        else if (state == 4)
        {
            double r = diam / 2.0;
            double frac = (h < 2 * r) ? (h - r) / r : 1.0;
            double A = -r*r*(acos(frac) - pow(1.0 - frac*frac, .5)*frac - 3.141592653589);
            ans[0] = c_1*A*pow(2 * g*h, .5);
        }
        else printf("Error in dam evaluation. Invalid state %u.\n", state);
    }
}


//Type 22
//Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
//The numbering is:	0   1   2  3   4   5   6  7   8          9	  10	  11   12     13      14        15  16   17
//Order of global_params: lambda_1,lambda_2,S_0,v_g
//The numbering is:         0        1       2   3
int dam_check2(
    double *y, unsigned int num_dof,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    QVSData *qvs,
    bool has_dam,
    void *user)
{
    if (!has_dam)
        return 0;

    double H_spill = params[10];
    double H_max = params[11];
    double S_max = params[12];
    double alpha = params[13];
    double diam = params[14];
    double S = y[1];
    double h = H_max * pow(S / S_max, alpha);

    if (h < diam)		return 4;
    if (h <= H_spill)	return 1;
    if (h <= H_max)		return 2;
    return 3;
}

// Type 23 ************************************************

//Type 23
//Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau
//The numbering is:	0   1   2   3  4   5   6  7    8
//Order of global_params: lambda_1,lambda_2,S_0,v_g
//The numbering is:         0        1       2   3
void nodam_rain_hillslope3(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;
    double lambda_1 = global_params[0];
    double A_h = params[2];
    double RC = params[3];
    double k2 = params[6];
    double k3 = params[7];
    double invtau = params[8];

    double S = y_i[1];
    double Ss = y_i[2];
    double Sg = y_i[3];

    //double q = invtau*S;
    double qm = invtau*pow(S, 1.0 / (1.0 - lambda_1));

    ans[1] = k2*Ss + k3*Sg - qm;
    for (i = 0; i<num_parents; i++)	ans[1] += y_p[i * dim] * 60.0;
    //ans[1] = ans[1];

    ans[2] = RC*forcing_values[0] * A_h*(.001 / 60.0) - k2*Ss;
    ans[3] = (1.0 - RC)*forcing_values[0] * A_h*(.001 / 60.0) - k3*Sg;
    //ans[3] = 0.0;
}


//Type 23
//Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
//The numbering is:	0   1   2  3   4   5   6  7   8          9	  10	  11   12     13      14        15  16   17
//Order of global_params: lambda_1,lambda_2,S_0,v_g
//The numbering is:         0        1       2   3
void dam_rain_hillslope3(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;
    double lambda_1 = global_params[0];
    double A_h = params[2];
    double RC = params[3];
    double k2 = params[6];
    double k3 = params[7];
    double invtau = params[8];
    double orifice_area = params[9];
    double H_spill = params[10];
    double H_max = params[11];
    double S_max = params[12];
    double alpha = params[13];
    double diam = params[14];
    double c_1 = params[15];
    double c_2 = params[16];
    double L_spill = params[17];
    double g = 9.81;

    double qm;
    double S = y_i[1];
    double Ss = y_i[2];
    double Sg = y_i[3];
    double h = H_max * pow(S / S_max, alpha);
    double diff = (h - H_spill >= 0) ? h - H_spill : 0.0;

    if (state == 0)
        //q = invtau*S;
        qm = invtau*pow(S, 1.0 / (1.0 - lambda_1));
    else if (state == 1)
        qm = c_1*orifice_area*pow(2 * g*h, .5) * 60.0;
    //qm = 0.0;
    else if (state == 2)
        qm = c_1*orifice_area*pow(2 * g*h, .5) * 60.0 + c_2*L_spill*pow(diff, 1.5) * 60.0;
    else if (state == 3)
        qm = c_1*orifice_area*pow(2 * g*h, .5) * 60.0 + c_2*L_spill*pow(diff, 1.5) * 60.0 + invtau*pow(S - S_max, 1.0 / (1.0 - lambda_1));
    else if (state == 4)
    {
        //if(h < 1e-6)	qm = 0.0;
        qm = c_2 * 2 * sqrt(h*diam - h*h) * pow(h, 1.5) * 60.0;
    }
    else if (state == 5)
        qm = c_2 * diam * pow(h, 1.5) * 60.0;
    else printf("Error in function evaluation. Invalid state %u.\n", state);

    ans[1] = k2*Ss + k3*Sg - qm;
    for (i = 0; i<num_parents; i++)	ans[1] += y_p[i * dim] * 60.0;
    //ans[1] = ans[1];

    ans[2] = RC*forcing_values[0] * A_h*(.001 / 60.0) - k2*Ss;
    ans[3] = (1.0 - RC)*forcing_values[0] * A_h*(.001 / 60.0) - k3*Sg;
}


//Type 23
//Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
//The numbering is:	0   1   2  3   4   5   6  7   8          9	  10	  11   12     13      14        15  16   17
//Order of global_params: lambda_1,lambda_2,S_0,v_g
//The numbering is:         0        1       2   3
void dam_q3(const double * const y_i, unsigned int num_dof, const double * const global_params, const double * const params, const QVSData * const qvs, int state, void* user, double *ans)
{
    double lambda_1 = global_params[0];
    double invtau = params[8];
    double S = (y_i[1] < 0.0) ? 0.0 : y_i[1];

    if (state == 0)
        //ans[0] = invtau*S;
        ans[0] = invtau / 60.0*pow(S, 1.0 / (1.0 - lambda_1));
    else
    {
        double orifice_area = params[9];
        double H_spill = params[10];
        double H_max = params[11];
        double S_max = params[12];
        double alpha = params[13];
        double diam = params[14];
        double c_1 = params[15];
        double c_2 = params[16];
        double L_spill = params[17];
        double g = 9.81;

        double h = H_max * pow(S / S_max, alpha);
        double diff = (h - H_spill >= 0) ? h - H_spill : 0.0;

        if (state == 1)
            ans[0] = c_1*orifice_area*pow(2 * g*h, .5);
        //ans[0] = 0.0;
        else if (state == 2)
            ans[0] = c_1*orifice_area*pow(2 * g*h, .5) + c_2*L_spill*pow(diff, 1.5);
        else if (state == 3)
            ans[0] = c_1*orifice_area*pow(2 * g*h, .5) + c_2*L_spill*pow(diff, 1.5) + invtau / 60.0*pow(S - S_max, 1.0 / (1.0 - lambda_1));
        else if (state == 4)
        {
            //if(h < 1e-6)	ans[0] = 0.0;
            ans[0] = c_2 * 2 * sqrt(h*diam - h*h) * pow(h, 1.5);
        }
        else if (state == 5)
            ans[0] = c_2 * diam * pow(h, 1.5);
        else printf("Error in dam evaluation. Invalid state %u.\n", state);
    }
}


//Type 23
//Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
//The numbering is:	0   1   2  3   4   5   6  7   8          9	  10	  11   12     13      14        15  16   17
//Order of global_params: lambda_1,lambda_2,S_0,v_g
//The numbering is:         0        1       2   3
int dam_check3(
    double *y, unsigned int num_dof,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    QVSData *qvs,
    bool has_dam,
    void *user)
{
    if (!has_dam)
        return 0;

    double H_spill = params[10];
    double H_max = params[11];
    double S_max = params[12];
    double alpha = params[13];
    double diam = params[14];
    double S = y[1];
    double h = H_max * pow(S / S_max, alpha);

    if (h < 0.5*diam)	return 4;
    if (h < diam)		return 5;
    if (h <= H_spill)	return 1;
    if (h <= H_max)		return 2;
    return 3;
}


// Type 40 ************************************************

//Type 40
//Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau
//The numbering is:	0   1   2   3  4   5   6  7    8
//Order of global_params: lambda_1,lambda_2,S_0,v_g
//The numbering is:         0        1       2   3
void nodam_rain_hillslope_qsv(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;
    double lambda_1 = global_params[0];
    double A_h = params[2];
    double RC = params[3];
    double k2 = params[6];
    double k3 = params[7];
    double invtau = params[8];

    double S = y_i[1];
    double Ss = y_i[2];
    double Sg = y_i[3];

    double qm = invtau*pow(S, 1.0 / (1.0 - lambda_1));

    ans[1] = k2*Ss + k3*Sg - qm;
    for (i = 0; i<num_parents; i++)	ans[1] += y_p[i * dim] * 60.0;

    ans[2] = RC*forcing_values[0] * A_h*(.001 / 60.0) - k2*Ss;
    ans[3] = (1.0 - RC)*forcing_values[0] * A_h*(.001 / 60.0) - k3*Sg;
}


//Type 40
//Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau
//The numbering is:	0   1   2  3   4   5   6  7   8
//Order of global_params: lambda_1,lambda_2,S_0,v_g
//The numbering is:         0        1       2   3
void dam_rain_hillslope_qsv(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;
    double lambda_1 = global_params[0];
    double A_h = params[2];
    double RC = params[3];
    double k2 = params[6];
    double k3 = params[7];
    double invtau = params[8];

    double qm, q1, q2, S1, S2, S_max, q_max;
    double S = y_i[1];
    double Ss = y_i[2];
    double Sg = y_i[3];

    if (state == -1)
        qm = invtau*pow(S, 1.0 / (1.0 - lambda_1));
    else if (state == (int)qvs->n_values - 1)
    {
        S_max = qvs->points[qvs->n_values - 1][0];
        q_max = qvs->points[qvs->n_values - 1][1];
        qm = q_max*60.0 + invtau*pow(max(S - S_max, 0.0), 1.0 / (1.0 - lambda_1));
        //qm = q_max * 60.0;
    }
    else
    {
        q2 = qvs->points[state + 1][1];
        q1 = qvs->points[state][1];
        S2 = qvs->points[state + 1][0];
        S1 = qvs->points[state][0];
        qm = ((q2 - q1) / (S2 - S1) * (S - S1) + q1) * 60.0;
    }

    ans[1] = k2*Ss + k3*Sg - qm;
    for (i = 0; i<num_parents; i++)	ans[1] += y_p[i * dim] * 60.0;

    ans[2] = RC*forcing_values[0] * A_h*(.001 / 60.0) - k2*Ss;
    ans[3] = (1.0 - RC)*forcing_values[0] * A_h*(.001 / 60.0) - k3*Sg;
}


//Type 40
//Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau
//The numbering is:	0   1   2  3   4   5   6  7   8
//Order of global_params: lambda_1,lambda_2,S_0,v_g
//The numbering is:         0        1       2   3
void dam_q_qvs(const double * const y_i, unsigned int num_dof, const double * const global_params, const double * const params, const QVSData * const qvs, int state, void* user, double *ans)
{
    double S_max, q_max;
    double lambda_1 = global_params[0];
    double invtau = params[8];
    double S = (y_i[1] < 0.0) ? 0.0 : y_i[1];
    double q1, q2, S1, S2;

    if (state == -1)
        ans[0] = invtau / 60.0*pow(S, 1.0 / (1.0 - lambda_1));
    else if (state == ((int)(qvs->n_values) - 1))
    {
        S_max = qvs->points[qvs->n_values - 1][0];
        q_max = qvs->points[qvs->n_values - 1][1];
        ans[0] = q_max + invtau / 60.0*pow(max(S - S_max, 0.0), 1.0 / (1.0 - lambda_1));
        //ans[0] = q_max;
    }
    else
    {
        q2 = qvs->points[state + 1][1];
        q1 = qvs->points[state][1];
        S2 = qvs->points[state + 1][0];
        S1 = qvs->points[state][0];
        ans[0] = (q2 - q1) / (S2 - S1) * (S - S1) + q1;
    }
}


// Rodica's Models **************************

//Type 0
//Function for simple river system
//Calculates the flow using simple parameters, using only the flow q.
//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
//The numbering is:        0      1        2     3   4   5
//This uses the units and functions from September 18, 2011 document
//y_i[0] = q
//void simple_river(double t,double *y_i,double *y_p,unsigned short num_parents, unsigned int max_dim,double *global_params,double *forcing_values,double *params,int state,void* user,double *ans)
void simple_river(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;

    double q = y_i[0];
    double lambda_1 = global_params[1];
    double invtau = params[12];

    //Flux equation (y_i[0])
    ans[0] = -q;
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];
}

//Type 0
//This assumes there is only 1 variable being passed link to link.
//double NormJx_simple_river(double *y_i,double *global_params,double *params)
void Jsimple_river(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, double *ans)
{
    ans[0] = params[12] * pow(y_i[0], global_params[1]);
}

//Type 0
//Jacobian of the simple river system******************
//y_i has q_i (y_i[0])
//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11   12     13     14  15  16  17  18  19
//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
//The numbering is:        0      1        2     3   4   5
//This uses the units and functions from September 18, 2011 document
void Jsimple(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, double *ans)
{
    unsigned short i;
    double q = y_i[0];
    double lambda_1 = global_params[1];
    double invtau = params[12];

    ans[0] = 0.0;
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = -(lambda_1 + 1.0) * invtau * pow(q, lambda_1) + lambda_1 * invtau * pow(q, lambda_1 - 1.0) * ans[0];
}

//Type 1/3
//Function for the simple river system with rainfall
//Calculates the flow using simple parameters, using only the flow q.
//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
//The numbering is:        0      1        2     3   4   5
//This uses the units and functions from September 18, 2011 document
//y_i[0] = q, y_i[1] = s
void river_rainfall(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;

    double q = y_i[0];
    double s = y_i[1];

    double invtau = params[12];
    double c_1 = params[14];
    double c_3 = params[16];
    double c_4 = params[17];
    double lambda_1 = global_params[1];

    double s53 = pow(s, 5.0 / 3.0);

    //Flux equation (y_i[0])
    ans[0] = -q + c_1 * s53;
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];

    //Ponded water equation (y_i[1])
    ans[1] = c_3 * forcing_values[0] - c_4 * s53;
}

//Type 2
//Function for simple river system with hillslope runoff
//Calculates the flow using simple parameters, using only the flow q.
//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
//The numbering is:        0      1        2     3   4   5
//This uses the units and functions from September 18, 2011 document
//y_i[0] = q, y_i[1] = s
void simple_hillslope(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;

    double q = y_i[0];
    double s = y_i[1];

    double invtau = params[12];
    double c_1 = params[14];
    double c_4 = params[17];
    double lambda_1 = global_params[1];

    double s53 = pow(s, 5.0 / 3.0);

    //Flux equation (y_i[0])
    ans[0] = -q + c_1 * s53;
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];

    //Ponded water equation (y_i[1])
    ans[1] = -c_4 * s53;
}

//Type 4
//Function for simple river system with soil properties******************
//Calculates the flow using the soil parameters.
//y_i has q_i (y_i[0]), s (y_i[1]), a (y_i[2]), v (y_i[3])
//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
//The numbering is:        0      1        2     3   4   5
//This uses the units and functions from September 18, 2011 document  
void simple_soil(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;

    double q = y_i[0];
    double s = y_i[1];
    double a = y_i[2];
    double v = y_i[3];

    double invtau = params[12];
    double epsilon = params[13];
    double c_1 = params[14];
    double c_4 = params[17];
    double c_5 = params[18];
    double c_6 = params[19];
    double lambda_1 = global_params[1];

    double s53 = pow(s, 5.0 / 3.0);
    double exp_term, big_frac;

    if (a < 1.0)
    {
        exp_term = c_6 * (1.0 - a) * v * exp(-(1.0 - a - v) / (1.0 - a));
        //big_frac = c_5 * s * epsilon*epsilon * (1.0-a) * (1.0-a-v)*(1.0-a-v)/((s*(1.0-a) + epsilon*(1.0 - a - v))*(s*(1.0-a) + 1.0 - a - v));
        big_frac = c_5 * s * sq(epsilon) * (1.0 - a) * sq(1.0 - a - v) / sq(s*(1.0 - a) + epsilon*(1.0 - a - v));
    }
    else
    {
        exp_term = 0.0;
        big_frac = 0.0;
    }

    //Flux equation (y_i[0])
    ans[0] = -q + c_1 * s53;
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1)*ans[0];

    //Ponded water equation (y_i[1])
    ans[1] = -c_4 * s53 - big_frac + c_6 * a * a;

    //Impermeable area (y_i[2])
    //ans[2] = -c_6 * a * a + c_6 * (1.0 - a) + exp_term;
    ans[2] = (-c_6 * a * a + exp_term) / epsilon;

    //Soil Volumetric Water Content (y_i[3])
    //ans[3] = big_frac - exp_term;
    ans[3] = (big_frac - exp_term) / epsilon;
}

//Type 5
//Function for simple river system with soil properties and rainfall******************
//Calculates the flow using the soil parameters.
//y_i has q_i (y_i[0]), s (y_i[1]), a (y_i[2]), v (y_i[3])
//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11   12     13     14  15  16  17  18  19
//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
//The numbering is:        0      1        2     3   4   5
//This uses the units and functions from September 18, 2011 document
void soil_rainfall(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;

    double q = y_i[0];
    double s = y_i[1];
    double a = y_i[2];
    double v = y_i[3];

    double invtau = params[12];
    double epsilon = params[13];
    double c_1 = params[14];
    double c_3 = params[16];
    double c_4 = params[17];
    double c_5 = params[18];
    double c_6 = params[19];
    double lambda_1 = global_params[1];

    double s53 = pow(s, 5.0 / 3.0);
    double exp_term, big_frac;

    if (a < 1.0)
    {
        exp_term = c_6 * (1.0 - a) * v * exp(-(1.0 - a - v) / (1.0 - a));
        //big_frac = c_5 * s * epsilon*epsilon * (1.0-a) * (1.0-a-v)*(1.0-a-v)/((s*(1.0-a) + epsilon*(1.0 - a - v))*(s*(1.0-a) + 1.0 - a - v));
        big_frac = c_5 * s * sq(epsilon) * (1.0 - a) * sq(1.0 - a - v) / sq(s*(1.0 - a) + epsilon*(1.0 - a - v));
    }
    else
    {
        exp_term = 0.0;
        big_frac = 0.0;
    }

    //Flux equation (y_i[0])
    ans[0] = -q + c_1 * s53;
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1)*ans[0];

    //Ponded water equation (y_i[1])
    ans[1] = c_3 * forcing_values[0] - c_4 * s53 - big_frac + c_6 * a * a;

    //Impermeable area (y_i[2])
    //ans[2] = -c_6 * a * a + c_6 * (1.0 - a) + exp_term;
    ans[2] = (-c_6 * a * a + exp_term) / epsilon;

    //Soil Volumetric Water Content (y_i[3])
    //ans[3] = big_frac - exp_term;
    ans[3] = (big_frac - exp_term) / epsilon;
}

//Type 6
//RHS function for discharge, surface ponding, saturated zone, unsaturated zone, and vegetation
//Includes rainfall
//This uses units and functions from Rodica's 1st (published) paper
//Order of parameters: L_i,A_h,A_i,h_b,K_sat,K_sp,d_0,d_1,d_2,invbeta,alpha_N,alpha_soil,a_res,v_res,invtau,gamma,c_1,c_2,c_3,c_4
//The numbering is:     0   1   2   3   4     5    6   7   8     9      10       11       12    13     14    15   16  17  18  19
//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r
//The numbering is:        0      1        2     3   4
void qsav_rainfall(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;

    //State variables
    double q = y_i[0];
    double s_p = y_i[1];
    double a = y_i[2];
    double v = y_i[3];

    //Parameters and constants
    double h_b = params[3];
    double d_0 = params[6];
    double d_1 = params[7];
    double d_2 = params[8];
    double invbeta = params[9];
    double a_res = params[12];
    double v_res = params[13];
    double invtau = params[14];
    double gamma = params[15];
    double c_1 = params[16];
    double c_2 = params[17];
    double c_3 = params[18];
    double c_4 = params[19];
    double lambda_1 = global_params[1];

    //Calculations
    double a_adj = a - a_res;
    double v_adj = v - v_res;
    double av_adjsum = a_adj + v_adj;
    double polyterm = d_0 * v_adj + d_1 * v_adj * sq(a_adj) + d_2 * sq(a_adj);
    double expterm = exp(a_adj / c_4);

    //Evaporation
    /*
    double c_evap = 1.0/60.0;
    if(t < 8000.0)		c_evap *= 0.002;
    else if(t < 15000.0)	c_evap *= 0.0005;
    else if(t < 30000.0)	c_evap *= 0.0001;
    else			c_evap = 0.0;

    double evap;
    if(t < 5000.0)		evap = 0.62;
    else if(t < 18000.0)	evap = 0.48;
    else if(t < 25000.0)	evap = 0.5;
    else if(t < 30000.0)	evap = 0.48;
    else if(t < 35000.0)	evap = 0.51;
    else if(t < 40000.0)	evap = 0.44;
    else			evap = 0.3;
    */
    double c_evap = 0.0;
    double evap = 0.0;

    //Flux equation
    ans[0] = -q + gamma * c_1 * s_p * av_adjsum + gamma * c_2 * a_adj * expterm;
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];

    //Ponded water equation
    //ans[1] = (1.0 - evap) * c_3 * forcing_values[0] - c_1 * s_p * av_adjsum - c_1 * s_p * (h_b - av_adjsum);
    ans[1] = (1.0 - evap) * c_3 * forcing_values[0] - c_1 * h_b * s_p;

    //Saturated zone equation (a)
    ans[2] = invbeta * (polyterm - c_2 * a_adj * expterm - c_evap * a_adj);

    //Unsaturated zone equation (v)
    ans[3] = invbeta * (c_1 * s_p * (h_b - av_adjsum) - polyterm);
}



//Type 15
//Function for the simple river system with rainfall
//Calculates the flow using simple parameters, using only the flow q.
//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
//Order of global_params: v_r,lambda_1,lambda_2,v_h,A_r,RC
//The numbering is:        0      1        2     3   4   5
//This uses the units and functions from September 18, 2011 document
//y_i[0] = q, y_i[1] = s_p
void river_rainfall_adjusted(
    double t,
    const double * const y_i, unsigned int dim,
    const double * const y_p, unsigned short num_parents, unsigned int max_dim,
    const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;

    double q = y_i[0];
    double s_p = y_i[1];

    double L = params[0];
    double invtau = params[12];
    double c_1 = params[14];
    double c_3 = params[16];
    double c_4 = params[17];
    double lambda_1 = global_params[1];

    //double q_pl = pow(0.5*L*s_p/(L+s_p),2.0/3.0);
    double q_pl = s_p;

    //Flux equation (y_i[0])
    ans[0] = -q + c_1 * q_pl;
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];

    //Ponded water equation (y_i[1])
    ans[1] = c_3 * forcing_values[0] - c_4 * q_pl;
    //ans[1] = c_3 * (forcing_values[0] + 100.0*sin(t/5.0)) - c_4 * q_pl;
}


//Type 315
//Function for simple river system with data assimilation.
//Calculates the flow using simple parameters, using only the flow q.
//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
//Order of global_params: v_r,lambda_1,lambda_2,v_h,A_r,RC
//The numbering is:        0      1        2     3   4   5
//This uses the units and functions from September 18, 2011 document
//y_i[0] = q, y_i[1] = s, followed by N entries for the variational equation
void assim_river_rainfall_adjusted(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    printf("Error: assim models do not exist...\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
    /*
    unsigned int i,j;
    unsigned int dim = ans.dim;
    unsigned int offset = 2;		//!!!! This needs to be num_dense, but without variational eqs !!!!
    unsigned int parent_offset;
    unsigned int problem_dim = 2;
    unsigned int all_states = (dim-offset)/problem_dim;
    unsigned int loc = iparams.ve[0];
    double inflow = 0.0;

    double q = y_i[0];
    double s_p = y_i[1];

    double L = params[0];
    double invtau = params[12];
    double c_1 = params[14];
    double c_3 = params[16];
    double c_4 = params[17];
    double lambda_1 = global_params[1];

    double q_to_lambda_1 = pow(q,lambda_1);
    //double q_to_lambda_1_m1 = (q > 1e-12) ? q_to_lambda_1 / q : 1e-12;
    double q_to_lambda_1_m1 = (q > 1e-12) ? q_to_lambda_1 / q : pow(1e-12,lambda_1 - 1.0);
    double deriv_qpl = 1.0;

    double q_pl = s_p;

    //Flux equation (y_i[0])
    ans[0] = -q + c_1 * q_pl;
    for(i=0;i<num_parents;i++)
    inflow += y_p[i * dim];
    ans[0] = invtau * q_to_lambda_1 * (inflow + ans[0]);

    //Ponded water equation (y_i[1])
    ans[1] = c_3 * forcing_values[0] - c_4 * q_pl;
    //ans[1] = c_3 * ( max(forcing_values[0] + 20.0*sin(t/5.0),0.0)) - c_4 * q_pl;

    //!!!! Pull if statements out of loops (should just need two cases total) !!!!
    //!!!! A lot of terms get repeated !!!!

    //Eqs for variational equations
    for(i=offset;i<dim;i++)	ans[i] = 0.0;

    //s variable from this link
    ans[offset] = -c_4*deriv_qpl*y_i[offset];

    //q variables from this link
    //	if(lambda_1 > 1e-12 && (inflow) > 1e-12)
    ans[offset + 1] = (lambda_1 * invtau * q_to_lambda_1_m1 * (inflow + c_1*s_p) - (lambda_1 + 1) * invtau * q_to_lambda_1) * y_i[offset + 1];
    //	else
    //		ans[offset + 1] = -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i[offset + 1];

    //	if(lambda_1 > 1e-12 && (inflow) > 1e-12)
    ans[offset + 2] = (lambda_1 * invtau * q_to_lambda_1_m1 * (inflow + c_1*s_p) - (lambda_1 + 1) * invtau * q_to_lambda_1) * y_i[offset + 2] + invtau*c_1*q_to_lambda_1*deriv_qpl * y_i[offset];
    //	else
    //		ans[offset + 2] = -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i[offset + 2] + invtau*c_1*deriv_qpl*y_i[offset];

    //Adjust offset
    offset += 3;

    //Variables from parents
    for(i=0;i<num_parents;i++)
    {
    parent_offset = 1 + problem_dim;

    for(j=0;j<numupstream[i];j++)
    {
    ans[offset] = invtau * q_to_lambda_1 * y_p[i].ve[parent_offset];
    //			if(lambda_1 > 1e-12 && (inflow) > 1e-12)
    ans[offset] += (lambda_1 * invtau * q_to_lambda_1_m1 * (inflow + c_1*s_p) - (lambda_1 + 1) * invtau * q_to_lambda_1) * y_i[offset];
    //			else
    //				ans[offset] += -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i[offset];

    ans[offset + 1] = invtau * q_to_lambda_1 * y_p[i].ve[parent_offset + 1];
    //			if(lambda_1 > 1e-12 && (inflow) > 1e-12)
    ans[offset + 1] += (lambda_1 * invtau * q_to_lambda_1_m1 * (inflow + c_1*s_p) - (lambda_1 + 1) * invtau * q_to_lambda_1) * y_i[offset + 1];
    //			else
    //				ans[offset + 1] += -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i[offset + 1];

    offset += 2;
    parent_offset += 2;
    }
    }
    */
}



//Type 30
//Function for river system with soil properties and rainfall******************
//This uses Luciana's model
//Calculates the flow using the soil parameters.
//y_i has q (y_i[0]), s_p (y_i[1]), h_w (y_i[2]), theta (y_i[3])
//Order of parameters: L_h,A_h,A_up,H_b,H_h,max_inf_rate,K_SAT,S_H,n_vh, b, c, d,K_Q,V_T,c_1,c_2,c_3,c_4,c_5,c_6,c_7
//The numbering is:     0   1   2    3   4       5         6    7   8    9 10 11  12  13  14  15  16  17  18  19  20
//Order of global_params: v_0,lambda_1,lambda_2,Q_r,A_r,K_T,C_r,e_pot
//The numbering is:        0      1        2     3   4   5   6   7
/*
void lcuencas_soilrain(double t,double *y_i,double *y_p,unsigned short num_parents, unsigned int max_dim,double *global_params,double *forcing_values,QVSData *qvs,double *params,int state,void* user,double *ans)
{
unsigned int i;

//States
double q = y_i[0];		//[m^3/s]
double s_p = y_i[1];	//[mm]
double h_w = y_i[2];	//[m]
double theta = y_i[3];	//[%]

//Global params
//double v_0 = global_params[0];
double lambda_1 = global_params[1];
//double lambda_2 = global_params[2];
//double A_r = global_params[4];
double K_T = global_params[5];
double e_pot = global_params[7];

//Local params
double A_h = params[1];
double H_b = params[3];
double H_h = params[4];
double maxinf = params[5];
//double K_SAT = params[6];
double K_SAT = (H_h > 1e-12) ? params[6] : 0.0;
double b = params[9];
double c = params[10];
double d = params[11];
double K_Q = params[12];
double V_T = params[13];
double c_1 = params[14];
double c_2 = params[15];
double c_3 = params[16];
double c_4 = params[17];
double c_5 = params[18];
double c_6 = params[19];
double c_7 = params[20];

//Other
//double H_T = H_b + H_h;
double h_rel = h_w - H_b * 1e-3;
//double H_relmax = H_T - H_b;
double H_relmax = (H_h > 1e-12) ? H_h : 1e-6;
double hoverH = h_rel / H_relmax;
double a_IoverA_H = b * hoverH + c * sq(hoverH) + d * sq(hoverH) * hoverH;
if(a_IoverA_H < 0.0)	a_IoverA_H = 1e-6;
if(a_IoverA_H > 0.9999)	a_IoverA_H = 0.9999;
double a_PoverA_H = 1.0 - a_IoverA_H;
double v_H = c_1 * pow(s_p * 1e-3,2.0/3.0);
if(v_H > 350.0)	v_H = 350.0;
double v_ssat = V_T * a_IoverA_H;
double v_sunsat = V_T - v_ssat;
double d_soil = (1.0 - theta) * H_b;
double RC = (s_p <= 0.0) ? 0.0 : s_p*(s_p + 2*d_soil) / sq(s_p + d_soil);
double deriv_a_I = c_7 * (b + 2*c*hoverH + 3*d*sq(hoverH));
if(deriv_a_I <= (1e-5*1e6)) deriv_a_I = (100.0*1e6);		//Should use this
//if(deriv_a_I <= (1e-0)) deriv_a_I = 1e-0;
double K_UNSAT = K_SAT * exp(5.0 * (theta - 1.0));
double infrate = (s_p > maxinf) ? infrate = maxinf : s_p * K_T;
//double infrate = s_p*K_T;
//if(infrate > maxinf)	infrate = maxinf;


//Evaporation
//double C_p = (e_pot > 0.0) ? s_p / e_pot : 0.0;
double D_unsat = 1e-3 * v_sunsat * theta / A_h;
double D_sat = 1e-3 * v_ssat / A_h;

double C_p,C_unsat,C_sat,C_T;
if(e_pot > 0.0)
{
C_p = s_p / e_pot;
C_unsat = D_unsat / e_pot;
C_sat = D_sat / e_pot;
C_T = C_p + C_unsat + C_sat;
}
else
{
C_p = 0.0;
C_unsat = 0.0;
C_sat = 0.0;
C_T = 0.0;
}


//double C_unsat = D_unsat / H_b * a_PoverA_H;
//double C_sat = D_sat/H_b * a_IoverA_H;
//double C_s = 0.0;
//C_T = C_s + C_unsat + C_sat;

double Corr_evap = (C_T > 1.0) ? 1.0/C_T : 1.0;
double e_p = Corr_evap * C_p * e_pot;
double e_sat = Corr_evap * C_sat * e_pot;
double e_unsat = Corr_evap * C_unsat * e_pot;


//Fluxes
double q_pl = c_2 * (a_PoverA_H * RC + a_IoverA_H) * v_H * s_p;
double q_sl = c_3 * v_ssat * a_IoverA_H;
double snow_M = 0.8 * 0.0;	//Assume no snow for now
double q_pu = a_PoverA_H * infrate * (1.0 - RC);
double q_us = 1e3 * K_UNSAT * theta * a_PoverA_H;

//!!!! This might be cleaner if H_b is converted into m !!!!

//Channel transport
ans[0] = -q;
for(i=0;i<num_parents;i++)	ans[0] += y_p[i * dim];
ans[0] = K_Q * pow(q,lambda_1) * (ans[0] + c_4 * (q_pl + q_sl + snow_M));

//Land surface storage
ans[1] = (1.0/60.0) * (forcing_values[0] - q_pl - q_pu - e_p);

//Water table
ans[2] = (H_b < 1e-6) ? 0.0 : c_5 / deriv_a_I * (q_us - q_sl - e_sat);

//Soil volumetric water content
//ans[3] = (v_sunsat < 1e-6) ? 0.0 : c_6 / v_sunsat * (q_pu - q_us - e_unsat + theta*(q_us - q_sl - e_sat));
//ans[3] = (v_sunsat < 1e-12) ? -100.0 : c_6 / v_sunsat * (q_pu - q_us - e_unsat + theta*(q_us - q_sl - e_sat));
ans[3] = c_6 / v_sunsat * (q_pu - q_us - e_unsat + theta*(q_us - q_sl - e_sat));
}
*/

void lcuencas_soilrain(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned int i;

    //States
    double q = y_i[0];
    double s_p = y_i[1];
    double h_w = y_i[2];
    double theta = y_i[3];

    //Global params
    //double v_0 = global_params[0];
    double lambda_1 = global_params[1];
    //double lambda_2 = global_params[2];
    //double A_r = global_params[4];
    double K_T = global_params[5];
    double e_pot = forcing_values[1];
    //double e_pot = global_params[7];

    //Local params
    //double L_h = params[0];
    double A_h = params[1];
    //double A_up = params[2];
    double H_b = params[3];
    double H_h = params[4];
    double maxinf = params[5];
    double K_SAT = params[6];
    //double S_H = params[7];
    //double n_vh = params[8];
    double b = params[9];
    double c = params[10];
    double d = params[11];
    double K_Q = params[12];
    double V_T = params[13];
    double c_1 = params[14];
    double c_2 = params[15];
    double c_3 = params[16];
    double c_4 = params[17];
    double c_5 = params[18];
    double c_6 = params[19];
    double c_7 = params[20];

    //Other
    //double H_T = H_b + H_h;
    double h_rel = h_w - H_b * 1e-3;
    //if(h_rel <= 0)	h_rel = 0.01;
    double H_relmax = (H_h > 1e-12) ? H_h : 1e-6;
    double hoverH = h_rel / H_relmax;
    double a_IoverA_H = b * hoverH + c * sq(hoverH) + d * sq(hoverH) * hoverH;
    if (a_IoverA_H < 0.0)	a_IoverA_H = 1e-6;
    //if(a_IoverA_H <= 0)	a_IoverA_H = 1e-8 / A_h;
    double a_PoverA_H = 1.0 - a_IoverA_H;
    double v_H = c_1 * pow(s_p * 1e-3, 2.0 / 3.0);
    if (v_H > 350.0)	v_H = 350.0;
    double v_ssat = V_T * a_IoverA_H;
    double v_sunsat = V_T - v_ssat;
    double d_soil = (1.0 - theta) * H_b;
    double RC = (s_p <= 0.0) ? 0.0 : s_p*(s_p + 2 * d_soil) / sq(s_p + d_soil);
    double deriv_a_I = c_7 * (b + 2 * c*hoverH + 3 * d*sq(hoverH));
    //if(deriv_a_I <= 0.0)	deriv_a_I = 100.0;
    if (deriv_a_I <= 1e-5) deriv_a_I = 100.0;
    double K_UNSAT = K_SAT * exp(5.0 * (theta - 1.0));
    double infrate = (s_p > maxinf) ? infrate = maxinf : s_p * K_T;

    //Evaporation
    //double e_pot = 0.0;	//!!!! This should be an input function !!!!
    double C_p = (e_pot > 0.0) ? s_p / e_pot : 0.0;
    double D_unsat = 1e-3 * v_sunsat * theta / A_h;
    double D_sat = 1e-3 * v_ssat / A_h;
    double C_unsat = D_unsat / H_b * a_PoverA_H;
    double C_sat = D_sat / H_b * a_IoverA_H;
    //double C_s = 0.0;		//!!!! Don't know what to use here !!!!
    //double C_T = C_s + C_unsat + C_sat;
    double C_T = C_p + C_unsat + C_sat;
    double Corr_evap = (C_T > 1.0) ? 1.0 / C_T : 1.0;
    double e_p = Corr_evap * C_p * e_pot;
    double e_sat = Corr_evap * C_sat * e_pot;
    double e_unsat = Corr_evap * C_unsat * e_pot;

    //Fluxes
    double q_pl = c_2 * (a_PoverA_H * RC + a_IoverA_H) * v_H * s_p;
    double q_sl = c_3 * v_ssat * a_IoverA_H;
    double snow_M = 0.8 * 0.0;	//Assume no snow for now
    double q_pu = a_PoverA_H * infrate * (1.0 - RC);
    double q_us = 1e3 * K_UNSAT * theta * a_PoverA_H;

    //!!!! This might be cleaner if H_b is converted into m !!!!

    //Channel transport
    ans[0] = -q;
    for (i = 0; i<num_parents; i++)	ans[0] += y_p[i * dim];
    ans[0] = K_Q * pow(q, lambda_1) * (ans[0] + c_4 * (q_pl + q_sl + snow_M));

    //Land surface storage
    ans[1] = (1.0 / 60.0) * (forcing_values[0] - q_pl - q_pu - e_p);

    //Water table
    ans[2] = (H_b < 1e-6) ? 0.0 : c_5 / deriv_a_I * (q_us - q_sl - e_sat);

    //Soil volumetric water content
    ans[3] = (v_sunsat < 1e-6) ? 0.0 : c_6 / v_sunsat * (q_pu - q_us - e_unsat + theta*(q_us - q_sl - e_sat));
    /*
    dump.ve[0] = e_p + e_sat + e_unsat;
    dump.ve[1] = q_pu;
    dump.ve[2] = q_pl;
    dump.ve[3] = q_sl;
    */
}


//Type 4
//Jacobian of the simple river system with soil properties function******************
//y_i has q_i (y_i[0]), s (y_i[1]), a (y_i[2]), v (y_i[3])
//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11   12     13     14  15  16  17  18  19
//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
//The numbering is:        0      1        2     3   4   5
//This uses the units and functions from September 18, 2011 document
void Jsimple_soil(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, double *ans)
{
    unsigned short i;

    double q = y_i[0];
    double s = y_i[1];
    double a = y_i[2];
    double v = y_i[3];

    double invtau = params[12];
    double epsilon = params[13];
    double c_1 = params[14];
    double c_4 = params[17];
    double c_5 = params[18];
    double c_6 = params[19];
    double lambda_1 = global_params[1];

    double qlambda_1 = pow(q, lambda_1);
    double s23 = pow(s, 2.0 / 3.0);
    double expterm = (v == 1.0 - a) ? 0.0 : exp(-(1 - a - v) / (1 - a));
    double denom = (s*(1 - a) + epsilon*(1 - a - v))*(s*(1 - a) + epsilon*(1 - a - v))*(s*(1 - a) + epsilon*(1 - a - v));

    //Calculate total inflow
    ans[0 * dim + 0] = 0.0;
    for (i = 0; i<num_parents; i++)
        ans[0 * dim + 0] += y_p[i * dim];

    //Calculate the jacobian
    ans[0 * dim + 0] = qlambda_1 * (lambda_1 * invtau * (ans[0 * dim + 0] - q + c_1*s23*s) / q - invtau);
    //ans[0 * dim + 0] = qlambda_1 * ( lambda_1 * invtau * (ans[0 * dim + 0] - q + c_1*s23*s + 100.0) / q - invtau );
    ans[0 * dim + 1] = 5 / 3 * c_1 * invtau * qlambda_1 * s23;
    ans[0 * dim + 2] = 0.0;
    ans[0 * dim + 3] = 0.0;

    ans[1 * dim + 0] = 0.0;
    if (denom == 0.0)
    {
        ans[1 * dim + 1] = c_5 * (1.0 - a);
        ans[1 * dim + 2] = 0.0;
        ans[1 * dim + 3] = 0.0;
    }
    else
    {
        ans[1 * dim + 1] = -5.0 / 3.0*c_4*s23 - c_5*(1 - a)*(epsilon*(1 - a - v))*(epsilon*(1 - a - v))*(-s*(1 - a) + epsilon*(1 - a - v)) / denom;
        ans[1 * dim + 2] = c_5*epsilon*epsilon*s*(1 - a - v)*((1 - a)*s*(1 + v - a) + epsilon*(1 - a - v)*(1 - a - v)) / denom;
        ans[1 * dim + 3] = 2 * c_5*epsilon*epsilon*s*s*(1 - a)*(1 - a)*(1 - a - v) / denom;
    }

    ans[2 * dim + 0] = 0.0;
    ans[2 * dim + 1] = 0.0;
    ans[2 * dim + 2] = (v == 1.0 - a) ? -2 * c_6*a / epsilon : c_6 / ((1 - a)*epsilon) * (v*(v - 1)*expterm + a*v*expterm) + 2 * c_6*a / epsilon;
    ans[2 * dim + 3] = c_6 / epsilon*(1 - a + v)*expterm;

    ans[3 * dim + 0] = 0.0;
    if (denom == 0)
    {
        ans[3 * dim + 1] = -c_5 / epsilon * (1 - a);
        ans[3 * dim + 2] = (a == 1.0) ? 0.0 : c_6 * v / epsilon * (1 - a - v) / (1 - a) * expterm;
        ans[3 * dim + 3] = -c_6 / epsilon * (1 + v - a) * expterm;
    }
    else
    {
        ans[3 * dim + 1] = c_5*epsilon*(1 - a)*(1 - a - v)*(1 - a - v)*(s - a*s - epsilon*(1 - a - v)) / denom;
        ans[3 * dim + 2] = (1 - a - v) / epsilon * (c_6 / (1 - a)*v*expterm + c_5*epsilon*epsilon*s*((1 - a)*s*(a - v - a) - epsilon*(1 - a - v)*(1 - a - v)) / denom);
        ans[3 * dim + 3] = -c_6 / epsilon*(1 + v - a)*expterm + 2 * c_5*epsilon*s*s*(1 - a)*(1 - a)*(-1 + a + v) / denom;
    }
}

//Type 101
void Robertson(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    double y1 = y_i[0];
    double y2 = y_i[1];
    double y3 = y_i[2];

    ans[0] = -.04 * y1 + 1e4 * y2 * y3;
    ans[1] = .04 * y1 - 1e4 * y2 * y3 - 3e7 * y2 * y2;
    ans[2] = 3e7 * y2 * y2;
}

//Type 101
void JRobertson(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, double *ans)
{
    double y2 = y_i[1];
    double y3 = y_i[2];

    ans[0 * dim + 0] = -.04;
    ans[0 * dim + 1] = 1e4 * y3;
    ans[0 * dim + 2] = 1e4 * y2;

    ans[1 * dim + 0] = .04;
    ans[1 * dim + 1] = -1e4 * y3 - 6e7 * y2;
    ans[1 * dim + 2] = -1e4 * y2;

    ans[2 * dim + 0] = 0.0;
    ans[2 * dim + 1] = 6e7 * y2;
    ans[2 * dim + 2] = 0.0;
}

//Type 105
//This is the model used in the Asynch summary paper
//Function for the simple river system with rainfall
//Calculates the flow using simple parameters, using only the flow q.
//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,c_1,c_2,c_3
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13  14  15 
//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
//The numbering is:        0      1        2     3   4   5
//This uses the units and functions from September 18, 2011 document
//y_i[0] = q, y_i[1] = s
void river_rainfall_summary(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;

    double q = y_i[0];
    double s = y_i[1];

    double invtau = params[12];
    double c_1 = params[13];
    double c_2 = params[14];
    double c_3 = params[15];
    double lambda_1 = global_params[1];

    double s53 = pow(s, 5.0 / 3.0);

    //Flux equation (y_i[0])
    ans[0] = -q + c_1 * s53;
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];
}
    //Ponded water equation (y_i[1])
//The numbering is:	0   1   2    3     4
//Order of global_params: v_0,lambda_1,lambda_2
//The numbering is:        0      1        2
void Tiling(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;

    double lambda_1 = global_params[1];
    double A_h = params[2];	//[m^2]
    double invtau = params[3];	//[1/min]
    double c_1 = params[4];

    double hillslope_flux = forcing_values[0] * c_1;	//[m^3/s]

    double q = y_i[0];		//[m^3/s]

    //Discharge
    ans[0] = -q + hillslope_flux;
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];
}
