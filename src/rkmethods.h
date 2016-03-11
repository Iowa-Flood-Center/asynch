#ifndef RKMETHODS_H
#define RKMETHODS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mathmethods.h"
#include "structs.h"
#include "problems.h"
#include "riversys.h"
#include "sort.h"
//#include "cblas.h"
//#include "clapack.h"
#include "definetype.h"

extern int my_rank;
extern int np;

double InitialStepSize(double t,Link* link_i,UnivVars* GlobalVars,TempStorage* workspace);
void BackupParent(Link* currentp,double time,VEC** backup,UnivVars* GlobalVars);
void store_k(VEC** full_k,VEC** k,unsigned int s,unsigned int* dense_indices,unsigned int num_dense);

//RKSolver methods
//int ExplicitRKSolverLinear(Link* link_i,UnivVars* GlobalVars,int* assignments,FILE* outputfile,TempStorage* workspace);
int ExplicitRKSolver(Link* link_i,UnivVars* GlobalVars,int* assignments,short int print_flag,FILE* outputfile,ConnData* conninfo,Forcing** forcings,TempStorage* workspace);
int ExplicitRKIndex1SolverDam(Link* link_i,UnivVars* GlobalVars,int* assignments,short int print_flag,FILE* outputfile,ConnData* conninfo,Forcing** forcings,TempStorage* workspace);
int ExplicitRKIndex1Solver(Link* link_i,UnivVars* GlobalVars,int* assignments,short int print_flag,FILE* outputfile,ConnData* conninfo,Forcing** forcings,TempStorage* workspace);
int ExplicitRKSolverDiscont(Link* link_i,UnivVars* GlobalVars,int* assignments,short int print_flag,FILE* outputfile,ConnData* conninfo,Forcing** forcings,TempStorage* workspace);
int RadauRKSolver(Link* link_i,UnivVars* GlobalVars,int* assignments,short int print_flag,FILE* outputfile,ConnData* conninfo,Forcing** forcings,TempStorage* workspace);
int ExplicitRKSolver_DataAssim(Link* link_i,UnivVars* GlobalVars,int* assignments,short int print_flag,FILE* outputfile,ConnData* conninfo,Forcing** forcings,TempStorage* workspace);

//Forced solution methods
int ForcedSolutionSolver(Link* link_i,UnivVars* GlobalVars,int* assignments,short int print_flag,FILE* outputfile,ConnData* conninfo,Forcing** forcings,TempStorage* workspace);

//RK methods data
RKMethod* RKDense3_2();
void RKDense3_2_b(double theta,VEC* b);
void RKDense3_2_bderiv(double theta,VEC* b);
RKMethod* TheRKDense4_3();
void TheRKDense4_3_b(double theta,VEC* b);
RKMethod* DOPRI5_dense();
void DOPRI5_b(double theta,VEC* b);
void DOPRI5_bderiv(double theta,VEC* b);
RKMethod* RadauIIA3_dense();
void RadauIIA3_b(double theta,VEC* b);

#endif

