#ifndef SYSTEM_H
#define SYSTEM_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "mathmethods.h"
#include "structs.h"

extern int my_rank;
extern int np;

//Discontinuity list
unsigned int Insert_Discontinuity(double time, unsigned int start, unsigned int end, unsigned int* count, unsigned int size, double* array, unsigned int id);
void Insert_SendDiscontinuity(double time, unsigned int order, unsigned int* count, unsigned int size, double* array, unsigned int*order_array, unsigned int id);

//Destructors
void Destroy_Link(Link* link_i, unsigned int list_length, int rkd_flag, Forcing* forcings, GlobalVars* GlobalVars);
void ForcingData_Free(ForcingData** forcing_buff);
void Destroy_RKMethod(RKMethod* method);
void Destroy_ErrorData(ErrorData* error);
void Destroy_Node(RKSolutionNode* node, unsigned short int s);
void Destroy_List(RKSolutionList* list, unsigned int list_length);
void Destroy_Workspace(TempStorage* workspace, unsigned short int s, unsigned short int max_parents);
void Destroy_UnivVars(GlobalVars* GlobalVars);

//RKSolutionList methods
RKSolutionNode* New_Step(RKSolutionList* list);
void Undo_Step(RKSolutionList* list);
void Remove_Head_Node(RKSolutionList* list);
RKSolutionList* Create_List(VEC y0, double t0, int dim, unsigned int dense_dim, unsigned short int s, unsigned int list_length);

//Workspace methods
TempStorage* Create_Workspace(unsigned int dim, unsigned short int s, unsigned short int max_parents);

#endif
