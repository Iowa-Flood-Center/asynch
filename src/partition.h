#ifndef PARTITION_H
#define PARTITION_H

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include "structs.h"
#include "system.h"
#include "sort.h"
//#include "metis.h"

extern int my_rank;
extern int np;

int* Partition_System_By_Leaves(Link** sys,unsigned int N,Link** leaves,unsigned int numleaves,unsigned int** my_sys,unsigned int* my_N,TransData* my_data,short int *getting);
int* Partition_System_By_Leaves_2(Link** sys,unsigned int N,Link** leaves,unsigned int numleaves,unsigned int** my_sys,unsigned int* my_N,TransData* my_data,short int *getting);
//int* Partition_METIS_Traditional(Link** sys,unsigned int N,Link** leaves,unsigned int numleaves,unsigned int** my_sys,unsigned int* my_N,TransData* my_data,short int *getting,UnivVars* GlobalVars);
//int* Partition_METIS_RainChanges(Link** sys,unsigned int N,Link** leaves,unsigned int numleaves,unsigned int** my_sys,unsigned int* my_N,TransData* my_data,short int *getting,UnivVars* GlobalVars);
//int* Partition_METIS_RainVolume(Link** sys,unsigned int N,Link** leaves,unsigned int numleaves,unsigned int** my_sys,unsigned int* my_N,TransData* my_data,short int *getting,UnivVars* GlobalVars);

#endif
