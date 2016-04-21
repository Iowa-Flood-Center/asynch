#ifndef SORT_H
#define SORT_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "structs.h"

//Merge sort routines
void merge_sort(Link** array,unsigned int size,unsigned int idx);
void merge_sort_distance(Link** sys,unsigned int* array,unsigned int size);
void merge_sort_ids(unsigned int** array,unsigned int size);
void merge_sort_1D(unsigned int* array,unsigned int size);

//Binary search routines
unsigned int find_link_by_id(unsigned int id,Link** sys,unsigned int N);
unsigned int find_link_by_idtoloc(unsigned int id,unsigned int** id_to_loc,unsigned int N);

#endif //SORT_H
