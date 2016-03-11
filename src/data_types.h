#ifndef DATA_TYPES_H
#define DATA_TYPES_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ASYNCH_BAD_TYPE -1
#define ASYNCH_CHAR 0
#define ASYNCH_SHORT 1
#define ASYNCH_INT 2
#define ASYNCH_FLOAT 3
#define ASYNCH_DOUBLE 4

typedef struct data_types
{
	int (**is_equal)(void*,void*);
} data_types;

data_types* Init_DataTypes();
void Free_DataTypes(data_types** dt_info);

//Comparison functions ***********************************************************************************************
int IsEqualInt(void* a,void* b);
int IsEqualShort(void* a,void* b);
int IsEqualDouble(void* a,void* b);
int IsEqualFloat(void* a,void* b);
int IsEqualChar(void* a,void* b);

#endif

