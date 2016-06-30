#ifndef DATA_TYPES_H
#define DATA_TYPES_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#define ASYNCH_BAD_TYPE -1
#define ASYNCH_CHAR 0
#define ASYNCH_SHORT 1
#define ASYNCH_INT 2
#define ASYNCH_FLOAT 3
#define ASYNCH_DOUBLE 4
#define ASYNCH_NUM_DATA_TYPES 5


typedef struct DataTypes
{
	int (*is_equal[ASYNCH_NUM_DATA_TYPES])(void*,void*);
} DataTypes;

void Init_DataTypes(DataTypes* dt_info);
void Free_DataTypes(DataTypes* dt_info);

//Comparison functions ***********************************************************************************************
int IsEqualInt(void* a,void* b);
int IsEqualShort(void* a,void* b);
int IsEqualDouble(void* a,void* b);
int IsEqualFloat(void* a,void* b);
int IsEqualChar(void* a,void* b);

#endif

