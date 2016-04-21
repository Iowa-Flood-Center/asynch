#include "data_types.h"

void Init_DataTypes(DataTypes* dt_info)
{
	dt_info->is_equal[ASYNCH_CHAR] = &IsEqualChar;
	dt_info->is_equal[ASYNCH_SHORT] = &IsEqualShort;
	dt_info->is_equal[ASYNCH_INT] = &IsEqualInt;
	dt_info->is_equal[ASYNCH_FLOAT] = &IsEqualFloat;
	dt_info->is_equal[ASYNCH_DOUBLE] = &IsEqualDouble;
}

void Free_DataTypes(DataTypes* dt_info)
{
    for (int i = 0; i < ASYNCH_NUM_DATA_TYPES; i++)
        dt_info->is_equal[i] = NULL;
}

//Comparison functions ***********************************************************************************************
int IsEqualInt(void* a,void* b)
{
	return *(int*)a == *(int*)b;
}

int IsEqualShort(void* a,void* b)
{
	return *(short int*)a == *(short int*)b;
}

int IsEqualDouble(void* a,void* b)
{
	return (fabs( *(double*)a - *(double*)b ) < 1e-10);
}

int IsEqualFloat(void* a,void* b)
{
	return (fabs( *(float*)a == *(float*)b ) < 1e-10);
}

int IsEqualChar(void* a,void* b)
{
	return *(char*)a == *(char*)b;
}


