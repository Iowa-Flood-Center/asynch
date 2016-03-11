#include "data_types.h"

data_types* Init_DataTypes()
{
	unsigned int i,num_data_types = 5;
	data_types* dt_info = (data_types*) malloc(sizeof(data_types));

	dt_info->is_equal = (int (**)(void*,void*)) malloc(num_data_types*sizeof(int (*)(void*,void*)));
	dt_info->is_equal[ASYNCH_CHAR] = &IsEqualChar;
	dt_info->is_equal[ASYNCH_SHORT] = &IsEqualShort;
	dt_info->is_equal[ASYNCH_INT] = &IsEqualInt;
	dt_info->is_equal[ASYNCH_FLOAT] = &IsEqualFloat;
	dt_info->is_equal[ASYNCH_DOUBLE] = &IsEqualDouble;

	return dt_info;
}

void Free_DataTypes(data_types** dt_info)
{
	if(*dt_info)
	{
		free((*dt_info)->is_equal);
		free(*dt_info);
		*dt_info = NULL;
	}
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


