#include "modeloutputs.h"

void SetOutputFunctions(char* outputname,char* specifier,unsigned int* states_used,unsigned int* num_states_used,short int* output_size,short int* output_type,int (**output_i)(double,VEC*,VEC*,VEC*,int,void*),double (**output_d)(double,VEC*,VEC*,VEC*,int,void*))
{
	if(strcmp(outputname,"Time") == 0)
	{
		*output_type = ASYNCH_DOUBLE;
		*output_d = &Output_Time;
	}
	else if(strcmp(outputname,"State0") == 0)
	{
		states_used[(*num_states_used)++] = 0;
		*output_type = ASYNCH_DOUBLE;
		*output_d = &Output_State0;
	}
	else if(strcmp(outputname,"State1") == 0)
	{
		states_used[(*num_states_used)++] = 1;
		*output_type = ASYNCH_DOUBLE;
		*output_d = &Output_State1;
	}
	else if(strcmp(outputname,"State2") == 0)
	{
		states_used[(*num_states_used)++] = 2;
		*output_type = ASYNCH_DOUBLE;
		*output_d = &Output_State2;
	}
	else if(strcmp(outputname,"State3") == 0)
	{
		states_used[(*num_states_used)++] = 3;
		*output_type = ASYNCH_DOUBLE;
		*output_d = &Output_State3;
	}
	else if(strcmp(outputname,"State4") == 0)
	{
		states_used[(*num_states_used)++] = 4;
		*output_type = ASYNCH_DOUBLE;
		*output_d = &Output_State4;
	}
	else if(strcmp(outputname,"State5") == 0)
	{
		states_used[(*num_states_used)++] = 5;
		*output_type = ASYNCH_DOUBLE;
		*output_d = &Output_State5;
	}
	else if(strcmp(outputname,"State6") == 0)
	{
		states_used[(*num_states_used)++] = 6;
		*output_type = ASYNCH_DOUBLE;
		*output_d = &Output_State6;
	}
	else if(strcmp(outputname,"State7") == 0)
	{
		states_used[(*num_states_used)++] = 7;
		*output_type = ASYNCH_DOUBLE;
		*output_d = &Output_State7;
	}
	else if(strcmp(outputname,"TimeI") == 0)
	{
		*output_type = ASYNCH_INT;
		*output_i = &Output_Time_Int;
	}
/*
	else if(strcmp(outputname,"LinkID") == 0)
	{
		*output_type = ASYNCH_INT;
		*output_i = &Output_Linkid;
	}
*/
	else	//Must be a user defined output
	{
		*output_type = -1;
		return;
		//printf("[%i]: Error setting output: bad output function name %s.\n",my_rank,outputname);
		//MPI_Abort(MPI_COMM_WORLD,1);
	}

	*output_size = GetByteSize(*output_type);
	GetSpecifier(specifier,*output_type);
}

void SetPeakflowOutputFunctions(char* outputname,void (**peak_output)(unsigned int,double,VEC*,VEC*,VEC*,double,unsigned int,void*,char*))
{
	if(strcmp(outputname,"Classic") == 0)
		*peak_output = &OutputPeakflow_Classic_Format;
	else if(strcmp(outputname,"Forecast") == 0)
		*peak_output = &OutputPeakflow_Forecast_Format;
	else	//Must be a user defined output
	{
		*peak_output = NULL;
		//printf("[%i]: Bad peakflow function %s.\n",my_rank,outputname);
		//MPI_Abort(MPI_COMM_WORLD,1);
	}
}

//ASYNCH_CHAR = 0
//ASYNCH_SHORT = 1
//ASYNCH_INT = 2
//ASYNCH_FLOAT = 3
//ASYNCH_DOUBLE = 4
short int GetByteSize(short int type)
{
	const short int data_sizes[] = {sizeof(char), sizeof(short int), sizeof(int), sizeof(float), sizeof(double)};
	return data_sizes[type];
}

void GetSpecifier(char* specifier,short int type)
{
	const char *data_specifiers[5] = {"%c", "%hi", "%i", "%.6e", "%.12e"};
	memcpy(specifier,data_specifiers[type],strlen(data_specifiers[type])+1);
}

unsigned int CalcTotalOutputSize(UnivVars* GlobalVars)
{
	unsigned int i,total = 0;

	for(i=0;i<GlobalVars->num_print;i++)
		total += GlobalVars->output_sizes[i];

	return total;
}

//Returns 0 if an output is not set, 1 if all outputs are good to go.
int OutputsSet(UnivVars* GlobalVars)
{
	unsigned int i;

	for(i=0;i<GlobalVars->num_print;i++)
	{
		switch(GlobalVars->output_types[i])
		{
			case ASYNCH_BAD_TYPE:
				return 0;
			case ASYNCH_DOUBLE:
				if(!(GlobalVars->outputs_d[i]))	return 0;
				break;
			case ASYNCH_INT:
				if(!(GlobalVars->outputs_i[i]))	return 0;
				break;
			default:
				printf("[%i]: Error: Bad output type (%hi) encountered while checking if outputs set.\n",my_rank,GlobalVars->output_types[i]);
				MPI_Abort(MPI_COMM_WORLD,1);
		}
	}

	return 1;
}

//Returns 0 if peakflow output is not set, 1 if it is.
int PeakflowOutputsSet(UnivVars* GlobalVars)
{
	if(GlobalVars->peakflow_function_name && !GlobalVars->peakflow_output)	return 0;
	return 1;
}


//Output functions ***************************************************************************************************

double Output_Time(double t,VEC* y_i,VEC* global_params,VEC* params,int state,void* user)
{
	return t;
}

double Output_State0(double t,VEC* y_i,VEC* global_params,VEC* params,int state,void* user)
{
	return y_i->ve[0];
}

double Output_State1(double t,VEC* y_i,VEC* global_params,VEC* params,int state,void* user)
{
	return y_i->ve[1];
}

double Output_State2(double t,VEC* y_i,VEC* global_params,VEC* params,int state,void* user)
{
	return y_i->ve[2];
}

double Output_State3(double t,VEC* y_i,VEC* global_params,VEC* params,int state,void* user)
{
	return y_i->ve[3];
}

double Output_State4(double t,VEC* y_i,VEC* global_params,VEC* params,int state,void* user)
{
	return y_i->ve[4];
}

double Output_State5(double t,VEC* y_i,VEC* global_params,VEC* params,int state,void* user)
{
	return y_i->ve[5];
}

double Output_State6(double t,VEC* y_i,VEC* global_params,VEC* params,int state,void* user)
{
	return y_i->ve[6];
}

double Output_State7(double t,VEC* y_i,VEC* global_params,VEC* params,int state,void* user)
{
	return y_i->ve[7];
}

int Output_Time_Int(double t,VEC* y_i,VEC* global_params,VEC* params,int state,void* user)
{
	return (int) (round(t) + 0.1);
}


//Peakflow output functions***********************************************************************************

void OutputPeakflow_Classic_Format(unsigned int ID,double peak_time,VEC* peak_value,VEC* params,VEC* global_params,double conversion,unsigned int area_idx,void* user,char* buffer)
{
	sprintf(buffer,"%u %.4f %.8f %.8f\n",ID,conversion*params->ve[area_idx],peak_time,peak_value->ve[0]);
}

void OutputPeakflow_Forecast_Format(unsigned int ID,double peak_time,VEC* peak_value,VEC* params,VEC* global_params,double conversion,unsigned int area_idx,void* user,char* buffer)
{
	unsigned int offset = *(unsigned int*)user;
	//sprintf(buffer,"%u %u %.12e %u NULL\n",ID,offset + (unsigned int)(peak_time*60 + .1),peak_value->ve[0],offset);
	sprintf(buffer,"%u %u %.12e %u\n",ID,offset + (unsigned int)(peak_time*60 + .1),peak_value->ve[0],offset);
}

