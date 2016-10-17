#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <stdio.h>
#include <string.h>

#include "io.h"

//Creates an io object
InOut* BuildIO(UnivVars* GlobalVars)
{
    InOut* output_data = (InOut*) malloc(sizeof(InOut));

	//Temporary Calculations
	output_data->PrepareTempOutput = &PrepareTempFiles;

	//Prepare Final Time Series Output
	if(GlobalVars->hydros_loc_flag == 3)	output_data->PrepareOutput = &PrepareDatabaseTable;
	else					output_data->PrepareOutput = NULL;

	//Prepare Peakflow Output
	if(GlobalVars->peaks_loc_flag == 1)	output_data->PreparePeakflowOutput = &PreparePeakFlowFiles;
	else					output_data->PreparePeakflowOutput = NULL;

	//Create Final Time Series Output
	if(GlobalVars->hydros_loc_flag == 1 || GlobalVars->hydros_loc_flag == 2 || GlobalVars->hydros_loc_flag == 4)
		output_data->CreateOutput = &Process_Data;
	else if(GlobalVars->hydros_loc_flag == 3)
		output_data->CreateOutput = &UploadHydrosDB;
	else
		output_data->CreateOutput = NULL;

	//Create Peakflow Output
	if(GlobalVars->peaks_loc_flag == 1)
		output_data->CreatePeakflowOutput = &DumpPeakFlowData;
	else if(GlobalVars->peaks_loc_flag == 2)
		output_data->CreatePeakflowOutput = &UploadPeakFlowData;
	else
		output_data->CreatePeakflowOutput = NULL;

	//Set data dump routines
    if (GlobalVars->dump_loc_flag == 1)
        output_data->CreateSnapShot = &DataDump2;
    else if (GlobalVars->dump_loc_flag == 2)
        output_data->CreateSnapShot = &UploadDBDataDump;
    else if ((GlobalVars->dump_loc_flag == 3) || (GlobalVars->dump_loc_flag == 4))
        output_data->CreateSnapShot = &DataDumpH5;
	else
		output_data->CreateSnapShot = NULL;

	return output_data;
}


//Reads a .dbc file and creates a corresponding database connection.
//This does NOT connect to the database.
//Returns NULL if there was an error.
ConnData* ReadDBC(char* filename,unsigned int string_size)
{
	unsigned int i = 0,j = 0,errorcode = 0;
	ConnData* conninfo = NULL;
	char c,*connectstring = (char*) malloc(string_size*sizeof(char));

	if(my_rank == 0)
	{
		FILE* input = fopen(filename,"r");

		if(!input)
		{
			printf("Error opening .dbc file %s.\n",filename);
			goto finish;
		}
		if(CheckWinFormat(input))
		{
			printf("Error: File %s appears to be in Windows format. Try converting to unix format using 'dos2unix' at the command line.\n",filename);
			fclose(input);
			goto finish;
		}

		//Read connection information
		//Currently, this expects 4 things like:
		//dbname=rm_model host=s-iihr58.iihr.uiowa.edu user=***REMOVED*** password=my_pass
		fgets(connectstring,string_size,input);
		conninfo = CreateConnData(connectstring);
		if(!conninfo)
		{
			fclose(input);
			goto finish;
		}

		//Get number of queries
		if (fscanf(input,"%u",&(conninfo->num_queries)) == EOF)
        {
            printf("[%i]: Error: failed to parse file.\n", my_rank);
            exit(EXIT_FAILURE);
        }
		conninfo->queries = (char**) malloc(conninfo->num_queries*sizeof(char*));
		for(i=0;i<conninfo->num_queries;i++)
			conninfo->queries[i] = (char*) malloc(ASYNCH_MAX_QUERY_SIZE*sizeof(char));

		//Get queries. They are delineated by a ;
		for(j=0;j<conninfo->num_queries;j++)
		{
			//Get rid of whitespace
			c = fgetc(input);
			while(c != EOF && (c == ' ' || c == '\n' || c == '\t'))	c = fgetc(input);
			if(c == EOF)
			{
				printf("[%i]: Warning: did not see %u queries in %s.\n",my_rank,conninfo->num_queries,filename);
				break;
			}

			//Read in query
			for(i=0;i<ASYNCH_MAX_QUERY_SIZE-2 && c != ';' && c != EOF;i++)
			{
				conninfo->queries[j][i] = c;
				c = fgetc(input);
			}

			//Check for problems and put stuff on the end
			if(i == ASYNCH_MAX_QUERY_SIZE)
				printf("[%i]: Warning: query %u is too long in %s.\n",my_rank,j,filename);
			else if(c == EOF)
			{
				printf("[%i]: Warning: did not see %u queries in %s.\n",my_rank,conninfo->num_queries,filename);
				break;
			}
			else
			{
				conninfo->queries[j][i] = ';';
				conninfo->queries[j][i+1] = '\0';
			}
		}

		fclose(input);

		//Get string length for other procs
		j = strlen(connectstring) + 1;
	}

	//Check if an error occurred
	finish:
	MPI_Bcast(&errorcode,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);

	//Transfer info from .dbc file
	if(!errorcode)
	{
		MPI_Bcast(&j,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
		MPI_Bcast(connectstring,j,MPI_CHAR,0,MPI_COMM_WORLD);
		if(my_rank != 0)
			conninfo = CreateConnData(connectstring);
		MPI_Bcast(&(conninfo->num_queries),1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
		if(my_rank == 0)
		{
			for(i=0;i<conninfo->num_queries;i++)
			{
				j = strlen(conninfo->queries[i]) + 1;
				MPI_Bcast(&j,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
				MPI_Bcast(conninfo->queries[i],j,MPI_CHAR,0,MPI_COMM_WORLD);
			}
		}
		else
		{
			conninfo->queries = (char**) malloc(conninfo->num_queries*sizeof(char*));
			for(i=0;i<conninfo->num_queries;i++)
			{
				conninfo->queries[i] = (char*) malloc(ASYNCH_MAX_QUERY_SIZE*sizeof(char));
				MPI_Bcast(&j,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
				MPI_Bcast(conninfo->queries[i],j,MPI_CHAR,0,MPI_COMM_WORLD);
			}
		}

	}

	free(connectstring);
	return conninfo;
}

//!!!! This really sucks. Is there any way to improve it? !!!!
//Writes a value to an ASCII file
void WriteValue(FILE* outputfile,char* specifier,char* data_storage,short int data_type,char* delim)
{
	switch(data_type)
	{
		case ASYNCH_DOUBLE:
			fprintf(outputfile,specifier,*(double*)data_storage);
			break;
		case ASYNCH_INT:
			fprintf(outputfile,specifier,*(int*)data_storage);
			break;
		case ASYNCH_FLOAT:
			fprintf(outputfile,specifier,*(float*)data_storage);
			break;
		case ASYNCH_SHORT:
			fprintf(outputfile,specifier,*(short int*)data_storage);
			break;
		case ASYNCH_CHAR:
			fprintf(outputfile,specifier,*(char*)data_storage);
			break;
		default:
			printf("[%i]: Error: Writing bad value to an ascii file (%hi).\n",my_rank,data_type);
			MPI_Abort(MPI_COMM_WORLD,1);
	}

	fprintf(outputfile,delim);
}

//void WriteStep(double t,VEC* y,UnivVars* GlobalVars,VEC* params,unsigned int state,FILE* outputfile,void* user,fpos_t* pos)
unsigned int WriteStep(double t,VEC y,UnivVars* GlobalVars,VEC params,unsigned int state,FILE* outputfile,void* user,long int* pos_offset)
{
	unsigned int i;

	double output_d;
	//float output_f;
	//short int output_s;
	int output_i;
	//char output_c;
	long int total_written = 0;

	//Set file to current position
	//fsetpos(outputfile,pos);
	if(pos_offset)	fseek(outputfile,*pos_offset,SEEK_SET);

	//Write the step
	for(i=0;i<GlobalVars->num_print;i++)
	{
		switch(GlobalVars->output_types[i])	//!!!! Get rid of this. Try char[] and output_sizes. !!!!
		{
			case ASYNCH_DOUBLE:
				output_d = (GlobalVars->outputs_d[i])(t,y,GlobalVars->global_params,params,state,user);
				fwrite(&output_d,sizeof(double),1,outputfile);
				break;
			case ASYNCH_INT:
				output_i = (GlobalVars->outputs_i[i])(t,y,GlobalVars->global_params,params,state,user);
				fwrite(&output_i,sizeof(int),1,outputfile);
				break;
			default:
				printf("[%i]: Error: Invalid output %s (%hi).\n",my_rank,GlobalVars->output_specifiers[i],GlobalVars->output_types[i]);
				MPI_Abort(MPI_COMM_WORLD,1);
		}
		total_written += GlobalVars->output_sizes[i];
		//if(pos_offset)	*pos_offset += GlobalVars->output_sizes[i];
	}

	if(pos_offset)	*pos_offset += total_written;
	return total_written;
}

unsigned int CatBinaryToString(char* submission,char* specifier,char* data_storage,short int data_type,char* delim)
{
	unsigned int written = 0;

	switch(data_type)
	{
		case ASYNCH_DOUBLE:
			written = sprintf(submission,specifier,*(double*)data_storage);
			break;
		case ASYNCH_INT:
			written = sprintf(submission,specifier,*(int*)data_storage);
			break;
		default:
			printf("[%i]: Error: Writing bad value to an ascii file (%hi).\n",my_rank,data_type);
			MPI_Abort(MPI_COMM_WORLD,1);
	}

	sprintf(&(submission[written++]),delim);

	return written;
}


