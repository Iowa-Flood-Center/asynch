#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <io.h>
#include <processdata.h>

//Creates an OutputFunc object
void OutputFunc_Init(
    unsigned short hydros_loc_flag,
    unsigned short peaks_loc_flag,
    unsigned short dump_loc_flag,
    OutputFunc* output_func)
{
    //Temporary Calculations
    output_func->PrepareTempOutput = &PrepareTempFiles;

    //Prepare Final Time Series Output
    if (hydros_loc_flag == 3)
    {
#if defined(HAVE_POSTGRESQL)

        output_func->PrepareOutput = &PrepareDatabaseTable;

#else //HAVE_POSTGRESQL

        if (my_rank == 0)	printf("Error: Asynch was build without PostgreSQL support.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

#endif //HAVE_POSTGRESQL
    }
    else
        output_func->PrepareOutput = NULL;

    //Prepare Peakflow Output
    if (peaks_loc_flag == 1)
        output_func->PreparePeakflowOutput = &PreparePeakFlowFiles;
    else
        output_func->PreparePeakflowOutput = NULL;

    //Create Final Time Series Output
    if (hydros_loc_flag == 1 || hydros_loc_flag == 2 || hydros_loc_flag == 4 || hydros_loc_flag == 5 || hydros_loc_flag == 6)
        output_func->CreateOutput = &DumpTimeSerieFile;
    else if (hydros_loc_flag == 3)
    {
#if defined(HAVE_POSTGRESQL)

        output_func->CreateOutput = &DumpTimeSerieDB;

#else //HAVE_POSTGRESQL

        if (my_rank == 0)	printf("Error: Asynch was build without PostgreSQL support.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

#endif //HAVE_POSTGRESQL
    }       
    else
        output_func->CreateOutput = NULL;

    //Create Peakflow Output
    if (peaks_loc_flag == 1)
        output_func->CreatePeakflowOutput = &DumpPeakFlowText;
    else if (peaks_loc_flag == 2)
    {
#if defined(HAVE_POSTGRESQL)

        output_func->CreatePeakflowOutput = &DumpPeakFlowDB;

#else //HAVE_POSTGRESQL

        if (my_rank == 0)	printf("Error: Asynch was build without PostgreSQL support.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

#endif //HAVE_POSTGRESQL
    }        
    else
        output_func->CreatePeakflowOutput = NULL;

    //Set data dump routines
    if (dump_loc_flag == 1)
        output_func->CreateSnapShot = &DumpStateText;
    else if (dump_loc_flag == 2)
    {
#if defined(HAVE_POSTGRESQL)

        output_func->CreateSnapShot = &DumpStateDB;

#else //HAVE_POSTGRESQL

        if (my_rank == 0)	printf("Error: Asynch was build without PostgreSQL support.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

#endif //HAVE_POSTGRESQL
    }
    else if ((dump_loc_flag == 3) || (dump_loc_flag == 4))
        output_func->CreateSnapShot = &DumpStateH5;
    else
        output_func->CreateSnapShot = NULL;
}

//!!!! This really sucks. Is there any way to improve it? !!!!
//Writes a value to an ASCII file
void WriteValue(FILE* outputfile, const char* specifier, char* data_storage, short int data_type, char* delim)
{
    switch (data_type)
    {
    case ASYNCH_DOUBLE:
        fprintf(outputfile, specifier, *(double*)data_storage);
        break;
    case ASYNCH_INT:
        fprintf(outputfile, specifier, *(int*)data_storage);
        break;
    case ASYNCH_FLOAT:
        fprintf(outputfile, specifier, *(float*)data_storage);
        break;
    case ASYNCH_SHORT:
        fprintf(outputfile, specifier, *(short int*)data_storage);
        break;
    case ASYNCH_CHAR:
        fprintf(outputfile, specifier, *(char*)data_storage);
        break;
    default:
        printf("[%i]: Error: Writing bad value to an ascii file (%hi).\n", my_rank, data_type);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    fprintf(outputfile, "%s", delim);
}

unsigned int WriteStep(Output *output, unsigned int num_outputs, FILE* outputfile, unsigned int id, double t, double *y, unsigned int num_dof, long int* pos_offset)
{
    unsigned int i;

    long int total_written = 0;

    //Set file to current position
    //fsetpos(outputfile,pos);
    if (pos_offset)
        fseek(outputfile, *pos_offset, SEEK_SET);

    //Write the step
    for (i = 0; i < num_outputs; i++)
    {
        switch (output[i].type)	//!!!! Get rid of this. Try char[] and output_sizes. !!!!
        {
        case ASYNCH_INT:
        {
            int output_i = output[i].callback.out_int(id, t, y, num_dof);
            fwrite(&output_i, sizeof(int), 1, outputfile);
            break;
        }
        case ASYNCH_DOUBLE:
        {
            double output_d = output[i].callback.out_double(id, t, y, num_dof);
            fwrite(&output_d, sizeof(double), 1, outputfile);
            break;
        }
        case ASYNCH_FLOAT:
        {
            float output_f = output[i].callback.out_float(id, t, y, num_dof);
            fwrite(&output_f, sizeof(float), 1, outputfile);
            break;
        }
        default:
            printf("[%i]: Error: Invalid output %s (%i).\n", my_rank, output[i].specifier, output[i].type);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        total_written += output[i].size;
        //if(pos_offset)	*pos_offset += globals->output_sizes[i];
    }

    if (pos_offset)
        *pos_offset += total_written;
    
    return total_written;
}


unsigned int CatBinaryToString(char* submission, const char* specifier, void* data_storage, short int data_type, char* delim)
{
    unsigned int written = 0;

    switch (data_type)
    {
    case ASYNCH_DOUBLE:
        written = sprintf(submission, specifier, *(double*)data_storage);
        break;
    case ASYNCH_FLOAT:
        written = sprintf(submission, specifier, *(float*)data_storage);
        break;
    case ASYNCH_SHORT:
        written = sprintf(submission, specifier, *(short*)data_storage);
        break;
    case ASYNCH_INT:
        written = sprintf(submission, specifier, *(int*)data_storage);
        break;
    default:
        printf("[%i]: Error: Writing bad value to an ascii file (%hi).\n", my_rank, data_type);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    sprintf(&(submission[written++]), "%s", delim);

    return written;
}


