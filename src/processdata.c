#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#if defined(HAVE_UNISTD_H)
#include <unistd.h>
#endif

#if !defined(_MSC_VER)
#define ASYNCH_SLEEP sleep
#else
#include <windows.h>
#define ASYNCH_SLEEP Sleep
#endif

#include <mpi.h>

#if defined(HAVE_POSTGRESQL)
#include <libpq-fe.h>
#endif

#if defined(HAVE_HDF5)
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

#include <minmax.h>
#include <outputs.h>
#include <io.h>
#include <processdata.h>
#include <blas.h>
#include <db.h>

/// Returns the H5 Type that match the Asynch type
static hid_t Get_H5_Type(enum AsynchTypes type)
{
    static hid_t h5_types[ASYNCH_NUM_DATA_TYPES];
    h5_types[0] = H5T_NATIVE_CHAR;
    h5_types[1] = H5T_NATIVE_SHORT;
    h5_types[2] = H5T_NATIVE_INT;
    h5_types[3] = H5T_NATIVE_FLOAT;
    h5_types[4] = H5T_NATIVE_DOUBLE;

    return h5_types[type];
}

//Reads the results stored in temporary files and outputs them conveniently to a new file.
//Use this for parallel implementations.
//There could be problems if there is lots of data in the temporary files. Improve on this in the future.
//Link* sys: The data for the riversystem.
//int* my_sys: The index in sys of all links assigned to this process.
//int my_N: The number of links assigned to this process.
//int save_size: Number of links for which data is to be written to file. This is the grand total for all processes.
//int my_save_size: Number of links for which data is to be written to file. This is the number for just this process.
//Return value = 0 means everything is good.
//1 means a database related error.
//2 means a file system related error.
//!!!! Is additional_temp needed? !!!!
int DumpTimeSerieFile(Link* sys, GlobalVars* globals, unsigned int N, unsigned int* save_list, unsigned int save_size, unsigned int my_save_size, const Lookup * const id_to_loc, int* assignments, char* additional_temp, char* additional_out, ConnData* conninfo, FILE** my_tempfile)
{
    char filename[ASYNCH_MAX_PATH_LENGTH];

    //Close the temp file, if open
    if (my_tempfile && *my_tempfile)
    {
        fclose(*my_tempfile);
        *my_tempfile = NULL;
    }

    //for (i = 0; i < size; i++)
    //    data_storage[i] = 0;

    //Find total size of a line in the temp files
    unsigned int line_size = CalcTotalOutputSize(globals);

    if (globals->hydros_loc_flag == 1)	//.dat
    {
        DumpTimeSerieDatFile(sys, globals, N, save_list, save_size, my_save_size, id_to_loc, assignments, additional_temp, additional_out);
    }
    else if (globals->hydros_loc_flag == 2)	//.csv
    {
        DumpTimeSerieCsvFile(sys, globals, N, save_list, save_size, my_save_size, id_to_loc, assignments, additional_temp, additional_out);
    }
#if defined(ASYNCH_HAVE_RADEK_PATENTED_COMPACT_BINARY_FORMAT_THAT_NO_ONE_ELSE_CAN_READ)
    else if (globals->hydros_loc_flag == 4)	//Radek's patented compact binary files
    {
        size_t io_val;
        char outputfilename_index[ASYNCH_MAX_PATH_LENGTH];

        //Create output .rad and .irad files
        //if(my_rank == 0)
        {
            if (globals->print_par_flag == 1)
            {
                if (!additional_out)
                    sprintf(outputfilename, "%s", globals->hydros_loc_filename);
                else
                    sprintf(outputfilename, "%s_%s", globals->hydros_loc_filename, additional_out);
                for (i = 0; i < globals->num_global_params; i++)
                {
                    sprintf(filenamespace, "_%.4e", globals->global_params.ve[i]);
                    strcat(outputfilename, filenamespace);
                }

                sprintf(outputfilename_index, "%s_%i.irad", outputfilename, my_rank);
                sprintf(filenamespace, "_%i.rad", my_rank);
                strcat(outputfilename, filenamespace);
            }
            else
            {
                if (!additional_out)
                {
                    sprintf(outputfilename_index, "%s_%i.irad", globals->hydros_loc_filename, my_rank);
                    sprintf(outputfilename, "%s_%i.rad", globals->hydros_loc_filename, my_rank);
                }
                else
                {
                    sprintf(outputfilename_index, "%s_%s_%i.irad", globals->hydros_loc_filename, additional_out, my_rank);
                    sprintf(outputfilename, "%s_%s_%i.rad", globals->hydros_loc_filename, additional_out, my_rank);
                }
            }

            outputfile = fopen(outputfilename, "wb");
            if (!outputfile)
            {
                printf("[%i]: Error creating outputfile %s.\n", my_rank, outputfilename);
                return 2;
            }

            outputfile_index = fopen(outputfilename_index, "wb");
            if (!outputfile_index)
            {
                printf("[%i]: Error creating index outputfile %s.\n", my_rank, outputfilename_index);
                fclose(outputfile);
                return 2;
            }

            //Write index file
            for (i = 0; i < save_size; i++)
            {
                loc = find_link_by_idtoloc(save_list[i], id_to_loc, N);
                if (assignments[loc] == my_rank)
                {
                    io_val = fwrite(&(save_list[i]), sizeof(unsigned int), 1, outputfile_index);
                    while (io_val != 1)
                    {
                        printf("[%i]: Error writing id %u to .irad file. Trying again...\n", my_rank, save_list[i]);
                        ASYNCH_SLEEP(2);
                        io_val = fwrite(&(save_list[i]), sizeof(unsigned int), 1, outputfile_index);
                    }
                }
            }
            //fwrite(save_list,sizeof(unsigned int),save_size,outputfile_index);
            fclose(outputfile_index);
        }

        if (my_save_size)
        {
            //Open temporary file
            if (!additional_temp)
                sprintf(filename, "%s", globals->temp_filename);
            else
                sprintf(filename, "%s_%s", globals->temp_filename, additional_temp);				//!!!! When is additional_temp useful? Get rid of it? Global params?? !!!!
            inputfile = fopen(filename, "rb");
            if (!inputfile)
            {
                printf("\n[%i]: Error opening inputfile %s.\n", my_rank, filename);
                return 2;
            }

            //Move data to final output
            for (i = 0; i < save_size; i++)
            {
                loc = find_link_by_idtoloc(save_list[i], id_to_loc, N);
                proc = assignments[loc];
                current = &sys[loc];

                if (proc == my_rank)
                {
                    //Find the link in the temp file inputfile
                    fread(&id, sizeof(unsigned int), 1, inputfile);
                    fread(&total_spaces, sizeof(unsigned int), 1, inputfile);
                    while (id != save_list[i] && !feof(inputfile))
                    {
                        jump_size = (long long int) total_spaces * (long long int) line_size;	//Watch overflows!
                        fseek(inputfile, jump_size, SEEK_CUR);
                        fread(&id, sizeof(unsigned int), 1, inputfile);
                        fread(&total_spaces, sizeof(unsigned int), 1, inputfile);
                    }
                    if (feof(inputfile))
                    {
                        printf("\n[%i]: Error: could not find id %u in temp file %s.\n", my_rank, save_list[i], filename);
                        return 2;
                    }

                    //Read data in the temp file
                    //if(my_rank == 0)
                    {
                        for (k = 0; k < current->disk_iterations; k++)
                        {
                            for (m = 0; m < dim; m++)
                            {
                                fread(data_storage, globals->output_sizes[m], 1, inputfile);
                                io_val = fwrite(data_storage, globals->output_sizes[m], 1, outputfile);
                                while (io_val != 1)
                                {
                                    printf("[%i]: Error writing data for id %u to .rad file. Trying again...\n", my_rank, save_list[i]);
                                    ASYNCH_SLEEP(2);
                                    io_val = fwrite(data_storage, globals->output_sizes[m], 1, outputfile);
                                }
                            }
                        }
                    }

                    //Skip over the last unused space. This is done so that the file does not need to be rewound.
                    for (k = 0; k < total_spaces - current->disk_iterations; k++)
                        for (j = 0; j < dim; j++)
                            fread(data_storage, globals->output_sizes[j], 1, inputfile);
                }
            }
        }

        //Cleanup
        if (inputfile)	fclose(inputfile);
        if (outputfile)	fclose(outputfile);
    }
#endif //ASYNCH_HAVE_RADEK_PATENTED_COMPACT_BINARY_FORMAT_THAT_NO_ONE_ELSE_CAN_READ
    else if (globals->hydros_loc_flag == 5)	//.h5 packet
    {
        DumpTimeSerieH5File(sys, globals, N, save_list, save_size, my_save_size, id_to_loc, assignments, additional_temp, additional_out);
    }
    else if (globals->hydros_loc_flag == 6)	//.h5 array
    {
        DumpTimeSerieNcFile(sys, globals, N, save_list, save_size, my_save_size, id_to_loc, assignments, additional_temp, additional_out);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //Reopen the tempfile
    if (my_tempfile && my_save_size > 0)
    {
        if (additional_temp != NULL)
            sprintf(filename, "%s_%s", globals->temp_filename, additional_temp);
        else
            sprintf(filename, "%s", globals->temp_filename);
        *my_tempfile = fopen(filename, "r+b");
        if (!*my_tempfile)
            printf("[%i]: Error reopening temp file %s.\n", my_rank, filename);
    }

    return 0;
}


int DumpTimeSerieDatFile(Link* sys, GlobalVars* globals, unsigned int N, unsigned int* save_list, unsigned int save_size, unsigned int my_save_size, const Lookup * const id_to_loc, int* assignments, char* additional_temp, char* additional_out)
{
    unsigned int size = 16;
    char filename[ASYNCH_MAX_PATH_LENGTH], filenamespace[ASYNCH_MAX_PATH_LENGTH], output_filename[ASYNCH_MAX_PATH_LENGTH];
    char data_storage[16];
    FILE *inputfile = NULL, *outputfile = NULL;

    memset(data_storage, 0, 16);

    //Find total size of a line in the temp files
    unsigned int line_size = CalcTotalOutputSize(globals);

    //Create output .dat file
    if (my_rank == 0)
    {
        if (globals->print_par_flag == 1)
        {
            if (!additional_out)
                sprintf(output_filename, "%s", globals->hydros_loc_filename);
            else
                sprintf(output_filename, "%s_%s", globals->hydros_loc_filename, additional_out);
            for (unsigned int i = 0; i < globals->num_global_params; i++)
            {
                sprintf(filenamespace, "_%.4e", globals->global_params[i]);
                strcat(output_filename, filenamespace);
            }
            sprintf(filenamespace, ".dat");
            strcat(output_filename, filenamespace);
        }
        else
            if (!additional_out)
                sprintf(output_filename, "%s.dat", globals->hydros_loc_filename);
            else
                sprintf(output_filename, "%s_%s.dat", globals->hydros_loc_filename, additional_out);
        outputfile = fopen(output_filename, "w");
        if (!outputfile)
        {
            printf("[%i]: Error opening outputfile %s.\n", my_rank, output_filename);
            return 2;
        }

        //Header
        fprintf(outputfile, "%i\n%i\n", save_size, globals->num_outputs);
    }

    //Open temporary files
    if (my_save_size)
    {
        if (!additional_temp)
            sprintf(filename, "%s", globals->temp_filename);
        else
            sprintf(filename, "%s_%s", globals->temp_filename, additional_temp);				//!!!! When is additional_temp useful? Get rid of it? Global params?? !!!!
        inputfile = fopen(filename, "rb");
        if (!inputfile)
        {
            printf("\n[%i]: Error opening inputfile %s.\n", my_rank, filename);
            return 2;
        }
    }

    //Move data to final output
    for (unsigned int i = 0; i < save_size; i++)
    {
        unsigned int loc = find_link_by_idtoloc(save_list[i], id_to_loc, N);
        int proc = assignments[loc];
        Link* current = &sys[loc];

        if (proc == my_rank)
        {
            //Find the link in the temp file inputfile
            unsigned int counter = 0;	//index of link in my_save_list of its process
            unsigned int id, total_spaces;
            fread(&id, sizeof(unsigned int), 1, inputfile);
            fread(&total_spaces, sizeof(unsigned int), 1, inputfile);
            while (id != save_list[i] && !feof(inputfile))
            {
                long jump_size = (long long int) total_spaces * (long long int) line_size;	//Watch overflows!
                fseek(inputfile, jump_size, SEEK_CUR);
                counter++;
                fread(&id, sizeof(unsigned int), 1, inputfile);
                fread(&total_spaces, sizeof(unsigned int), 1, inputfile);
            }
            if (feof(inputfile))
            {
                printf("\n[%i]: Error: could not find id %u in temp file %s.\n", my_rank, save_list[i], filename);
                return 2;
            }

            //Write id and number of steps
            if (my_rank == 0)
                fprintf(outputfile, "\n%u %u\n", save_list[i], current->disk_iterations);
            else
                MPI_Send(&(current->disk_iterations), 1, MPI_UNSIGNED, 0, save_list[i], MPI_COMM_WORLD);

            //Read data in the temp file
            if (my_rank == 0)
            {
                for (unsigned int k = 0; k < current->disk_iterations; k++)
                {
                    for (unsigned int m = 0; m < globals->num_outputs; m++)
                    {
                        const Output * const out = &globals->outputs[m];

                        fread(data_storage, out->size, 1, inputfile);
                        WriteValue(outputfile, out->specifier, data_storage, out->type, " ");
                    }
                    fprintf(outputfile, "\n");
                }
            }
            else
            {
                for (unsigned int k = 0; k < current->disk_iterations; k++)
                {
                    for (unsigned int m = 0; m < globals->num_outputs; m++)
                    {
                        const Output * const out = &globals->outputs[m];

                        fread(data_storage, out->size, 1, inputfile);
                        MPI_Ssend(&data_storage, size, MPI_CHAR, 0, save_list[i], MPI_COMM_WORLD);
                    }
                }
            }

            //Skip over the last unused space. This is done so that the file does not need to be rewound.
            for (unsigned int k = 0; k < total_spaces - current->disk_iterations; k++)
                for (unsigned int j = 0; j < globals->num_outputs; j++)
                {
                    const Output * const out = &globals->outputs[j];

                    fread(data_storage, out->size, 1, inputfile);
                }
        }
        else if (my_rank == 0)
        {
            //Write to file
            MPI_Recv(&(current->disk_iterations), 1, MPI_UNSIGNED, proc, save_list[i], MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            fprintf(outputfile, "\n%u %u\n", save_list[i], current->disk_iterations);

            for (unsigned int k = 0; k < current->disk_iterations; k++)
            {
                for (unsigned int m = 0; m < globals->num_outputs; m++)
                {
                    const Output * const out = &globals->outputs[m];

                    MPI_Recv(data_storage, size, MPI_CHAR, proc, save_list[i], MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    WriteValue(outputfile, out->specifier, data_storage, out->type, " ");
                }
                fprintf(outputfile, "\n");
            }
        }
    }

    //Cleanup
    if (inputfile)	fclose(inputfile);
    if (outputfile)	fclose(outputfile);

    if (my_rank == 0)
        printf("\nResults written to file %s.\n", output_filename);

    return 0;
}


int DumpTimeSerieCsvFile(Link* sys, GlobalVars* globals, unsigned int N, unsigned int* save_list, unsigned int save_size, unsigned int my_save_size, const Lookup * const id_to_loc, int* assignments, char* additional_temp, char* additional_out)
{
    unsigned int size = 16;
    char filename[ASYNCH_MAX_PATH_LENGTH], filenamespace[ASYNCH_MAX_PATH_LENGTH], output_filename[ASYNCH_MAX_PATH_LENGTH];
    char data_storage[16];
    FILE *inputfile = NULL, *outputfile = NULL;

    memset(data_storage, 0, 16);

    //Find total size of a line in the temp files
    unsigned int line_size = CalcTotalOutputSize(globals);

    //Create output file
    if (my_rank == 0)
    {
        if (globals->print_par_flag == 1)
        {
            if (!additional_out)
                sprintf(output_filename, "%s", globals->hydros_loc_filename);
            else
                sprintf(output_filename, "%s_%s", globals->hydros_loc_filename, additional_out);
            for (unsigned int i = 0; i < globals->num_global_params; i++)
            {
                sprintf(filenamespace, "_%.4e", globals->global_params[i]);
                strcat(output_filename, filenamespace);
            }
            sprintf(filenamespace, ".csv");
            strcat(output_filename, filenamespace);
        }
        else
        {
            if (!additional_out)
                sprintf(output_filename, "%s.csv", globals->hydros_loc_filename);
            else
                sprintf(output_filename, "%s_%s.csv", globals->hydros_loc_filename, additional_out);
        }
        outputfile = fopen(output_filename, "w");
        if (!outputfile)
        {
            printf("\n[%i]: Error opening outputfile %s.\n", my_rank, output_filename);
            return 2;
        }
    }

    //Open input files
    if (my_save_size)
    {
        if (!additional_temp)
            sprintf(filename, "%s", globals->temp_filename);
        else
            sprintf(filename, "%s_%s", globals->temp_filename, additional_temp);
        inputfile = fopen(filename, "rb");
        if (!inputfile)
        {
            printf("\n[%i]: Error opening inputfile %s.\n", my_rank, filename);
            return 2;
        }
    }

    //Initializations
    unsigned int *space_counter = (unsigned int*)calloc(save_size, sizeof(unsigned int));
    unsigned int my_max_disk = 0;
    for (unsigned int j = 0; j < save_size; j++)
    {
        unsigned int loc = find_link_by_idtoloc(save_list[j], id_to_loc, N);
        Link *current = &sys[loc];

        if (assignments[loc] == my_rank)
        {
            if (my_max_disk < current->disk_iterations)	my_max_disk = current->disk_iterations;
            if (my_rank != 0)
                MPI_Send(&(current->disk_iterations), 1, MPI_INT, 0, save_list[j], MPI_COMM_WORLD);
        }
        else if (my_rank == 0)
            MPI_Recv(&(current->disk_iterations), 1, MPI_INT, assignments[loc], save_list[j], MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    unsigned int max_disk;
    MPI_Allreduce(&my_max_disk, &max_disk, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
    fpos_t *positions = (fpos_t*)malloc(save_size * sizeof(fpos_t));

    //Find the starting position of each link in each file
    //Note: We assume the ids in the temp files are ordered by save_list
    for (unsigned int j = 0; j < save_size; j++)
    {
        unsigned int loc = find_link_by_idtoloc(save_list[j], id_to_loc, N);	//!!!! Ugh... !!!!

        if (assignments[loc] == my_rank)
        {
            fseek(inputfile, sizeof(unsigned int), SEEK_CUR);	//Skip ID
            unsigned int total_spaces;
            fread(&total_spaces, sizeof(unsigned int), 1, inputfile);	//Grab total spaces for this link
            fgetpos(inputfile, &(positions[j]));		//Here is the place to start for this link
            long jump_size = total_spaces * line_size;	//Watch overflows!
            fseek(inputfile, jump_size, SEEK_CUR);	//Skip to next link
        }
    }

    //Make the .csv header
    if (my_rank == 0)
    {
        for (unsigned int i = 0; i < save_size; i++)
        {
            fprintf(outputfile, "Link %u", save_list[i]);
            for (unsigned int k = 0; k < globals->num_outputs; k++)	fprintf(outputfile, " ,");
        }
        fprintf(outputfile, "\n");

        for (unsigned int i = 0; i < save_size; i++)
        {
            unsigned int loc = find_link_by_idtoloc(save_list[i], id_to_loc, N);	//!!!! Ugh... !!!!
            for (unsigned int k = 0; k < globals->num_outputs; k++)
                fprintf(outputfile, "Output_%u,", k);	//!!!! Use names. What if skipping some? !!!!
        }
        fprintf(outputfile, "\n");
    }

    //Make the .csv body
    if (my_rank == 0)
    {
        for (unsigned int m = 0; m < max_disk; m++)
        {
            for (unsigned int i = 0; i < save_size; i++)
            {
                unsigned int loc = find_link_by_idtoloc(save_list[i], id_to_loc, N);	//!!!! Ugh... !!!!
                Link *current = &sys[loc];

                if (space_counter[i] > current->disk_iterations)	//This link is done, leave blanks
                    for (unsigned int k = 0; k < globals->num_outputs; k++)	fprintf(outputfile, ",");
                else
                {
                    if (assignments[loc] == 0)
                    {
                        fsetpos(inputfile, &(positions[i]));
                        for (unsigned int l = 0; l < globals->num_outputs; l++)
                        {
                            const Output * const out = &globals->outputs[l];
                            fread(data_storage, out->size, 1, inputfile);
                            WriteValue(outputfile, out->specifier, data_storage, out->type, ",");
                        }
                        fgetpos(inputfile, &(positions[i]));
                        (space_counter[i])++;
                    }
                    else
                    {
                        for (unsigned int l = 0; l < globals->num_outputs; l++)
                        {
                            const Output * const out = &globals->outputs[l];
                            MPI_Recv(&data_storage, 16, MPI_CHAR, assignments[loc], save_list[i], MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                            WriteValue(outputfile, out->specifier, data_storage, out->type, ",");
                        }
                        (space_counter[i])++;
                    }
                }
            }
            fprintf(outputfile, "\n");
        }
    }
    else
    {
        for (unsigned int m = 0; m < max_disk; m++)
        {
            for (unsigned int i = 0; i < save_size; i++)
            {
                unsigned int loc = find_link_by_idtoloc(save_list[i], id_to_loc, N);	//!!!! Ugh... !!!!
                Link *current = &sys[loc];

                if (assignments[loc] == my_rank && space_counter[i] <= current->disk_iterations)
                {
                    fsetpos(inputfile, &(positions[i]));
                    for (unsigned int l = 0; l < globals->num_outputs; l++)
                    {
                        fread(data_storage, globals->outputs[l].size, 1, inputfile);
                        MPI_Send(&data_storage, 16, MPI_CHAR, 0, save_list[i], MPI_COMM_WORLD);
                    }
                    fgetpos(inputfile, &(positions[i]));
                    (space_counter[i])++;
                }
            }
        }
    }

    //Clean up
    if (outputfile)	fclose(outputfile);
    if (inputfile)	fclose(inputfile);
    free(positions);
    free(space_counter);

    if (my_rank == 0)
        printf("\nResults written to file %s.\n", output_filename);

    return 0;
}


int DumpTimeSerieH5File(Link* sys, GlobalVars* globals, unsigned int N, unsigned int* save_list, unsigned int save_size, unsigned int my_save_size, const Lookup * const id_to_loc, int* assignments, char* additional_temp, char* additional_out)
{
    unsigned int size = 16;
    char filename[ASYNCH_MAX_PATH_LENGTH], filenamespace[ASYNCH_MAX_PATH_LENGTH], output_filename[ASYNCH_MAX_PATH_LENGTH];
    char *data_storage;
    FILE *inputfile = NULL;
    hid_t file_id;
    hid_t compound_id;
    hid_t packet_file_id;
    
    const hsize_t chunk_size = 512; // Chunk size, in number of table entries per chunk
    const int compression = 5;      // Compression level, a value of 0 through 9.

    //Find total size of a line in the temp files
    unsigned int line_size = CalcTotalOutputSize(globals);

    data_storage = malloc(chunk_size * line_size);

    //Create output file
    if (my_rank == 0)
    {
        if (globals->print_par_flag == 1)
        {
            if (!additional_out)
                sprintf(output_filename, "%s", globals->hydros_loc_filename);
            else
                sprintf(output_filename, "%s_%s", globals->hydros_loc_filename, additional_out);
            for (unsigned int i = 0; i < globals->num_global_params; i++)
            {
                sprintf(filenamespace, "_%.4e", globals->global_params[i]);
                strcat(output_filename, filenamespace);
            }
            sprintf(filenamespace, ".h5");
            strcat(output_filename, filenamespace);
        }
        else
        {
            if (!additional_out)
                sprintf(output_filename, "%s.h5", globals->hydros_loc_filename);
            else
                sprintf(output_filename, "%s_%s.h5", globals->hydros_loc_filename, additional_out);
        }
        
        //Create output file
        file_id = H5Fcreate(output_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (file_id < 0)
        {
            printf("Error: could not open h5 file %s.\n", output_filename);
            return 2;
        }

        //Set attributes
        unsigned short type = globals->model_uid;
        unsigned int issue_time = (unsigned int) globals->begin_time;
        H5LTset_attribute_string(file_id, "/", "version", PACKAGE_VERSION);
        H5LTset_attribute_ushort(file_id, "/", "model", &type, 1);
        H5LTset_attribute_uint(file_id, "/", "issue_time", &issue_time, 1);

        //Create compound type
        compound_id = H5Tcreate(H5T_COMPOUND, line_size);
        size_t offset = 0;
        for (unsigned int i = 0 ; i < globals->num_outputs; i++)
        {
            const Output * const out = &globals->outputs[i];
            herr_t status = H5Tinsert(compound_id, out->name, offset, Get_H5_Type(out->type));
            offset += out->size;
        }

        // Create packet file
        packet_file_id = H5PTcreate_fl(file_id, "outputs", compound_id, chunk_size, compression);
        if (packet_file_id < 0)
        {
            printf("Error: could not initialize h5 file %s.\n", output_filename);
            return 2;
        }
    }

    //Open input files
    if (my_save_size)
    {
        if (!additional_temp)
            sprintf(filename, "%s", globals->temp_filename);
        else
            sprintf(filename, "%s_%s", globals->temp_filename, additional_temp);
        inputfile = fopen(filename, "rb");
        if (!inputfile)
        {
            printf("\n[%i]: Error opening inputfile %s.\n", my_rank, filename);
            return 2;
        }
    }

    //Move data to final output
    for (unsigned int i = 0; i < save_size; i++)
    {
        unsigned int loc = find_link_by_idtoloc(save_list[i], id_to_loc, N);
        int proc = assignments[loc];
        Link* current = &sys[loc];

        if (proc == my_rank)
        {
            //Find the link in the temp file inputfile
            unsigned int counter = 0;	//index of link in my_save_list of its process
            unsigned int id, total_spaces;
            fread(&id, sizeof(unsigned int), 1, inputfile);
            fread(&total_spaces, sizeof(unsigned int), 1, inputfile);
            while (id != save_list[i] && !feof(inputfile))
            {
                long jump_size = (long) total_spaces * (long) line_size;	//Watch overflows!
                fseek(inputfile, jump_size, SEEK_CUR);
                counter++;
                fread(&id, sizeof(unsigned int), 1, inputfile);
                fread(&total_spaces, sizeof(unsigned int), 1, inputfile);
            }
            if (feof(inputfile))
            {
                printf("\n[%i]: Error: could not find id %u in temp file %s.\n", my_rank, save_list[i], filename);
                return 2;
            }

            //Read data in the temp file
            if (my_rank == 0)
            {
                for (size_t k = 0; k < current->disk_iterations; k += chunk_size)
                {
                    size_t reminder = min(chunk_size, current->disk_iterations - k);
                    size_t num_read = fread(data_storage, line_size, reminder, inputfile);

                    herr_t ret = H5PTappend(packet_file_id, num_read, data_storage);
                }
            }
            else
            {
                MPI_Ssend(&(current->disk_iterations), 1, MPI_UNSIGNED, 0, save_list[i], MPI_COMM_WORLD);

                for (hsize_t k = 0; k < current->disk_iterations; k += chunk_size)
                {
                    size_t reminder = min(chunk_size, current->disk_iterations - k);
                    unsigned int num_read = (int) fread(data_storage, line_size, reminder, inputfile);

                    MPI_Ssend(&num_read, 1, MPI_UNSIGNED, 0, save_list[i], MPI_COMM_WORLD);
                    MPI_Ssend(data_storage, num_read * line_size, MPI_CHAR, 0, save_list[i], MPI_COMM_WORLD);
                }
            }

            //Skip over the last unused space. This is done so that the file does not need to be rewound.
            for (unsigned int k = 0; k < total_spaces - current->disk_iterations; k++)
                for (unsigned int j = 0; j < globals->num_outputs; j++)
                    fread(data_storage, globals->outputs[j].size, 1, inputfile);
        }
        else if (my_rank == 0)
        {
            //Write to file
            MPI_Recv(&(current->disk_iterations), 1, MPI_UNSIGNED, proc, save_list[i], MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for (hsize_t k = 0; k < current->disk_iterations; k += chunk_size)
            {
                unsigned int num_read;

                MPI_Recv(&num_read, 1, MPI_UNSIGNED, proc, save_list[i], MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(data_storage, num_read * line_size, MPI_CHAR, proc, save_list[i], MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                H5PTappend(packet_file_id, num_read, data_storage);
            }
        }
    }

    //Cleanup
    if (inputfile)
        fclose(inputfile);
    
    if (my_rank == 0)
    {
        H5PTclose(packet_file_id);
        H5Fclose(file_id);
        H5Tclose(compound_id);
    }

    if (my_rank == 0)
        printf("\nResults written to file %s.\n", output_filename);

    return 0;
}


int DumpTimeSerieNcFile(Link* sys, GlobalVars* globals, unsigned int N, unsigned int* save_list, unsigned int save_size, unsigned int my_save_size, const Lookup * const id_to_loc, int* assignments, char* additional_temp, char* additional_out)
{
    unsigned int size = 16;
    char filename[ASYNCH_MAX_PATH_LENGTH], filenamespace[ASYNCH_MAX_PATH_LENGTH], output_filename[ASYNCH_MAX_PATH_LENGTH];
    char *data_storage;
    FILE *inputfile = NULL;
    hid_t file_id;
    hid_t dataset_id;
    hid_t mem_dataspace_id, file_dataspace_id;
    
    hsize_t start[3];  // Start of hyperslab
    hsize_t count[3];  // Block count

    const hsize_t chunk_size = 512; // Chunk size, in number of table entries per chunk
    const int compression = 5;      // Compression level, a value of 0 through 9.

    //Find total size of a line in the temp files
    unsigned int line_size = CalcTotalOutputSize(globals);

    data_storage = malloc(chunk_size * line_size);

    //Create output file
    if (my_rank == 0)
    {
        if (globals->print_par_flag == 1)
        {
            if (!additional_out)
                sprintf(output_filename, "%s", globals->hydros_loc_filename);
            else
                sprintf(output_filename, "%s_%s", globals->hydros_loc_filename, additional_out);
            for (unsigned int i = 0; i < globals->num_global_params; i++)
            {
                sprintf(filenamespace, "_%.4e", globals->global_params[i]);
                strcat(output_filename, filenamespace);
            }
            sprintf(filenamespace, ".h5");
            strcat(output_filename, filenamespace);
        }
        else
        {
            if (!additional_out)
                sprintf(output_filename, "%s.h5", globals->hydros_loc_filename);
            else
                sprintf(output_filename, "%s_%s.h5", globals->hydros_loc_filename, additional_out);
        }

        //Create output file
        file_id = H5Fcreate(output_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (file_id < 0)
        {
            printf("Error: could not open h5 file %s.\n", output_filename);
            return 2;
        }

        //Set attributes
        unsigned short type = globals->model_uid;
        unsigned int issue_time = (unsigned int)globals->begin_time;
        H5LTset_attribute_string(file_id, "/", "version", PACKAGE_VERSION);
        H5LTset_attribute_ushort(file_id, "/", "model", &type, 1);
        H5LTset_attribute_uint(file_id, "/", "issue_time", &issue_time, 1);

        //Set the link_id axis
        hsize_t num_link = save_size;
        H5LTmake_dataset_int(file_id, "link_id", 1, &num_link, save_list);

        //Make it a dimension scale
        hid_t link_id_ds = H5Dopen(file_id, "link_id", H5P_DEFAULT);
        H5DSset_scale(link_id_ds, "link_id");

        //Set the time axis
        hsize_t num_timestep = (int) (globals->maxtime / globals->print_time) + 1;
        int *buffer = malloc(num_timestep * sizeof(int));
        for (unsigned int i = 0; i < num_timestep; i++)
        {
            buffer[i] = (int) (globals->begin_time + (int)(i * globals->print_time * 60));
        }
        H5LTmake_dataset_int(file_id, "time", 1, &num_timestep, buffer);
        free(buffer);
        
        //Make it a dimension scale
        hid_t time_ds = H5Dopen(file_id, "time", H5P_DEFAULT);
        H5DSset_scale(time_ds, "time");

        //Check if all the outputs are of type float (limitation of the current implementation)
        for (unsigned int i = 0; i < globals->num_outputs; i++)
        {
            const Output * const out = &globals->outputs[i];
            if (out->type != ASYNCH_FLOAT)
            {
                printf("Error: All output types should be float (got %i for %s).\n", out->type, out->name);
                MPI_Abort(MPI_COMM_WORLD, 0);
            }
        }

        //Set the output axis
        hsize_t num_outputs = globals->num_outputs;
        buffer = malloc(globals->num_outputs * sizeof(int));
        for (unsigned int i = 0; i < num_outputs; i++)
        {
            buffer[i] = i;
        }
        H5LTmake_dataset_int(file_id, "output", 1, &num_outputs, buffer);
        free(buffer);

        //Make it a dimension scale
        hid_t output_ds = H5Dopen(file_id, "output", H5P_DEFAULT);
        H5DSset_scale(output_ds, "output");

        //Create the memory data space
        hsize_t dims[2];  // Block count
        dims[0] = num_timestep;
        dims[1] = num_outputs;
        mem_dataspace_id = H5Screate_simple(2, dims, NULL);

        // Describe the size of the array and create the data space for fixed size dataset.
        hsize_t dimsf[3];
        dimsf[0] = save_size;
        dimsf[1] = num_timestep;
        dimsf[2] = globals->num_outputs;
        file_dataspace_id = H5Screate_simple(3, dimsf, NULL);

        // Create a new dataset within the file using defined dataspace and
        // datatype and default dataset creation properties.
        dataset_id = H5Dcreate(file_id, "outputs", H5T_NATIVE_FLOAT, file_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (dataset_id < 0)
        {
            printf("Error: could not initialize h5 file %s.\n", output_filename);
            return 2;
        }

        H5DSattach_scale(dataset_id, link_id_ds, 0);
        H5DSattach_scale(dataset_id, time_ds, 1);
        H5DSattach_scale(dataset_id, output_ds, 2);

        H5Dclose(link_id_ds);
        H5Dclose(time_ds);
        H5Dclose(output_ds);
    }

    //Open input files
    if (my_save_size)
    {
        if (!additional_temp)
            sprintf(filename, "%s", globals->temp_filename);
        else
            sprintf(filename, "%s_%s", globals->temp_filename, additional_temp);
        inputfile = fopen(filename, "rb");
        if (!inputfile)
        {
            printf("\n[%i]: Error opening inputfile %s.\n", my_rank, filename);
            return 2;
        }
    }

    //Move data to final output
    for (unsigned int i = 0; i < save_size; i++)
    {
        unsigned int loc = find_link_by_idtoloc(save_list[i], id_to_loc, N);
        int proc = assignments[loc];
        Link* current = &sys[loc];

        if (proc == my_rank)
        {
            //Find the link in the temp file inputfile
            unsigned int counter = 0;	//index of link in my_save_list of its process
            unsigned int id, total_spaces;
            fread(&id, sizeof(unsigned int), 1, inputfile);
            fread(&total_spaces, sizeof(unsigned int), 1, inputfile);
            while (id != save_list[i] && !feof(inputfile))
            {
                long jump_size = (long)total_spaces * (long)line_size;	//Watch overflows!
                fseek(inputfile, jump_size, SEEK_CUR);
                counter++;
                fread(&id, sizeof(unsigned int), 1, inputfile);
                fread(&total_spaces, sizeof(unsigned int), 1, inputfile);
            }
            if (feof(inputfile))
            {
                printf("\n[%i]: Error: could not find id %u in temp file %s.\n", my_rank, save_list[i], filename);
                return 2;
            }

            //Read data in the temp file
            if (my_rank == 0)
            {
                for (size_t k = 0; k < current->disk_iterations; k += chunk_size)
                {
                    size_t reminder = min(chunk_size, current->disk_iterations - k);
                    size_t num_read = fread(data_storage, line_size, reminder, inputfile);
                    assert(num_read <= (globals->maxtime / globals->print_time) + 1);

                    start[0] = i;
                    start[1] = k;
                    start[2] = 0;
                    //stride[0] = 1;
                    //stride[1] = 1;
                    count[0] = 1;
                    count[1] = num_read;
                    count[2] = globals->num_outputs;
                    //block[0] = 1;
                    //block[1] = 1;
                    H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL);

                    herr_t ret = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, mem_dataspace_id, file_dataspace_id, H5P_DEFAULT, data_storage);
                }
            }
            else
            {
                MPI_Ssend(&(current->disk_iterations), 1, MPI_UNSIGNED, 0, save_list[i], MPI_COMM_WORLD);

                for (hsize_t k = 0; k < current->disk_iterations; k += chunk_size)
                {
                    size_t reminder = min(chunk_size, current->disk_iterations - k);
                    unsigned int num_read = (unsigned int)fread(data_storage, line_size, reminder, inputfile);

                    MPI_Ssend(&num_read, 1, MPI_UNSIGNED, 0, save_list[i], MPI_COMM_WORLD);
                    MPI_Ssend(data_storage, num_read * line_size, MPI_CHAR, 0, save_list[i], MPI_COMM_WORLD);
                }
            }

            //Skip over the last unused space. This is done so that the file does not need to be rewound.
            for (unsigned int k = 0; k < total_spaces - current->disk_iterations; k++)
                for (unsigned int j = 0; j < globals->num_outputs; j++)
                    fread(data_storage, globals->outputs[j].size, 1, inputfile);
        }
        else if (my_rank == 0)
        {
            //Write to file
            MPI_Recv(&(current->disk_iterations), 1, MPI_UNSIGNED, proc, save_list[i], MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for (hsize_t k = 0; k < current->disk_iterations; k += chunk_size)
            {
                unsigned int num_read;

                MPI_Recv(&num_read, 1, MPI_UNSIGNED, proc, save_list[i], MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(data_storage, num_read * line_size, MPI_CHAR, proc, save_list[i], MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                start[0] = i;
                start[1] = k;
                start[2] = 0;
                //stride[0] = 1;
                //stride[1] = 1;
                count[0] = 1;
                count[1] = num_read;
                count[2] = globals->num_outputs;
                //block[0] = 1;
                //block[1] = 1;
                H5Sselect_hyperslab(file_dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL);

                herr_t ret = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, mem_dataspace_id, file_dataspace_id, H5P_DEFAULT, data_storage);
            }
        }
    }

    //Cleanup
    if (inputfile)
        fclose(inputfile);

    if (my_rank == 0)
    {
        H5Sclose(file_dataspace_id);
        H5Sclose(mem_dataspace_id);
        H5Dclose(dataset_id);
        H5Fclose(file_id);
    }

    if (my_rank == 0)
        printf("\nResults written to file %s.\n", output_filename);

    return 0;
}


#if defined(HAVE_POSTGRESQL)

void PrepareDatabaseTable(GlobalVars* globals, ConnData* conninfo)
{
    unsigned int num_cols, i;
    PGresult* res;
    char* query;

    if (conninfo == NULL)	return;

    if (my_rank == 0)
    {
        ConnectPGDB(conninfo);
        query = conninfo->query;

        //Create table
        if (conninfo->num_queries > 1)
        {
            sprintf(query, conninfo->queries[1], globals->hydro_table);
            res = PQexec(conninfo->conn, query);
            CheckResError(res, "deleting hydrograph table");
            PQclear(res);
        }
        sprintf(query, conninfo->queries[0], globals->hydro_table);
        res = PQexec(conninfo->conn, query);
        CheckResError(res, "creating hydrograph table");
        PQclear(res);

        //Check if the table name contains a schema
        char* sep = strchr(globals->hydro_table, '.');
        if (sep != NULL)
        {
            *sep = '\0';
            const char* schema = globals->hydro_table;
            const char* table = sep + 1;
            sprintf(query, "SELECT data_type FROM information_schema.columns WHERE table_schema='%s' AND table_name='%s' ORDER BY ordinal_position;", schema, table);
            *sep = '.';
        }
        else
            sprintf(query, "SELECT data_type FROM information_schema.columns WHERE table_name='%s' ORDER BY ordinal_position;", globals->hydro_table);

        //Make sure table is consistent with outputs
        res = PQexec(conninfo->conn, query);
        CheckResError(res, "obtaining data types from hydrograph table");

        num_cols = PQntuples(res);
        if (num_cols > globals->num_outputs)
        {
            printf("[%i]: Error: need more outputs in .gbl file.\n", my_rank);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else if (num_cols < globals->num_outputs)
        {
            printf("[%i]: Error: need more columns in table %s. Got %i, expected %i.\n", my_rank, globals->hydro_table, num_cols, globals->num_outputs);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        for (i = 0; i < num_cols; i++)
        {
            if (globals->outputs[i].type == ASYNCH_BAD_TYPE)
            {
                printf("[%i]: Error: output %i must be set to prepare database tables.\n", my_rank, i);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            else if (strcmp(PQgetvalue(res, i, 0), "double precision") == 0)
            {
                if (globals->outputs[i].type != ASYNCH_DOUBLE)
                {
                    printf("[%i]: Error: Output %i is of type double precision in output table, but should not be. %i\n", my_rank, i, globals->outputs[i].type);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
            }
            else if (strcmp(PQgetvalue(res, i, 0), "integer") == 0)
            {
                if (globals->outputs[i].type != ASYNCH_INT)
                {
                    printf("[%i]: Error: Output %i is of type integer in output table, but should not be. %i\n", my_rank, i, globals->outputs[i].type);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
            }
            else if (strcmp(PQgetvalue(res, i, 0), "real") == 0)
            {
                if (globals->outputs[i].type != ASYNCH_FLOAT)
                {
                    printf("[%i]: Error: Output %i is of type single precision in output table, but should not be. %i\n", my_rank, i, globals->outputs[i].type);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
            }
            else if (strcmp(PQgetvalue(res, i, 0), "short integer") == 0)
            {
                if (globals->outputs[i].type != ASYNCH_SHORT)
                {
                    printf("[%i]: Error: Output %i is of type short integer in output table, but should not be. %i\n", my_rank, i, globals->outputs[i].type);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
            }
            else if (strcmp(PQgetvalue(res, i, 0), "character") == 0)
            {
                if (globals->outputs[i].type != ASYNCH_CHAR)
                {
                    printf("[%i]: Error: Output %i is of type character in output table, but should not be. %i\n", my_rank, i, globals->outputs[i].type);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
            }
            else
            {
                printf("[%i]: Error: Bad datatype for output %i while preparing database tables. (%s)\n", my_rank, i, PQgetvalue(res, i, 0));
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }

        PQclear(res);
        DisconnectPGDB(conninfo);
    }
}



//Assumes temp files.
//Return value = 0 means everything is good.
//1 means a database related error.
//2 means a file system related error.
int DumpTimeSerieDB(Link* sys, GlobalVars* globals, unsigned int N, unsigned int* save_list, unsigned int save_size, unsigned int my_save_size, const Lookup * const id_to_loc, int* assignments, char* additional_temp, char* additional_out, ConnData* conninfo, FILE** my_tempfile)
{
    int my_result = 0, result = 0, return_val = 0;
    //int i, k, nbytes, my_result = 0, result = 0, return_val = 0;
    //unsigned int j, loc, id, total_spaces;
    char filename[ASYNCH_MAX_PATH_LENGTH], temptablename[ASYNCH_MAX_PATH_LENGTH];
    //char* submission;
    char data_storage[16];
    //unsigned int dim = globals->num_print;
    FILE *inputfile;
    //PGresult *res;
    //Link* current;
    //long long int jump_size;

    //Close the temp file, if open
    if (my_tempfile && *my_tempfile)
    {
        fclose(*my_tempfile);
        *my_tempfile = NULL;
    }

    //Assumes 24 chars max per double, plus delims, plus link id. So %.8e should work.
    char* submission = (char*)malloc((globals->num_outputs * 24 + globals->num_outputs) * sizeof(char));

    //Find total size of a line in the temp files
    unsigned int line_size = CalcTotalOutputSize(globals);

    //Open temporary file
    if (!additional_temp)
        sprintf(filename, "%s", globals->temp_filename);
    else
        sprintf(filename, "%s_%s", globals->temp_filename, additional_temp);
    inputfile = fopen(filename, "rb");

    if (my_save_size > 0)	//!!!! This wasn't here before. But I think it should be... !!!!
    {
        if (!inputfile)
        {
            printf("\n[%i]: Error opening inputfile %s.\n", my_rank, filename);
            my_result = 1;
        }
        else
            my_result = 0;
    }

    //Check if an error occurred
    MPI_Allreduce(&my_result, &result, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
    if (result)
    {
        return_val = 2;
        goto finished;	//Don't bother uploading anymore
    }

    if (my_rank == 0)
    {
        //sprintf(temptablename,"tmp_%s",globals->hydro_table);
        sprintf(temptablename, "%s_tmp", globals->hydro_table);
        return_val = 0;
        ConnectPGDB(conninfo);

        //Delete temporary table, if it exists
        sprintf(conninfo->query, "DROP TABLE IF EXISTS %s;", temptablename);
        PGresult *res = PQexec(conninfo->conn, conninfo->query);
        result = CheckResError(res, "dropping temporary output table");
        PQclear(res);

        //Create temporary table
        sprintf(conninfo->query, conninfo->queries[0], temptablename);
        res = PQexec(conninfo->conn, conninfo->query);
        result = CheckResError(res, "creating temporary output table");
        PQclear(res);

        //Tell database to prepare for copying
        sprintf(conninfo->query, "COPY %s FROM STDIN WITH DELIMITER ',';", temptablename);
        res = PQexec(conninfo->conn, conninfo->query);
        result = CheckResState(res, PGRES_COPY_IN);
        PQclear(res);
    }

    //Check if an error occurred
    MPI_Bcast(&result, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (result)
    {
        return_val = 1;
        goto finished;
    }

    //Upload data
    for (unsigned int i = 0; i < save_size; i++)
    {
        unsigned int loc = find_link_by_idtoloc(save_list[i], id_to_loc, N);
        Link *current = &sys[loc];

        if (my_rank == 0)
        {
            //Find the link in the temp file inputfile
            if (assignments[loc] == my_rank)
            {
                unsigned int id, total_spaces;
                fread(&id, sizeof(unsigned int), 1, inputfile);
                fread(&total_spaces, sizeof(unsigned int), 1, inputfile);
                while (id != save_list[i] && !feof(inputfile))
                {
                    long jump_size = (long)total_spaces * (long)line_size;	//Watch overflows!
                    fseek(inputfile, jump_size, SEEK_CUR);
                    fread(&id, sizeof(unsigned int), 1, inputfile);
                    fread(&total_spaces, sizeof(unsigned int), 1, inputfile);
                }
                if (feof(inputfile))
                {
                    printf("\n[%i]: Error: could not find id %u in temp file %s.\n", my_rank, save_list[i], filename);
                    return_val = 2;
                    goto finished;	//Don't bother uploading anymore
                }

                //Now read in the data in the temp file, and submit them to the database
                for (unsigned int k = 0; k < current->disk_iterations; k++)
                {
                    unsigned int nbytes = 0;
                    unsigned int j;
                    for (j = 0; j < globals->num_outputs - 1; j++)
                    {
                        const Output * const out = &globals->outputs[j];
                        fread(data_storage, out->size, 1, inputfile);
                        nbytes += CatBinaryToString(&(submission[nbytes]), out->specifier, data_storage, out->type, ",");
                    }
                    const Output * const out = &globals->outputs[j];
                    fread(data_storage, out->size, 1, inputfile);
                    nbytes += CatBinaryToString(&(submission[nbytes]), out->specifier, data_storage, out->type, "\n");
                    result = PQputCopyData(conninfo->conn, submission, nbytes);
                }

                //Skip over the last unused space. This is done so that the file does not need to be rewound.
                for (unsigned int k = 0; k < total_spaces - current->disk_iterations; k++)
                    for (unsigned int j = 0; j < globals->num_outputs; j++)
                        fread(data_storage, globals->outputs[j].size, 1, inputfile);
            }
            else
            {
                //Get disk_iterations
                MPI_Recv(&(current->disk_iterations), 1, MPI_UNSIGNED, assignments[loc], loc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                //Now read in the data in the temp file, and submit them to the database
                for (unsigned int k = 0; k < current->disk_iterations; k++)
                {
                    unsigned int nbytes;
                    MPI_Recv(&nbytes, 1, MPI_UNSIGNED, assignments[loc], loc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(submission, nbytes, MPI_CHAR, assignments[loc], loc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    result = PQputCopyData(conninfo->conn, submission, nbytes);
                }
            }

            //Check if an error occurred.
            MPI_Bcast(&result, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if (result != 1)
            {
                printf("Error: copy returned result %i.\n", result);
                return_val = 1;
                goto finished;	//Don't bother uploading anymore
            }
        }
        else
        {
            //Find the link in the temp file inputfile
            if (assignments[loc] == my_rank)
            {
                unsigned int id, total_spaces;
                fread(&id, sizeof(unsigned int), 1, inputfile);
                fread(&total_spaces, sizeof(unsigned int), 1, inputfile);

                while (id != save_list[i] && !feof(inputfile))
                {
                    long jump_size = total_spaces * line_size;	//Watch overflows!
                    fseek(inputfile, jump_size, SEEK_CUR);
                    fread(&id, sizeof(unsigned int), 1, inputfile);
                    fread(&total_spaces, sizeof(unsigned int), 1, inputfile);
                }
                if (feof(inputfile))
                {
                    printf("\n[%i]: Error: could not find id %u in temp file %s.\n", my_rank, save_list[i], filename);
                    return_val = 2;
                    goto finished;	//Don't bother uploading anymore
                }

                //Send disk_iterations
                MPI_Send(&(current->disk_iterations), 1, MPI_UNSIGNED, 0, loc, MPI_COMM_WORLD);

                //Now read in the data in the temp file, and submit them to the database
                for (unsigned int k = 0; k < current->disk_iterations; k++)
                {
                    unsigned int nbytes = 0;
                    unsigned int j;
                    for (j = 0; j < globals->num_outputs - 1; j++)
                    {
                        const Output * const out = &globals->outputs[j];
                        fread(data_storage, out->size, 1, inputfile);
                        nbytes += CatBinaryToString(&(submission[nbytes]), out->specifier, data_storage, out->type, ",");
                    }
                    const Output * const out = &globals->outputs[j];
                    fread(data_storage, out->size, 1, inputfile);
                    nbytes += CatBinaryToString(&(submission[nbytes]), out->specifier, data_storage, out->type, "\n");

                    MPI_Send(&nbytes, 1, MPI_UNSIGNED, 0, loc, MPI_COMM_WORLD);
                    MPI_Send(submission, nbytes, MPI_CHAR, 0, loc, MPI_COMM_WORLD);
                }

                //Skip over the last unused space. This is done so that the file does not need to be rewound.
                for (unsigned int k = 0; k < total_spaces - current->disk_iterations; k++)
                    for (unsigned int j = 0; j < globals->num_outputs; j++)
                        fread(data_storage, globals->outputs[j].size, 1, inputfile);
            }

            //Check if an error occurred
            MPI_Bcast(&result, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if (result != 1)
            {
                return_val = 1;
                goto finished;	//Don't bother uploading anymore
            }
        }
    }

finished:

    if (my_rank == 0)
    {
        //Tell database everything is uploaded
        int res = PQputCopyEnd(conninfo->conn, NULL);
        if (res != 1)
        {
            printf("Returned %i while closing copy to hydrograph database.\n", res);
            return_val = 1;
        }

        //If the temporary table was loaded successfully, inserted it into the main table
        if (!return_val)
        {
            sprintf(conninfo->query, "INSERT INTO %s (SELECT * FROM %s);", globals->hydro_table, temptablename);
            PGresult *res = PQexec(conninfo->conn, conninfo->query);
            CheckResError(res, "inserting temporary output table to final table");
            PQclear(res);

            //Delete temporary table, if it exists
            sprintf(conninfo->query, "DROP TABLE IF EXISTS %s;", temptablename);
            res = PQexec(conninfo->conn, conninfo->query);
            CheckResError(res, "dropping temporary output table");
            PQclear(res);
        }

        //Clean up
        DisconnectPGDB(conninfo);

        if (!return_val)	printf("\nResults written to table %s.\n", globals->hydro_table);
    }

    MPI_Bcast(&return_val, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //Reopen the tempfile
    if (my_save_size > 0)
    {
        if (additional_temp != NULL)
            sprintf(filename, "%s_%s", globals->temp_filename, additional_temp);
        else
            sprintf(filename, "%s", globals->temp_filename);
        *my_tempfile = fopen(filename, "r+b");
        if (!*my_tempfile)
            printf("[%i]: Error reopening temp file %s.\n", my_rank, filename);
    }

    if (inputfile)	fclose(inputfile);
    free(submission);

    return return_val;
}

#endif //HAVE_POSTGRESQL

//Creates a peakflow file.
//Make sure GlobalVars has been initialized first.
int PreparePeakFlowFiles(GlobalVars* globals, unsigned int peaksave_size)
{
    unsigned int i;

    if (!globals)
    {
        printf("[%i]: Error: A .gbl file must be read before preparing peakflow files.\n", my_rank);
        return 1;
    }

    if (peaksave_size)
    {
        //Put together the output filename string
        char outputfilename[ASYNCH_MAX_PATH_LENGTH], filename[ASYNCH_MAX_PATH_LENGTH];
        sprintf(outputfilename, "%s", globals->peaks_loc_filename);
        if (globals->print_par_flag == 1)
        {
            for (i = 0; i < globals->num_global_params; i++)
            {
                sprintf(filename, "_%.4e", globals->global_params[i]);
                strcat(outputfilename, filename);
            }
        }

        //Setup peakflow file
        globals->peakfilename = (char*)malloc(globals->string_size * sizeof(char));

        //Setup the output peakflows files
        sprintf(globals->peakfilename, "%s.pea", outputfilename);
    }

    return 0;
}


//Writes to disk the current peakflow information.
int DumpPeakFlowText(Link* sys, GlobalVars* globals, unsigned int N, int* assignments, unsigned int* peaksave_list, unsigned int peaksave_size, const Lookup * const id_to_loc, ConnData* conninfo)
{
    FILE* peakfile = NULL;
    char buffer[256];
    double conversion = (globals->convertarea_flag) ? 1e-6 : 1.0;

    if (peaksave_size)
    {
        if (!(globals->peakfilename))
        {
            if (my_rank == 0)	printf("Error: Cannot write peakflow data to disk. Peakflow file not prepared.\n");
            return 1;
        }

        //Create file
        unsigned int error = 0;
        if (my_rank == 0)
        {
            peakfile = fopen(globals->peakfilename, "w");
            if (!peakfile)
            {
                error = 1;
                printf("Error: Cannot open peakflow file %s.\n", globals->peakfilename);
            }
            else
                fprintf(peakfile, "%i\n%i\n\n", peaksave_size, globals->model_uid);
        }
        
        MPI_Bcast(&error, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        if (error)
            return 1;

        //Write data to file
        if (my_rank == 0)
        {
            for (unsigned int i = 0; i < peaksave_size; i++)
            {
                unsigned int loc = find_link_by_idtoloc(peaksave_list[i], id_to_loc, N);
                Link *current = &sys[loc];

                if (assignments[loc] == my_rank)
                    globals->peakflow_output(current->ID, current->peak_time, current->peak_value, current->params, globals->global_params, conversion, globals->area_idx, current->peakoutput_user, buffer);
                else
                {
                    unsigned int length;
                    MPI_Recv(&length, 1, MPI_UNSIGNED, assignments[loc], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(buffer, length + 1, MPI_CHAR, assignments[loc], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }

                fprintf(peakfile, "%s", buffer);
            }

            fclose(peakfile);
            printf("Peakflows written to file %s.\n", globals->peakfilename);
        }
        else
        {
            for (unsigned int i = 0; i < peaksave_size; i++)
            {
                unsigned int loc = find_link_by_idtoloc(peaksave_list[i], id_to_loc, N);
                Link *current = &sys[loc];

                if (assignments[loc] == my_rank)
                {
                    //TODO
                    globals->peakflow_output(current->ID, current->peak_time, current->peak_value, current->params, globals->global_params, conversion, globals->area_idx, current->peakoutput_user, buffer);
                    unsigned int length = (unsigned int)strlen(buffer);
                    MPI_Send(&length, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD);
                    MPI_Send(buffer, length + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
                }
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    return 0;
}

#if defined(HAVE_POSTGRESQL)

//Uploads the current peakflow information to a database.
int DumpPeakFlowDB(Link* sys, GlobalVars* globals, unsigned int N, int* assignments, unsigned int* peaksave_list, unsigned int peaksave_size, const Lookup * const id_to_loc, ConnData* conninfo)
{
    unsigned int i, loc, length, result, return_val = 0, error = 0;
    char temptablename[ASYNCH_MAX_QUERY_LENGTH], buffer[256];
    Link* current;
    PGresult *res;
    double conversion = (globals->convertarea_flag) ? 1e-6 : 1.0;

    if (conninfo->num_queries < 1)
        printf("[%i]: Warning: I think you need a file to create a peakflow table...\n", my_rank);

    if (peaksave_size)
    {
        if (my_rank == 0)
        {
            //Prepare temporary table name
            //sprintf(temptablename,"tmp_%s",globals->peak_table);
            sprintf(temptablename, "%s_tmp", globals->peak_table);

            error |= ConnectPGDB(conninfo);

            //Delete and create table. This should NOT be done in an init routine, in case there is an error.
            if (!error && conninfo->num_queries > 1)
            {
                sprintf(conninfo->query, conninfo->queries[1], globals->peak_table);
                res = PQexec(conninfo->conn, conninfo->query);
                error |= CheckResError(res, "deleting peakflow table");
                PQclear(res);
            }
            //if(conninfo->num_queries > 0)
            if (!error)
            {
                //Delete temporary table
                sprintf(conninfo->query, "DROP TABLE IF EXISTS %s;", temptablename);
                res = PQexec(conninfo->conn, conninfo->query);
                error = CheckResError(res, "dropping temporary peakflow table");
                PQclear(res);

                //Create temporary table
                sprintf(conninfo->query, conninfo->queries[0], temptablename);
                res = PQexec(conninfo->conn, conninfo->query);
                error = CheckResError(res, "creating temporary peakflow table");
                PQclear(res);

                //Create final peakflow table
                sprintf(conninfo->query, conninfo->queries[0], globals->peak_table);
                res = PQexec(conninfo->conn, conninfo->query);
                //error |= CheckResError(res,"creating peakflow table");
                CheckResError(res, "creating peakflow table");
                PQclear(res);
            }

            if (error)	DisconnectPGDB(conninfo);
        }

        //Return if proc 0 encountered an error
        MPI_Bcast(&error, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        if (error)	return 1;

        //Upload data
        if (my_rank == 0)
        {
            //Tell database to prepare for copying
            sprintf(conninfo->query, "COPY %s FROM STDIN WITH DELIMITER ' ';", temptablename);
            res = PQexec(conninfo->conn, conninfo->query);
            error = CheckResState(res, PGRES_COPY_IN);
            PQclear(res);

            for (i = 0; i < peaksave_size; i++)
            {
                loc = find_link_by_idtoloc(peaksave_list[i], id_to_loc, N);
                current = &sys[loc];

                if (assignments[loc] == my_rank)
                {
                    //TODO
                    globals->peakflow_output(current->ID, current->peak_time, current->peak_value, current->params, globals->global_params, conversion, globals->area_idx, current->peakoutput_user, buffer);
                    length = (unsigned int)strlen(buffer);
                }
                else
                {
                    MPI_Recv(&length, 1, MPI_UNSIGNED, assignments[loc], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(buffer, length + 1, MPI_CHAR, assignments[loc], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }

                result = PQputCopyData(conninfo->conn, buffer, length);
                if (result != 1)
                {
                    printf("[%i]: Returned %i while copying to peakflow database.\n", my_rank, result);
                    return_val = 1;
                    break;
                }
            }

            //Finish copy
            result = PQputCopyEnd(conninfo->conn, NULL);
            if (result != 1)
            {
                printf("[%i]: Returned %i while closing copy to peakflow database.\n", my_rank, result);
                return_val = 1;
            }

            //Insert temporary table, if no errors have occurred
            if (!return_val)
            {
                sprintf(conninfo->query, "INSERT INTO %s (SELECT * FROM %s);", globals->peak_table, temptablename);
                res = PQexec(conninfo->conn, conninfo->query);
                if (CheckResError(res, "inserting peakflow temporary table into permanent table"))
                    return_val = 1;
                PQclear(res);
            }

            //Delete temporary table
            sprintf(conninfo->query, "DROP TABLE IF EXISTS %s;", temptablename);
            res = PQexec(conninfo->conn, conninfo->query);
            if (CheckResError(res, "dropping temporary peakflow table"))
                return_val = 1;
            PQclear(res);

            DisconnectPGDB(conninfo);
        }
        else
        {
            for (i = 0; i < peaksave_size; i++)
            {
                loc = find_link_by_idtoloc(peaksave_list[i], id_to_loc, N);
                current = &sys[loc];

                if (assignments[loc] == my_rank)
                {
                    //TODO
                    globals->peakflow_output(current->ID, current->peak_time, current->peak_value, current->params, globals->global_params, conversion, globals->area_idx, current->peakoutput_user, buffer);
                    length = (unsigned int)strlen(buffer);
                    MPI_Send(&length, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD);
                    MPI_Send(buffer, length + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
                }
            }
        }

        //Make sure everyone knows if an error occured
        MPI_Bcast(&return_val, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

        if (my_rank == 0 && !return_val)
            printf("Peakflows written to table %s.\n", globals->peak_table);
    }

    return return_val;
}

#endif //HAVE_POSTGRESQL


//This is the same as DataDump, but allows for a generic name
int DumpStateText(Link* sys, unsigned int N, int* assignments, GlobalVars* globals, char* preface, ConnData* conninfo)
{
    unsigned int i, j;
    FILE* output;
    //char filename[256];
    //unsigned int dim = globals->dim;
    double buffer[ASYNCH_MAX_DIM];

    if (my_rank == 0)	//Creating the file
    {
        //sprintf(filename,"%s/%s.rec",globals->dump_location,name);
        //if(name)	sprintf(filename,"%s%s.rec",globals->dump_loc_filename,name);
        //else		sprintf(filename,"%s.rec",globals->dump_loc_filename);
        //output = fopen(filename,"w");
        output = fopen(globals->dump_loc_filename, "w");
        if (output == NULL)
        {
            printf("[%i]: Error opening file %s.\n", my_rank, globals->dump_loc_filename);
            i = 1;
        }
        else	i = 0;
        MPI_Bcast(&i, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        if (i)	return 1;

        fprintf(output, "%hu\n%u\n0.0\n\n", globals->model_uid, N);

        for (i = 0; i < N; i++)
        {
            assert(sys[i].dim <= ASYNCH_MAX_DIM);
            if (assignments[i] != 0)
                MPI_Recv(buffer, sys[i].dim, MPI_DOUBLE, assignments[i], i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            else
                memcpy(buffer, sys[i].my->list.tail->y_approx, sys[i].dim * sizeof(double));

            fprintf(output, "%u\n", sys[i].ID);
            for (j = 0; j < sys[i].dim; j++)	fprintf(output, "%.6e ", buffer[j]);
            fprintf(output, "\n");
        }

        fclose(output);
    }
    else			//Sending data to proc 0
    {
        //Check for error
        MPI_Bcast(&i, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        if (i)	return 1;

        for (i = 0; i < N; i++)
        {
            if (assignments[i] == my_rank)
            {
                memcpy(buffer, sys[i].my->list.tail->y_approx, sys[i].dim * sizeof(double));
                MPI_Send(buffer, sys[i].dim, MPI_DOUBLE, 0, i, MPI_COMM_WORLD);
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    return 0;
}

int DumpStateH5(Link* sys, unsigned int N, int* assignments, GlobalVars* globals, char* suffix, ConnData* conninfo)
{
    unsigned int i;
    unsigned int res = 0;

    const hsize_t chunk_size = 512;   // Chunk size, in number of table entries per chunk
    const int compression = 5;      // Compression level, a value of 0 through 9.

    if (my_rank == 0)
    {
        // Creating the file
        unsigned int unix_time = 0;
        char dump_loc_filename[ASYNCH_MAX_PATH_LENGTH];
        if (globals->dump_loc_flag == 4)
            {
                char timestamp[ASYNCH_MAX_TIMESTAMP_LENGTH + 1];
            unix_time = (unsigned int)(globals->begin_time + globals->t * 60);
            sprintf(timestamp, "%u", unix_time);
            snprintf(dump_loc_filename, ASYNCH_MAX_PATH_LENGTH, "%s_%s.h5", globals->dump_loc_filename, timestamp);
        }
        else
        {
            unix_time = (unsigned int)globals->end_time;
            snprintf(dump_loc_filename, ASYNCH_MAX_PATH_LENGTH, "%s", globals->dump_loc_filename);
        }


        hid_t file_id = H5Fcreate(dump_loc_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (file_id < 0)
        {
            printf("Error: could not open h5 file %s.\n", dump_loc_filename);
            res = 1;
            MPI_Bcast(&res, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
            return res;
        }
        else
            MPI_Bcast(&res, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

        //Set attributes
        unsigned short type = globals->model_uid;
        H5LTset_attribute_string(file_id, "/", "version", PACKAGE_VERSION);
        H5LTset_attribute_ushort(file_id, "/", "model", &type, 1);
        H5LTset_attribute_uint(file_id, "/", "unix_time", &unix_time, 1);

        // Find the first link that belongs to this process
        i = 0;
        while (assignments[i] != my_rank)
            i++;

        //Assume that every links have the same dimension
        unsigned int dim = sys[i].dim;

        size_t line_size = sizeof(unsigned int) + dim * sizeof(double);

        //Create compound type
        hid_t compound_id = H5Tcreate(H5T_COMPOUND, line_size);
        herr_t status = H5Tinsert(compound_id, "link_id", 0, H5T_NATIVE_UINT);

        size_t offset = sizeof(unsigned int);
        for (unsigned int i = 0; i < dim; i++)
        {
            char name[10];
            sprintf(name, "state_%d", i);

            herr_t status = H5Tinsert(compound_id, name, offset, H5T_NATIVE_DOUBLE);
            offset += sizeof(double);
        }

        // Create packet file
        hid_t packet_file_id = H5PTcreate_fl(file_id, "snapshot", compound_id, chunk_size, compression);
        if (packet_file_id < 0)
        {
            printf("Error: could not initialize h5 packet file %s.\n", dump_loc_filename);
            return 2;
        }
        
        char *data_storage = malloc(N * line_size);

        for (i = 0; i < N; i++)
        {
            char *data = data_storage + i * line_size;

            // Copy link id
            memcpy(data, &sys[i].ID, sizeof(unsigned int));
            data += sizeof(unsigned int);

            // Copy states
            assert(sys[i].dim >= dim);
            if (assignments[i] != 0)
                MPI_Recv(data, sys[i].dim, MPI_DOUBLE, assignments[i], i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            else
                memcpy(data, sys[i].my->list.tail->y_approx, sys[i].dim * sizeof(double));
        }

        // Append to the packet table
        herr_t ret = H5PTappend(packet_file_id, N, data_storage);

        //Clean up
        H5PTclose(packet_file_id);
        H5Fclose(file_id);
        H5Tclose(compound_id);
    }
    else			//Sending data to proc 0
    {
        //Check for error while opening the file
        MPI_Bcast(&res, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        if (res)
            return 1;

        for (i = 0; i < N; i++)
        {
            if (assignments[i] == my_rank)
                MPI_Ssend(sys[i].my->list.tail->y_approx, sys[i].dim, MPI_DOUBLE, 0, i, MPI_COMM_WORLD);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    return 0;
}

#if defined(HAVE_POSTGRESQL)

int DumpStateDB(Link* sys, unsigned int N, int* assignments, GlobalVars* globals, char* preface, ConnData* conninfo)
{
    unsigned int i, j, init_length;
    char query[ASYNCH_MAX_QUERY_LENGTH], temptablename[ASYNCH_MAX_QUERY_LENGTH];
    unsigned int size = globals->max_dim*(16 + 1) + 16 + 1 + 16 + 1;
    char submission[ASYNCH_MAX_DIM * (16 + 1) + 16 + 1 + 16 + 1];	//Assumes 16 bytes for each double
    static short int first_call = 0;
    PGresult* res;
    int result, error;

    memset(submission, 0, size);

    //Since forecast_time does not change, just put them in submission once
    if (preface)
    {
        strcpy(submission, preface);
        init_length = (unsigned int)strlen(preface);
        submission[init_length++] = ',';
    }
    else
        init_length = 0;

    //Make a snapshot
    if (my_rank == 0)
    {
        ConnectPGDB(conninfo);
        error = 0;

        if (!first_call)	//!!!! Should this be done in an init routine? Does the file need an init routine? !!!!
        {
            //Trash table
            if (conninfo->num_queries > 1)
            {
                sprintf(query, conninfo->queries[1], globals->dump_table);
                res = PQexec(conninfo->conn, query);
                error |= CheckResError(res, "applying query 1 to snapshot table");
                PQclear(res);
            }

            //Create table
            if (conninfo->num_queries > 0)
            {
                sprintf(query, conninfo->queries[0], globals->dump_table);
                res = PQexec(conninfo->conn, query);
                error |= CheckResError(res, "applying query 0 to snapshot table");
                PQclear(res);
            }

            //Don't create the table again
            first_call = 1;
        }

        if (conninfo->num_queries > 2)	//For truncating the table
        {
            sprintf(query, conninfo->queries[2], globals->dump_table);
            res = PQexec(conninfo->conn, query);
            error |= CheckResError(res, "applying query 2 to snapshot table");
            PQclear(res);
        }

        //Create temp table
        //sprintf(temptablename,"tmp_%s",globals->dump_table);
        sprintf(temptablename, "%s_tmp", globals->dump_table);

        sprintf(conninfo->query, "DROP TABLE IF EXISTS %s;", temptablename);
        res = PQexec(conninfo->conn, conninfo->query);
        CheckResError(res, "dropping temporary snapshot table");
        PQclear(res);

        sprintf(conninfo->query, conninfo->queries[0], temptablename);
        res = PQexec(conninfo->conn, conninfo->query);
        CheckResError(res, "creating temporary snapshot table");
        PQclear(res);

        //Tell database to prepare for copying
        if (!error)
        {
            sprintf(conninfo->query, "COPY %s FROM STDIN WITH DELIMITER ',' NULL AS 'NULL';", temptablename);
            //sprintf(conninfo->query,"COPY %s FROM STDIN WITH DELIMITER ',' NULL AS 'NULL';",globals->dump_table);
            res = PQexec(conninfo->conn, conninfo->query);
            error |= CheckResState(res, PGRES_COPY_IN);
            PQclear(res);
        }

        //Tell all procs if an error occured
        MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);

        //Upload data
        if (!error)
        {
            for (i = 0; i < N; i++)
            {
                int nbytes = init_length;
                if (assignments[i] == my_rank)
                {
                    nbytes += CatBinaryToString(&(submission[nbytes]), "%i", (char*)&(sys[i].ID), ASYNCH_INT, ",");
                    for (j = 0; j < sys[i].dim; j++)
                        nbytes += CatBinaryToString(&(submission[nbytes]), "%.6e", sys[i].my->list.tail->y_approx + j, ASYNCH_DOUBLE, ",");
                    submission[nbytes - 1] = '\n';
                }
                else
                {
                    MPI_Recv(&(submission[init_length]), size - init_length, MPI_CHAR, assignments[i], sys[i].ID, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    nbytes += (int)strlen(&(submission[init_length]));
                }

                result = PQputCopyData(conninfo->conn, submission, nbytes);
                if (result != 1)
                {
                    printf("[%i]: Error: copy returned result %i.\n", my_rank, result);
                    error = 1;
                }
            }
        }

        //Finish copy
        result = PQputCopyEnd(conninfo->conn, NULL);
        if (result != 1)
        {
            printf("[%i]: Returned %i while closing copy to peakflow database.\n", my_rank, result);
            error = 1;
        }

        //If no error, insert temporary table into permanent table
        if (!error)
        {
            sprintf(conninfo->query, "INSERT INTO %s (SELECT * FROM %s);", globals->dump_table, temptablename);
            res = PQexec(conninfo->conn, conninfo->query);
            if (CheckResError(res, "inserting temp snapshot table into permanent table"))
                error = 1;
            PQclear(res);
        }

        //Delete temporary table
        sprintf(conninfo->query, "DROP TABLE IF EXISTS %s;", temptablename);
        res = PQexec(conninfo->conn, conninfo->query);
        CheckResError(res, "dropping temporary snapshot table");
        PQclear(res);

        //Disconnect
        DisconnectPGDB(conninfo);
    }
    else			//Sending data to proc 0
    {
        //See if proc 0 encountered an error
        MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (!error)
        {
            for (i = 0; i < N; i++)
            {
                if (assignments[i] == my_rank)
                {
                    size_t nbytes = init_length;
                    nbytes += CatBinaryToString(&(submission[nbytes]), "%i", (char*)&(sys[i].ID), ASYNCH_INT, ",");
                    for (j = 0; j < sys[i].dim; j++)
                        nbytes += CatBinaryToString(&(submission[nbytes]), "%.6e", sys[i].my->list.tail->y_approx + j, ASYNCH_DOUBLE, ",");
                    submission[nbytes - 1] = '\n';
                    MPI_Send(&(submission[init_length]), size - init_length, MPI_CHAR, 0, sys[i].ID, MPI_COMM_WORLD);
                }
            }
        }
    }

    //See if proc 0 encountered an error
    MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);

    return error;
}

#endif //HAVE_POSTGRESQL

FILE* PrepareTempFiles(Link* sys, unsigned int N, int* assignments, GlobalVars* globals, unsigned int* save_list, unsigned int save_size, unsigned int my_save_size, char* additional, const Lookup * const id_to_loc)
{
    FILE* outputfile = NULL;
    char filename[ASYNCH_MAX_PATH_LENGTH];
    long int current_pos = SEEK_SET;

    const int buffer_size = 1024;
    char buffer[1024];
    memset(buffer, 0, buffer_size);

    //Setup temporary output data file
    if (my_save_size > 0)
    {
        //dummy_value = (double*) malloc((1+globals->num_print)*sizeof(double));
        //for(i=0;i<1+globals->num_print;i++)	dummy_value[i] = 0.0;
        if (additional != NULL)
            sprintf(filename, "%s_%s", globals->temp_filename, additional);
        //sprintf(filename,"%s_%s_%.3i",globals->temp_filename,additional,my_rank);
        else
            sprintf(filename, "%s", globals->temp_filename);
        //sprintf(filename,"%s_%.3i",globals->temp_filename,my_rank);
        outputfile = fopen(filename, "w+b");
        if (outputfile == NULL)
        {
            printf("[%i]: Error: Could not create file %s.\n", my_rank, filename);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        //if(globals->assim_flag == 1)	start = 0;
        //else				start = 1;
        unsigned int start = 0;

        for (unsigned int i = 0; i < save_size; i++)
        {
            unsigned int loc = find_link_by_idtoloc(save_list[i], id_to_loc, N);

            if (assignments[loc] == my_rank)
            {
                Link* current = &sys[loc];
                fwrite(&(current->ID), sizeof(unsigned int), 1, outputfile);

                //Calculate how many steps should be stored
                current->expected_file_vals = (unsigned int)rint(globals->maxtime / current->print_time) + 2;
                fwrite(&(current->expected_file_vals), sizeof(unsigned int), 1, outputfile);
                //fgetpos(outputfile,&(current->pos));
                current_pos += 2 * sizeof(unsigned int);
                current->pos_offset = current_pos;

                unsigned int line_size = 0;
                for (unsigned int j = 0; j < globals->num_outputs; j++)
                    line_size += globals->outputs[j].size;

                long  offset = line_size * current->expected_file_vals;
                while (offset)
                {
                    unsigned int num_bytes = min(buffer_size, offset);
                    fwrite(buffer, sizeof(char), num_bytes, outputfile);
                    offset -= num_bytes;
                }                
                
                //Grow the file (fseek is not agood idea since it keep blowing the cache)
                //int ret = fseek(outputfile, offset, SEEK_CUR);
                
                current_pos += line_size * current->expected_file_vals;
            }
        }

        //Add a few padding bytes to the end of the file.
        //This is to fix an issue with having the temp files open while by proc p while proc 0 reads them.
        //On Helium (at least), the last few bytes are not readable by 0 until p closes the file.
        //for(i=0;i<4;i++)	fwrite(&dummy_t,sizeof(double),1,outputfile);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    return outputfile;
}

//Deletes this process's temporary file.
//Returns 0 if file deleted, errno if an error occurred.
int RemoveTemporaryFiles(GlobalVars* globals, unsigned int my_save_size, char* additional_temp)
{
    int ret_val = 0;
    char filename[ASYNCH_MAX_PATH_LENGTH];

    if (my_save_size > 0)
    {
        //Open the temp file
        if (!additional_temp)
            sprintf(filename, "%s", globals->temp_filename);
        //sprintf(filename,"%s_%.3i",globals->temp_filename,my_rank);
        else
            sprintf(filename, "%s_%s", globals->temp_filename, additional_temp);
        //sprintf(filename,"%s_%s_%.3i",globals->temp_filename,additional_temp,my_rank);

        ret_val = remove(filename);
        if (ret_val == -1)
            ret_val = errno;
    }

    return ret_val;
}

//Resets all temp files to the beginning of each link. This will allow all data in the temp files to be overwritten.
//Returns 2 if there is an error, 1 if a warning, 0 if good, -1 if no tempfile is open.
int ResetTempFiles(double set_time, Link* sys, unsigned int N, FILE* tempfile, GlobalVars* globals, unsigned int my_save_size, const Lookup * const id_to_loc)
{
    unsigned int i, id;
    long int current_pos = SEEK_SET;
    Link* current;

    //If nothing to write, don't reset!
    if (my_save_size == 0)	return 0;

    //Make sure file is open (this may not be an error necessarily)
    if (!tempfile)
    {
        //printf("[%i]: Error: Cannot reset temporary file. A temporary file is not open.\n",my_rank);
        return -1;
    }

    //Restart the file
    rewind(tempfile);

    //Compute the number of bytes in for the component_idx
    unsigned int line_size = CalcTotalOutputSize(globals);

    //Find set_time at each link
    for (i = 0; i < my_save_size; i++)
    {
        //Get linkid from file
        fread(&id, sizeof(unsigned int), 1, tempfile);
        fseek(tempfile, sizeof(unsigned int), SEEK_CUR);	//Skip over number of stored values (I don't think it's needed here)
        current = &sys[find_link_by_idtoloc(id, id_to_loc, N)];
        //current_pos = fgetpos(tempfile,&(current->pos));
        current_pos += 2 * sizeof(unsigned int);
        current->pos_offset = current_pos;

        //Set file position
        //fgetpos(tempfile,&(current->pos));
        current->disk_iterations = 0;
        current->next_save = set_time;		//!!!! This forces the print times to match up with the assimilation times !!!!

        //Get to next link in file
        fseek(tempfile, (current->expected_file_vals)*line_size, SEEK_CUR);
        current_pos += (current->expected_file_vals)*line_size;
    }

    return 0;
}


//Sets all temp files to set_value. The set_value is in the component_idx spot. All data that has been stored later can be overwritten.
//set_time is the time that is used for the next_save times.
//Returns 2 if there is an error, 1 if a warning, 0 if good, -1 if no file is open.
int SetTempFiles(double set_time, void* set_value, enum AsynchTypes data_type, unsigned int component_idx, Link* sys, unsigned int N, FILE* tempfile, GlobalVars* globals, unsigned int my_save_size, const Lookup * const id_to_loc)
{
    unsigned int i, j, id;
    int warning = 0;
    long int current_pos = SEEK_SET;
    //double current_time;
    void* current_value;
    Link* current;

    //printf("!!!! Note: SetTempFiles has not been tested. Be careful... !!!!\n");

    if (!tempfile)
    {
        //printf("[%i]: Error: Cannot reset temporary file. A temporary file is not open.\n",my_rank);
        return -1;
    }

    if (component_idx >= globals->num_outputs)
    {
        printf("[%i]: Error: Cannot reset file to index %u. Only %u values are being outputed.\n", my_rank, component_idx, globals->num_outputs);
        return 2;
    }

    //Restart the file
    rewind(tempfile);

    //Compute the number of bytes in for the component_idx
    unsigned int line_size = CalcTotalOutputSize(globals);
    unsigned int start_offset = 0;
    for (i = 0; i < component_idx; i++)
        start_offset += globals->outputs[i].size;
    
    current_value = malloc(globals->outputs[component_idx].size);
    
    //unsigned int tempy,obtained;
    //int got;
    //double tempy_d;
    
    //Find set_time at each link
    for (i = 0; i < my_save_size; i++)
    {
        //Get linkid from file
        fread(&id, sizeof(unsigned int), 1, tempfile);
        /*
        printf("read id %u (%u)\n",id,obtained);
        if(obtained == 0)
        {
        printf("Error function: %i %i\n",feof(tempfile),ferror(tempfile));
        }
        */
        fseek(tempfile, sizeof(unsigned int), SEEK_CUR);	//Skip over number of stored values (I don't think it's needed here)
/*
obtained = fread(&tempy,sizeof(unsigned int),1,tempfile);
printf("read steps %u (%u)\n",tempy,obtained);
*/
        current = &sys[find_link_by_idtoloc(id, id_to_loc, N)];
        current_pos += 2 * sizeof(unsigned int);
        //printf("Got id = %u, stored %u, disk_iters = %u, start offset = %u\n",id,tempy,current->disk_iterations,start_offset);
        /*
        {
        obtained = fread(&tempy_d,sizeof(double),1,tempfile);
        printf("read %e (%u)\n",tempy_d,obtained);
        obtained = fread(&tempy_d,sizeof(double),1,tempfile);
        printf("read %e (%u)\n",tempy_d,obtained);
        obtained = fread(&tempy_d,sizeof(double),1,tempfile);
        printf("read %e (%u)\n",tempy_d,obtained);
        obtained = fread(&tempy_d,sizeof(double),1,tempfile);
        printf("read %e (%u)\n",tempy_d,obtained);
        obtained = fread(&tempy_d,sizeof(double),1,tempfile);
        printf("read %e (%u)\n",tempy_d,obtained);
        obtained = fread(&tempy_d,sizeof(double),1,tempfile);
        printf("read %e (%u)\n",tempy_d,obtained);
        fseek(tempfile,-6*sizeof(double),SEEK_CUR);
        }
        */
        //Move to component index
        fseek(tempfile, start_offset, SEEK_CUR);

        for (j = 0; j < current->disk_iterations; j++)
        {
            fread(current_value, globals->outputs[component_idx].size, 1, tempfile);
            /*
            if(id == 0)
            {
            printf("Comparing %e to %e. (%u) Got %i.\n",*(double*)set_value,*(double*)current_value,obtained,dt_info->is_equal[data_type](set_value,current_value));
            printf("i = %i j = %i iters = %u\n",i,j,current->disk_iterations);
            //getchar();
            }
            */
            if (DataTypes_IsEqual(data_type, set_value, current_value))
            {
                /*
                if(id == 0)
                {
                printf("breaking...\n");
                printf("i = %i j = %i iters = %u\n",i,j,current->disk_iterations);
                //getchar();
                }
                */
                //fseek(tempfile,-start_offset,SEEK_CUR);	//Backup to start of step
                fseek(tempfile, (long int)-(int)(start_offset + globals->outputs[component_idx].size), SEEK_CUR);	//Backup to start of step
/*
if(id == 0)
{
long int val = (long int) -(int)(start_offset+globals->output_sizes[component_idx]);
printf("Went back %li, got = %i, ferror = %i, feof = %i\n",val,got,ferror(tempfile),feof(tempfile));
obtained = fread(&tempy_d,sizeof(double),1,tempfile);
printf("checky: read %e (%u %i)\n",tempy_d,obtained,feof(tempfile));
obtained = fread(&tempy_d,sizeof(double),1,tempfile);
printf("checky: read %e (%u %i)\n",tempy_d,obtained,feof(tempfile));
fseek(tempfile,-2*sizeof(double),SEEK_CUR);
}
*/
                break;
            }

            fseek(tempfile, line_size - globals->outputs[component_idx].size, SEEK_CUR);	//Skip to next step
            current_pos += line_size;
        }

        if (j == current->disk_iterations)
        {
            printf("[%i]: Warning: Cannot set file position for link %u.\n", my_rank, id);
            warning = 1;
            //continue;
/*
if(id == 0)
{
printf("Failed on %e\n",*(double*)set_value);
printf("i = %i j = %i iters = %u\n",i,j,current->disk_iterations);
//getchar();
}
*/
        }

        //Set file position
        //fgetpos(tempfile,&(current->pos));
        current->pos_offset = current_pos;
        current->disk_iterations = j;
        current->next_save = set_time;		//!!!! This forces the print times to match up with the assimilation times !!!!

/*
if(id == 0)
{
printf("i = %i j = %i line_size = %u expected = %u shift = %i pos_offset = %i\n",i,j,line_size,current->expected_file_vals,(current->expected_file_vals-j)*line_size,current_pos);
//getchar();
}
*/
//Get to next link in file
        fseek(tempfile, (current->expected_file_vals - j)*line_size, SEEK_CUR);
        current_pos += (current->expected_file_vals - j)*line_size;
        //printf("Seek gave %i, current_pos = %i\n",got,current_pos);
    }

    free(current_value);

    return warning;
}

//!!!! This assumes the first value in the file is time. !!!!
/*
int SetTempFiles(double set_time,Link** sys,unsigned int N,FILE* tempfile,UnivVars* GlobalVars,unsigned int my_save_size,const Lookup * const id_to_loc)
{
    unsigned int i,j,id;
    int current_pos,warning = 0;
    double current_time;
    Link* current;

    if(!tempfile)
    {
        printf("[%i]: Error: Cannot reset temporary file. A temporary file is not open.\n",my_rank);
        return 2;
    }

    //Restart the file
    rewind(tempfile);

    //Find set_time at each link
    for(i=0;i<my_save_size;i++)
    {
        //Get linkid from file
        fread(&id,sizeof(unsigned int),1,tempfile);
        fseek(tempfile,sizeof(unsigned int),SEEK_CUR);	//Skip over number of stored values (I don't think it's needed here)
        current = sys[find_link_by_idtoloc(id,id_to_loc,N)];
        current_pos = fgetpos(tempfile,&(current->pos));
        for(j=0;j<current->disk_iterations;j++)
        {
            fread(&current_time,sizeof(double),1,tempfile);
            if( fabs(set_time - current_time) < 1e-10 )		//!!!! Is this a good bound? Maybe use a relative error. !!!!
            {
                fseek(tempfile,-sizeof(double),SEEK_CUR);	//Backup one double
                break;
            }

            //fseek(tempfile,globals->num_print*sizeof(double),SEEK_CUR);	//Skip ahead to next timestep
            fseek(tempfile,(globals->num_print-1)*sizeof(double),SEEK_CUR);	//Skip ahead to next timestep
        }

        if(j == current->disk_iterations)
        {
            printf("[%i]: Warning: Cannot set file position for link %u to time %f.\n",my_rank,id,set_time);
            warning = 1;
            continue;
        }

        //Set file position
        fgetpos(tempfile,&(current->pos));
        current->disk_iterations = j;
        current->next_save = set_time;		//!!!! This forces the print times to match up with the assimilation times !!!!

        //Get to next link in file
        //fseek(tempfile,(current->expected_file_vals-j)*(1+globals->num_print)*sizeof(double),SEEK_CUR);
        fseek(tempfile,(current->expected_file_vals-j)*globals->num_print*sizeof(double),SEEK_CUR);
    }

    return warning;
}
*/

//Read in a .rec file from disk and loads it into the intial condition for sys (tail).
//This does NOT set the current time at each link.
void LoadRecoveryFile(char* filename, Link* sys, unsigned int N, unsigned int my_N, unsigned int* assignments, GlobalVars* globals)
{
    MPI_Barrier(MPI_COMM_WORLD);

    if (my_rank == 0)
    {
        FILE *input = fopen(filename, "r");
        if (!input)
        {
            printf("[%i]: Error opening recovery file %s.\n", my_rank, filename);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        unsigned int read_type, read_N;
        fscanf(input, "%u %u %*f", &read_type, &read_N);
        if (N != read_N)		//!!!! Ignoring read_type for now (190 vs 19) !!!!
        {
            printf("[%i]: Error reading recovery file: bad model type (%i) or wrong number of links (%i).\n", my_rank, read_type, read_N);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        unsigned int max_dim = 0;
        for (unsigned int i = 0; i < N; i++)
            max_dim = max(max_dim, sys[i].dim);

        double *buffer = malloc(max_dim * sizeof(double));
        for (unsigned int i = 0; i < N; i++)
        {
            unsigned int id;
            fscanf(input, "%u", &id);
            if (sys[i].ID != id)
            {
                printf("[%i]: Error reading recovery file: bad link id (%i); expected %i.\n", my_rank, id, sys[i].ID);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }

            for (unsigned int j = 0; j < sys[i].dim; j++)
                fscanf(input, "%lf", buffer + j);

            if (assignments[i] == my_rank)
                dcopy(buffer, sys[i].my->list.tail->y_approx, 0, sys[i].dim);
            else
            {
                MPI_Send(&i, 1, MPI_UNSIGNED, assignments[i], 0, MPI_COMM_WORLD);
                MPI_Send(buffer, sys[i].dim, MPI_DOUBLE, assignments[i], 0, MPI_COMM_WORLD);
            }
        }
        free(buffer);
    }
    else
    {
        for (unsigned int i = 0; i < N; i++)
        {
            unsigned int j;
            MPI_Recv(&j, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(sys[j].my->list.tail->y_approx, sys[j].dim, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

//Rewrites the previous step written at link_i with the current step. The current time and state is used.
//!!!! This assumes we are writting into files. !!!!
//Returns 0 if all is well
//Returns 1 if there is no previous iteration to overwrite
//Returns 2 a step as been previously written, but it is not the expected number of bytes
int overwrite_last_step(Link* link_i, GlobalVars *globals, FILE* outputfile)
{
    long step_byte_size = 0;

    //Check that something has actually been written for this link
    if (link_i->disk_iterations == 0)	return 1;

    //Backup a step in the file
    for (unsigned int i = 0; i < globals->num_outputs; i++)
        step_byte_size += globals->outputs[i].size;
    if (link_i->pos_offset < step_byte_size)	return 2;
    link_i->pos_offset -= step_byte_size;

    //Write the current step
    WriteStep(globals->outputs, globals->num_outputs, outputfile, link_i->ID, link_i->last_t, link_i->my->list.tail->y_approx, link_i->dim, &(link_i->pos_offset));
    return 0;
}


