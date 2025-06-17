#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#if defined(HAVE_POSTGRESQL)
#include <libpq-fe.h>
#endif

#if defined(HAVE_LIBZ)
#include <zlib.h>
#endif


#include <globals.h>
#include <compression.h>
#include <comm.h>
#include <date_manip.h>
#include <db.h>
#include <sort.h>
#include <forcings_io.h>


//This reads in a set of binary files for the rainfall at each link.
//Assumes the file is full of floats. Assumes no IDs are in the file and that IDs are consecutive starting from 0
//Link** sys: An array of links.
//int N: The number of links in sys.
//int my_N: The number of links assigned to this process.
//UnivVars* GlobalVars: Contains all the information that is shared by every link in the system.
//int* my_sys: Array of links assigned to this process (value is location in sys array).
//int* assignments (set by this method): Will be an array with N entries. assignments[i] will the process link sys[i] is assigned to.
//char strfilename[]: String of the filename for the rain file to read (NOT .str files). Filenames should be indexed.
//unsigned int first: The index of the file to read first.
//unsigned int last: The index of the file to read last.
//double t_0: The time at which the first file starts.
//double increment: The amount of time between consecutively indexed files.
//int** id_to_loc (set by this method): Will be an array with N rows and 2 columns, sorted by first col. First col is a link id and second is
//				the location of the id in sys.
//unsigned int max_files: The maximum number of files to be read.
int Create_Rain_Data_Par(
    Link *sys, unsigned int N,
    Link **my_sys, unsigned int my_N,
    const GlobalVars * const globals,
    int* assignments,
    char strfilename[],
    unsigned int first, unsigned int last,
    double t_0, double increment,
    Forcing* forcing, const Lookup * const id_to_loc, unsigned int max_files, unsigned int forcing_idx)
{
    unsigned int i, j, curr_idx;
    unsigned int k;
    Link* current;
    float forcing_buffer;
    unsigned int holder;
    char filename[ASYNCH_MAX_PATH_LENGTH];
    FILE* stormdata = NULL;
    unsigned int numfiles = last - first + 1;

    //This is a time larger than any time in which the integrator is expected to get
    double ceil_time = 1e300;
    if (my_sys[0]->last_t > ceil_time * 0.1)
        printf("[%i]: Warning: integrator time is extremely large (about %e). Loss of precision may occur.\n", my_rank, my_sys[0]->last_t);

    for (i = 0; i < my_N; i++)
    {
        my_sys[i]->my->forcing_data[forcing_idx].data = malloc((numfiles + 1) * sizeof(DataPoint));
        my_sys[i]->my->forcing_data[forcing_idx].num_points = numfiles + 1;
    }

    //Read through the files.
    for (k = 0; k < numfiles; k++)
    {
        //sprintf(filename,"%s%i",strfilename,first+k);		//This is the usual one to use
        sprintf(filename, "%s%i", strfilename, first + k);
        //sprintf(filename,"%srain%i",strfilename,first+k);
        //sprintf(filename,"%sfile-%i",strfilename,first+k);
        stormdata = fopen(filename, "r");
        if (stormdata == NULL)	printf("[%i]: Error opening file %s\n", my_rank, filename);

        for (i = 0; i < N; i++)
        {
            //Find the location of ID i (should be i by assumption)
            //curr_idx = id_to_loc[i][1];
            curr_idx = i;

            //Read the index
            if (assignments[curr_idx] == my_rank)	//If information about the link is needed on this process
            {
                //Read in the storm data for this link
                fread(&forcing_buffer, sizeof(float), 1, stormdata);

                //This assumes different endianness
                holder = *(unsigned int*)&forcing_buffer;	//Pointers are fun
                holder = (((holder & 0x0000ffff) << 16) | ((holder & 0xffff0000) >> 16));
                holder = (((holder & 0x00ff00ff) << 8) | ((holder & 0xff00ff00) >> 8));
                forcing_buffer = *(float*)&holder;

                //Store the data
                sys[curr_idx].my->forcing_data[forcing_idx].data[k].time = (float) (t_0 + k*increment);
                sys[curr_idx].my->forcing_data[forcing_idx].data[k].value = forcing_buffer;
            }
            else		//This link data is not needed on this process.
                fread(&forcing_buffer, sizeof(float), 1, stormdata);
        }

        fclose(stormdata);
    }

    if (my_rank == 0)
        printf("Read %i binary files.\n", numfiles);

    //Add in terms for no rainfall if max_files > numfiles
    for (i = 0; i < my_N; i++)
    {
        curr_idx = my_sys[i]->location;
        for (j = numfiles; j < max_files; j++)
        {
            sys[curr_idx].my->forcing_data[forcing_idx].data[j].time = sys[curr_idx].my->forcing_data[forcing_idx].data[j - 1].time + .0001;
            sys[curr_idx].my->forcing_data[forcing_idx].data[j].value = 0.0;
        }
    }

    //Add a ceiling term
    for (i = 0; i < my_N; i++)
    {
        curr_idx = my_sys[i]->location;
        //sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[max_files].time = GlobalVars->maxtime + 1.0;
        sys[curr_idx].my->forcing_data[forcing_idx].data[max_files].time = ceil_time;
        sys[curr_idx].my->forcing_data[forcing_idx].data[max_files].value = -1.0;
    }

    //Calculate the first rain change time and set rain_value
    for (i = 0; i < my_N; i++)
    {
        current = my_sys[i];
        forcing_buffer = current->my->forcing_data[forcing_idx].data[0].value;
        current->my->forcing_values[forcing_idx] = forcing_buffer;
        current->my->forcing_indices[forcing_idx] = 0;

        for (j = 1; j < current->my->forcing_data[forcing_idx].num_points; j++)
        {
            if (fabs(forcing_buffer - current->my->forcing_data[forcing_idx].data[j].value) > 1e-14)
            {
                current->my->forcing_change_times[forcing_idx] = current->my->forcing_data[forcing_idx].data[j].time;
                break;
            }
        }
        if (j == current->my->forcing_data[forcing_idx].num_points)
        {
            current->my->forcing_change_times[forcing_idx] = current->my->forcing_data[forcing_idx].data[j - 1].time;
        }
    }

    return 0;
}

//This reads in a set of irregular binary files for the rainfall at each link.
//Assumes the file is full of floats. Assumes no IDs are in the file and that IDs are consecutive starting from 0
//Link** sys: An array of links.
//int N: The number of links in sys.
//int my_N: The number of links assigned to this process.
//UnivVars* GlobalVars: Contains all the information that is shared by every link in the system.
//int* my_sys: Array of links assigned to this process (value is location in sys array).
//int* assignments (set by this method): Will be an array with N entries. assignments[i] will the process link sys[i] is assigned to.
//char strfilename[]: String of the filename for the rain file to read (NOT .str files). Filenames should be indexed.
//unsigned int first: The index of the file to read first.
//unsigned int last: The index of the file to read last.
//double t_0: The time at which the first file starts.
//double increment: The amount of time between consecutively indexed files.
//int** id_to_loc (set by this method): Will be an array with N rows and 2 columns, sorted by first col. First col is a link id and second is
//				the location of the id in sys.
//unsigned int max_files: The maximum number of files to be read.
int Create_Rain_Data_Par_IBin(
    Link *sys, unsigned int N,
    Link **my_sys, unsigned int my_N,
    const GlobalVars * const globals,
    int* assignments,
    char strfilename[],
    unsigned int first, unsigned int last,
    double t_0, double increment,
    Forcing* forcing, const Lookup * const id_to_loc, unsigned int max_files, unsigned int forcing_idx)
{
    unsigned int i, j, curr_idx, k, linkid_buffer;
    unsigned int total_entry_index, total_entry_counts, entry_counts, entry_index;
    Link* current;
    float forcing_buffer;
    unsigned int holder;
    char filename[ASYNCH_MAX_PATH_LENGTH];
    FILE* stormdata = NULL;
    unsigned int numfiles = ((last - first) / (forcing->file_time * 60));
    int* ibin_link_id;
    float* ibin_rain_intens;
    int** ibin_link_id_set;
    float** ibin_rain_intens_set;
    unsigned int* entry_count_set;
    unsigned int* timestamp_set;
    unsigned int* timestamps;
    unsigned int c_timestamp;
    unsigned int next_exp_time;
    int received_time;

    // do something
    unsigned int *total_times = (unsigned int*)calloc(my_N, sizeof(unsigned int));
    for (i = 0; i < my_N; i++) {
        total_times[i] = my_sys[i]->my->forcing_data[forcing_idx].num_points;
        my_sys[i]->my->forcing_data[forcing_idx].num_points = 1;
    }

    // This is a time larger than any time in which the integrator is expected to get
    double ceil_time = 1e300;
    if (my_sys[0]->last_t > ceil_time * 0.1)
        printf("[%i]: Warning: integrator time is extremely large (about %e). Loss of precision may occur.\n", my_rank, my_sys[0]->last_t);

    if (my_rank == 0) {
        // main process: read the files and broadcast the data

        // Open memory for file buffers
        if (numfiles > 0) {
            ibin_link_id_set = malloc(numfiles * sizeof(int*));
            ibin_rain_intens_set = malloc(numfiles * sizeof(float*));
            entry_count_set = malloc(numfiles * sizeof(unsigned int));
            timestamp_set = malloc(numfiles * sizeof(unsigned int));
        }

        // main process open the files and save them locally
        for (k = 0; k < numfiles; k++) {
            // open file and read the number of entries
            c_timestamp = first + (k*(forcing->file_time*60));
            sprintf(filename, "%s%i", strfilename, c_timestamp);  // define file name
            timestamp_set[k] = c_timestamp;
            stormdata = fopen(filename, "rb");
            if (stormdata == NULL) {
                entry_counts = 0;
            }
            else {
                fread(&entry_counts, sizeof(int), 1, stormdata);
            }

            // allocate memory for those entries
            entry_count_set[k] = entry_counts;
            ibin_link_id_set[k] = malloc(entry_counts * sizeof(unsigned int));
            ibin_rain_intens_set[k] = malloc(entry_counts * sizeof(float));

            // read memory and fill the temporary data
            if (entry_counts > 0) {
                for (entry_index = 0; entry_index < entry_counts; entry_index++) {
                    fread(&linkid_buffer, sizeof(unsigned int), 1, stormdata);
                    fread(&forcing_buffer, sizeof(float), 1, stormdata);

                    ibin_link_id_set[k][entry_index] = linkid_buffer;
                    ibin_rain_intens_set[k][entry_index] = forcing_buffer;
                }
            }
            if (stormdata != NULL) {
                fclose(stormdata);
            }
        }

        // count total number of entries and broadcast it
        total_entry_counts = 0;
        for (k = 0; k < numfiles; k++) {
            total_entry_counts += entry_count_set[k];
        }
        MPI_Bcast(&total_entry_counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // process if there is something to send
        if (total_entry_counts > 0) {
            // allocates final total space
            ibin_link_id = malloc(total_entry_counts * sizeof(unsigned int));
            ibin_rain_intens = malloc(total_entry_counts * sizeof(float));
            timestamps = malloc(total_entry_counts * sizeof(unsigned int));

            // reorganizes all data readen into a single buffer vector
            total_entry_index = 0;
            for (k = 0; k < numfiles; k++) {
                for (entry_index = 0; entry_index < entry_count_set[k]; entry_index++) {
                    ibin_link_id[total_entry_index] = ibin_link_id_set[k][entry_index];
                    ibin_rain_intens[total_entry_index] = ibin_rain_intens_set[k][entry_index];
                    timestamps[total_entry_index] = timestamp_set[k];
                    total_entry_index += 1;
                }
            }

            // broadcast-send single buffers
            MPI_Bcast(ibin_link_id, total_entry_counts, MPI_FLOAT, 0, MPI_COMM_WORLD);
            MPI_Bcast(ibin_rain_intens, total_entry_counts, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(timestamps, total_entry_counts, MPI_INT, 0, MPI_COMM_WORLD);
        }

        // free memory of "(...)_set" variables
        for (k = 0; k < numfiles; k++) {
            free(ibin_link_id_set[k]);
            free(ibin_rain_intens_set[k]);
        }
        
        if (numfiles > 0) {
            free(ibin_link_id_set);
            free(ibin_rain_intens_set);
            free(timestamp_set);
        }
    }
    else
    {
        // other processes: receive broadcasted data
        MPI_Bcast(&total_entry_counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // process if there is something to send
        if (total_entry_counts > 0) {
            // allocates final total space
            ibin_link_id = malloc(total_entry_counts * sizeof(unsigned int));
            ibin_rain_intens = malloc(total_entry_counts * sizeof(float));
            timestamps = malloc(total_entry_counts * sizeof(unsigned int));

            // broadcast-receive single buffers
            MPI_Bcast(ibin_link_id, total_entry_counts, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(ibin_rain_intens, total_entry_counts, MPI_FLOAT, 0, MPI_COMM_WORLD);
            MPI_Bcast(timestamps, total_entry_counts, MPI_INT, 0, MPI_COMM_WORLD);
        }
    } // end file reading

    MPI_Barrier(MPI_COMM_WORLD);

    //Setup initial time in rainfall data with 0 rain
    for (i = 0; i < my_N; i++)
    {
        my_sys[i]->my->forcing_data[forcing_idx].data[0].time = (int)(first - forcing->first_file) / 60.0;  // in minutes from 0
        my_sys[i]->my->forcing_data[forcing_idx].data[0].value = 0.0;
    }

    // import data from buffer
    for (i = 0; i < total_entry_counts; i++)
    {
        curr_idx = find_link_by_idtoloc(ibin_link_id[i], id_to_loc, N);

        if (curr_idx < N && assignments[curr_idx] == my_rank)
        {
            k = sys[curr_idx].my->forcing_data[forcing_idx].num_points;
            received_time = (timestamps[i] - forcing->first_file) / 60.0;  // minutes - time just read

            if (received_time >(int) (sys[curr_idx].my->forcing_data[forcing_idx].data[k - 1].time + 0.01)) // if just read something after last timestamp in the timeseries
            {
                next_exp_time = (int)(sys[curr_idx].my->forcing_data[forcing_idx].data[k - 1].time) + (unsigned int)(forcing->file_time);
                if ((received_time <= next_exp_time) || (sys[curr_idx].my->forcing_data[forcing_idx].data[k - 1].value == 0.0))
                {
                    // add the read value in the end of the time series
                    sys[curr_idx].my->forcing_data[forcing_idx].data[k].time = received_time;
                    sys[curr_idx].my->forcing_data[forcing_idx].data[k].value = ibin_rain_intens[i];
                    sys[curr_idx].my->forcing_data[forcing_idx].num_points += 1;
                }
                else	//Add an intermediate 0 rainfall data
                {
                    // add the read value in the end of the time series
                    sys[curr_idx].my->forcing_data[forcing_idx].data[k].time = next_exp_time;
                    sys[curr_idx].my->forcing_data[forcing_idx].data[k].value = 0.0;
                    sys[curr_idx].my->forcing_data[forcing_idx].data[k + 1].time = received_time;
                    sys[curr_idx].my->forcing_data[forcing_idx].data[k + 1].value = ibin_rain_intens[i];
                    sys[curr_idx].my->forcing_data[forcing_idx].num_points += 2;
                }
            }
            else if (received_time <= (int)(sys[curr_idx].my->forcing_data[forcing_idx].data[k - 1].time + 0.01))	//If the initial rate needs to be reset
            {
                sys[curr_idx].my->forcing_data[forcing_idx].data[k - 1].time = received_time;
                sys[curr_idx].my->forcing_data[forcing_idx].data[k - 1].value = ibin_rain_intens[i];
            }
            else
            {
                printf("!!!! Uh oh... !!!!\n");
                printf("!!!! i = %i k = %i received = %i unix_time = %i raindb_start = %i\n", i, k, received_time, timestamps[i], forcing->raindb_start_time);
            }

        }
    }

    //Add ceiling terms
    for (i = 0; i < my_N; i++)
    {
        curr_idx = my_sys[i]->location;
        k = sys[curr_idx].my->forcing_data[forcing_idx].num_points;
        if (sys[curr_idx].my->forcing_data[forcing_idx].data[k - 1].value == 0.0)	//No rain, add just a ceiling
        {
            //sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k].time = maxtime * (1.1) + 1.0;
            sys[curr_idx].my->forcing_data[forcing_idx].data[k].time = ceil_time;
            sys[curr_idx].my->forcing_data[forcing_idx].data[k].value = -1.0;
        }
        else	//Add a 0.0, and a ceiling
        {
            sys[curr_idx].my->forcing_data[forcing_idx].data[k].time = sys[curr_idx].my->forcing_data[forcing_idx].data[k - 1].time + forcing->file_time;
            sys[curr_idx].my->forcing_data[forcing_idx].data[k].value = 0.0;

            sys[curr_idx].my->forcing_data[forcing_idx].data[k + 1].time = ceil_time;
            sys[curr_idx].my->forcing_data[forcing_idx].data[k + 1].value = -1.0;
        }
    }

    for (i = 0; i < my_N; i++)	my_sys[i]->my->forcing_data[forcing_idx].num_points = total_times[i];

    //Calculate the first rain change time and set rain_value
    for (i = 0; i < my_N; i++)
    {

        current = my_sys[i];
        forcing_buffer = current->my->forcing_data[forcing_idx].data[0].value; // why this? data[0] or data[1]?
        current->my->forcing_values[forcing_idx] = forcing_buffer;
        current->my->forcing_indices[forcing_idx] = 0;

        for (j = 1; j < current->my->forcing_data[forcing_idx].num_points; j++)
        {
            if (fabs(forcing_buffer - current->my->forcing_data[forcing_idx].data[j].value) > 1e-14)
            {
                current->my->forcing_change_times[forcing_idx] = current->my->forcing_data[forcing_idx].data[j].time;
                break;
            }
        }
        if (j == current->my->forcing_data[forcing_idx].num_points)
        {
            current->my->forcing_change_times[forcing_idx] = current->my->forcing_data[forcing_idx].data[j - 1].time;
        }
    }

    return 0;
}

#define ASYNCH_BUFFER_SIZE sizeof(unsigned int) + sizeof(float)

/// This reads in a set of gzip compressed binary files for the rainfall at each link.
/// Assumes the file is full of floats. Assumes no IDs are in the file and that IDs are consecutive starting from 0
///
/// \param sys An array of links.
/// \param N   The number of links in sys.
/// \param my_sys An array of pointer to links.
/// \param my_N The number of links assigned to this process.
/// \param globals  Contains all the information that is shared by every link in the system.
/// \param assignments (set by this method): Will be an array with N entries. assignments[i] will the process link sys[i] is assigned to.
/// \param strfilename: String of the filename for the rain file to read (NOT .str files). Filenames should be indexed.
/// \param first: The index of the file to read first.
/// \param last: The index of the file to read last.
/// \param t_0: The time at which the first file starts.
/// \param increment: The amount of time between consecutively indexed files.
/// \param id_to_loc (set by this method): Will be an array with N rows and 2 columns, sorted by first col. First col is a link id and second is
///				the location of the id in sys.
//unsigned int max_files: The maximum number of files to be read.
int Create_Rain_Data_GZ(
    Link *sys, unsigned int N,
    Link **my_sys, unsigned int my_N,
    const GlobalVars * const globals, int* assignments, char strfilename[], unsigned int first, unsigned int last, double t_0, double increment,
    Forcing* forcing,
    const Lookup * const id_to_loc,
    unsigned int max_files,
    unsigned int forcing_idx)
{
    unsigned int i, j, curr_idx;
    unsigned int k;
    Link* current;
    float rainfall_buffer;
    unsigned int holder;
    char filename[128];
    FILE* stormdata = NULL;
    unsigned int numfiles = last - first + 1;
    FILE* compfile = NULL;
    char transferbuffer[ASYNCH_BUFFER_SIZE];

    //This is a time larger than any time in which the integrator is expected to get
    double ceil_time = 1e300;
    if (my_sys[0]->last_t > ceil_time*0.1)
        printf("[%i]: Warning: integrator time is extremely large (about %e). Loss of precision may occur.\n", my_rank, my_sys[0]->last_t);

    //Check that space for rain data has been allocated.
    for (i = 0; i < my_N; i++)
    {
        my_sys[i]->my->forcing_data[forcing_idx].data = malloc((numfiles + 1) * sizeof(DataPoint));
        my_sys[i]->my->forcing_data[forcing_idx].num_points = numfiles + 1;
    }

    //Read through the files.
    MPI_Barrier(MPI_COMM_WORLD);
    for (k = 0; k < numfiles; k++)
    {
        if (my_rank == 0)
        {
            //For compressed file
            sprintf(filename, "%s%i.gz", strfilename, first + k);
            compfile = fopen(filename, "r");
            if (!compfile)
            {
                printf("[%i]: Error opening file %s.\n", my_rank, filename);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }

            sprintf(filename, "%s%i_%i", strfilename, first + k, my_rank);
            stormdata = fopen(filename, "w");
            if (!stormdata)	printf("[%i]: Error opening file %s.\n", my_rank, filename);

            int ret = uncompress_gzfile(compfile, stormdata);
            if (ret != Z_OK)		zerr(ret);
            fclose(compfile);
            fclose(stormdata);

            sprintf(filename, "%s%i_%i", strfilename, first + k, my_rank);
            stormdata = fopen(filename, "r");
            if (!stormdata)	printf("[%i]: Error opening file %s.\n", my_rank, filename);

            for (i = 0; i < N; i++)
            {
                //Find the location of ID i (should be i by assumption)
                //curr_idx = id_to_loc[i][1];
                curr_idx = i;

                fread(&rainfall_buffer, sizeof(float), 1, stormdata);

                if (assignments[curr_idx] == my_rank)	//Data needed on proc 0; don't send
                {
                    //This assumes the files have a different endianness from the system
                    holder = *(unsigned int*)&rainfall_buffer;	//Pointers are fun!
                    holder = (((holder & 0x0000ffff) << 16) | ((holder & 0xffff0000) >> 16));
                    holder = (((holder & 0x00ff00ff) << 8) | ((holder & 0xff00ff00) >> 8));
                    rainfall_buffer = *(float*)&holder;

                    //Store the data
                    sys[curr_idx].my->forcing_data[forcing_idx].data[k].time = t_0 + k*increment;
                    sys[curr_idx].my->forcing_data[forcing_idx].data[k].value = rainfall_buffer;
                }
                else	//Send it to the correct proc
                {
                    int pos = 0;
                    MPI_Pack(&curr_idx, 1, MPI_INT, transferbuffer, ASYNCH_BUFFER_SIZE, &pos, MPI_COMM_WORLD);
                    MPI_Pack(&rainfall_buffer, 1, MPI_FLOAT, transferbuffer, ASYNCH_BUFFER_SIZE, &pos, MPI_COMM_WORLD);
                    MPI_Send(transferbuffer, ASYNCH_BUFFER_SIZE, MPI_PACKED, assignments[curr_idx], N, MPI_COMM_WORLD);		//N = tag, as typical communication should not send this many links.
                }
            }

            fclose(stormdata);
            remove(filename);
        }
        else
        {
            //MPI_Status status;
            for (i = 0; i < my_N; i++)
            {
                int pos = 0;
                MPI_Recv(transferbuffer, ASYNCH_BUFFER_SIZE, MPI_PACKED, 0, N, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Unpack(transferbuffer, ASYNCH_BUFFER_SIZE, &pos, &curr_idx, 1, MPI_INT, MPI_COMM_WORLD);
                MPI_Unpack(transferbuffer, ASYNCH_BUFFER_SIZE, &pos, &rainfall_buffer, 1, MPI_FLOAT, MPI_COMM_WORLD);

                //This assumes the files have a different endianness from the system
                holder = *(unsigned int*)&rainfall_buffer;	//Pointers are fun!
                holder = (((holder & 0x0000ffff) << 16) | ((holder & 0xffff0000) >> 16));
                holder = (((holder & 0x00ff00ff) << 8) | ((holder & 0xff00ff00) >> 8));
                rainfall_buffer = *(float*)&holder;

                //Store the data
                sys[curr_idx].my->forcing_data[forcing_idx].data[k].time = t_0 + k*increment;
                sys[curr_idx].my->forcing_data[forcing_idx].data[k].value = rainfall_buffer;
            }
        }
    }

    if (my_rank == 0)
        printf("Read %i binary files.\n", numfiles);

    //Add in terms for no rainfall if max_files > numfiles
    for (i = 0; i < my_N; i++)
    {
        curr_idx = my_sys[i]->location;
        for (j = numfiles; j < max_files; j++)
        {
            sys[curr_idx].my->forcing_data[forcing_idx].data[j].time = sys[curr_idx].my->forcing_data[forcing_idx].data[j - 1].time + .0001;
            sys[curr_idx].my->forcing_data[forcing_idx].data[j].value = 0.0;
        }
    }

    //Add a ceiling term
    for (i = 0; i < my_N; i++)
    {
        curr_idx = my_sys[i]->location;
        //sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[max_files].time = GlobalVars->maxtime + 1.0;
        sys[curr_idx].my->forcing_data[forcing_idx].data[max_files].time = ceil_time;
        sys[curr_idx].my->forcing_data[forcing_idx].data[max_files].value = -1.0;
    }

    //Calculate the first rain change time and set rain_value
    for (i = 0; i < my_N; i++)
    {
        current = my_sys[i];
        rainfall_buffer = current->my->forcing_data[forcing_idx].data[0].value;
        current->my->forcing_values[forcing_idx] = rainfall_buffer;
        current->my->forcing_indices[forcing_idx] = 0;

        for (j = 1; j < current->my->forcing_data[forcing_idx].num_points; j++)
        {
            if (rainfall_buffer != current->my->forcing_data[forcing_idx].data[j].value)
            {
                current->my->forcing_change_times[forcing_idx] = current->my->forcing_data[forcing_idx].data[j].time;
                break;
            }
        }
        if (j == current->my->forcing_data[forcing_idx].num_points)
            current->my->forcing_change_times[forcing_idx] = current->my->forcing_data[forcing_idx].data[j - 1].time;
    }

    return 0;
}


//This reads in a set of binary files for the rainfall at each link.
//The data is given as intensities per grid cell.
//Link** sys: An array of links.
//int N: The number of links in sys.
//int my_N: The number of links assigned to this process.
//UnivVars* GlobalVars: Contains all the information that is shared by every link in the system.
//int* my_sys: Array of links assigned to this process (value is location in sys array).
//int* assignments (set by this method): Will be an array with N entries. assignments[i] will the process link sys[i] is assigned to.
//char strfilename[]: String of the filename for the rain file to read (NOT .str files). Filenames should be indexed.
//unsigned int first: The index of the file to read first.
//unsigned int last: The index of the file to read last.
//double t_0: The time at which the first file starts.
//double increment: The amount of time between consecutively indexed files.
//int** id_to_loc (set by this method): Will be an array with N rows and 2 columns, sorted by first col. First col is a link id and second is
//				the location of the id in sys.
//unsigned int max_files: The maximum number of files to be read.

int Create_Rain_Data_Grid(
    Link *sys, unsigned int N,
    Link **my_sys, unsigned int my_N,
    const GlobalVars * const globals, int* assignments, char strfilename[], unsigned int first, unsigned int last, double t_0, double increment, Forcing* forcing, const Lookup * const id_to_loc, unsigned int max_files, unsigned int forcing_idx)
{
    unsigned int i, j, curr_idx, k, holder, endianness, cell;
    short unsigned int intensity;
    Link* current;
    char filename[128];
    float forcing_buffer;
    FILE* stormdata = NULL;
    unsigned int numfiles = last - first + 1;
    size_t result;

    //This is a time larger than any time in which the integrator is expected to get
    double ceil_time = 1e300;
    if (my_sys[0]->last_t > ceil_time*0.1)
        printf("[%i]: Warning: integrator time is extremely large (about %e). Loss of precision may occur.\n", my_rank, my_sys[0]->last_t);

    //Check that space for rain data has been allocated.
    for (i = 0; i < my_N; i++)
    {
        my_sys[i]->my->forcing_data[forcing_idx].data = malloc((numfiles + 1) * sizeof(DataPoint));
        my_sys[i]->my->forcing_data[forcing_idx].num_points = numfiles + 1;
    }

    //Read through the files.
    for (k = 0; k < numfiles; k++)
    {
        if (my_rank == 0)
        {
            sprintf(filename, "%s%i", strfilename, first + k);
            stormdata = fopen(filename, "r");
            if (stormdata)
            {
                for (i = 0; i < forcing->num_cells; i++)
                    forcing->received[i] = 0;

                //Check endianness
                fread(&i, sizeof(unsigned int), 1, stormdata);
                if (i == 0x1)			endianness = 0;
                else if (i == 0x80000000)	endianness = 1;
                else
                {
                    printf("Error: Cannot read endianness flag in binary file %s.\n", filename);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }

                //Read file
                if (endianness)
                {
                    while (!feof(stormdata))
                    {
                        //Read intensity
                        result = fread(&cell, sizeof(unsigned int), 1, stormdata);
                        if (!result)	break;
                        fread(&intensity, sizeof(short unsigned int), 1, stormdata);

                        //Swap byte order
                        holder = (((cell & 0x0000ffff) << 16) | ((cell & 0xffff0000) >> 16));
                        cell = (((holder & 0x00ff00ff) << 8) | ((holder & 0xff00ff00) >> 8));
                        intensity = (((intensity & 0x00ff00ff) << 8) | ((intensity & 0xff00ff00) >> 8));

                        if (cell < forcing->num_cells)
                        {
                            if (forcing->received[cell])
                                printf("Warning: Received multiple intensities for cell %u in file %s.\n", cell, filename);
                            forcing->received[cell] = 1;
                            forcing->intensities[cell] = (float)(intensity * forcing->factor);
                        }
                        else{
                            printf("Warning: bad grid cell id in file %s.\n", filename);
                            printf("Read cell id %u num_cells is %u.\n", cell,forcing->num_cells);
                        }
                            
                    }
                }
                else
                {
                    while (!feof(stormdata))
                    {
                        //Read intensity
                        result = fread(&cell, sizeof(unsigned int), 1, stormdata);
                        if (!result)	break;
                        fread(&intensity, sizeof(short unsigned int), 1, stormdata);

                        if (cell < forcing->num_cells)
                        {
                            if (forcing->received[cell])
                                printf("Warning: Received multiple intensities for cell %u in file %s.\n", cell, filename);
                            forcing->received[cell] = 1;
                            forcing->intensities[cell] = (float)(intensity * forcing->factor);
                        }
                        else{
                            printf("Warning: bad grid cell id in file %s.\n", filename);
                            printf("Read cell id %u num_cells is %u.\n", cell,forcing->num_cells);
                        }
                    }
                }

                fclose(stormdata);

                //Store 0's for remaining cells
                for (i = 0; i < forcing->num_cells; i++)
                    if (!forcing->received[i])	forcing->intensities[i] = 0.0;
            }
            else	//No file, no rain
            {
                for (i = 0; i < forcing->num_cells; i++)
                    forcing->intensities[i] = 0.0;
            }
        }

        MPI_Bcast(forcing->intensities, forcing->num_cells, MPI_FLOAT, 0, MPI_COMM_WORLD);

        //Load the data
        for (cell = 0; cell < forcing->num_cells; cell++)
        {
            for (i = 0; i < forcing->num_links_in_grid[cell]; i++)	//!!!! Assuming only links on this proc !!!!
            {
                curr_idx = forcing->grid_to_linkid[cell][i];
                assert(assignments[curr_idx] == my_rank);
                
                sys[curr_idx].my->forcing_data[forcing_idx].data[k].time = t_0 + k*increment;
                sys[curr_idx].my->forcing_data[forcing_idx].data[k].value = forcing->intensities[cell];
            }
        }
    }

    //Add in terms for no rainfall if max_files > numfiles
    for (i = 0; i < my_N; i++)
    {
        curr_idx = my_sys[i]->location;
        for (j = numfiles; j < max_files; j++)
        {
            sys[curr_idx].my->forcing_data[forcing_idx].data[j].time = sys[curr_idx].my->forcing_data[forcing_idx].data[j - 1].time + .0001;
            sys[curr_idx].my->forcing_data[forcing_idx].data[j].value = 0.0;
        }
    }

    //Add a ceiling term
    for (i = 0; i < my_N; i++)
    {
        curr_idx = my_sys[i]->location;
        //sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[max_files].time = GlobalVars->maxtime + 1.0;
        sys[curr_idx].my->forcing_data[forcing_idx].data[max_files].time = ceil_time;
        sys[curr_idx].my->forcing_data[forcing_idx].data[max_files].value = -1.0;
    }

    //Calculate the first rain change time and set rain_value
    for (i = 0; i < my_N; i++)
    {
        current = my_sys[i];
        forcing_buffer = current->my->forcing_data[forcing_idx].data[0].value;
        current->my->forcing_values[forcing_idx] = forcing_buffer;
        current->my->forcing_indices[forcing_idx] = 0;

        for (j = 1; j < current->my->forcing_data[forcing_idx].num_points; j++)
        {
            if (fabs(forcing_buffer - current->my->forcing_data[forcing_idx].data[j].value) > 1e-14)
            {
                current->my->forcing_change_times[forcing_idx] = current->my->forcing_data[forcing_idx].data[j].time;
                break;
            }
        }
        if (j == current->my->forcing_data[forcing_idx].num_points)
        {
            current->my->forcing_change_times[forcing_idx] = current->my->forcing_data[forcing_idx].data[j - 1].time;
        }
    }

    return 0;
}

//This reads in rainfall data at each link from an SQL database.
//Link** sys: An array of links.
//int N: The number of links in sys.
//int my_N: The number of links assigned to this process.
//UnivVars* GlobalVars: Contains all the information that is shared by every link in the system.
//int* my_sys: Array of links assigned to this process (value is location in sys array).
//int* assignments (set by this method): Will be an array with N entries. assignments[i] will the process link sys[i] is assigned to.
//char strfilename[]: String of the filename for the rain file to read (NOT .str files). Filenames should be indexed.
//unsigned int first: The index of the file to read first.
//unsigned int last: The index of the file to read last.
//double t_0: The time at which the first file starts.
//double increment: The amount of time between consecutively indexed files.
//int** id_to_loc (set by this method): Will be an array with N rows and 2 columns, sorted by first col. First col is a link id and second is
//				the location of the id in sys.
//unsigned int max_files: The maximum number of files to be read.
int Create_Rain_Database(
    Link *sys, unsigned int N,
    Link **my_sys, unsigned int my_N,
    const GlobalVars * const globals, int* assignments, ConnData *conninfo, unsigned int first, unsigned int last, Forcing* forcing, const Lookup * const id_to_loc, double maxtime, unsigned int forcing_idx)
{
#if defined(HAVE_POSTGRESQL)
    unsigned int i, j, k, curr_idx, tuple_count;
    Link* current;
    float forcing_buffer;
    int received_time;
    char* query = conninfo->query;
    PGresult *res;
    unsigned int *db_unix_time = NULL, *db_link_id = NULL;
    float *db_rain_intens = NULL;

    unsigned int *total_times = (unsigned int*)calloc(my_N, sizeof(unsigned int));
    for (i = 0; i < my_N; i++)
    {
        total_times[i] = my_sys[i]->my->forcing_data[forcing_idx].num_points;
        my_sys[i]->my->forcing_data[forcing_idx].num_points = 1;
    }

    //This is a time larger than any time in which the integrator is expected to get
    double ceil_time = 1e300;
    static short int gave_warning = 0;
    if (my_sys[0]->last_t > ceil_time*0.1 && !gave_warning)
    {
        gave_warning = 1;
        printf("[%i]: Warning: integrator time is extremely large (about %e). Loss of precision may occur.\n", my_rank, my_sys[i]->last_t);
    }

    /*
        //!!!! Fix this (or don't?) !!!!
        total_times = my_sys[0]->my->forcing_data[forcing_idx]->n_times;
        for(i=0;i<my_N;i++)	my_sys[i]->my->forcing_data[forcing_idx]->n_times = 1;
    */

    //GlobalVars->outletlink = 318213;

        //Query the database
    if (my_rank == 0)
    {
        //Connect to the database
        ConnectPGDB(conninfo);

        if (globals->outletlink == 0)
            sprintf(query, conninfo->queries[0], first, last);
        else
            sprintf(query, conninfo->queries[1], globals->outletlink, first, last);
        //printf("*************************\n");
        //printf("First = %u Last = %u t = %e increment = %u\n",first,last,my_sys[0]->last_t,forcing->increment);
        //printf("*************************\n");
        //printf("Gmax = %e maxtime = %e\n",GlobalVars->maxtime,maxtime);
        //printf("query: %s\n",query);
        //printf("*************************\n");
        res = PQexec(conninfo->conn, query);
        CheckResError(res, "downloading rainfall data");
        tuple_count = PQntuples(res);
        //printf("Received %u intensities.\n", tuple_count);
        
        //Disconnect
        DisconnectPGDB(conninfo);

        //Allocate space
        MPI_Bcast(&tuple_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
        db_unix_time = malloc(tuple_count * sizeof(unsigned int));
        db_rain_intens = malloc(tuple_count * sizeof(float));
        db_link_id = malloc(tuple_count * sizeof(unsigned int));

        //Load up the buffers
        for (i = 0; i < tuple_count; i++)
        {
            db_unix_time[i] = atoi(PQgetvalue(res, i, 0));
            db_rain_intens[i] = (float) atof(PQgetvalue(res, i, 1));
            db_link_id[i] = atoi(PQgetvalue(res, i, 2));
        }

        //Broadcast the data
        MPI_Bcast(db_unix_time, tuple_count, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(db_rain_intens, tuple_count, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Bcast(db_link_id, tuple_count, MPI_INT, 0, MPI_COMM_WORLD);

        //Clean up
        PQclear(res);
    }
    else
    {
        //Allocate space
        MPI_Bcast(&tuple_count, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (tuple_count)
        {
            db_unix_time = malloc(tuple_count * sizeof(unsigned int));
            db_rain_intens = malloc(tuple_count * sizeof(float));
            db_link_id = malloc(tuple_count * sizeof(unsigned int));

            //Receive the data
            MPI_Bcast(db_unix_time, tuple_count, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(db_rain_intens, tuple_count, MPI_FLOAT, 0, MPI_COMM_WORLD);
            MPI_Bcast(db_link_id, tuple_count, MPI_INT, 0, MPI_COMM_WORLD);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //double stop = time(NULL);
    //if(my_rank == 0)	printf("Total time to get rain: %f\n%u %u\n",difftime(stop,start),first,last);

        //Setup initial time in rainfall data with 0 rain
    for (i = 0; i < my_N; i++)
    {
        my_sys[i]->my->forcing_data[forcing_idx].data[0].time = (int)(first - forcing->raindb_start_time) / 60.0;
        my_sys[i]->my->forcing_data[forcing_idx].data[0].value = 0.0;
    }
    //printf("Started with first = %u, raindb_start = %u, %f\n",first,forcing->raindb_start_time,(int)(first - forcing->raindb_start_time)/60.0);
        //Setup the data received
    for (i = 0; i < tuple_count; i++)
    {
        //Find the location of ID i (ids should be numbered from 2 to whatever)
        //!!!! This will need to be a search at some point !!!!
/*
        if(GlobalVars->outletlink == 0)
        {
            curr_idx = id_to_loc[db_link_id[i]-2][1];
            //curr_idx = id_to_loc[atoi(PQgetvalue(res,i,2))-2][1];
            if(sys[curr_idx].ID != db_link_id[i])
                printf("Indices do not match %u %u\n",sys[curr_idx].ID,db_link_id[i]);
            //if(sys[curr_idx].ID != atoi(PQgetvalue(res,i,2)))
                //printf("Indices do not match %u %u\n",sys[curr_idx].ID,atoi(PQgetvalue(res,i,2)));
        }
        else
*/
        curr_idx = find_link_by_idtoloc(db_link_id[i], id_to_loc, N);
        //curr_idx = find_link_by_idtoloc(atoi(PQgetvalue(res,i,2)),id_to_loc,N);

        if (curr_idx < N && assignments[curr_idx] == my_rank)
        {
            k = sys[curr_idx].my->forcing_data[forcing_idx].num_points;
            received_time = db_unix_time[i] - forcing->raindb_start_time;	//In seconds
//printf("Got k = %i received_time = %i in secs = %i in mins = %f\n",k,received_time,(int)(sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k-1].time * 60.0 + 0.01),sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k-1].time);
            if (received_time > (int) (sys[curr_idx].my->forcing_data[forcing_idx].data[k - 1].time * 60.0 + 0.01))
            {
                if (received_time <= (int)(sys[curr_idx].my->forcing_data[forcing_idx].data[k - 1].time * 60.0) + (unsigned int)(forcing->file_time*60.0) || sys[curr_idx].my->forcing_data[forcing_idx].data[k - 1].value == 0.0)
                {
                    //printf("stored ID = %u i = %i k = %i received = %i unix_time = %i raindb_start = %i\n",sys[curr_idx].ID,i,k,received_time,db_unix_time[i],forcing->raindb_start_time);
                    sys[curr_idx].my->forcing_data[forcing_idx].data[k].time = received_time / 60.0;
                    sys[curr_idx].my->forcing_data[forcing_idx].data[k].value = db_rain_intens[i];
                    (sys[curr_idx].my->forcing_data[forcing_idx].num_points)++;
                    //printf("k = %i (%f %f)  Also: received = %i next = %u\n",k,sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k].time,sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k].value,received_time,(int) (sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k-1].time*60.0) + (unsigned int) (forcing->file_time*60.0));
                }
                else	//Add a 0 rainfall data
                {
                    //printf("stored 0 i = %i k = %i received = %i unix_time = %i raindb_start = %i\n",i,k,received_time,db_unix_time[i],forcing->raindb_start_time);
                    sys[curr_idx].my->forcing_data[forcing_idx].data[k].time = sys[curr_idx].my->forcing_data[forcing_idx].data[k - 1].time + forcing->file_time;
                    sys[curr_idx].my->forcing_data[forcing_idx].data[k].value = 0.0;
                    sys[curr_idx].my->forcing_data[forcing_idx].data[k + 1].time = received_time / 60.0;
                    sys[curr_idx].my->forcing_data[forcing_idx].data[k + 1].value = db_rain_intens[i];
                    sys[curr_idx].my->forcing_data[forcing_idx].num_points += 2;
                    //printf("(%f %f)\n",sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k].time,sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k].value);
                    //printf("(%f %f)\n",sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k+1].time,sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k+1].value);
                }
            }
            else if (received_time <= (int)(sys[curr_idx].my->forcing_data[forcing_idx].data[k - 1].time * 60.0 + 0.01))	//If the initial rate needs to be reset
            {
                //printf("stored init i = %i k = %i received = %i unix_time = %i raindb_start = %i\n",i,k,received_time,db_unix_time[i],forcing->raindb_start_time);
                sys[curr_idx].my->forcing_data[forcing_idx].data[k - 1].time = received_time / 60.0;
                sys[curr_idx].my->forcing_data[forcing_idx].data[k - 1].value = db_rain_intens[i];
                //printf("k = %i Init (%f %f)\n",k,sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k-1].time,sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k-1].value);
            }
            else
            {
                printf("!!!! Uh oh... !!!!\n");
                printf("!!!! i = %i k = %i received = %i unix_time = %i raindb_start = %i\n", i, k, received_time, db_unix_time[i], forcing->raindb_start_time);
            }
        }
    }
    //printf("Got %u\n",sys[0]->forcing_data[forcing_idx]->n_times);
        //Add ceiling terms
    for (i = 0; i < my_N; i++)
    {
        curr_idx = my_sys[i]->location;
        k = sys[curr_idx].my->forcing_data[forcing_idx].num_points;
        if (sys[curr_idx].my->forcing_data[forcing_idx].data[k - 1].value == 0.0)	//No rain, add just a ceiling
        {
            //sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k].time = maxtime * (1.1) + 1.0;
            sys[curr_idx].my->forcing_data[forcing_idx].data[k].time = ceil_time;
            sys[curr_idx].my->forcing_data[forcing_idx].data[k].value = -1.0;
        }
        else	//Add a 0.0, and a ceiling
        {
            sys[curr_idx].my->forcing_data[forcing_idx].data[k].time = sys[curr_idx].my->forcing_data[forcing_idx].data[k - 1].time + forcing->file_time;
            sys[curr_idx].my->forcing_data[forcing_idx].data[k].value = 0.0;
            //sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k+1].time = maxtime * (1.1) + 1.0;
            sys[curr_idx].my->forcing_data[forcing_idx].data[k + 1].time = ceil_time;
            sys[curr_idx].my->forcing_data[forcing_idx].data[k + 1].value = -1.0;
        }
    }

    //Reset n_times  !!!! Fix (well, this might be ok to do) !!!!
    //for(i=0;i<my_N;i++)	my_sys[i]->my->forcing_data[forcing_idx]->n_times = total_times;
    for (i = 0; i < my_N; i++)	my_sys[i]->my->forcing_data[forcing_idx].num_points = total_times[i];

    //Calculate the first rain change time and set rain_value
    for (i = 0; i < my_N; i++)
    {
        current = my_sys[i];

        //Get to the time block that corresponds to the current time, set forcing value
        for (j = 1; j < current->my->forcing_data[forcing_idx].num_points; j++)
        {
            if (current->last_t < current->my->forcing_data[forcing_idx].data[j].time - 1e-12)
                break;
        }

        forcing_buffer = current->my->forcing_data[forcing_idx].data[j - 1].value;
        current->my->forcing_values[forcing_idx] = forcing_buffer;
        current->my->forcing_indices[forcing_idx] = j - 1;

        for (; j < current->my->forcing_data[forcing_idx].num_points; j++)
        {
            if (fabs(forcing_buffer - current->my->forcing_data[forcing_idx].data[j].value) > 1e-12)
            {
                current->my->forcing_change_times[forcing_idx] = current->my->forcing_data[forcing_idx].data[j].time;
                break;
            }
        }
        if (j == current->my->forcing_data[forcing_idx].num_points)
            current->my->forcing_change_times[forcing_idx] = current->my->forcing_data[forcing_idx].data[j - 1].time;
    }

    //Clean up
    free(total_times);
    free(db_unix_time);
    free(db_link_id);
    free(db_rain_intens);

#else //HAVE_POSTGRESQL

    if (my_rank == 0)	printf("Error: Asynch was build without PostgreSQL support.\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

#endif //HAVE_POSTGRESQL

    return 0;
}


//This reads in rainfall data at each link from an SQL database. The timestamps are assumed to be irregularly spaced.
//Link** sys: An array of links.
//int N: The number of links in sys.
//int my_N: The number of links assigned to this process.
//UnivVars* GlobalVars: Contains all the information that is shared by every link in the system.
//int* my_sys: Array of links assigned to this process (value is location in sys array).
//int* assignments (set by this method): Will be an array with N entries. assignments[i] will the process link sys[i] is assigned to.
//char strfilename[]: String of the filename for the rain file to read (NOT .str files). Filenames should be indexed.
//unsigned int first: The index of the file to read first.
//unsigned int last: The index of the file to read last.
//double t_0: The time at which the first file starts.
//double increment: The amount of time between consecutively indexed files.
//int** id_to_loc (set by this method): Will be an array with N rows and 2 columns, sorted by first col. First col is a link id and second is
//				the location of the id in sys.
//unsigned int max_files: The maximum number of files to be read.
int Create_Rain_Database_Irregular(
    Link *sys, unsigned int N,
    Link **my_sys, unsigned int my_N,
    const GlobalVars * const globals, int* assignments, ConnData *conninfo, unsigned int first, unsigned int last, Forcing* forcing, const Lookup * const id_to_loc, double maxtime, unsigned int forcing_idx)
{
#if defined(HAVE_POSTGRESQL)
    unsigned int i, j, k, curr_idx, tuple_count, current_timestamp;
    Link* current;
    float forcing_buffer;
    int received_time;
    char* query = conninfo->query;
    PGresult *res;
    unsigned int *db_unix_time, *db_link_id;
    float *db_rain_intens;
    unsigned int *actual_timestamps, num_actual_timestamps, max_timestamps = forcing->increment;	//The maximum number of intensities to get for each link

/*
    unsigned int *total_times = (unsigned int*) calloc(my_N,sizeof(unsigned int));
    for(i=0;i<my_N;i++)
    {
        if(my_sys[i]->my->forcing_data[forcing_idx])
        {
            total_times[i] = my_sys[i]->my->forcing_data[forcing_idx]->n_times;
            my_sys[i]->my->forcing_data[forcing_idx]->n_times = 1;
        }
    }
*/

//This is a time larger than any time in which the integrator is expected to get
    double ceil_time = 1e300;
    static short int gave_warning = 0;
    if (my_sys[0]->last_t > ceil_time*0.1 && !gave_warning)
    {
        gave_warning = 1;
        printf("[%i]: Warning: integrator time is extremely large (about %e). Loss of precision may occur.\n", my_rank, my_sys[0]->last_t);
    }

    //Query the database
    if (my_rank == 0)
    {
        //Connect to the database
        ConnectPGDB(conninfo);

        //Download timestamps
        sprintf(query, conninfo->queries[3], first, max_timestamps);
        res = PQexec(conninfo->conn, query);
        CheckResError(res, "downloading rainfall timestamps");
        num_actual_timestamps = PQntuples(res);

        //Unpack and broadcast the timestamps
        actual_timestamps = (unsigned int*)malloc(num_actual_timestamps * sizeof(unsigned int));
        for (i = 0; i < num_actual_timestamps; i++)
            actual_timestamps[i] = atoi(PQgetvalue(res, i, 0));
        MPI_Bcast(&num_actual_timestamps, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        MPI_Bcast(actual_timestamps, num_actual_timestamps, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        last = actual_timestamps[num_actual_timestamps - 1];
        PQclear(res);

        //Download intensities
        if (globals->outletlink == 0)
            sprintf(query, conninfo->queries[0], first, last);
        else
            sprintf(query, conninfo->queries[1], globals->outletlink, first, last);
        //printf("*************************\n");
        //printf("First = %u Last = %u t = %e increment = %u\n",first,last,my_sys[0]->last_t,forcing->increment);
        //printf("*************************\n");
        //printf("Gmax = %e maxtime = %e\n",GlobalVars->maxtime,maxtime);
        //printf("query: %s\n",query);
        //printf("*************************\n");
        res = PQexec(conninfo->conn, query);
        CheckResError(res, "downloading rainfall data");
        tuple_count = PQntuples(res);
        //printf("Received %u intensities.\n", tuple_count);
        
        //Disconnect
        DisconnectPGDB(conninfo);

        //Allocate space
        MPI_Bcast(&tuple_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
        db_unix_time = malloc(tuple_count * sizeof(unsigned int));
        db_rain_intens = malloc(tuple_count * sizeof(float));
        db_link_id = malloc(tuple_count * sizeof(unsigned int));

        //Load up the buffers
        for (i = 0; i < tuple_count; i++)
        {
            db_unix_time[i] = atoi(PQgetvalue(res, i, 0));
            db_rain_intens[i] = (float) atof(PQgetvalue(res, i, 1));
            db_link_id[i] = atoi(PQgetvalue(res, i, 2));
        }

        //Broadcast the data
        MPI_Bcast(db_unix_time, tuple_count, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(db_rain_intens, tuple_count, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Bcast(db_link_id, tuple_count, MPI_INT, 0, MPI_COMM_WORLD);

        //Clean up
        PQclear(res);
    }
    else
    {
        MPI_Bcast(&num_actual_timestamps, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        actual_timestamps = (unsigned int*)malloc(num_actual_timestamps * sizeof(unsigned int));
        MPI_Bcast(actual_timestamps, num_actual_timestamps, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        last = actual_timestamps[num_actual_timestamps - 1];

        //Allocate space
        MPI_Bcast(&tuple_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
        db_unix_time = malloc(tuple_count * sizeof(unsigned int));
        db_rain_intens = malloc(tuple_count * sizeof(float));
        db_link_id = malloc(tuple_count * sizeof(unsigned int));

        //Receive the data
        MPI_Bcast(db_unix_time, tuple_count, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(db_rain_intens, tuple_count, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Bcast(db_link_id, tuple_count, MPI_INT, 0, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //double stop = time(NULL);
    //if(my_rank == 0)	printf("Total time to get rain: %f\n%u %u\n",difftime(stop,start),first,last);

    /*
    printf("+++++++++\n");
    for(i=0;i<num_actual_timestamps;i++)
    {
    printf("%u\n",actual_timestamps[i]);
    }
    printf("+++++++++\n");
    */

    //Set the times and zero out the intensities (some of these might change later)
    for (i = 0; i < my_N; i++)
    {
        for (j = 0; j < num_actual_timestamps; j++)
        {
            my_sys[i]->my->forcing_data[forcing_idx].data[j].time = (double)(actual_timestamps[j] - forcing->raindb_start_time) / 60.0;
            my_sys[i]->my->forcing_data[forcing_idx].data[j].value = 0.0;
        }
    }

    //Setup the forcing values
    i = 0;
    for (j = 0; j < num_actual_timestamps; j++)
    {
        current_timestamp = actual_timestamps[j];
        for (; i < tuple_count && db_unix_time[i] <= current_timestamp; i++)
        {
            received_time = db_unix_time[i] - forcing->raindb_start_time;	//In seconds
            curr_idx = find_link_by_idtoloc(db_link_id[i], id_to_loc, N);
            if (curr_idx < N && assignments[curr_idx] == my_rank)
            {
                sys[curr_idx].my->forcing_data[forcing_idx].data[j].time = received_time / 60.0;
                sys[curr_idx].my->forcing_data[forcing_idx].data[j].value = db_rain_intens[i];
                //if(sys[curr_idx].ID == 456117)
                //printf("j = %u received_time = %u intensity = %f\n",j,received_time,sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[j].value);
            }
            //else
            //printf("!!!! Huh, wtf? !!!!\n");
        }
    }

    /*
        current_timestamp = actual_timestamps[0];
        for(i=0;i < tuple_count && current_timestamp <= db_unix_time[i];i++)
        {
            received_time = db_unix_time[i] - forcing->raindb_start_time;	//In seconds
            //if(received_time > (int) (sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k-1].time * 60.0 + 0.01))	break;
            curr_idx = find_link_by_idtoloc(db_link_id[i],id_to_loc,N);
            if(curr_idx < N && assignments[curr_idx] == my_rank)
            {
                //k = sys[curr_idx].my->forcing_data[forcing_idx]->n_times;
                //if(received_time <= (int) (sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k-1].time * 60.0 + 0.01))	//If the initial rate needs to be reset
                {
                    sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[0].time = received_time / 60.0;
                    sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[0].value = db_rain_intens[i];
                }
            }
        }

        //!!!! Uh, any difference between this and j=0? !!!!
        for(j=1;j<num_actual_timestamps;j++)
        {
            current_timestamp = actual_timestamps[j];
            for(;i < tuple_count && current_timestamp <= db_unix_time[i];i++)
            {
                received_time = db_unix_time[i] - forcing->raindb_start_time;	//In seconds
                curr_idx = find_link_by_idtoloc(db_link_id[i],id_to_loc,N);

                if(curr_idx < N && assignments[curr_idx] == my_rank)
                {
                    sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[j].time = received_time / 60.0;
                    sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[j].value = db_rain_intens[i];
                }
            }
        }
    */

    /*
        //Setup the data received
        for(i=0;i<tuple_count;i++)
        {
            //Find the location of ID i (ids should be numbered from 2 to whatever)
            curr_idx = find_link_by_idtoloc(db_link_id[i],id_to_loc,N);

            if(curr_idx < N && assignments[curr_idx] == my_rank)
            {
                k = sys[curr_idx].my->forcing_data[forcing_idx]->n_times;
                received_time = db_unix_time[i] - forcing->raindb_start_time;	//In seconds

                if(received_time > (int) (sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k-1].time * 60.0 + 0.01))
                {
                    if( received_time <= (int) (sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k-1].time*60.0) + (unsigned int) (forcing->file_time*60.0) || sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k-1].value == 0.0)
                    {
                        sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k].time = received_time / 60.0;
                        sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k].value = db_rain_intens[i];
                        (sys[curr_idx].my->forcing_data[forcing_idx]->n_times)++;
                    }
                    else	//Add a 0 rainfall data before adding this data
                    {
                        sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k].time = sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k-1].time + forcing->file_time;
                        sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k].value = 0.0;
                        sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k+1].time = received_time / 60.0;
                        sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k+1].value = db_rain_intens[i];
                        sys[curr_idx].my->forcing_data[forcing_idx]->n_times += 2;
                    }
                }
                else //if(received_time <= (int) (sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k-1].time * 60.0 + 0.01))	//If the initial rate needs to be reset
                {
                    sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k-1].time = received_time / 60.0;
                    sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k-1].value = db_rain_intens[i];
                }
            }
        }
    */
    //Add ceiling terms
    for (i = 0; i < my_N; i++)
    {
        curr_idx = my_sys[i]->location;
        //k = sys[curr_idx].my->forcing_data[forcing_idx]->n_times;
        k = num_actual_timestamps;
        //if(sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k-1].value == 0.0)	//No rain, add just a ceiling
        {
            sys[curr_idx].my->forcing_data[forcing_idx].data[k].time = ceil_time;
            sys[curr_idx].my->forcing_data[forcing_idx].data[k].value = -1.0;
        }
        /*
                else	//Add a 0.0, and a ceiling
                {
                    sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k].time = sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k-1].time + forcing->file_time;
                    sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k].value = 0.0;
                    //sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k+1].time = maxtime * (1.1) + 1.0;
                    sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k+1].time = ceil_time;
                    sys[curr_idx].my->forcing_data[forcing_idx]->rainfall[k+1].value = -1.0;
                }
        */
    }

    //Reset n_times  !!!! Fix (well, this might be ok to do) !!!!
    //for(i=0;i<my_N;i++)	my_sys[i]->my->forcing_data[forcing_idx]->n_times = total_times[i];

    //Calculate the first rain change time and set rain_value
    for (i = 0; i < my_N; i++)
    {
        current = my_sys[i];

        //Get to the time block that corresponds to the current time, set forcing value
        for (j = 1; j < current->my->forcing_data[forcing_idx].num_points; j++)
        {
            if (current->last_t < current->my->forcing_data[forcing_idx].data[j].time - 1e-12)
                break;
        }

        forcing_buffer = current->my->forcing_data[forcing_idx].data[j - 1].value;
        current->my->forcing_values[forcing_idx] = forcing_buffer;
        current->my->forcing_indices[forcing_idx] = j - 1;

        for (; j < current->my->forcing_data[forcing_idx].num_points; j++)
        {
            if (fabs(forcing_buffer - current->my->forcing_data[forcing_idx].data[j].value) > 1e-12)
            {
                current->my->forcing_change_times[forcing_idx] = current->my->forcing_data[forcing_idx].data[j].time;
                break;
            }
        }
        if (j == current->my->forcing_data[forcing_idx].num_points)
            current->my->forcing_change_times[forcing_idx] = current->my->forcing_data[forcing_idx].data[j - 1].time;
    }

    //Clean up
    //free(total_times);
    free(actual_timestamps);
    free(db_unix_time);
    free(db_link_id);
    free(db_rain_intens);

    return last;

#else //HAVE_POSTGRESQL

    if (my_rank == 0)	printf("Error: Asynch was build without PostgreSQL support.\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

    return 0;

#endif //HAVE_POSTGRESQL
}


//Sets a forcing based on monthly data.
//Assumes the rates are already set.
double CreateForcing_Monthly(
    Link *sys, unsigned int N,
    Link **my_sys, unsigned int my_N,
    const GlobalVars * const globals, TimeSerie* global_forcings, unsigned int forcing_idx, struct tm *current_time, time_t first_time, time_t last_time, double t_0)
{
    unsigned int num_months = 12;
    int i;
    int month_0, current_year;
    char buffer[4];
    double t = t_0;

    //Find the current month
    strftime(buffer, 4, "%m", current_time);
    month_0 = atoi(buffer) - 1;
    if (month_0 < 0 || month_0 > num_months)
    {
        printf("[%i]: Error: Bad month %i.\n", my_rank, month_0);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    //Set the (local) times for the current month and previous months
    global_forcings->data[month_0].time = t_0;
    for (i = month_0 - 1; i > -1; i--)
        global_forcings->data[i].time = global_forcings->data[i + 1].time - 1.0;
    current_year = current_time->tm_year + 1900;

    //Find days until next month
    int days_this_month = days_in_month(month_0, current_year);
    int curr_day = current_time->tm_mday;
    t += (days_this_month - curr_day) * (24.0*60.0);	//Convert to mins

    //Add in hours and mins and seconds
    t += (double)((23 - current_time->tm_hour) * 60 + (59 - current_time->tm_min)) + (60.0 - current_time->tm_sec) / 60.0;

    //Save t
    global_forcings->data[month_0 + 1].time = t;

    //Set the (local) times for each future month
    for (i = month_0 + 2; i <= num_months; i++)	//This should set the ceiling term
    {
        t += days_in_month(i - 1, current_year) * (24.0*60.0);
        global_forcings->data[i].time = t;
    }

    //Set ceiling term
    global_forcings->data[num_months].value = 0.0;

    //Check if this data goes past the last time
    double final_time = (last_time - first_time) / 60.0;
    if (t > final_time)	//Need to set 0s past final_time
    {
        i = num_months - 1;
        global_forcings->data[i].time = t + 1.0;
        global_forcings->data[i].value = 0.0;

        for (i -= 1; global_forcings->data[i].time > final_time; i--)
        {
            global_forcings->data[i].time = final_time;
            global_forcings->data[i].value = 0.0;
            global_forcings->data[i + 1].time = t + 1.0;
            global_forcings->data[i + 1].value = 0.0;
        }
    }

    //Set the current forcing value at each link
    for (i = 0; i < my_N; i++)
    {
        Link *current = my_sys[i];
        current->my->forcing_values[forcing_idx] = global_forcings->data[month_0].value;
        current->my->forcing_indices[forcing_idx] = month_0;
        current->my->forcing_change_times[forcing_idx] = global_forcings->data[month_0 + 1].time;
    }

    return t;
}

