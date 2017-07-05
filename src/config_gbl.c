#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(HAVE_UNISTD_H)
#include <unistd.h>
#endif

#if defined(_MSC_VER)
#include <process.h>
#endif

#include <structs.h>
#include <globals.h>
#include <config_gbl.h>
#include <db.h>
#include <outputs.h>
#include <io.h>
#include <models/definitions.h>


#if defined(_MSC_VER)
#   define timegm _mkgmtime
#endif


/// Put a vector of global_params onto the end of a filename.
int AttachParameters(
    char* filename,
    unsigned int max_size,
    double *global_params, unsigned int num_global_params,
    unsigned int string_size);

int CheckWinFormat(FILE* file);
int CheckFilenameExtension(char *filename, char *extension);

void ReadLineFromTextFile(FILE* globalfile, char *line_buffer, unsigned int size);

int ReadLineError(int valsread, int valswant, char *message);

void ReadDBC(char* filename, ConnData* const conninfo);

int RemoveSuffix(char *filename, const char *suffix);


GlobalVars* Read_Global_Data(
    char *globalfilename,
    ErrorData *errors,
    Forcing *forcings,
    ConnData *db_connections,
    char *rkdfilename,
    AsynchModel const *model,
    void *external)
{
    unsigned int i, total, written;
    int flag, valsread;
    char endmark;
    char line_buffer[ASYNCH_MAX_LINE_LENGTH];
    const unsigned int line_buffer_len = ASYNCH_MAX_LINE_LENGTH;

    GlobalVars* globals = (GlobalVars*)malloc(sizeof(GlobalVars));
    memset(globals, 0, sizeof(GlobalVars));

    globals->string_size = ASYNCH_MAX_PATH_LENGTH;
    char db_filename[ASYNCH_MAX_PATH_LENGTH];
    globals->query_size = ASYNCH_MAX_QUERY_LENGTH;

    FILE* globalfile = NULL;

    if (my_rank == 0)
    {
        globalfile = fopen(globalfilename, "r");
        if (globalfile == NULL)
        {
            printf("Error: Global file %s was not found.\n", globalfilename);
            return NULL;
        }

        if (CheckWinFormat(globalfile))
        {
            printf("Error: File %s appears to be in Windows format. Try converting to unix format using 'dos2unix' at the command line.\n", globalfilename);
            return NULL;
        }
    }

    globals->rain_filename = NULL;

    //Grab the model uid
    unsigned short uid;
    ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
    valsread = sscanf(line_buffer, "%hu", &uid);
    if (ReadLineError(valsread, 1, "model type"))	return NULL;

    globals->model_uid = uid;

    //model = GetModel(uid);
    //if (model == NULL)
    //{
    //    printf("Error: Model %hu not defined.\n", uid);
    //    return NULL;
    //}

    //Grab the begin and end time
    struct tm begin_tm;
    memset(&begin_tm, 0, sizeof(struct tm));
    
    ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
    valsread = sscanf(line_buffer, "%d-%d-%d %d:%d", &begin_tm.tm_year, &begin_tm.tm_mon, &begin_tm.tm_mday, &begin_tm.tm_hour, &begin_tm.tm_min);
    if (valsread == 5) 
    {
        begin_tm.tm_year = begin_tm.tm_year - 1900;
        begin_tm.tm_mon = begin_tm.tm_mon - 1;
        globals->begin_time = timegm(&begin_tm);
    }
    else
    {
        int begin_time;
        valsread = sscanf(line_buffer, "%d", &begin_time);
        if (ReadLineError(valsread, 1, "begin YYYY-MM-DD HH:MM || unix_time"))	return NULL;
        globals->begin_time = begin_time;
    }
        
    struct tm end_tm;
    memset(&end_tm, 0, sizeof(struct tm));

    ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
    valsread = sscanf(line_buffer, "%d-%d-%d %d:%d", &end_tm.tm_year, &end_tm.tm_mon, &end_tm.tm_mday, &end_tm.tm_hour, &end_tm.tm_min);
    if (valsread == 5)
    {
        end_tm.tm_year = end_tm.tm_year - 1900;
        end_tm.tm_mon = end_tm.tm_mon - 1;
        globals->end_time = timegm(&end_tm);
    }
    else
    {
        int end_time;
        valsread = sscanf(line_buffer, "%d", &end_time);
        if (ReadLineError(valsread, 1, "end YYYY-MM-DD HH:MM || unix_time"))	return NULL;
        globals->end_time = end_time;
    }

    globals->maxtime = (double)(globals->end_time - globals->begin_time) / 60.0;
    if (globals->maxtime <= 0.0)
    {
        printf("Error: Simulation period invalid (begin >= end)\n");
        return NULL;
    }

    //Grab the output filename info
    ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
    valsread = sscanf(line_buffer, "%hu", &(globals->print_par_flag));
    if (ReadLineError(valsread, 1, "to print filename parameters"))	return NULL;

    //Grab components to print
    ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
    valsread = sscanf(line_buffer, "%u", &globals->num_outputs);
    if (ReadLineError(valsread, 1, "number of indices to print"))	return NULL;
    globals->outputs = malloc(globals->num_outputs * sizeof(Output));
    for (i = 0; i < globals->num_outputs; i++)
    {
        globals->outputs[i].name = (char*)malloc(ASYNCH_MAX_SYMBOL_LENGTH * sizeof(char));
        ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
        valsread = sscanf(line_buffer, "%s", globals->outputs[i].name);
        if (ReadLineError(valsread, 1, "a component to print"))	return NULL;
    }

    //Peakflow function
    globals->peakflow_function_name = (char*)malloc(ASYNCH_MAX_SYMBOL_LENGTH * sizeof(char));
    ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
    valsread = sscanf(line_buffer, "%s", globals->peakflow_function_name);
    if (ReadLineError(valsread, 1, "peakflow function name"))	return NULL;
    SetPeakflowOutputFunctions(globals->peakflow_function_name, &(globals->peakflow_output));

    //Grab the global parameters
    ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
    valsread = sscanf(line_buffer, "%u%n", &globals->num_global_params, &total);
    if (ReadLineError(valsread, 1, "number of global parameters"))	return NULL;
    globals->global_params = malloc(globals->num_global_params * sizeof(double));
    //if (my_rank == 0 && globals->num_global_params != model->num_global_params)
    //{
    //    printf("[%i]: Error: Got %u global params in the .gbl file. Expected %u for model %u.\n", my_rank, globals->num_global_params, model->num_global_params, globals->model_uid);
    //    MPI_Abort(MPI_COMM_WORLD, 1);
    //}
    for (i = 0; i < globals->num_global_params; i++)
    {
        valsread = sscanf(&(line_buffer[total]), "%lf%n", &globals->global_params[i], &written);
        if (ReadLineError(valsread, 1, "a global parameter"))	return NULL;
        total += written;
    }

    //Set dim and other sizes
    if (model)
        model->set_param_sizes(globals, external);
    else
        SetParamSizes(globals, external);

    //Find the states needed for printing
    globals->num_states_for_printing = 0;
    globals->print_indices = (unsigned int*)calloc(globals->num_outputs, sizeof(unsigned int));
    for (i = 0; i < globals->num_outputs; i++)
        SetDefaultOutputFunctions(globals->outputs[i].name, &globals->outputs[i], globals->print_indices, &globals->num_states_for_printing);

    globals->print_indices = (unsigned int*)realloc(globals->print_indices, globals->num_states_for_printing * sizeof(unsigned int));

    //Grab the stored steps limits
    ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
    valsread = sscanf(line_buffer, "%u %i %u", &(globals->iter_limit), &(globals->max_transfer_steps), &(globals->discont_size));
    if (ReadLineError(valsread, 3, "steps stored, steps transfered, and discontinuity buffer size"))	return NULL;

    //Grab the topology data filename
    globals->outletlink = 0;
    ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
    valsread = sscanf(line_buffer, "%hu", &(globals->rvr_flag));
    if (ReadLineError(valsread, 1, "topology data flag"))	return NULL;
    if (globals->rvr_flag == 0)
    {
        globals->rvr_filename = (char*)malloc(ASYNCH_MAX_PATH_LENGTH * sizeof(char));
        valsread = sscanf(line_buffer, "%*u %s", globals->rvr_filename);
        if (ReadLineError(valsread, 1, "filename for topology data"))	return NULL;
        if (!CheckFilenameExtension(globals->rvr_filename, ".rvr"))	return NULL;
    }
    else	//Reading from database
    {
        valsread = sscanf(line_buffer, "%*u %u %s", &(globals->outletlink), db_filename);
        if (ReadLineError(valsread, 2, "link id of downstream link for topology data or .dbc for topology"))	return NULL;
        if (!CheckFilenameExtension(db_filename, ".dbc"))	return NULL;
        globals->rvr_filename = NULL;
        ReadDBC(db_filename, &db_connections[ASYNCH_DB_LOC_TOPO]);
    }

    //Grab the parameter data filename
    ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
    valsread = sscanf(line_buffer, "%hu", &(globals->prm_flag));
    if (ReadLineError(valsread, 1, "parameter flag"))	return NULL;
    if (globals->prm_flag == 0)
    {
        globals->prm_filename = (char*)malloc(ASYNCH_MAX_PATH_LENGTH * sizeof(char));
        valsread = sscanf(line_buffer, "%*u %s", globals->prm_filename);
        if (ReadLineError(valsread, 1, ".prm filename"))	return NULL;
        if (!CheckFilenameExtension(globals->prm_filename, ".prm"))	return NULL;
    }
    else
    {
        valsread = sscanf(line_buffer, "%*u %s", db_filename);
        if (ReadLineError(valsread, 1, ".dbc for parameters"))	return NULL;
        if (!CheckFilenameExtension(db_filename, ".dbc"))	return NULL;
        globals->prm_filename = NULL;
        ReadDBC(db_filename, &db_connections[ASYNCH_DB_LOC_PARAMS]);
    }

    //Grab the initial data file
    globals->init_filename = NULL;
    ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
    valsread = sscanf(line_buffer, "%hu", &(globals->init_flag));
    if (ReadLineError(valsread, 1, "initial data flag"))	return NULL;
    if (globals->init_flag == 0 || globals->init_flag == 1 || globals->init_flag == 2 || globals->init_flag == 4)
    {
        globals->init_filename = (char*)malloc(ASYNCH_MAX_PATH_LENGTH * sizeof(char));
        valsread = sscanf(line_buffer, "%*u %s", globals->init_filename);
        if (ReadLineError(valsread, 1, "initial data flag"))	return NULL;
        if (globals->init_flag == 0 && !CheckFilenameExtension(globals->init_filename, ".ini")) return NULL;
        if (globals->init_flag == 1 && !CheckFilenameExtension(globals->init_filename, ".uini")) return NULL;
        if (globals->init_flag == 2 && !CheckFilenameExtension(globals->init_filename, ".rec")) return NULL;
        if (globals->init_flag == 4 && !CheckFilenameExtension(globals->init_filename, ".h5")) return NULL;
    }
    else if (globals->init_flag == 3)
    {
        valsread = sscanf(line_buffer, "%*u %s %u", db_filename, &(globals->init_timestamp));
        if (ReadLineError(valsread, 1, ".dbc for parameters"))	return NULL;
        if (!CheckFilenameExtension(db_filename, ".dbc"))	return NULL;
        ReadDBC(db_filename, &db_connections[ASYNCH_DB_LOC_INIT]);
    }

    //Grab number of forcings
    unsigned int got_forcings;
    ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
    valsread = sscanf(line_buffer, "%u", &got_forcings);
    if (ReadLineError(valsread, 1, "rainfall flag"))	return NULL;
    if (got_forcings < globals->num_forcings && my_rank == 0)
    {
        printf("[%i]: Error: Got %u forcings in the .gbl file. Expected %u for model %u.\n", my_rank, got_forcings, globals->num_forcings, globals->model_uid);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (got_forcings > globals->num_forcings && my_rank == 0)
    {
        printf("[%i]: Warning: Got %u forcings in the .gbl file. Expected %u for model %u.\n", my_rank, got_forcings, globals->num_forcings, globals->model_uid);
        globals->num_forcings = got_forcings;
    }

    //Grab the forcing parameters
    //0 for no rain, 1 for .str file, 2 for binary files, 3 for database, 4 for uniform rain (.ustr)
    globals->hydro_table = globals->peak_table = NULL;
    for (i = 0; i < globals->num_forcings; i++)
    {
        ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
        valsread = sscanf(line_buffer, "%hi", &(forcings[i].flag));
        if (ReadLineError(valsread, 1, "forcings flag"))	return NULL;

        if (forcings[i].flag == 1 || forcings[i].flag == 2 || forcings[i].flag == 4 || forcings[i].flag == 6 || forcings[i].flag == 8)
        {
            forcings[i].filename = (char*)malloc(ASYNCH_MAX_PATH_LENGTH * sizeof(char));
            valsread = sscanf(line_buffer, "%*i %s", forcings[i].filename);
            if (ReadLineError(valsread, 1, "forcing data filename"))	return NULL;
            if (forcings[i].flag == 1 && !CheckFilenameExtension(forcings[i].filename, ".str"))	return NULL;
            if (forcings[i].flag == 4 && !CheckFilenameExtension(forcings[i].filename, ".ustr"))	return NULL;

            if (forcings[i].flag == 2 || forcings[i].flag == 6)
            {
                ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
                valsread = sscanf(line_buffer, "%u %lf %u %u", &(forcings[i].increment), &(forcings[i].file_time), &(forcings[i].first_file), &(forcings[i].last_file));
                if (ReadLineError(valsread, 4, "time increment, file time, first file, and last file"))	return NULL;
            }
            else if (forcings[i].flag == 8)
            {
                ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
                valsread = sscanf(line_buffer, "%u %u %u", &(forcings[i].increment), &(forcings[i].first_file), &(forcings[i].last_file));
                if (ReadLineError(valsread, 3, "time increment, first file, and last file"))	return NULL;
            }
        }
        else if (forcings[i].flag == 3)	//Database
        {
            valsread = sscanf(line_buffer, "%*i %s", db_filename);
            if (ReadLineError(valsread, 1, ".dbc for forcing"))	return NULL;
            if (!CheckFilenameExtension(db_filename, ".dbc"))	return NULL;
            ReadDBC(db_filename, &db_connections[ASYNCH_DB_LOC_FORCING_START + i]);

            ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
            valsread = sscanf(line_buffer, "%u %lf %u %u", &(forcings[i].increment), &(forcings[i].file_time), &(forcings[i].first_file), &(forcings[i].last_file));
            if (ReadLineError(valsread, 4, "chunk size, time resolution, beginning unix time, and ending unix time"))	return NULL;
            forcings[i].raindb_start_time = forcings[i].first_file;
        }
        else if (forcings[i].flag == 9)	//Database with irregular timesteps
        {
            valsread = sscanf(line_buffer, "%*i %s", db_filename);
            if (ReadLineError(valsread, 1, ".dbc for forcing"))	return NULL;
            if (!CheckFilenameExtension(db_filename, ".dbc"))	return NULL;
            ReadDBC(db_filename, &db_connections[ASYNCH_DB_LOC_FORCING_START + i]);

            ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
            valsread = sscanf(line_buffer, "%u %u %u", &(forcings[i].increment), &(forcings[i].first_file), &(forcings[i].last_file));
            if (ReadLineError(valsread, 3, "time increment, file time, first file, and last file"))	return NULL;
            forcings[i].raindb_start_time = forcings[i].first_file;

            forcings[i].lastused_first_file = forcings[i].lastused_last_file = 0;
            forcings[i].next_timestamp = forcings[i].first_file;
            forcings[i].number_timesteps = 0;
        }
        else if (forcings[i].flag == 7)	//Recurring
        {
            forcings[i].filename = (char*)malloc(ASYNCH_MAX_PATH_LENGTH * sizeof(char));
            valsread = sscanf(line_buffer, "%*i %s", forcings[i].filename);
            if (ReadLineError(valsread, 1, "recurring rainfall filename"))	return NULL;
            if (!CheckFilenameExtension(forcings[i].filename, ".mon"))	return NULL;

            ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
            valsread = sscanf(line_buffer, "%u %u", &(forcings[i].first_file), &(forcings[i].last_file));
            if (ReadLineError(valsread, 2, "first time, and last time"))	return NULL;
            forcings[i].raindb_start_time = forcings[i].first_file;
        }
        else if (forcings[i].flag == 0) //No forcing
        {
            forcings[i].filename = NULL;
        }
        else
        {
            printf("[%i]: Error reading %s: Invalid forcing flag %i.\n", my_rank, globalfilename, forcings[i].flag);
            return NULL;
        }
    }

    //Grab the dam filename
    ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
    valsread = sscanf(line_buffer, "%hu", &(globals->dam_flag));
    if (ReadLineError(valsread, 1, "dam flag"))	return NULL;

    if (globals->dam_flag)
    {
        globals->dam_filename = (char*)malloc(ASYNCH_MAX_PATH_LENGTH * sizeof(char));
        valsread = sscanf(line_buffer, "%*u %s", globals->dam_filename);
        if (ReadLineError(valsread, 1, "filename for dam info"))	return NULL;
        if (globals->dam_flag == 1 && !CheckFilenameExtension(globals->dam_filename, ".dam"))	return NULL;
        if (globals->dam_flag == 2 && !CheckFilenameExtension(globals->dam_filename, ".qvs"))	return NULL;

        if (globals->dam_flag == 3)
        {
            if (!CheckFilenameExtension(globals->dam_filename, ".dbc"))	return NULL;
            ReadDBC(globals->dam_filename, &db_connections[ASYNCH_DB_LOC_QVS]);
        }
    }
    else
        globals->dam_filename = NULL;

    //Get the link ids where reservoirs exist
    ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
    valsread = sscanf(line_buffer, "%hu", &(globals->res_flag));
    if (ReadLineError(valsread, 1, "res flag"))	return NULL;

    if (globals->res_flag)
    {
        if (globals->res_flag == 1)
        {
            globals->rsv_filename = (char*)malloc(ASYNCH_MAX_PATH_LENGTH * sizeof(char));
            valsread = sscanf(line_buffer, "%*u %s %hi", globals->rsv_filename, &(globals->res_forcing_idx));
            if (ReadLineError(valsread, 2, ".rsv filename"))	return NULL;
            if (!CheckFilenameExtension(globals->rsv_filename, ".rsv"))	return NULL;
        }
        else	//Flag is 2
        {
            globals->rsv_filename = (char*)malloc(ASYNCH_MAX_PATH_LENGTH * sizeof(char));
            valsread = sscanf(line_buffer, "%*u %s %hi", globals->rsv_filename, &(globals->res_forcing_idx));
            if (ReadLineError(valsread, 2, ".dbc for reservoirs"))	return NULL;
            if (!CheckFilenameExtension(globals->rsv_filename, ".dbc"))	return NULL;
            ReadDBC(globals->rsv_filename, &db_connections[ASYNCH_DB_LOC_RSV]);
        }

        if (globals->res_forcing_idx >= globals->num_forcings)
        {
            if (my_rank == 0)	printf("Bad forcing index for a reservoir feed (%hi). Only %i forcings available.\n", globals->res_forcing_idx, globals->num_forcings);
            return NULL;
        }
    }
    else
    {
        globals->rsv_filename = NULL;
        globals->res_forcing_idx = -1;
    }

    //Grab where to write the timeseries
    ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
    valsread = sscanf(line_buffer, "%hu", &(globals->hydros_loc_flag));
    if (ReadLineError(valsread, 1, "timeseries location"))	return NULL;

    globals->hydros_loc_filename = NULL;
    globals->hydro_table = NULL;

    if (globals->hydros_loc_flag == 1 || globals->hydros_loc_flag == 2 || globals->hydros_loc_flag == 4 || globals->hydros_loc_flag == 5 || globals->hydros_loc_flag == 6)
    {
        globals->hydros_loc_filename = (char*)malloc(ASYNCH_MAX_PATH_LENGTH * sizeof(char));
        valsread = sscanf(line_buffer, "%*u %lf %s", &(globals->print_time), globals->hydros_loc_filename);
        if (ReadLineError(valsread, 2, "hydrographs location"))	return NULL;
        if (globals->hydros_loc_flag == 1 && !CheckFilenameExtension(globals->hydros_loc_filename, ".dat"))	return NULL;
        if (globals->hydros_loc_flag == 2 && !CheckFilenameExtension(globals->hydros_loc_filename, ".csv"))	return NULL;
        if (globals->hydros_loc_flag == 4 && !CheckFilenameExtension(globals->hydros_loc_filename, ".rad"))	return NULL;
        if (globals->hydros_loc_flag == 5 && !CheckFilenameExtension(globals->hydros_loc_filename, ".h5"))	return NULL;
        if (globals->hydros_loc_flag == 6 && !CheckFilenameExtension(globals->hydros_loc_filename, ".h5"))	return NULL;
        //globals->output_flag = (globals->hydros_loc_flag == 1) ? 0 : 1;

        if (globals->hydros_loc_flag == 1)	RemoveSuffix(globals->hydros_loc_filename, ".dat");
        else if (globals->hydros_loc_flag == 2)	RemoveSuffix(globals->hydros_loc_filename, ".csv");
        else if (globals->hydros_loc_flag == 4)	RemoveSuffix(globals->hydros_loc_filename, ".rad");
        else if (globals->hydros_loc_flag == 5)	RemoveSuffix(globals->hydros_loc_filename, ".h5");
        else if (globals->hydros_loc_flag == 6)	RemoveSuffix(globals->hydros_loc_filename, ".h5");
    }
    else if (globals->hydros_loc_flag == 3)
    {
        globals->hydros_loc_filename = (char*)malloc(ASYNCH_MAX_PATH_LENGTH * sizeof(char));
        globals->hydro_table = (char*)malloc(ASYNCH_MAX_SYMBOL_LENGTH * sizeof(char));
        valsread = sscanf(line_buffer, "%*u %lf %s %s", &(globals->print_time), globals->hydros_loc_filename, globals->hydro_table);
        if (ReadLineError(valsread, 3, "hydrographs location"))	return NULL;
        if (!CheckFilenameExtension(globals->hydros_loc_filename, ".dbc"))	return NULL;
        ReadDBC(globals->hydros_loc_filename, &db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);
    }

    //Grab where to write the peakflow data
    ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
    valsread = sscanf(line_buffer, "%hu", &(globals->peaks_loc_flag));
    if (ReadLineError(valsread, 1, "peakflow location"))	return NULL;

    if (globals->peaks_loc_flag == 0)
    {
        globals->peaks_loc_filename = NULL;
    }
    else if (globals->peaks_loc_flag == 1)
    {
        globals->peaks_loc_filename = (char*)malloc(ASYNCH_MAX_PATH_LENGTH * sizeof(char));
        valsread = sscanf(line_buffer, "%*u %s", globals->peaks_loc_filename);
        if (ReadLineError(valsread, 1, "peakflow location"))	return NULL;
        if (!CheckFilenameExtension(globals->peaks_loc_filename, ".pea"))	return NULL;

        RemoveSuffix(globals->peaks_loc_filename, ".pea");
    }
    else if (globals->peaks_loc_flag == 2)
    {
        globals->peaks_loc_filename = (char*)malloc(ASYNCH_MAX_PATH_LENGTH * sizeof(char));
        globals->peak_table = (char*)malloc(ASYNCH_MAX_SYMBOL_LENGTH * sizeof(char));
        valsread = sscanf(line_buffer, "%*u %s %s", globals->peaks_loc_filename, globals->peak_table);
        if (ReadLineError(valsread, 2, "peakflow location"))	return NULL;
        if (!CheckFilenameExtension(globals->peaks_loc_filename, ".dbc"))	return NULL;
        ReadDBC(globals->peaks_loc_filename, &db_connections[ASYNCH_DB_LOC_PEAK_OUTPUT]);
    }

    //Grab the .sav files
    ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
    valsread = sscanf(line_buffer, "%hu", &(globals->hydrosave_flag));
    if (ReadLineError(valsread, 1, "hydrographs save flag"))	return NULL;

    if (globals->hydrosave_flag == 1)
    {
        globals->hydrosave_filename = (char*)malloc(ASYNCH_MAX_PATH_LENGTH * sizeof(char));
        valsread = sscanf(line_buffer, "%*u %s", globals->hydrosave_filename);
        if (ReadLineError(valsread, 1, "hydrographs .sav filename"))	return NULL;
        if (!CheckFilenameExtension(globals->hydrosave_filename, ".sav"))	return NULL;
    }
    else if (globals->hydrosave_flag == 2)
    {
        valsread = sscanf(line_buffer, "%*u %s", db_filename);
        if (ReadLineError(valsread, 1, ".dbc for hydrograph save ids"))	return NULL;
        if (!CheckFilenameExtension(db_filename, ".dbc"))	return NULL;
        globals->hydrosave_filename = NULL;
        ReadDBC(db_filename, &db_connections[ASYNCH_DB_LOC_HYDROSAVE]);
    }
    else
        globals->hydrosave_filename = NULL;

    ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
    valsread = sscanf(line_buffer, "%hu", &(globals->peaksave_flag));
    if (ReadLineError(valsread, 1, "peakflows save flag"))	return NULL;

    if (globals->peaksave_flag == 1)
    {
        globals->peaksave_filename = (char*)malloc(ASYNCH_MAX_PATH_LENGTH * sizeof(char));
        valsread = sscanf(line_buffer, "%*u %s", globals->peaksave_filename);
        if (ReadLineError(valsread, 1, "peakflows .sav filename"))	return NULL;
        if (!CheckFilenameExtension(globals->peaksave_filename, ".sav"))	return NULL;
    }
    else if (globals->peaksave_flag == 2)
    {
        valsread = sscanf(line_buffer, "%*u %s", db_filename);
        if (ReadLineError(valsread, 1, ".dbc for peakflow save ids"))	return NULL;
        if (!CheckFilenameExtension(db_filename, ".dbc"))	return NULL;
        globals->peaksave_filename = NULL;
        ReadDBC(db_filename, &db_connections[ASYNCH_DB_LOC_PEAKSAVE]);
    }
    else
        globals->peaksave_filename = NULL;
    globals->peakfilename = NULL;

    //Grab data dump info
    ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
    valsread = sscanf(line_buffer, "%hu", &(globals->dump_loc_flag));
    if (ReadLineError(valsread, 1, "snapshot save flag"))	return NULL;

    globals->dump_loc_filename = NULL;
    globals->dump_table = NULL;

    if (globals->dump_loc_flag == 1 || globals->dump_loc_flag == 3)
    {
        globals->dump_loc_filename = (char*)malloc(ASYNCH_MAX_PATH_LENGTH * sizeof(char));
        valsread = sscanf(line_buffer, "%*u %s", globals->dump_loc_filename);
        if (ReadLineError(valsread, 1, "snapshot filename"))	return NULL;
        if ((globals->dump_loc_flag == 1) && !CheckFilenameExtension(globals->dump_loc_filename, ".rec"))
            return NULL;
        if ((globals->dump_loc_flag == 3) && !CheckFilenameExtension(globals->dump_loc_filename, ".h5"))
            return NULL;
    }
    else if (globals->dump_loc_flag == 2)
    {
        globals->dump_table = (char*)malloc(ASYNCH_MAX_SYMBOL_LENGTH * sizeof(char));
        valsread = sscanf(line_buffer, "%*u %s %s", db_filename, globals->dump_table);
        if (ReadLineError(valsread, 2, ".dbc for snapshots"))	return NULL;
        if (!CheckFilenameExtension(db_filename, ".dbc"))	return NULL;
        //globals->dump_loc_filename = NULL;
        ReadDBC(db_filename, &db_connections[ASYNCH_DB_LOC_SNAPSHOT_OUTPUT]);
    }
    else if (globals->dump_loc_flag == 4)
    {
        globals->dump_loc_filename = (char*)malloc(ASYNCH_MAX_PATH_LENGTH * sizeof(char));
        valsread = sscanf(line_buffer, "%*u %lf %s", &globals->dump_time, globals->dump_loc_filename);
        if (ReadLineError(valsread, 2, "snapshot time and filename"))	return NULL;
        if (!CheckFilenameExtension(globals->dump_loc_filename, ".h5"))
            return NULL;

        //Strip the extension
        char *ext = strrchr(globals->dump_loc_filename, '.');
        unsigned int filename_len = (unsigned int)(ext - globals->dump_loc_filename);
        globals->dump_loc_filename[filename_len] = '\0';
    }

    //Grab folder locations
    //globals->results_folder = (char*) malloc(string_size*sizeof(char));
    globals->temp_filename = (char*)malloc(ASYNCH_MAX_PATH_LENGTH * sizeof(char));
    //ReadLineFromTextFile(globalfile,line_buffer,line_buffer_len,string_size);
    //valsread = sscanf(line_buffer,"%s",globals->results_folder);
    //if(ReadLineError(valsread,1,"results folder"))	return NULL;
    ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
    valsread = sscanf(line_buffer, "%s", globals->temp_filename);
    if (ReadLineError(valsread, 1, "scratch work folder"))	return NULL;

    if (globals->print_par_flag)	//!!!! Is this needed? Why bother? !!!!
    {
        if (AttachParameters(globals->temp_filename, ASYNCH_MAX_PATH_LENGTH, globals->global_params, globals->num_global_params, ASYNCH_MAX_PATH_LENGTH))
        {
            printf("[%i]: Error attaching global parameters to temporary filenames.\n", my_rank);
            return NULL;
        }
    }

    sprintf(db_filename, "_%i_%i", getpid(), my_rank);
    strcat(globals->temp_filename, db_filename);

    //Grab adapative data
    ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
    valsread = sscanf(line_buffer, "%lf %lf %lf", &errors->facmin, &errors->facmax, &errors->fac);
    if (ReadLineError(valsread, 3, "facmin, facmax, fac"))	return NULL;

    //Read in the flag for the error tolerances
    ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
    valsread = sscanf(line_buffer, "%hi", &(globals->rkd_flag));
    if (ReadLineError(valsread, 1, "error tolerance flag"))	return NULL;

    //Set some parameters
    globals->max_rk_stages = 0;
    globals->max_parents = 0;

    if (globals->rkd_flag == 0)	//Error data is found in the universal file
    {
        rkdfilename[0] = '\0';
        ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
        valsread = sscanf(line_buffer, "%u", &flag);
        if (ReadLineError(valsread, 1, "RK method index"))	return NULL;
        globals->method = flag;

        //Count the number of states given
        ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
        total = 0;
        int error_dim = -1;
        double tempy;
        do
        {
            error_dim++;
            valsread = sscanf(&(line_buffer[total]), "%lf%n", &tempy, &written);	//Need to actually read to trigger valsread
            total += written;
        } while (valsread > 0);

        //Reserve memory
        errors->abstol = calloc(error_dim, sizeof(double));
        errors->reltol = calloc(error_dim, sizeof(double));
        errors->abstol_dense = calloc(error_dim, sizeof(double));
        errors->reltol_dense = calloc(error_dim, sizeof(double));

        //ReadLineFromTextFile(globalfile,line_buffer,line_buffer_len,string_size);
        total = 0;
        for (int i = 0; i < error_dim; i++)	//Note: Don't read the line from disk again
        {
            valsread = sscanf(&(line_buffer[total]), "%lf%n", &errors->abstol[i], &written);
            //if(ReadLineError(valsread,1,"an abstol component"))	return NULL;
            total += written;
        }

        ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
        total = 0;
        for (int i = 0; i < error_dim; i++)
        {
            valsread = sscanf(&(line_buffer[total]), "%lf%n", &errors->reltol[i], &written);
            if (ReadLineError(valsread, 1, "a reltol component"))	return NULL;
            total += written;
        }

        ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
        total = 0;
        for (int i = 0; i < error_dim; i++)
        {
            valsread = sscanf((&line_buffer[total]), "%lf%n", &errors->abstol_dense[i], &written);
            if (ReadLineError(valsread, 1, "an abstol dense output component"))	return NULL;
            total += written;
        }

        ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
        total = 0;
        for (int i = 0; i < error_dim; i++)
        {
            valsread = sscanf(&(line_buffer[total]), "%lf%n", &errors->reltol_dense[i], &written);
            if (ReadLineError(valsread, 1, "a reltol dense output component"))	return NULL;
            total += written;
        }
    }
    else if (globals->rkd_flag == 1)	//Error data is in an .rkd file
    {
        globals->method = -1;
        errors->abstol = NULL;
        errors->reltol = NULL;
        errors->abstol_dense = NULL;
        errors->reltol_dense = NULL;

        //ReadLineFromTextFile(globalfile,line_buffer,line_buffer_len);
        valsread = sscanf(line_buffer, "%*i %s", rkdfilename);
        if (ReadLineError(valsread, 1, ".rkd filename"))
            return NULL;
        if (!CheckFilenameExtension(rkdfilename, ".rkd"))
            return NULL;
    }
    else
    {
        printf("Error: bad flag for error tolerance. Got %hu.\n", globals->rkd_flag);
        return NULL;
    }

    //Check for end mark
    ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
    sscanf(line_buffer, "%c", &endmark);
    if (endmark != '#')
    {
        printf("Error: an ending # not seen in %s file. Got %c.\n", globalfilename, endmark);
        return NULL;
    }

    //Setup an io object
    OutputFunc_Init(
        globals->hydros_loc_flag,
        globals->peaks_loc_flag,
        globals->dump_loc_flag,
        &globals->output_func);

    //Clean up
    if (my_rank == 0)
        fclose(globalfile);

    return globals;
}


//Put a vector of global_params onto the end of a filename.
//Returns 1 if filename is not long enough to support this.
//Retruns 0 if the parameters are attached.
int AttachParameters(char* filename, unsigned int max_size, double *global_params, unsigned int num_global_params, unsigned int string_size)
{
    unsigned int i, count, total = 0;
    char *buffer;
    size_t length = strlen(filename);

    //TODO use snprintf
    buffer = (char *)malloc(string_size);

    for (i = 0; i < num_global_params; i++)
    {
        sprintf(buffer, "_%.4e%n", global_params[i], &count);
        total += count;
        if (count + 1 > string_size)	return 1;
        if (total + 1 > max_size)		return 1;
        strcat(filename, buffer);
    }

    free(buffer);

    return 0;
}



void ReadLineFromTextFile(FILE *globalfile, char *line_buffer, unsigned int size)
{
    unsigned int line_buffer_length;
    if (my_rank == 0)
    {
        line_buffer[0] = '%';
        while (!feof(globalfile) && (line_buffer[0] == '%' || line_buffer[0] == '\n'))
            fgets(line_buffer, size, globalfile);

        line_buffer_length = (unsigned int)strlen(line_buffer);
        if (line_buffer_length > size - 1)
            printf("Warning: %zu %u Line in .gbl file may be too long. Read in the long line:\n\"%s\"\n", strlen(line_buffer), size, line_buffer);
    }

    MPI_Bcast(&line_buffer_length, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(line_buffer, line_buffer_length + 1, MPI_CHAR, 0, MPI_COMM_WORLD);
}


int ReadLineError(int valsread, int valswant, char *message)
{
    if (valsread < valswant)
    {
        if (my_rank == 0)	printf("Error: Did not get a value from .gbl file for %s. %i\n", message, valsread);
        return 1;
    }
    return 0;
}


//Reads a .dbc file and creates a corresponding database connection.
//This does NOT connect to the database.
//Returns NULL if there was an error.
void ReadDBC(char *filename, ConnData* const conninfo)
{
    bool res = true;
    unsigned int i = 0, j = 0;
    char connstring[ASYNCH_MAX_CONNSTRING_LENGTH];
    char c;

    //if(my_rank == 0)
    //{
    FILE* input = fopen(filename, "r");

    if (!input)
    {
        printf("Error opening .dbc file %s.\n", filename);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    if (CheckWinFormat(input))
    {
        printf("Error: File %s appears to be in Windows format. Try converting to unix format using 'dos2unix' at the command line.\n", filename);
        fclose(input);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    //Read connection information
    //Currently, this expects 4 things like:
    fgets(connstring, ASYNCH_MAX_CONNSTRING_LENGTH, input);
    ConnData_Init(conninfo, connstring);
    if (!conninfo)
    {
        fclose(input);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    //Get number of queries
    if (fscanf(input, "%u", &(conninfo->num_queries)) == EOF)
    {
        printf("[%i]: Error: failed to parse file.\n", my_rank);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    for (i = 0; i < conninfo->num_queries; i++)
        conninfo->queries[i] = (char*)malloc(ASYNCH_MAX_QUERY_LENGTH * sizeof(char));

    //Get queries. They are delineated by a ;
    for (j = 0; j < conninfo->num_queries; j++)
    {
        //Get rid of whitespace
        c = fgetc(input);
        while (c != EOF && (c == ' ' || c == '\n' || c == '\t'))	c = fgetc(input);
        if (c == EOF)
        {
            printf("[%i]: Warning: did not see %u queries in %s.\n", my_rank, conninfo->num_queries, filename);
            break;
        }

        //Read in query
        for (i = 0; i < ASYNCH_MAX_QUERY_LENGTH - 2 && c != ';' && c != EOF; i++)
        {
            conninfo->queries[j][i] = c;
            c = fgetc(input);
        }

        //Check for problems and put stuff on the end
        if (i == ASYNCH_MAX_QUERY_LENGTH)
            printf("[%i]: Warning: query %u is too long in %s.\n", my_rank, j, filename);
        else if (c == EOF)
        {
            printf("[%i]: Warning: did not see %u queries in %s.\n", my_rank, conninfo->num_queries, filename);
            break;
        }
        else
        {
            conninfo->queries[j][i] = ';';
            conninfo->queries[j][i + 1] = '\0';
        }
    }

    fclose(input);

    ////Get string length for other procs
    //j = (unsigned int) strlen(connstring) + 1;
    //}

    ////Check if an error occurred
    //finish:
    //MPI_Bcast(&has_error,1,MPI_C_BOOL,0,MPI_COMM_WORLD);

    ////Transfer info from .dbc file
    //if(!errorcode)
    //{
    //	MPI_Bcast(&j,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
    //	MPI_Bcast(connstring,j,MPI_CHAR,0,MPI_COMM_WORLD);
    //	if(my_rank != 0)
    //		ConnData_Init(conninfo, connstring);
    //	MPI_Bcast(&(conninfo->num_queries),1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
    //	if(my_rank == 0)
    //	{
    //		for(i=0;i<conninfo->num_queries;i++)
    //		{
    //			j = (unsigned int) strlen(conninfo->queries[i]) + 1;
    //			MPI_Bcast(&j,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
    //			MPI_Bcast(conninfo->queries[i],j,MPI_CHAR,0,MPI_COMM_WORLD);
    //		}
    //	}
    //	else
    //	{
    //		conninfo->queries = (char**) malloc(conninfo->num_queries*sizeof(char*));
    //		for(i=0;i<conninfo->num_queries;i++)
    //		{
    //			conninfo->queries[i] = (char*) malloc(ASYNCH_MAX_QUERY_LENGTH*sizeof(char));
    //			MPI_Bcast(&j,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
    //			MPI_Bcast(conninfo->queries[i],j,MPI_CHAR,0,MPI_COMM_WORLD);
    //		}
    //	}
    //}


}

//Checks that the filename ends with the approriate extension.
//Returns 1 if the extension is correct, 0 if not.
int CheckFilenameExtension(char *filename, char *extension)
{
    size_t l_filename, l_extension, i;

    if (!extension)	return 1;
    if (!filename)
    {
        if (my_rank == 0) printf("Error: no filename given.\n");
        return 0;
    }

    l_filename = strlen(filename);
    l_extension = strlen(extension);

    if (l_filename < l_extension)
    {
        if (my_rank == 0) printf("Error: filename %s is of wrong file type (expected %s).\n", filename, extension);
        return 0;
    }
    for (i = 0; i < l_extension; i++)
    {
        if (extension[i] != filename[l_filename - l_extension + i])
        {
            if (my_rank == 0) printf("Error: filename %s is of wrong file type (expected %s).\n", filename, extension);
            return 0;
        }
    }

    return 1;
}


//Removes a suffix from filename, if present.
//Returns 1 if suffix removed
//0 if not (not present)
int RemoveSuffix(char *filename, const char *suffix)
{
    size_t filename_length = strlen(filename);
    size_t suffix_length = strlen(suffix);

    if (suffix_length > filename_length)	return 0;

    char *dot = strrchr(filename, '.');
    if (!dot || dot == filename || strcmp(dot, suffix) != 0)
        return 0;

    *dot = '\0';
    return 1;
}
