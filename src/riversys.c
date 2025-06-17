#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <assert.h>
#include <math.h>
#include <memory.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#if defined(HAVE_UNISTD_H)
#include <unistd.h>
#endif

#if defined(_MSC_VER)
#include <process.h>
#endif

#if defined(HAVE_POSTGRESQL)
#include <libpq-fe.h>
#endif

#if defined(HAVE_HDF5)
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

#include <data_types.h>
#include <rkmethods.h>
#include <rksteppers.h>
#include <partition.h>
#include <db.h>
#include <comm.h>
#include <io.h>
#include <forcings.h>
#include <forcings_io.h>
#include <outputs.h>
//#include <builtin.h>
//#include <vector_mpi.h>
#include <minmax.h>
#include <blas.h>
#include <models/definitions.h>

#include <riversys.h>


//Read topo data and build the network.
//Also creates id_to_loc.
void Create_River_Network(GlobalVars* globals, Link** system, unsigned int* N, Lookup** id_to_loc, ConnData* db_connections)
{
    FILE* riverdata = NULL;
    unsigned int *link_ids = NULL;
    unsigned int *dbres_link_id = NULL;
    unsigned int *dbres_parent = NULL;
    unsigned int sizeres = 0;
    unsigned int i, j;
    unsigned int max_children = 10;
    unsigned int *num_parents = NULL;
    unsigned int **loc_to_children = NULL;
    unsigned int *loc_to_children_array = NULL;
    unsigned int curr_loc;

    *system = NULL;
    *N = 0;

    if (globals->rvr_flag == 0)	//Read topo data from file
    {
        if (my_rank == 0)
        {
            riverdata = fopen(globals->rvr_filename, "r");
            if (!riverdata)
            {
                if (my_rank == 0)	printf("Error: file %s not found for .rvr file.\n", globals->rvr_filename);
                *N = 0;
                return;
            }
            if (CheckWinFormat(riverdata))
            {
                printf("Error: File %s appears to be in Windows format. Try converting to unix format using 'dos2unix' at the command line.\n", globals->rvr_filename);
                fclose(riverdata);
                return;
            }

            fscanf(riverdata, "%u", N);
            link_ids = (unsigned int*)malloc(*N * sizeof(unsigned int));
            loc_to_children_array = (unsigned int*)calloc(*N*max_children, sizeof(unsigned int));
            loc_to_children = (unsigned int**)malloc(*N * sizeof(unsigned int*));	//This holds the ID of the children
            for (i = 0; i < *N; i++)	
                loc_to_children[i] = &(loc_to_children_array[i*max_children]);
            num_parents = (unsigned int*)malloc(*N * sizeof(unsigned int));

            for (i = 0; i < *N; i++)
            {
                fscanf(riverdata, "%u %u", &(link_ids[i]), &(num_parents[i]));
                for (j = 0; j < num_parents[i]; j++)
                {
                    unsigned int id;
                    fscanf(riverdata, "%u", &id);
                    loc_to_children[i][j] = id;
                }

                if (num_parents[i] > max_children)
                {
                    printf("Error: assumed no link has more than %u parents, but link %u has %u.\n", max_children, link_ids[i], num_parents[i]);
                    printf("If this is not an error in the input data, modify max_children in riversys.c, function Create_River_Network.\n");
                    return;
                }
            }

            fclose(riverdata);

            //Broadcast data
            MPI_Bcast(N, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
            MPI_Bcast(link_ids, *N, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
            MPI_Bcast(loc_to_children_array, *N*max_children, MPI_UNSIGNED, 0, MPI_COMM_WORLD);	//Yeah, this is probably not the most efficient...
            MPI_Bcast(num_parents, *N, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        }
        else
        {
            MPI_Bcast(N, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
            link_ids = (unsigned int*)malloc(*N * sizeof(unsigned int));
            loc_to_children_array = (unsigned int*)malloc(*N*max_children * sizeof(unsigned int));
            loc_to_children = (unsigned int**)malloc(*N * sizeof(unsigned int*));	//This holds the ID of the children
            for (i = 0; i < *N; i++)	
                loc_to_children[i] = &(loc_to_children_array[i*max_children]);
            num_parents = (unsigned int*)malloc(*N * sizeof(unsigned int));
            MPI_Bcast(link_ids, *N, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
            MPI_Bcast(loc_to_children_array, *N*max_children, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
            MPI_Bcast(num_parents, *N, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        }
    }
    else if (globals->rvr_flag == 1)	//Download topo data from database
    {

#if defined(HAVE_POSTGRESQL)
        PGresult *mainres, *res;
        int db;

        if (my_rank == 0)
        {
            if (globals->outletlink == 0)	//Grab entire network
            {
                //Note: parent_id on the SQL database is like the child id here
                db = ConnectPGDB(&db_connections[ASYNCH_DB_LOC_TOPO]);
                if (db)
                {
                    printf("[%i]: Error connecting to the topology database.\n", my_rank);
                    return;
                }
                mainres = PQexec(db_connections[ASYNCH_DB_LOC_TOPO].conn, db_connections[ASYNCH_DB_LOC_TOPO].queries[0]);	//For all link ids
                res = PQexec(db_connections[ASYNCH_DB_LOC_TOPO].conn, db_connections[ASYNCH_DB_LOC_TOPO].queries[1]);		//For parent data
                if (CheckResError(res, "querying connectivity") || CheckResError(mainres, "querying DEM data"))
                    return;
                DisconnectPGDB(&db_connections[ASYNCH_DB_LOC_TOPO]);

                *N = PQntuples(mainres);

                //Get the list of link ids
                link_ids = (unsigned int*)malloc(*N * sizeof(unsigned int));
                for (i = 0; i < *N; i++)
                    link_ids[i] = atoi(PQgetvalue(mainres, i, 0));
                PQclear(mainres);

                sizeres = PQntuples(res);
            }
            else	//Grab a sub basin (Note: all link ids (except the outlet) is in the second column of res)
            {
                db = ConnectPGDB(&db_connections[ASYNCH_DB_LOC_TOPO]);
                if (db)
                {
                    printf("[%i]: Error connecting to the topology database.\n", my_rank);
                    return;
                }

                //Make the queries
                //Be careful to not overload the database
                sprintf(db_connections[ASYNCH_DB_LOC_TOPO].query, db_connections[ASYNCH_DB_LOC_TOPO].queries[2], globals->outletlink);	//For parent data
                res = PQexec(db_connections[ASYNCH_DB_LOC_TOPO].conn, db_connections[ASYNCH_DB_LOC_TOPO].query);
                if (CheckResError(res, "querying connectivity"))
                    return;
                DisconnectPGDB(&db_connections[ASYNCH_DB_LOC_TOPO]);

                *N = PQntuples(res) + 1;

                //Get the list of link ids
                link_ids = (unsigned int*)malloc(*N * sizeof(unsigned int));
                for (i = 0; i < *N - 1; i++)	
                    link_ids[i] = atoi(PQgetvalue(res, i, 1));
                link_ids[i] = globals->outletlink;
                merge_sort_1D(link_ids, *N);

                sizeres = *N - 1;
            }

            dbres_link_id = (unsigned int*)malloc(sizeres * sizeof(unsigned int));
            dbres_parent = (unsigned int*)malloc(sizeres * sizeof(unsigned int));

            for (i = 0; i < sizeres; i++)
            {
                dbres_link_id[i] = atoi(PQgetvalue(res, i, 0));
                dbres_parent[i] = atoi(PQgetvalue(res, i, 1));
            }

            PQclear(res);

            //Modify the data format
            loc_to_children_array = (unsigned int*)calloc(*N*max_children, sizeof(unsigned int));
            loc_to_children = (unsigned int**)malloc(*N * sizeof(unsigned int*));	//This holds the IDs of the parents. Really needs a better name...
            for (i = 0; i < *N; i++)	
                loc_to_children[i] = &(loc_to_children_array[i*max_children]);
            num_parents = (unsigned int*)malloc(*N * sizeof(unsigned int));

            curr_loc = 0;
            for (i = 0; i < sizeres; i += j)
            {
                //Select the next link
                unsigned int id = dbres_link_id[i];

                if (id == link_ids[curr_loc])	//Extract parents
                {
                    //Count the parents
                    for (j = 0; i + j < sizeres; j++)
                        if (dbres_link_id[i + j] != id)	break;

                    //Set the number of parents
                    num_parents[curr_loc] = j;

                    //Set the information to find the child links
                    for (unsigned int k = 0; k < j; k++)
                        loc_to_children[curr_loc][k] = dbres_parent[i + k];
                }
                else	//No parents
                {
                    num_parents[curr_loc] = 0;
                    j = 0;
                }

                //Set the next location
                curr_loc++;
            }
            for (; curr_loc < *N; curr_loc++)	//Finish any leaves at the end of link_ids
                num_parents[curr_loc] = 0;

            //Send sizes and data to other processes
            MPI_Bcast(N, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(link_ids, *N, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
            MPI_Bcast(loc_to_children_array, *N*max_children, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
            MPI_Bcast(num_parents, *N, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

            //Clean up
            free(dbres_link_id);
            free(dbres_parent);
        }
        else
        {
            //Receive data
            MPI_Bcast(N, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
            link_ids = (unsigned int*)malloc(*N * sizeof(unsigned int));
            loc_to_children_array = (unsigned int*)calloc(*N*max_children, sizeof(unsigned int));
            loc_to_children = (unsigned int**)malloc(*N * sizeof(unsigned int*));	//This holds the ID of the children
            for (i = 0; i < *N; i++)	loc_to_children[i] = &(loc_to_children_array[i*max_children]);
            num_parents = (unsigned int*)calloc(*N, sizeof(unsigned int));
            MPI_Bcast(link_ids, *N, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
            MPI_Bcast(loc_to_children_array, *N*max_children, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
            MPI_Bcast(num_parents, *N, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        }

#else //HAVE_POSTGRESQL

        if (my_rank == 0)	printf("Error: Asynch was build without PostgreSQL support.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

#endif //HAVE_POSTGRESQL
    }
    else
    {
        if (my_rank == 0)	printf("Error: Bad topology flag %hi in .gbl file.\n", globals->rvr_flag);
        *N = 0;
        return;
    }

    //Make a list of ids and locations, sorted by id
    *id_to_loc = malloc(*N * sizeof(Lookup));
    for (i = 0; i < *N; i++)
    {
        (*id_to_loc)[i].id = link_ids[i];
        (*id_to_loc)[i].loc = i;
    }
    merge_sort_by_ids(*id_to_loc, *N);

    //Check for crapiness
    if (my_rank == 0 && *N < (unsigned int)np)
        printf("\nWarning: using more processes (%i) than links (%u).\n", np, *N);

    //Allocate some space for the network
    Link* sys = (Link*)calloc(*N, sizeof(Link));
    *system = sys;

    //Build the network
    globals->max_parents = 0;
    for (i = 0; i < *N; i++)
    {
        sys[i].my = NULL;

        sys[i].location = i;
        sys[i].ID = link_ids[i];

        //Set the parents
        sys[i].num_parents = num_parents[i];
        sys[i].parents = (Link**)calloc(sys[i].num_parents, sizeof(Link*));
        sys[i].child = NULL;
        globals->max_parents = max(globals->max_parents, num_parents[i]);

        //Set a few other data
        sys[i].params = NULL;
        sys[i].discont_end = globals->discont_size - 1;

        sys[i].output_user = NULL;
        //system[i].peakoutput_user = NULL;
        //system[i].differential = NULL;
        //system[i].params = v_init(0);
        //system[i].is_dam = 0;
        //system[i].method = NULL;
        //system[i].error_data = NULL;
        //system[i].qvs = NULL;
        //system[i].discont = NULL;
        //system[i].discont_start = 0;
        //system[i].discont_end = globals->discont_size-1;
        //system[i].discont_count = 0;
        //system[i].discont_send = NULL;
        //system[i].discont_order_send = NULL;
        //system[i].discont_send_count = 0;
        //system[i].res = 0;
        //system[i].num_dense = 0;
        //system[i].dense_indices = NULL;
        //system[i].dim = 0;
        //system[i].diff_start = 0;
        //system[i].no_ini_start = 0;
        //system[i].disk_iterations = 0;
        //system[i].forcing_data = NULL;
        //system[i].forcing_values = NULL;
        //system[i].forcing_change_times = NULL;
        //system[i].forcing_indices = NULL;
        //system[i].save_flag = 0;
        //system[i].peak_flag = 0;
        //system[i].user = NULL;
        //system[i].equations = NULL;
    }

    //Setup the child and parent information
    for (i = 0; i < *N; i++)
    {
        for (j = 0; j < sys[i].num_parents; j++)
        {
            curr_loc = find_link_by_idtoloc(loc_to_children[i][j], *id_to_loc, *N);
            if (curr_loc > *N)
            {
                if (my_rank == 0)
                    printf("Error: Invalid id in topology data (%u).\n", loc_to_children[i][j]);
                    printf("N: (%u) curr_loc (%u).\n", *N, curr_loc);
                *N = 0;
                return;
            }
            sys[i].parents[j] = &sys[curr_loc];
            sys[curr_loc].child = &sys[i];
        }
    }

    //Set an outletlink id. This only sets one outlet, and only if rvr_flag is not set.
    if (globals->prm_flag == 1 && globals->rvr_flag == 0)
    {
        for (i = 0; i < *N; i++)
            if (sys[i].child == NULL)	break;
        globals->outletlink = sys[i].ID;
    }

    //Clean up
    if (loc_to_children)
        free(loc_to_children);
    if (loc_to_children_array)
        free(loc_to_children_array);
    if (link_ids)
        free(link_ids);
    if (num_parents)
        free(num_parents);
}



//Read in the local paramters for the network.
//Returns 1 if there is an error, 0 otherwise.
//If load_all == 1, then the parameters for every link are available on every proc.
//If load_all == 0, then the parameters are only available for links assigned to this proc.
int Load_Local_Parameters(
    Link *system, unsigned int N,
    Link ** my_sys, unsigned int my_N,
    int* assignments, short int* getting, const Lookup * const id_to_loc,
    const GlobalVars * const globals,
    ConnData* db_connections,
    AsynchModel* model,
    void* external)
{
    unsigned int *db_link_id, curr_loc;
    double *db_params_array, **db_params;
    FILE* paramdata;

    //Error checking
    if (!assignments || !getting)
    {
        if (my_rank == 0)
        {
            printf("Error loading link parameters: network partitioning must occur before loading parameters.\n");
            printf("(Hint: either partition the network before loading parameters, or load parameters at every link.)\n");
        }
        return 1;
    }

    //Allocate space
    db_link_id = (unsigned int*)malloc(N * sizeof(unsigned int));
    db_params_array = (double*)malloc(N * globals->num_disk_params * sizeof(double));
    db_params = (double**)malloc(N * sizeof(double*));
    for (unsigned int i = 0; i < N; i++)
        db_params[i] = &(db_params_array[i * globals->num_disk_params]);

    //Read parameters
    if (my_rank == 0)
    {
        if (globals->prm_flag == 0)
        {
            paramdata = fopen(globals->prm_filename, "r");
            if (paramdata == NULL)
            {
                printf("Error: file %s not found for .prm file.\n", globals->prm_filename);
                return 1;
            }
            if (CheckWinFormat(paramdata))
            {
                printf("Error: File %s appears to be in Windows format. Try converting to unix format using 'dos2unix' at the command line.\n", globals->prm_filename);
                fclose(paramdata);
                return 1;
            }

            unsigned int n;
            fscanf(paramdata, "%u", &n);
            if (n != N)
            {
                printf("Error: expected %u links in parameter file. Got %u.\n", N, n);
                return 1;
            }

            for (unsigned int i = 0; i < N; i++)
            {
                fscanf(paramdata, "%u", &(db_link_id[i]));
                for (unsigned int j = 0; j < globals->num_disk_params; j++)
                {
                    if (fscanf(paramdata, "%lf", &(db_params[i][j])) == 0)
                    {
                        printf("Error reading from parameter file %s.\n", globals->prm_filename);
                        return 1;
                    }
                }
            }

            fclose(paramdata);
        }
        else if (globals->prm_flag == 1)
        {
#if defined(HAVE_POSTGRESQL)
            int db;
            PGresult *res;

            if (globals->outletlink == 0)	//Grab entire network
            {
                db = ConnectPGDB(&db_connections[ASYNCH_DB_LOC_PARAMS]);
                if (db)
                {
                    printf("[%i]: Error connecting to the parameter database.\n", my_rank);
                    return 1;
                }
                res = PQexec(db_connections[ASYNCH_DB_LOC_PARAMS].conn, db_connections[ASYNCH_DB_LOC_PARAMS].queries[0]);
                if (CheckResError(res, "querying DEM data"))	return 1;
                DisconnectPGDB(&db_connections[ASYNCH_DB_LOC_PARAMS]);
            }
            else	//Grab a sub basin
            {
                db = ConnectPGDB(&db_connections[ASYNCH_DB_LOC_PARAMS]);
                if (db)
                {
                    printf("[%i]: Error connecting to the parameter database.\n", my_rank);
                    return 1;
                }

                //Make the queries
                sprintf(db_connections[ASYNCH_DB_LOC_PARAMS].query, db_connections[ASYNCH_DB_LOC_PARAMS].queries[1], globals->outletlink);
                res = PQexec(db_connections[ASYNCH_DB_LOC_PARAMS].conn, db_connections[ASYNCH_DB_LOC_PARAMS].query);
                if (CheckResError(res, "querying DEM data"))	return 1;
                DisconnectPGDB(&db_connections[ASYNCH_DB_LOC_PARAMS]);
            }

            unsigned int n = PQntuples(res);
            if (n != N)
            {
                printf("Error processing link parameters: Got %u, expected %u.\n(Hint: make sure your topology and parameter sources have the same number of links.)\n", n, N);
                return 1;
            }

            //Load buffers
            for (unsigned int i = 0; i < N; i++)
                db_link_id[i] = atoi(PQgetvalue(res, i, 0));
            for (unsigned int i = 0; i < N; i++)
                for (unsigned int j = 0; j < globals->num_disk_params; j++)
                    db_params[i][j] = atof(PQgetvalue(res, i, 1 + j));

            //Cleanup
            PQclear(res);

#else //HAVE_POSTGRESQL

            if (my_rank == 0)	printf("Error: Asynch was build without PostgreSQL support.\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

#endif //HAVE_POSTGRESQL
        }
    }

    //Broadcast data
    MPI_Bcast(db_link_id, N, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(db_params_array, N * globals->num_disk_params, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (unsigned int i = 0; i < N; i++)
    {
        curr_loc = find_link_by_idtoloc(db_link_id[i], id_to_loc, N);
        if (curr_loc > N)
        {
            if (my_rank == 0)	printf("Error: link id %u appears in the link parameters, but not in the topology data.\n(Hint: make sure your topology and parameter sources are correct.)\n", db_link_id[i]);
            return 1;
        }
        else
        {
            if (assignments[curr_loc] == my_rank || getting[curr_loc])
            {
                system[curr_loc].num_params = globals->num_params;
                system[curr_loc].params = malloc(globals->num_params * sizeof(double));
                for (unsigned int j = 0; j < globals->num_disk_params; j++)
                    system[curr_loc].params[j] = db_params[i][j];

                if (model)
                    model->convert(system[curr_loc].params, globals->model_uid, external);
                else
                    ConvertParams(system[curr_loc].params, globals->model_uid, external);
            }
        }
    }

    //Clean up
    free(db_link_id);
    free(db_params_array);
    free(db_params);

    return 0;
}



//Partitions the network amongst different MPI processes.
//!!!! Perhaps the leaves info could be moved deeper? How about errors here? !!!!
int Partition_Network(
    Link *system, unsigned int N,
    const GlobalVars * const globals,
    Link*** my_sys, unsigned int* my_N, int** assignments,
    TransData** my_data, short int** getting, AsynchModel* model)
{
    Link *current, *prev;

    //Perform a DFS to sort the leaves
    Link** stack = malloc(N * sizeof(Link*)); //Holds the index in system
    int stack_size = 0;

    Link** leaves = calloc(N, sizeof(Link*));
    unsigned int leaves_size = 0;
    unsigned short int num_parents;

    for (unsigned int j = 0; j < N; j++)	//!!!! Iterate over a list of roots? No need to sort then... !!!!
    {
        if (system[j].child == NULL)
        {
            stack[0] = &system[j];
            stack_size = 1;
            while (stack_size > 0)
            {
                current = stack[stack_size - 1];	//Top of stack
                num_parents = current->num_parents;

                if (num_parents == 0)
                {
                    stack_size--;
                    leaves[leaves_size] = current;
                    leaves_size++;
                }
                else
                {
                    //If current is not a leaf, replace it with its parents
                    for (unsigned int i = 0; i < num_parents; i++)
                    {
                        stack[stack_size - 1 + i] = current->parents[num_parents - 1 - i];
                        //stack[stack_size - 1 + i]->child = current;
                    }
                    stack_size += num_parents - 1;
                }
            }
        }
        //else
        //	break;
    }

    //!!!! Why is this needed? !!!!

    //for(i=0;i<N;i++)	system[i].distance = 0;

    //Calculate the distance and number of upstream links for each link
    for (unsigned int i = 0; i < leaves_size; i++)
    {
        prev = leaves[i];
        prev->distance = 1;
        for (current = prev->child; current != NULL; current = current->child)
        {
            if (current->distance > prev->distance + 1)
                break;
            else
                current->distance = prev->distance + 1;
            prev = current;
        }
    }

    //Partition the system and assign the links
    *my_data = Initialize_TransData();
    *getting = (short int*)malloc(N * sizeof(short int));
    if (model && model->partition)
        *assignments = model->partition(system, N, leaves, leaves_size, my_sys, my_N, *my_data, *getting);
    else
    {
        *assignments = Partition_System_By_Leaves(system, N, leaves, leaves_size, my_sys, my_N, *my_data, *getting);
        //*assignments = Partition_System_By_Leaves_2(system,N,leaves,leaves_size,my_sys,my_N,*my_data,*getting);
        //*assignments = Partition_METIS_Traditional(system,N,leaves,leaves_size,my_sys,my_N,*my_data,*getting,GlobalVars);
        //*assignments = Partition_METIS_RainChanges(system,N,leaves,leaves_size,my_sys,my_N,*my_data,*getting,GlobalVars);
        //*assignments = Partition_METIS_RainVolume(system,N,leaves,leaves_size,my_sys,my_N,*my_data,*getting,GlobalVars);	//!!!! Requires params for all links !!!!
    }

    //Clean up
    free(stack);
    free(leaves);
    //free(upstream_order);

    //Allocate the link data according to the partition
    for (unsigned int i = 0; i < N; i++)
        if ((*assignments)[i] == my_rank || (*getting)[i])
        {
            system[i].my = malloc(sizeof(LinkData));
            memset(system[i].my, 0, sizeof(LinkData));
        }

    return 0;
}


//Reads numerical error tolerances. Builds RK methods.
//!!!! I'm not really sure how to handle specifying the dimension here. Should the rkd file allow a variable number of tols? !!!!
int Build_RKData(
    Link *system, unsigned int N,
    Link **my_sys, unsigned int my_N,
    char rk_filename[],
    int* assignments, short int* getting,
    GlobalVars *globals,
    ErrorData* error_data,
    RKMethod** methods,
    unsigned int* num_methods)
{
    FILE* rkdata;
    double *filedata_abs, *filedata_rel, *filedata_abs_dense, *filedata_rel_dense;

    //Build all the RKMethods
    static RKMethod rk_methods[4];
    RKDense3_2(&rk_methods[0]);
    TheRKDense4_3(&rk_methods[1]);
    DOPRI5_dense(&rk_methods[2]);
    RadauIIA3_dense(&rk_methods[3]);

    *methods = rk_methods;
    *num_methods = 4;

    globals->max_localorder = rk_methods[0].localorder;
    globals->max_rk_stages = rk_methods[0].num_stages;
    for (unsigned int i = 1; i < *num_methods; i++)
    {
        globals->max_localorder = (globals->max_localorder < rk_methods[i].localorder) ? rk_methods[i].localorder : globals->max_localorder;
        globals->max_rk_stages = (globals->max_rk_stages > rk_methods[i].num_stages) ? globals->max_rk_stages : rk_methods[i].num_stages;
        //!!!! Note: Use a +1 for Radau solver? !!!!
    }

    if (rk_filename[0] != '\0')
    {
        unsigned int *link_ids = (unsigned int*)malloc(N * sizeof(unsigned int));
        unsigned int *rk_methods_idx;
        unsigned int num_states;

        if (my_rank == 0)
        {
            rkdata = fopen(rk_filename, "r");
            if (rkdata == NULL)
            {
                printf("Error: file %s not found for .rkd file\n", rk_filename);
                return 1;
            }
            if (CheckWinFormat(rkdata))
            {
                printf("Error: File %s appears to be in Windows format. Try converting to unix format using 'dos2unix' at the command line.\n", rk_filename);
                fclose(rkdata);
                return 1;
            }

            unsigned int n;
            fscanf(rkdata, "%u %u", &n, &num_states);
            if (n != N)
            {
                printf("Error: the number of links in the rkd file differ from the number in the topology data (Got %u, expected %u).\n", n, N);
                return 1;
            }

            MPI_Bcast(&num_states, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
            filedata_abs = (double*)malloc(N*num_states * sizeof(double));
            filedata_rel = (double*)malloc(N*num_states * sizeof(double));
            filedata_abs_dense = (double*)malloc(N*num_states * sizeof(double));
            filedata_rel_dense = (double*)malloc(N*num_states * sizeof(double));
            rk_methods_idx = (unsigned int*)malloc(num_states * sizeof(unsigned int));

            //Read the file
            for (unsigned int i = 0; i < N; i++)
            {
                if (fscanf(rkdata, "%u", &(link_ids[i])) == 0)
                {
                    printf("Error reading .rkd file: Not enough links in file (expected %u, got %u).\n", N, i);
                    return 1;
                }
                for (unsigned int j = 0; j < num_states; i++)	fscanf(rkdata, "%lf", &(filedata_abs[i*num_states + j]));
                for (unsigned int j = 0; j < num_states; i++)	fscanf(rkdata, "%lf", &(filedata_rel[i*num_states + j]));
                for (unsigned int j = 0; j < num_states; i++)	fscanf(rkdata, "%lf", &(filedata_abs_dense[i*num_states + j]));
                for (unsigned int j = 0; j < num_states; i++)	fscanf(rkdata, "%lf", &(filedata_rel_dense[i*num_states + j]));
                fscanf(rkdata, "%u", &(rk_methods_idx[i]));
            }
        }
        else
        {
            MPI_Bcast(&num_states, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
            filedata_abs = (double*)malloc(N*num_states * sizeof(double));
            filedata_rel = (double*)malloc(N*num_states * sizeof(double));
            filedata_abs_dense = (double*)malloc(N*num_states * sizeof(double));
            filedata_rel_dense = (double*)malloc(N*num_states * sizeof(double));
            rk_methods_idx = (unsigned int*)malloc(num_states * sizeof(unsigned int));
        }

        //Broadcast data
        MPI_Bcast(link_ids, N, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        MPI_Bcast(filedata_abs, N*num_states, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(filedata_rel, N*num_states, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(filedata_abs_dense, N*num_states, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(filedata_rel_dense, N*num_states, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(methods, N, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

        //Construct error data at each link
        for (unsigned int i = 0; i < N; i++)
        {
            Link *current = &system[i];
            if (assignments[i] == my_rank || getting[i])
            {
                current->my->error_data->abstol = calloc(num_states, sizeof(double));
                current->my->error_data->reltol = calloc(num_states, sizeof(double));
                current->my->error_data->abstol_dense = calloc(num_states, sizeof(double));
                current->my->error_data->reltol_dense = calloc(num_states, sizeof(double));
                current->my->error_data->facmax = error_data->facmax;
                current->my->error_data->facmin = error_data->facmin;
                current->my->error_data->fac = error_data->fac;

                for (unsigned int j = 0; j < num_states; j++)
                {
                    current->my->error_data->abstol[j] = filedata_abs[i*num_states + j];
                    current->my->error_data->reltol[j] = filedata_rel[i*num_states + j];
                    current->my->error_data->abstol_dense[j] = filedata_abs_dense[i*num_states + j];
                    current->my->error_data->reltol_dense[j] = filedata_rel_dense[i*num_states + j];
                }
                current->method = &rk_methods[rk_methods_idx[i]];
            }
        }
    }
    else
    {
        for (unsigned int i = 0; i < N; i++)
        {
            Link *current = &system[i];
            if (assignments[i] == my_rank || getting[i])
            {
                current->my->error_data = error_data;
                current->method = &rk_methods[globals->method];
            }
        }
    }

    return 0;
}

//Runs the init routine for the model. Also performs precalculations.
//Returns 0 if everything is ok, 1 if an error occurred.
int Initialize_Model(
    Link* system, unsigned int N,
    Link **my_sys, unsigned int my_N,
    int* assignments, short int* getting,
    GlobalVars* globals,
    AsynchModel* model, void* external)
{
    unsigned int i, j, max_dim = 0;
    int my_error_code = 0, error_code;

    for (i = 0; i < N; i++)
    {
        if (assignments[i] == my_rank || getting[i])
        {
            if (model)
            {
                model->routines(&system[i], globals->model_uid, system[i].method->exp_imp, system[i].has_dam, external);
                model->precalculations(&system[i], globals->global_params, system[i].params, system[i].has_dam, external);
            }
            else
            {
                InitRoutines(&system[i], globals->model_uid, system[i].method->exp_imp, system[i].has_dam, external);
                Precalculations(&system[i], globals->global_params, globals->num_global_params, system[i].params, globals->num_disk_params, globals->num_params, system[i].has_dam, globals->model_uid, external);
            }

            max_dim = (max_dim < system[i].dim) ? system[i].dim : max_dim;

            //Be sure the problem dimension and number of error tolerances are compatible
            //if (assignments[i] == my_rank)
            //{
            //    smallest_dim = system[i].my->error_data.abstol.dim;
            //    if (smallest_dim > system[i].my->error_data.reltol.dim)
            //        smallest_dim = system[i].my->error_data.reltol.dim;
            //    if (smallest_dim > system[i].my->error_data.abstol_dense.dim)
            //        smallest_dim = system[i].my->error_data.abstol_dense.dim;
            //    if (smallest_dim > system[i].my->error_data.reltol_dense.dim)
            //        smallest_dim = system[i].my->error_data.reltol_dense.dim;
            //    if (globals->min_error_tolerances > smallest_dim)
            //    {
            //        printf("[%i] Error: link id %u does not have enough error tolerances (got %u, expected %u)\n", my_rank, system[i].ID, smallest_dim, globals->min_error_tolerances);
            //        my_error_code = 1;
            //        break;
            //    }
            //}
        }
    }

    //Check if an error occurred
    MPI_Allreduce(&my_error_code, &error_code, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
    if (error_code)	return 1;

    //Make sure all procs know how large the problem is everywhere
    MPI_Allreduce(&max_dim, &(globals->max_dim), 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
    for (i = 0; i < N; i++)
        MPI_Bcast(&(system[i].dim), 1, MPI_UNSIGNED, assignments[i], MPI_COMM_WORLD);

    //Mix dense_indices with print_indices
    unsigned int loc, num_to_add;
    unsigned int* states_to_add = (unsigned int*)malloc(globals->num_states_for_printing * sizeof(unsigned int));	//!!!! Only mix if save_flag set !!!!
    for (loc = 0; loc < N; loc++)
    {
        if (assignments[loc] == my_rank || getting[loc])
        {
            Link* current = &system[loc];
            num_to_add = 0;

            for (i = 0; i < globals->num_states_for_printing; i++)
            {
                if (globals->print_indices[i] > current->dim)	continue;	//State is not present at this link
                for (j = 0; j < current->num_dense; j++)
                    if (globals->print_indices[i] == current->dense_indices[j])	break;
                if (j == current->num_dense)
                    states_to_add[num_to_add++] = globals->print_indices[i];
            }

            if (num_to_add)
            {
                current->dense_indices = (unsigned int*)realloc(current->dense_indices, (current->num_dense + num_to_add) * sizeof(unsigned int));
                for (i = 0; i < num_to_add; i++)
                    current->dense_indices[i + current->num_dense] = states_to_add[i];
                current->num_dense += num_to_add;
                merge_sort_1D(current->dense_indices, current->num_dense);
            }
        }
    }

    free(states_to_add);

    //Make sure all procs know the number of dense states at each link
    for (i = 0; i < N; i++)
        MPI_Bcast(&(system[i].num_dense), 1, MPI_UNSIGNED, assignments[i], MPI_COMM_WORLD);

    return 0;
}


static int Load_Initial_Conditions_Ini(
    Link *system, unsigned int N,
    int* assignments, short int* getting, const Lookup * const id_to_loc,
    GlobalVars* globals,
    ConnData* db_connections,
    AsynchModel* model,
    void* external)
{
    unsigned int id, loc, no_ini_start, diff_start = 0, dim;
    FILE* initdata = NULL;
    short int *who_needs = NULL;
    short int my_need;

    unsigned int max_dim = 0;
    for (unsigned int i = 0; i < N; i++)
        max_dim = max(max_dim, system[i].dim);

    //Proc 0 reads the file and sends the data to the other procs
    if (my_rank == 0)
    {
        initdata = fopen(globals->init_filename, "r");
        if (!initdata)
        {
            printf("Error: file %s not found for .ini file.\n", globals->init_filename);
            return 1;
        }
        if (CheckWinFormat(initdata))
        {
            printf("Error: File %s appears to be in Windows format. Try converting to unix format using 'dos2unix' at the command line.\n", globals->init_filename);
            fclose(initdata);
            return 1;
        }

        unsigned int n;
        fscanf(initdata, "%*i %u %lf", &n, &(globals->t_0));	//Read model type, number of links, init time

        if (n != N)
        {
            printf("Error: the number of links in %s differs from the number in the topology data. (Got %u, expected %u)\n", globals->init_filename, n, N);
            return 1;
        }

        //Broadcast initial time
        MPI_Bcast(&(globals->t_0), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        double *y_0 = malloc(max_dim * sizeof(double));

        //Read the .ini file
        who_needs = (short int*)malloc(np * sizeof(short int));
        for (unsigned int i = 0; i < N; i++)
        {
            //Send current location
            fscanf(initdata, "%u", &id);
            loc = find_link_by_idtoloc(id, id_to_loc, N);
            if (loc > N)
            {
                printf("Error: link id %u in initial condition file, but not in network.\n", id);
                return 1;
            }
            MPI_Bcast(&loc, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

            //See who needs info about this link.
            //0 means the proc doesn't need it, 1 means link is assigned to proc, 2 means the link is a ghost to the proc.
            my_need = (assignments[loc] == my_rank) ? 1 : ((getting[loc]) ? 2 : 0);
            MPI_Gather(&my_need, 1, MPI_SHORT, who_needs, 1, MPI_SHORT, 0, MPI_COMM_WORLD);
            if (my_need == 1)
            {
                no_ini_start = system[loc].no_ini_start;
                diff_start = system[loc].diff_start;
                dim = system[loc].dim;
            }
            else
            {
                MPI_Recv(&no_ini_start, 1, MPI_UNSIGNED, assignments[loc], 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&diff_start, 1, MPI_UNSIGNED, assignments[loc], 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&dim, 1, MPI_UNSIGNED, assignments[loc], 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            //Read init data
            for (unsigned int j = diff_start; j < no_ini_start; j++)
            {
                if (0 == fscanf(initdata, "%lf", y_0 + j))
                {
                    printf("Error: not enough states in .ini file.\n");
                    return 1;
                }
            }

            //Send data to assigned proc and getting proc
            //if (assignments[loc] == my_rank || getting[loc])
            if (system[loc].my)
            {
                if (model && model->initialize_eqs)
                    system[i].state = model->initialize_eqs(
                        globals->global_params, globals->num_global_params,
                        system[loc].params, globals->num_params,
                        y_0, system[loc].dim,
                        system[loc].user);
                else // load from ini
                    system[i].state = ReadInitData(
                        globals->global_params, globals->num_global_params,
                        system[loc].params, globals->num_params,
                        system[loc].qvs, system[loc].has_dam, y_0, system[loc].dim, globals->model_uid, diff_start, no_ini_start, system[loc].user, external);

                Init_List(&system[loc].my->list, globals->t_0, y_0, system[loc].dim, system[loc].num_dense, system[loc].method->num_stages, globals->iter_limit);
                system[loc].my->list.head->state = system[loc].state;
                system[loc].last_t = globals->t_0;
            }

            if (assignments[loc] != my_rank)
                MPI_Send(y_0 + diff_start, no_ini_start - diff_start, MPI_DOUBLE, assignments[loc], 2, MPI_COMM_WORLD);
            if (!(getting[loc]))
            {
                int j;
                for (j = 0; j < np; j++)
                    if (who_needs[j] == 2)
                        break;

                if (j < np)
                    MPI_Send(y_0 + diff_start, no_ini_start - diff_start, MPI_DOUBLE, (int)j, 2, MPI_COMM_WORLD);
            }
        }

        //Clean up
        free(y_0);
        fclose(initdata);
        free(who_needs);
    }
    else
    {
        double *y_0 = malloc(max_dim * sizeof(double));

        //Get initial time
        MPI_Bcast(&(globals->t_0), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        who_needs = (short int*)malloc(np * sizeof(short int));
        for (unsigned int i = 0; i < N; i++)
        {
            //Get link location
            MPI_Bcast(&loc, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

            //Is data needed for this link assigned at this proc?
            my_need = (assignments[loc] == my_rank) ? 1 : ((getting[loc]) ? 2 : 0);
            MPI_Gather(&my_need, 1, MPI_SHORT, who_needs, 1, MPI_SHORT, 0, MPI_COMM_WORLD);

            if (my_need)
            {
                no_ini_start = system[loc].no_ini_start;
                diff_start = system[loc].diff_start;
                dim = system[loc].dim;

                if (assignments[loc] == my_rank)
                {
                    MPI_Send(&no_ini_start, 1, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD);
                    MPI_Send(&diff_start, 1, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD);
                    MPI_Send(&dim, 1, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD);	//!!!! Actually, this might be available everywhere now !!!!
                }

                MPI_Recv(y_0 + diff_start, no_ini_start - diff_start, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                if (model && model->initialize_eqs)
                    system[i].state = model->initialize_eqs(
                        globals->global_params, globals->num_global_params,
                        system[loc].params, globals->num_params,
                        y_0, system[loc].dim,
                        system[loc].user);
                else // also from ini
                    system[i].state = ReadInitData(
                        globals->global_params, globals->num_global_params,
                        system[loc].params, globals->num_params,
                        system[loc].qvs, system[loc].has_dam, y_0, system[loc].dim, globals->model_uid, diff_start, no_ini_start, system[loc].user, external);
                
                Init_List(&system[loc].my->list, globals->t_0, y_0, system[loc].dim, system[loc].num_dense, system[loc].method->num_stages, globals->iter_limit);
                system[loc].my->list.head->state = system[loc].state;
                system[loc].last_t = globals->t_0;
            }
        }

        //Clean up
        free(y_0);
        free(who_needs);
    }

    return 0;
}

static int Load_Initial_Conditions_Uini(
    Link *system, unsigned int N,
    int* assignments, short int* getting, const Lookup * const id_to_loc,
    GlobalVars* globals,
    ConnData* db_connections,
    AsynchModel* model,
    void* external)
{
    FILE* initdata = NULL;
    short int *who_needs = NULL;
    double *y_0 = NULL;

    //Proc 0 reads the initial conds, and send them to the other procs
    if (my_rank == 0)
    {
        initdata = fopen(globals->init_filename, "r");
        if (!initdata)
        {
            printf("Error: file %s not found for .uini file.\n", globals->init_filename);
            return 1;
        }
        if (CheckWinFormat(initdata))
        {
            printf("Error: File %s appears to be in Windows format. Try converting to unix format using 'dos2unix' at the command line.\n", globals->init_filename);
            fclose(initdata);
            return 1;
        }

        fscanf(initdata, "%*i %lf", &(globals->t_0));	//Read model type, init time
    }

    //Broadcast the initial time
    MPI_Bcast(&(globals->t_0), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //Get number of values to read from disk (and error checking)
    unsigned int max_dim = 0;
    unsigned int no_ini_start, diff_start = 0;
    unsigned int i = 0;
    for (i = 0; i < N; i++)
    {
        if (assignments[i] == my_rank)
        {
            no_ini_start = system[i].no_ini_start;
            diff_start = system[i].diff_start;
            max_dim = max(max_dim, system[i].dim);
            break;
        }
    }
    for (; i < N; i++)
    {
        //if(assignments[i] == my_rank && y_0.dim != system[i].dim)
        if (assignments[i] == my_rank && (no_ini_start != system[i].no_ini_start || diff_start != system[i].diff_start))
        {
            printf("[%i]: Error: model type %u does not support .uini files (because a variable number of states must be specified for the initial conditions).\n", my_rank, globals->model_uid);
            return 1;
        }
    }

    //no_ini_start = system[loc].no_ini_start;
    //diff_start = system[loc].diff_start;
    double *y_0_backup = malloc((no_ini_start - diff_start) * sizeof(double));
    //y_0_backup->dim = no_ini_start - diff_start;
    //y_0_backup.ve = (double*) calloc(y_0_backup->dim,sizeof(double));

    if (my_rank == 0)
    {
        //for(i=diff_start;i<no_ini_start;i++)
        for (i = 0; i < no_ini_start - diff_start; i++)
        {
            if (fscanf(initdata, "%lf", y_0_backup + i) == 0)
            {
                printf("Error reading .uini file: Not enough initial states.\n");
                return 1;
            }
        }

        //Done with file, so close it
        fclose(initdata);
    }

    MPI_Bcast(y_0_backup, no_ini_start - diff_start, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //Store the data
    for (unsigned int i = 0; i < N; i++)
    {
        if (assignments[i] == my_rank || getting[i])
        {
            y_0 = realloc(y_0, system[i].dim * sizeof(double));
            for (unsigned int j = diff_start; j < no_ini_start; j++)
                y_0[j] = y_0_backup[j - diff_start];

            if (model && model->initialize_eqs)
                system[i].state = model->initialize_eqs(
                    globals->global_params, globals->num_global_params,
                    system[i].params, globals->num_params,
                    y_0, system[i].dim,
                    system[i].user);
            else // from uini
                system[i].state = ReadInitData(
                    globals->global_params, globals->num_global_params,
                    system[i].params, globals->num_params,
                    system[i].qvs, system[i].has_dam, y_0, system[i].dim, globals->model_uid, diff_start, no_ini_start, system[i].user, external);

            Init_List(&system[i].my->list, globals->t_0, y_0, system[i].dim, system[i].num_dense, system[i].method->num_stages, globals->iter_limit);
            system[i].my->list.head->state = system[i].state;
            system[i].last_t = globals->t_0;
            //v_copy(y_0_backup,y_0);
        }
    }

    //Clean up
    if (y_0)
        free(y_0);
    
    free(y_0_backup);

    return 0;
}

static int Load_Initial_Conditions_Rec(
    Link *system, unsigned int N,
    int* assignments, short int* getting, const Lookup * const id_to_loc,
    GlobalVars* globals,
    ConnData* db_connections,
    AsynchModel* model,
    void* external)
{
    unsigned int id, dim;
    FILE* initdata = NULL;
    short int *who_needs = NULL;
    short int my_need;
    double *y_0 = NULL;

    //Proc 0 reads the file and sends the data to the other procs
    if (my_rank == 0)
    {
        initdata = fopen(globals->init_filename, "r");
        if (!initdata)
        {
            printf("Error: file %s not found for .rec file.\n", globals->init_filename);
            return 1;
        }
        if (CheckWinFormat(initdata))
        {
            printf("Error: File %s appears to be in Windows format. Try converting to unix format using 'dos2unix' at the command line.\n", globals->init_filename);
            fclose(initdata);
            return 1;
        }

        unsigned int n;
        fscanf(initdata, "%*i %u %lf", &n, &(globals->t_0));	//Read model type, number of links, init time

        if (n != N)
        {
            printf("Error: the number of links in %s differs from the number in the topology data. (Got %u, expected %u)\n", globals->init_filename, n, N);
            return 1;
        }

        //Broadcast the initial time
        MPI_Bcast(&(globals->t_0), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        //Read the .rec file
        who_needs = (short int*)malloc(np * sizeof(short int));
        for (unsigned int i = 0; i < N; i++)
        {
            //Send current location
            fscanf(initdata, "%u", &id);
            unsigned int loc = find_link_by_idtoloc(id, id_to_loc, N);
            if (loc > N)
            {
                printf("Error: link id %u in initial condition file, but not in network.\n", id);
                return 1;
            }
            MPI_Bcast(&loc, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

            //See who needs info about this link.
            //0 means the proc doesn't need it, 1 means link is assigned to proc, 2 means the link is a ghost to the proc.
            my_need = (assignments[loc] == my_rank) ? 1 : ((getting[loc]) ? 2 : 0);
            MPI_Gather(&my_need, 1, MPI_SHORT, who_needs, 1, MPI_SHORT, 0, MPI_COMM_WORLD);
            if (my_need == 1)
                dim = system[loc].dim;
            else
                MPI_Recv(&dim, 1, MPI_UNSIGNED, assignments[loc], 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            //Read init data
            y_0 = realloc(y_0, dim * sizeof(double));
            for (unsigned int j = 0; j < dim; j++)
            {
                if (0 == fscanf(initdata, "%lf", y_0 + j))
                {
                    printf("Error: not enough states in .rec file.\n");
                    return 1;
                }
            }

            //Send data to assigned proc and getting proc
            if (assignments[loc] == my_rank || getting[loc])
            {
                if (system[loc].check_state)
                    system[loc].state = system[loc].check_state(y_0, system[loc].dim, globals->global_params, globals->num_global_params, system[loc].params, system[loc].num_params, system[loc].qvs, system[loc].state, system[loc].user);
                
                Init_List(&system[loc].my->list, globals->t_0, y_0, system[loc].dim, system[loc].num_dense, system[loc].method->num_stages, globals->iter_limit);
                system[loc].my->list.head->state = system[loc].state;
                system[loc].last_t = globals->t_0;
            }

            if (assignments[loc] != my_rank)
                MPI_Send(y_0, dim, MPI_DOUBLE, assignments[loc], 2, MPI_COMM_WORLD);

            if (!(getting[loc]))
            {
                int j;
                for (j = 0; j < np; j++)
                    if (who_needs[j] == 2)
                        break;

                if (j < np)
                    MPI_Send(y_0, dim, MPI_DOUBLE, (int)j, 2, MPI_COMM_WORLD);
            }
        }

        //Clean up
        fclose(initdata);
        free(who_needs);
    }
    else
    {
        //Get the initial time
        MPI_Bcast(&(globals->t_0), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        for (unsigned int i = 0; i < N; i++)
        {
            //Get link location
            unsigned int loc;
            MPI_Bcast(&loc, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

            //Is data needed for this link assigned at this proc?
            my_need = (assignments[loc] == my_rank) ? 1 : ((getting[loc]) ? 2 : 0);
            MPI_Gather(&my_need, 1, MPI_SHORT, who_needs, 1, MPI_SHORT, 0, MPI_COMM_WORLD);

            if (my_need)
            {
                dim = system[loc].dim;

                if (assignments[loc] == my_rank)
                    MPI_Send(&dim, 1, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD);

                y_0 = realloc(y_0, dim * sizeof(double));
                MPI_Recv(y_0, dim, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                if (system[loc].check_state)
                    system[loc].state = system[loc].check_state(y_0, dim, globals->global_params, globals->num_global_params, system[loc].params, system[loc].num_params, system[loc].qvs, system[loc].state, system[loc].user);
                
                Init_List(&system[loc].my->list, globals->t_0, y_0, dim, system[loc].num_dense, system[loc].method->num_stages, globals->iter_limit);
                system[loc].my->list.head->state = system[loc].state;
                system[loc].last_t = globals->t_0;
            }
        }
    }

    //Clean up
    if (y_0)
        free(y_0);

    return 0;
}


static int Load_Initial_Conditions_Dbc(
    Link *system, unsigned int N,
    int* assignments, short int* getting, const Lookup * const id_to_loc,
    GlobalVars* globals,
    ConnData* db_connections,
    AsynchModel* model,
    void* external)
{

#if defined(HAVE_POSTGRESQL)

    unsigned int loc;
    short int *who_needs = NULL;
    short int my_need;
    double *y_0 = NULL;
    PGresult *res;

    //!!!! Note: this assumes the database is like a .rec, with each state given. It also !!!!
    //!!!! assumes the same number of states at each link. !!!!

    //Set t_0 (I'm not really sure what else to do here...)
    globals->t_0 = 0.0;

    if (my_rank == 0)
    {
        //Download data
        if (ConnectPGDB(&db_connections[ASYNCH_DB_LOC_INIT]))
        {
            printf("Error connecting to database for init conditions.\n");
            return 1;
        }
        sprintf(db_connections[ASYNCH_DB_LOC_INIT].query, db_connections[ASYNCH_DB_LOC_INIT].queries[0], globals->init_timestamp);
        res = PQexec(db_connections[ASYNCH_DB_LOC_INIT].conn, db_connections[ASYNCH_DB_LOC_INIT].query);
        if (CheckResError(res, "downloading init data"))	return 1;

        if (PQntuples(res) != N)
        {
            printf("Error downloading init data. Got %i conditions, expected %u.\n", PQntuples(res), N);
            return 1;
        }

        //Get dim
        unsigned int dim = PQnfields(res) - 1;
        y_0 = realloc(y_0, dim * sizeof(double));
        MPI_Bcast(&dim, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        who_needs = (short int*)malloc(np * sizeof(short int));

        //Read data
        for (unsigned int i = 0; i < N; i++)
        {
            loc = find_link_by_idtoloc(atoi(PQgetvalue(res, i, 0)), id_to_loc, N);
            for (unsigned int j = 0; j < dim; j++)
                y_0[j] = atof(PQgetvalue(res, i, j + 1));
            MPI_Bcast(&loc, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

            //See who needs info about this link.
            //0 means the proc doesn't need it, 1 means link is assigned to proc, 2 means the link is a ghost to the proc.
            my_need = (assignments[loc] == my_rank) ? 1 : ((getting[loc]) ? 2 : 0);
            MPI_Gather(&my_need, 1, MPI_SHORT, who_needs, 1, MPI_SHORT, 0, MPI_COMM_WORLD);

            //Send the data
            if (assignments[loc] == my_rank || getting[loc])
            {
                if (system[loc].check_state != NULL)
                    system[loc].state = system[loc].check_state(y_0, dim, globals->global_params, globals->num_global_params, system[loc].params, system[loc].num_params, system[loc].qvs, system[loc].state, system[loc].user);

                Init_List(&system[loc].my->list, globals->t_0, y_0, dim, system[loc].num_dense, system[loc].method->num_stages, globals->iter_limit);
                system[loc].my->list.head->state = system[loc].state;
                system[loc].last_t = globals->t_0;
            }

            if (assignments[loc] != my_rank)
                MPI_Send(y_0, dim, MPI_DOUBLE, assignments[loc], 2, MPI_COMM_WORLD);
            if (!(getting[loc]))
            {
                int j;
                for (j = 0; j < np; j++)
                    if (who_needs[j] == 2)
                        break;

                if (j < np)
                    MPI_Send(y_0, dim, MPI_DOUBLE, (int)j, 2, MPI_COMM_WORLD);
            }
        }

        //Clean up
        PQclear(res);
        DisconnectPGDB(&db_connections[ASYNCH_DB_LOC_INIT]);
        free(y_0);
        free(who_needs);
    }
    else
    {
        //Get dim
        unsigned int dim;
        MPI_Bcast(&dim, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        y_0 = realloc(y_0, dim * sizeof(double));

        for (unsigned int i = 0; i < N; i++)
        {
            //Get link location
            MPI_Bcast(&loc, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

            //Is data needed for this link assigned at this proc?
            my_need = (assignments[loc] == my_rank) ? 1 : ((getting[loc]) ? 2 : 0);
            MPI_Gather(&my_need, 1, MPI_SHORT, who_needs, 1, MPI_SHORT, 0, MPI_COMM_WORLD);

            if (my_need)
            {
                MPI_Recv(y_0, dim, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                if (system[loc].check_state != NULL)
                    system[loc].state = system[loc].check_state(y_0, dim, globals->global_params, globals->num_global_params, system[loc].params, system[loc].num_params, system[loc].qvs, system[loc].state, system[loc].user);

                Init_List(&system[loc].my->list, globals->t_0, y_0, dim, system[loc].num_dense, system[loc].method->num_stages, globals->iter_limit);
                system[loc].my->list.head->state = system[loc].state;
                system[loc].last_t = globals->t_0;
            }
        }

        //Clean up
        free(y_0);
    }

#else //HAVE_POSTGRESQL

    if (my_rank == 0)	printf("Error: Asynch was build without PostgreSQL support.\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

#endif //HAVE_POSTGRESQL

    return 0;
}


static int Load_Initial_Conditions_H5(
    Link *system, unsigned int N,
    int* assignments, short int* getting, const Lookup * const id_to_loc,
    GlobalVars* globals,
    ConnData* db_connections,
    AsynchModel* model,
    void* external)
{
    unsigned char who_needs[ASYNCH_MAX_NUMBER_OF_PROCESS];
    unsigned char my_need;

    //Proc 0 reads the file and sends the data to the other procs
    if (my_rank == 0)
    {
        hid_t file_id = H5Fopen(globals->init_filename, H5F_ACC_RDONLY, H5P_DEFAULT);
        if (!file_id)
        {
            printf("Error: file %s not found for .h5 file.\n", globals->init_filename);
            return 1;
        }

        hid_t packet_file_id = H5PTopen(file_id, "/snapshot");
        if (packet_file_id < 0)
        {
            printf("Error: could not initialize h5 packet file %s.\n", globals->init_filename);
            return 2;
        }

        hsize_t num_packets;
        H5PTget_num_packets(packet_file_id, &num_packets);

        if (num_packets != N)
        {
            printf("Error: the number of links in %s differs from the number in the topology data. (Got %llu, expected %u)\n", globals->init_filename, num_packets, N);
            return 1;
        }

        //Assume that every links have the same dimension
        unsigned int dim = system[0].dim;

        //Read model type, init time
        unsigned short type;
        H5LTget_attribute_ushort(file_id, "/", "model", &type);

        if (type != globals->model_uid)
        {
            printf("Error: model type do no match. (Got %hu, expected %hu)\n", globals->model_uid, type);
            return 1;
        }

        //Broadcast the initial time
        //MPI_Bcast(&(globals->t_0), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        size_t line_size = sizeof(unsigned int) + dim * sizeof(double);
        char *buffer = malloc(line_size);

        //Read the .h5 file
        for (unsigned int i = 0; i < N; i++)
        {
            H5PTget_next(packet_file_id, 1, buffer);

            //Send current location
            unsigned int id = *((unsigned int *)buffer);
            unsigned int loc = find_link_by_idtoloc(id, id_to_loc, N);
            if (loc > N)
            {
                printf("Error: link id %u in initial condition file, but not in network.\n", id);
                return 1;
            }
            MPI_Bcast(&loc, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

            //See who needs info about this link.
            //0 means the proc doesn't need it, 1 means link is assigned to proc, 2 means the link is a ghost to the proc.
            my_need = (assignments[loc] == my_rank) ? 1 : ((getting[loc]) ? 2 : 0);
            MPI_Gather(&my_need, 1, MPI_UNSIGNED_CHAR, who_needs, 1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

            unsigned int dim;
            if (my_need == 1)
                dim = system[loc].dim;
            else
                MPI_Recv(&dim, 1, MPI_UNSIGNED, assignments[loc], 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            //Read init data
            double *y_0 = (double *)(buffer + sizeof(unsigned int));

            //Send data to assigned proc and getting proc
            if (assignments[loc] == my_rank || getting[loc])
            {
                if (system[loc].check_state)
                    system[loc].state = system[loc].check_state(y_0, dim, globals->global_params, globals->num_global_params, system[loc].params, system[loc].num_params, system[loc].qvs, system[loc].state, system[loc].user);

                Init_List(&system[loc].my->list, globals->t_0, y_0, dim, system[loc].num_dense, system[loc].method->num_stages, globals->iter_limit);
                system[loc].my->list.head->state = system[loc].state;
                system[loc].last_t = globals->t_0;
            }

            if (assignments[loc] != my_rank)
                MPI_Send(y_0, dim, MPI_DOUBLE, assignments[loc], 2, MPI_COMM_WORLD);
            if (!(getting[loc]))
            {
                int j;
                for (j = 0; j < np; j++)
                    if (who_needs[j] == 2)
                        break;

                if (j < np)
                    MPI_Send(y_0, dim, MPI_DOUBLE, (int)j, 2, MPI_COMM_WORLD);
            }
        }

        //Clean up
        free(buffer);
        H5PTclose(packet_file_id);
        H5Fclose(file_id);
    }
    else
    {
        ////Get the initial time
        //MPI_Bcast(&(globals->t_0), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        double *y_0 = NULL;

        for (unsigned int i = 0; i < N; i++)
        {
            //Get link location
            unsigned int loc;
            MPI_Bcast(&loc, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

            //Is data needed for this link assigned at this proc?
            my_need = (assignments[loc] == my_rank) ? 1 : ((getting[loc]) ? 2 : 0);
            MPI_Gather(&my_need, 1, MPI_UNSIGNED_CHAR, who_needs, 1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

            if (my_need)
            {
                unsigned int dim = system[loc].dim;
                
                if (assignments[loc] == my_rank)
                    MPI_Send(&dim, 1, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD);

                y_0 = realloc(y_0, dim * sizeof(double));
                MPI_Recv(y_0, dim, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                if (system[loc].check_state)
                    system[loc].state = system[loc].check_state(y_0, dim, globals->global_params, globals->num_global_params, system[loc].params, system[loc].num_params, system[loc].qvs, system[loc].state, system[loc].user);

                Init_List(&system[loc].my->list, globals->t_0, y_0, dim, system[loc].num_dense, system[loc].method->num_stages, globals->iter_limit);
                system[loc].my->list.head->state = system[loc].state;
                system[loc].last_t = globals->t_0;
            }
        }

        //Clean up
        free(y_0);
    }

    return 0;
}


//Loads the initial conditions.
//Initial state (0 = .ini, 1 = .uini, 2 = .rec, 3 = .dbc)
int Load_Initial_Conditions(
    Link *system, unsigned int N,
    int* assignments, short int* getting, const Lookup * const id_to_loc,
    GlobalVars* globals,
    ConnData* db_connections,
    AsynchModel* model,
    void* external)
{
    int res = 0;

    switch (globals->init_flag)
    {
        //.ini
    case 0:
        res = Load_Initial_Conditions_Ini(system, N, assignments, getting, id_to_loc, globals, db_connections, model, external);
        break;

        //.uini
    case 1:
        res = Load_Initial_Conditions_Uini(system, N, assignments, getting, id_to_loc, globals, db_connections, model, external);
        break;

        //.rec
    case 2:
        res = Load_Initial_Conditions_Rec(system, N, assignments, getting, id_to_loc, globals, db_connections, model, external);
        break;

        //.dbc
    case 3:
        res = Load_Initial_Conditions_Dbc(system, N, assignments, getting, id_to_loc, globals, db_connections, model, external);
        break;

        //.h5
    case 4:
        res = Load_Initial_Conditions_H5(system, N, assignments, getting, id_to_loc, globals, db_connections, model, external);
        break;

    default:
        printf("Error: invalid intial condition file type.\n");
        res = -1;
    }

    return res;
}


//Loads the forcing data specified in the global file.
int Load_Forcings(
    Link *system, unsigned int N,
    Link **my_sys, unsigned int my_N,
    int* assignments, short int* getting, unsigned int* res_list, unsigned int res_size, const Lookup * const id_to_loc,
    const GlobalVars * const globals,
    Forcing* forcings,
    ConnData* db_connections)
{
    // Create an MPI type for struct DataPoint
    const int nitems = 2;
    int          blocklengths[2] = { 1,1 };
    MPI_Datatype types[2] = { MPI_DOUBLE, MPI_FLOAT };
    MPI_Datatype mpi_datapoint_type;
    MPI_Aint     offsets[2];

    offsets[0] = offsetof(DataPoint, time);
    offsets[1] = offsetof(DataPoint, value);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_datapoint_type);
    MPI_Type_commit(&mpi_datapoint_type);

    //Reserve space for forcing data
    for (unsigned int i = 0; i < my_N; i++)
    {
        Link *current = my_sys[i];
        current->my->forcing_data = malloc(globals->num_forcings * sizeof(TimeSerie));
        current->my->forcing_values = calloc(globals->num_forcings, sizeof(double));
        current->my->forcing_change_times = calloc(globals->num_forcings, sizeof(double));
        current->my->forcing_indices = malloc(globals->num_forcings * sizeof(double));
    }

    //Setup forcings. Read uniform forcing data and open .str files. Also initialize rainfall from database.
    for (unsigned int l = 0; l < globals->num_forcings; l++)
    {
        forcings[l].maxtime = globals->t_0;
        forcings[l].iteration = 0;
        forcings[l].active = 1;

        //Go through each possible flag
        if (forcings[l].flag == 0)	//0 forcing
        {
            //Set routines
            forcings[l].GetPasses = &PassesOther;
            forcings[l].GetNextForcing = &NextForcingOther;

            //Setup buffers at each link
            for (unsigned int i = 0; i < N; i++)
            {
                if (assignments[i] == my_rank)
                {
                    system[i].my->forcing_data[l].data = NULL;
                    system[i].my->forcing_data[l].num_points = 0;
                    //v_set(system[i].forcing_values, l, 0.0;
                    system[i].my->forcing_change_times[l] = globals->maxtime + 1.0;
                }
            }
        }
        else if (forcings[l].flag == 1)	//Storm file
        {
            unsigned int limit, loc;
            FILE* forcingfile = NULL;

            //Set routines
            forcings[l].GetPasses = &PassesOther;
            forcings[l].GetNextForcing = &NextForcingOther;

            if (my_rank == 0)
            {
                //Open .str file
                forcingfile = fopen(forcings[l].filename, "r");
                if (!forcingfile)
                {
                    printf("Error: cannot open forcing file %s.\n", forcings[l].filename);
                    return 1;
                }
                if (CheckWinFormat(forcingfile))
                {
                    printf("Error: File %s appears to be in Windows format. Try converting to unix format using 'dos2unix' at the command line.\n", forcings[l].filename);
                    fclose(forcingfile);
                    return 1;
                }
                
                fscanf(forcingfile, "%u", &limit);
                if (limit != N && (!(globals->res_flag) || l != globals->res_forcing_idx))
                {
                    printf("Error: Number of links in .str file differs from number of links in network (%u vs %u).\n", limit, N);
                    return 1;
                }
            }

            //Get total number of links with data
            MPI_Bcast(&limit, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

            //Setup buffers at each link
            for (unsigned int i = 0; i < limit; i++)
            {
                DataPoint *buffer = NULL;
                unsigned int m;
                if (my_rank == 0)
                {
                    //Get location
                    unsigned int id;
                    fscanf(forcingfile, "%i", &id);
                    loc = find_link_by_idtoloc(id, id_to_loc, N);
                    if (loc >= N)
                    {
                        printf("Error: forcing data provided for link id %u in forcing %u, but link id is not in the network.\n", id, l);
                        return 1;
                    }

                    //Read values
                    fscanf(forcingfile, "%i", &m);
                    m++;		//Increase this by one to add a "ceiling" term
                    buffer = realloc(buffer, m * sizeof(DataPoint));
                    for (unsigned int j = 0; j < m - 1; j++)
                        fscanf(forcingfile, "%lf %f", &(buffer[j].time), &(buffer[j].value));
                    buffer[m - 1].time = globals->maxtime;
                    buffer[m - 1].value = -1.0f;

                    //Send data
                    MPI_Bcast(&loc, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
                    if (assignments[loc] != my_rank)
                    {
                        MPI_Send(&m, 1, MPI_UNSIGNED, assignments[loc], 1, MPI_COMM_WORLD);
                        MPI_Send(buffer, m, mpi_datapoint_type, assignments[loc], 1, MPI_COMM_WORLD);
                    }
                }
                else
                {
                    MPI_Bcast(&loc, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
                    if (assignments[loc] == my_rank)
                    {
                        MPI_Recv(&m, 1, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        buffer = realloc(buffer, m * sizeof(DataPoint));
                        MPI_Recv(buffer, m, mpi_datapoint_type, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }

                if (assignments[loc] == my_rank)
                {
                    //TimeSerie* forcing_data = malloc(sizeof(TimeSerie));
                    //system[loc].my->forcing_data[l] = forcing_data;
                    TimeSerie *forcing_data = &system[loc].my->forcing_data[l];

                    if (!(globals->res_flag) || !(l == globals->res_forcing_idx) || system[loc].has_res)
                    {
                        forcing_data->data = malloc(m * sizeof(DataPoint));
                        forcing_data->num_points = m;

                        //Read in the storm data for this link
                        for (unsigned int j = 0; j < m; j++)
                            forcing_data->data[j] = buffer[j];

                        double rainfall_buffer = forcing_data->data[0].value;
                        system[loc].my->forcing_values[l] = rainfall_buffer;
                        system[loc].my->forcing_indices[l] = 0;
                        unsigned int j;
                        for (j = 1; j < forcing_data->num_points; j++)
                        {
                            if (fabs(forcing_data->data[j].value - rainfall_buffer) > 1e-8)
                            {
                                system[loc].my->forcing_change_times[l] = forcing_data->data[j].time;
                                break;
                            }
                        }
                        if (j == forcing_data->num_points)
                        {
                            system[loc].my->forcing_change_times[l] = forcing_data->data[j - 1].time;
                            system[loc].my->forcing_indices[l] = j - 1;
                        }
                    }
                    else	//No reservoir here
                    {
                        unsigned int m = 2;	//Init value (assumed 0.0)

                        forcing_data->data = malloc(m * sizeof(DataPoint));
                        forcing_data->num_points = m;

                        forcing_data->data[0].time = globals->t_0;
                        forcing_data->data[0].value = 0.0;
                        forcing_data->data[1].time = globals->maxtime;
                        forcing_data->data[1].value = -1.0;

                        double rainfall_buffer = forcing_data->data[0].value;
                        system[loc].my->forcing_values[l] = rainfall_buffer;
                        system[loc].my->forcing_indices[l] = 0;
                        unsigned int j;
                        for (j = 1; j < forcing_data->num_points; j++)
                        {
                            if (fabs(forcing_data->data[j].value - rainfall_buffer) > 1e-8)
                            {
                                system[loc].my->forcing_change_times[l] = forcing_data->data[j].time;
                                break;
                            }
                        }
                        if (j == forcing_data->num_points)
                            system[loc].my->forcing_change_times[l] = forcing_data->data[j - 1].time;
                    }
                } //if(assigned to me)
            }

            //Allocate space for links without a reservoir
            if (globals->res_flag && l == globals->res_forcing_idx)
            {
                for (unsigned int i = 0; i < N; i++)
                {
                    if (assignments[i] == my_rank && !(system[i].my->forcing_data[l].data))
                    {
                        unsigned int m = 2;	//Init value (assumed 0.0)

                        TimeSerie* forcing_data = &system[loc].my->forcing_data[l];

                        forcing_data->data = malloc(m * sizeof(DataPoint));
                        forcing_data->num_points = m;

                        forcing_data->data[0].time = globals->t_0;
                        forcing_data->data[0].value = 0.0;
                        forcing_data->data[1].time = globals->maxtime + 3.0;
                        forcing_data->data[1].value = -1.0;

                        double rainfall_buffer = forcing_data->data[0].value;
                        system[i].my->forcing_values[l] = rainfall_buffer;
                        system[i].my->forcing_indices[l] = 0;
                        unsigned int j;
                        for (j = 1; j < forcing_data->num_points; j++)
                        {
                            if (fabs(forcing_data->data[j].value - rainfall_buffer) > 1e-8)
                            {
                                system[i].my->forcing_change_times[l] = forcing_data->data[j].time;
                                break;
                            }
                        }
                        if (j == forcing_data->num_points)
                            system[i].my->forcing_change_times[l] = forcing_data->data[j - 1].time;
                    }
                }
            }

            //Clean up
            if (forcingfile)
                fclose(forcingfile);
        }
        else if (forcings[l].flag == 2)	//Binary files
        {
            //Set routines
            forcings[l].GetPasses = &PassesBinaryFiles;
            forcings[l].GetNextForcing = &NextForcingBinaryFiles;

            //Setup buffers at each link
            //!!!! This should really be improved... !!!!
            for (unsigned int i = 0; i < N; i++)
                if (assignments[i] == my_rank)
                    system[i].my->forcing_data[l].data = NULL;
        }
        else if (forcings[l].flag == 5)	//Irregular Binary files
        {
            //Set routines
            forcings[l].GetPasses = &PassesIrregularBinaryFiles;
            forcings[l].GetNextForcing = &NextForcingIrregularBinaryFiles;

            //Setup buffers at each link
            //!!!! This should really be improved... it was just copied from above !!!!
            for (unsigned int i = 0; i < N; i++)
                if (assignments[i] == my_rank)
                    system[i].my->forcing_data[l].data = NULL;

            //Allocate memory as it is done with databases
            for (unsigned int i = 0; i < N; i++)
            {
                if (assignments[i] == my_rank)
                {
                    TimeSerie* forcing_data = &system[i].my->forcing_data[l];

                    if (!(globals->res_flag) || !(l == globals->res_forcing_idx) || system[i].has_res)
                    {
                        unsigned int m = forcings[l].increment + 4;	//+1 for init, +1 for ceiling, +2 for when init time doesn't line up with file_time
                        forcing_data->data = malloc(m * sizeof(DataPoint));
                        if (!forcing_data->data) {
                            printf("Reached memory limit on link %u out of %u. Aborting.\n", i, N);
                            exit(1);
                        }
                        forcing_data->num_points = m;

                        forcing_data->data[0].time = globals->t_0;
                        system[i].my->forcing_values[l] = 0.0;
                        system[i].my->forcing_change_times[l] = fabs(globals->t_0 + globals->maxtime) + 10.0;	//Just pick something away from t_0, and positive
                    }
                    else	//Reservoir, so allocate only a little memory
                    {
                        unsigned int m = 4;	//+1 for init, +1 for ceiling, +2 for when init time doesn't line up with file_time
                        forcing_data->data = malloc(m * sizeof(DataPoint));
                        forcing_data->num_points = m;

                        forcing_data->data[0].time = globals->t_0;
                        system[i].my->forcing_values[l] = 0.0;
                        system[i].my->forcing_change_times[l] = fabs(globals->t_0 + globals->maxtime) + 10.0;	//Just pick something away from t_0, and positive
                    }
                }
            }
        }
        else if (forcings[l].flag == 6)	//GZ binary files
        {
            //Set routines
            forcings[l].GetPasses = &PassesBinaryFiles;
            forcings[l].GetNextForcing = &NextForcingGZBinaryFiles;

            //Setup buffers at each link
            //!!!! This should really be improved... !!!!
            for (unsigned int i = 0; i < N; i++)
                if (assignments[i] == my_rank)
                    system[i].my->forcing_data[l].data = NULL;
        }
        else if (forcings[l].flag == 8)	//Grid cell based files
        {
            //Set routines
            forcings[l].GetPasses = &PassesBinaryFiles;
            forcings[l].GetNextForcing = &NextForcingGridCell;

            //Setup buffers at each link
            //!!!! This should really be improved... !!!!
            unsigned int i;
            for (i = 0; i < N; i++)
                if (assignments[i] == my_rank)
                    system[i].my->forcing_data[l].data = NULL;

            //Read the index file
            if (my_rank == 0)
            {
                FILE* forcingfile = fopen(forcings[l].filename, "r");
                if (!forcingfile)
                {
                    printf("Error: forcing file %s not found.\n", forcings[l].filename);
                    return 1;
                }

                char tempspace1[ASYNCH_MAX_PATH_LENGTH];
                char tempspace2[ASYNCH_MAX_PATH_LENGTH];
                forcings[l].lookup_filename = (char*)malloc(ASYNCH_MAX_PATH_LENGTH * sizeof(char));
                forcings[l].fileident = (char*)malloc(ASYNCH_MAX_PATH_LENGTH * sizeof(char));
                FindPath(forcings[l].filename, forcings[l].fileident);
                unsigned int length = (unsigned int)strlen(forcings[l].fileident);
                strcpy(forcings[l].lookup_filename, forcings[l].fileident);
                //fscanf(forcingfile,"%lf %lf %u %s %s",&(forcings[l].file_time),&(forcings[l].factor),&(forcings[l].num_cells),&(forcings[l].fileident[i]),&(forcings[l].lookup_filename[i]));
                fscanf(forcingfile, "%lf %lf %u %s %s", &(forcings[l].file_time), &(forcings[l].factor), &(forcings[l].num_cells), tempspace1, tempspace2);
                fclose(forcingfile);

                //Add path from index file?
                //if (tempspace1[0] == '/')	//Use absolute path
                    strcpy(forcings[l].fileident, tempspace1);
                //else		//Relative path
                //    strcpy(&(forcings[l].fileident[i]), tempspace1);
                //if (tempspace2[0] == '/')	//Use absolute path
                    strcpy(forcings[l].lookup_filename, tempspace2);
                //else		//Relative path
                //    strcpy(&(forcings[l].lookup_filename[i]), tempspace2);

                //Process the lookup file
                forcingfile = fopen(forcings[l].lookup_filename, "r");
                if (!forcingfile)
                {
                    printf("Error: lookup file %s not found.\n", forcings[l].lookup_filename);
                    return 1;
                }
                forcings[l].grid_to_linkid = malloc(forcings[l].num_cells * sizeof(unsigned int*));
                forcings[l].num_links_in_grid = calloc(forcings[l].num_cells, sizeof(unsigned int));

                while (!feof(forcingfile))	//Count the number of links in each grid cell
                {
                    unsigned int j;
                    if (!fscanf(forcingfile, "%*u %u", &j))	break;
                    (forcings[l].num_links_in_grid[j])++;
                }

                rewind(forcingfile);
                for (unsigned int i = 0; i < forcings[l].num_cells; i++)
                    forcings[l].grid_to_linkid[i] = (unsigned int*)malloc(forcings[l].num_links_in_grid[i] * sizeof(unsigned int));
                unsigned int* counters = (unsigned int*)calloc(forcings[l].num_cells, sizeof(unsigned int));

                while (!feof(forcingfile))	//Read in the link ids
                {
                    unsigned int id, j;
                    if (!fscanf(forcingfile, "%u %u", &id, &j))	break;
                    forcings[l].grid_to_linkid[j][counters[j]++] = id;
                }

                fclose(forcingfile);
                free(counters);

                //Check if a grid cell file actually exists
                unsigned int i;
                for (i = forcings[l].first_file; i < forcings[l].last_file; i++)
                {
                    sprintf(tempspace2, "%s%u", forcings[l].fileident, i);
                    forcingfile = fopen(tempspace2, "rb");
                    if (forcingfile)
                    {
                        fclose(forcingfile);
                        break;
                    }
                }
                if (i == forcings[l].last_file)
                    printf("Warning: No forcing files found at %s for forcing %u. Be sure this is the correct directory.\n", forcings[l].fileident, l);
            }

            //Broadcast data
            MPI_Bcast(&(forcings[l].file_time), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(&(forcings[l].factor), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(&(forcings[l].num_cells), 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
            if (my_rank != 0)
            {
                forcings[l].grid_to_linkid = (unsigned int**)malloc(forcings[l].num_cells * sizeof(unsigned int*));
                forcings[l].num_links_in_grid = (unsigned int*)malloc(forcings[l].num_cells * sizeof(unsigned int*));
            }
            MPI_Bcast(forcings[l].num_links_in_grid, forcings[l].num_cells, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
            if (my_rank != 0)
            {
                for (unsigned int i = 0; i < forcings[l].num_cells; i++)
                    forcings[l].grid_to_linkid[i] = (unsigned int*)malloc(forcings[l].num_links_in_grid[i] * sizeof(unsigned int));
            }
            for (unsigned int i = 0; i < forcings[l].num_cells; i++)
                MPI_Bcast(forcings[l].grid_to_linkid[i], forcings[l].num_links_in_grid[i], MPI_UNSIGNED, 0, MPI_COMM_WORLD);

            forcings[l].received = (char*)malloc(forcings[l].num_cells * sizeof(char));
            forcings[l].intensities = (float*)malloc(forcings[l].num_cells * sizeof(float));

            //Remove from grid_to_linkid all links not on this proc
            for (i = 0; i < forcings[l].num_cells; i++)
            {
                unsigned int drop = 0;
                for (unsigned int j = 0; j < forcings[l].num_links_in_grid[i]; j++)
                {
                    if (assignments[forcings[l].grid_to_linkid[i][j]] != my_rank)
                        drop++;
                    else
                        forcings[l].grid_to_linkid[i][j - drop] = forcings[l].grid_to_linkid[i][j];
                }
                forcings[l].num_links_in_grid[i] -= drop;
                forcings[l].grid_to_linkid[i] = (unsigned int*)realloc(forcings[l].grid_to_linkid[i], forcings[l].num_links_in_grid[i] * sizeof(unsigned int));
            }
        }
        else if (forcings[l].flag == 3 || forcings[l].flag == 9)	//Database (could be irregular timesteps)
        {
#if defined(HAVE_POSTGRESQL)
            PGresult *res;

            //Set routines
            if (forcings[l].flag == 3)
            {
                forcings[l].GetPasses = &PassesDatabase;
                forcings[l].GetNextForcing = &NextForcingDatabase;
            }
            else
            {
                forcings[l].GetPasses = &PassesDatabase_Irregular;
                forcings[l].GetNextForcing = &NextForcingDatabase_Irregular;
            }

            //Get good_timestamp
            unsigned int good_time, is_null;
            char* query = db_connections[ASYNCH_DB_LOC_FORCING_START + l].query;
            if (my_rank == 0)
            {
                ConnectPGDB(&db_connections[ASYNCH_DB_LOC_FORCING_START + l]);
                sprintf(query, db_connections[ASYNCH_DB_LOC_FORCING_START + l].queries[2], forcings[l].first_file);
                res = PQexec(db_connections[ASYNCH_DB_LOC_FORCING_START + l].conn, query);
                if (CheckResError(res, "querying a valid forcing time"))	return 1;
                is_null = PQgetisnull(res, 0, 0);
                if (is_null)	printf("Warning: forcing %u may have no data...\n", l);
                else		good_time = atoi(PQgetvalue(res, 0, 0));
                PQclear(res);
                DisconnectPGDB(&db_connections[ASYNCH_DB_LOC_FORCING_START + l]);
            }

            MPI_Bcast(&is_null, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if (is_null)
                forcings[l].good_timestamp = forcings[l].raindb_start_time;
            else
            {
                MPI_Bcast(&good_time, 1, MPI_INT, 0, MPI_COMM_WORLD);
                forcings[l].good_timestamp = good_time;
            }

            if (forcings[l].flag == 9)	//Download some extra data if the first_file is not on a timestamp with a forcing
                forcings[l].next_timestamp = forcings[l].good_timestamp;	//!!!! Could probably do something similar for flag 3 !!!!

                                                                            //Allocate space
            for (unsigned int i = 0; i < N; i++)
            {
                if (assignments[i] == my_rank)
                {
                     TimeSerie* forcing_data = &system[i].my->forcing_data[l];

                    if (!(globals->res_flag) || !(l == globals->res_forcing_idx) || system[i].has_res)
                    {
                        unsigned int m = forcings[l].increment + 4;	//+1 for init, +1 for ceiling, +2 for when init time doesn't line up with file_time
                        forcing_data->data = malloc(m * sizeof(DataPoint));
                        forcing_data->num_points = m;

                        if (i == 1000) {
                            printf("For link %u we are allocating %u.\n", i, m);
                        }

                        forcing_data->data[0].time = globals->t_0;
                        system[i].my->forcing_values[l] = 0.0;
                        system[i].my->forcing_change_times[l] = fabs(globals->t_0 + globals->maxtime) + 10.0;	//Just pick something away from t_0, and positive
                    }
                    else	//Reservoir, so allocate only a little memory
                    {
                        unsigned int m = 4;	//+1 for init, +1 for ceiling, +2 for when init time doesn't line up with file_time
                        forcing_data->data = malloc(m * sizeof(DataPoint));
                        forcing_data->num_points = m;
                        
                        forcing_data->data[0].time = globals->t_0;
                        system[i].my->forcing_values[l] = 0.0;
                        system[i].my->forcing_change_times[l] = fabs(globals->t_0 + globals->maxtime) + 10.0;	//Just pick something away from t_0, and positive
                    }
                }
            }

#else //HAVE_POSTGRESQL

            if (my_rank == 0)	printf("Error: Asynch was build without PostgreSQL support.\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

#endif //HAVE_POSTGRESQL

        }
        else if (forcings[l].flag == 4)	//.ustr
        {
            DataPoint *buffer = NULL;
            double univ_forcing_change_time[ASYNCH_MAX_NUM_FORCINGS];
            unsigned int m;

            //Set routines
            forcings[l].GetPasses = &PassesOther;
            forcings[l].GetNextForcing = &NextForcingOther;

            //Read uniform data
            if (my_rank == 0)
            {
                FILE* forcingfile = fopen(forcings[l].filename, "r");
                if (!forcingfile)
                {
                    printf("[%i]: Error: cannot open uniform forcing file %s.\n", my_rank, forcings[l].filename);
                    return 1;
                }
                if (CheckWinFormat(forcingfile))
                {
                    printf("Error: File %s appears to be in Windows format. Try converting to unix format using 'dos2unix' at the command line.\n", forcings[l].filename);
                    fclose(forcingfile);
                    return 1;
                }

                //Read the number of times for the rain data for this link                
                fscanf(forcingfile, "%i", &m);
                m++;		//Increase this by one to add a "ceiling" term
                buffer = realloc(buffer, m * sizeof(DataPoint));

                //Read the data
                for (unsigned int j = 0; j < m - 1; j++)
                    fscanf(forcingfile, "%lf %f", &buffer[j].time, &buffer[j].value);
                
                //TODO WTF?
                buffer[m - 1].time = globals->maxtime + 3.0;
                buffer[m - 1].value = -1.0;

                fclose(forcingfile);
            }

            MPI_Bcast(&m, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
            if (my_rank != 0)
                buffer = realloc(buffer, m * sizeof(DataPoint));
            MPI_Bcast(buffer, m, mpi_datapoint_type, 0, MPI_COMM_WORLD);

            //Create a global forcing object
            forcings[l].global_forcing.data = malloc(m * sizeof(DataPoint));
            forcings[l].global_forcing.num_points = m;

            for (unsigned int j = 0; j < m; j++)
                forcings[l].global_forcing.data[j] = buffer[j];

            double rainfall_buffer = forcings[l].global_forcing.data[0].value;
            unsigned int j;
            for (j = 1; j < forcings[l].global_forcing.num_points; j++)
            {
                if (rainfall_buffer != forcings[l].global_forcing.data[j].value)
                {
                    univ_forcing_change_time[l] = forcings[l].global_forcing.data[j].time;
                    break;
                }
            }
            if (j == forcings[l].global_forcing.num_points)
                univ_forcing_change_time[l] = forcings[l].global_forcing.data[j - 1].time;

            //Load the links
            for (unsigned int i = 0; i < N; i++)
            {
                if (assignments[i] == my_rank)
                {
                    system[i].my->forcing_data[l] = forcings[l].global_forcing;
                    system[i].my->forcing_values[l] = system[i].my->forcing_data[l].data[0].value;
                    system[i].my->forcing_change_times[l] = univ_forcing_change_time[l];
                    system[i].my->forcing_indices[l] = 0;
                }
            }

            if (buffer)
                free(buffer);
        }
        else if (forcings[l].flag == 7)	//monthly recurring files
        {
            double univ_forcing_change_time[ASYNCH_MAX_NUM_FORCINGS];

            //Set routines
            forcings[l].GetPasses = &PassesRecurring;
            forcings[l].GetNextForcing = &NextForcingRecurring;

            unsigned int num_months = 12;
            float *buffer = malloc(num_months * sizeof(float));

            //Read monthly file
            if (my_rank == 0)
            {
                FILE* forcingfile = fopen(forcings[l].filename, "r");
                if (!forcingfile)
                {
                    printf("Error: cannot open uniform forcing file %s.\n", forcings[l].filename);
                    return 1;
                }
                if (CheckWinFormat(forcingfile))
                {
                    printf("Error: File %s appears to be in Windows format. Try converting to unix format using 'dos2unix' at the command line.\n", forcings[l].filename);
                    fclose(forcingfile);
                    return 1;
                }

                //Read the data
                for (unsigned int j = 0; j < num_months; j++)
                    fscanf(forcingfile, "%f", &(buffer[j]));

                fclose(forcingfile);
            }

            MPI_Bcast(buffer, num_months, MPI_FLOAT, 0, MPI_COMM_WORLD);

            //Create a global forcing object
            forcings[l].global_forcing.data = malloc((num_months + 1) * sizeof(DataPoint));
            forcings[l].global_forcing.num_points = num_months + 1;
            for (unsigned int j = 0; j < num_months; j++)
                forcings[l].global_forcing.data[j].value = buffer[j];
            forcings[l].global_forcing.data[num_months].value = -1.0;

            //Find the starting month, and use the next month for the change time
            time_t start_time_t = forcings[l].first_file;
            struct tm *start_time = gmtime(&start_time_t);

            (start_time->tm_mon)++;
            start_time->tm_mday = 1;
            start_time->tm_hour = 0;
            start_time->tm_min = 0;
            start_time->tm_sec = 0;
            time_t next_month = mktime(start_time);
            univ_forcing_change_time[l] = difftime(next_month, forcings[l].first_file) / 60.0;

            //Load the links
            for (unsigned int i = 0; i < N; i++)
            {
                if (assignments[i] == my_rank)
                {
                    system[i].my->forcing_data[l] = forcings[l].global_forcing;
                    system[i].my->forcing_values[l] = system[i].my->forcing_data[l].data[0].value;
                    system[i].my->forcing_change_times[l] = univ_forcing_change_time[l];
                    system[i].my->forcing_indices[l] = 0;
                }
            }

            if (buffer)
                free(buffer);
        }
        else
        {
            if (my_rank == 0)	printf("Error: Bad forcing flag for forcing %u.\n", forcings[l].flag);
            return 1;
        }

    }

    //Find links with state forcing
    if (globals->res_flag)
    {
        //Download forcing data	!!!! Not sure if this is the way to go. Maybe separate function? !!!!
        forcings[globals->res_forcing_idx].passes = 1;
        unsigned int start_iteration = forcings[globals->res_forcing_idx].iteration;
        forcings[globals->res_forcing_idx].GetNextForcing(system, N, my_sys, my_N, assignments, globals, &forcings[globals->res_forcing_idx], db_connections, id_to_loc, globals->res_forcing_idx);
        forcings[globals->res_forcing_idx].iteration = start_iteration;	//Keep this the same

        //Setup init condition at links with forcing
        //Exchange_InitState_At_Forced(system, N, assignments, getting, res_list, res_size, id_to_loc, globals);
    }

    //Clean up
    MPI_Type_free(&mpi_datapoint_type);
    
    return 0;
}


//Reads any is_dam sources. Sets up appropriate buffers. Also sets flag for reservoirs.
int Load_Dams(
    Link *system, unsigned int N,
    Link **my_sys, unsigned int my_N,
    int* assignments, short int* getting, const Lookup * const id_to_loc, GlobalVars* globals, ErrorData* errors, ConnData* db_connections, unsigned int** res_list, unsigned int* res_size, unsigned int* my_res_size)
{

    unsigned int i, j, k, m, num_dams, id, size, num_values;
    Link* current;
    FILE* damfile = NULL;
    double* buffer = NULL;

    //Read is_dam file
    if (globals->uses_dam && globals->dam_flag == 1)	//.dam file
    {
		/*

        //VEC* buffer;
        size = globals->dam_params_size - the_model.num_params;
        buffer = (double*)malloc(size * sizeof(double));

        //Setup needs array. This is the procs that have getting set to 1.
        int* needs = NULL;
        if (my_rank == 0) needs = (int*)malloc(N * sizeof(int));
        int not_needed = -1;
        for (i = 0; i < N; i++)
        {
            if (getting[i])	MPI_Reduce(&my_rank, &(needs[i]), 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
            else		MPI_Reduce(&not_needed, &(needs[i]), 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        }

        if (my_rank == 0)
        {
            damfile = fopen(globals->dam_filename, "r");
            if (damfile == NULL)
            {
                printf("Error: Dam file %s not found.\n", globals->dam_filename);
                return 1;
            }
            if (CheckWinFormat(damfile))
            {
                printf("Error: File %s appears to be in Windows format. Try converting to unix format using 'dos2unix' at the command line.\n", globals->dam_filename);
                fclose(damfile);
                return 1;
            }

            fscanf(damfile, "%u", &num_dams);
            MPI_Bcast(&num_dams, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

            for (i = 0; i < num_dams; i++)
            {
                fscanf(damfile, "%u", &id);
                m = find_link_by_idtoloc(id, id_to_loc, N);
                if (m > N)
                {
                    printf("Error: Link id %u from .dam file %s not present in network.\n", id, globals->dam_filename);
                    return 1;
                }

                for (j = 0; j < size; j++)
                    fscanf(damfile, "%lf", &(buffer[j]));

                MPI_Bcast(&m, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

                //Either store the parameters or send them
                if (my_rank == assignments[m] || my_rank == needs[m])
                {
                    v_resize(&system[m].params, globals->dam_params_size);

                    for (j = globals->num_params; j < globals->dam_params_size; j++)
                        v_set(system[m].params, j, buffer[j - globals->num_params]);

                    system[m].is_dam = 1;
                }

                if (my_rank != assignments[m])
                    MPI_Send(buffer, size, MPI_DOUBLE, assignments[m], 1, MPI_COMM_WORLD);

                if (needs[m] > 0)
                    MPI_Send(buffer, size, MPI_DOUBLE, needs[m], 1, MPI_COMM_WORLD);
            }

            fclose(damfile);
            free(needs);
        }
        else
        {
            MPI_Bcast(&num_dams, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

            for (i = 0; i < num_dams; i++)
            {
                MPI_Bcast(&m, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
                if (my_rank == assignments[m] || getting[m])
                {
                    MPI_Recv(buffer, size, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                    v_resize(&system[m].params, globals->dam_params_size);

                    for (j = globals->num_params; j < globals->dam_params_size; j++)
                        v_set(system[m].params, j, buffer[j - globals->num_params]);
                    system[m].is_dam = 1;
                }
            }
        }

        //Clean up
        if (buffer)	free(buffer);

        //Set error tolerance
        if (globals->rkd_flag == 0)
        {
            globals->discont_tol = errors->abstol[0];
            if (globals->discont_tol < 1e-12 && my_rank == 0)
                printf("Warning: Discontinuity tolerance has been set to %e.\n", globals->discont_tol);
        }
        else
        {
            globals->discont_tol = 1e-8;
            if (my_rank == 0)	printf("Notice: Discontinuity tolerance has been set to %e.\n", globals->discont_tol);
        }
		*/
    }
    else if (globals->uses_dam && globals->dam_flag == 2)	//.qvs file
    {
        //Setup needs array. This is the procs that have getting set to 1.
        int* needs = NULL;
        if (my_rank == 0) needs = (int*)malloc(N * sizeof(int));
        int not_needed = -1;
        for (i = 0; i < N; i++)
        {
            if (getting[i])	MPI_Reduce(&my_rank, &(needs[i]), 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
            else		MPI_Reduce(&not_needed, &(needs[i]), 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        }

        if (my_rank == 0)
        {
            damfile = fopen(globals->dam_filename, "r");
            if (damfile == NULL)
            {
                printf("Error: Dam file %s not found.\n", globals->dam_filename);
                return 1;
            }
            if (CheckWinFormat(damfile))
            {
                printf("Error: File %s appears to be in Windows format. Try converting to unix format using 'dos2unix' at the command line.\n", globals->dam_filename);
                fclose(damfile);
                return 1;
            }

            fscanf(damfile, "%u", &num_dams);
            MPI_Bcast(&num_dams, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

            for (i = 0; i < num_dams; i++)
            {
                //Get link id and number of values
                fscanf(damfile, "%u %u", &id, &num_values);
                m = find_link_by_idtoloc(id, id_to_loc, N);
                if (m > N)
                {
                    printf("Error: Link id %u from .qvs file %s not present in network.\n", id, globals->dam_filename);
                    return 1;
                }

                //Read q vs s data
                buffer = realloc(buffer, 2 * num_values * sizeof(double));
                for (j = 0; j < num_values; j++)
                    fscanf(damfile, "%lf %lf", &(buffer[2 * j]), &(buffer[2 * j + 1]));

                //Send location
                MPI_Bcast(&m, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

                //Either store the parameters or send them
                if (my_rank == assignments[m] || getting[m])
                {
                    system[m].has_dam = 1;
                    system[m].qvs = (QVSData*)malloc(sizeof(QVSData));
                    system[m].qvs->n_values = num_values;
                    system[m].qvs->points_array = (double*)malloc(2 * num_values * sizeof(double));
                    system[m].qvs->points = (double**)malloc(num_values * sizeof(double*));
                    for (j = 0; j < num_values; j++)	system[m].qvs->points[j] = &(system[m].qvs->points_array[2 * j]);

                    for (j = 0; j < num_values; j++)
                    {
                        // system[m].qvs->points[j].time = buffer[2 * j];
                        // system[m].qvs->points[j].value = buffer[2 * j + 1];
						system[m].qvs->points[j][0] = buffer[2 * j];     // time
						system[m].qvs->points[j][1] = buffer[2 * j + 1]; // value
                    }
                }

                if (my_rank != assignments[m])
                {
                    MPI_Send(&num_values, 1, MPI_UNSIGNED, assignments[m], 1, MPI_COMM_WORLD);
                    MPI_Send(buffer, 2 * num_values, MPI_DOUBLE, assignments[m], 1, MPI_COMM_WORLD);
                }

                if (needs[m] > 0)
                {
                    MPI_Send(&num_values, 1, MPI_UNSIGNED, needs[m], 1, MPI_COMM_WORLD);
                    MPI_Send(buffer, 2 * num_values, MPI_DOUBLE, needs[m], 1, MPI_COMM_WORLD);
                }
            }

            fclose(damfile);
            free(needs);
        }
        else
        {
            MPI_Bcast(&num_dams, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

            for (i = 0; i < num_dams; i++)
            {
                MPI_Bcast(&m, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
                if (my_rank == assignments[m] || getting[m])
                {
                    MPI_Recv(&num_values, 1, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    buffer = (double*)realloc(buffer, 2 * num_values * sizeof(double));
                    MPI_Recv(buffer, 2 * num_values, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                    system[m].has_dam = 1;
                    system[m].qvs = (QVSData*)malloc(sizeof(QVSData));
                    system[m].qvs->n_values = num_values;
                    system[m].qvs->points_array = (double*)malloc(2 * num_values * sizeof(double));
                    system[m].qvs->points = (double**)malloc(num_values * sizeof(double*));
                    for (j = 0; j < num_values; j++)	system[m].qvs->points[j] = &(system[m].qvs->points_array[2 * j]);

                    for (j = 0; j < num_values; j++)
                    {
                        // system[m].qvs->points[j].time = buffer[2 * j];
                        // system[m].qvs->points[j].value = buffer[2 * j + 1];
						system[m].qvs->points[j][0] = buffer[2 * j];      // time
						system[m].qvs->points[j][1] = buffer[2 * j + 1];  // value
                    }
                }
            }
        }

        //Clean up
        if (buffer)	free(buffer);

        //Set error tolerance
        if (globals->rkd_flag == 0)
        {
            globals->discont_tol = errors->abstol[0];
            if (globals->discont_tol < 1e-12 && my_rank == 0)
                printf("Warning: Discontinuity tolerance has been set to %e.\n", globals->discont_tol);
        }
        else
        {
            globals->discont_tol = 1e-8;
            if (my_rank == 0)	printf("Notice: Discontinuity tolerance has been set to %e.\n", globals->discont_tol);
        }
    }
    else if (globals->uses_dam && globals->dam_flag == 3)	//database connection
    {
#if defined(HAVE_POSTGRESQL)
        PGresult *res;
        Link* current;
        unsigned int num_pts, curr_loc;
        num_dams = 0;
        short int procs_sending_to[ASYNCH_MAX_NUMBER_OF_PROCESS], mine;
        double* array_holder;

        //Connect to the database and download data
        if (my_rank == 0)
        {
            ConnectPGDB(&db_connections[ASYNCH_DB_LOC_QVS]);
            res = PQexec(db_connections[ASYNCH_DB_LOC_QVS].conn, db_connections[ASYNCH_DB_LOC_QVS].queries[0]);
            if (CheckResError(res, "querying qvs relations"))	return 1;
            num_pts = PQntuples(res);

            i = 0;
            while (i < num_pts)
            {
                //Get link id and location
                id = atoi(PQgetvalue(res, i, 0));
                curr_loc = find_link_by_idtoloc(id, id_to_loc, N);

                //Count number of qvs values for current
                for (j = i; j < num_pts && id == atoi(PQgetvalue(res, j, 0)); j++);
                num_values = j - i;

                //Extract the data points
                array_holder = (double*)malloc(2 * num_values * sizeof(double));
                for (j = 0; j < num_values; j++)
                {
                    array_holder[2 * j] = atof(PQgetvalue(res, i + j, 1));
                    array_holder[2 * j + 1] = atof(PQgetvalue(res, i + j, 2));
                }

                //Check for an error real quick
                for (j = 1; j < num_values; j++)
                {
                    if (array_holder[2 * (j - 1)] > array_holder[2 * j] || array_holder[2 * (j - 1) + 1] > array_holder[2 * j + 1])
                    {
                        printf("[%i]: Bad storage or discharge values found at link id %u. Check that the data is sorted correctly. (%u)\n", my_rank, id, j);
                        //break;
                        free(array_holder);
                        return 1;
                    }
                }

                if (curr_loc < N)
                {
                    //Tell everyone what link has the current is_dam
                    MPI_Bcast(&curr_loc, 1, MPI_INT, 0, MPI_COMM_WORLD);
                    mine = (my_rank == assignments[curr_loc] || getting[curr_loc]);
                    MPI_Gather(&mine, 1, MPI_SHORT, procs_sending_to, 1, MPI_SHORT, 0, MPI_COMM_WORLD);

                    //Send the data to whoever needs it
                    for (int j = 1; j < np; j++)
                    {
                        if (procs_sending_to[j])
                        {
                            MPI_Send(&num_values, 1, MPI_UNSIGNED, j, 0, MPI_COMM_WORLD);
                            MPI_Send(array_holder, 2 * num_values, MPI_DOUBLE, j, 0, MPI_COMM_WORLD);
                        }
                    }

                    //Check if proc 0 needs the data
                    if (mine)
                    {
                        current = &system[curr_loc];
                        current->has_dam = 1;
                        current->qvs = (QVSData*)malloc(sizeof(QVSData));
                        current->qvs->n_values = num_values;
                        current->qvs->points_array = array_holder;
                        current->qvs->points = (double**)malloc(num_values * sizeof(double*));
                        for (j = 0; j < num_values; j++)	current->qvs->points[j] = &(current->qvs->points_array[2 * j]);
                    }
                    else
                        free(array_holder);
                }
                else
                    free(array_holder);

                i += num_values;
            }

            //Finish up
            PQclear(res);
            DisconnectPGDB(&db_connections[ASYNCH_DB_LOC_QVS]);
            curr_loc = -1;
            MPI_Bcast(&curr_loc, 1, MPI_INT, 0, MPI_COMM_WORLD);
        }
        else	//Receive is_dam data
        {
            curr_loc = 0;
            MPI_Bcast(&curr_loc, 1, MPI_INT, 0, MPI_COMM_WORLD);

            while ((int)curr_loc != -1)
            {
                //Check if I need the data for this link
                mine = (my_rank == assignments[curr_loc] || getting[curr_loc]);
                MPI_Gather(&mine, 1, MPI_SHORT, procs_sending_to, 1, MPI_SHORT, 0, MPI_COMM_WORLD);

                //Receive data
                if (mine)
                {
                    MPI_Recv(&num_values, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    array_holder = (double*)malloc(2 * num_values * sizeof(double));
                    MPI_Recv(array_holder, 2 * num_values, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    current = &system[curr_loc];
                    current->has_dam = 1;
                    current->qvs = (QVSData*)malloc(sizeof(QVSData));
                    current->qvs->n_values = num_values;
                    current->qvs->points_array = array_holder;
                    current->qvs->points = (double**)malloc(num_values * sizeof(double*));
                    for (j = 0; j < num_values; j++)	current->qvs->points[j] = &(current->qvs->points_array[2 * j]);
                }

                //Check next signal
                MPI_Bcast(&curr_loc, 1, MPI_INT, 0, MPI_COMM_WORLD);
            }
        }

        //Set error tolerance
        if (globals->rkd_flag == 0)
        {
            globals->discont_tol = errors->abstol[0];
            if (globals->discont_tol < 1e-12 && my_rank == 0)
                printf("Warning: Discontinuity tolerance has been set to %e.\n", globals->discont_tol);
        }
        else
        {
            globals->discont_tol = 1e-8;
            if (my_rank == 0)	printf("Notice: Discontinuity tolerance has been set to %e.\n", globals->discont_tol);
        }

#else //HAVE_POSTGRESQL

        if (my_rank == 0)	printf("Error: Asynch was build without PostgreSQL support.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

#endif //HAVE_POSTGRESQL
    }
    else	//Some other type of discontinuity (or none at all)
    {
        //Set error tolerance
        if (globals->rkd_flag == 0)
        {
            globals->discont_tol = errors->abstol[0];
            if (globals->discont_tol < 1e-12 && my_rank == 0)
                printf("Warning: Discontinuity tolerance has been set to %e.\n", globals->discont_tol);
        }
        else
        {
            globals->discont_tol = 1e-8;
            //if(my_rank == 0)	printf("Notice: Discontinuity tolerance has been set to %e.\n",globals->discont_tol);
        }
    }

    //Check where state forcings should be set
    *my_res_size = 0;
    if (globals->res_flag)
    {
        if (Create_SAV_Data(globals->rsv_filename, system, N, res_list, res_size, &db_connections[ASYNCH_DB_LOC_RSV], globals->res_flag))
            return 1;

        //Setup links with forcing
        for (j = 0; j < *res_size; j++)
        {
            k = find_link_by_idtoloc((*res_list)[j], id_to_loc, N);

            if (k == N + 1)
            {
                if (my_rank == 0)
                    printf("Warning: Reservoir forcing requested for a link NOT in the network (link id = %u).\n", (*res_list)[j]);

                //Shift ids in peak save list to spot k, and decrement peaksave_size
                for (k = j + 1; k < *res_size; k++)
                    (*res_list)[k - 1] = (*res_list)[k];
                j--;
                (*res_size)--;
            }
            else if (assignments[k] == my_rank || getting[k])
            {
                (*my_res_size)++;
                current = &system[k];
                current->has_res = 1;
            }
        }
    }
    else
    {
        (*res_list) = NULL;
        *res_size = 0;
    }

    return 0;
}

//Calculates an estimate for the initial stepsizes at each link.
//Returns 0 if initial conditions set. 1 if a forcing is not ready.
int CalculateInitialStepSizes(
    Link **my_sys, unsigned int my_N,
    const GlobalVars * const globals,
    Workspace* workspace,
    short int watch_forcings)
{
    unsigned int i, j;

    for (i = 0; i < my_N; i++)
    {
        my_sys[i]->h = InitialStepSize(globals->t_0, my_sys[i], globals, workspace);
        assert(my_sys[i]->h > 0.0);
    }

    if (watch_forcings)
    {
        for (i = 0; i < my_N; i++)
        {
            for (j = 0; j < globals->num_forcings; j++)
                my_sys[i]->h = min(my_sys[i]->h, my_sys[i]->my->forcing_change_times[j] - my_sys[i]->last_t);                
        }
    }

    return 0;
}

//Reads the link ids for saving data.
//Sets in GlobalVars: save_list, save_size, peaksave_list, peaksave_size.
//Also sets: my_save_size.
int BuildSaveLists(
    Link* system, unsigned int N,
    Link **my_sys, unsigned int my_N,
    int* assignments, const Lookup * const id_to_loc,
    GlobalVars* globals, unsigned int** save_list, unsigned int* save_size, unsigned int* my_save_size, unsigned int** peaksave_list, unsigned int* peaksave_size, unsigned int* my_peaksave_size, ConnData* db_connections)
{
    unsigned int j, k;
    Link* current;

    if (!save_list || !peaksave_list)
    {
        if (my_rank == 0)	printf("Error: id list for output data cannot be null.\n");
        return 1;
    }

    //For time series
    *my_save_size = 0;
    if (globals->hydros_loc_flag)
    {
        if (Create_SAV_Data(globals->hydrosave_filename, system, N, save_list, save_size, &db_connections[ASYNCH_DB_LOC_HYDROSAVE], globals->hydrosave_flag))
            return 1;

        for (j = 0; j < *save_size; j++)
        {
            k = find_link_by_idtoloc((*save_list)[j], id_to_loc, N);

            if (k == N + 1)
            {
                if (my_rank == 0)
                    printf("Warning: Time series output requested for a link NOT in the network (link id = %u).\n", (*save_list)[j]);

                //Shift ids in save list to spot k, and decrement save_size
                for (k = j + 1; k < *save_size; k++)
                    (*save_list)[k - 1] = (*save_list)[k];
                j--;
                (*save_size)--;
            }
            else if (assignments[k] == my_rank)
            {
                (*my_save_size)++;
                current = &system[k];
                current->save_flag = 1;
                current->next_save = globals->t_0;
                if (globals->print_time > 0.0)
                    current->print_time = globals->print_time;
                else
                    current->print_time = sqrt(current->params[globals->area_idx] * 0.1);
            }
        }
    }
    else
    {
        (*save_list) = NULL;
        *save_size = 0;
    }

    //For peakflows
    *my_peaksave_size = 0;
    if (globals->peaks_loc_flag)
    {
        if (Create_SAV_Data(globals->peaksave_filename, system, N, peaksave_list, peaksave_size, &db_connections[ASYNCH_DB_LOC_PEAKSAVE], globals->peaksave_flag))
            return 1;

        for (j = 0; j < *peaksave_size; j++)
        {
            k = find_link_by_idtoloc((*peaksave_list)[j], id_to_loc, N);

            if (k == N + 1)
            {
                if (my_rank == 0)
                    printf("Warning: Peakflow output requested for a link NOT in the network (link id = %u).\n", (*peaksave_list)[j]);

                //Shift ids in peak save list to spot k, and decrement peaksave_size
                for (k = j + 1; k < *peaksave_size; k++)
                    (*peaksave_list)[k - 1] = (*peaksave_list)[k];
                j--;
                (*peaksave_size)--;
            }
            else if (assignments[k] == my_rank)
            {
                (*my_peaksave_size)++;
                current = &system[k];
                current->peak_flag = 1;
            }
        }
    }
    else
    {
        (*peaksave_list) = NULL;
        *peaksave_size = 0;
    }

    return 0;
}

//Finalizes the system for calculations. Basically, lots of small miscellaneous initializations are made.
int FinalizeSystem(
    Link* system, unsigned int N,
    Link **my_sys, unsigned int my_N,
    int* assignments, short int* getting, const Lookup * const id_to_loc, TransData* my_data, GlobalVars* globals, ConnData* db_connections, Workspace* workspace)
{
    unsigned int i, j;
    int ii;
    Link* current;

    //Initialize workspace
    Create_Workspace(workspace, globals->max_dim, globals->max_rk_stages, globals->max_parents);

    //Need space for nodes, number of iterations, discontinuities
    //Data: ( size(double)*(max_rk_stages*max_dim + max_dim + time)*# steps to transfer + size(int)(for stage)*# steps to transfer + size(int)*(location + # steps to transfer) ) * # of sending links
    //Upstream: + size(int) * (location + # of iterations) * # of receiving links
    //Discontinuities: + (size(int) + size(double)*discont_size + size(int)*discont_size) * # of sending links
    //unsigned int bytes1 = ( (sizeof(double)*(globals->max_rk_stages*2 + 2 + 1) + sizeof(int) )*globals->max_transfer_steps + sizeof(int)*2);
    unsigned int bytes2 = 2 * sizeof(int);
    unsigned int bytes3 = sizeof(int) + (sizeof(int) + sizeof(double))*globals->discont_size;
    for (ii = 0; ii < np; ii++)
    {
        my_data->send_buffer_size[ii] = bytes2 * my_data->receive_size[ii] + bytes3 * my_data->send_size[ii];
        for (j = 0; j < my_data->send_size[ii]; j++)
        {
            current = my_data->send_data[ii][j];
            my_data->send_buffer_size[ii] += (sizeof(double)*(current->method->num_stages * current->num_dense + current->dim + 1) + sizeof(int))*globals->max_transfer_steps + sizeof(int) * 2;
        }

        my_data->receive_buffer_size[ii] = bytes2 * my_data->send_size[ii] + bytes3 * my_data->receive_size[ii];
        for (j = 0; j < my_data->receive_size[ii]; j++)
        {
            current = my_data->receive_data[ii][j];
            my_data->receive_buffer_size[ii] += (sizeof(double)*(current->method->num_stages * current->num_dense + current->dim + 1) + sizeof(int))*globals->max_transfer_steps + sizeof(int) * 2;
        }

        //(*my_data)->send_buffer_size[ii] = bytes1 * (*my_data)->send_size[ii] + bytes2 * (*my_data)->receive_size[ii] + bytes3 * (*my_data)->send_size[ii];
        //(*my_data)->receive_buffer_size[ii] = bytes1 * (*my_data)->receive_size[ii] + bytes2*(*my_data)->send_size[ii] + bytes3 * (*my_data)->receive_size[ii];

        if (my_data->send_buffer_size[ii])	my_data->send_buffer[ii] = (char*)malloc(my_data->send_buffer_size[ii] * sizeof(char));
        else					my_data->send_buffer[ii] = NULL;

        if (my_data->receive_buffer_size[ii])	my_data->receive_buffer[ii] = (char*)malloc(my_data->receive_buffer_size[ii] * sizeof(char));
        else					my_data->receive_buffer[ii] = NULL;
    }

    //Do initializations
    for (i = 0; i < N; i++)
    {
        if (assignments[i] == my_rank || getting[i])
        {
            //Discontinuity information
            if (system[i].num_parents)
                system[i].discont = (double*)malloc(globals->discont_size * sizeof(double));
            if (system[i].child && my_rank != assignments[system[i].child->location])
            {
                system[i].discont_send = (double*)malloc(globals->discont_size * sizeof(double));
                system[i].discont_order_send = (unsigned int*)malloc(globals->discont_size * sizeof(unsigned int));
                system[i].discont_send_count = 0;
            }

            //Setup most of the remaining data
            //system[i].last_t = globals->t_0;
            //system[i].print_time = globals->t_0;
            system[i].current_iterations = 1;
            if (!(system[i].save_flag))	system[i].next_save = globals->t_0 - 1.0;
            system[i].peak_time = globals->t_0;
            //system[i].save_flag = 0;
            //system[i].peak_flag = 0;
            system[i].peak_value = malloc(system[i].dim * sizeof(double));
            dcopy(system[i].my->list.head->y_approx, system[i].peak_value, 0, system[i].dim);
            if (system[i].num_parents)	system[i].ready = 0;
            else				system[i].ready = 1;
            system[i].iters_removed = 0;
            system[i].steps_on_diff_proc = 1;	//Note: This won't be used if link isn't stored on another proc
            system[i].rejected = 0;

#if defined (ASYNCH_HAVE_IMPLICIT_SOLVER)
            system[i].last_eta = 1e10;
            system[i].compute_J = 1;
            system[i].compute_LU = 1;

            if (system[i].method->exp_imp == 1)
            {
                system[i].JMatrix = m_get(system[i].dim, system[i].dim);
                system[i].CoefMat = m_get(globals->max_s*system[i].dim, globals->max_s*system[i].dim);
                system[i].Z_i = malloc(globals->max_s * sizeof(VEC*));
                for (j = 0; j < globals->max_s; j++)
                    system[i].Z_i[j] = v_init(system[i].dim);
                system[i].sol_diff = v_init(system[i].dim);
                system[i].h_old = -1.0;
                system[i].value_old = -1.0;
            }
            else
            {
                system[i].JMatrix = m_get(0, 0);
                system[i].CoefMat = m_get(0, 0);
                system[i].Z_i = NULL;
                system[i].sol_diff = v_init(0);
            }
#endif
        }
        else	//If link not assigned to this process, and not receiving any information from this link.
        {
            system[i].method = NULL;
            system[i].peak_value = NULL;
            system[i].params = NULL;

#if defined (ASYNCH_HAVE_IMPLICIT_SOLVER)
            system[i].JMatrix = m_get(0, 0);
            system[i].CoefMat = m_get(0, 0);
            system[i].Z_i = NULL;
            system[i].sol_diff = v_init(0);
#endif
        }
    }

    return 0;
}


//Reads in the contents of a .sav file.
//char filename[]: The filename of the .sav file.
//int N: The number of links in the river system.
//int** save_list (set by this method): Array of link ids where data will be written.
//int* size (set by this method): Will be the number of links for which data must be written to disk (number of links in the .sav file).
//Returns 1 if an error occured. 0 otherwise.
int Create_SAV_Data(char filename[], Link* sys, unsigned int N, unsigned int** save_list, unsigned int* size, ConnData *conninfo, unsigned short int flag)
{
    unsigned int i, id, error = 0, list_size = N;
    FILE* save_file = NULL;
    *size = 0;
    *save_list = NULL;

    //.sav file
    if (flag == 1)
    {
        if (my_rank == 0)
        {
            //Open save file
            save_file = fopen(filename, "r");
            if (!save_file)
            {
                printf("Error opening .sav file %s.\n", filename);
                error = 1;
            }
            else
            {
                if (CheckWinFormat(save_file))
                {
                    printf("Error: File %s appears to be in Windows format. Try converting to unix format using 'dos2unix' at the command line.\n", filename);
                    fclose(save_file);
                    return 1;
                }

                *save_list = (unsigned int*)malloc(list_size * sizeof(unsigned int));

                //Read the save file
                while (fscanf(save_file, "%u", &id) != EOF)
                {
                    (*save_list)[*size] = id;
                    (*size)++;
                    if ((*size) == list_size)
                    {
                        list_size *= 2;
                        *save_list = (unsigned int*)realloc(*save_list, list_size * sizeof(unsigned int));
                    }
                }

                //Close save file
                fclose(save_file);
            }
        }

        //Send data
        MPI_Bcast(&error, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        if (error)	return 1;
        MPI_Bcast(size, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        *save_list = realloc(*save_list, *size * sizeof(unsigned int));
        MPI_Bcast(*save_list, *size, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    }
    else if (flag == 2)	//Grab from database
    {
#if defined(HAVE_POSTGRESQL)

        //char* query = conninfo.query;
        PGresult *res;

        if (my_rank == 0)
        {
            ConnectPGDB(conninfo);
            sprintf(conninfo->query, "%s", conninfo->queries[0]);
            res = PQexec(conninfo->conn, conninfo->query);
            error = CheckResError(res, "locating links with sensors");

            if (!error)
            {
                *size = PQntuples(res);
                *save_list = malloc(*size * sizeof(unsigned int));
                for (i = 0; i < *size; i++)	(*save_list)[i] = atoi(PQgetvalue(res, i, 0));
            }

            PQclear(res);
            DisconnectPGDB(conninfo);
        }

        MPI_Bcast(&error, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        if (error)	return 1;
        MPI_Bcast(size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (my_rank != 0)	*save_list = malloc(*size * sizeof(unsigned int));
        MPI_Bcast(*save_list, *size, MPI_INT, 0, MPI_COMM_WORLD);

#else //HAVE_POSTGRESQL

        if (my_rank == 0)	printf("Error: Asynch was build without PostgreSQL support.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

#endif //HAVE_POSTGRESQL
    }
    else if (flag == 3)	//All links
    {
        *size = N;
        *save_list = malloc(*size * sizeof(unsigned int));
        for (i = 0; i < N; i++)	(*save_list)[i] = sys[i].ID;
    }

    return 0;
}


//Returns 1 if the ascii file is in Windows format. 0 otherwise.
int CheckWinFormat(FILE* file)
{
    char buffer[256];
    fpos_t pos;
    fgetpos(file, &pos);

    if (!file)	return 0;
    fgets(buffer, 256, file);
    if (buffer[strlen(buffer) - 2] == 13)	//CR
    {
        fsetpos(file, &pos);
        return 1;
    }
    fsetpos(file, &pos);
    return 0;
}

//Finds the path in filename. The path is copied into variable path.
//Returns 0 if the path is is extracted successfully.
//Returns 1 if filename does not contain a path (path variable set to empty).
//Returns 2 if filename is just a path (the path variable is still set).
int FindPath(char* filename, char* path)
{
    size_t len = strlen(filename);
    char holder;

    if (len == 0)
    {
        path[0] = '\0';
        return 1;
    }

    if (filename[len - 1] == '/')
    {
        strcpy(path, filename);
        return 2;
    }

    for (int i = (int)(len - 2); i >= 0; i--)
    {
        if (filename[i] == '/')
        {
            holder = filename[i + 1];
            filename[i + 1] = '\0';
            strcpy(path, filename);
            filename[i + 1] = holder;
            return 0;
        }
    }

    path[0] = '\0';
    return 1;
}

//Finds the filename from the full path. The filename is copied into variable filename.
//Returns 0 if the filename is is extracted successfully.
//Returns 1 if fullpath is just a path (the filename variable is set to empty).
int FindFilename(char* fullpath, char* filename)
{
    int i;
    size_t len = strlen(fullpath);

    if (len == 0 || fullpath[len - 1] == '/')
    {
        filename[0] = '\0';
        return 1;
    }

    for (i = (int)(len - 2); i >= 0; i--)
    {
        if (fullpath[i] == '/')
            break;
    }

    i++;
    sprintf(filename, "%s", &(fullpath[i]));
    return 0;
}
