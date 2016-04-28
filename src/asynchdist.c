#include <stdio.h>
#include <stdbool.h>

#if !defined(_MSC_VER)
#include <unistd.h>
#endif

#include "mpi.h"
#include "asynch_interface.h"
#include "optparse.h"

#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

// Global variables
int my_rank = 0;
int np = 0;

int Output_Linkid(double t,VEC y_i,VEC global_params,VEC params,int state,void* user);
void Set_Output_User_LinkID(asynchsolver* asynch);

//Print to stdout only for process of rank 0
int print_out(const char* format, ...)
{
    int res = 0;
    if (my_rank == 0)
    {
        va_list args;
        va_start(args, format);
        res = vprintf(format, args);
        va_end(args);
    }

    return res;
}

//Print to stderr only for process of rank 0
int print_err(const char* format, ...)
{
    int res = 0;
    if (my_rank == 0)
    {
        va_list args;
        va_start(args, format);
        res = vfprintf(stderr, format, args);
        va_end(args);
    }

    return res;
}

//Make sure we finalize MPI
void asynch_onexit(void)
{
    int flag;
    MPI_Finalized(&flag);
    if (!flag)
        MPI_Finalize();
}

int main(int argc,char* argv[])
{
    int res;

    //Initialize MPI stuff
    res = MPI_Init(&argc, &argv);
    if (res == MPI_SUCCESS)
        atexit(asynch_onexit);
    else
    {
        print_err("Failed to initialize MPI");
        exit(EXIT_FAILURE);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    //Command line options
    bool debug = false;
    bool help = false;
    bool version = false;    

    //Parse command line
    struct optparse options;
    optparse_init(&options, argv);
    struct optparse_long longopts[] = {
        { "debug", 'd', OPTPARSE_NONE },
        { "help", 'h', OPTPARSE_NONE },
        { "version", 'v', OPTPARSE_NONE },
        { 0 }
    };
    int option;
    while ((option = optparse_long(&options, longopts, NULL)) != -1) {
        switch (option) {
        case 'd':
            debug = true;
            break;
        case 'h':
            help = true;
            break;
        case 'v':
            version = true;
            break;
        case '?':
            print_err("%s: %s\n", argv[0], options.errmsg);
            exit(EXIT_FAILURE);
        }
    }

    if (version) print_out("This is %s\n", PACKAGE_STRING);
    if (help)
    {
        print_out("Usage: asynch <global file>\n", PACKAGE_STRING);
        print_out(
"  -d [--debug]   : Wait for the user input at the begining of the program (useful" \
"                   for attaching a debugger)\n" \
"  -v [--version] : Print the current version of ASYNCH\n");
        exit(EXIT_SUCCESS);
    }
    if (version || help) exit(EXIT_SUCCESS);

	//Parse remaining arguments
    char *global_filename = optparse_arg(&options);
    if (global_filename == NULL && my_rank == 0)
    {
        print_err("Command line parameter required:  A universal variable file (.gbl).\n");
        exit(EXIT_FAILURE);
    }

    if (debug)
    {
        //Disable stdout buffering
        setvbuf(stdout, NULL, _IONBF, 0);

        //When the program first starts to execute, at the very beginning of our program, we 
        //ask the user to type some sort of input to simply stall the application until start your
        //"Attach to Process" and you can attach to all the different threads in your program.
        if (my_rank == 0)
        {
            printf("You may now attach the debugger then press enter.\n");
            //fflush(stdout);
            getchar();
        }

        MPI_Barrier(MPI_COMM_WORLD); // All threads will wait here until you give thread 0 an input
    }

	//Declare variables
    double start, stop;
	double total_time;
	asynchsolver* asynch;

    print_out("\nBeginning initialization...\n*****************************\n");
	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();

	//Init asynch object and the river network
	asynch = Asynch_Init(MPI_COMM_WORLD,&argc,&argv);
	print_out("Reading global file...\n");
	Asynch_Parse_GBL(asynch,argv[1]);
	print_out("Loading network...\n");
	Asynch_Load_Network(asynch);
	print_out("Partitioning network...\n");
	Asynch_Partition_Network(asynch);
	print_out("Loading parameters...\n");
	Asynch_Load_Network_Parameters(asynch,0);
	print_out("Reading dam and reservoir data...\n");
	Asynch_Load_Dams(asynch);
	print_out("Setting up numerical error data...\n");
	Asynch_Load_Numerical_Error_Data(asynch);
	print_out("Initializing model...\n");
	Asynch_Initialize_Model(asynch);
	print_out("Loading initial conditions...\n");
	Asynch_Load_Initial_Conditions(asynch);
	print_out("Loading forcings...\n");
	Asynch_Load_Forcings(asynch);
	print_out("Loading output data information...\n");
	Asynch_Load_Save_Lists(asynch);
	print_out("Finalizing network...\n");
	Asynch_Finalize_Network(asynch);
	print_out("Calculating initial step sizes...\n");
	Asynch_Calculate_Step_Sizes(asynch);

	if(my_rank == 0)
	{
		printf("\nModel type is %u.\nGlobal parameters are:\n",asynch->GlobalVars->type);
		Print_Vector(asynch->GlobalVars->global_params);
		printf("\n");
	}

	//Setup output for link id, if needed
	int id_setup = Asynch_Check_Output(asynch,"LinkID");
	if(id_setup != -1)
	{
		Set_Output_User_LinkID(asynch);
		Asynch_Set_Output(asynch,"LinkID",ASYNCH_INT,&Output_Linkid,NULL,0);
	}

	//Prepare output files
	Asynch_Prepare_Temp_Files(asynch);
	Asynch_Write_Current_Step(asynch);		//!!!! Wow, this sucks. Is there a way to get rid of it? !!!!
	Asynch_Prepare_Peakflow_Output(asynch);
	Asynch_Prepare_Output(asynch);

	//Make sure everyone is good before getting down to it...
	printf("Process %i (%i total) is good to go with %i links.\n",my_rank,np,asynch->my_N);
	ASYNCH_SLEEP(1);
	MPI_Barrier(MPI_COMM_WORLD);

	if(my_rank == 0)
	{
		stop = MPI_Wtime();
		total_time = stop - start;
		printf("Finished initialization. Total time for initialization: %f\n\n\nComputing solution at each link...\n************************************\n", stop - start);
	}
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);

	//Perform the calculations
    start = MPI_Wtime();
	Asynch_Advance(asynch,1);
	MPI_Barrier(MPI_COMM_WORLD);
    stop = MPI_Wtime();

	//Out information
	total_time += stop - start;
	print_out("\nComputations complete. Total time for calculations: %f\n", stop - start);
	if(asynch->sys[asynch->my_sys[0]]->c == NULL)
	{
		printf("[%i]: The solution at ID %i at time %.12f is\n",my_rank,asynch->sys[asynch->my_sys[0]]->ID,asynch->sys[asynch->my_sys[0]]->last_t);
		Print_Vector(asynch->sys[asynch->my_sys[0]]->list->tail->y_approx);
	}

	//Take a snapshot
	Asynch_Take_System_Snapshot(asynch,NULL);

	//Create output files
	Asynch_Create_Output(asynch,NULL);
	Asynch_Create_Peakflows_Output(asynch);

	//Clean up
	Asynch_Delete_Temporary_Files(asynch);
	Asynch_Free(asynch);

	return 0;
}



int Output_Linkid(double t,VEC y_i,VEC global_params,VEC params,int state,void* user)
{
	return ((Link*)user)->ID;
}

//!!!! Gross, but not sure how else to handle this. Maybe with a lot of interface functions? !!!!
void Set_Output_User_LinkID(asynchsolver* asynch)
{
	unsigned int i,my_N = asynch->my_N,*my_sys = asynch->my_sys;
	Link** sys = asynch->sys;

	for(i=0;i<my_N;i++)
		sys[my_sys[i]]->output_user = (void*) sys[my_sys[i]];
}




