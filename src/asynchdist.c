#include <stdio.h>

#if !defined(_MSC_VER)
#include <unistd.h>
#endif

#include "mpi.h"
#include "asynch_interface.h"

int my_rank;
int np;

int Output_Linkid(double t,VEC* y_i,VEC* global_params,VEC* params,int state,void* user);
void Set_Output_User_LinkID(asynchsolver* asynch);

int main(int argc,char* argv[])
{
	//Initialize MPI stuff
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&np);

	//Parse input
	if(argc < 2)
	{
		if(my_rank == 0)
		{
			printf("Command line parameter required:  A universal variable file (.gbl).\n");
			printf("\n");
		}
		MPI_Finalize();
		return 1;
	}

#if !defined(NDEBUG)
    //Disable stdout buffering
    setvbuf(stdout, NULL, _IONBF, 0);

    //When the program first starts to execute, at the very beginning of our program, we 
    //ask the user to type some sort of inpu to simply stall the application until start your
    //"Attach to Process" and you can attach to all the different threads in your program.
    if (my_rank == 0)
    {
        printf("You may now attach the debugger then press enter.\n");
        //fflush(stdout);
        getchar();
    }

    MPI_Barrier(MPI_COMM_WORLD); // All threads will wait here until you give thread 0 an input
#endif

	//Declare variables
	time_t start,stop;
	double total_time;
	asynchsolver* asynch;

	if(my_rank == 0)
		printf("\nBeginning initialization...\n*****************************\n");
	MPI_Barrier(MPI_COMM_WORLD);
	start = time(NULL);

	//Init asynch object and the river network
	asynch = Asynch_Init(MPI_COMM_WORLD,&argc,&argv);
	if(my_rank == 0)	printf("Reading global file...\n");
	Asynch_Parse_GBL(asynch,argv[1]);
	if(my_rank == 0)	printf("Loading network...\n");
	Asynch_Load_Network(asynch);
	if(my_rank == 0)	printf("Partitioning network...\n");
	Asynch_Partition_Network(asynch);
	if(my_rank == 0)	printf("Loading parameters...\n");
	Asynch_Load_Network_Parameters(asynch,0);
	if(my_rank == 0)	printf("Reading dam and reservoir data...\n");
	Asynch_Load_Dams(asynch);
	if(my_rank == 0)	printf("Setting up numerical error data...\n");
	Asynch_Load_Numerical_Error_Data(asynch);
	if(my_rank == 0)	printf("Initializing model...\n");
	Asynch_Initialize_Model(asynch);
	if(my_rank == 0)	printf("Loading initial conditions...\n");
	Asynch_Load_Initial_Conditions(asynch);
	if(my_rank == 0)	printf("Loading forcings...\n");
	Asynch_Load_Forcings(asynch);
	if(my_rank == 0)	printf("Loading output data information...\n");
	Asynch_Load_Save_Lists(asynch);
	if(my_rank == 0)	printf("Finalizing network...\n");
	Asynch_Finalize_Network(asynch);
	if(my_rank == 0)	printf("Calculating initial step sizes...\n");
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
		Asynch_Set_Output(asynch,"LinkID",ASYNCH_INT,(void (*)(double,VEC*,VEC*,VEC*,int,void*)) &Output_Linkid,NULL,0);
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
		stop = time(NULL);
		total_time = difftime(stop,start);
		printf("Finished initialization. Total time: %f\n\n\nComputing solution at each link...\n************************************\n",total_time);
	}
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);

	//Perform the calculations
	time(&start);
	Asynch_Advance(asynch,1);
	MPI_Barrier(MPI_COMM_WORLD);
	time(&stop);

	//Out information
	total_time += difftime(stop,start);
	if(my_rank == 0)	printf("\nComputations complete. Total time for calculations: %f\n",difftime(stop,start));
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
	MPI_Finalize();
	return 0;
}



int Output_Linkid(double t,VEC* y_i,VEC* global_params,VEC* params,int state,void* user)
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




