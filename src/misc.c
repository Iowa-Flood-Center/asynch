#include "misc.h"

//Calculates the order of each link in sys
//Assumes order and complete have N spaces already reserved.
//This could be fairly easily incorporated into riversys.c to be more efficient, if needed.
void CalcHortonOrder(Link** sys,unsigned int N,unsigned int* order,unsigned short int* complete)
{
	unsigned int i,j,loc,parentsval;
	Link *current,*root = NULL;

	Link** stack = malloc(N * sizeof(Link*));
	int stack_size = 0;
	Link** leaves = malloc(N * sizeof(Link*));
	unsigned int leaves_size = 0;
	unsigned short int numparents;

	for(i=0;i<N;i++)
	{
		order[i] = 0;
		complete[i] = 0;
	}

	//Find the root
	for(i=0;i<N;i++)
	{
		if(sys[i]->c == NULL)
		{
			root = sys[i];
			break;
		}
	}

	//Find the leaves
	//Note: this assumes only one root. If there are multiple roots, this needs a for loop.
	stack[0] = root;
	stack_size = 1;
	while(stack_size > 0)
	{
		current = stack[stack_size-1];	//Top of stack
		numparents = current->numparents;

		if(numparents == 0)
		{
			stack_size--;
			leaves[leaves_size] = current;
			leaves_size++;
		}
		else
		{
			//If current is not a leaf, replace it with it's parents
			for(i=0;i<numparents;i++)
			{
				stack[stack_size - 1 + i] = current->parents[numparents - 1 - i];
				stack[stack_size - 1 + i]->c = current;
			}
			stack_size += numparents - 1;
		}
	}

	//Calculate order
	for(i=0;i<leaves_size;i++)
	{
		current = leaves[i];
		order[current->location] = 1;
		//complete[current->location] = 1;	//Exclude order 1 links

		while(current->c != NULL)
		{
			current = current->c;
			loc = current->location;
			numparents = current->numparents;
			for(j=0;j<numparents;j++)
				order[loc] = (order[loc] > order[current->parents[j]->location]) ? order[loc] : order[current->parents[j]->location];

			//Check if current is a complete ordered link
			parentsval = 0;
			for(j=0;j<numparents;j++)
				if(order[loc] == order[current->parents[j]->location])	parentsval++;
			if(parentsval == numparents)
			{
				order[loc]++;
				complete[loc] = 1;
			}
		}
	}

	//for(i=0;i<N;i++)
	//	printf("id = %u  order = %u  complete = %hu\n",sys[i]->ID,order[i],complete[i]);

	free(stack);
	free(leaves);
}

//Generate a .str file for use with complete ordered links
void CreateStrComplete(Link** sys,unsigned int N)
{
	unsigned int i,j,changes,rain_order;
	FILE* output = fopen("Cedar30Complete9.str","w");
	if(!output)
	{
		printf("Error creating file in CreateStrComplete.\n");
		return;
	}
	fprintf(output,"%u\n\n",N);

	//Get the order of each link
	printf("Calculating orders...\n");
	unsigned int* order = malloc(N*sizeof(unsigned int));
	unsigned short int* complete = malloc(N*sizeof(unsigned short int));
	CalcHortonOrder(sys,N,order,complete);

	//Create .str file
	changes = 100;
	rain_order = 9;
	printf("Creating .str file...\n");

	for(i=0;i<N;i++)
	{
		if(complete[i] == 1 && order[i] == rain_order)
		{
			fprintf(output,"%u\n%u\n",sys[i]->ID,changes);
			for(j=0;j<changes-1;j++)
				fprintf(output,"%f %f\n",20.0*j,5.0*j);
			for(j=0;j<changes-1;j++)
				fprintf(output,"%f %f\n",20.0*(changes+j),5.0*j);
			for(j=0;j<changes-1;j++)
				fprintf(output,"%f %f\n",20.0*(2*changes+j),5.0*j);
			fprintf(output,"%f %f\n",20.0*3*changes,0.0);
		}
		else
		{
			fprintf(output,"%u\n%u\n",sys[i]->ID,1);
			fprintf(output,"%f %f\n",0.0,0.0);
		}
		fprintf(output,"\n");
	}

	free(order);
	free(complete);
	fclose(output);
}

//Generates a graph for use by the standalone programs from METIS
void CreateGraph(Link** sys,unsigned int N)
{
	unsigned int i,j;
	int offset = sys[0]->ID + 2;
	FILE* output = fopen("Cedar30.gra","w");
	if(!output)
	{
		printf("Error opening file for graph.\n");
		return;
	}

	fprintf(output,"%u %u\n",N,N-1);

	for(i=0;i<N;i++)
	{
		if(sys[i]->c != NULL)	fprintf(output,"%u ",sys[i]->c->location + offset);
		for(j=0;j<sys[i]->numparents;j++)
			fprintf(output,"%u ",sys[i]->parents[j]->location + offset);
		fprintf(output,"\n");
	}

	fclose(output);
}


//Generates a graph for use by the standalone programs from METIS but takes rainfall into account
/*
void CreateGraphRain(Link** sys,unsigned int N,UnivVars* GlobalVars)
{
	unsigned int i,j,holder;
	char filename[256];
	int offset = -sys[0]->ID + 2;
	unsigned int* total = calloc(N,sizeof(unsigned int));
	float* last = malloc(N*sizeof(float));
	float rainfall_buffer;
	FILE* rainfile;
	FILE* output = fopen("Cedar30Complete1.gra","w");
	if(!output)
	{
		printf("Error opening file for graph.\n");
		return;
	}

	printf("Calculating the number of rainfall changes...\n");

	if(GlobalVars->rain_flag == 2)
	{
		//Calculate the number of rainfall changes that occur
		unsigned int totalfiles = GlobalVars->last_file - GlobalVars->first_file + 1;
		for(i=0;i<totalfiles;i++)
		{
			sprintf(filename,"%srain%i",GlobalVars->rain_filename,GlobalVars->first_file+i);
			rainfile = fopen(filename,"r");
			if(!rainfile)
			{
				printf("Error opening file %s.\n",filename);
				return;
			}

			for(j=0;j<N;j++)
			{
				//Read in the storm data for this link
				fread(&rainfall_buffer,sizeof(float),1,rainfile);

				//This assumes different endianness
				holder = *(unsigned int*) &rainfall_buffer;
				holder = (((holder & 0x0000ffff)<<16) | ((holder & 0xffff0000)>>16));
				holder = (((holder & 0x00ff00ff)<<8) | ((holder & 0xff00ff00)>>8));
				rainfall_buffer = *(float*) &holder;

				if(i == 0 || !(.999 * rainfall_buffer <= last[j] && last[j] <= 1.001 * rainfall_buffer) )
				{
					last[j] = rainfall_buffer;
					total[j]++;
				}
			}

			fclose(rainfile);
		}
	}
	else if(GlobalVars->rain_flag == 1)
	{
		unsigned int id,numtimes;
		rainfile = fopen(GlobalVars->rain_filename,"r");
		if(!rainfile)
		{
			printf("Error opening file %s.\n",GlobalVars->rain_filename);
			return;
		}

		fscanf(rainfile,"%*u");	//N

		for(i=0;i<N;i++)
		{
			fscanf(rainfile,"%u\n%u",&id,&numtimes);
			total[i] = numtimes;
			for(j=0;j<numtimes;j++)
				fscanf(rainfile,"%*f %*f");
		}
	}
	else
	{
		printf("Invalid rain flag (%hu) in CreateGraphRain.\n",GlobalVars->rain_flag);
		return;
	}

	printf("Writing graph to file...\n");

	//Make the graph
	fprintf(output,"%u %u %u%u%u\n",N,N-1,0,1,0);
	for(i=0;i<N;i++)
	{
		fprintf(output,"%u ",total[i]);
		if(sys[i]->c != NULL)	fprintf(output,"%u ",sys[i]->c->location + offset);
		for(j=0;j<sys[i]->numparents;j++)
			fprintf(output,"%u ",sys[i]->parents[j]->location + offset);
		fprintf(output,"\n");
	}

	fclose(output);
	free(total);
	free(last);
}
*/

/*
int* Partition_System_File(Link** sys,unsigned int N,Link** leaves,unsigned int numleaves,unsigned int** my_sys,unsigned int* my_N,unsigned int* my_max_nodes,TransData* my_data,short int *getting)
{
	unsigned int i,j,start_index,end_index,extras,partition,loc;
	unsigned int nodes_per_proc = numleaves/np;	//Number of leaves assigned to each process (except the last)
	char filename[256];

	start_index = nodes_per_proc * my_rank;
	if(my_rank == np-1)		end_index = numleaves;
	else				end_index = nodes_per_proc * (my_rank + 1);
//	*my_N = end_index - start_index;
	*my_N = 0;
	*my_max_nodes = N - numleaves + nodes_per_proc;
	*my_sys = (unsigned int*) malloc(*my_max_nodes * sizeof(unsigned int));	//The indices of this processes links (in sys)
	for(i=0;i<*my_max_nodes;i++)	(*my_sys)[i] = -1;
	for(i=start_index;i<end_index;i++)	(*my_sys)[i-start_index] = leaves[i]->location;
	for(i=0;i<N;i++)	getting[i] = 0;

	//Initialize assignments
	int* assignments = (int*) malloc(N*sizeof(int));
	for(i=0;i<N;i++)	assignments[i] = -1;

	//Open a file
	sprintf(filename,"testgraph.gra.part.%i",np);
	FILE* graph = fopen(filename,"r");
	if(!graph)	printf("Error reading graph file %s.\n",filename);

	//Read in the partitioning data
	//Note: The links should be listed in the graph file in the same order as in the .rvr file
	for(i=0;i<N;i++)
	{
		fscanf(graph,"%i",&(assignments[i]));
		if(assignments[i] == my_rank)
		{
			(*my_sys)[*my_N] = i;
			(*my_N)++;
		}
	}

	//Close file
	fclose(graph);

	//Set the getting array and determine number of sending and receiving links
	for(i=0;i<*my_N;i++)
	{
		//Receiving
		for(j=0;j<sys[(*my_sys)[i]]->numparents;j++)
		{
			loc = sys[(*my_sys)[i]]->parents[j]->location;
			if(assignments[loc] != my_rank)
			{
				getting[loc] = 1;
				my_data->receive_size[assignments[loc]]++;
			}
		}

		//Sending
		if(sys[(*my_sys)[i]]->c != NULL)
		{
			loc = sys[(*my_sys)[i]]->c->location;
			if(assignments[loc] != my_rank)
				my_data->send_size[assignments[loc]]++;
		}
	}

	//Reorder my_sys so that the links with lower numbering are towards the beginning
	merge_sort_distance(sys,*my_sys,*my_N);

	//Allocate space in my_data for recieving and sending
	for(j=0;j<np;j++)
	{
		my_data->receive_data[j] = (Link**) malloc(my_data->receive_size[j] * sizeof(Link*));
		my_data->send_data[j] = (Link**) malloc(my_data->send_size[j] * sizeof(Link*));
	}

	//Set the receive_data and send_data arrays
	int* current_receive_size = (int*) calloc(np,sizeof(int));
	int* current_send_size = (int*) calloc(np,sizeof(int));
	for(i=0;i<*my_N;i++)
	{
		//Receiving
		for(j=0;j<sys[(*my_sys)[i]]->numparents;j++)
		{
			loc = sys[(*my_sys)[i]]->parents[j]->location;
			if(assignments[loc] != my_rank)
			{
				my_data->receive_data[assignments[loc]][current_receive_size[assignments[loc]]] = sys[loc];
				current_receive_size[assignments[loc]]++;
			}
		}

		//Sending
		if(sys[(*my_sys)[i]]->c != NULL)
		{
			loc = sys[(*my_sys)[i]]->c->location;
			if(assignments[loc] != my_rank)
			{
				my_data->send_data[assignments[loc]][current_send_size[assignments[loc]]] = sys[(*my_sys)[i]];
				current_send_size[assignments[loc]]++;
			}
		}
	}

	free(current_receive_size);
	free(current_send_size);

	return assignments;
}
*/

//Calculates the width function for sys. This assumes only one outlet.
void CalculateWidth(Link** sys,unsigned int N)
{
	Link** queue = (Link**) malloc(N*sizeof(Link*));
	unsigned int* width = (unsigned int*) calloc(N,sizeof(unsigned int));
	unsigned int* width_I = (unsigned int*) calloc(N,sizeof(unsigned int));
	unsigned int* width_E = (unsigned int*) calloc(N,sizeof(unsigned int));
	unsigned int* distance = (unsigned int*) calloc(N,sizeof(unsigned int));
	unsigned int i,j,curr_idx;
	Link* current;
	FILE* output;

	//Find the root
	for(i=0;i<N;i++)
		if(sys[i]->c == NULL)	break;

	//Find the distance of each link to the root
	queue[0] = sys[i];
	distance[0] = 1;
	width[0] = 1;
	width_I[0] = 1;
	curr_idx = 1;

	for(i=0;i<N;i++)
	{
		current = queue[i];
		for(j=0;j<current->numparents;j++)
		{
			queue[curr_idx] = current->parents[j];
			distance[curr_idx] = distance[i] + 1;
			width[distance[curr_idx]-1]++;
			if(current->parents[j]->numparents == 0)	width_E[distance[curr_idx]-1]++;
			else						width_I[distance[curr_idx]-1]++;
			curr_idx++;
		}
	}

	//Find the largest size
	for(j=0;j<N;j++)
		if(width[j] == 0)	break;

	//Write results to file
	output = fopen("width","w");
	if(output == NULL)	printf("Error opening file for width function.\n");
	for(i=0;i<j;i++)
		fprintf(output,"%u\n",width[i]);
	fclose(output);

	output = fopen("widthE","w");
	if(output == NULL)	printf("Error opening file for widthE function.\n");
	for(i=0;i<j;i++)
		fprintf(output,"%u\n",width_E[i]);
	fclose(output);

	output = fopen("widthI","w");
	if(output == NULL)	printf("Error opening file for width function.\n");
	for(i=0;i<j;i++)
		fprintf(output,"%u\n",width_I[i]);
	fclose(output);

	//Cleanup
	free(width);
	free(width_I);
	free(width_E);
	free(queue);
	free(distance);

	//Compute averages
	double ave_l = 0.0;
	double ave_l_E = 0.0;
	double ave_l_I = 0.0;
	double ave_a = 0.0;
	unsigned int count_E = 0;
	unsigned int count_I = 0;
	for(i=0;i<N;i++)
	{
		ave_l += sys[i]->params->ve[1];
		ave_a += sys[i]->params->ve[2];
		if(sys[i]->numparents == 0)
		{
			ave_l_E += sys[i]->params->ve[1];
			count_E++;
		}
		else
		{
			ave_l_I += sys[i]->params->ve[1];
			count_I++;
		}
	}
	ave_l /= N;
	ave_a /= N;
	ave_l_E /= count_E;
	ave_l_I /= count_I;

	printf("\nAverage hillslope area: %.12f\n",ave_a);
	printf("Average length: %.12f\n",ave_l);
	printf("Average length external: %.12f\n",ave_l_E);
	printf("Average length internal: %.12f\n\n",ave_l_I);

	printf("Width function computed and stored to disk.\n");
	getchar();
}

//Calculates the evaporation rates for model 30
void EvapRates(double t,VEC* y_i,VEC* params,VEC* global_params,VEC* ans)
{
	//States
	double q = y_i->ve[0];
	double s_p = y_i->ve[1];
	double h_w = y_i->ve[2];
	double theta = y_i->ve[3];

	//Global params
	//double v_0 = global_params->ve[0];
	double lambda_1 = global_params->ve[1];
	//double lambda_2 = global_params->ve[2];
	//double A_r = global_params->ve[4];
	double K_T = global_params->ve[5];
	double e_pot = global_params->ve[7];

	//Local params
	//double L_h = params->ve[0];
	double A_h = params->ve[1];
	//double A_up = params->ve[2];
	double H_b = params->ve[3];
	double H_h = params->ve[4];
	double maxinf = params->ve[5];
	double K_SAT = params->ve[6];
	//double S_H = params->ve[7];
	//double n_vh = params->ve[8];
	double b = params->ve[9];
	double c = params->ve[10];
	double d = params->ve[11];
	double K_Q = params->ve[12];
	double V_T = params->ve[13];
	double c_1 = params->ve[14];
	double c_2 = params->ve[15];
	double c_3 = params->ve[16];
	double c_4 = params->ve[17];
	double c_5 = params->ve[18];
	double c_6 = params->ve[19];
	double c_7 = params->ve[20];


	//Other
	//double H_T = H_b + H_h;
	double h_rel = h_w - H_b * 1e-3;
	//double H_relmax = H_T - H_b;
	double H_relmax = (H_h > 1e-12) ? H_h : 1e-6;
	double hoverH = h_rel / H_relmax;
	double a_IoverA_H = b * hoverH + c * sq(hoverH) + d * sq(hoverH) * hoverH;
	if(a_IoverA_H < 0.0)	a_IoverA_H = 1e-6;
	if(a_IoverA_H > 0.9999)	a_IoverA_H = 0.9999;
	double a_PoverA_H = 1.0 - a_IoverA_H;
	double v_H = c_1 * pow(s_p * 1e-3,2.0/3.0);
	if(v_H > 350.0)	v_H = 350.0;
	double v_ssat = V_T * a_IoverA_H;
	double v_sunsat = V_T - v_ssat;
	double d_soil = (1.0 - theta) * H_b;
	double RC = (s_p <= 0.0) ? 0.0 : s_p*(s_p + 2*d_soil) / sq(s_p + d_soil);
	double deriv_a_I = c_7 * (b + 2*c*hoverH + 3*d*sq(hoverH));
	if(deriv_a_I <= (1e-0)) deriv_a_I = 1e-0;

	//Evaporation
	double D_unsat = 1e-3 * v_sunsat * theta / A_h;
	double D_sat = 1e-3 * v_ssat / A_h;

	double C_p,C_unsat,C_sat,C_T;
	if(e_pot > 0.0)
	{
		C_p = s_p / e_pot;
		C_unsat = D_unsat / e_pot;
		C_sat = D_sat / e_pot;
		C_T = C_p + C_unsat + C_sat;
	}
	else
	{
		C_p = 0.0;
		C_unsat = 0.0;
		C_sat = 0.0;
		C_T = 0.0;
	}

	double Corr_evap = (C_T > 1.0) ? 1.0/C_T : 1.0;
	double e_p = Corr_evap * C_p * e_pot;
	double e_sat = Corr_evap * C_sat * e_pot;
	double e_unsat = Corr_evap * C_unsat * e_pot;

	ans->ve[0] = e_p + e_sat + e_unsat;
	ans->ve[1] = e_p;
	ans->ve[2] = e_sat;
	ans->ve[3] = e_unsat;
}


//Calculates the evaporation rates for model 30
void EvapRates_linear(double t,VEC* y_i,VEC* params,VEC* global_params,VEC* ans)
{
	double q = y_i->ve[0];		//[m^3/s]
	double s_p = y_i->ve[1];	//[m]
	double s_a = y_i->ve[2];	//[m]

	double e_pot = global_params->ve[6];

	//Evaporation
	double C_p,C_a,C_T;
	if(e_pot > 0.0)
	{
		C_p = s_p / e_pot;
		C_a = s_a / e_pot;
		C_T = C_p + C_a;
	}
	else
	{
		C_p = 0.0;
		C_a = 0.0;
		C_T = 0.0;
	}

	double Corr_evap = (C_T > 1.0) ? 1.0/C_T : 1.0;
	double e_p = Corr_evap * C_p * e_pot;
	double e_a = Corr_evap * C_a * e_pot;

	ans->ve[0] = e_p + e_a;
	ans->ve[1] = e_p;
	ans->ve[2] = e_a;
}

