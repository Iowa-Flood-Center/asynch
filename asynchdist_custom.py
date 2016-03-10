from asynch_py.asynch_interface import *
import sys
import timeit


#User defined routines ******************************************************************************
@ASYNCH_SETPARAMSIZES_DATATYPE
def SetParamSizes_MyModel(GlobalVars_ptr,lib_ptr):
	lib = cast(lib_ptr,py_object).value
	GlobalVars = GlobalVars_ptr.contents

	GlobalVars.uses_dam = 0
	GlobalVars.params_size = 8
	GlobalVars.dam_params_size = 0
	GlobalVars.area_idx = 0
	GlobalVars.areah_idx = 2
	GlobalVars.disk_params = 3
	GlobalVars.convertarea_flag = 0
	GlobalVars.num_forcings = 3
	GlobalVars.min_error_tolerances = 6

#	GlobalVars.dim = GlobalVars.problem_dim = 6
#	GlobalVars.template_flag = 0
#	GlobalVars.assim_flag = 0
#	GlobalVars.diff_start = 0
#	GlobalVars.no_ini_start = 3
#	GlobalVars.uses_dam = 0
#	GlobalVars.params_size = 8
#	GlobalVars.iparams_size = 0
#	GlobalVars.dam_params_size = 0
#	GlobalVars.area_idx = 0
#	GlobalVars.areah_idx = 2
#	GlobalVars.disk_params = 3
#	GlobalVars.num_dense = 2
#	GlobalVars.convertarea_flag = 0
#	GlobalVars.num_forcings = 3

#	GlobalVars.dense_indices = lib.Allocate_CUINT_Array(GlobalVars.num_dense)
#	GlobalVars.dense_indices[0] = 0	#Discharge
#	GlobalVars.dense_indices[1] = 5	#Subsurface


@ASYNCH_CONVERT_DATATYPE
def ConvertParams_MyModel(params,model_type,lib_ptr):
	params.contents.ve[1] *= 1000	#L: km -> m
	params.contents.ve[2] *= 1e6	#A_h: km^2 -> m^2

@ASYNCH_ROUTINES_DATATYPE
def InitRoutines_MyModel(link_ptr,model_type,exp_imp,dam,lib_ptr):
	lib = cast(lib_ptr,py_object).value
	link = link_ptr.contents

	link.dim = 6;
	link.no_ini_start = 3;
	link.diff_start = 0;

	link.num_dense = 2;
	link.dense_indices = lib.Allocate_CUINT_Array(link.num_dense)
	link.dense_indices[0] = 0;
	link.dense_indices[1] = 5;

	if link.res:
		link.f = LinearHillslope_Reservoirs_MyModel
		#link.f = cast(lib.LinearHillslope_Reservoirs_extras,ASYNCH_F_DATATYPE)
		link.RKSolver = cast(lib.ForcedSolutionSolver,ASYNCH_RKSOLVER_DATATYPE)
	else:
		link.f = LinearHillslope_MyModel
		#link.f = cast(lib.LinearHillslope_MonthlyEvap_extras,ASYNCH_F_DATATYPE)
		link.RKSolver = cast(lib.ExplicitRKSolver,ASYNCH_RKSOLVER_DATATYPE)
	link.alg = cast(None,ASYNCH_ALG_DATATYPE)
	link.state_check = cast(None,ASYNCH_STATECHECK_DATATYPE)
	link.CheckConsistency = cast(lib.CheckConsistency_Nonzero_AllStates_q,ASYNCH_CONSISTENCY_DATATYPE)

@ASYNCH_PRECALCULATIONS_DATATYPE
def Precalculations_MyModel(link_i,global_params,params,disk_params,params_size,dam,model_type,lib_ptr):
	#Order of parameters: A_i,L_i,A_h,k2,k3,invtau,c_1,c_2
	#The numbering is:	0   1   2   3  4    5    6   7
	#Order of global_params: v_r,lambda_1,lambda_2,RC,v_h,v_g (,v_B)
	#The numbering is:        0      1        2     3  4   5     6
	vals = params.contents.ve
	global_vals = global_params.contents.ve
	A_i = vals[0]
	L_i = vals[1]
	A_h = vals[2]
	v_r = global_vals[0]
	lambda_1 = global_vals[1]
	lambda_2 = global_vals[2]
	RC = global_vals[3]
	v_h = global_vals[4]
	v_g = global_vals[5]

	vals[3] = v_h * L_i / A_h * 60.0	#[1/min]  k2
	vals[4] = v_g * L_i / A_h * 60.0	#[1/min]  k3
	vals[5] = 60.0*v_r*A_i**lambda_2 / ((1.0-lambda_1)*L_i)	#[1/min]  invtau
	vals[6] = RC*(0.001/60.0)		#(mm/hr->m/min)  c_1
	vals[7] = (1.0-RC)*(0.001/60.0)		#(mm/hr->m/min)  c_2


@ASYNCH_INITIALIZEEQS_DATATYPE
def ReadInitData_MyModel(global_params,params,qvs,dam,y_0,model_type,diff_start,no_ini_start,user_ptr,lib_ptr):
	init_conds = y_0.contents.ve

	#For this type, the extra states need to be set (3,4,5)
	init_conds[3] = 0.0
	init_conds[4] = 0.0
	init_conds[5] = init_conds[0]

	return 0

@ASYNCH_F_DATATYPE
def LinearHillslope_MyModel(t,y_i_ptr,y_p,numparents,global_params_ptr,forcing_values,qvs_ptr,params_ptr,state,user_ptr,ans_ptr):

	#Unpack the pointers
	y_i = y_i_ptr.contents.ve
	params = params_ptr.contents.ve
	global_params = global_params_ptr.contents.ve
	ans = ans_ptr.contents.ve

	#Unpack the parameters
	lambda_1 = global_params[1]
	v_B = global_params[6]
	L = params[1]
	A_h = params[2]
	k2 = params[3]
	k3 = params[4]
	invtau = params[5]
	c_1 = params[6]
	c_2 = params[7]

	#Unpack the system states
	q = y_i[0]	#[m^3/s]
	s_p = y_i[1]	#[m]
	s_a = y_i[2]	#[m]
	q_b = y_i[5]	#[m^3/s]

	#Calculate fluxes
	q_pl = k2 * s_p
	q_al = k3 * s_a

	#Evaporation
	e_pot = forcing_values[1] * (1e-3/(30.0*24.0*60.0))	#[mm/month] -> [m/min]

	if e_pot > 0.0:
		C_p = s_p / e_pot
		C_a = s_a / e_pot
		C_T = C_p + C_a
	else:
		C_p = 0.0
		C_a = 0.0
		C_T = 0.0

	if C_T > 1.0:	Corr_evap = 1.0/C_T
	else:		Corr_evap = 1.0

	e_p = Corr_evap * C_p * e_pot
	e_a = Corr_evap * C_a * e_pot

	#Discharge
	ans[0] = -q + (q_pl + q_al) * A_h/60.0
	for i in range(0,numparents):
		ans[0] += y_p[i].contents.ve[0]
	ans[0] = invtau * q**lambda_1 * ans[0]

	#Hillslope
	ans[1] = forcing_values[0]*c_1 - q_pl - e_p
	ans[2] = forcing_values[0]*c_2 - q_al - e_a

	#Additional states
	ans[3] = forcing_values[0]*c_1
	ans[4] = q_pl
	ans[5] = q_al * A_h - q_b * 60.0
	for i in range(0,numparents):
		ans[5] += y_p[i].contents.ve[5] * 60.0
	ans[5] *= v_B/L


@ASYNCH_F_DATATYPE
def LinearHillslope_Reservoirs_MyModel(t,y_i_ptr,y_p_ptr,numparents,global_params_ptr,forcing_values_ptr,qvs_ptr,params_ptr,state,user_ptr,ans_ptr):
	ans = ans_ptr.contents.ve
	#forcing_values = forcing_values_ptr.contents
	ans[0] = forcing_values_ptr[2]
	ans[1] = 0.0
	ans[2] = 0.0
	ans[3] = 0.0
	ans[4] = 0.0
	ans[5] = 0.0


#Routines for custom output **************************************************************
def Set_Output_User_LinkID(asynch):
	my_N = asynch.Get_Local_Number_Links()
	for i in range(my_N):
		linkid = asynch.Get_Local_LinkID(i)
		asynch.Set_Size_Local_OutputUser_Data(i,sys.getsizeof(linkid))
		asynch.Copy_Local_OutputUser_Data(i,linkid)

@ASYNCH_OUTPUT_INT_DATATYPE
def Output_Linkid(t,y_i,global_params,params,iparams,state,user):
	return cast(user,py_object).value


#Routines for custom peakflow output ****************************************************
@ASYNCH_PEAKOUTPUT_DATATYPE
def NewClassicPeakflow(ID,peak_time,peak_value,params,global_params,conversion,area_idx,user,outputbuffer):
	array = (c_char*1028).from_address(addressof(outputbuffer.contents))
	to_write = str(2*ID)+' ' + str(conversion*params.contents.ve[area_idx])+' ' + str(peak_time)+' ' + str(peak_value.contents.ve[0])+'\n\0'
	array[:len(to_write)] = to_write


#Main program *******************************************************************

#Parse command line arguments
numargs = len(sys.argv)
if numargs != 2:
	print 'Need an input .gbl file'
	sys.exit(1)

#Initialize the solver object
asynch = asynchsolver()
comm = asynch.comm
my_rank = asynch.my_rank
np = asynch.np

#Load a custom model
asynch.Custom_Model(SetParamSizes_MyModel,ConvertParams_MyModel,InitRoutines_MyModel,Precalculations_MyModel,ReadInitData_MyModel)

#Build the network
if my_rank == 0:	print 'Parsing .gbl...'
asynch.Parse_GBL(sys.argv[1])
if my_rank == 0:	print "Loading network..."
asynch.Load_Network()
if my_rank == 0:	print "Partitioning network..."
asynch.Partition_Network()
if my_rank == 0:	print "Loading parameters..."
asynch.Load_Network_Parameters(False)
if my_rank == 0:	print "Reading dam and reservoir data..."
asynch.Load_Dams()
if my_rank == 0:	print "Setting up numerical error data..."
asynch.Load_Numerical_Error_Data()
if my_rank == 0:	print "Initializing model..."
asynch.Initialize_Model()
if my_rank == 0:	print "Loading initial conditions..."
asynch.Load_Initial_Conditions()
if my_rank == 0:	print "Loading forcings..."
asynch.Load_Forcings()
if my_rank == 0:	print "Loading output data information..."
asynch.Load_Save_Lists()
if my_rank == 0:	print "Finalizing network..."
asynch.Finalize_Network()
if my_rank == 0:	print "Calculating initial step sizes..."
asynch.Calculate_Step_Sizes()


#Check for LinkID in output
printing_linkid = asynch.Check_Output('LinkID')
if printing_linkid == 0:
	Set_Output_User_LinkID(asynch)
	asynch.Set_Output('LinkID',ASYNCH_INT,Output_Linkid,None)

#Check for NewClassic peakflows
printing_newpeakflow = asynch.Check_Peakflow_Output('NewClassic')
if printing_newpeakflow == 0:
	asynch.Set_Peakflow_Output('NewClassic',NewClassicPeakflow)

#Prepare outputs
if my_rank == 0:	print 'Preparing output files...'
asynch.Prepare_Temp_Files()
asynch.Write_Current_Step()
asynch.Prepare_Peakflow_Output()
asynch.Prepare_Output()

#Ready to go
N = asynch.Get_Number_Links()
my_N = asynch.Get_Local_Number_Links()
print '[',my_rank,'/',np,']:','''I'm ready to go with''',my_N,'/',N,'links.'

#Advance solver
comm.Barrier()
start = timeit.default_timer()
asynch.Advance(True)
comm.Barrier()
stop = timeit.default_timer()
if my_rank == 0:	print 'Calculations done! Total time',stop-start,'seconds'

#Take a snapshot
asynch.Take_System_Snapshot(None)

#Create output files
asynch.Create_Output(None)
asynch.Create_Peakflows_Output()

#Cleanup
if printing_linkid == 0:
	asynch.Free_OutputUser_Data()




