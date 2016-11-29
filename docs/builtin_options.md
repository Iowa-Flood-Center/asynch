# Built-In Options

Through a global file, many options are selected, including which model to use, what time series to output, and how to calculate peakflows Each of these can be customized by the user However, several precreated options do exist This section provides a list of these options

## Built-In Models

Although numerous models exist, a quick list of the most frequently used models is given in Figure 5 A complete list can be found in the source codes problems.c and defnetype.c These models can be activated by selecting the corresponding model type in the input global file See Section 6 1 1.

Model Type | Description | States
- | - | -
21 | Linear reservoir model, with dams | q, S, Ss, Sg
40 | Linear reservoir model, with qvs dams | q, S, Ss, Sg
190 | IFC model with constant runoff | q, sp, ss
191 | IFC model with constant runoff | q, sp, ss, sprecip, Vr, qb
252 | IFC toplayer model | q, sp, st, ss
253 | IFC toplayer model, with reservoirs | q, sp, st, ss
254 | IFC toplayer model, with reservoirs | q, sp, st, ss, sprecip, Vr, qb

Figure 5: Some models built-in to ASYNCH

## Built-In Output Time Series

Figure 6 contains the names and a description of built-in output time series These outputs are defined in the source file `modeloutputs.c`. Up to seven states can be outputted with the built-in output time series. In addition to these, users can create their own time series outputs See Section 9 3.

Output Name |  Description  
- | -
Time | Simulation time
TimeI | Simulation time, truncated to an integer
State0 | State 0 of the model
State1 | State 1 of the model
State6 | State 6 of the model

Figure 6: Built-in output time series

## Built-In Peakflow Functions

Two built-in peakflow functions exist: *Classic* and *Forecast*. The two are described in Figure 7. The peak discharges are the largest values obtained in the state with index 0 in the state vectors. The time to peak for the Classic function is given in simulation time For Forecast, the time to peak is measured in unix time. The time period output is a parameter that can be altered by user programs to provide additional output information.

Function Name |  Outputs
- | -
Classic | Link ID, upstream area, time to peak, peak discharge
Forecast | Link ID, time to peak, peak discharge, time period

Figure 7: Built-in peakflow functions

## Built-In Runge-Kutta Methods

The ASYNCH solver is based upon using Runge-Kutta methods at the link level. These methods are selected either in the input global file or in a Runge-Kutta data file by the RK index in Figure 8.

RK index | Name | Local order / Dense order
- | - | -
0 | Kutta's Method | 3 / 2
1 | The RK Method | 4 / 3
2 | Dormand and Prince's Method | 5 / 4
3 | RadauII 3A | 3 / 2

Figure 8: Built-in Runge-Kutta methods

The application of these methods is done through the `RKSolver` routine in the `UnivVars` structure. This is set with a call to the `InitRoutines` method See Section 8 1 3. Several choices exist for the `RKSolver`, see Figure 9. Some solvers are only appropriate if the model uses ODEs, while others support DAEs. Similarly, some methods support discontinuity states, while others do not. Currently, only one method is equipped to handle stif ODEs. Certainly, the routine `ExplicitRKIndex1SolverDam` could be used to solve any problem. However, using a more appropriate solver is signifcantly more efficient.

Name | DAEs | Discontinuities | Stif
- | - | - | -
ExplicitRKSolver | No | Yes | Yes
ExplicitRKIndex1SolverDam | No | No | No
ExplicitRKIndex1Solver | Yes | No | Yes
ExplicitRKSolverDiscont | No | No | No
RadauRKSolver | No | No | Yes

Figure 9: Built-in Runge-Kutta solvers
