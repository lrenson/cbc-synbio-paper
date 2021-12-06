Here there is a quick guide on how to use the Matlab2019b code.

MAIN_simulation_CBC.m is the runnable code. It contains all the variables 
that the user might want to might want to modify (such as what type of controller to be used,
how many points to be collected, sampling time and so on). All variables are commented. 
This code calls CBC_inner_loop for each reference value (which can as well be defined by the user).


CBC_inner_loop.m contains the  control based continuation routine for a single reference step.
It can either call the MPC optimization algorithm (MPC_optimization_block.m) or the PID algorithm (PIDCtrl.m)
which works as a simple proportional control. 

plots_cbc.m allows the user to reproduce Control based continuation data seen as well in the manuscript.
The user will need to upload a '.mat' file and then run the firts section of the code ('COLORS') 
and then an appropriate section depending on the type of data loaded (deterministic, stochastic, bsim etc.)

ParametersEstimation_fun.m allows the user to try the parameters estimation process. Main data resulting 
from the estimation are contained in 'allDATA_toplot_scurves.mat' and 'allDATA_toplot_scurvesWITHconstraintsonIPTG.mat'.
To plot those, refer to the plot function above.

BruteForceCurve.m is a code to generate bruteforce data.

LugagneToggleODE.m, LugagneToggle.m, LugagnePropAndD.m and LugagneParameters.m are all functions 
related to the toggle switch models (ODEs and SDEs). They cannot run on their own.

Optimizer.m is called by ParametersEstimation_fun.m for the estimation of the toggle switch parameters.

state_space_model_from_DATA_HIGH.mat and identified_systems_and_kalman.mat contain data used by
the MPC for the optimization process (such as the linear model associated to the toggle switch and
the kalman filter to reconstruct the states).

SDESolver.m contains the solver for SDEs.