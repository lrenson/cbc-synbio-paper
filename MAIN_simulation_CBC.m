clear all;
close all;
clc;
strname = 'INSERT_A_NAME'
%% *** COMMON PARAMETERS ***
MPC_flag = false; %true= MPC strategy; false = Proportional control
Parameters.MPC_flag = MPC_flag;

ODE_flag = true; %true= ODE system; false = SDE system
Parameters.ODE_flag = ODE_flag;
Parameters.sim.endflag = false; %if false, all steady states are evaluated after a fixed time.
                                %if true, the reference steps after a variable time, depending
                                %on whether the system is at steady state or not.
Parameters.sim.nrep = 1; %1,5 or 10; number of repetitions, in case of a stochastic simulation
Parameters.sim.npts = 30; %or 30, 10, 12; number of points to be caught by the CBC algorithm
% Parameters.sim.reference =linspace(100,1500,Parameters.sim.npts);%65-1200,1500 Higher span for P controller
Parameters.sim.reference =fliplr(linspace(0,1500,Parameters.sim.npts ));%0-1200,0-1500
% Initial_x0 = [37; 0.8; 2237; 61 ];
initIPTG = 1; %initial IPTG
initATC = 25; %initial aTc
Initial_x0 = [ 2.55166047363230 38.7108543679906 102.155003051775 1196.05604522200]';
[~,xtemp] = ode45(@(t,x)LugagneToggle(x,[25 initIPTG]),[0 5*12*6],Initial_x0);
ya_init= xtemp(end,4);
Initial_x0 = xtemp(end,:)';

Parameters.time.Ts=5;%sampling time
Parameters.time.h=60;%minutes x h

%tolerances: linear fit
%     Parameters.sim.tol_u = 0.001;%0.1/0.0001
%     Parameters.sim.tol_e = 0.1;%0.1
Parameters.sim.tol_u = 0.01;%/0.01;/0.0001
Parameters.sim.tol_e = 0.1;%0.1

Parameters.sim.solvertime= 0.01; %step for the SDEsolver
Parameters.sim.Par=LugagneParameters(); %Parameters of the Toggle Switch model
Parameters.sim.sigma = 100; % strenght percentage of the noise signal

if Parameters.sim.endflag
%  variable time
Parameters.time.Tend = 5;%10 h
Parameters.sim.timetoavg=12;% 12(xTs minutes from the end of each simulation)
Parameters.sim.n_wait = 12*3 ;%-->80

else
Parameters.time.Tend = 10;%12 h
Parameters.sim.timetoavg=12;% 12?(xTs minutes from the end of each simulation)
Parameters.sim.n_wait = 12*3 ;%-->80
end

%% reference points and repetitions

npts = Parameters.sim.npts;
nrep = Parameters.sim.nrep;
ref_ampl= Parameters.sim.reference; 

%% ______________________________________________  %%

if MPC_flag==1
    %% -----같같� MPC 같같�-----
    % % Identified model and kalman
    load('state_space_model_from_DATA_HIGH.mat');%identified system (continuous time)
    load('kalman_systems.mat');%system with kalman filter
    ID_SYS.kalm_sys=Kalman_SYS_HIGH;
    ID_SYS.system_td=c2d(SYS_TC_HIGHER,Parameters.time.Ts);% control horizon
    x_hat2 = ID_SYS.system_td.x0;
    for k=1:4
    x_hat2 = ID_SYS.kalm_sys.a*x_hat2+ID_SYS.kalm_sys.b*[Initial_x0(3);Initial_x0(4);initATC;initIPTG];
    end

    %% Sim Parameters
    
    Parameters.sim.x_hat = x_hat2;
    Parameters.sim.N = 3;%prediction horizon for the MPC
    Parameters.sim.ctrlbandwidth = 0.3; %control bandwidth tolerance in the MPC optimization routine 0.3
    Parameters.ctrl.last_ctrl = 0; %Last control action, initial value
    
    %% CELLS STRUCT DATA
    RESULTS.YA = {};
    RESULTS.BP = {};
    RESULTS.J = {};%cost function vector
    RESULTS.RMS = {};%vector to save the invasiveness
    RESULTS.U = {};%Control action vector
    RESULTS.e_TetR = {};%error
    RESULTS.derr_TetR= {};%error derivative
    RESULTS.du_IPTG = {};%control derivative
    RESULTS.x_esti_LacI = {};% estimated state
    RESULTS.x_esti_TetR = {};% estimated state
    RESULTS.x_esti_mrnaLacI = {};% estimated state
    RESULTS.x_esti_mrnaTetR = {};% estimated state
    RESULTS.x_real_LacI = {};%real state
    RESULTS.x_real_TetR = {};%real state
    RESULTS.Endtime = {};%Total time of each point of the algorithm (in samples)
    
    stack.error = {};
    stack.derr = {};
    stack.du = {};
else
    %% ---P---
    
    Parameters.ctrl.Kp = 0.0008*2;%*5 %0.005,0.001; min: 0.0005/ working det 0.003 10 pts/ 0.0008 working stoch
    Parameters.ctrl.Kd = 0;
    Parameters.ctrl.Ki = 0;
    
    Parameters.ctrl.past_err = 0 ;
    Parameters.ctrl.I = 0 ;
    Parameters.ctrl.min_u = 0 ;
    Parameters.ctrl.max_u = 1 ;
    
    %% CELLS STRUCT DATA
    RESULTS.YA = {};
    RESULTS.BP = {};
    RESULTS.RMS = {};%vector to save the invasiveness
    RESULTS.U = {};%Control action vector
    RESULTS.e_TetR = {};%error
    RESULTS.derr_TetR= {};%error derivative
    RESULTS.du_IPTG = {};%control derivative
    RESULTS.x_real_LacI = {};%real state
    RESULTS.x_real_TetR = {};%real state
    RESULTS.Endtime = {};%Total time of each point of the algorithm (in samples)
    
end
%% Main loop

for i= 1:nrep
    
    %Initial parameters in case of multiple simulations
    Parameters.sim.x0 = Initial_x0;
    Parameters.Init_input.INPUT_act=initIPTG;%IPTG
    Parameters.Init_input.fixed_act=initATC;%Atc
    Parameters.Init_input.IPTG0=initIPTG;
    Parameters.Init_input.Atc0=initATC;
    
    
    if MPC_flag==1
        Parameters.sim.x_hat = x_hat2;%ID_SYS.system_td.x0;
        XESTI = zeros(4,1);
        J_vect = [];% Cost function vector
        x_esti_LacI = [];% estimated state
        x_esti_TetR = [];% estimated state
        x_esti_mrnaLacI = [];% estimated state
        x_esti_mrnaTetR = [];% estimated state
        
    else 
        ID_SYS = 0;
   
    end
    
    XOUT = zeros(4,1);
    
    yA = [ya_init];%System Output
    bp = [Parameters.Init_input.INPUT_act];%bifurcation parameter(IPTG)
    rmsvect = [];%vector to save the invasiveness
    U = []; %control
    e_TetR = [];%error
    derr_TetR = [];%error derivative
    du_IPTG = [];%control derivative
    x_real_LacI = [];%real state
    x_real_TetR = [];%real state
    LAST_U = []; %last control action sample
    Endtime = []; %time associated to each reference step

     
    for j = 1:npts
        
        Parameters.sim.iter = j;
        if MPC_flag==1
            
            if j==1
                Parameters.ctrl.last_c = -1;
            else
                Parameters.ctrl.last_c = CONTROL.mean_u;
            end
        end
        %INNER CBC LOOP: called for each new reference-steady state point
        [ STATES, CONTROL, Parameters, endtime, Optimizer,Internal_concentration] = CBC_inner_loop( ID_SYS, Parameters , ref_ampl(j));
            
        
        if MPC_flag==1
            
            stack.error{i,j} = Optimizer.error;
            stack.derr{i,j} = Optimizer.derr ;
            stack.du{i,j} = Optimizer.du;
            J_vect = [J_vect Optimizer.J];
            x_esti_mrnaLacI = [x_esti_mrnaLacI STATES.x_esti(1,1:endtime)];
            x_esti_mrnaTetR = [x_esti_mrnaTetR STATES.x_esti(2,1:endtime)];
            x_esti_LacI = [x_esti_LacI STATES.x_esti(3,1:endtime)];
            x_esti_TetR = [x_esti_TetR STATES.x_esti(4,1:endtime)];
            XESTI(:,1) = STATES.x_esti(:,endtime);

        end
        yA =  [yA STATES.mean_y_out(2)] ;%update the output with last sample
        bp = [bp CONTROL.mean_u];%update bif parameter with the new sample
        rmsvect = [rmsvect CONTROL.inv];

        e_TetR = [e_TetR STATES.err_stack(1,1:endtime)];
        derr_TetR = [derr_TetR STATES.err_stack(2,1:endtime)];
        du_IPTG = [du_IPTG STATES.err_stack(3,1:endtime)];
        U = [U CONTROL.u(1:endtime)];
        x_real_LacI = [x_real_LacI STATES.x_real(3,1:endtime)];
        x_real_TetR = [x_real_TetR STATES.x_real(4,1:endtime)];
        Endtime = [Endtime endtime];
        
        Parameters.Init_input.Atc0=Internal_concentration(1);
        Parameters.Init_input.IPTG0=Internal_concentration(2);
        
        LAST_U = [LAST_U CONTROL.u(endtime)];
        XOUT(:,1) = STATES.x_real(:,endtime);

    end
    
    if MPC_flag==1
        RESULTS.J{i} = J_vect;
        RESULTS.x_esti_LacI{i} = x_esti_LacI;% estimated state
        RESULTS.x_esti_TetR{i} = x_esti_TetR;% estimated state
        RESULTS.x_esti_mrnaLacI{i} = x_esti_mrnaLacI;% estimated state
        RESULTS.x_esti_mrnaTetR{i} = x_esti_mrnaTetR;% estimated state
        
    end
    RESULTS.YA{i} = yA;
    RESULTS.BP{i} = bp;
    RESULTS.RMS{i} = rmsvect;%vector to save the invasiveness
    RESULTS.U{i} = U;
    RESULTS.e_TetR{i} = e_TetR;%error
    RESULTS.derr_TetR{i} = derr_TetR;%error derivative
    RESULTS.du_IPTG{i} = du_IPTG;%control derivative
    
    RESULTS.x_real_LacI{i} = x_real_LacI;%real state
    RESULTS.x_real_TetR{i} = x_real_TetR;%real state
    RESULTS.Endtime{i} = Endtime;

end
save(strcat(strname,'.mat'))

%% COMPUTE SLOPE of error and control
% clear SLOPE Err_temp_overNpts U_temp_overNpts Err_temp U_temp slope_e slope_u 
% SLOPE.e = {};
% SLOPE.u = {};

% endtime
timetoavg = Parameters.sim.timetoavg;
for k = 1:1%Parameters.sim.nrep
    temp_slope_e = [];
    temp_slope_u = [];
    Err_temp_overNpts = [ RESULTS.e_TetR{k} ];
    U_temp_overNpts = [ RESULTS.U{k} ];
    Err_temp = [];
    U_temp = [];
    %variable endtime    
    endtime = [ RESULTS.Endtime{k} ];
    for i = 1:Parameters.sim.npts
        %variable endtime
        if i==1
            Err_temp = [Err_temp_overNpts(1:sum(endtime(1:i)))];
            U_temp = [U_temp_overNpts(1:sum(endtime(1:i)))];
        else
            Err_temp = [Err_temp_overNpts(sum(endtime(1:i-1))+1:sum(endtime(1:i)))];
            U_temp = [U_temp_overNpts(sum(endtime(1:i-1))+1:sum(endtime(1:i)))];
        end

        [slope_e,slope_u] = Slope(Err_temp,U_temp,timetoavg);
        temp_slope_e = [temp_slope_e slope_e];
        temp_slope_u = [temp_slope_u slope_u];
        SLOPE(i,:) = [slope_e slope_u];
        
%         Err_temp = [Err_temp_overNpts(((i-1)*endtime)+1:(i*endtime))];
%         U_temp = [U_temp_overNpts(((i-1)*endtime)+1:(i*endtime))];
%         [slope_e,slope_u] = Slope(Err_temp,U_temp,timetoavg);
%         temp_slope_e = [temp_slope_e slope_e];
%         temp_slope_u = [temp_slope_u slope_u];
    end  
%     SLOPE.e{k} = temp_slope_e;
%     SLOPE.u{k} = temp_slope_u;
end
YA=RESULTS.YA;
BP=RESULTS.BP;

% save(strcat('Data ',strname,'.mat'),'YA','BP','SLOPE')