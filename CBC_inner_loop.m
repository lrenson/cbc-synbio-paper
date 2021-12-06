
function [ STATES, CONTROL, Par, endtime, Optimizer, Internal_concentration] = CBC_inner_loop( ID_SYS, Par , ref_ampl)

Ts=Par.time.Ts;%sampling time
h=Par.time.h;%minutes x h
Tend=Par.time.Tend;%final time (h)

timetoavg = Par.sim.timetoavg;%number of samples we use to average the control action
wait_time = Par.sim.n_wait;%time waited before starting to compute the SLOPE
INPUT_act = Par.Init_input.INPUT_act;%IPTG
fixed_act = Par.Init_input.fixed_act;%Atc
%initial concentrations
IPTG0=Par.Init_input.IPTG0;
Atc0=Par.Init_input.Atc0;

if Par.MPC_flag==1
    Cost_val = [];
    
    Optimizer.error = [];
    Optimizer.derr = [];
    Optimizer.du = [];
    
    x_esti = zeros(4,Tend*60/Ts+1);% estimated state
    
    %% Identified Discrete time model TS for the MPC
    system_td = ID_SYS.system_td;
    
    A=system_td.A;
    B=system_td.B;
    C=system_td.C;
    D=system_td.D;
    
    sys.a = A; sys.b = B; sys.c = C; sys.d = D;
    x_hat = Par.sim.x_hat;
    
    %% Identified Kalman Predictor
    kalm_sys = ID_SYS.kalm_sys;
end
% vectors 
e = zeros(1,Tend*h/Ts+1);%error
derr = zeros(1,Tend*h/Ts+1);%derivative of the error
du = zeros(1,Tend*h/Ts+1);%derivative of the control action
x_real = zeros(4,Tend*h/Ts+1);%real state
y_LacI = zeros(1,Tend*h/Ts+1);
y_TetR = zeros(1,Tend*h/Ts+1);

x0 = Par.sim.x0;%initial state
ampl=ref_ampl;

time_dim = ((Tend*h)/Ts);
%% =========================================================================

disp('EXPERIMENT INITIALIZATION');
reference=ampl*ones(1,time_dim);
disp('Initialization Concluded. Reference for the Control Loop evaluated.');
%========================================================================
j=1;%index for the elements of the reference signal
% Continuous time system used for the simulation of the plant
Par.sim.j=j;

y_LacI(j) = x0(3);
y_TetR(j)= x0(4);

%error first iteration
e(j) =(reference(j)-y_TetR(j));

if Par.MPC_flag==1
x_esti(:,j) = x_hat;
end
x_real(:,j) = x0;
endflag=1;

%% STARTING POINT OF THE LOOP
while (endflag)
    clc;
    disp(['Control iteration number: ',num2str(j)]);
    disp('STATE: EVALUATION OF CONTROL ACTION u(k).');
    disp(['error offset: ',num2str(e(j))])
    disp(['reference value: ',num2str(reference(j))])
    disp(['Control input: ',num2str(INPUT_act(end))])
    
    %% STATE ESTIMATION
    
    if Par.MPC_flag==1

        %% MPC INPUT Computation
        [input_act, cost_val, stack] = MPC_optimization_block(x_esti(:,j), reference, j, INPUT_act(end), sys, Par.sim, Par.ctrl, Ts);
        Optimizer.error = [Optimizer.error stack.err];
        Optimizer.du = [Optimizer.du stack.d_u];
        Cost_val = [Cost_val; cost_val];
        
        x_hat=kalm_sys.a*x_hat+kalm_sys.b*[y_LacI(j);y_TetR(j);fixed_act;INPUT_act(end)];%INPUT_act(j)
        x_esti(:,j+1)=x_hat;
        
    else
        %% PID INPUT Computation
        [input_act, Par.ctrl] = PIDCtrl(reference(j),y_TetR(j), Par.ctrl, Par.time.Ts);
        
    end
     
    fprintf('Control input=%.3f \n',input_act);
    disp(['Control input: ',num2str(input_act)])
    INPUT_act = [INPUT_act input_act];% we need to save it as many times as TC
    disp('STATE: EVALUATION OF THE PLANT STATE AND ERROR');
    %% PLANT SIMULATION
    
    if (Par.ODE_flag==1)%deterministic
        [~,xout]=ode45(@(t,x)LugagneToggle(x,[fixed_act INPUT_act(end)]),[0 Ts],x_real(:,j));
        Atco = 0; IPTGo = 0;
        xout = xout';
    else%stochastic
        [~, xout, Atco, IPTGo]  = SDESolver(Par.sim, x_real(:,j),[0 Ts], Atc0(end), IPTG0(end), INPUT_act(end), fixed_act);
    end
    x_real(:,j+1)=xout(:,end);% x_esti(:,j+1) = xout(:,end);
    
    y_LacI(j+1)=xout(3,end); %LacI
    y_TetR(j+1)=xout(4,end); %TetR
    Atc0=[Atc0 Atco];
    IPTG0=[IPTG0 IPTGo];
    
    %% derivative CTRL ESTIMATION
    if j==1
        du(j)=(INPUT_act(end)-INPUT_act(1))/Ts;
    else
        du(j)=((INPUT_act(end)-INPUT_act(end-1))/Ts);
    end
    %% derivative ERR ESTIMATION
    if j==1
        derr(j)=0;%save error from prev iter
    else
        derr(j)=((e(j)-e(j-1))/Ts);
    end
    %% error 
    if j==time_dim
        e(j+1)=(reference(j)-y_TetR(j+1));
    else
        e(j+1)=(reference(j+1)-y_TetR(j+1));
    end
    %% steady state check
    if (j >= wait_time)
        [endflag] = steadystate_check(INPUT_act, e, j, endflag, Par.sim);
    end
    endtime = j;
    j=j+1;
    Par.sim.j=j;
    if (j==time_dim)
        endflag=0;
    end

end
%% OUTPUT VARIABLES

STATES.x_real = x_real;
err_stack =[e;derr;du];
STATES.err_stack = err_stack;
STATES.mean_y_out= [ mean(y_LacI(endtime-timetoavg:endtime)); mean(y_TetR(endtime-timetoavg:endtime))];
CONTROL.mean_u= mean(INPUT_act(endtime-timetoavg:endtime));
CONTROL.inv=rms(INPUT_act(endtime-timetoavg:endtime)-CONTROL.mean_u);%/mean_u
CONTROL.u = INPUT_act(1:endtime);
Optimizer.J = 0;
Internal_concentration=[Atc0(endtime) IPTG0(endtime)];
Par.Init_input.INPUT_act = INPUT_act(endtime);
Par.sim.x0 = x_real(:,endtime);

if Par.MPC_flag==1
    STATES.x_esti = x_esti;
    Optimizer.J = mean (Cost_val);
    Par.sim.x_hat = x_esti(:,endtime);
end

end


%% STEADY STATE CHECK

function [endflag] = steadystate_check(INPUT_act,e,j, endflag, sim)
timetoavg = sim.timetoavg;

%% CASE SLOPE computation
tempErr = e( j-timetoavg:j );
tempctrl = INPUT_act( j-timetoavg:j );
x_u = (0:timetoavg)';
y_u = tempctrl';
X_u = [ones(length(x_u),1) x_u];
B_u = X_u\y_u;
slope_u = abs(B_u(2));
x_e = (0:timetoavg)';
y_e = tempErr';
X_e = [ones(length(x_e),1) x_e];
B_e = X_e\y_e;
slope_e = abs(B_e(2));
if (slope_e<=sim.tol_e)&&(slope_u<=sim.tol_u)
    endflag=not(sim.endflag);
end

end

%% PID CONTROLLER
function [u, ctrl] = PIDCtrl(ref, y, ctrl, Ts)
%NO ANTI WIND-UP ACTION!!!!

% Compute errors
err = ref-y ; %error
derr = (err-ctrl.past_err)/Ts ;

% Compute control signal
P=ctrl.Kp*err;%take the gains from the external parameters
I=ctrl.I+ctrl.Ki*(err*Ts);
D=ctrl.Kd*derr;
u = P+I+D ;
% saturation
if u <= ctrl.min_u
    u = ctrl.min_u ;
end
if u >= ctrl.max_u
    u = ctrl.max_u ;
end
% Store errors
ctrl.past_err = err ;
ctrl.I = I ;
end

