function [u_x_esti] = MPC_BSIM(output, ref, IPTG0, aTc0, x_hat_init, u_previous_pt)
LacI=output(1);
TetR=output(2);

% disp(output)
% disp(ref);
% disp(IPTG0);
% disp(aTc0)
% disp(x_hat_init);
% disp(u_previous_pt);

if size(x_hat_init,2)>size(x_hat_init,1)
    x_hat_init=x_hat_init';
end

Par.time.Ts=5;%sampling time
Par.time.h=60;%minutes x h
Par.time.Tend=10;
Ts=Par.time.Ts;%sampling time

IPTG = IPTG0;%init IPTG
aTc = aTc0;%init Atc


%% Identified Discrete time model TS for the MPC
load('state_space_model_from_DATA_HIGH.mat');%identified system (continuous time)
load('kalman_systems.mat');%system with kalman filter
ID_SYS.kalm_sys=Kalman_SYS_HIGH;
ID_SYS.system_td=c2d(SYS_TC_HIGHER,Par.time.Ts);% control horizon

system_td = ID_SYS.system_td;

A=system_td.A;
B=system_td.B;
C=system_td.C;
D=system_td.D;

sys.a = A; sys.b = B; sys.c = C; sys.d = D;

% Identified Kalman Predictor
kalm_sys = ID_SYS.kalm_sys;


%%

y_LacI = LacI;
y_TetR = TetR;
reference=ref;

%% STATE ESTIMATION


x_hat=kalm_sys.a*x_hat_init+kalm_sys.b*[y_LacI;y_TetR;aTc;IPTG];
x_esti=x_hat;

%% INPUT Evaluation through optimization
[input_act, ~ ] = state1_MPC_discrete(x_hat, reference, sys, u_previous_pt);
%Cost_val = [Cost_val; cost_val];

% fprintf('Control input=%.3f \n',input_act);
disp(['Control input: ',num2str(input_act)])

u=input_act;

u_x_esti=[u;x_esti];

end


function [u_opt, cost_val]= state1_MPC_discrete(x_est, ref, system_td, u_previous_pt)

%% Sim parameters
ctrlbandwidth = 0.3;
N = 3;
u_previous_pt = u_previous_pt;
%% lower and upper bounds 
% % CONSTRAINTS ON U
if u_previous_pt<0
     lb = 0*ones(1,N);
     ub = 1*ones(1,N);
else
    lb = (1-ctrlbandwidth)*u_previous_pt*ones(1,N);
    ub = min([(1+ctrlbandwidth)*u_previous_pt ,1])*ones(1,N);
end

     
 [a, cost_val] = ga(@(a) costfunction(a, N, x_est, ref, system_td),N,[],[],[],[],lb,ub,[]);
 
%  [~, stack] = costfunction(a, N, x_est, ref, system_td);
 
 u_opt=a(1);
 
end

function [cost, stack] = costfunction(a, N, x_e, ref, system_td)
%COST FUNCTION

%%%%% discrete model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=system_td.a;
B=system_td.b;
C=system_td.c;
D=system_td.d;
alpha = 1;

cost = 0;
X = x_e;
err = zeros(1,N);

for k=1:N
      
    u=a(k);%a(k)
    
    XP = A*(X)+B*[25;u];
    Y = C*X;
    
%     YP = C*XP;
    X = XP;
    
    err(k)=(ref-Y(2));
    
    cost = cost+(alpha)*(N-k+1)*((err(k))^2);
    
if nargout>1
    stack.err = err;
end
end
end
