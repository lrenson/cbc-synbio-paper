function [u_opt, cost_val, stack]= MPC_optimization_block(x_est, ref, j, last_ctrl, system_td, sim, ctrl, Ts)
%state1_MPC_discrete   Discrete Model Predictive Control.
%   [U_OPT, COST_VAL, STACK] = MPC_optimization_block(X0,REF,J,LAST_CTRL,SYSTEM_TD, SIM, CTRL, TS)  
%   computes an array of optimal control actions [u_opt1 u_opt2...], the
%   cost function value COST_VAL and a STACK with the predicted error,
%   error derivative and control derivative.
%   It returns only the first action [u_opt1] of the optimal
%   control sequence. X0 is the array of state, REF is the setpoint
%   value for the control problem. J is the current iteration of the
%   algorithm, LAST_CTRL is the last value of the control input and 
%   SYSTEM_TD contains the matrices A,B,C, D to evaluate the cost function.
%   CTRL contains previous coleccted point control action, to compute constraints on u.
%   SIM contains some simulation parameters and TS is the CBC sampling
%   time.

%% Sim parameters
ctrlbandwidth = sim.ctrlbandwidth;
N = sim.N;
u_previous_pt = ctrl.last_c;
%% lower and upper bounds 
% % CONSTRAINTS ON U
if u_previous_pt<0
     lb = 0*ones(1,N);
     ub = 1*ones(1,N);
else
    lb = (1-ctrlbandwidth)*u_previous_pt*ones(1,N);
    ub = min([(1+ctrlbandwidth)*u_previous_pt ,1])*ones(1,N);
end
% % CASE NO CONSTRAINTS
% lb = 0*ones(1,N);
% ub = 1*ones(1,N);

 if length(ref)<(j+N+1)
     ref = [ref,ref(end)*ones(1,N+1)];
 end
     
 [a, cost_val] = ga(@(a) costfunction(a, N, Ts, x_est, ref((j+1):(j+N+1)), last_ctrl, system_td),N,[],[],[],[],lb,ub,[]);
 [~, stack] = costfunction(a, N, Ts, x_est, ref((j+1):(j+N+1)), last_ctrl, system_td);
 
 u_opt=a(1);
 
end

%%
function [cost, stack] = costfunction(a,N,Ts,x_e,ref,last_ctrl, system_td)
%COST FUNCTION

%%%%% discrete model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=system_td.a;
B=system_td.b;
C=system_td.c;
D=system_td.d;

cost = 0;
X = x_e;
u0=last_ctrl;
err = zeros(1,N);
d_err = zeros(1,N);
d_u = zeros(1,N);

for k=1:N
      
    u=a(k);%a(k)
    
    XP = A*(X)+B*[25;u];
    Y = C*X;
    
    YP = C*XP;
    X = XP;
    
    err(k)=(ref(k)-Y(2));
    d_err(k)=((ref(k+1)-YP(2)-ref(k)+Y(2))/Ts);
    
    d_u(k)=(u-u0)/Ts;
    u0=u;
   
    cost = cost+(N-k+1)*((err(k))^2);
    
if nargout>1
    stack.err = err;
    stack.d_err = d_err;
    stack.d_u = d_u;
end
end
end