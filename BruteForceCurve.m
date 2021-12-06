clc;clear all;close all;
%Parameter initialization
    sim.solvertime = 0.05;
    Tend = 3*24*60;% 3days
    timetoavg = 60/sim.solvertime;%samples
    tspan = [0 Tend];
    nSamples = 30*10;% Must be even
    nRepetition = 1;
    ext_Atc = 25;
    sim.Par = LugagneParameters();
    
    %% ODE CASE:


    % x0 = [37    0.8    2237    61]; %First Initial Condition
    x0= [0 0 0 0]';
    range=[0 1];
    ctrl_vectd=linspace(range(1),range(2),nSamples/2);
    output_temp = zeros(nRepetition,1);
    state_temp = zeros(nRepetition,4);
    ATC_IPTG_temp = zeros(nRepetition,2);
    ATC_IPTG_vect_f = zeros(nSamples/2,2);
    ATC_IPTG_vect_b = zeros(nSamples/2,2);
    output_vect_f = zeros(nSamples/2,2);
    output_vect_b = zeros(nSamples/2,2);
    output_all_b = zeros(nSamples/2,nRepetition);
    output_all_f = zeros(nSamples/2,nRepetition);
    
    [t,state_out]=ode45(@(t,x)LugagneToggle(x,[ext_Atc 0.1]),[0 Tend],x0);
    x0 = state_out(end,:)';
    
    
    %% Simulation loop (forward)
    for samp=1:nSamples/2
        ctrl_val=ctrl_vectd(samp);
        for rep=1:nRepetition%questa parte era per il codice stocastico
            [t,xout]=ode45(@(t,x)LugagneToggle(x,[ext_Atc ctrl_val]),[0 Tend],x0);
            xout=xout';
            output_temp(rep)=xout(4,end);
            state_temp(rep,:)=xout(:,end);
%             ATC_IPTG_temp(rep,:)=[Atco IPTGo];
        end
        output_all_f(samp,:)=output_temp;
        output_vect_f(samp,:)=[min(output_temp) max(output_temp)];
        ATC_IPTG_vect_f(samp+1,:) = mean(ATC_IPTG_temp,1);
        x0=mean(state_temp,1);
        x0=x0';
        fprintf('Sample number %d \n',samp);
    end
    
    
    [t,state_out]=ode45(@(t,x)LugagneToggle(x,[ext_Atc 1]),[0 Tend],x0);
    x0=state_out(end,:)';
    ctrl_vectu=fliplr(ctrl_vectd);%+.05;
    %% Simulation loop (backward)
    for samp=1:nSamples/2
        ctrl_val=ctrl_vectu(samp);
        for rep=1:nRepetition
            [t,xout]=ode45(@(t,x)LugagneToggle(x,[ext_Atc ctrl_val]),[0 Tend],x0);
            xout=xout';
            output_temp(rep)=xout(4,end);
            state_temp(rep,:)=xout(:,end);
%             ATC_IPTG_temp(rep,:)=[Atco IPTGo];
        end
        output_all_b(samp,:)=output_temp;
        output_vect_b(samp,:)=[min(output_temp) max(output_temp)];
%         ATC_IPTG_vect_b(samp+1,:) = mean(ATC_IPTG_temp,1);
        x0=mean(state_temp,1);
        x0=x0';
        fprintf('Sample number %d \n',samp);
    end
      
    control_vect=[ctrl_vectd ctrl_vectu];
    output_vect=[output_vect_f;output_vect_b];
    %% Plotting results
    figure;
    for samp=1:nSamples/2
        hold on;
        if output_vect(samp,1)==output_vect(samp,2)
            plot(control_vect(samp),output_vect(samp,1),'b.');
        else
            line([control_vect(samp) control_vect(samp)],output_vect(samp,:),'color','b','LineWidth',2);
        end
        
    end
    for samp=nSamples/2:nSamples
        hold on;
        if output_vect(samp,1)==output_vect(samp,2)
            plot(control_vect(samp),output_vect(samp,1),'r.');
        else
            line([control_vect(samp) control_vect(samp)],output_vect(samp,:),'color','r','LineWidth',2);
        end
    end
%     pause;
%     %% save variables
%     var_name=strcat('output_TetR_min_max_values_and_control_IPTG',num2str(sim.sigma),'%noise.mat');
%     save(var_name,'output_vect','control_vect');
%     var_name = strcat('all_points_',num2str(sim.sigma),'%noise.mat');
%     save(var_name,'output_all_f','output_all_b') ;
    %% save figures
    fig = gcf; % current figure
    ax = fig.CurrentAxes;% current axes
    ax.FontSize = 16;
    ax.TickDir = 'out';
    ax.TickLength = [0.02 0.02];
    ax.YLim = [0 1500];
    ax.XLim = [0 1.2];
    ax.XLim = [0 1200];
    ax.Title.String = 'Stochastic Bruteforce diagram';
    xlabel('IPTG');ylabel('TetR');
%     figename=strcat('Stochastic_Bruteforce_diagram_noiseAmplitude',num2str(sim.sigma),'%');
%     savefig(fig,figename);
    

    
    
%% LOOP SDE CASE
for iter=1:1

    if iter==1
        sim.sigma=100;
    elseif iter==2
        sim.sigma=70;
    elseif iter==3
        sim.sigma=50;
    end
   
    % x0 = [37    0.8    2237    61]; %First Initial Condition
    x0= [0 0 0 0]';
    range=[0 1];
    ctrl_vectd=linspace(range(1),range(2),nSamples/2);
    output_temp = zeros(nRepetition,1);
    state_temp = zeros(nRepetition,4);
    ATC_IPTG_temp = zeros(nRepetition,2);
    ATC_IPTG_vect_f = zeros(nSamples/2,2);
    ATC_IPTG_vect_b = zeros(nSamples/2,2);
    output_vect_f = zeros(nSamples/2,2);
    output_vect_b = zeros(nSamples/2,2);
    output_all_b = zeros(nSamples/2,nRepetition);
    output_all_f = zeros(nSamples/2,nRepetition);
    
    [t,state_out]=ode45(@(t,x)LugagneToggle(x,[ext_Atc 0.1]),[0 Tend],x0);
    x0 = state_out(end,:)';
    
    
    %% Simulation loop (forward)
    for samp=1:nSamples/2
        ctrl_val=ctrl_vectd(samp);
        for rep=1:nRepetition%questa parte era per il codice stocastico
            [tout, xout, Atco, IPTGo] = SDESolver(sim, x0, tspan, ATC_IPTG_vect_f(samp,1), ATC_IPTG_vect_f(samp,2), ctrl_val, ext_Atc);
            output_temp(rep)= mean(xout(4,end-timetoavg:end));%mean
            state_temp(rep,:)=xout(:,end);
            ATC_IPTG_temp(rep,:)=[Atco IPTGo];
        end
        output_all_f(samp,:)=output_temp;
        output_vect_f(samp,:)=[min(output_temp) max(output_temp)];
        ATC_IPTG_vect_f(samp,:) = mean(ATC_IPTG_temp,1);%samp+1
        x0=mean(state_temp,1);
        x0=x0';
        fprintf('Sample number %d \n',samp);
    end
    
    
    [t,state_out]=ode45(@(t,x)LugagneToggle(x,[ext_Atc 1]),[0 Tend],x0);
    x0=state_out(end,:)';
    ctrl_vectu=fliplr(ctrl_vectd);%+.05;
    %% Simulation loop (backward)
    for samp=1:nSamples/2
        ctrl_val=ctrl_vectu(samp);
        for rep=1:nRepetition
            [tout, xout, Atco, IPTGo] = SDESolver(sim, x0, tspan, ATC_IPTG_vect_b(samp,1), ATC_IPTG_vect_b(samp,2), ctrl_val, ext_Atc);
            output_temp(rep)= mean(xout(4,end-timetoavg:end));;
            state_temp(rep,:)= xout(:,end);
            ATC_IPTG_temp(rep,:)=[Atco IPTGo];
        end
        output_all_b(samp,:)=output_temp;
        output_vect_b(samp,:)=[min(output_temp) max(output_temp)];
        ATC_IPTG_vect_b(samp,:) = mean(ATC_IPTG_temp,1);
        x0=mean(state_temp,1);
        x0=x0';
        fprintf('Sample number %d \n',samp);
    end
      
    control_vect=[ctrl_vectd ctrl_vectu];
    output_vect=[output_vect_f;output_vect_b];
    %% Plotting results
    figure;
    for samp=1:nSamples/2
        hold on;
        if output_vect(samp,1)==output_vect(samp,2)
            plot(control_vect(samp),output_vect(samp,1),'b.');
        else
            line([control_vect(samp) control_vect(samp)],output_vect(samp,:),'color','b','LineWidth',2);
        end
        
    end
    for samp=nSamples/2:nSamples
        hold on;
        if output_vect(samp,1)==output_vect(samp,2)
            plot(control_vect(samp),output_vect(samp,1),'r.');
        else
            line([control_vect(samp) control_vect(samp)],output_vect(samp,:),'color','r','LineWidth',2);
        end
    end
%     pause;
%     %% save variables
%     var_name=strcat('output_TetR_min_max_values_and_control_IPTG',num2str(sim.sigma),'%noise.mat');
%     save(var_name,'output_vect','control_vect');
%     var_name = strcat('all_points_',num2str(sim.sigma),'%noise.mat');
%     save(var_name,'output_all_f','output_all_b') ;
    %% save figures
    fig = gcf; % current figure
    ax = fig.CurrentAxes;% current axes
    ax.FontSize = 16;
    ax.TickDir = 'out';
    ax.TickLength = [0.02 0.02];
    ax.YLim = [0 1500];
    ax.XLim = [0 1.2];
%     ax.XLim = [0 1200];
%     ax.Title.String = 'Stochastic Bruteforce diagram';
    xlabel('IPTG');ylabel('TetR');
   
end
