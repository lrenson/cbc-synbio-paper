%PLOTS 
%% COLORS:
cyan = rgb('cyan');
blue = rgb('blue');
red = rgb('red');
gray = rgb('DarkGray');
orange = rgb('orange');
green = rgb('darkgreen');
teal = rgb('teal');

openfig('coco_bif_curve')
lin=findobj(gca,'Type', 'line');
x_bif = lin(1).XData;
y_bif=lin(1).YData;
x_coco=lin(2).XData;
y_coco=lin(2).YData;
close

%% DETERMINISTIC PLOT

% load('.\MATLAB DATA\MPC ODE 30 pts fix ref time tend10_IPTGdecreasing.mat')
% load('.\MATLAB DATA\P ODE 30 pts fix ref time tend10_IPTGdecreasing.mat')

clear yA_overNrep bp_overNrep
for k=1:nrep
yA_overNrep(k,:) = RESULTS.YA{k}; %TetR
bp_overNrep(k,:) = RESULTS.BP{k};
end

mean_tetR = mean(yA_overNrep,1);%mean of each column
mean_bp = mean(bp_overNrep,1);

subplot(1,2,1)
time=Parameters.time;
temp_x_real_TetR = [];
temp_u = [];
t = 0:time.Ts:length(RESULTS.x_real_TetR{1,1})*time.Ts-time.Ts;
for i=1:nrep
temp_x_real_TetR(i,:) = RESULTS.x_real_TetR{1,i};
temp_u(i,:) = RESULTS.U{1,i};
end
% figure()
plot(x_coco, y_coco, '--', 'LineWidth', 2, 'Color', gray);
hold on;
plot(x_bif, y_bif, '*', 'LineWidth', 6, 'Color', red);

x_real_TetR = temp_x_real_TetR;
u = temp_u;

plot(u,x_real_TetR,':','linewidth',2, 'Color', red);
plot(mean_bp, mean_tetR, '.', 'MarkerSize', 20, 'Color', teal);


xlabel('IPTG [mM]')
ylabel('TetR [a.u.]')
set(gca,'Fontsize',12,'FontWeight','bold');


time=Parameters.time;
% figure(2)
subplot(2,2,2)

for i=1:nrep
t = 0:time.Ts:length(RESULTS.x_real_TetR{1,i})*time.Ts-time.Ts;
plot(t/60,RESULTS.x_real_TetR{1,i},'linewidth',2);
hold on;
end 

temp=cell2mat(RESULTS.Endtime(1,1));
tt=0;
for i=1:length(Endtime)
   tt = [tt sum(Endtime(1:i))];
end
if (length(ref_ampl)~=length(tt))
ref_ampl2=[ref_ampl ref_ampl(end)]
else 
    ref_ampl2 = ref_ampl;
end
Endtime*5/60

stairs(tt*5/60,ref_ampl2,'linewidth',2)

hold on
% xlim([0 t(end)/60]);
% ylim([0 max(ref_ampl2)])
ylabel('TetR [a.u.]')
xlabel('time[h]')
set(gca,'Fontsize',12,'FontWeight','bold');


% figure(3)
subplot(2,2,4)
for i=1:nrep
t = 0:time.Ts:length(RESULTS.U{1,i})*time.Ts-time.Ts;
plot(t/60,RESULTS.U{1,i},'Color',orange,'linewidth',2);
hold on;
end

ylim([0 1])
ylabel('IPTG [mM]')
xlabel('time[h]')
set(gca,'Fontsize',12,'FontWeight','bold');

set(gcf,'OuterPosition',[15 15 1400 700])
% clear all
%% STOCHASTIC from Matlab _ fixed time reference

% load('.\MATLAB DATA\MPC SDE 30 pts 10 rep fix ref time tend10_IPTGdecreasing.mat')
% load('.\MATLAB DATA\P SDE 30 pts 10 rep fix ref time tend10 kp016_IPTGdecreasing.mat')

CTRL = [];
TETR = [];
for i=1:10%length(RESULTS.BP)
    CTRL = [CTRL; (RESULTS.BP{i})];
end

for i=1:10%length(RESULTS.YA)
    TETR = [TETR; (RESULTS.YA{i})];
end
% 

subplot(1,2,1)
ct_vec=reshape(CTRL,size(CTRL,1)*size(CTRL,2),1);
tet_vec=reshape(TETR,size(TETR,1)*size(TETR,2),1);
N = hist3([ct_vec tet_vec],[20 20]);
N_pcolor = N';
N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
xl = linspace(min(ct_vec),max(ct_vec),size(N_pcolor,2)); % Columns of N_pcolor
yl = linspace(min(tet_vec),max(tet_vec),size(N_pcolor,1)); % Rows of N_pcolor

h = pcolor(xl,yl,N_pcolor);
colormap('hsv') % Change color scheme 
colorbar % Display colorbar
h.ZData = -max(N_pcolor(:))*ones(size(N_pcolor));
ax = gca;
ax.ZTick(ax.ZTick < 0) = [];
% ax.YLim=[min(tet_vec)-5 max(tet_vec)+5];
ax.YLim=[min(tet_vec)-5 1350];
ax.XLim=[min([0 min(ct_vec)-0.01]) max(ct_vec)+0.01*max(ct_vec)];
title('Density Map');
shading interp;
hold on
openfig('coco_bif_curve')
lin=findobj(gca,'Type', 'line');
x_bif = lin(1).XData;
y_bif=lin(1).YData;
x_coco=lin(2).XData;
y_coco=lin(2).YData;
close
z2=10*ones(length(x_coco),length(x_coco));

hold on
plot3(x_coco,y_coco,z2, '--', 'LineWidth', 2, 'Color', '[0.65,0.65,0.65]');
hold on
xlabel('IPTG [mM]')
ylabel('TetR [a.u.]')
set(gca,'Fontsize',12,'FontWeight','bold');

time=Parameters.time;

subplot(2,2,2)
for i=1:1
t = 0:time.Ts:length(RESULTS.x_real_TetR{1,i})*time.Ts-time.Ts;
plot(t/60,RESULTS.x_real_TetR{1,i},'linewidth',2,'Color',blue);
hold on;
end 

temp=cell2mat(RESULTS.Endtime(1,1));
tt=0;
for i=1:length(Endtime)
   tt = [tt sum(Endtime(1:i))];
end
if (length(ref_ampl)~=length(tt))
ref_ampl2=[ref_ampl ref_ampl(end)]
else 
    ref_ampl2 = ref_ampl;
end
Endtime*5/60
stairs(tt*5/60,ref_ampl2,'linewidth',2,'Color',gray)
ax = gca;
ax.YLim = [0 2000]

legend('TetR','TetR*', 'Location','best')
ylabel('TetR/TetR* [a.u.]')
xlabel('time[h]')
set(gca,'Fontsize',12,'FontWeight','bold');

subplot(2,2,4)
for i=1:1
t = 0:time.Ts:length(RESULTS.U{1,i})*time.Ts-time.Ts;
plot(t/60,RESULTS.U{1,i},'Color',orange,'linewidth',2);

hold on;
end

ylim([0 1])
ylabel('IPTG [mM]')
xlabel('time[h]')
set(gca,'Fontsize',12,'FontWeight','bold');

annotation('textbox', [0.05, 0.98, 0, 0], 'string', 'A','FontSize',12,'FontWeight','bold')
annotation('textbox', [0.5, 0.98, 0, 0], 'string', 'B','FontSize',12,'FontWeight','bold')
annotation('textbox', [0.5, 0.51, 0, 0], 'string', 'C','FontSize',12,'FontWeight','bold')

% set(0, 'DefaultFigureRenderer', 'painters');

%% STOCH from BSIM _ fixed time reference

% load('.\MATLAB DATA\BSIM_MPC_300pts.mat')
% load('.\MATLAB DATA\BSIM_P_300pts_dataforstochastic plot.mat')

cyan = rgb('cyan');
blue = rgb('blue');
red = rgb('red');
gray = rgb('DarkGray');
orange = rgb('orange');
green = rgb('darkgreen');
teal = rgb('teal');
orangered = rgb('OrangeRed');
darkred = rgb('DarkRed');

openfig('coco_bif_curve')
lin=findobj(gca,'Type', 'line');
x_bif = lin(1).XData;
y_bif=lin(1).YData;
x_coco=lin(2).XData;
y_coco=lin(2).YData;
close

ct_vec=bP';
tet_vec=yA';

subplot(1,2,1)

N = hist3([ct_vec tet_vec],[20 20]);
N_pcolor = N';
N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
xl = linspace(min(ct_vec),max(ct_vec),size(N_pcolor,2)); % Columns of N_pcolor
yl = linspace(min(tet_vec),max(tet_vec),size(N_pcolor,1)); % Rows of N_pcolor

h = pcolor(xl,yl,N_pcolor);
colormap('hsv') % Change color scheme 
colorbar % Display colorbar
h.ZData = -max(N_pcolor(:))*ones(size(N_pcolor));
ax = gca;
ax.ZTick(ax.ZTick < 0) = [];
% ax.YLim=[min(tet_vec)-5 max(tet_vec)+5];
ax.YLim=[min(tet_vec)-5 1350];
ax.XLim=[min([0 min(ct_vec)-0.01]) max(ct_vec)+0.01*max(ct_vec)];
title('Density Map');
shading interp;
hold on
openfig('coco_bif_curve')
lin=findobj(gca,'Type', 'line');
x_bif = lin(1).XData;
y_bif=lin(1).YData;
x_coco=lin(2).XData;
y_coco=lin(2).YData;
close
z2=10*ones(length(x_coco),length(x_coco));

hold on
plot3(x_coco,y_coco,z2, '--', 'LineWidth', 2, 'Color', '[0.65,0.65,0.65]');
hold on
xlabel('IPTG [mM]')
ylabel('TetR [a.u.]')
set(gca,'Fontsize',12,'FontWeight','bold');

%%%
% ref=control(:,3);
ref=refs;
t = t(not(ref==1))
av_TetR = av_TetR(not(ref==1))
IPTG = IPTG(not(ref==1))
ref = ref(not(ref==1))

subplot(2,2,2)
plot(t/60,av_TetR,'linewidth',2,'Color',blue);
hold on;
plot(t/60,ref,'linewidth',2,'Color',gray);
ax = gca;
ax.YLim=[0 2000];
ax.XLim=[0 t(end)/60];
legend('TetR Signal','Reference Signal', 'Location','best')

ylabel('TetR [a.u.]')
xlabel('time[h]')
set(gca,'Fontsize',12,'FontWeight','bold');

xt = ax.XTick;
xticks([xt fix(t(end)/60)]);

subplot(2,2,4)
plot(t/60,IPTG,'Color',orange,'linewidth',2);

ax = gca;
ax.YLim=[0 1];

ax.XLim=[0 t(end)/60];
xt = ax.XTick;
xticks([xt fix(t(end)/60)]);

ylabel('IPTG [mM]')
xlabel('time[h]')
set(gca,'Fontsize',12,'FontWeight','bold');

annotation('textbox', [0.05, 0.98, 0, 0], 'string', 'A','FontSize',12,'FontWeight','bold')
annotation('textbox', [0.5, 0.98, 0, 0], 'string', 'B','FontSize',12,'FontWeight','bold')
annotation('textbox', [0.5, 0.51, 0, 0], 'string', 'C','FontSize',12,'FontWeight','bold')

%% STOCHASTIC VARIABLE TIME_minimum experimental lenght MATLAB

% load('.\MATLAB DATA\SDE P 12 pts 10 rep backward kp 0016.mat')
% load('.\MATLAB DATA\MPC SDE 12 pts 10 rep var ref time tend10_IPTGdecreasing.mat')

cyan = rgb('cyan');
blue = rgb('blue');
red = rgb('red');
gray = rgb('DarkGray');
orange = rgb('orange');
green = rgb('darkgreen');
teal = rgb('teal');
orangered = rgb('OrangeRed');
darkred = rgb('DarkRed');

openfig('coco_bif_curve')
lin=findobj(gca,'Type', 'line');
x_bif = lin(1).XData;
y_bif=lin(1).YData;
x_coco=lin(2).XData;
y_coco=lin(2).YData;
close

CTRL = [];
TETR = [];
for i=1:10%length(RESULTS.BP)
    CTRL = [CTRL; (RESULTS.BP{i})];
end

for i=1:10%length(RESULTS.YA)
    TETR = [TETR; (RESULTS.YA{i})];
end

subplot(1,2,1)
ct_vec=reshape(CTRL,size(CTRL,1)*size(CTRL,2),1);
tet_vec=reshape(TETR,size(TETR,1)*size(TETR,2),1);
N = hist3([ct_vec tet_vec],[20 20]);
N_pcolor = N';
N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
xl = linspace(min(ct_vec),max(ct_vec),size(N_pcolor,2)); % Columns of N_pcolor
yl = linspace(min(tet_vec),max(tet_vec),size(N_pcolor,1)); % Rows of N_pcolor

h = pcolor(xl,yl,N_pcolor);
colormap('hsv') % Change color scheme 
colorbar % Display colorbar
h.ZData = -max(N_pcolor(:))*ones(size(N_pcolor));
ax = gca;
ax.ZTick(ax.ZTick < 0) = [];
% ax.YLim=[min(tet_vec)-5 max(tet_vec)+5];
ax.YLim=[min(tet_vec)-5 1350];
ax.XLim=[min([0 min(ct_vec)-0.01]) max(ct_vec)+0.01*max(ct_vec)];
title('Density Map');
shading interp;
hold on
openfig('coco_bif_curve')
lin=findobj(gca,'Type', 'line');
x_bif = lin(1).XData;
y_bif=lin(1).YData;
x_coco=lin(2).XData;
y_coco=lin(2).YData;
close
z2=10*ones(length(x_coco),length(x_coco));

hold on
plot3(x_coco,y_coco,z2, '--', 'LineWidth', 2, 'Color', '[0.65,0.65,0.65]');
hold on
xlabel('IPTG [mM]')
ylabel('TetR [a.u.]')
set(gca,'Fontsize',12,'FontWeight','bold');

time=Parameters.time;

mintime=sum(RESULTS.Endtime{1,1});
maxtime=sum(RESULTS.Endtime{1,1});
indexmin=1;
indexmax=1;
for i=1:10%length(RESULTS.Endtime)
    tempsum=sum(RESULTS.Endtime{1,i});
    if tempsum<mintime
        indexmin=i;
        mintime=tempsum;
    elseif tempsum>=maxtime
        indexmax=i;
        maxtime=tempsum;
    end
end

subplot(2,2,2)

t = 0:time.Ts:length(RESULTS.x_real_TetR{1,indexmin})*time.Ts-time.Ts;
plot(t/60,RESULTS.x_real_TetR{1,indexmin},'linewidth',2,'Color',blue);

hold on
temp=cell2mat(RESULTS.Endtime(1,indexmin));
tt=0;
for i=1:length(RESULTS.Endtime{1,indexmin})
   tt = [tt sum(RESULTS.Endtime{1,indexmin}(1:i))];
end
if (length(ref_ampl)~=length(tt))
ref_ampl2=[ref_ampl ref_ampl(end)]
else 
    ref_ampl2 = ref_ampl;
end
stairs(tt*5/60,ref_ampl2,'linewidth',2,'Color',gray)
ax = gca;

ax.YLim=[0 2000];
ax.XLim=[0 sum(RESULTS.Endtime{1,indexmin}*5/60)];
xt = ax.XTick;
xticks([xt fix(sum(RESULTS.Endtime{1,indexmin}*5/60))]);
% ax.XLim=[0 300];
legend('TetR','TetR*', 'Location','best')
% xlim([0 sum(RESULTS.Endtime{1,indexmax}*5/60)])
% ylim([0 2000])
ylabel('TetR/TetR* [a.u.]')
xlabel('time[h]')
set(gca,'Fontsize',12,'FontWeight','bold');

subplot(2,2,4)

t = 0:time.Ts:length(RESULTS.U{1,indexmin})*time.Ts-time.Ts;
plot(t/60,RESULTS.U{1,indexmin},'Color',orange,'linewidth',2);

ax = gca;
ax.YLim=[0 1];
% ax.XLim=[0 300];

ax.XLim=[0 sum(RESULTS.Endtime{1,indexmin}*5/60)];
xt = ax.XTick;
xticks([xt fix(sum(RESULTS.Endtime{1,indexmin}*5/60))]);
% xlim([0 sum(RESULTS.Endtime{1,indexmax}*5/60)])
% ylim([0 1])
ylabel('IPTG [mM]')
xlabel('time[h]')
set(gca,'Fontsize',12,'FontWeight','bold');


annotation('textbox', [0.05, 0.98, 0, 0], 'string', 'A','FontSize',12,'FontWeight','bold')
annotation('textbox', [0.5, 0.98, 0, 0], 'string', 'B','FontSize',12,'FontWeight','bold')
annotation('textbox', [0.5, 0.51, 0, 0], 'string', 'C','FontSize',12,'FontWeight','bold')

%% STOCHASTIC VARIABLE TIME_ BSIM

% load('.\MATLAB DATA\BSIM_P_120pts_varendtime_1500maxref_alldataforstochplot.mat')
% load('.\MATLAB DATA\BSIM_MPC_120pts_varendtime_alldataforstochplot.mat')

cyan = rgb('cyan');
blue = rgb('blue');
red = rgb('red');
gray = rgb('DarkGray');
orange = rgb('orange');
green = rgb('darkgreen');
teal = rgb('teal');
orangered = rgb('OrangeRed');
darkred = rgb('DarkRed');

openfig('coco_bif_curve')
lin=findobj(gca,'Type', 'line');
x_bif = lin(1).XData;
y_bif=lin(1).YData;
x_coco=lin(2).XData;
y_coco=lin(2).YData;
close

ct_vec=bP';
tet_vec=yA';

subplot(1,2,1)

N = hist3([ct_vec tet_vec],[20 20]);
N_pcolor = N';
N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
xl = linspace(min(ct_vec),max(ct_vec),size(N_pcolor,2)); % Columns of N_pcolor
yl = linspace(min(tet_vec),max(tet_vec),size(N_pcolor,1)); % Rows of N_pcolor

h = pcolor(xl,yl,N_pcolor);
colormap('hsv') % Change color scheme 
colorbar % Display colorbar
h.ZData = -max(N_pcolor(:))*ones(size(N_pcolor));
ax = gca;
ax.ZTick(ax.ZTick < 0) = [];
% ax.YLim=[min(tet_vec)-5 max(tet_vec)+5];
ax.YLim=[min(tet_vec)-5 1350];
ax.XLim=[min([0 min(ct_vec)-0.01]) max(ct_vec)+0.01*max(ct_vec)];
title('Density Map');
shading interp;
hold on
openfig('coco_bif_curve')
lin=findobj(gca,'Type', 'line');
x_bif = lin(1).XData;
y_bif=lin(1).YData;
x_coco=lin(2).XData;
y_coco=lin(2).YData;
close
z2=10*ones(length(x_coco),length(x_coco));

hold on
plot3(x_coco,y_coco,z2, '--', 'LineWidth', 2, 'Color', '[0.65,0.65,0.65]');
hold on
xlabel('IPTG [mM]')
ylabel('TetR [a.u.]')
set(gca,'Fontsize',12,'FontWeight','bold');

ref=control(:,3);
t = t(not(ref==1))
av_TetR = av_TetR(not(ref==1))
IPTG = IPTG(not(ref==1))
ref = ref(not(ref==1))
% ref=refs;
subplot(2,2,2)
plot(t/60,av_TetR,'linewidth',2,'Color',blue);
hold on;
plot(t/60,ref,'linewidth',2,'Color',gray);
ax = gca;
ax.YLim=[0 2000];
ax.XLim=[0 t(end)/60];
legend('TetR Signal','Reference Signal', 'Location','best')

ylabel('TetR [a.u.]')
xlabel('time[h]')
set(gca,'Fontsize',12,'FontWeight','bold');

xt = ax.XTick;
xticks([xt fix(t(end)/60)]);

subplot(2,2,4)
plot(t/60,IPTG,'Color',orange,'linewidth',2);

ax = gca;
ax.YLim=[0 1];

ax.XLim=[0 t(end)/60];
xt = ax.XTick;
xticks([xt fix(t(end)/60)]);

ylabel('IPTG [mM]')
xlabel('time[h]')
set(gca,'Fontsize',12,'FontWeight','bold');


annotation('textbox', [0.05, 0.98, 0, 0], 'string', 'A','FontSize',12,'FontWeight','bold')
annotation('textbox', [0.5, 0.98, 0, 0], 'string', 'B','FontSize',12,'FontWeight','bold')
annotation('textbox', [0.5, 0.51, 0, 0], 'string', 'C','FontSize',12,'FontWeight','bold')
%%
set(gcf,'OuterPosition',[15 15 1400 700])
set(0, 'DefaultFigureRenderer', 'painters');
saveas(gcf,strcat('.\figfiles\',strname),'fig');saveas(gcf,strcat('.\figfiles\',strname),'png');
saveas(gcf,strcat('.\figfiles\',strname),'svg');

%% COMPUTED BIFURCATION CUIRVES OUT OF THE ESTIMATED PARAMETERS
% load('.\MATLAB DATA\allDATA_toplot_scurves.mat')
% load('.\MATLAB DATA\allDATA_toplot_scurvesWITHconstraintsonIPTG.mat')

subplot(2,3,1)
hold on
plot(measured_BR_ODE, xdata_BR_ODE, '.', 'MarkerSize', 6, 'Color', blue);

hold on;
plot(x_coco, y_coco, '--', 'LineWidth', 2, 'Color', gray);
plot(IP2_BR_ODE,TetR_vect_BR_ODE, 'LineWidth', 2, 'Color', blue);

fig = gcf; % current figure
ax = fig.CurrentAxes;% current axes
ax.YLim=[0 1];
% ax.XLim=[0 300];
ax.TickLength = [0.02 0.02];
ax.YLim = [0 1500];
ax.XLim = [0 1.5];
ylabel('TetR');
set(ax,'Fontsize',12,'FontWeight','bold');
title('Parameter Sweeps- Matlab ODE')



subplot(2,3,2)
plot(measured_BR_SDE, xdata_BR_SDE, '.', 'MarkerSize', 6, 'Color', blue);
hold on;
plot(x_coco, y_coco, '--', 'LineWidth', 2, 'Color', gray);
plot(IP2_BR_SDE,TetR_vect_BR_SDE, 'LineWidth', 2, 'Color', blue);

fig = gcf; % current figure
ax = fig.CurrentAxes;% current axes
ax.FontSize = 10;
% ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.YLim = [0 1500];
ax.XLim = [0 1.5];
set(ax,'Fontsize',12,'FontWeight','bold');
title('Parameter Sweeps- Matlab SDE')


subplot(2,3,3)
plot(measured_SWEEPS_BSIM, xdata_SWEEPS_BSIM, '.', 'MarkerSize', 6, 'Color', blue);
hold on;
plot(x_coco, y_coco, '--', 'LineWidth', 2, 'Color', gray);
plot(IP2_SWEEPS_BSIM,TetR_vect_SWEEPS_BSIM, 'LineWidth', 2, 'Color', blue);

fig = gcf; % current figure
ax = fig.CurrentAxes;% current axes
ax.FontSize = 10;
% ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.YLim = [0 1500];
ax.XLim = [0 1.5];

title('Parameter Sweeps- BSIM')
set(ax,'Fontsize',12,'FontWeight','bold');
legend('Measured','IPTG_{BF}=f(TetR,\theta)','location','best')


subplot(2,3,4)
hold on
% plot(measured_BR_ODE, xdata_BR_ODE, '.', 'MarkerSize', 6, 'Color', blue);
plot(measured_MPCODE, xdata_MPCODE, '.', 'MarkerSize', 6, 'Color', red);
plot(measured_PODE, xdata_PODE, '.', 'MarkerSize', 6, 'Color', green);
hold on;
plot(x_coco, y_coco, '--', 'LineWidth', 2, 'Color', gray);
% plot(IP2_BR_ODE,TetR_vect_BR_ODE, 'LineWidth', 2, 'Color', blue);
plot(IP2_MPCODE,TetR_vect_MPCODE, 'LineWidth', 2, 'Color', red);
plot(IP2_PODE,TetR_vect_PODE, 'LineWidth', 2, 'Color', green);
title('CBC Matlab ODE')

% legend('MPC measured','P measured','Numerical Continuation','IPTG_{MPC}=f(TetR,\theta)','IPTG_{P}=f(TetR,\theta)')
fig = gcf; % current figure
ax = fig.CurrentAxes;% current axes
ax.FontSize = 10;
% ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.YLim = [0 1500];
ax.XLim = [0 1.5];
xlabel('IPTG');
set(ax,'Fontsize',12,'FontWeight','bold');


subplot(2,3,5)
hold on
% plot(measured_BR_SDE, xdata_BR_SDE, '.', 'MarkerSize', 6, 'Color', blue);
plot(measured_MPCSDE, xdata_MPCSDE, '.', 'MarkerSize', 6, 'Color', red);
plot(measured_PSDE, xdata_PSDE, '.', 'MarkerSize', 6, 'Color', green);
hold on;
plot(x_coco, y_coco, '--', 'LineWidth', 2, 'Color', gray);
% plot(IP2_BR_SDE,TetR_vect_BR_SDE, 'LineWidth', 2, 'Color', blue);
plot(IP2_MPCSDE,TetR_vect_MPCSDE, 'LineWidth', 2, 'Color', red);
plot(IP2_PSDE,TetR_vect_PSDE, 'LineWidth', 2, 'Color', green);
title('CBC Matlab SDE')

% legend('MPC measured','P measured','Numerical Continuation','IPTG_{MPC}=f(TetR,\theta)','IPTG_{P}=f(TetR,\theta)')
fig = gcf; % current figure
ax = fig.CurrentAxes;% current axes
ax.FontSize = 10;
% ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.YLim = [0 1500];
ax.XLim = [0 1.5];
xlabel('IPTG');
set(ax,'Fontsize',12,'FontWeight','bold');


subplot(2,3,6)
hold on
plot(measured_BSIM_MPC, xdata_BSIM_MPC, '.', 'MarkerSize', 6, 'Color', red);
plot(measured_BSIM_P, xdata_BSIM_P, '.', 'MarkerSize', 6, 'Color', green);
hold on;
plot(x_coco, y_coco, '--', 'LineWidth', 2, 'Color', gray);
plot(IP2_BSIM_MPC,TetR_vect_BSIM_MPC, 'LineWidth', 2, 'Color', red);
plot(IP2_BSIM_P,TetR_vect_BSIM_P, 'LineWidth', 2, 'Color', green);
title('BSIM SDE')
legend('Measured MPC','Measured P','Numerical Continuation','IPTG_{MPC}=f(TetR,\theta)','IPTG_{P}=f(TetR,\theta)')
% legend('Measured MPC','Numerical Continuation','Calibrated Model')
fig = gcf; % current figure
ax = fig.CurrentAxes;% current axes
ax.FontSize = 10;
% ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.YLim = [0 1500];
ax.XLim = [0 1.5];
xlabel('IPTG');ylabel('TetR');
set(ax,'Fontsize',12,'FontWeight','bold');

annotation('textbox', [0.1, 0.98, 0, 0], 'string', 'A','FontSize',12,'FontWeight','bold')
annotation('textbox', [0.38, 0.98, 0, 0], 'string', 'B','FontSize',12,'FontWeight','bold')
annotation('textbox', [0.66, 0.98, 0, 0], 'string', 'C','FontSize',12,'FontWeight','bold')
annotation('textbox', [0.1, 0.51, 0, 0], 'string', 'D','FontSize',12,'FontWeight','bold')
annotation('textbox', [0.38, 0.51, 0, 0], 'string', 'E','FontSize',12,'FontWeight','bold')
annotation('textbox', [0.66, 0.51, 0, 0], 'string', 'F','FontSize',12,'FontWeight','bold')
set(gcf,'OuterPosition',[15 15 1500 800])