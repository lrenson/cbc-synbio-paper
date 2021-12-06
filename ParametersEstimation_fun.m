clear all;
close all;
clc;

%% Parameter Evaluation
%% NUMERICAL SOL
%original parameters
klm0=3.20e-2;

klm=8.30;

thetaAtc=11.65;

etaAtc=2.00;

thetaTet=30.00;

etaTet=2.00;

glm=1.386e-1;

ktm0=1.19e-1;

ktm=2.06;   

thetaIptg=9.06e-2;

etaIptg=2.00;

thetaLac=31.94;

etaLac=2.00;

gtm=1.386e-1;

klp=9.726e-1;

glp=1.65e-2;

ktp=1.170;

gtp=1.65e-2;


%%
% % clear all
% % clc
% load('.\bruteforceODE30points.mat','xdata','measured')
% load('.\bruteforceSDEnew.mat','xdata','xdata')
% load('.\MATLAB DATA\BSIM_SWEEPS_300pts.mat')
% load('.\MATLAB DATA\MPC ODE 30 pts fix ref time tend10_IPTGdecreasing.mat', 'RESULTS')

%TO MAKE THE CODE WORK WE WANT 'measured' to contain measured IPTG values
%and 'xdata' to contain measured TetR data


% measured = [];
% xdata = [];

%from Matlab data
% for i=1:length(RESULTS.BP)
%     measured = [measured RESULTS.BP{1,i}(2:end)];
%     xdata = [xdata RESULTS.YA{1,i}(2:end)];
% end

%from BSIM data
% measured = control_vect;
% xdata = output_vect(:,1)';


% GA
[Params, cost_val] = Optimizer(14, measured, xdata);

%%
klm0 = Params(1);

klm = Params(2);

thetaAtc = Params(3);

etaAtc = Params(4);

thetaTet = Params(5);

etaTet = Params(6);

ktm0 = Params(7);

ktm = Params(8);   

thetaIptg = Params(9);

etaIptg = Params(10);

thetaLac = Params(11);

etaLac = Params(12);

klp = Params(13);

ktp = Params(14);


%fixed
aTc= 25;

gtm=1.386e-1;
glm=1.386e-1;
glp=1.65e-2;
gtp=1.65e-2;


siz = 100;
TetR_vect = linspace(min(xdata),max(xdata),siz);
TetR = TetR_vect;
IP2 =((-(TetR ./ (1 + (aTc ./ thetaAtc) .^ etaAtc) ./ thetaTet) .^ etaTet .* glm .* glp .* thetaLac - glp .* glm .* thetaLac + (-(gtm .* gtp .* TetR - ktm .* ktp - ktm0 .* ktp) ./ (gtm .* gtp .* TetR - ktm0 .* ktp)) .^ (-0.1e1 ./ etaLac) .* klm0 .* klp .* (TetR ./ (1 + (aTc ./ thetaAtc) .^ etaAtc) ./ thetaTet) .^ etaTet + (-(gtm .* gtp .* TetR - ktm .* ktp - ktm0 .* ktp) ./ (gtm .* gtp .* TetR - ktm0 .* ktp)) .^ (-0.1e1 ./ etaLac) .* klm .* klp + (-(gtm .* gtp .* TetR - ktm .* ktp - ktm0 .* ktp) ./ (gtm .* gtp .* TetR - ktm0 .* ktp)) .^ (-0.1e1 ./ etaLac) .* klm0 .* klp) ./ glp ./ glm ./ thetaLac ./ (1 + (TetR ./ (1 + (aTc ./ thetaAtc) .^ etaAtc) ./ thetaTet) .^ etaTet)) .^ (0.1e1 ./ etaIptg) .* thetaIptg;
LacI = klp .* (klm0 .* (TetR ./ (1 + (aTc ./ thetaAtc) .^ etaAtc) ./ thetaTet) .^ etaTet + klm + klm0) ./ glp ./ glm ./ (1 + (TetR ./ (1 + (aTc ./ thetaAtc) .^ etaAtc) ./ thetaTet) .^ etaTet);

%delete imag sol
% jj=1;
% for s=1:length(IP2)
%     
%     if isreal(IP2(s))
%         IPTG_r(jj)=IP2(s);
%         TetR_r(jj)=TetR(s);
%         jj=jj+1;
%     end
%    
% end
% 
% IP2=IPTG_r;
% TetR=TetR_r;

% end

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
figure()
hold on
plot(measured, xdata, '.', 'MarkerSize', 8, 'Color', red);%8
hold on;
plot(x_coco, y_coco, '--', 'LineWidth', 2, 'Color', gray);
plot(IP2,TetR, 'LineWidth', 2, 'Color', red);
plot(measured, xdata, '.', 'MarkerSize', 10, 'Color', blue);%8

xlim([0 1])
