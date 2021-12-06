function [Params, cost_val]= Optimizer(n, measured, xdata)
%% Sim parameters
N = n;
TetR_vect = xdata;
%% lower and upper bounds 
% % CONSTRAINTS ON U
lb = [3.20e-2-0.5*3.20e-2, 8.30-0.5*8.30, 11.65-0.5*11.65, 2.00-0.5*2.00,...
    30.00-0.5*30.00, 2.00-0.5*2.00, 1.19e-1-0.5*1.19e-1,...
    2.06-0.5*2.06, 9.06e-2-0.5*9.06e-2, 2.00-0.5*2.00, 31.94-0.5*31.94,...
    2.00-0.5*2.00, 9.726e-1-0.5*9.726e-1, 1.170-0.5*1.170];
ub = [3.20e-2+0.5*3.20e-2, 8.30+0.5*8.30, 11.65+0.5*11.65, 2.00+0.5*2.00,...
    30.00+0.5*30.00, 2.00+0.5*2.00, 1.19e-1+0.5*1.19e-1,...
    2.06+0.5*2.06, 9.06e-2+0.5*9.06e-2, 2.00+0.5*2.00, 31.94+0.5*31.94,...
    2.00+0.5*2.00, 9.726e-1+0.5*9.726e-1, 1.170+0.5*1.170];

  opt = optimoptions('ga','Display','iter');
  opt.UseParallel = true;
  opt.MaxGenerations = 500*2;%500*N;
  opt.FunctionTolerance = 1e-7;
 [a_opt,fval,exitflag,output,population,scores] = ga(@(a) costfunction(a, measured, TetR_vect),N,[],[],[],[],lb,ub,[],opt);

 cost = costfunction(a_opt, measured,TetR_vect);
 Params = a_opt;
 cost_val=cost;
 
end

%%
function [cost] = costfunction(a, measured, xdata)
%COST FUNCTION
%%%%% discrete model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
klm0 = a(1);

klm = a(2);

thetaAtc = a(3);

etaAtc = a(4);

thetaTet = a(5);

etaTet = a(6);

ktm0 = a(7);

ktm = a(8);   

thetaIptg = a(9);

etaIptg = a(10);

thetaLac = a(11);

etaLac = a(12);

klp = a(13);

ktp = a(14);

cost = 0;

aTc= 25;

gtm=1.386e-1;
glm=1.386e-1;
glp=1.65e-2;
gtp=1.65e-2;

TetR = xdata;

mRNAt = gtp .* TetR ./ ktp;
mRNAl = (klm0 .* (TetR ./ (1 + (aTc ./ thetaAtc) .^ etaAtc) ./ thetaTet) .^ etaTet + klm + klm0) ./ glm ./ (1 + (TetR ./ (1 + (aTc ./ thetaAtc) .^ etaAtc) ./ thetaTet) .^ etaTet);
LacI = klp .* (klm0 .* (TetR ./ (1 + (aTc ./ thetaAtc) .^ etaAtc) ./ thetaTet) .^ etaTet + klm + klm0) ./ glp ./ glm ./ (1 + (TetR ./ (1 + (aTc ./ thetaAtc) .^ etaAtc) ./ thetaTet) .^ etaTet);
IPTG = ((-(TetR ./ (1 + (aTc ./ thetaAtc) .^ etaAtc) ./ thetaTet) .^ etaTet .* glm .* glp .* thetaLac - glp .* glm .* thetaLac + (-(gtm .* gtp .* TetR - ktm .* ktp - ktm0 .* ktp) ./ (gtm .* gtp .* TetR - ktm0 .* ktp)) .^ (-0.1e1 ./ etaLac) .* klm0 .* klp .* (TetR ./ (1 + (aTc ./ thetaAtc) .^ etaAtc) ./ thetaTet) .^ etaTet + (-(gtm .* gtp .* TetR - ktm .* ktp - ktm0 .* ktp) ./ (gtm .* gtp .* TetR - ktm0 .* ktp)) .^ (-0.1e1 ./ etaLac) .* klm .* klp + (-(gtm .* gtp .* TetR - ktm .* ktp - ktm0 .* ktp) ./ (gtm .* gtp .* TetR - ktm0 .* ktp)) .^ (-0.1e1 ./ etaLac) .* klm0 .* klp) ./ glp ./ glm ./ thetaLac ./ (1 + (TetR ./ (1 + (aTc ./ thetaAtc) .^ etaAtc) ./ thetaTet) .^ etaTet)) .^ (0.1e1 ./ etaIptg) .* thetaIptg;

%%
%%NEW LINES

% [IPTGr,TetRr]=real_IPTG(IPTG,TetR);
% TetRr=TetRr./max(xdata);
% 
% TetR=xdata/max(xdata);
% 
% for l=1:length(measured)
%    
%     cost=cost+min_distance(IPTGr,TetRr,measured(l),TetR(l));
%     
% end

%%END NEW LINES
%%

if (isreal(IPTG) && isreal(LacI) && isreal(mRNAl) && isreal(mRNAt) && all(IPTG>=0) && all(LacI>=0) && all(mRNAl>=0) && all(mRNAt>=0))
    for l=1:length(measured)
        
        cost=cost+min_distance(IPTG,TetR,measured(l),TetR(l));
        
    end
else
    cost = inf;
end

end



%% NEW fcn


function [IPTG_r,TetR_r] = real_IPTG(IPTG,TetR)

jj=1;
IPTG_r=[];

for s=1:length(IPTG)
       
    if isreal(IPTG(s))
        IPTG_r(jj)=IPTG(s);
        TetR_r(jj)=TetR(s);
        jj=jj+1;
    end
      
end

if isempty(IPTG_r)
   IPTG_r=50;
   TetR_r=50;
end

end


function [md] = min_distance(IPTG_r,TetR_r,m_i,m_t)

dist_m=inf;
for kk=1:length(IPTG_r)
    
    dist=norm([IPTG_r(kk) TetR_r(kk)]-[m_i m_t]);
    if dist<dist_m
           dist_m = dist; 
%            ind=kk;
    end
    
    md=dist_m;

end

end
