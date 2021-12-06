function [IPTG] = IPTG_fun(Params,TetR)
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

glp= Params(14);

ktp = Params(15);

gtp = Params(16);

%fixed
aTc= 25;
gtm=1.386e-1;
glm=1.386e-1;

%% ANALYTICAL FUNCTIONS COMPUTED WITH MAPLE:

mRNAt = gtp * TetR / ktp;
mRNAl = (klm0 * (TetR / (1 + (aTc / thetaAtc) ^ etaAtc) / thetaTet) ^ etaTet + klm + klm0) / glm / (1 + (TetR / (1 + (aTc / thetaAtc) ^ etaAtc) / thetaTet) ^ etaTet); 
LacI = klp * (klm0 * (TetR / (1 + (aTc / thetaAtc) ^ etaAtc) / thetaTet) ^ etaTet + klm + klm0) / glp / glm / (1 + (TetR / (1 + (aTc / thetaAtc) ^ etaAtc) / thetaTet) ^ etaTet);
IPTG = ((-(TetR / (1 + (aTc / thetaAtc) ^ etaAtc) / thetaTet) ^ etaTet * glm * glp * thetaLac - glp * glm * thetaLac + (-(gtm * gtp * TetR - ktm * ktp - ktm0 * ktp) / (gtm * gtp * TetR - ktm0 * ktp)) ^ (-0.1e1 / etaLac) * klm0 * klp * (TetR / (1 + (aTc / thetaAtc) ^ etaAtc) / thetaTet) ^ etaTet + (-(gtm * gtp * TetR - ktm * ktp - ktm0 * ktp) / (gtm * gtp * TetR - ktm0 * ktp)) ^ (-0.1e1 / etaLac) * klm * klp + (-(gtm * gtp * TetR - ktm * ktp - ktm0 * ktp) / (gtm * gtp * TetR - ktm0 * ktp)) ^ (-0.1e1 / etaLac) * klm0 * klp) / glp / glm / thetaLac / (1 + (TetR / (1 + (aTc / thetaAtc) ^ etaAtc) / thetaTet) ^ etaTet)) ^ (0.1e1 / etaIptg) * thetaIptg;

end