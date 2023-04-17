function [pl,ql,pr,qr] = CaCO3bc(xl,ul,xr,ur,t)
global F_CaCO3 v_burial Bioturb alpha_bioturb

NPP2 = F_CaCO3 * 1E-4 * 1E6; % umol/cm2/year
v_burial1 = v_burial(1,1);  %cm/year
Bioturb1 = Bioturb(1,1);%.*alpha_bioturb;

pl = 1E-3 * ul * v_burial1 - NPP2; % umol/cm2/year
ql = - Bioturb1 * 1E-3;
pr = 0;
qr = 1;

end