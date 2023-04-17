
function [pl,ql,pr,qr] = organicbc(xl,ul,xr,ur,t)
global NPP v_burial poros rho Bioturb BE alpha_bioturb

NPP1 = BE * NPP * 1E-4; %gram/cm2/year
v_burial1 = v_burial(1,1);  %cm/year
poros1 = poros(1,1); 
A1 = rho * (1-poros1);
Bioturb1 = Bioturb(1,1);%.*alpha_bioturb;

pl = A1 * ul * v_burial1 - NPP1;
ql = - Bioturb1 * A1;
pr = 0;
qr = 1;

end

