function [c,f,s] = pde_NH4(x,t,u,dudx)
global DO2 xmesh tspan RC Alpha_Bioirrig NH4init v_burial_Fluid 
global z_sed DOC alpha_DOC k_DOC Knitrif N_C_ratio Oxygen


v_burial_f = interp1(z_sed,v_burial_Fluid,x);
Alpha_Bioirrig_1 = interp1(z_sed,Alpha_Bioirrig,x);

[x1,t1]=meshgrid(x,t);
A1_DOC = interp2(xmesh,tspan,alpha_DOC.*(DOC./(DOC+k_DOC)),x1,t1);
RC1 = interp2(xmesh,tspan,RC,x1,t1);
O2 = interp2(xmesh,tspan,Oxygen,x1,t1);

c = 1;

f = DO2.*dudx;

S2 = - v_burial_f.* dudx + N_C_ratio.* RC1.*A1_DOC.* 1E9 - ...
      (Knitrif.*u(1).*O2) + (Alpha_Bioirrig_1.*(NH4init-u(1))); % umol/l/year 

s = S2;  %umol/l/year

end