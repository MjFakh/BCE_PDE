function [c,f,s] = pde_HS(x,t,u,dudx)
global k_SO4 xmesh tspan RC Oxygen k_O2 DH2S v_burial_Fluid HSinit
global Alpha_Bioirrig z_sed DOC k_DOC alpha_DOC Iron Kreox kFeS Sulfate

v_burial_f = interp1(z_sed,v_burial_Fluid,x);
Alpha_Bioirrig_1 = interp1(z_sed,Alpha_Bioirrig,x);

[x1,t1]=meshgrid(x,t);
A1_DOC = interp2(xmesh,tspan,alpha_DOC.*(DOC./(DOC+k_DOC)),x1,t1);
RC1 = interp2(xmesh,tspan,RC,x1,t1);
SO4_1 = interp2(xmesh,tspan,Sulfate,x1,t1);
O2 = interp2(xmesh,tspan,Oxygen,x1,t1);
Fe = interp2(xmesh,tspan,Iron,x1,t1);
Inh = (k_O2./(O2+k_O2));

c = 1;

f = DH2S.*dudx;

S2 = - v_burial_f.* dudx + 0.5 * RC1.*A1_DOC.*Inh.* (SO4_1/(SO4_1+k_SO4)) * 1E9 - (Kreox.*u(1).*O2) ...
     - (kFeS.*u(1).*Fe) + (Alpha_Bioirrig_1.*(HSinit-u(1))); % umol/l/year

s = S2;  %umol/l/year

end

