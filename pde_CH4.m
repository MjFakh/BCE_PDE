function [c,f,s] = pde_CH4(x,t,u,dudx)
global k_SO4 xmesh tspan RC Oxygen k_O2 DCH4 v_burial_Fluid CH4init k1_AOM k1_CH4O2
global Alpha_Bioirrig z_sed DOC k_DOC alpha_DOC k_CH4O2 Sulfate k_AOM

v_burial_f = interp1(z_sed,v_burial_Fluid,x);
Alpha_Bioirrig_1 = interp1(z_sed,Alpha_Bioirrig,x);

[x1,t1]=meshgrid(x,t);
A1_DOC = interp2(xmesh,tspan,alpha_DOC.*(DOC./(DOC+k_DOC)),x1,t1);
RC1 = interp2(xmesh,tspan,RC,x1,t1);
SO4_1 = interp2(xmesh,tspan,Sulfate,x1,t1);
O2 = interp2(xmesh,tspan,Oxygen,x1,t1);
Inh = (k_O2./(O2+k_O2));
Inh1 = (k_SO4/(SO4_1+k_SO4));

c = 1;

f = DCH4.*dudx;

S2 = - v_burial_f.* dudx + RC1.*A1_DOC.*Inh.*Inh1.* 1E9 - (k_AOM.*u(1).*(SO4_1/(SO4_1+k1_AOM))) ...
     - (k_CH4O2.*u(1).*(O2/(O2+k1_CH4O2))) + (Alpha_Bioirrig_1.*(CH4init-u(1))); % umol/l/year

s = S2;  %umol/l/year

end