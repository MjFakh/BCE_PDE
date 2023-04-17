function [c,f,s] = pdefun3(x,t,u,dudx)
global DHCO3 xmesh tspan Alpha_Bioirrig R_respi H2CO3init v_burial_Fluid z_sed

v_burial_f = interp1(z_sed,v_burial_Fluid,x);
Alpha_Bioirrig_1 = interp1(z_sed,Alpha_Bioirrig,x);

[x1,t1]=meshgrid(x,t);
R_respi1 = interp2(xmesh,tspan,R_respi,x1,t1);

c = 1;

f = DHCO3.*dudx;


S2 = - v_burial_f.* dudx + R_respi1 + (Alpha_Bioirrig_1.*(H2CO3init-u(1))); % umol/l/year
s = S2;  %umol/l/year

end