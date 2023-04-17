function [c,f,s] = pdefun2(x,t,u,dudx)
global DHCO3 xmesh tspan Alpha_Bioirrig R_SRR HCO3init RC
global v_burial_Fluid z_sed R1_carb R_HS_Ox R_respi alpha_DOC k_DOC DOC R1_carb_form R1_carb_disso

v_burial_f = interp1(z_sed,v_burial_Fluid,x);
Alpha_Bioirrig_1 = interp1(z_sed,Alpha_Bioirrig,x);

[x1,t1]=meshgrid(x,t);
%R_SRR1 = interp2(xmesh,tspan,R_SRR,x1,t1);R_respi
R_respi1 = interp2(xmesh,tspan,R_respi,x1,t1);
R_HS_Ox1 = interp2(xmesh,tspan,R_HS_Ox,x1,t1);
RC1 = interp2(xmesh,tspan,RC,x1,t1);
R1_carb1 = interp2(xmesh,tspan,R1_carb_form+R1_carb_disso,x1,t1);
A1_DOC = interp2(xmesh,tspan,alpha_DOC.*(DOC./(DOC+k_DOC)),x1,t1);

c = 1;

f = DHCO3.*dudx;

S2 = - v_burial_f.* dudx + A1_DOC.* RC1.* 1E9 - R_respi1 ...
     - 2.*R1_carb1 - 2.*R_HS_Ox1 + (Alpha_Bioirrig_1.*(HCO3init-u(1))); % umol/l/year  
s = S2;  %umol/l/year 

end


