function [c,f,s] = pdefun4(x,t,u,dudx)
global DHCO3 xmesh tspan RC Alpha_Bioirrig DICinit R1_carb v_burial_Fluid z_sed alpha_DOC DOC k_DOC
global R1_carb_form R1_carb_disso R_HS_Ox

v_burial_f = interp1(z_sed,v_burial_Fluid,x);
Alpha_Bioirrig_1 = interp1(z_sed,Alpha_Bioirrig,x);

[x1,t1]=meshgrid(x,t);
RC1 = interp2(xmesh,tspan,RC,x1,t1);
R1_carb1 = interp2(xmesh,tspan,R1_carb_form+R1_carb_disso,x1,t1);
A1_DOC = interp2(xmesh,tspan,alpha_DOC.*(DOC./(DOC+k_DOC)),x1,t1);
R_HS_Ox1 = interp2(xmesh,tspan,R_HS_Ox,x1,t1);


c = 1;

f = DHCO3.*dudx;

S2 = - v_burial_f.* dudx + A1_DOC.* RC1.* 1E9 - R1_carb1 + (Alpha_Bioirrig_1.*(DICinit-u(1))); % umol/l/year + CO2_Mineral.* R_dissol1
s = S2;  %umol/l/year

end