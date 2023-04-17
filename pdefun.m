
function [c,f,s] = pdefun(x,t,u,dudx)
global DO2 k_O2 xmesh tspan RC Alpha_Bioirrig O2init v_burial_Fluid 
global O2_root z_sed DOC alpha_DOC k_DOC Kreox Sulfide

v_burial_f = interp1(z_sed,v_burial_Fluid,x);
O2_root_PDE = interp1(xmesh,O2_root,x);
Alpha_Bioirrig_1 = interp1(z_sed,Alpha_Bioirrig,x);

[x1,t1]=meshgrid(x,t);
A1_DOC = interp2(xmesh,tspan,alpha_DOC.*(DOC./(DOC+k_DOC)),x1,t1);
RC1 =  interp2(xmesh,tspan,RC,x1,t1);
H2S = interp2(xmesh,tspan,Sulfide,x1,t1);

c = 1;

f = DO2.*dudx;

S2 = - v_burial_f.* dudx - A1_DOC * RC1 * (u(1)/(u(1)+k_O2)) * 1E9 ...
     - (Kreox.*u(1).*H2S) + (Alpha_Bioirrig_1.*(O2init-u(1))) + O2_root_PDE; % umol/l/year
s = S2;  %umol/l/year

end






