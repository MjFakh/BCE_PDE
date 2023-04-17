function [c,f,s] = pde_Fe(x,t,u,dudx)
global xmesh tspan RC Oxygen k_O2 DH2S v_burial_Fluid Feinit FeooH kFeOx
global Alpha_Bioirrig z_sed DOC k_DOC alpha_DOC Sulfide kFeS KFEMonod

v_burial_f = interp1(z_sed,v_burial_Fluid,x);
FeOx = interp1(z_sed,FeooH,x);
Alpha_Bioirrig_1 = interp1(z_sed,Alpha_Bioirrig,x);

[x1,t1]=meshgrid(x,t);
A1_DOC = interp2(xmesh,tspan,alpha_DOC.*(DOC./(DOC+k_DOC)),x1,t1);
RC1 = interp2(xmesh,tspan,RC,x1,t1);
O2 = interp2(xmesh,tspan,Oxygen,x1,t1);
H2S = interp2(xmesh,tspan,Sulfide,x1,t1);
Inh = (k_O2./(O2+k_O2));

c = 1;

f = DH2S.*dudx;

S2 = - v_burial_f.* dudx + 4 * RC1.*A1_DOC.*Inh.* (FeOx./(FeOx+KFEMonod)) * 1E9 ...
     - (kFeS.*u(1).*H2S) - (kFeOx.*u(1).*O2) + (Alpha_Bioirrig_1.*(Feinit-u(1))); % umol/l/year

s = S2;  %umol/l/year

end

