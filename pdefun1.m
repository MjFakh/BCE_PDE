function [c,f,s] = pdefun1(x,t,u,dudx)
global k_SO4 xmesh tspan RC Oxygen k_O2 DSO4 v_burial_Fluid SO4init 
global Alpha_Bioirrig z_sed DOC k_DOC alpha_DOC Sulfide Kreox

v_burial_f = interp1(z_sed,v_burial_Fluid,x);
Alpha_Bioirrig_1 = interp1(z_sed,Alpha_Bioirrig,x);

[x1,t1]=meshgrid(x,t);
A1_DOC = interp2(xmesh,tspan,alpha_DOC.*(DOC./(DOC+k_DOC)),x1,t1);
RC1 = interp2(xmesh,tspan,RC,x1,t1);
O2 = interp2(xmesh,tspan,Oxygen,x1,t1);
Sulfide_1 = interp2(xmesh,tspan,Sulfide,x1,t1);
Inh = (k_O2./(O2+k_O2));

c = 1;

f = DSO4.*dudx;

S2 = - v_burial_f.* dudx - 0.5 * RC1.*A1_DOC.*Inh.* (u(1)/(u(1)+k_SO4)) * 1E9 + ...
      (Kreox.*Sulfide_1.*O2) + (Alpha_Bioirrig_1.*(SO4init-u(1))); % umol/l/year 
s = S2;  %umol/l/year

end