function [c,f,s] = pde_NO3(x,t,u,dudx)
global DO2 k_NO3 xmesh tspan RC Alpha_Bioirrig NO3init v_burial_Fluid 
global z_sed DOC alpha_DOC k_DOC Knitrif Oxygen NH4 k_O2


v_burial_f = interp1(z_sed,v_burial_Fluid,x);
Alpha_Bioirrig_1 = interp1(z_sed,Alpha_Bioirrig,x);

[x1,t1]=meshgrid(x,t);
A1_DOC = interp2(xmesh,tspan,alpha_DOC.*(DOC./(DOC+k_DOC)),x1,t1);
RC1 = interp2(xmesh,tspan,RC,x1,t1);
NH4_1 = interp2(xmesh,tspan,NH4,x1,t1);
O2 = interp2(xmesh,tspan,Oxygen,x1,t1);
Inh = (k_O2./(O2+k_O2));

c = 1;

f = DO2.*dudx;

S2 = - v_burial_f.* dudx - 0.8.* RC1.*A1_DOC.*Inh.* (u(1)/(u(1)+k_NO3)) * 1E9 + ...
      (Knitrif.*NH4_1.*O2) + (Alpha_Bioirrig_1.*(NO3init-u(1))); % umol/l/year 

s = S2;  %umol/l/year

end