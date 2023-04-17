function [c,f,s] = DOCfun(x,t,u,dudx)
global DDOC k_DOC xmesh tspan RC Alpha_Bioirrig DOCinit v_burial_Fluid DOC_root z_sed alpha_DOC

v_burial_f = interp1(z_sed,v_burial_Fluid,x);
DOC_root_PDE = interp1(xmesh,DOC_root,x);
Alpha_Bioirrig_1 = interp1(z_sed,Alpha_Bioirrig,x);

[x1,t1]=meshgrid(x,t);
RC1 = interp2(xmesh,tspan,RC,x1,t1);

c = 1;

f = DDOC.*dudx;

S2 = - v_burial_f.* dudx + RC1 * 1E9 - alpha_DOC * RC1 * (u(1)/(u(1)+k_DOC)) * 1E9 + ...
     (Alpha_Bioirrig_1.*(DOCinit-u(1))) + DOC_root_PDE; % umol/l/year
s = S2;  %umol/l/year

end