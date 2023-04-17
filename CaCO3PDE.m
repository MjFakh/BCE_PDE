function [c,f,s] = CaCO3PDE(x,t,u,dudx)
global v_burial z_sed Bioturb R1_carb xmesh tspan

[x1,t1]=meshgrid(x,t);
R1_carb1 = interp2(xmesh,tspan,R1_carb,x1,t1); %umol/l/year
v_burial_1 = interp1(z_sed,v_burial,x);
Db = interp1(z_sed,Bioturb,x);

c = 1;
f = Db.*dudx;

s = - v_burial_1.* dudx + R1_carb1;  %umol/l/year

end