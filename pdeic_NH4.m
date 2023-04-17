function u0 = pdeic_NH4(x)
global NH4_IC xmesh

NH4_IC1 = interp1(xmesh,NH4_IC,x);
u0=NH4_IC1;

end