function u0 = CaCO3IC(x)
global CaCO3_IC xmesh

CaCO3_IC1 = interp1(xmesh,CaCO3_IC,x);


u0 = CaCO3_IC1;

end