function u0 = pdeic_HS(x)
global HS_IC xmesh

 Sulfate_IC1 = interp1(xmesh,HS_IC,x);

 u0=Sulfate_IC1;

end