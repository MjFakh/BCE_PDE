function u0 = pdeic1(x)
global Sulfate_IC xmesh

 Sulfate_IC1 = interp1(xmesh,Sulfate_IC,x);

 u0=Sulfate_IC1;

end