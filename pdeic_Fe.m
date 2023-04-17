function u0 = pdeic_Fe(x)
global Fe_IC xmesh

 Sulfate_IC1 = interp1(xmesh,Fe_IC,x);

 u0=Sulfate_IC1;

end