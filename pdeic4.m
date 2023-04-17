function u0 = pdeic4(x)
global DIC_IC xmesh

DIC_IC1 = interp1(xmesh,DIC_IC,x);

u0=DIC_IC1;

end