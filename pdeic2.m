function u0 = pdeic2(x)
global Alk_IC xmesh

Alk_IC1 = interp1(xmesh,Alk_IC,x);

u0=Alk_IC1;

end