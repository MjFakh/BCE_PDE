function u0 = pdeic3(x)
global H2CO3_IC xmesh

H2CO3_IC1 = interp1(xmesh,H2CO3_IC,x);

u0=H2CO3_IC1;

end