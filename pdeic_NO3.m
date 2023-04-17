function u0 = pdeic_NO3(x)
global NO3_IC xmesh

NO3_IC1 = interp1(xmesh,NO3_IC,x);
u0=NO3_IC1;

end