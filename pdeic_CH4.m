function u0 = pdeic_CH4(x)
global CH4_IC xmesh

 CH4_IC1 = interp1(xmesh,CH4_IC,x);

 u0=CH4_IC1;

end