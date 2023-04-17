function u0 = organicIC(x)
global Organic_IC xmesh

Organic_IC1 = interp1(xmesh,Organic_IC,x);


u0 = Organic_IC1;

end
