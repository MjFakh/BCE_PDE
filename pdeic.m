% function u0 = pdeic(x)
% 
% u0=[0; 0; 0];
% 
% end


function u0 = pdeic(x)
global Oxygen_IC xmesh

Oxygen_IC1 = interp1(xmesh,Oxygen_IC,x);
u0=Oxygen_IC1;

end