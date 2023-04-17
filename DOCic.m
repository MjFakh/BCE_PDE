% function u0 = pdeic(x)
% 
% u0=[0; 0; 0];
% 
% end


function u0 = DOCic(x)
global DOC_IC xmesh

DOC_IC1 = interp1(xmesh,DOC_IC,x);
u0=DOC_IC1;

end