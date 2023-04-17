% function [pl,ql,pr,qr] = pdebc(xl,ul,xr,ur,t)
% global O2init SO4init HSinit
% 
% pl = [ul(1)-O2init; ul(2)-SO4init; ul(3)-HSinit];
% ql = [0; 0; 0];
% pr = [0; 0; 0];
% qr = [1; 1; 1];
% 
% end


function [pl,ql,pr,qr] = DOCbc(xl,ul,xr,ur,t)
global DOCinit

pl = ul(1)-DOCinit;
ql = 0;
pr = 0;
qr = 1;

end

