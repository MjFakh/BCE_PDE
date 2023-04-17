function [pl,ql,pr,qr] = pdebc_HS(xl,ul,xr,ur,t)
global  Feinit

pl = ul(1)-Feinit;
ql = 0;
pr = 0;
qr = 1;

end