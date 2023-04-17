function [pl,ql,pr,qr] = pdebc_HS(xl,ul,xr,ur,t)
global  HSinit

pl = ul(1)-HSinit;
ql = 0;
pr = 0;
qr = 1;

end