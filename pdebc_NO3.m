function [pl,ql,pr,qr] = pdebc_NO3(xl,ul,xr,ur,t)
global NO3init

pl = ul(1)-NO3init;
ql = 0;
pr = 0;
qr = 1;

end