function [pl,ql,pr,qr] = pdebc_NH4(xl,ul,xr,ur,t)
global NH4init

pl = ul(1)-NH4init;
ql = 0;
pr = 0;
qr = 1;

end