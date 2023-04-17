function [pl,ql,pr,qr] = pdebc2(xl,ul,xr,ur,t)
global HCO3init

pl = ul(1)-HCO3init;
ql = 0;
pr = 0;
qr = 1;

end