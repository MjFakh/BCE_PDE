function [pl,ql,pr,qr] = pdebc3(xl,ul,xr,ur,t)
global H2CO3init

pl = ul(1)-H2CO3init;
ql = 0;
pr = 0;
qr = 1;

end