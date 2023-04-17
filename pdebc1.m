function [pl,ql,pr,qr] = pdebc1(xl,ul,xr,ur,t)
global  SO4init

pl = ul(1)-SO4init;
ql = 0;
pr = 0;
qr = 1;

end