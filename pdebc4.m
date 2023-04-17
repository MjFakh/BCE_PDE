function [pl,ql,pr,qr] = pdebc4(xl,ul,xr,ur,t)
global DICinit

pl = ul(1)-DICinit;
ql = 0;
pr = 0;
qr = 1;

end