function [pl,ql,pr,qr] = pdebc_CH4(xl,ul,xr,ur,t)
global  CH4init

pl = ul(1)-CH4init;
ql = 0;
pr = 0;
qr = 1;

end