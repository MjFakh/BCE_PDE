function [c,f,s] = organicPDE(x,t,u,dudx)
global k_sed v_burial z_sed Bioturb

v_burial_1 = interp1(z_sed,v_burial,x);
k_sed1 = interp1(z_sed,k_sed,x);
Db = interp1(z_sed,Bioturb,x);

c = 1;
f = Db.*dudx;
RC = k_sed1.*u; %k_sed1.*u.*rho.*((1-phi)/phi)*12; % molCorg/cm3sed/yr mineralization rate

s = - v_burial_1.* dudx - RC;  %g/gDw/year

end


% % Seagrass POC release 
% POC_root_1 = 10; % POC flux release in seagrass root zone (mmol/m2/day; Eldridge & MorserMarine 2000)
% depth_rootzone = 10; % seagrass root length (cm)
% z_root = 0:0.1:depth_rootzone;
% mu_root = 5; % value for the center of rootzone in normal distribution
% sigma_root = 1; % sigma for normal distribution of flux in the rootzone
% POC_root_2 = POC_root_1 * 1E-4 * normpdf(z_root,mu_root,sigma_root);  % mmol/cm2/day
% POC_root_22 = interp1(z_root,POC_root_2,xmesh);
% poros_root = interp1(z_sed,poros,xmesh);
% POC_root = ((12*365)/1000).*(1./(rho * (1-poros_root))).*(POC_root_22./(xmesh(1,2)-xmesh(1,1))); % (gr/grDw/year)