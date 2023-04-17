%% Sediment Diagenesis Model for Blue Carbon Ecosystems

clearvars -except M_ICs & kk & ALK_DIC_output

% ----------------------------- INPUT PARAMETERS ---------------------------

global v_burial z_sed KFe_HS Oxygen
global k_sed k_O2 DSO4 DH2S DO2 k_SO4 Kreox Iron_conc Bioturb Calcium DHCO3 HCO3init NH4init NO3init
global O2init SO4init HSinit C_organic xmesh tspan rho poros RC Alpha_Bioirrig  alpha_bioturb
global R_respi R_SRR Ksp_ca k_calcite DICinit R1_carb CO3_1 BE
global Organic_IC Oxygen_IC Sulfate_IC DIC_IC H2CO3_IC Alk_IC DOC_IC HS_IC  Fe_IC NO3_IC NH4_IC v_burial_Fluid 
global O2_root DOC_root DDOC k_DOC Calcium_activity NPP DOCinit DOC alpha_DOC CaCO3_IC F_CaCO3
global Sulfide kFeS FeooH Iron KFEMonod Feinit Sulfate R_HS_Ox CaCO3 kFeOx R1_carb_form R1_carb_disso
global CH4init CH4_IC k1_AOM k1_CH4O2 DCH4 k_AOM k_CH4O2 k_NO3 Knitrif NH4 N_C_ratio R_NH4_Ox

Bioirrig_top = 80; %1/yr
Bioirrig_bottom = 0; %1/yr
Bioirrig_scale = 0.75; %1/yr
Lbottom = 10; %cm  ''Maximum depth''
t_final = 100; %year  ''Target Year''
Bioturbtop = 2; %cm2/yr bioturbation coefficient at top
Bioturbbottom = 0; %cm2/yr bioturbation coefficient at bottom
bioturbscale = 1.5; %cm - depth scale for decrease in bioturbation
vbottom = 0.5;  %cm/year
vbottom_fluid = 0.2;  %cm/year
porosbottom = 0.6; %porosity
porostop = 0.9; 
rho = 2.73; %gramDw/cm3Dw
porosscale = 3; %cm - depth scale for decrease in porosity
ageinit = 0.1; % initial age of organic matter at the sediment water interface
age_root = 5; 
k_O2 = 2; % Oxygen half-saturation constant (uM)
k_DOC = 200; % Oxygen half-saturation constant (uM)
k_SO4 = 20; % Sulfate half-saturation constant (uM)
k_NO3 = 2; %
KFEMonod = 20; %umol/g - Monod constant for FeOOH reduction
DSO4 = 300; % cm2/yr diffusion coefficient sulfate
DH2S = 300; % cm2/yr diffusion coefficient sulfide
DO2 = 300;   % cm2/yr diffusion coefficient oxygen
DDOC = 100;   % cm2/yr diffusion coefficient oxygen
DHCO3 = 300;   % cm2/yr diffusion coefficient oxygen
DCH4 = 300;   % cm2/yr diffusion coefficient oxygen
KFe_HS = 100;  % 1/umol/l/year ''Iron Sulfide oxidation rate constant''
Iron_conc = 50; % Iron concentration (uM) 
Calcium_1 = 1E3.*[10.7 10.7 10.73 11.4 11.5 11.1 11.7 11.8 11.8]; % calcium concentration
z_calcium = [0 0.98 1.89 3.86 4.92 3.18 6.51 8.10 10];
Calcium_activity = 0.6; % calcium concentration
O2init = 220; % uM - O2 concentration at SWI
SO4init = 28000; % uM - SO4 concentration at SWI
HSinit = 0;  % uM - H2S concentration at SWI
CH4init = 0;  % uM - H2S concentration at SWI
Feinit = 0;  % uM - Fe concentration at SWI
NO3init = 10; %uM
NH4init = 0; %uM
Ksp_ca = 2000000; % uM2
k_calcite = 10; % umol/l/year
k_calcite_dis = 0.005; % yr-1
DOCinit = 1000;     % [DOC] [umol/kg]
DICinit = 2000;     % [DIC] [umol/kg]
HCO3init = 2400; %2400
Sed_rate = 1;  % sediment accumulation rate (gram/cm2/year)
BE = 0.2;    % burial efficiency of organic from the water column model
NPP = 1000;      % Net Primary Production (gram/m2/year)
Depth_WC = 200; % Water column depth
grain_size = 1; % um
K_dissolution_WC = 0.01; % mol/m2/year rate of dissolution in the water column (Hartmann)
r_a_bas = [1 2 5 10 30 50 100 200 500 1000]; % grain size matrix
alpha_bioturb = (2.375/Bioturbtop)*(0.1/vbottom);
alpha_DOC = 1;
F_CaCO3 = 1;  % Flux of CaCO3 to sediment water interface
Kreox = 500;  % 1/umol/l/year
kFeS = 0.2;
kFeOx = 100;  % 1/umol/l/year
k1_AOM = 50;  %uM
k1_CH4O2 = 2; %uM
k_AOM = 0.01; % l/year
k_CH4O2 = 0.01; % l/year
Knitrif = 1000; %1/umol/l/year
N_C_ratio = 0.15;

% Space and time resolution

Time_resolution = 100; % MaxTime
Space_resolution = 10*Lbottom; % MaxDepth
ALK = zeros(Time_resolution,Space_resolution);

% --------------- Temperature dependency parameter ------------------------

T_seawater = 25;
T_seawater_ref = 25;
Q10 = 2.5;
K_temp = (Q10.^((T_seawater-T_seawater_ref)./10));

% ---------------Calculating initial depth profiles for input parameters ----------------

n=3001;
MaxDepth = Lbottom;
z_sed=linspace(0,MaxDepth,n);
dz_sed=MaxDepth/(n-1);

poros = porosbottom + (porostop-porosbottom)*exp(-z_sed/porosscale);
Bioturb = Bioturbbottom + (Bioturbtop-Bioturbbottom)*exp(-z_sed/bioturbscale);
Alpha_Bioirrig = Bioirrig_bottom + (Bioirrig_top-Bioirrig_bottom)*exp(-z_sed/Bioirrig_scale);
v_burial = vbottom*(1-porosbottom)./(1-poros); %cm/year
v_burial_Fluid = vbottom_fluid*(1+porosbottom)./(1+poros); %cm/year

age = ageinit + cumsum(dz_sed./v_burial);

k1_org = 10.^(-0.977*log10(age) - 0.312); % more reactive
k2_org = 10.^(-0.8*log10(age) - 1.1); 
k_sed_nonbio = 10.^(-0.95*log10(age) - 0.81);    % organic C reactivity Middelburg

alpha_org = (k_sed_nonbio-k2_org)./(k1_org-k2_org);

xmesh_org = linspace(0,Lbottom,2*Lbottom);
Probability_mixing = Bioturb./Bioturbtop;
alpha_org1 = interp1(z_sed,alpha_org,xmesh_org);
Probability_mixing1 = interp1(z_sed,Probability_mixing,xmesh_org);

for i=1:size(alpha_org1,2)  
    for j=1:10*size(alpha_org1,2)
        alpha_org_11(1,j) = alpha_org1(randi(numel(alpha_org1)));  
    end 
    alpha_mean1(1,i) = mean(alpha_org_11,'all');
    alpha_org2(1,i) = Probability_mixing1(1,i).*alpha_mean1(1,i) + (1-Probability_mixing1(1,i))*alpha_org1(1,i);  
end

alpha_org_final = interp1(xmesh_org,alpha_org2,z_sed);
k_sed = alpha_org_final.*k1_org + (1-alpha_org_final).*k2_org;

% ------------------- Determining with or without BCE -------------


IC_PDEs = 0;    % Swithcing betweeen without and with blue carbon coditions. Zero value assumes no blue carbon whereas ''one'' 
                % correspond to with blue carbon and uses the initial
                % conditions solved from no blue carbon when the solution of PDEs reached steady state

if IC_PDEs == 0     
     
  Organic_IC = zeros(Space_resolution,1);
  Oxygen_IC = zeros(Space_resolution,1);
  DOC_IC = ones(Space_resolution,1);
  Sulfate_IC = ones(Space_resolution,1);
  HS_IC = ones(Space_resolution,1);
  Fe_IC = ones(Space_resolution,1);
  DIC_IC = ones(Space_resolution,1);
  H2CO3_IC = zeros(Space_resolution,1);
  Alk_IC = ones(Space_resolution,1);
  CaCO3_IC = ones(Space_resolution,1);
  CH4_IC = ones(Space_resolution,1);
  NO3_IC = zeros(Space_resolution,1);
  NH4_IC = zeros(Space_resolution,1);

  DOC_root_1 = 0; % POC flux release in seagrass root zone
  O2_root_1  = 0; % POC flux release in seagrass root zone 
  POC_root_1 = 0;  % POC flux release in seagrass root zone 
  
else
    
  Organic_IC = M_ICs(:,1)';
  Oxygen_IC = M_ICs(:,2)';
  Sulfate_IC = M_ICs(:,3)';
  DIC_IC = M_ICs(:,4)';
  H2CO3_IC = M_ICs(:,5)';
  Alk_IC = M_ICs(:,6)';
  DOC_IC = M_ICs(:,7)';
  CaCO3_IC = M_ICs(:,8)';
  HS_IC = M_ICs(:,9)';
  Fe_IC = M_ICs(:,10)';
  CH4_IC = M_ICs(:,11)';
  NH4_IC = M_ICs(:,12)';
  NO3_IC = M_ICs(:,13)';

  DOC_root_1 = 1000;    % DOC flux release in seagrass root zone (mmol/m2/day; Eldridge & MorserMarine 2000)
  O2_root_1  = 500;     % O2 flux release in seagrass root zone (mmol/m2/day; Eldridge & MorserMarine 2000)
  POC_root_1 = 0.08;    % POC flux release in seagrass root zone (mol/m2/day; Eldridge & MorserMarine 2000)

  NPP = 1000;   %   % Net Primary Production (gram/m2/year)
    
end

% ----------------------- Initial Concentrations -----------------

pH_top = carbonate(HCO3init,DICinit); % bottom water pH based on DIC and ALK top boundary
pH_initial = pH_top*ones(Time_resolution,Space_resolution);
xmesh = linspace(0,Lbottom,Space_resolution);
CO3_top = Carb_CO3(HCO3init,DICinit); % bottom water pH based on DIC and ALK top boundary
CO3_1 = CO3_top*ones(Time_resolution,Space_resolution);
CaCO3 = 1E5.*ones(Time_resolution,Space_resolution);
Sulfide = zeros(Time_resolution,Space_resolution); %intial value for sulfide
Calcium = interp1(z_calcium,Calcium_1,xmesh);

R1_carb_form = zeros(Time_resolution,Space_resolution);
R1_carb_disso = zeros(Time_resolution,Space_resolution);

% ------------------- BCE root fluxes ------------------------------------

% DOC release 

xmesh_BCE = linspace(0,Lbottom,3*Lbottom);
depth_rootzone = 10; % seagrass root length (cm)
z_root = 0:0.1:depth_rootzone;
mu_root = 4; % value for the center of rootzone in normal distribution
sigma_root = 6; % sigma for normal distribution of flux in the rootzone
DOC_root_2 = DOC_root_1 * 1E-4 * normpdf(z_root,mu_root,sigma_root);  % mmol/cm2/year
DOC_root_22 = interp1(z_root,DOC_root_2,xmesh_BCE);
poros_root = interp1(z_sed,poros,xmesh_BCE);

for i_root = 1:(3*Lbottom)

    if depth_rootzone >= xmesh_BCE(1,i_root) 

DOC_root3(1,i_root) = 1E6 * (DOC_root_22(1,i_root)./(xmesh_BCE(1,2)-xmesh_BCE(1,1))); % (umol/l/year)

    else 

DOC_root3(1,i_root) = 0;

    end

end

DOC_root = interp1(xmesh_BCE,DOC_root3,xmesh);

% O2 release 

muO2_root = 4; % value for the center of rootzone in normal distribution
sigmaO2_root = 6; % sigma for normal distribution of flux in the rootzone
O2_root_2 = O2_root_1 * 1E-4 * normpdf(z_root,muO2_root,sigmaO2_root);  % mmol/cm2/year
O2_root_22 = interp1(z_root,O2_root_2,xmesh_BCE);

for i_root = 1:(3*Lbottom)
    if depth_rootzone >= xmesh_BCE(1,i_root) 
O2_root3(1,i_root) = 1E6 * (O2_root_22(1,i_root)./(xmesh_BCE(1,2)-xmesh_BCE(1,1))); % (umol/l/year)
    else 
O2_root3(1,i_root) = 0;
    end
end

O2_root = interp1(xmesh_BCE,O2_root3,xmesh);

% POC release 

muPOC_root = 4; % value for the center of rootzone in normal distribution
sigmaPOC_root = 6; % sigma for normal distribution of flux in the rootzone

POC_root_2 = POC_root_1 * normpdf(z_root,muPOC_root,sigmaPOC_root);  % mmol/cm2/year
POC_root = interp1(z_root,POC_root_2,xmesh_BCE);

k_sed_root = 10.^(-0.95*log10(age_root) - 0.8); % more reactive
poros2 = interp1(z_sed,poros,xmesh);


for i_root = 1:(3*Lbottom)
    if depth_rootzone >= xmesh_BCE(1,i_root) 
RC_root3(1,i_root) = k_sed_root .*POC_root(1,i_root).*rho.*((1-poros_root(1,i_root))./(poros_root(1,i_root).*12)); % (umol/l/year)
    else 
RC_root3(1,i_root) = 0;
    end
end

RC_root = interp1(xmesh_BCE,RC_root3,xmesh);

% -----------------------------------------------------------------

K_converge = 1;
iteration = 1;
n_converge = 1;
iteration_tolerance = 0.2;

while abs(K_converge) > iteration_tolerance % for iteration_1=1:2 
            
% ------------- ORGANIC MATTER DEGRADATION (POC) --------------------------

hold on

m = 0;
xmesh = linspace(0,Lbottom,Space_resolution);
tspan = linspace(0,t_final,Time_resolution);


% Solving PDE

sol = pdepe(m,@organicPDE, @organicIC,@organicbc,xmesh,tspan);

C_organic = sol(:,:,1);
xmesh_oxygen = linspace(0,Lbottom,2000);


for i = 1:Time_resolution

BEsed_org(i,:) = C_organic (i,Space_resolution)./C_organic(i,1);  % Burial Efficiency of Organic

end

k_sed2 = interp1(z_sed,k_sed,xmesh);
poros2 = interp1(z_sed,poros,xmesh);

RC = K_temp.*(RC_root + (k_sed2.*C_organic.*rho.*((1-poros2)./(poros2.*12)))); % molCorg/cm3/yr mineralization rate

% ------------- ORGANIC MATTER DEGRADATION (DOC) --------------------------

hold on

m = 0;
xmesh = linspace(0,Lbottom,Space_resolution);
tspan = linspace(0,t_final,Time_resolution);


% Solving PDE

sol = pdepe(m,@DOCfun, @DOCic,@DOCbc,xmesh,tspan);

C_DOC = sol(:,:,1);
DOC = C_DOC;

% ---------------------------- OXYGEN -------------------------------------

% Solving PDE

m = 0; 
xmesh = linspace(0,Lbottom,Space_resolution);
tspan = linspace(0,t_final,Time_resolution);
sol = pdepe(m,@pdefun,@pdeic,@pdebc,xmesh,tspan);

%  Calculating concentrations

C_O2 = sol(:,:,1);

for i = 1:Time_resolution
    
   for  j=1:Space_resolution

     if C_O2(i,j) < 1

     C_O2(i,j) = 0;
    
     end

   end 
end

Oxygen = C_O2;
R_respi = alpha_DOC.* RC.* (DOC./(DOC+k_DOC)).* (Oxygen./(Oxygen+k_O2)).*1E9; %rate of aerobic respiration umol/l/year

% ---------------------------- AMMONIUM -----------------------------------

% Solving PDE

m = 0; 
xmesh = linspace(0,Lbottom,Space_resolution);
tspan = linspace(0,t_final,Time_resolution);
sol = pdepe(m,@pde_NH4,@pdeic_NH4,@pdebc_NH4,xmesh,tspan);
C_NH4 = sol(:,:,1);

for i = 1:Time_resolution
    
   for  j=1:Space_resolution

     if C_NH4(i,j) < 0

     C_NH4(i,j) = 0;
    
     end

   end 
end

NH4 = C_NH4;
R_NH4_Ox = NH4.* Oxygen.*Knitrif; %rate of sulfide reduction umol/l/year

% ---------------------------- NITROGEN -----------------------------------

% Solving PDE

m = 0; 
xmesh = linspace(0,Lbottom,Space_resolution);
tspan = linspace(0,t_final,Time_resolution);
sol = pdepe(m,@pde_NO3,@pdeic_NO3,@pdebc_NO3,xmesh,tspan);
C_NO3 = sol(:,:,1);

for i = 1:Time_resolution
    
   for  j=1:Space_resolution

     if C_NO3(i,j) < 0

     C_NO3(i,j) = 0;
    
     end

   end 
end

NO3 = C_NO3;

% ---------------------------- OXYGEN PENTRATION DEPTH --------------------

Oxygen_1 = interp1(xmesh,Oxygen(Time_resolution,:),z_sed);
count_OPD = 0;

for i=1:n
    if Oxygen_1(1,i) < 1

        count_OPD = count_OPD + 1;
        OPD_1(1,count_OPD) = z_sed(1,i);
        num_OPD1(1,count_OPD) = i;

    end
end

if min(Oxygen_1) > 1

    OPD_1 = Lbottom;
    num_OPD1 = n;
end

OPD = min(OPD_1);
num_OPD = min(num_OPD1);

mm_count = 0;
for i=1:n
    Ironoxy(1,i) = 0.1*((1/(1+exp(z_sed(1,i)-OPD)))+2*exp(-((z_sed(1,i)-OPD)^2)/2));
    mm_count=mm_count+1;
    FeOx(1,mm_count)=Ironoxy(1,i);
end
FeooH =  FeOx;

% -------------------------- IRON ----------------------------------------

m = 0; 
xmesh = linspace(0,Lbottom,Space_resolution);
tspan = linspace(0,t_final,Time_resolution);
sol = pdepe(m,@pde_Fe,@pdeic_Fe,@pdebc_Fe,xmesh,tspan);

%  Calculating concentration

C_Fe = sol(:,:,1);
Iron = C_Fe;

% ------------------------ SULFATE ---------------------------------------

% Solving PDE

m = 0; 
xmesh = linspace(0,Lbottom,Space_resolution);
tspan = linspace(0,t_final,Time_resolution);
sol = pdepe(m,@pdefun1,@pdeic1,@pdebc1,xmesh,tspan);

%  Calculating concentrations and rate

C_SO4 = sol(:,:,1);
Sulfate = C_SO4;
Inhib = (k_O2./(Oxygen+k_O2)); % inhibition term for sulfate reduction by oxic respiration
R_SRR = alpha_DOC.* RC.* (DOC./(DOC+k_DOC)).*Inhib.* (Sulfate./(Sulfate+k_SO4)).*1E9; %rate of sulfate reduction umol/l/year

% ------------------------- SULFIDE ---------------------------------------

m = 0; 
xmesh = linspace(0,Lbottom,Space_resolution);
tspan = linspace(0,t_final,Time_resolution);
sol = pdepe(m,@pde_HS,@pdeic_HS,@pdebc_HS,xmesh,tspan);

%  Calculating concentration

C_HS = sol(:,:,1);

for i = 1:Time_resolution
    
   for  j=1:Space_resolution

     if C_HS(i,j) < 0.5

     C_HS(i,j) = 0;
    
     end

   end 
end

Sulfide = C_HS;

R_HS_Ox = Sulfide.* Oxygen.*Kreox; %rate of sulfide reduction umol/l/year
R_HS_Fe = Sulfide.* Iron.*kFeS; %rate of sulfide reduction umol/l/year

% ------------------------------- CaCO3 -----------------------------------

sigma_carb = (Calcium_activity.*Calcium.* CO3_1)./Ksp_ca - 1;  


for i=1:Time_resolution
    
    for j=1:Space_resolution
    
    if sigma_carb(i,j) > 0
         n_power_CaCO3 = 1.76;
         K_sigma = abs(sigma_carb(i,j)).^n_power_CaCO3;
         R1_carb_form(i,j) = K_sigma.* k_calcite; %umol/l/year
    else 
         R1_carb_form(i,j) = 0;
    end


    if sigma_carb(i,j) < 0

        if -0.2 < sigma_carb(i,j) && sigma_carb(i,j) < 0

        n_power_CaCO3 = 0.11;
        k_calcite_dis = 0.005;
        K_sigma = abs(sigma_carb(i,j)).^n_power_CaCO3;
        R1_carb_disso(i,j) = - K_sigma.* k_calcite_dis.*CaCO3(i,j); %umol/l/year 

        else 

        n_power_CaCO3 = 4;
        k_calcite_dis = 0.1;
        K_sigma = abs(sigma_carb(i,j)).^n_power_CaCO3;
        R1_carb_disso(i,j) = - K_sigma.* k_calcite_dis.*CaCO3(i,j); %umol/l/year 

        end

 
    else 
        R1_carb_disso(i,j) = 0;
    end

    end

end


m = 0;

R1_carb  = R1_carb_form + R1_carb_disso;

xmesh = linspace(0,Lbottom,Space_resolution);
tspan = linspace(0,t_final,Time_resolution);

% Solving PDE

sol = pdepe(m,@CaCO3PDE, @CaCO3IC,@CaCO3bc,xmesh,tspan);

C_CaCO3 = sol(:,:,1); % umol/l

for i = 1:Time_resolution
    
   for  j=1:Space_resolution

     if C_CaCO3(i,j) < 0

     C_CaCO3(i,j) = 0;
    
     end

   end 
end

CaCO3 = C_CaCO3;

% ------------------------- METHANE ---------------------------------------

m = 0; 
xmesh = linspace(0,Lbottom,Space_resolution);
tspan = linspace(0,t_final,Time_resolution);
sol = pdepe(m,@pde_CH4,@pdeic_CH4,@pdebc_CH4,xmesh,tspan);

%  Calculating concentration

C_CH4 = sol(:,:,1);

Methane = C_CH4;

      for i = 1:Time_resolution
    
          F_diff_CH4 (1,i) = DCH4.*((C_CH4(i,2) - C_CH4(i,1))./(xmesh(1,2)-xmesh(1,1)))*1E-3;
    
      end

R_CH4 = alpha_DOC.* RC.* (DOC./(DOC+k_DOC)).*Inhib.* (k_SO4./(Sulfate+k_SO4)).*1E9;
R_OM =  (k_AOM.*C_CH4.*(Sulfate./(Sulfate+k1_AOM)));%
R_AOM  = (k_CH4O2.*C_CH4.*(Oxygen./(Oxygen+k1_CH4O2)));

% ------------------------------- DIC -------------------------------------

m = 0;
xmesh = linspace(0,Lbottom,Space_resolution);
tspan = linspace(0,t_final,Time_resolution);
sol = pdepe(m,@pdefun4, @pdeic4,@pdebc4,xmesh,tspan);

C_DIC = sol(:,:,1);


       for i = 1:Time_resolution
    
       
          F_diff_DIC (1,i) = DHCO3.*((C_DIC(i,2) - C_DIC(i,1))./(xmesh(1,2)-xmesh(1,1)))*1E-3;
    
    
       end

% ------------------------------ ALKALINITY ------------------------------

m = 0;
xmesh = linspace(0,Lbottom,Space_resolution);
tspan = linspace(0,t_final,Time_resolution);
sol = pdepe(m,@pdefun2, @pdeic2,@pdebc2,xmesh,tspan);

C_alka = sol(:,:,1);

ALK = C_alka;

      for i = 1:Time_resolution
    
          F_diff (1,i) = DHCO3.*((C_alka(i,2) - C_alka(i,1))./(xmesh(1,2)-xmesh(1,1)))*1E-3;
    
      end
     
      
% ------------------------------- pH and H2CO3 (2 for 6 calculation) --------------------------------------

for i=1:Time_resolution
    
    for j=1:Space_resolution

           pH_1(i,j) = carbonate(ALK(i,j),C_DIC(i,j));
           C_H2CO3(i,j) = carb_acid(ALK(i,j),C_DIC(i,j));
           CO3_1(i,j) = Carb_CO3(ALK(i,j),C_DIC(i,j));

    end

end

pH = pH_1;

% -------------------------- Convergence coefficient ----------------------

pH_iteration(iteration,:) = pH(Time_resolution,:);
alpha_converge(1,iteration) = F_diff(1,Time_resolution);
alpha_converge_1(1,n_converge) = alpha_converge(1,iteration);


CO3_loop(1,iteration) = CO3_1(end);


     if n_converge > 1 %&& alpha_converge(1,iteration) > 0
        
       K_converge = (alpha_converge_1(1,n_converge)-alpha_converge_1(1,n_converge-1))/alpha_converge_1(1,n_converge-1);
  
     end


     iteration = iteration + 1;
     n_converge = n_converge + 1;

     
end  % iteration ends here 


% Storing steady-state solutions as initial conditions for PDEs
M_ICs = [C_organic(Time_resolution,:)' C_O2(Time_resolution,:)' C_SO4(Time_resolution,:)' C_DIC(Time_resolution,:)' ...
         C_H2CO3(Time_resolution,:)' ALK(Time_resolution,:)' C_DOC(Time_resolution,:)' CaCO3(Time_resolution,:)'...
         Sulfide(Time_resolution,:)' C_Fe(Time_resolution,:)' C_CH4(Time_resolution,:)' NH4(Time_resolution,:)' NO3(Time_resolution,:)'];
     
ALK_DIC_R = F_diff./F_diff_DIC;


% -------------------------------------------------------------------------
% ----------------------------- PLOTS -------------------------------------

n_plot = 7; % number of plots in each row
m_plot = 2; % number of total rows


% Organic

subplot(m_plot,n_plot,1);

title('Organic (%gDw)')

for i = 1:Time_resolution
    
    C_organic11(i,:)=interp1(xmesh,C_organic(i,:),xmesh_oxygen);
        
    hold on

end

POC_root_plot = interp1(xmesh_BCE,POC_root,xmesh_oxygen);
plot((C_organic11(Time_resolution,:) + POC_root_plot).*100,xmesh_oxygen,'lineWidth',2); axis ij
xlim([0 6])
box on


% Oxygen

subplot(m_plot,n_plot,2);
title('Oxygen (\muM)')


for i = 1:Time_resolution
    
    C_O21(i,:)=interp1(xmesh,C_O2(i,:),xmesh_oxygen);
     
    hold on

end


plot(C_O21(Time_resolution,:),xmesh_oxygen,'lineWidth',2); axis ij
box on


% Sulfate

subplot(m_plot,n_plot,3);
title('Sulfate (\muM)')


for i = 1:Time_resolution
    
    C_SO41(i,:)=interp1(xmesh,C_SO4(i,:),xmesh_oxygen);
        
    hold on

end

plot(C_SO41(Time_resolution,:),xmesh_oxygen,'lineWidth',2); axis ij
box on

% Sulfide

subplot(m_plot,n_plot,4);
title('Sulfide (\muM)')


for i = 1:Time_resolution
    
    C_H2S1(i,:)=interp1(xmesh,Sulfide(i,:),xmesh_oxygen);
        
    hold on

end

plot(C_H2S1(Time_resolution,:),xmesh_oxygen,'lineWidth',2); axis ij
box on

% Iron

subplot(m_plot,n_plot,5);
title('Fe^{2+} (\muM)')


for i = 1:Time_resolution
    
    C_Fe1(i,:)=interp1(xmesh,Iron(i,:),xmesh_oxygen);
        
    hold on

end

plot(C_Fe1(Time_resolution,:),xmesh_oxygen,'lineWidth',2); axis ij
box on

% DIC

subplot(m_plot,n_plot,6);
title('DIC (\muM)')


for i = 1:Time_resolution
    
    
    C_DIC1(i,:)=interp1(xmesh,C_DIC(i,:),xmesh_oxygen);
    
    F_diff_DIC (1,i) = DHCO3.*((C_DIC(i,2) - C_DIC(i,1))./(xmesh(1,2)-xmesh(1,1)))*1E-3;
    
    hold on

    
end


plot(C_DIC1(Time_resolution,:),xmesh_oxygen,'lineWidth',2); axis ij
box on


% ALK

subplot(m_plot,n_plot,7);
title('ALK (\muM)')


for i = 1:Time_resolution
    
    C_alka1(i,:)=interp1(xmesh,C_alka(i,:),xmesh_oxygen);
    F_diff (1,i) = DHCO3.*((C_alka(i,2) - C_alka(i,1))./(xmesh(1,2)-xmesh(1,1)))*1E-3;
    F_diff_bottom (1,i) = DHCO3.*((C_alka(i,size(xmesh,2)) - C_alka(i,size(xmesh,2)-1))./(xmesh(1,2)-xmesh(1,1)))*1E-3;
        
    hold on

end

plot(C_alka1(Time_resolution,:),xmesh_oxygen,'lineWidth',2); axis ij
box on


% Carbonic Acid

subplot(m_plot,n_plot,8);
title('Carb Acid (\muM)')


for i = 1:Time_resolution
    
    
    C_H2CO31(i,:)=interp1(xmesh,C_H2CO3(i,:),xmesh_oxygen);
        
    hold on

    
end

plot(C_H2CO31(Time_resolution,:),xmesh_oxygen,'lineWidth',2); axis ij
box on


% pH
subplot(m_plot,n_plot,9);
title('pH')


for i = 1:Time_resolution
    
    
    pH1(i,:)=interp1(xmesh,pH(i,:),xmesh_oxygen);
        
    hold on

    
end

plot(pH1(Time_resolution,:),xmesh_oxygen,'lineWidth',2); axis ij
box on


dz_integ = xmesh(1,2)-xmesh(1,1);
R_SRR_integ = cumsum(R_SRR(Time_resolution,:).*dz_integ.*1E-3);
R1_carb_integ = cumsum(R1_carb(Time_resolution,:).*dz_integ.*1E-3);


% Organic degradation rate

subplot(m_plot,n_plot,10);


for i = 1:Time_resolution
    
    RC1(i,:) = interp1(xmesh,RC(i,:),xmesh_oxygen);
            
end

plot(RC1(Time_resolution,:).* 1E9,xmesh_oxygen,'lineWidth',2); axis ij  %umol/l/year
title('Mineralization Rate (\mumol/l/year)')
box on


% DIC rate

subplot(m_plot,n_plot,11);


R_DIC = RC.* 1E9 - R1_carb; %umol/l/year

for i = 1:Time_resolution
    
    R_DIC1(i,:) = interp1(xmesh,R_DIC(i,:),xmesh_oxygen);
        
    
end

plot(R_DIC1(Time_resolution,:),xmesh_oxygen,'lineWidth',2); axis ij  %umol/l/year
title('DIC rate (\mumol/l/year)')
box on


% ALK rate

subplot(m_plot,n_plot,12);


R_ALK = R_SRR - R1_carb; %umol/l/year


for i = 1:Time_resolution
    
    R_ALK1(i,:) = interp1(xmesh,R_ALK(i,:),xmesh_oxygen);
        
    
end


for i = 1:Time_resolution
    
    R_ALK1_integ(i,:) = (cumsum(R_ALK(i,:).*dz_integ)).*1E-3; %umol/cm2/year
    R_ALK1_integ1(i,:) = interp1(xmesh,R_ALK1_integ(i,:),xmesh_oxygen);
    R_ALK1_integ1_final(1,i) = R_ALK1_integ1(end);
       
end

plot(R_ALK1(Time_resolution,:),xmesh_oxygen,'lineWidth',2); axis ij  %umol/l/year
title('ALK rate (\mumol/l/year)')
box on


% Aerobic respiration and sulfate reduction rates

subplot(m_plot,n_plot,13);
plot(R_SRR(Time_resolution,:),xmesh,R_respi(Time_resolution,:),xmesh,'lineWidth',2); axis ij %umol/l/year
title('Rate (\mumol/l/year)')
legend('Sulfate Red','Aerobic Resp');


