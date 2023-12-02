%%   Dynamical analysis for atmospheric oxygen

% -------------------------------------------------------------------------------------------------------------- 
% -------------------------------------------------------------------------------------------------------------- 
% -------------------------------------------------------------------------------------------------------------- 
% ------------------------------------------------ NOTES ------------------------------------------------------- 
% -------------------------------------------------------------------------------------------------------------- 
% -------------------------------------------------------------------------------------------------------------- 
% -------------------------------------------------------------------------------------------------------------- 
% --------------------- A Stochastic dynamic analysis of atmospheric oxygen level ------------------------------
% --------- The code calculates the range of steady state solutions for the atmospheric oxygen mass balance ---- 
% -------------------------------------------------------------------------------------------------------------- 
% -------------------------------------------------------------------------------------------------------------- 
% -------------------------------------------------------------------------------------------------------------- 
% -------------------------------------------------------------------------------------------------------------- 
% -------------------------------------------------------------------------------------------------------------- 

clear all
alpha_land = 1;
count_pCH4 = 1;

% -------  SET PARAMETERS -----------

for i=1:1401

power_pO2 = -13:0.01:1;
pO2_modern = 210000;                                     % modern atmospheric oxygen level (ppm)
pO2 = pO2_modern.*10.^(power_pO2);

Weathering_pyrite_MC = [1500 1700 2000 2200 2500];       % Range of pyrite weathering  
Weathering_pyrite  = (Weathering_pyrite_MC(randi(numel(Weathering_pyrite_MC)))).*1E4;

% -------  Parameters for Net Primary Production on Land -----------

pCO2_half = 183.6;
pCO2_min = 10;
pCO2_1_MC = [100 200 280 300 500 600 700 1000];  % range of pCO2 during Phanerozoic (ppm)
pCO2_1  = (pCO2_1_MC(randi(numel(pCO2_1_MC)))); 

Temp_MC = [15 16 17 18 19 20 25];
Temp  = (Temp_MC(randi(numel(Temp_MC)))); % temperature range during Phanerozoic (oC)

alpha_CO2 = (pCO2_1-pCO2_min)./(pCO2_half + (pCO2_1-pCO2_min));
alpha_temperature = 1 - ((Temp - 25)./25).^2;
alpha_oxygen = 1.5 - 0.5.*((pO2./pO2_modern));

for ff = 1:1401

if alpha_oxygen(1,ff) < 0.0001

    alpha_oxygen(1,ff) = 0.0001;

end

end

Fire_feedback = max(586.2.*(pO2.*1E-6) - 122.102,0);

for ff = 1:1401

if Fire_feedback(1,ff) > 24.45

    Fire_feedback(1,ff) = 24.45;

end

end


alpha_NPP_land = interp1([-13 -5 -3 -2 -1 -0.5 -0.8 -0.9 0 1],[0 0 0 0 0 0 0 0 1 1],power_pO2);            

NPP_land_MC = [1 2 5 7];                                 
NPP_land_1  = (NPP_land_MC(randi(numel(NPP_land_MC))));

k_fire = 20;                % Fire feedback coefficient

alpha_fire = k_fire./(k_fire - 1 + Fire_feedback);

alpha_plant_feedback = alpha_temperature.*alpha_oxygen.*alpha_fire.*alpha_CO2;

NPP_land  = 2.*alpha_NPP_land.*NPP_land_1.*alpha_plant_feedback.*1E4.*1E3; % Calculation of land Net Primary Production
Pyrite_weathering_threshold_MC = [5E-5 7E-4 1E-3 5E-3];
Pyrite_weathering_threshold   = (Pyrite_weathering_threshold_MC(randi(numel(Pyrite_weathering_threshold_MC)))).*pO2_modern;       % atm

NPP_MC = [10 20 30 40 50 60 70 80 90 100];   % Range of modern marine net primary production (grC/m2/yr)
NPP = (NPP_MC(randi(numel(NPP_MC)))).* (1/12).*1E-9.* 1E4;     

AREAtot = 3.8*1E14;

% ------ burial efficiency and fraction of anoxia -----------

% ---- fraction of anoxia ---------

f_anoxia_MC_1 = 1;
f_anoxia_MC_2 = [0.85 0.9 0.95 0.99];
f_anoxia_MC_3 = [0.4 0.5 0.6];
f_anoxia_MC_4 = [0.2 0.3 0.4];
f_anoxia_MC_5 = [0.1 0.15 0.2 0.25];
f_anoxia_MC_6 = [0.01 0.03 0.05 0.07 0.1];
f_anoxia_MC_7 = [0.001 0.003 0.005];
f_anoxia_MC_8 = [0.001 0.003 0.005];

f_anoxia  = interp1([-13 -5 -3 -2 -1 -0.5 0 1],[f_anoxia_MC_1(randi(numel(f_anoxia_MC_1))) f_anoxia_MC_2(randi(numel(f_anoxia_MC_2))) f_anoxia_MC_3(randi(numel(f_anoxia_MC_3))) ...
                                           f_anoxia_MC_4(randi(numel(f_anoxia_MC_4))) f_anoxia_MC_5(randi(numel(f_anoxia_MC_5))) f_anoxia_MC_6(randi(numel(f_anoxia_MC_6))) ...
                                           f_anoxia_MC_7(randi(numel(f_anoxia_MC_7))) f_anoxia_MC_8(randi(numel(f_anoxia_MC_8)))],power_pO2); %0.01;% * (Oxygen_modern./Oxy1conc(W));

% ---- burial efficiency ---------

BE_1 = [0.6 0.7 0.8];
BE_2 = [0.4 0.5 0.6];
BE_3 = [0.2 0.3 0.4];
BE_4 = [0.05 0.1 0.2];
BE_5 = [0.03 0.05 0.1];
BE_6 = [0.01 0.03 0.05];
BE_7 = [0.001 0.003 0.005 0.007];
BE_8 = [0.001 0.003 0.005 0.007];
BE_modern  = BE_7(randi(numel(BE_7)));

BE  = interp1([-13 -5 -3 -2 -1 -0.5 0 1],[BE_1(randi(numel(BE_1))) BE_2(randi(numel(BE_2))) BE_3(randi(numel(BE_3))) ...
                                               BE_4(randi(numel(BE_4))) BE_5(randi(numel(BE_5))) BE_6(randi(numel(BE_6))) ...
                                               BE_modern BE_modern],power_pO2); 


% -------- productivity factor -------------

fprod_1 = 0;
fprod_2 = [0.0003 0.001 0.002];
fprod_3 = [0.006 0.007 0.01 0.02];
fprod_4 = [0.01 0.02 0.03];
fprod_5 = [0.03 0.05 0.07];
fprod_6 = [0.05 0.1 0.2 0.3 0.5];
fprod_7 = 1;
fprod_8 = 1;

fprod_MC_1 = fprod_1(randi(numel(fprod_1))).*BE_1(randi(numel(BE_1)));
fprod_MC_2 = fprod_2(randi(numel(fprod_2))).*BE_2(randi(numel(BE_2)));
fprod_MC_3 = fprod_3(randi(numel(fprod_3))).*BE_3(randi(numel(BE_3)));
fprod_MC_4 = fprod_4(randi(numel(fprod_4))).*BE_4(randi(numel(BE_4)));
fprod_MC_5 = fprod_5(randi(numel(fprod_5))).*BE_5(randi(numel(BE_5)));
fprod_MC_6 = fprod_6(randi(numel(fprod_6))).*BE_6(randi(numel(BE_6)));
fprod_MC_7 = fprod_7(randi(numel(fprod_7))).*BE_modern;
fprod_MC_8 = fprod_8(randi(numel(fprod_8))).*BE_modern;

fprod  = interp1([-13 -5 -3 -2 -1 -0.5 0 1],[fprod_MC_1(randi(numel(fprod_MC_1))) fprod_MC_2(randi(numel(fprod_MC_2))) fprod_MC_3(randi(numel(fprod_MC_3))) ...
                                           fprod_MC_4(randi(numel(fprod_MC_4))) fprod_MC_5(randi(numel(fprod_MC_5))) fprod_MC_6(randi(numel(fprod_MC_6))) ...
                                           fprod_MC_7(randi(numel(fprod_MC_7))) fprod_MC_8(randi(numel(fprod_MC_8)))],power_pO2); %0.01;% * (Oxygen_modern./Oxy1conc(W));


NPP2 = NPP.* (BE_modern./BE);
NPP1 = NPP2;
F_anoxic_factor = [0.01 0.03 0.05 0.07 0.1];
F_anoxic = F_anoxic_factor(randi(numel(F_anoxic_factor))).*(1-BE).* NPP1.* f_anoxia.*AREAtot;   % rate of marine oxygen sink

% --------- Marine Organic burial and Methane Production ------------

alpha_CH4_1 = [0.05 0.07 0.1 0.2 0.5];
alpha_CH4   = alpha_CH4_1(randi(numel(alpha_CH4_1)));

F_NPP   = NPP.*fprod.*AREAtot;                                % flux of organic burial
F_NPP_CH4   = (1-BE).*f_anoxia.* alpha_CH4.* NPP1.* AREAtot;  % biogenic methane flux

% ----- Calculating Rates of Methane Oxidation -------

mu = 1.773*10^20;       % number of moles in 1 atmosphere

a_1 = 0.002998345;      % parameters for the kinetics of atmospheric methane oxidation
a_2 = -0.165030964;
a_3 = 3.221048922;
a_4 = -25.757487116;
a_5 = 70.985147970;

si = log10(mu.*pO2.*1E-6);

psi = 10.^(a_1.*si.^4 + a_2.*si.^3 + a_3.*si.^2 + a_4.*si.^1 + a_5);   % constant value for rate of atmospheric methane oxidation (Goldblatt et al. 2006)

% --------------------- pCH4 -------------------------

alpha_pCH4_1 = [1E-19 1E-18];
alpha_pCH4   = alpha_pCH4_1(randi(numel(alpha_pCH4_1)));
k_pCH4_1 = [1E-2 1E-1 1];
k_pCH4 = k_pCH4_1(randi(numel(k_pCH4_1))).*1E-6;
max_pCH4_1 = [100 300 500 800 1000];
max_pCH4 = max_pCH4_1(randi(numel(max_pCH4_1))); 
R_H2_MC = [1 2 3 4 5 6 7 8 9 10];
pCH4 = (max_pCH4.*(k_pCH4./(pO2./210000 + k_pCH4))).*mu.*1E-6;

% ----------------------------------------------------

R_CH4O2 = (psi.*mu.*pO2.*((pCH4).^0.7)).*1E4.*1E-9.*1E-6.*alpha_pCH4;    %  rate of methane oxidation
R_H2 = R_H2_MC(randi(numel(R_H2_MC))).*1E-5.*pCH4.*1E4.*1E-9;            %  rate of hydrogen escape

% ------------ Pyrite and organic weathering ----------

FW_pyr = (Weathering_pyrite).* (pO2./(pO2 + Pyrite_weathering_threshold)); % Oxidative weathering, pyrite
FW_oc_MC = [1 1.5 2 3 3.5 4 4.5 5 7 10];
FW_oc =  (FW_oc_MC(randi(numel(FW_oc_MC)))).*1E12.*1E-9.*1E4.* (pO2./(pO2 + Pyrite_weathering_threshold)).* ...
         ((NPP1.*fprod)/(NPP.*BE_modern));    % Oxidative weathering, organic 

% ------------- metamorphism/volcanic reductant input -----------

F_volc_MC_1 = [1 2 3 4 5];
F_volc_MC_2 = [1 2 3 4 5];
F_volc_MC_3 = [1 2 3 4 5];
F_volc_MC_4 = [0.5 1 2];
F_volc_MC_5 = [0.5 1 2];
F_volc_MC_6 = [0.5 1 2];
F_volc_MC_7 = [0.01 0.05 0.1];
F_volc_MC_8 = [0.01 0.05 0.1];
F_volc_MC_9 = [0.01 0.05 0.1];

F_redvolc_1    = interp1([-13 -10 -7 -5 -3 -1 -0.5 0 1],[F_volc_MC_1(randi(numel(F_volc_MC_1))) F_volc_MC_2(randi(numel(F_volc_MC_2))) ...
                                                   F_volc_MC_3(randi(numel(F_volc_MC_3))) F_volc_MC_4(randi(numel(F_volc_MC_4))) ...
                                                   F_volc_MC_5(randi(numel(F_volc_MC_5))) F_volc_MC_6(randi(numel(F_volc_MC_6))) ...
                                                   F_volc_MC_7(randi(numel(F_volc_MC_7))) F_volc_MC_8(randi(numel(F_volc_MC_8))) ...
                                                   F_volc_MC_9(randi(numel(F_volc_MC_9)))],power_pO2);

F_redvolc = F_redvolc_1.*1E12.*1E-9.*1E4;

% ------------- pO2 mass balance ------------------

dpO2dt = F_NPP - F_anoxic + R_H2 - R_CH4O2 - F_redvolc - (15/8).* FW_pyr - FW_oc + (alpha_land.*NPP_land);

dpO2dt_MC(i,:)    = dpO2dt;
F_NPP_MC(i,:)     = F_NPP;
F_anoxic_MC(i,:)  = F_anoxic;
R_CH4O2_MC(i,:)   = R_CH4O2;
F_H2_MC(i,:)      = R_H2;
FW_pyr_MC(i,:)    = (15/8).* FW_pyr;
FW_oc_MC1(i,:)    = FW_oc;
NPP_land_MC1(i,:) = alpha_land.*NPP_land;
pCH4_MC(i,:)      = pCH4;
F_NPP_CH4_MC(i,:) = F_NPP_CH4;

end


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% ------------------------- END OF CALCULATION ----------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% ------------------------- PlOTTING THE RESULTS --------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

color_1 = [0.04,0.33,0.52];

% pO2     

subplot(2,1,1);

hold on

fanChart(pO2./210000, dpO2dt_MC'.*1E-4.*1E-3,'median', 0:50:50, ...
    'alpha', 0, 'colormap', {'shadesOfColor', [0,0,0]});

h = findobj(gca,'Type','line');
medianline_1=get(h,'Ydata');

medianline = movmean(medianline_1,90);

hold off

plot(pO2./210000,medianline,'lineWidth',2,'Color',color_1);

set(gca, 'FontName', 'Arial','FontSize',14,'lineWidth',2);
set(gca,'xscale','log')
xlabel ('{\itp}O_2 (PAL)' );
ylabel ('d({\itp}O_2)/dt (Tmol/year)');
xlim([1E-10 2.5]);
ylim([-4 4]);
xticks([1E-12 1E-10 1E-8 1E-6 1E-4 1E-2 1E0])
xticklabels({'10^{-12}','10^{-10}','10^{-8}','10^{-6}','10^{-4}','10^{-2}','10^{0}'});

hold on
plot(pO2./210000,zeros(1,size(pO2,2)),'--k','lineWidth',2);

box on
grid off
ax.LineWidth = 2;

