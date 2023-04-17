% -------------------------------------------------------------------------
%
% THIS IS A SIMPLE SCRIPT FOR SOLVING THE CARBONATE SYSTEM, BASED ON THE
% ROUTINES IN ZEEBE + WOLF-GLADROW [2001].
%
% Lines 15-24 define physical and chemical constants. Note that we have to
%    specify a total boron concentration (B_T)
%
% Lines 26-252 just calculate all of the equilibrium constants (I added the
%    solubility product of calcite as a bonus). These are all TEMPERATURE 
%    and PRESSURE corrected. You could simplify things by ignoring
%    pressure, but I see no reason not to include this.
%
% Lines 255-284 are the carbonate system solver. It takes in [DIC], [ALK], 
%    and the equilibrium constants and produces a solution for the full
%    carbonate system [pH, pOH, H2CO3, HCO3, CO3, CO2g, BOH4, BOH3]
%
% -------------------------------------------------------------------------

function pH = carbonate(x,y)

% x is ALK and y is DIC
% -------------------------------------------------------------------------
% some constants
% -------------------------------------------------------------------------
  B_T = 400;                     % total boron [umol/kg]
  T   = 25;                      % bottom water temp [degC]
  T_K = T + 273.15;              % converted temp [K]
  P   = 500;                     % pressure at depth [atm]
  S   = 35;                      % salinity [ppt]
%   DIC = 2100;                    % [DIC] [umol/kg]
%   Alk = 2300;                    % [ALK] [umol/kg]

% -------------------------------------------------------------------------
% solubility of CO2 in seawater
% -------------------------------------------------------------------------
% Note: units are mol/kg*atm (same as µmol/kg*ppm)
  K0 = exp( ...
      (9345.17/T_K)-60.2409+23.3585*log(T_K/100.) ...
      +S*( 0.023517-0.00023656*T_K  +  ...
           0.00000047036*T_K*T_K ));  
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% first dissociation constant of carbonic acid
% -------------------------------------------------------------------------
% Note: calculated in mol/kg, converted at end of function
  RGAS    = 8.314510;             % J mol-1 deg-1 (perfect Gas)  
  R       = 83.131;               % mol bar deg-1 
                                  % conversion cm3 -> m3          *1.e-6
                                  %            bar -> Pa = N m-2  *1.e+5
                                  %                => *1.e-1 or *1/10
  
  tmp1    = 2.83655 - 2307.1266 ./ T_K - 1.5529413 .* log(T_K);
  tmp2    =         - (0.20760841 + 4.0484 ./ T_K) .* sqrt(S);
  tmp3    =         + 0.08468345 .* S - 0.00654208 .* S .* sqrt(S);   
  tmp4    =         + log(1 - 0.001005 .* S);
  lnK1roy = tmp1 + tmp2 + tmp3 + tmp4;
  K1      = exp(lnK1roy);
  
% correct for pressure:
    % index: 
    %       1: K1 
    %       2: K2
    %       3: Kb
    %       4: Kw
    %       5: Kspc
    %       6: Kspa

  a0 = -[25.5    15.82  29.48   25.60  48.76   46.];
  a1 =  [0.1271 -0.0219 0.1622  0.2324 0.5304  0.5304];
  a2 =  [0.0     0.0    2.608  -3.6246 0.0     0.0] * 1.e-3;
  b0 = -[3.08   -1.13   2.84    5.13   11.76   11.76] * 1.e-3;
  b1 =  [0.0877 -0.1475 0.0     0.0794 0.3692  0.3692] * 1.e-3;
  b2 =  [0.0     0.0    0.0     0.0    0.0     0.0];

  for ipc=1:length(a0)
    deltav(ipc)    =  a0(ipc) + a1(ipc).*T + a2(ipc).*T.*T;
    deltak(ipc)    = (b0(ipc) + b1(ipc).*T + b2(ipc).*T.*T);  
    p.lnkpok0(ipc) = -(deltav(ipc)./(R.*T_K)).*P + (0.5*deltak(ipc)./(R.*T_K)).*P.*P;
  end

  K1      = K1*exp(p.lnkpok0(1));

% convert units (µmol/kg):
  K1      = K1*1.e6;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% second dissociation constant of carbonic acid
% -------------------------------------------------------------------------
% Note: calculated in mol/kg, converted at end of function
  RGAS = 8.314510;             % J mol-1 deg-1 (perfect Gas)  
  R    = 83.131;               % mol bar deg-1 
                               % conversion cm3 -> m3          *1.e-6
                               %            bar -> Pa = N m-2  *1.e+5
                               %                => *1.e-1 or *1/10
  
  tmp1    = -9.226508 - 3351.6106 ./ T_K - 0.2005743 .* log(T_K);
  tmp2    = (-0.106901773 - 23.9722 ./ T_K) .* sqrt(S);
  tmp3    = 0.1130822 .* S - 0.00846934 .* S.^1.5 + log(1 - 0.001005 * S);
  lnK2roy = tmp1 + tmp2 + tmp3;
  K2      = exp(lnK2roy);

% correct for pressure:
    % index: 
    %       1: K1 
    %       2: K2
    %       3: Kb
    %       4: Kw
    %       5: Kspc
    %       6: Kspa

  a0 = -[25.5    15.82  29.48   25.60  48.76   46.];
  a1 =  [0.1271 -0.0219 0.1622  0.2324 0.5304  0.5304];
  a2 =  [0.0     0.0    2.608  -3.6246 0.0     0.0] * 1.e-3;
  b0 = -[3.08   -1.13   2.84    5.13   11.76   11.76] * 1.e-3;
  b1 =  [0.0877 -0.1475 0.0     0.0794 0.3692  0.3692] * 1.e-3;
  b2 =  [0.0     0.0    0.0     0.0    0.0     0.0];

  for ipc=1:length(a0)
    deltav(ipc)    =  a0(ipc) + a1(ipc).*T + a2(ipc).*T.*T;
    deltak(ipc)    = (b0(ipc) + b1(ipc).*T + b2(ipc).*T.*T);  
    p.lnkpok0(ipc) = -(deltav(ipc)./(R.*T_K)).*P + (0.5*deltak(ipc)./(R.*T_K)).*P.*P;
  end

  K2      = K2*exp(p.lnkpok0(2));

% convert units (µmol/kg):
  K2      = K2*1.e6;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% dissociation constant for boric acid
% -------------------------------------------------------------------------
% Note: calculated in mol/kg, converted at end of function
  RGAS = 8.314510;             % J mol-1 deg-1 (perfect Gas)  
  R    = 83.131;               % mol bar deg-1 
                               % conversion cm3 -> m3          *1.e-6
                               %            bar -> Pa = N m-2  *1.e+5
                               %                => *1.e-1 or *1/10
  
  tmp1 =  (-8966.90-2890.53*sqrt(S)-77.942*S+1.728*S.^(3./2.)-0.0996*S.*S);
  tmp2 =   +148.0248+137.1942*sqrt(S)+1.62142*S;
  tmp3 = +(-24.4344-25.085*sqrt(S)-0.2474*S).*log(T_K);
  lnKb = tmp1 ./ T_K + tmp2 + tmp3 + 0.053105*sqrt(S).*T_K;
  Kb   = exp(lnKb);

% correct for pressure:
    % index: 
    %       1: K1 
    %       2: K2
    %       3: Kb
    %       4: Kw
    %       5: Kspc
    %       6: Kspa

  a0 = -[25.5    15.82  29.48   25.60  48.76   46.];
  a1 =  [0.1271 -0.0219 0.1622  0.2324 0.5304  0.5304];
  a2 =  [0.0     0.0    2.608  -3.6246 0.0     0.0] * 1.e-3;
  b0 = -[3.08   -1.13   2.84    5.13   11.76   11.76] * 1.e-3;
  b1 =  [0.0877 -0.1475 0.0     0.0794 0.3692  0.3692] * 1.e-3;
  b2 =  [0.0     0.0    0.0     0.0    0.0     0.0];

  for ipc=1:length(a0)
    deltav(ipc)    =  a0(ipc) + a1(ipc).*T + a2(ipc).*T.*T;
    deltak(ipc)    = (b0(ipc) + b1(ipc).*T + b2(ipc).*T.*T);  
    p.lnkpok0(ipc) = -(deltav(ipc)./(R.*T_K)).*P + (0.5*deltak(ipc)./(R.*T_K)).*P.*P;
  end

  Kb   = Kb*exp(p.lnkpok0(3));

% convert units (µmol/kg):
  Kb   = Kb*1.e6;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% ion product for water
% -------------------------------------------------------------------------
% Note: calculated in mol^2/kg^2, converted at end of function
  RGAS = 8.314510;             % J mol-1 deg-1 (perfect Gas)  
  R    = 83.131;               % mol bar deg-1 
                               % conversion cm3 -> m3          *1.e-6
                               %            bar -> Pa = N m-2  *1.e+5
                               %                => *1.e-1 or *1/10
  
  Kw = exp( 148.96502-13847.26/T_K-23.6521*log(T_K) ...
           +(118.67/T_K-5.977+1.0495*log(T_K))*sqrt(S)-0.01615*S );

% correct for pressure:
    % index: 
    %       1: K1 
    %       2: K2
    %       3: Kb
    %       4: Kw
    %       5: Kspc
    %       6: Kspa

  a0 = -[25.5    15.82  29.48   25.60  48.76   46.];
  a1 =  [0.1271 -0.0219 0.1622  0.2324 0.5304  0.5304];
  a2 =  [0.0     0.0    2.608  -3.6246 0.0     0.0] * 1.e-3;
  b0 = -[3.08   -1.13   2.84    5.13   11.76   11.76] * 1.e-3;
  b1 =  [0.0877 -0.1475 0.0     0.0794 0.3692  0.3692] * 1.e-3;
  b2 =  [0.0     0.0    0.0     0.0    0.0     0.0];

  for ipc=1:length(a0)
    deltav(ipc)    =  a0(ipc) + a1(ipc).*T + a2(ipc).*T.*T;
    deltak(ipc)    = (b0(ipc) + b1(ipc).*T + b2(ipc).*T.*T);  
    p.lnkpok0(ipc) = -(deltav(ipc)./(R.*T_K)).*P + (0.5*deltak(ipc)./(R.*T_K)).*P.*P;
  end

  Kw = Kw*exp(p.lnkpok0(4));

% convert units (µmol^2/kg^2):
  Kw = Kw*1.e12;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% solubility product constant for calcite
% -------------------------------------------------------------------------
% Note: calculated in mol^2/kg^2, converted at end of function
  RGAS = 8.314510;             % J mol-1 deg-1 (perfect Gas)  
  R    = 83.131;               % mol bar deg-1 
                               % conversion cm3 -> m3          *1.e-6
                               %            bar -> Pa = N m-2  *1.e+5
                               %                => *1.e-1 or *1/10

  tmp1      = -171.9065-0.077993.*T_K+2839.319./T_K+71.595.*log10(T_K);
  tmp2      = +(-0.77712+0.0028426.*T_K+178.34./T_K).*sqrt(S);
  tmp3      = -0.07711.*S+0.0041249.*S.^1.5;
  log10Kspc = tmp1 + tmp2 + tmp3;
  Kspc      = 10.^(log10Kspc);

% correct for pressure:
    % index: 
    %       1: K1 
    %       2: K2
    %       3: Kb
    %       4: Kw
    %       5: Kspc
    %       6: Kspa

  a0 = -[25.5    15.82  29.48   25.60  48.76   46.];
  a1 =  [0.1271 -0.0219 0.1622  0.2324 0.5304  0.5304];
  a2 =  [0.0     0.0    2.608  -3.6246 0.0     0.0] * 1.e-3;
  b0 = -[3.08   -1.13   2.84    5.13   11.76   11.76] * 1.e-3;
  b1 =  [0.0877 -0.1475 0.0     0.0794 0.3692  0.3692] * 1.e-3;
  b2 =  [0.0     0.0    0.0     0.0    0.0     0.0];

  for ipc=1:length(a0)
    deltav(ipc)    =  a0(ipc) + a1(ipc).*T + a2(ipc).*T.*T;
    deltak(ipc)    = (b0(ipc) + b1(ipc).*T + b2(ipc).*T.*T);  
    p.lnkpok0(ipc) = -(deltav(ipc)./(R.*T_K)).*P + (0.5*deltak(ipc)./(R.*T_K)).*P.*P;
  end

Kspc      = Kspc*exp(p.lnkpok0(5));

% convert units (µmol^2/kg^2):
out       = Kspc*1.e12;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% csys3 routine for solving carbonate system
% -------------------------------------------------------------------------
% calculate H+ abundance based on DIC and Alk:
% x is ALK and y is DIC


  p5  = -1.;
  p4  = -x-Kb-K1;
  p3  = y*K1-x*(Kb+K1)+Kb*B_T+Kw-Kb*K1-K1*K2;
  tmp = y*(Kb*K1+2.*K1*K2)-x*(Kb*K1+K1*K2)+Kb*B_T*K1;
  p2  = tmp+(Kw*Kb+Kw*K1-Kb*K1*K2);
  tmp = 2.*y*Kb*K1*K2-x*Kb*K1*K2+Kb*B_T*K1*K2;
  p1  = tmp+(+Kw*Kb*K1+Kw*K1*K2);
  p0  = Kw*Kb*K1*K2;
  p   = [p5 p4 p3 p2 p1 p0];
  r   = roots(p);
 
    
% calculate other carbonate system variables:
  H     = max(real(r));
  pH    = -log10(H)+6;
  OH    = Kw/H;
  H2CO3 = y / (1.+K1/H+K1*K2/H/H);
  HCO3  = y / (1+H/K1+K2/H);
  CO3   = y / (1+H/K2+H*H/K1/K2);
  CO2g  = H2CO3/K0;
  BOH4  = (B_T*Kb)/(Kb+H);
  BOH3  = (H*BOH4)/Kb;
  
  
end
% -------------------------------------------------------------------------