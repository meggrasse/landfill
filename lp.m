% this code calculates the profit for the different permutations of a 
% composite liner system. 
% the next step is to explore if increasing thickness could increase profit
% by decreasing leachate.

% the parameters to calculating `q` using the composite liner equation
% for each permuation of the composite liner system (high vs. low quality
% clay and good vs. bad contact]
%
% params:
% [n, c, a, k]
% n: number of holes
% c: factor per good contact/poor contact
% a: hole area (m^2)
% k: hydraulic conductivity of mineral layer
Q_params = [ 20 1.15 0.05 1*10^-8 
      5 0.21 0.01 1*10^-8 
      20 1.15 0.05 1*10^-9
      5 0.21 0.01 1*10^-9 
    ];

area_m2 = 10000;
% Q (m^3/s) for each lining option
Q = (Q_params(:,1).*Q_params(:,2).*1.^0.9.*Q_params(:,3).^0.1.*Q_params(:,4).^0.74);

seconds_per_year = 365*24*60*60;

% parameters to calculate cost. common sense assumption as LLDPE as the
% geomembrane (cheapest) and that composite liners require a protection
% layer.

% params:
% [protection, geomembrane, clay, leachate]
C_params = [ 4*area_m2 8*area_m2 20*area_m2 0 300*60*Q(1)*seconds_per_year
      4*area_m2 8*area_m2 20*area_m2 50*area_m2 300*60*Q(2)*seconds_per_year
      4*area_m2 8*area_m2 100*area_m2 0 300*60*Q(3)*seconds_per_year
      4*area_m2 8*area_m2 100*area_m2 50*area_m2 300*60*Q(4)*seconds_per_year
    ];

clay_thickness = 1;
C = C_params(:,1)+C_params(:,2)+C_params(:,3)*clay_thickness+C_params(:,4)+C_params(:,5);

% income assuming settlement to meet regulation
I=7.50*(30-1)/.8 * area_m2;

%profit
P=I-C

% quick expl. of mineral layer only
% thickness 1m
Q_min=10^-9*(1+1)/1;
% high quality clay, leakage
C_min=100*Q_min*(300*60*365*24*60*60);
% income per m2
I_min=7.50*(30-1)/.8;
P_min = I_min-C_min;

% next steps:
% 1. explore increasing thickness with additional layers for current
% optimal solution
% 2. add in drainage costs
% 3. add in capping costs




