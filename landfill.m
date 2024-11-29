% this code calculates the profit for the different permutations of a 
% composite liner system. 

% the parameters to calculating `Q` using the composite liner equation
% for each permuation of the composite liner system (high vs. low quality
% clay and good vs. bad contact]
%
% params:
% [n, c, a, k]
% n: number of holes
% c: factor per good contact/poor contact
% a: hole area (m^2)
% k: hydraulic conductivity of mineral layer
Q_params = [ 20 1.15 0.05 1*10^-8 %poor contact, poor mineral layer
      5 0.21 0.01 1*10^-8 %good contact, poor mineral layer
      20 1.15 0.05 1*10^-9 %poor contact, good mineral layer
      5 0.21 0.01 1*10^-9 %good contact, good mineral layer
    ];

area_m2 = 10000;
% Q (m^3/s) for each lining option.
% implements the `Q=n*0.21(h^0.9*a^0.1*k^0.74)` [poor contact] or
% `Q=n*0.21(h^0.9*a^0.1*k^0.74)` [good contact] G&B (1989) equation.
Q = (Q_params(:,1).*Q_params(:,2).*1.^0.9.*Q_params(:,3).^0.1.*Q_params(:,4).^0.74);

seconds_per_year = 365*24*60*60;

% parameters to calculate cost, following the same row structure as
% `Q_params` for the different options.
% 
% assuming that a composite liner must require a protection layer.
%
% params:
% [pl, gm, gml, cqa, l]
% pl: cost of protection layer (per m^2)
% gm: cost of geomembrane (per m^2) [using LLDPE for now as it's cheapest]
% ml: cost of mineral layer (per m^2)
% cqa: cost of CQA (per m^2)
% l: cost of pumping and treating leakage. calculated by multiplying
% quanitity (flow rate x time) by cost (per landfill area)

landfill_lifetime_yrs = 60;
C_params = [ 4 8 20 0 landfill_lifetime_yrs*seconds_per_year*Q(1)*300 %poor contact, poor mineral layer
      4 8 20 50 landfill_lifetime_yrs*seconds_per_year*Q(2)*300 %good contact, poor mineral layer
      4 8 100 0 landfill_lifetime_yrs*seconds_per_year*Q(3)*300 %poor contact, good mineral layer
      4 8 100 50 landfill_lifetime_yrs*seconds_per_year*Q(4)*300 %good contact, good mineral layer
    ];
% convert all costs to cost per landfill area 
C_params = [C_params(:,1:4).*area_m2, C_params(:,5)];

clay_thickness = 1;
% add up the costs for each option, multiplying by 1m of clay thickness
C = C_params(:,1)+C_params(:,2)+C_params(:,3)*clay_thickness+C_params(:,4)+C_params(:,5);

% income assuming settlement to meet regulation
I=7.50*(30-1)/.8 * area_m2;

%profit
P=I-C

% this is a quick exploration of just using a 1m "good quality" clay layer 
% with 10^-9 permeability; the mininmum needed to achieve the landfill 
% directive.

% calculate hydraulic gradient via `i=h+D/D` where h=1m and D(thickness)=1m
i_ml = 1+1/1;
% calcuate's darcy's law: q=k*i*a 
Q_good_ml=10^-9*i_ml*area_m2;
% cost: cost of mineral layer + leachate quanitity cost
C_good_ml=100+area_m2 + landfill_lifetime_yrs*seconds_per_year*Q_good_ml*300;
% income per m2
I_good_ml=7.50*(30-1)/.8 * area_m2;
P_good_ml = I_good_ml-C_good_ml

% next steps:
% 1. explore increasing thickness with additional layers for current
% optimal solution
% 2. add in drainage costs
% 3. add in capping costs




