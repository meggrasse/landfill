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
C_materials = C_params(:,1)+C_params(:,2)+C_params(:,3)*clay_thickness+C_params(:,4);
C_leakage = C_params(:,5);

disp("for a single liner");
disp("the cost of materials is :");
disp(C_materials);
disp("the cost of leakage is :");
disp(C_leakage);

C = C_materials + C_leakage;

disp("the total cost is :");
disp(C); 

% income assuming settlement to meet regulation
I=7.50*(30-1)/.8 * area_m2;

disp("the total income is :");
disp(I);

%profit
P=I-C

disp("thus, profit is :");
disp(P);

% calculate k
% thickness is 1m
d = 1;
% head is 1m
h = 1;
i = (h+d)/d;
k = Q/(i*area_m2);

disp("by the way, k for each option is:");
disp(k);

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

% this section evaluates a double liner system: 
%      - top: geomembrane
%      - bottom: composite liner
% if we're assuming a head on the bottom liner, we actually don't even need
% to calcuate a flow through the top.

% params:
% [n, c, h, a, k]
% n: number of holes
% c: factor per good contact/poor contact
% a: hole area (m^2)
% k: hydraulic conductivity of mineral layer
Q_params_bottom = [
    20 1.15 0.05 10^-8 % no CQA, expected head
    5 0.21 0.01 10^-8 % CQA, expected head
    20 1.15 0.05 10^-9 % no CQA, worst case head 
    5 0.21 0.01 10^-9 % CQA, worst case head
    ];
Q_bottom_worst_case = Q_params_bottom(:,1).*Q_params_bottom(:,2).*0.1.^0.9.*Q_params_bottom(:,3).^0.1.*Q_params_bottom(:,4).^0.74;
Q_bottom_expected = Q_params_bottom(:,1).*Q_params_bottom(:,2).*0.01.^0.9.*Q_params_bottom(:,3).^0.1.*Q_params_bottom(:,4).^0.74;

% parameters to calculate cost, following the same row structure as
% `Q_params` for the different options.
% 
% assuming that a composite liner must require a protection layer.
%
% params:
% [tl, idr, pl, gm, gml, cqa, l]
% tl: cost of top liner
% idr: cost of intermediary drainage layer (MG: i'm making this as 0.25 for
% now)
% pl: cost of protection layer (per m^2)
% gm: cost of geomembrane (per m^2) [using LLDPE for now as it's cheapest]
% ml: cost of mineral layer (per m^2)
% cqa: cost of CQA (per m^2)
% l: cost of pumping and treating leakage. calculated by multiplying
% quanitity (flow rate x time) by cost (per landfill area)

landfill_lifetime_yrs = 60;
C_params_double_expected = [ 8 10*.25 4 8 20 0 landfill_lifetime_yrs*seconds_per_year*Q_bottom_expected(1)*300 %poor contact, poor mineral layer
      8 10*.25 4 8 20 50 landfill_lifetime_yrs*seconds_per_year*Q_bottom_expected(2)*300 %good contact, poor mineral layer
      8 10*.25 4 8 100 0 landfill_lifetime_yrs*seconds_per_year*Q_bottom_expected(3)*300 %poor contact, good mineral layer
      8 10*.25 4 8 100 50 landfill_lifetime_yrs*seconds_per_year*Q_bottom_expected(4)*300 %good contact, good mineral layer
    ];
C_params_double_worst_case = [ 8 10*.25 4 8 20 0 landfill_lifetime_yrs*seconds_per_year*Q_bottom_worst_case(1)*300 %poor contact, poor mineral layer
      8 10*.25 4 8 20 50 landfill_lifetime_yrs*seconds_per_year*Q_bottom_worst_case(2)*300 %good contact, poor mineral layer
      8 10*.25 4 8 100 0 landfill_lifetime_yrs*seconds_per_year*Q_bottom_worst_case(3)*300 %poor contact, good mineral layer
      8 10*.25 4 8 100 50 landfill_lifetime_yrs*seconds_per_year*Q_bottom_worst_case(4)*300 %good contact, good mineral layer
    ];

% assumes by 1m of clay thickness
C_materials_double = sum(C_params_double_expected(:,1:end-1).*area_m2, 2);
C_leakage_double_expected = C_params_double_expected(:,end);
C_leakage_double_worst_case = C_params_double_worst_case(:,end);

disp("for a double liner");
disp("the cost of materials is :");
disp(C_materials_double);
disp("the cost of leakage is in the expected case is :");
disp(C_leakage_double_expected);
disp("and in the worst case is :");
disp(C_leakage_double_worst_case);

C_double_expected = C_materials_double + C_leakage_double_expected;
C_double_worst_case = C_materials_double + C_leakage_double_worst_case;

disp("the total cost is :");
disp(C_double_expected); 
disp("and in the worst case:");
disp(C_double_worst_case); 

% income assuming settlement to meet regulation
I=7.50*(30-1.25)/.8 * area_m2;

disp("the total income is :");
disp(I);

%profit
P_expected=I-C_double_expected
P_worst_case = I-C_double_worst_case

disp("thus, profit is :");
disp(P_expected);
disp("and in the worst case")
disp(P_worst_case);

% next steps:
% 1. try out GCL with our own clay?
% 2. add in drainage costs
% 3. add in capping costs

