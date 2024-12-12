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

disp("for the top or single liner");
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

% this section evaluates double liner systems: 
%     - a top geomembrane and a bottom composite liner
%     - a double composite liner
%
% without a derived equation to calculate the resulting flow, we've
% made the following assumptions:
%    - the flow from the top liner will comprise the head of the bottom
%    liner
%    - if the head of the bottom liner reaches 1m it will be drained by
%    the intermediary drainage system

% we assume the top liner system is good contact, poor quality soil because
% that yielded highest profit in section 1.

hours_per_year = 365*60*24;

% i'm assuming from part 3 if we have CQA we have less and smaller holes
% [n, a]
% n: number of holes
% a: hole area (m^2)
Q_params_geomembrane_top = [
    20 0.05 % low quality
    5 0.01 % high quality
    20 0.05 % low quality
    5 0.01 % high quality
    ];
% calculate with geomembrane. m^3/s
Q_geomembrane_top = Q_params_geomembrane_top(:,1).*0.6.*Q_params_geomembrane_top(:,2).*sqrt(2*9.81*1)

% follows the same row structure as `Q_params_top`.
% assuming GCL permeability of 1*10^-10--12 (c.f. Waste Treatment and
% Disposal p. 187).
%
% params:
% [n, c, a, k]
% n: number of holes
% c: factor per good contact/poor contact
% a: hole area (m^2)
% k: hydraulic conductivity of mineral layer
Q_params_composite_bottom = [ 20 1.15 0.05 1*10^-12 %poor contact, GCL
                    5 0.21 0.01 1*10^-12 %good contact, GCL
                    20 1.15 0.05 1*10^-12 %poor contact, GCL
                    5 0.21 0.01 1*10^-12 %good contact, GCL
    ];

% this is the head on the bottom liner taken at hourly increments over the
% lifetime of the landfill.
% initialized to 0s.
%head_bottom_hourly = zeros(4,hours_per_year*landfill_lifetime_yrs);
% this is the leakage of the bottom liner at each hour according to the
% dynamic head. in m^3/s.
% initialized to 0s.
q_bottom_composite = zeros(4,hours_per_year*landfill_lifetime_yrs);
q_bottom_geomembrane_composite = zeros(4,hours_per_year*landfill_lifetime_yrs);
% this is the head at the current time interval.
current_head_composite = zeros(4,1);
current_head_geomembrane = zeros(4,1);
% this is the volume of leakage that flowed out of the bottom liner during
% the last hour. this should be removed from the calculation of the head.
bottom_q_volume_composite = zeros(4,1);
bottom_q_volume_geomembrane = zeros(4,1);
% the volume of leakage out of the top liner in one hour.
top_q_volume_composite = zeros(4,1);
top_q_volume_geomembrane = zeros(4,1);
% this loop interates through each hour in the lifetime of the landfill to
% calculate the `q` of the bottom liner at that time.
for hour = 1:hours_per_year*landfill_lifetime_yrs
    % the volume of leakage out of the top liner in one hour.
    top_q_volume_composite = Q*60*60;
    top_q_volume_geomembrane = Q_geomembrane_top*60*60;
    % to calculate the current head we subtract what flowed out of the
    % bottom liner from the previous head and add what flowed from the top liner 
    current_head_composite = current_head_composite - bottom_q_volume_composite + top_q_volume_composite;
    current_head_geomembrane = current_head_geomembrane - bottom_q_volume_geomembrane + top_q_volume_geomembrane;
    % we're assuming the operator maintains a max of 1m leachate head MG IS
    % THIS A FAIR ASSUMPTION??
    % this is actually really important to do because the bottom liner is 
    % less permeable, the head head grows from flow out of the top liner.
    current_head_composite = min(current_head_composite, 1);
    % MG: for the geomembrane can it be even less?
    current_head_geomembrane = min(current_head_geomembrane, 1);
    %head_bottom_hourly(:,hour) = current_head_composite;
    bottom_q_volume_composite = Q_params_composite_bottom(:,1).*Q_params_composite_bottom(:,2).*current_head_composite.^0.9.*Q_params_composite_bottom(:,3).^0.1.*Q_params_composite_bottom(:,4).^0.74;
    bottom_q_volume_geomembrane = Q_params_composite_bottom(:,1).*Q_params_composite_bottom(:,2).*current_head_geomembrane.^0.9.*Q_params_composite_bottom(:,3).^0.1.*Q_params_composite_bottom(:,4).^0.74;
    q_bottom_composite(:,hour) = bottom_q_volume_composite;
    q_bottom_geomembrane_composite(:,hour) = bottom_q_volume_geomembrane;
end

% to calculate the total leakage over the lifetime of the landfill then
% would be to convert the q of the bottom liner from m^3/s to m^3/h and sum
% up for the lifetime of the landfill.
q_bottom_composite_per_hour = q_bottom_composite.*60*60;
q_dbl_composite_lifetime = sum(q_bottom_composite_per_hour,2);
q_bottom_geomembrane_composite_per_hour = q_bottom_geomembrane_composite.*60*60;
q_geomembrane_composite_lifetime = sum(q_bottom_geomembrane_composite_per_hour,2);


% similar utility as `C_params_top`; follows same row structure.
% assuming we don't have to pay for CQA  for good contact solutions 
% because we've paid for it for the top liner.
%
% params:
% [pl, gm, gml, l]
% pl: cost of protection layer (per m^2)
% dl: cost of drainage layer (per m^2)
% gm: cost of geomembrane (per m^2) [using LLDPE for now as it's cheapest]
% ml: cost of mineral layer (per m^2)
% l: cost of pumping and treating leakage. calculated by multiplying
% quanitity (flow rate x time) by cost (per landfill area)
C_params_bottom = [4 2 8 15 q_dbl_composite_lifetime(1)*300 % protection, seperator (MG: i'm assuming this is a geonet), LLDPE, GCL
                   4 2 8 15 q_dbl_composite_lifetime(2)*300 % protection, seperator (MG: i'm assuming this is a geonet), LLDPE, GCL
                   4 2 8 15 q_dbl_composite_lifetime(3)*300 % protection, seperator (MG: i'm assuming this is a geonet), LLDPE, GCL
                   4 2 8 15 q_dbl_composite_lifetime(4)*300 % protection, seperator (MG: i'm assuming this is a geonet), LLDPE, GCL
    ];

C_params_bottom_geomembrane = [4 2 8 15 q_geomembrane_composite_lifetime(1)*300 % protection, seperator (MG: i'm assuming this is a geonet), LLDPE, GCL
                   4 2 8 15 q_geomembrane_composite_lifetime(2)*300 % protection, seperator (MG: i'm assuming this is a geonet), LLDPE, GCL
                   4 2 8 15 q_geomembrane_composite_lifetime(3)*300 % protection, seperator (MG: i'm assuming this is a geonet), LLDPE, GCL
                   4 2 8 15 q_geomembrane_composite_lifetime(4)*300 % protection, seperator (MG: i'm assuming this is a geonet), LLDPE, GCL
    ];

% cost includes seperator geotextile (geonet) and LLDPE
C_params_geoembrane = [2 8
                   2 8
                   2 8
                   2 8
    ];

% convert all costs to cost per landfill area 
C_materials_bottom = sum(C_params_bottom(:,1:4).*area_m2,2);
C_leakage_bottom = C_params_bottom(:,5);
C_bottom = C_materials_bottom + C_leakage_bottom;

disp("for the composite bottom liner");
disp("the cost of materials is :");
disp(C_materials_bottom);
disp("the cost of leakage is :");
disp(C_leakage_bottom);
disp("the total cost of the bottom liner is :");
disp(C_bottom); 

% cost of top liner excluding the leakage cost
C_top = sum(C_params(:,1:end-1),2);
disp("the cost of the top liner without leakage cost :");
disp(C_top); 
C_dbl_liner = C_bottom + C_top;

disp("the total cost of the double lining system is :");
disp(C_dbl_liner); 

% income assuming settlement to meet regulation
% GCL is 1cm thick
I_dbl_liner=7.50*(30-1.01)/.8 * area_m2;

disp("the total income is :");
disp(I_dbl_liner);

%profit
P_dbl_liner = I_dbl_liner-C_dbl_liner

disp("thus, profit is :");
disp(P_dbl_liner);

% cost of top liner excluding the leakage cost
C_top = sum(C_params(:,1:end-1),2);
disp("the cost of the top liner without leakage cost :");
disp(C_top); 
C_dbl_liner = C_bottom + C_top;

C_materials_geomembrane= sum(C_params_geoembrane(:,1:2).*area_m2,2);
C_leakage_geomembrane_composite = C_params_bottom_geomembrane(:,5);

disp("for the geomembrane top liner");
disp("the cost of the geomembrane top liner without leakage cost :");
disp(C_materials_geomembrane);

disp("the cost of bottom composite materials is :");
disp(C_materials_bottom);
disp("the cost of leakage is :");
disp(C_leakage_geomembrane_composite);
C_geomembrane_composite_bottom = C_leakage_geomembrane_composite + C_materials_bottom;
disp("the total cost of the bottom liner is :");
disp(C_geomembrane_composite_bottom); 

C_geomembrane_composite = C_geomembrane_composite_bottom + C_materials_geomembrane;

disp("the total cost of the geomembrane and composite lining system is :");
disp(C_geomembrane_composite); 

% income assuming settlement to meet regulation
% GCL is 1cm thick
I_geomembrane_composite=7.50*(30-1.00)/.8 * area_m2;

disp("the total income is :");
disp(I_geomembrane_composite);

%profit
P_geomembrane_composite = I_geomembrane_composite-C_geomembrane_composite;
disp("thus, profit is :");
disp(P_geomembrane_composite);

%ah wait, this exploration isn't valid because i'm just calcuating a
%geomembrane on top of a GCL

% next steps:
% 1. try out other multiple liner systems in the textbook
% 2. add in drainage costs
% 3. add in capping costs

% possible errors:
% 1. is hole area right?
% 2. understand cost of lechate vs. cost of materials for each option

% oh, my assumption is broken because we end up paying for 


