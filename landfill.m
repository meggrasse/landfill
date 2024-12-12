% this code calculates the profit for the different permutations of a 
% basal liner system, including the drainage and capping system.

linerKind = struct('MINERAL_LINER', 1, 'SINGLE_COMPOSITE', 2, 'DOUBLE_LINER', 3);
mineralLinerPermeability = struct('LOW_PERMEABILITY_CLAY', 10^-9, 'SEMI_LOW_PERMEABILITY_CLAY', 10^-8, 'GCL', 10^-10);

LinerKind = [linerKind.MINERAL_LINER linerKind.MINERAL_LINER, linerKind.SINGLE_COMPOSITE, linerKind.SINGLE_COMPOSITE, linerKind.SINGLE_COMPOSITE, linerKind.SINGLE_COMPOSITE, linerKind.SINGLE_COMPOSITE, linerKind.SINGLE_COMPOSITE, linerKind.DOUBLE_LINER, linerKind.DOUBLE_LINER, linerKind.DOUBLE_LINER, linerKind.DOUBLE_LINER, linerKind.DOUBLE_LINER, linerKind.DOUBLE_LINER].';
CQA = [false false false true false true false true false true false true false true].';
MineralLinerPermeability = [mineralLinerPermeability.LOW_PERMEABILITY_CLAY, mineralLinerPermeability.GCL, mineralLinerPermeability.SEMI_LOW_PERMEABILITY_CLAY, mineralLinerPermeability.SEMI_LOW_PERMEABILITY_CLAY, mineralLinerPermeability.LOW_PERMEABILITY_CLAY, mineralLinerPermeability.LOW_PERMEABILITY_CLAY, mineralLinerPermeability.GCL, mineralLinerPermeability.GCL, mineralLinerPermeability.SEMI_LOW_PERMEABILITY_CLAY, mineralLinerPermeability.SEMI_LOW_PERMEABILITY_CLAY, mineralLinerPermeability.LOW_PERMEABILITY_CLAY, mineralLinerPermeability.LOW_PERMEABILITY_CLAY, mineralLinerPermeability.GCL, mineralLinerPermeability.GCL].';

% Part 3: Calculate Leakage

seconds_per_year = 365*24*60*60;
head = 1;
head_bottom_liner = 0.1; % worst case
thickness = 1;
leakage_rate = 0;
hole_count = 0;
contact_factor = 0;
hole_area = 0;
indexes = length(LinerKind);
LeakageRate = zeros(indexes,1);
for ii = 1:indexes
    leakage_rate = 0;
    if CQA(ii)
        hole_count = 5;
        contact_factor = 0.21;
        hole_area = 0.01;
    elseif MineralLinerPermeability(ii) == mineralLinerPermeability.GCL
        % GCL can assume good contact (G&P)
        hole_count = 20;
        contact_factor = 0.21;
        hole_area = 0.05;
    else
        hole_count = 20;
        contact_factor = 1.15;
        hole_area = 0.05;
    end

    switch LinerKind(ii)
        case linerKind.MINERAL_LINER
            hydraulic_gradient = (thickness+head)/thickness;
            % Darcy's Law
            leakage_rate = MineralLinerPermeability(ii)*hydraulic_gradient*area_m2;
        case linerKind.SINGLE_COMPOSITE
            % leakage through a single composite liner (Bonaparte et al., 1989)
            leakage_rate = hole_count*contact_factor*head^0.9*hole_area^0.1*MineralLinerPermeability(ii)^0.74;
        case linerKind.DOUBLE_LINER
            % leakage through a single composite liner (Bonaparte et al., 1989)
            leakage_rate = hole_count*contact_factor*head_bottom_liner^0.9*hole_area^0.1*MineralLinerPermeability(ii)^0.74;
        otherwise
            warning('Unexpected linerKind.')
    end
    LeakageRate(ii,1) = leakage_rate;
end

LifetimeLeakage = LeakageRate.*seconds_per_year*landfill_lifetime_yrs;

% TODO: add seperator geotextile in for intermediary drainage layer in basal system?
costs = struct('LLDPE', 8, 'TYRES', 10, 'PROTECTION_GEOTEXTILE', 4, 'LOW_PERMEABILITY_CLAY', 100, 'SEMI_LOW_PERMEABILITY_CLAY', 20, 'GCL', 15, 'CQA', 50, 'SEPARATOR_GEOTEXTILE', 2, 'DRAINAGE_GRAVEL', 50, 'RESTORATION_SOILS', 1);

MaterialCost = zeros(indexes,1);
current_material_cost = 0; % per m^2
for ii = 1:indexes
    current_material_cost = 0;
    if CQA(ii)
        current_material_cost = current_material_cost + costs.CQA;
    end

    switch LinerKind(ii)
        case linerKind.SINGLE_COMPOSITE
            current_material_cost = current_material_cost + costs.PROTECTION_GEOTEXTILE + costs.LLDPE;
        case linerKind.DOUBLE_LINER
            % tyre layer height of 10cm
            current_material_cost = current_material_cost + costs.LLDPE + costs.TYRES*0.1 + costs.PROTECTION_GEOTEXTILE + costs.LLDPE;
    end

    switch MineralLinerPermeability(ii)
        case mineralLinerPermeability.GCL
            current_material_cost = current_material_cost + costs.GCL;
        case mineralLinerPermeability.LOW_PERMEABILITY_CLAY
            current_material_cost = current_material_cost + costs.LOW_PERMEABILITY_CLAY;
        case mineralLinerPermeability.SEMI_LOW_PERMEABILITY_CLAY
            current_material_cost = current_material_cost + costs.SEMI_LOW_PERMEABILITY_CLAY;
    end
    current_material_cost = current_material_cost*area_m2;
    MaterialCost(ii) = current_material_cost;
end

LifetimeLeakageCost = LifetimeLeakage.*300;

thickness = 0;
head = 1;
hydraulic_gradient = 0;
Permeability = zeros(indexes,1);
for ii = 1:indexes
    hydraulic_gradient = 0;
    if LinerKind(ii) == linerKind.DOUBLE_LINER
        % includes intermediary drainage layer
        thickness = 1.1;
    else
        thickness = 1;
    end

    hydraulic_gradient = (thickness + head)/thickness;
    Permeability(ii) = LeakageRate(ii)/(hydraulic_gradient*area_m2);
end

total_cover_soil_cost = 134523;
drainage_cost = 170925;
capping_cost = costs.RESTORATION_SOILS*1 + costs.RESTORATION_SOILS*1 + costs.SEPARATOR_GEOTEXTILE + costs.DRAINAGE_GRAVEL*0.5 + costs.PROTECTION_GEOTEXTILE + costs.LLDPE + costs.PROTECTION_GEOTEXTILE + costs.DRAINAGE_GRAVEL*0.3 + costs.SEPARATOR_GEOTEXTILE;
capping_cost = capping_cost*area_m2;

available_volume = 143721.2;
income = available_volume*7.5;

TotalCost = MaterialCost+LifetimeLeakageCost+total_cover_soil_cost+drainage_cost+capping_cost;
Profit = income-TotalCost;

table(LinerKind, CQA, MineralLinerPermeability, LeakageRate, LifetimeLeakage, Permeability, MaterialCost, LifetimeLeakageCost, TotalCost, Profit)

% hold on
% Q_temp = [Q_good_ml; Q_gcl; Q; Q_bottom_worst_case]
% % need to fix calclulations
% P_temp = [P_good_ml; P_gcl; P; P_worst_case]
% for ii = 1:length(Q_temp)
% %     colormap = colormap winter;
% %     if ii > 2
% %         colormap = colormap("summer");
% %     end
% %     if ii > 8
% %          colormap = colormap("default");
% %     end
% %     color = colormap[ii];
%     scatter(Q_temp(ii),P_temp(ii))
% end
% legend("LP", "GCL", "GM+HP(no CQA)","GM+HP(CQA)","GM+LP(no CQA)","GM+LP(CQA)", "GM+GCL(no CQA)", "GM+GCL(CQA)", "GM+GM+HP(CQA)","GM+GM+LP(no CQA)","GM+GM+LP(CQA)", "GM+GM+GCL(no CQA)", "GM+GM+GCL(CQA)");

