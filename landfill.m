% this code calculates the profit for the different permutations of a 
% basal liner system, including the drainage and capping system.

linerKind = struct('MINERAL_LINER', 1, 'SINGLE_COMPOSITE', 2, 'DOUBLE_LINER_COMPOSITE', 3, 'GEOMEMBRANE', 4, 'DOUBLE_GEOMEMBRANE', 5);
mineralLinerPermeability = struct('LOW_PERMEABILITY_CLAY', 10^-9, 'SEMI_LOW_PERMEABILITY_CLAY', 10^-8, 'GCL', 10^-10);

LinerKind = [linerKind.MINERAL_LINER linerKind.MINERAL_LINER, linerKind.SINGLE_COMPOSITE, linerKind.SINGLE_COMPOSITE, linerKind.SINGLE_COMPOSITE, linerKind.SINGLE_COMPOSITE, linerKind.SINGLE_COMPOSITE, linerKind.SINGLE_COMPOSITE, linerKind.DOUBLE_LINER_COMPOSITE, linerKind.DOUBLE_LINER_COMPOSITE, linerKind.DOUBLE_LINER_COMPOSITE, linerKind.DOUBLE_LINER_COMPOSITE, linerKind.DOUBLE_LINER_COMPOSITE, linerKind.DOUBLE_LINER_COMPOSITE linerKind.GEOMEMBRANE linerKind.GEOMEMBRANE linerKind.DOUBLE_GEOMEMBRANE linerKind.DOUBLE_GEOMEMBRANE].';
CQA = [false false false true false true false true false true false true false true false true false true].';
MineralLinerPermeability = [mineralLinerPermeability.LOW_PERMEABILITY_CLAY, mineralLinerPermeability.GCL, mineralLinerPermeability.SEMI_LOW_PERMEABILITY_CLAY, mineralLinerPermeability.SEMI_LOW_PERMEABILITY_CLAY, mineralLinerPermeability.LOW_PERMEABILITY_CLAY, mineralLinerPermeability.LOW_PERMEABILITY_CLAY, mineralLinerPermeability.GCL, mineralLinerPermeability.GCL, mineralLinerPermeability.SEMI_LOW_PERMEABILITY_CLAY, mineralLinerPermeability.SEMI_LOW_PERMEABILITY_CLAY, mineralLinerPermeability.LOW_PERMEABILITY_CLAY, mineralLinerPermeability.LOW_PERMEABILITY_CLAY, mineralLinerPermeability.GCL, mineralLinerPermeability.GCL 0 0 0 0].';

AREA_M2 = 10000;
LANDFILL_LIFETIME_YRS = 60;
SECONDS_PER_DAY = 60*60*24;
SECONDS_PER_YEAR = 365*SECONDS_PER_DAY;
HEAD = 1;
indexes = length(LinerKind);

%% Part 3: Calculate leakage

% calculate i
thickness = 0;
hydraulic_gradient = 0;
HydraulicGradient = zeros(indexes,1)
for ii = 1:indexes
    hydraulic_gradient = 0;
    if LinerKind(ii) == linerKind.DOUBLE_LINER_COMPOSITE | LinerKind(ii) == linerKind.DOUBLE_GEOMEMBRANE
        % includes intermediary drainage layer
        thickness = 0.6;
    else
        thickness = 0.5;
    end

    hydraulic_gradient = (thickness + HEAD)/thickness;
    HydraulicGradient(ii) = hydraulic_gradient;
end

% calculate q
HEAD_BOTTOM_LINER = 0.1; % worst case head of 0.1m on bottom liner (Giroud & Bonaparte, 1989)
leakage_rate = 0;
hole_count = 0;
contact_factor = 0;
hole_area = 0;
LeakageRate = zeros(indexes,1);
for ii = 1:indexes
    leakage_rate = 0;
    if CQA(ii)
        hole_count = 5;
        contact_factor = 0.21;
        hole_area = 1*10^-5;
    elseif MineralLinerPermeability(ii) == mineralLinerPermeability.GCL
        % GCL can assume good contact (Giroud, 1997)
        hole_count = 20;
        contact_factor = 0.21;
        hole_area = 5*10^-5;
    else
        hole_count = 20;
        contact_factor = 1.15;
        hole_area = 5*10^-5;
    end

    switch LinerKind(ii)
        case linerKind.MINERAL_LINER
            % Darcy's Law
            leakage_rate = MineralLinerPermeability(ii)*HydraulicGradient(ii)*AREA_M2;
        case linerKind.SINGLE_COMPOSITE
            % leakage through a single composite liner (Bonaparte et al., 1989)
            leakage_rate = hole_count*contact_factor*HEAD^0.9*hole_area^0.1*MineralLinerPermeability(ii)^0.74;
        case linerKind.DOUBLE_LINER_COMPOSITE
            % leakage through a single composite liner (Bonaparte et al., 1989)
            leakage_rate = hole_count*contact_factor*HEAD_BOTTOM_LINER^0.9*hole_area^0.1*MineralLinerPermeability(ii)^0.74;
        case linerKind.GEOMEMBRANE
            leakage_rate = 0.6*hole_area*sqrt(2*9.81*HEAD);
        case linerKind.DOUBLE_GEOMEMBRANE
            leakage_rate = 0.6*hole_area*sqrt(2*9.81*HEAD_BOTTOM_LINER);
        otherwise
            warning('Unexpected linerKind.')
    end
    LeakageRate(ii,1) = leakage_rate;
end

% calculate k (for LD compliance)
Permeability = LeakageRate./(HydraulicGradient.*AREA_M2);

LeakageRateDay = LeakageRate.*SECONDS_PER_DAY;
LifetimeLeakage = LeakageRate.*SECONDS_PER_YEAR*LANDFILL_LIFETIME_YRS;

leakage = table(LinerKind, CQA, MineralLinerPermeability, Permeability, LeakageRate, LeakageRateDay, LifetimeLeakage)

%% Part 4: Calculate cost

costs = struct('LLDPE', 8, 'TYRES', 10, 'PROTECTION_GEOTEXTILE', 4, 'LOW_PERMEABILITY_CLAY', 100, 'SEMI_LOW_PERMEABILITY_CLAY', 20, 'GCL', 15, 'CQA', 50, 'SEPARATOR_GEOTEXTILE', 2, 'DRAINAGE_GRAVEL', 50, 'RESTORATION_SOILS', 1);

% calculate basal lining cost
LiningMaterialCost = zeros(indexes,1);
current_material_cost = 0; % per m^2
for ii = 1:indexes
    current_material_cost = 0;
    if CQA(ii)
        current_material_cost = current_material_cost + costs.CQA;
    end

    switch LinerKind(ii)
        case linerKind.SINGLE_COMPOSITE
            current_material_cost = current_material_cost + costs.LLDPE;
        case linerKind.DOUBLE_LINER_COMPOSITE
            % tyre layer height of 10cm
            current_material_cost = current_material_cost + costs.LLDPE + costs.TYRES*0.1 + costs.LLDPE;
        case linerKind.GEOMEMBRANE
            current_material_cost = current_material_cost + costs.LLDPE;
        case linerKind.DOUBLE_GEOMEMBRANE
            current_material_cost = current_material_cost + costs.LLDPE + costs.TYRES*0.1 + costs.LLDPE;
    end

    switch MineralLinerPermeability(ii)
        case mineralLinerPermeability.GCL
            current_material_cost = current_material_cost + costs.GCL;
        case mineralLinerPermeability.LOW_PERMEABILITY_CLAY
            current_material_cost = current_material_cost + costs.LOW_PERMEABILITY_CLAY;
        case mineralLinerPermeability.SEMI_LOW_PERMEABILITY_CLAY
            current_material_cost = current_material_cost + costs.SEMI_LOW_PERMEABILITY_CLAY;
    end
    current_material_cost = current_material_cost*AREA_M2;
    LiningMaterialCost(ii) = current_material_cost;
end

LifetimeLeakageCost = LifetimeLeakage.*300;

drainage_cost = 170925;
DrainageCost = ones(indexes, 1);
DrainageCost = DrainageCost.*drainage_cost;

capping_cost = costs.RESTORATION_SOILS*1 + costs.RESTORATION_SOILS*1 + costs.SEPARATOR_GEOTEXTILE + costs.DRAINAGE_GRAVEL*0.5 + costs.PROTECTION_GEOTEXTILE + costs.LLDPE + costs.PROTECTION_GEOTEXTILE + costs.DRAINAGE_GRAVEL*0.3 + costs.SEPARATOR_GEOTEXTILE;
capping_cost = capping_cost*AREA_M2;
CappingCost = ones(indexes, 1);
CappingCost = CappingCost.*capping_cost;

cover_soil_cost = 51000;
CoverSoilCost = ones(indexes, 1);
CoverSoilCost = CoverSoilCost.*cover_soil_cost;

TotalCost = LiningMaterialCost+LifetimeLeakageCost+CoverSoilCost+DrainageCost+CappingCost;

cost = table(LinerKind, CQA, MineralLinerPermeability, Permeability, LifetimeLeakageCost, LiningMaterialCost, CoverSoilCost, DrainageCost, CappingCost, TotalCost)

%% Part 5: Calculate profit

% calculate income
% c.f. https://imperiallondon.sharepoint.com/:x:/r/sites/Landfilldesigngroup4-CI/_layouts/15/Doc.aspx?sourcedoc=%7B12524882-6197-42F4-A1A9-1AC2A37C1969%7D&file=Landfill%20dimensions.xlsx&action=default&mobileredirect=true&DefaultItemOpen=1
available_volume = 0;
AvailableVolume = zeros(indexes,1);
CoverSoilCost = zeros(indexes, 1);
for ii = 1:indexes
    available_volume = 0;
    if LinerKind(ii) == linerKind.DOUBLE_LINER_COMPOSITE
        % liner thickness is 0.6
        available_volume = 243552;
    else
        available_volume = 244375;
    end
    AvailableVolume(ii) = available_volume;
end

Income = AvailableVolume*7.5;
Profit = Income-TotalCost;

profit = table(LinerKind, CQA, MineralLinerPermeability, TotalCost, AvailableVolume, Income, Profit)

permeability_table = table(LinerKind, CQA, MineralLinerPermeability, Permeability)

figure(1);
color = [0 0 0 0];
marker = "";
legend_labels = string(zeros(indexes, 1));
legend_label_index = 0;
for ii = 1:indexes
    if Permeability(ii) > 10^-9
        continue
    end
    legend_label = "";
    switch LinerKind(ii) 
        case linerKind.MINERAL_LINER
            color = [0.8500 0.3250 0.0980];
            legend_label = "ML: ";
        case linerKind.SINGLE_COMPOSITE
            color = [0.4660 0.6740 0.1880];
            legend_label = "SC: ";
        case linerKind.DOUBLE_LINER_COMPOSITE
            color = [0 0.4470 0.7410];
            legend_label = "DLC: ";
        case linerKind.GEOMEMBRANE
            color = [0.4940 0.1840 0.5560];
            legend_label = "GM:";
        case linerKind.DOUBLE_GEOMEMBRANE
            color = [0.6350 0.0780 0.1840];
            legend_label = "DGM:";
    end
    switch MineralLinerPermeability(ii)
        case mineralLinerPermeability.LOW_PERMEABILITY_CLAY
            marker = "o";
            legend_label = legend_label + "LP Clay";
        case  mineralLinerPermeability.SEMI_LOW_PERMEABILITY_CLAY
            marker = "square";
            legend_label = legend_label + "Semi-LP Clay";
        case  mineralLinerPermeability.GCL
            marker = "diamond";
            legend_label = legend_label + "GCL";
        otherwise
            marker = "hexagram";
    end

    if CQA(ii)
        scatter(LeakageRateDay(ii), Profit(ii), 400, color, "filled", marker, LineWidth=5);
        legend_label = legend_label + " - CQA";
    else
        scatter(LeakageRateDay(ii), Profit(ii), 400, color, marker, LineWidth=5);
    end
    legend_label_index = legend_label_index + 1;
    legend_labels(legend_label_index) = legend_label;
    hold on;
end
title("Basal Liner Leakage vs. Profit");
xlabel("Leakage Rate (m^3/d)");
ylabel("Total Profit (£)");
legend(legend_labels, 'Location','westoutside');

ax = gca;
ax.XAxisLocation = "origin";
ax.YAxisLocation = "origin";
ax.YAxis.Exponent = 0;
ax.XAxis.Exponent = 0;
ax.FontSize = 30;

hold off;