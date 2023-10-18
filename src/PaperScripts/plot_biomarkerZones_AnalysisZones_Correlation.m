clc; clear all; close all;

paths = getFilesPaths();

patTablesPath = strcat(paths.workspacePath, 'CharacterizationTables\FeatureCharacterizationTables\');
groupTablesPath = strcat(paths.workspacePath, 'CharacterizationTables\GroupCharacterizationTablesAvg\');
analysisTablesPath = strcat(paths.workspacePath, 'AnalysisResults\Avg\');mkdir(analysisTablesPath);

biomarkersList =  {'HFO', 'iesHFO', 'IES'};
featuresList = {'rate', 'maxAmpl', 'power', 'variance'};

channSelStr = 'FlexK';
normStr = ''; %'', '_Normalized'};
zoneFormation = ''; % '', '_NonExclusive';
prctlTh = 75;
outcomeTh = 90;
corrLimitP = 0.001;

colors = {[0 0.4470 0.7410];...	
[0.8500 0.3250 0.0980];...
[0.9290 0.6940 0.1250];...
[0.4940 0.1840 0.5560];...
[0.4660 0.6740 0.1880];...
[0.3010 0.7450 0.9330];...
[0.6350 0.0780 0.1840]};

% Tables to read
spreadSheetNamePre = strcat(groupTablesPath, 'GroupCharacterization_', channSelStr, '_Pre', normStr,'.xls');
spreadSheetNamePost = strcat(groupTablesPath, 'GroupCharacterization_', channSelStr, '_Post', normStr,'.xls');

% Tables to read
tableBiomCorrelAllPatients_SpreadSheetName = strcat(analysisTablesPath, 'allPatsBiomarkerZoneCorrelation', channSelStr, '_', normStr,zoneFormation,'.xls');
  
maxMCC = 0;
minMCC = 0;
for bm = biomarkersList
    biomarker = bm{1};
    biomarkerZoneCorrT = readtable(tableBiomCorrelAllPatients_SpreadSheetName, 'Sheet', biomarker);
    % build correlation matrix
    corrMatrix = table2cell(biomarkerZoneCorrT(:,3:end));
    corrMatrix = cell2mat(corrMatrix);
    corrMatrixStruct.(biomarker) = corrMatrix;
    minBM = min(corrMatrix,[],'all');
    maxBM = max(corrMatrix,[],'all');
    if minBM < minMCC
        minMCC = minBM;
    end
    
    if maxBM > maxMCC
        maxMCC = maxBM;
    end
end

figure(1)
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w');
t = tiledlayout(3,3,'TileSpacing','Compact');
title(t,{'Correlation between'; 'Features and Analysis Zones'}, 'FontSize', 20)

for bmi = 1:length(biomarkersList)
    biomarker = biomarkersList{bmi};
    biomarkerZoneCorrT = readtable(tableBiomCorrelAllPatients_SpreadSheetName, 'Sheet', biomarker);
    % define Labels
    featNames = {'Occ.Rate', 'Amplitude', 'Variance', 'Power'};
    %featNames = biomarkerZoneCorrT.Feature;
    zoneNames = biomarkerZoneCorrT.Properties.VariableNames(3:end);

    % build correlation matrix
    corrMatrix = corrMatrixStruct.(biomarker)';
    
    nexttile
    h = heatmap(featNames, zoneNames, corrMatrix, 'Colormap',summer, 'CellLabelColor','black')
    h.Title = biomarker;
    h.FontSize = 16;
    %h.Colormap = summer;
    h.CellLabelFormat = '%.2f'
    caxis([minMCC maxMCC])
    %imagesc(corrMatrix);
%     set(gca, 'XTick', 1:size(corrMatrix,2)); % center x-axis ticks on bins
%     set(gca, 'YTick', 1:size(corrMatrix,1)); % center y-axis ticks on bins
%     set(gca, 'XTickLabel', zoneNames); % set x-axis labels
%     set(gca, 'YTickLabel', featNames); % set y-axis labels
    %xtickangle(30);
    %ax = gca;
    %ax.FontSize = 18; 
    %ax.FontWeight = 'bold';
    %title({biomarker; 'Feature-Zones Correlation'}, 'FontSize', 20); % set title
    %title(biomarker, 'FontSize', 20); % set title
    %colormap('parula'); % Choose jet or any other color scheme
%     for ri = 1:size(corrMatrix,1)
%         for ci = 1:size(corrMatrix,2)
%             label = num2str(corrMatrix(ri, ci), '%.2f');
%             text(ci-0.1, ri, label,'FontSize', 16, 'FontWeight', 'bold');
%         end
%     end   
    %if bmi == 3
%         c = colorbar; %
%         c.FontSize = 18;
%         c.FontWeight = 'bold';
%         %c.Limits = [minMCC maxMCC];
%         caxis([minMCC maxMCC])
%         c.Label.String = 'MCC';
    %end
    nexttile
    nexttile
end

boxchartPath = strcat('F:\ForschungsProjekte\RFTC\MATLAB_ver2\AnalysisResults\Avg\Plotheas\'); mkdir(boxchartPath)
boxchartFN = strcat(boxchartPath, 'biomarkerZones_AnalysisZones_Correlation_', zoneFormation);
saveas(gcf,boxchartFN);
boxchartFN = strcat(boxchartFN, '.jpg');
saveas(gcf,boxchartFN);
close(); 



