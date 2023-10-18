clc; clear all; close all;
paths = getFilesPaths();
preFiles  = getPreFiles();%getAllFiles, getPreFiles, getPostFiles
postFiles  = getPostFiles();%getAllFiles, getPreFiles, getPostFiles
analysisPerZonesPath = 'F:\ForschungsProjekte\RFTC\MATLAB\DetectHFO\AnalysisPerZones_Exclusive_Boxplots\';

selDetector = 'MOSSDET'; 'Delphos'; 'MOSSDET'; 'MOSSDET_Depurated';
analysisPerZonesTablePath = strcat(analysisPerZonesPath, selDetector, '\');
features = {'OccRate', 'Power', 'Frequency'};
%features = {'OccRate', 'Power'};
%features = {'OccRate'};
normOptions = {'Normalized'}; %'Normalized'; % {'', 'Scaled', 'Normalized'};
freqBandConnList = { 'maxAllBands'};
thList = 0;
patsSelList = {'allPatients', 'improvedPatients', 'nonImprovedPatients'};%{'allPatients', 'improvedPatients', 'nonImprovedPatients'};
patsSelList = {'allPatients', 'improvedPatients', 'nonImprovedPatients'};
biomarkersList =  {'iesHFOVals', 'iaVals', 'allHFOVals', 'isolHFOVals', 'allRippleVals', 'iesRippleVals', 'isolRippleVals', 'allFR_Vals', 'iesFR_Vals', 'isolFR_Vals'};
biomarkersList =  {'iesHFOVals', 'allHFOVals', 'iaVals'};

params.analysisPerZonesTablePath = analysisPerZonesTablePath;
params.selDetector = selDetector;
params.connTh = 75;
params.eiTh = 75;
params.outcomeTh = 90;

normStrList = {'Normalized', 'noNorm'};

                
% Occ.Rate Delta
features = {'OccRate', 'Power', 'Frequency'};
testsCell = {'Feature', 'rftcSite', 'rftcConnected', 'rftcStructure', 'highEI', 'rftcLobe', 'rftcHemisphere'};
for normStrCell = normStrList
    normStr = normStrCell{1};
    predictionDataPath = strcat(paths.workspacePath, 'PredictionData\Prediction_Data_Structures_', normStr, '.mat');
    load(predictionDataPath)

    for featStr = features
    %     impPats_OccRateData = predictionData.('iesHFOVals').(normStr).('improvedPatients').('OccRate');
    %     nonImpPats_OccRateData = predictionData.('iesHFOVals').(normStr).('nonImprovedPatients').('OccRate');
    %     impPats_PowData = predictionData.('iesHFOVals').(normStr).('improvedPatients').('Power');
    %     nonImpPats_PowData = predictionData.('iesHFOVals').(normStr).('nonImprovedPatients').('Power');
    %     impPats_FreqData = predictionData.('iesHFOVals').(normStr).('improvedPatients').('Frequency');
    %     nonImpPats_FreqData = predictionData.('iesHFOVals').(normStr).('nonImprovedPatients').('Frequency');
    % 
    %     impPats_CombinedData = getCombinedMetric(impPats_OccRateData, impPats_PowData);
    %     nonImpPats_CombinedData = getCombinedMetric(nonImpPats_OccRateData, nonImpPats_PowData);
    % 
    %     impPatsData = impPats_PowData;
    %     nonImpPatsData = nonImpPats_PowData;

        impPatsData = predictionData.('iesHFOVals').(normStr).('improvedPatients').(featStr{1});
        nonImpPatsData = predictionData.('iesHFOVals').(normStr).('nonImprovedPatients').(featStr{1});

        rftcSite_p = ranksum(impPatsData.rftcDiff, nonImpPatsData.rftcDiff, 'tail', 'right');
        rftcConn_p = ranksum(impPatsData.rftcConn, nonImpPatsData.rftcConn, 'tail', 'right');
        rftcStruct_p = ranksum(impPatsData.rftcStruct, nonImpPatsData.rftcStruct, 'tail', 'right');
        highEI_p = ranksum(impPatsData.highEI, nonImpPatsData.highEI, 'tail', 'left');
        rftcLobe_p = ranksum(impPatsData.rftcLobe, nonImpPatsData.rftcLobe, 'tail', 'left');
        rftcHemis_p = ranksum(impPatsData.rftcHemis, nonImpPatsData.rftcHemis, 'tail', 'left');


        impPats_MeanDiff = cellfun(@mean, {impPatsData.rftcDiff, impPatsData.rftcConn, impPatsData.rftcStruct, impPatsData.highEI, impPatsData.rftcLobe, impPatsData.rftcHemis});
        nonImpPats_MeanDiff = cellfun(@mean, {nonImpPatsData.rftcDiff, nonImpPatsData.rftcConn, nonImpPatsData.rftcStruct, nonImpPatsData.highEI, nonImpPatsData.rftcLobe, nonImpPatsData.rftcHemis});
        meandDiffCell = {'rftcSite', 'rftcConnected', 'rftcStructure', 'highEI', 'rftcLobe', 'rftcHemisphere'};
        meandDiffCell = cat(1, meandDiffCell, num2cell(impPats_MeanDiff));
        meandDiffCell = cat(1, meandDiffCell, num2cell(nonImpPats_MeanDiff));

        testsCell = cat(1, testsCell, {strcat(featStr{1}, '_',normStr), rftcSite_p, rftcConn_p, rftcStruct_p, highEI_p, rftcLobe_p, rftcHemis_p})
    end
end

function combinedMetric = getCombinedMetric(occRateData, powData)
    combinedMetric.rftcDiff = occRateData.rftcDiff .* powData.rftcDiff;
    combinedMetric.rftcConn = occRateData.rftcConn .* powData.rftcConn;
    combinedMetric.rftcStruct = occRateData.rftcStruct .* powData.rftcStruct;
    combinedMetric.highEI = occRateData.highEI .* powData.highEI;
    combinedMetric.rftcLobe = occRateData.rftcLobe .* powData.rftcLobe;
    combinedMetric.rftcHemis = occRateData.rftcHemis .* powData.rftcHemis;
end
