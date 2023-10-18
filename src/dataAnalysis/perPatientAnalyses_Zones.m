clc; clear all; close all;
paths = getFilesPaths();
preFiles  = getPreFiles();%getAllFiles, getPreFiles, getPostFiles
postFiles  = getPostFiles();%getAllFiles, getPreFiles, getPostFiles
analysisPerZonesPath = 'F:\ForschungsProjekte\RFTC\MATLAB\DetectHFO\AnalysisPerZones\Patientwise\';

selDetector = 'MOSSDET'; 'Delphos'; 'MOSSDET'; 'MOSSDET_Depurated';
features = {'OccRate', 'Power', 'Frequency'};
normOptions = {'', 'Scaled', 'Normalized'};
%normOptions = {'Normalized'};
freqBandConnList = {'delta', 'theta', 'alpha', 'beta', 'gamma', 'highGamma', 'ripple', 'fr', 'maxAllBands', 'meanAllBands'};
freqBandConnList = {'maxAllBands'};
thList = 0:2;
thList = 0;
if strcmp(selDetector, 'Delphos')
    thList = 0;
end
patsSelList = {'allPatients', 'improvedPatients', 'nonImprovedPatients'};
connTh = 75;
eiTh = 75;
outcomeTh = 49;

%%Start analyses
analysisPerZonesTablePath = strcat(analysisPerZonesPath, selDetector, '\'); mkdir(analysisPerZonesTablePath);

%get_Same_BrainParcel_Channels_Analysis_PrePost(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh);
%get_Opposite_BrainParcel_Channels_Analysis_PrePost(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh);
%get_Same_Hemisphere_Channels_Analysis_PrePost(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh);
%get_Opposite_Hemisphere_Channels_Analysis_PrePost(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh);

%All Pats
%get_Biomarker_EI_Correlation(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh);
%getAllChannelsPairedAnalysis_PrePost(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh);
%getRFTC_Channels_PairedAnalysis_PrePost(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh);
%get_AllMinusRFTC_Channels_PairedAnalysis_PrePost(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh);
%get_HighEI_Channels_Analysis_PrePost(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh);
get_RFTC_Connected_Channels_Analysis_PrePost(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh);
get_RFTC_NonConnected_Channels_Analysis_PrePost(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh);


function get_Biomarker_EI_Correlation(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh)
    analysisType = {'Biomarker_EI_Corr'};
    
    analysis = {};
    threshold = [];
    nrPats = [];
    nrChanns = [];
    iaVals = [];
    allHFOVals = []; iesHFOVals = []; isolHFOVals = [];
    allRippleVals = []; iesRippleVals = []; isolRippleVals = [];
    allFR_Vals = []; iesFR_Vals = []; isolFR_Vals = [];
            
    for th = thList
        for ni = 1:length(normOptions)
            for fi = 1:length(features)
                for psi = 1:length(patsSelList)
                    patsSel = patsSelList{psi};
                    feature = features{fi};
                    normalization = normOptions{ni};
                    analysisSubType = feature;

                    selDetectorFolder = selDetector;
                    if th > 0
                        selDetectorFolder = strcat('MOSSDET_Depurated\Th', num2str(th));
                    end
                    groupAnalysisTablesFilePath = strcat(paths.workspacePath, 'GroupAnalysis_ChannelCharacterizationTables\', selDetectorFolder, '\');
                    preGroupTableFN = strcat(groupAnalysisTablesFilePath, normalization, 'GroupAnalysis_ChannelCharacterization_Pre', selDetector, '.xls');
                    postGroupTableFN = strcat(groupAnalysisTablesFilePath, normalization, 'GroupAnalysis_ChannelCharacterization_Post', selDetector, '.xls');
                    groupTablePre = readtable(preGroupTableFN, 'Sheet', feature);
                    groupTablePost = readtable(postGroupTableFN, 'Sheet', feature);

                    %% read normalized EI
                    preGroupTableFN = strcat(groupAnalysisTablesFilePath, '', 'GroupAnalysis_ChannelCharacterization_Pre', selDetector, '.xls');
                    postGroupTableFN = strcat(groupAnalysisTablesFilePath, '', 'GroupAnalysis_ChannelCharacterization_Post', selDetector, '.xls');
                    groupTablePreEI = readtable(preGroupTableFN, 'Sheet', feature);
                    groupTablePostEI = readtable(postGroupTableFN, 'Sheet', feature);
                    groupTablePre.eiVals = groupTablePreEI.eiVals;
                    groupTablePost.eiVals = groupTablePostEI.eiVals;
                    
                    %% Select outcomes
                    if strcmp(patsSel, 'allPatients')
                        analysisSubType = strcat(analysisSubType, '_', patsSel);
                    elseif strcmp(patsSel, 'improvedPatients')
                        analysisSubType = strcat(analysisSubType, '_', patsSel);
                        tableDelIdxs = groupTablePre.outcome <= outcomeTh;
                        groupTablePre(tableDelIdxs, :) = [];
                        tableDelIdxs = groupTablePost.outcome <= outcomeTh;
                        groupTablePost(tableDelIdxs, :) = [];
                    elseif strcmp(patsSel, 'nonImprovedPatients')
                        analysisSubType = strcat(analysisSubType, '_', patsSel);
                        tableDelIdxs = groupTablePre.outcome > outcomeTh;
                        groupTablePre(tableDelIdxs, :) = [];
                        tableDelIdxs = groupTablePost.outcome > outcomeTh;
                        groupTablePost(tableDelIdxs, :) = [];
                    end
                    
                    analysisSubType = strcat(analysisSubType, '_', normalization);
                    %%

                    analysis = cat(1, analysis, analysisSubType);
                    threshold = cat(1, threshold, th);
                    nrPats = cat(1, nrPats, length(unique(groupTablePost.patName)));
                    nrChanns = cat(1, nrChanns, length(groupTablePost.patName));
                    
                    patNameCol = groupTablePre.patName;
                    patNames = unique(patNameCol);
                    iaValsSum = 0;
                    allHFOValsSum = 0; iesHFOValsSum = 0; isolHFOValsSum = 0;
                    allRippleValsSum = 0; iesRippleValsSum = 0; isolRippleValsSum = 0;
                    allFR_ValsSum = 0; iesFR_ValsSum = 0; isolFR_ValsSum = 0;
                    nrPatsForAvg = length(patNames);
                    for pi = 1:nrPatsForAvg
                        patName = patNames{pi};
                        patSelIdx = ismember(patNameCol, patName);
                        perPatTablePre = groupTablePre(patSelIdx, :);
                        
                        iaValsSum = iaValsSum + getSigCorr(perPatTablePre.iaVals, perPatTablePre.eiVals, 'Type','Spearman');
                        allHFOValsSum = allHFOValsSum + getSigCorr(groupTablePre.allHFOVals, groupTablePre.eiVals, 'Type','Spearman');
                        iesHFOValsSum = iesHFOValsSum + getSigCorr(groupTablePre.iesHFOVals, groupTablePre.eiVals, 'Type','Spearman');
                        isolHFOValsSum = isolHFOValsSum + getSigCorr(groupTablePre.isolHFOVals, groupTablePre.eiVals, 'Type','Spearman');

                        allRippleValsSum = allRippleValsSum + getSigCorr(groupTablePre.allRippleVals, groupTablePre.eiVals, 'Type','Spearman');
                        iesRippleValsSum = iesRippleValsSum + getSigCorr(groupTablePre.iesRippleVals, groupTablePre.eiVals, 'Type','Spearman');
                        isolRippleValsSum = isolRippleValsSum + getSigCorr(groupTablePre.isolRippleVals, groupTablePre.eiVals, 'Type','Spearman');

                        allFR_ValsSum = allFR_ValsSum + getSigCorr(groupTablePre.allFR_Vals, groupTablePre.eiVals, 'Type','Spearman');
                        iesFR_ValsSum = iesFR_ValsSum + getSigCorr(groupTablePre.iesFR_Vals, groupTablePre.eiVals, 'Type','Spearman');
                        isolFR_ValsSum = isolFR_ValsSum + getSigCorr(groupTablePre.isolFR_Vals, groupTablePre.eiVals, 'Type','Spearman');
                    end
        
                    iaVals = cat(1, iaVals, iaValsSum / nrPatsForAvg*100);
                    allHFOVals = cat(1, allHFOVals, allHFOValsSum / nrPatsForAvg*100);
                    iesHFOVals = cat(1, iesHFOVals, iesHFOValsSum / nrPatsForAvg*100);
                    isolHFOVals = cat(1, isolHFOVals, isolHFOValsSum / nrPatsForAvg*100);
                    allRippleVals = cat(1, allRippleVals, allRippleValsSum / nrPatsForAvg*100);
                    iesRippleVals = cat(1, iesRippleVals, iesRippleValsSum / nrPatsForAvg*100);
                    isolRippleVals = cat(1, isolRippleVals, isolRippleValsSum / nrPatsForAvg*100);
                    allFR_Vals = cat(1, allFR_Vals, allFR_ValsSum / nrPatsForAvg*100);
                    iesFR_Vals = cat(1, iesFR_Vals, iesFR_ValsSum / nrPatsForAvg*100);
                    isolFR_Vals = cat(1, isolFR_Vals, isolFR_ValsSum / nrPatsForAvg*100);
                end
            end
        end
        iaVals = round(iaVals);
        allHFOVals = round(allHFOVals); iesHFOVals = round(iesHFOVals); isolHFOVals = round(isolHFOVals);
        allRippleVals = round(allRippleVals); iesRippleVals = round(iesRippleVals); isolRippleVals = round(isolRippleVals);
        allFR_Vals = round(allFR_Vals); iesFR_Vals = round(iesFR_Vals); isolFR_Vals = round(isolFR_Vals);
                    
        zonesAnalysisTable = table(analysis, threshold, nrPats, nrChanns, iaVals, allHFOVals, iesHFOVals, isolHFOVals, allRippleVals, iesRippleVals, isolRippleVals, allFR_Vals, iesFR_Vals, isolFR_Vals); 
        spaceT = zonesAnalysisTable(1,:);
        for ti = 1:length(spaceT.Properties.VariableNames)
            spaceT.(spaceT.Properties.VariableNames{ti}) = nan(1,1);
        end
        zonesAnalysisTableSpaces = [];
        for si = 0:(size(zonesAnalysisTable,1)/3)-1
            idx = 3*si;
            zonesAnalysisTableSpaces = cat(1, zonesAnalysisTableSpaces, zonesAnalysisTable(idx+1:idx+3,:), spaceT);
        end
        
        newSelDetector = selDetector;
        if th > 0
            analysisPerZonesTablePath(end) = [];
            analysisPerZonesTablePath = strcat(analysisPerZonesTablePath, '_Depurated\', 'Th', num2str(th), '\');mkdir(analysisPerZonesTablePath);
            newSelDetector = strcat(selDetector, '_Th', num2str(th));
        end
        analysisPerZonesTableFN = strcat(analysisPerZonesTablePath, 'Biomarker_EI_Correlation_Analysis_', newSelDetector,'.xls');
        writetable(zonesAnalysisTableSpaces, analysisPerZonesTableFN, 'Sheet', analysisType{1});
    end
end

function rho = getSigCorr(vecA, vecB, typeSet, corrType)
    [rho,pval] = corr(vecA, vecB, typeSet,corrType);
    if not(pval < 0.05)
        rho = 0;
    end
end


function getAllChannelsPairedAnalysis_PrePost(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh)
    analysisType = {'allChannels'};
    
    analysis = {};
    threshold = [];
    nrPats = [];
    nrChanns = [];
    iaVals = [];
    allHFOVals = []; iesHFOVals = []; isolHFOVals = [];
    allRippleVals = []; iesRippleVals = []; isolRippleVals = [];
    allFR_Vals = []; iesFR_Vals = []; isolFR_Vals = [];
            
    for th = thList
        for ni = 1:length(normOptions)
            for fi = 1:length(features)
                for psi = 1:length(patsSelList)
                    patsSel = patsSelList{psi};
                    feature = features{fi};
                    normalization = normOptions{ni};
                    analysisSubType = feature;
                    
                    selDetectorFolder = selDetector;
                    if th > 0
                        selDetectorFolder = strcat('MOSSDET_Depurated\Th', num2str(th));
                    end
                    groupAnalysisTablesFilePath = strcat(paths.workspacePath, 'GroupAnalysis_ChannelCharacterizationTables\', selDetectorFolder, '\');
                    preGroupTableFN = strcat(groupAnalysisTablesFilePath, normalization, 'GroupAnalysis_ChannelCharacterization_Pre', selDetector, '.xls');
                    postGroupTableFN = strcat(groupAnalysisTablesFilePath, normalization, 'GroupAnalysis_ChannelCharacterization_Post', selDetector, '.xls');
                    groupTablePre = readtable(preGroupTableFN, 'Sheet', feature);
                    groupTablePost = readtable(postGroupTableFN, 'Sheet', feature);

                    %% Select outcomes
                    if strcmp(patsSel, 'allPatients')
                        analysisSubType = strcat(analysisSubType, '_', patsSel);
                    elseif strcmp(patsSel, 'improvedPatients')
                        analysisSubType = strcat(analysisSubType, '_', patsSel);
                        tableDelIdxs = groupTablePre.outcome <= outcomeTh;
                        groupTablePre(tableDelIdxs, :) = [];
                        tableDelIdxs = groupTablePost.outcome <= outcomeTh;
                        groupTablePost(tableDelIdxs, :) = [];
                    elseif strcmp(patsSel, 'nonImprovedPatients')
                        analysisSubType = strcat(analysisSubType, '_', patsSel);
                        tableDelIdxs = groupTablePre.outcome > outcomeTh;
                        groupTablePre(tableDelIdxs, :) = [];
                        tableDelIdxs = groupTablePost.outcome > outcomeTh;
                        groupTablePost(tableDelIdxs, :) = [];
                    end
                    
                    analysisSubType = strcat(analysisSubType, '_', normalization);
                    %%

                    analysis = cat(1, analysis, analysisSubType);
                    threshold = cat(1, threshold, th);
                    nrPats = cat(1, nrPats, length(unique(groupTablePost.patName)));
                    nrChanns = cat(1, nrChanns, length(groupTablePost.patName));
                    
                    patNameCol = groupTablePre.patName;
                    patNames = unique(patNameCol);
                    iaValsSum = 0;
                    allHFOValsSum = 0; iesHFOValsSum = 0; isolHFOValsSum = 0;
                    allRippleValsSum = 0; iesRippleValsSum = 0; isolRippleValsSum = 0;
                    allFR_ValsSum = 0; iesFR_ValsSum = 0; isolFR_ValsSum = 0;
                    nrPatsForAvg = length(patNames);
                    pValLim = 0.05;
                    for pi = 1:nrPatsForAvg
                        patName = patNames{pi};
                        patSelIdx = ismember(patNameCol, patName);
                        perPatTablePre = groupTablePre(patSelIdx, :);
                        perPatTablePost = groupTablePost(patSelIdx, :);
                            
                        iaValsSum = iaValsSum + (signrank(perPatTablePre.iaVals, perPatTablePost.iaVals, 'tail', 'right') < pValLim);
                        allHFOValsSum = allHFOValsSum + (signrank(perPatTablePre.allHFOVals, perPatTablePost.allHFOVals, 'tail', 'right') < pValLim);
                        iesHFOValsSum = iesHFOValsSum + (signrank(perPatTablePre.iesHFOVals, perPatTablePost.iesHFOVals, 'tail', 'right') < pValLim);
                        isolHFOValsSum = isolHFOValsSum + (signrank(perPatTablePre.isolHFOVals, perPatTablePost.isolHFOVals, 'tail', 'right') < pValLim);

                        allRippleValsSum = allRippleValsSum + (signrank(perPatTablePre.allRippleVals, perPatTablePost.allRippleVals, 'tail', 'right') < pValLim);
                        iesRippleValsSum = iesRippleValsSum + (signrank(perPatTablePre.iesRippleVals, perPatTablePost.iesRippleVals, 'tail', 'right') < pValLim);
                        isolRippleValsSum = isolRippleValsSum + (signrank(perPatTablePre.isolRippleVals, perPatTablePost.isolRippleVals, 'tail', 'right') < pValLim);

                        allFR_ValsSum = allFR_ValsSum + (signrank(perPatTablePre.allFR_Vals, perPatTablePost.allFR_Vals, 'tail', 'right') < pValLim);
                        iesFR_ValsSum = iesFR_ValsSum + (signrank(perPatTablePre.iesFR_Vals, perPatTablePost.iesFR_Vals, 'tail', 'right') < pValLim);
                        isolFR_ValsSum = isolFR_ValsSum + (signrank(perPatTablePre.isolFR_Vals, perPatTablePost.isolFR_Vals, 'tail', 'right') < pValLim);
                    end

                    iaVals = cat(1, iaVals, iaValsSum / nrPatsForAvg*100);
                    allHFOVals = cat(1, allHFOVals, allHFOValsSum / nrPatsForAvg*100);
                    iesHFOVals = cat(1, iesHFOVals, iesHFOValsSum / nrPatsForAvg*100);
                    isolHFOVals = cat(1, isolHFOVals, isolHFOValsSum / nrPatsForAvg*100);
                    allRippleVals = cat(1, allRippleVals, allRippleValsSum / nrPatsForAvg*100);
                    iesRippleVals = cat(1, iesRippleVals, iesRippleValsSum / nrPatsForAvg*100);
                    isolRippleVals = cat(1, isolRippleVals, isolRippleValsSum / nrPatsForAvg*100);
                    allFR_Vals = cat(1, allFR_Vals, allFR_ValsSum / nrPatsForAvg*100);
                    iesFR_Vals = cat(1, iesFR_Vals, iesFR_ValsSum / nrPatsForAvg*100);
                    isolFR_Vals = cat(1, isolFR_Vals, isolFR_ValsSum / nrPatsForAvg*100);
                end
            end
        end
        iaVals = round(iaVals);
        allHFOVals = round(allHFOVals); iesHFOVals = round(iesHFOVals); isolHFOVals = round(isolHFOVals);
        allRippleVals = round(allRippleVals); iesRippleVals = round(iesRippleVals); isolRippleVals = round(isolRippleVals);
        allFR_Vals = round(allFR_Vals); iesFR_Vals = round(iesFR_Vals); isolFR_Vals = round(isolFR_Vals);
        
        zonesAnalysisTable = table(analysis, threshold, nrPats, nrChanns, iaVals, allHFOVals, iesHFOVals, isolHFOVals, allRippleVals, iesRippleVals, isolRippleVals, allFR_Vals, iesFR_Vals, isolFR_Vals);
        spaceT = zonesAnalysisTable(1,:);
        for ti = 1:length(spaceT.Properties.VariableNames)
            spaceT.(spaceT.Properties.VariableNames{ti}) = nan(1,1);
        end
        zonesAnalysisTableSpaces = [];
        for si = 0:(size(zonesAnalysisTable,1)/3)-1
            idx = 3*si;
            zonesAnalysisTableSpaces = cat(1, zonesAnalysisTableSpaces, zonesAnalysisTable(idx+1:idx+3,:), spaceT);
        end
        newSelDetector = selDetector;
        if th > 0
            analysisPerZonesTablePath(end) = [];
            analysisPerZonesTablePath = strcat(analysisPerZonesTablePath, '_Depurated\', 'Th', num2str(th), '\');mkdir(analysisPerZonesTablePath);
            newSelDetector = strcat(selDetector, '_Th', num2str(th));
        end
        analysisPerZonesTableFN = strcat(analysisPerZonesTablePath, 'AnalysisPerZones_', newSelDetector,'.xls');
        writetable(zonesAnalysisTableSpaces, analysisPerZonesTableFN, 'Sheet', analysisType{1});
    end
end

function getRFTC_Channels_PairedAnalysis_PrePost(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh)
    analysisType = {'RFTC_Channels_PrePost'};
    
    analysis = {};
    threshold = [];
    nrPats = [];
    nrChanns = [];
    iaVals = [];
    allHFOVals = []; iesHFOVals = []; isolHFOVals = [];
    allRippleVals = []; iesRippleVals = []; isolRippleVals = [];
    allFR_Vals = []; iesFR_Vals = []; isolFR_Vals = [];
            
    for th = thList
        for ni = 1:length(normOptions)
            for fi = 1:length(features)
                for psi = 1:length(patsSelList)
                    patsSel = patsSelList{psi};
                    feature = features{fi};
                    normalization = normOptions{ni};
                    analysisSubType = feature;

                    selDetectorFolder = selDetector;
                    if th > 0
                        selDetectorFolder = strcat('MOSSDET_Depurated\Th', num2str(th));
                    end
                    groupAnalysisTablesFilePath = strcat(paths.workspacePath, 'GroupAnalysis_ChannelCharacterizationTables\', selDetectorFolder, '\');
                    preGroupTableFN = strcat(groupAnalysisTablesFilePath, normalization, 'GroupAnalysis_ChannelCharacterization_Pre', selDetector, '.xls');
                    postGroupTableFN = strcat(groupAnalysisTablesFilePath, normalization, 'GroupAnalysis_ChannelCharacterization_Post', selDetector, '.xls');
                    groupTablePre = readtable(preGroupTableFN, 'Sheet', feature);
                    groupTablePost = readtable(postGroupTableFN, 'Sheet', feature);

                    %% Select channels to compare
                    tableDelIdxs = not(groupTablePre.rftcVals == 1);
                    groupTablePre(tableDelIdxs, :) = [];
                    tableDelIdxs = not(groupTablePost.rftcVals == 1);
                    groupTablePost(tableDelIdxs, :) = [];
                    
                    %% Select outcomes
                    if strcmp(patsSel, 'allPatients')
                        analysisSubType = strcat(analysisSubType, '_', patsSel);
                    elseif strcmp(patsSel, 'improvedPatients')
                        analysisSubType = strcat(analysisSubType, '_', patsSel);
                        tableDelIdxs = groupTablePre.outcome <= outcomeTh;
                        groupTablePre(tableDelIdxs, :) = [];
                        tableDelIdxs = groupTablePost.outcome <= outcomeTh;
                        groupTablePost(tableDelIdxs, :) = [];
                    elseif strcmp(patsSel, 'nonImprovedPatients')
                        analysisSubType = strcat(analysisSubType, '_', patsSel);
                        tableDelIdxs = groupTablePre.outcome > outcomeTh;
                        groupTablePre(tableDelIdxs, :) = [];
                        tableDelIdxs = groupTablePost.outcome > outcomeTh;
                        groupTablePost(tableDelIdxs, :) = [];
                    end
                    
                    analysisSubType = strcat(analysisSubType, '_', normalization);
                    %%

                    analysis = cat(1, analysis, analysisSubType);
                    threshold = cat(1, threshold, th);
                    nrPats = cat(1, nrPats, length(unique(groupTablePost.patName)));
                    nrChanns = cat(1, nrChanns, length(groupTablePost.patName));

                    patNameCol = groupTablePre.patName;
                    patNames = unique(patNameCol);
                    iaValsSum = 0;
                    allHFOValsSum = 0; iesHFOValsSum = 0; isolHFOValsSum = 0;
                    allRippleValsSum = 0; iesRippleValsSum = 0; isolRippleValsSum = 0;
                    allFR_ValsSum = 0; iesFR_ValsSum = 0; isolFR_ValsSum = 0;
                    nrPatsForAvg = length(patNames);
                    pValLim = 0.05;
                    for pi = 1:nrPatsForAvg
                        patName = patNames{pi};
                        patSelIdx = ismember(patNameCol, patName);
                        perPatTablePre = groupTablePre(patSelIdx, :);
                        perPatTablePost = groupTablePost(patSelIdx, :);
                            
                        iaValsSum = iaValsSum + (signrank(perPatTablePre.iaVals, perPatTablePost.iaVals, 'tail', 'right') < pValLim);
                        allHFOValsSum = allHFOValsSum + (signrank(perPatTablePre.allHFOVals, perPatTablePost.allHFOVals, 'tail', 'right') < pValLim);
                        iesHFOValsSum = iesHFOValsSum + (signrank(perPatTablePre.iesHFOVals, perPatTablePost.iesHFOVals, 'tail', 'right') < pValLim);
                        isolHFOValsSum = isolHFOValsSum + (signrank(perPatTablePre.isolHFOVals, perPatTablePost.isolHFOVals, 'tail', 'right') < pValLim);

                        allRippleValsSum = allRippleValsSum + (signrank(perPatTablePre.allRippleVals, perPatTablePost.allRippleVals, 'tail', 'right') < pValLim);
                        iesRippleValsSum = iesRippleValsSum + (signrank(perPatTablePre.iesRippleVals, perPatTablePost.iesRippleVals, 'tail', 'right') < pValLim);
                        isolRippleValsSum = isolRippleValsSum + (signrank(perPatTablePre.isolRippleVals, perPatTablePost.isolRippleVals, 'tail', 'right') < pValLim);

                        allFR_ValsSum = allFR_ValsSum + (signrank(perPatTablePre.allFR_Vals, perPatTablePost.allFR_Vals, 'tail', 'right') < pValLim);
                        iesFR_ValsSum = iesFR_ValsSum + (signrank(perPatTablePre.iesFR_Vals, perPatTablePost.iesFR_Vals, 'tail', 'right') < pValLim);
                        isolFR_ValsSum = isolFR_ValsSum + (signrank(perPatTablePre.isolFR_Vals, perPatTablePost.isolFR_Vals, 'tail', 'right') < pValLim);
                    end

                    iaVals = cat(1, iaVals, iaValsSum / nrPatsForAvg*100);
                    allHFOVals = cat(1, allHFOVals, allHFOValsSum / nrPatsForAvg*100);
                    iesHFOVals = cat(1, iesHFOVals, iesHFOValsSum / nrPatsForAvg*100);
                    isolHFOVals = cat(1, isolHFOVals, isolHFOValsSum / nrPatsForAvg*100);
                    allRippleVals = cat(1, allRippleVals, allRippleValsSum / nrPatsForAvg*100);
                    iesRippleVals = cat(1, iesRippleVals, iesRippleValsSum / nrPatsForAvg*100);
                    isolRippleVals = cat(1, isolRippleVals, isolRippleValsSum / nrPatsForAvg*100);
                    allFR_Vals = cat(1, allFR_Vals, allFR_ValsSum / nrPatsForAvg*100);
                    iesFR_Vals = cat(1, iesFR_Vals, iesFR_ValsSum / nrPatsForAvg*100);
                    isolFR_Vals = cat(1, isolFR_Vals, isolFR_ValsSum / nrPatsForAvg*100);
                end
            end
        end
        iaVals = round(iaVals);
        allHFOVals = round(allHFOVals); iesHFOVals = round(iesHFOVals); isolHFOVals = round(isolHFOVals);
        allRippleVals = round(allRippleVals); iesRippleVals = round(iesRippleVals); isolRippleVals = round(isolRippleVals);
        allFR_Vals = round(allFR_Vals); iesFR_Vals = round(iesFR_Vals); isolFR_Vals = round(isolFR_Vals);
        
        zonesAnalysisTable = table(analysis, threshold, nrPats, nrChanns, iaVals, allHFOVals, iesHFOVals, isolHFOVals, allRippleVals, iesRippleVals, isolRippleVals, allFR_Vals, iesFR_Vals, isolFR_Vals); 
        spaceT = zonesAnalysisTable(1,:);
        for ti = 1:length(spaceT.Properties.VariableNames)
            spaceT.(spaceT.Properties.VariableNames{ti}) = nan(1,1);
        end
        zonesAnalysisTableSpaces = [];
        for si = 0:(size(zonesAnalysisTable,1)/3)-1
            idx = 3*si;
            zonesAnalysisTableSpaces = cat(1, zonesAnalysisTableSpaces, zonesAnalysisTable(idx+1:idx+3,:), spaceT);
        end
        
        newSelDetector = selDetector;
        if th > 0
            analysisPerZonesTablePath(end) = [];
            analysisPerZonesTablePath = strcat(analysisPerZonesTablePath, '_Depurated\', 'Th', num2str(th), '\');mkdir(analysisPerZonesTablePath);
            newSelDetector = strcat(selDetector, '_Th', num2str(th));
        end
        analysisPerZonesTableFN = strcat(analysisPerZonesTablePath, 'AnalysisPerZones_', newSelDetector,'.xls');
        writetable(zonesAnalysisTableSpaces, analysisPerZonesTableFN, 'Sheet', analysisType{1});
    end
end

function get_AllMinusRFTC_Channels_PairedAnalysis_PrePost(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh)
    analysisType = {'AllMinusRFTC_Channels'};
    
    analysis = {};
    threshold = [];
    nrPats = [];
    nrChanns = [];
    iaVals = [];
    allHFOVals = []; iesHFOVals = []; isolHFOVals = [];
    allRippleVals = []; iesRippleVals = []; isolRippleVals = [];
    allFR_Vals = []; iesFR_Vals = []; isolFR_Vals = [];
            
    for th = thList
        for ni = 1:length(normOptions)
            for fi = 1:length(features)
                for psi = 1:length(patsSelList)
                     patsSel = patsSelList{psi};
                    feature = features{fi};
                    normalization = normOptions{ni};
                    analysisSubType = feature;

                    selDetectorFolder = selDetector;
                    if th > 0
                        selDetectorFolder = strcat('MOSSDET_Depurated\Th', num2str(th));
                    end
                    groupAnalysisTablesFilePath = strcat(paths.workspacePath, 'GroupAnalysis_ChannelCharacterizationTables\', selDetectorFolder, '\');
                    preGroupTableFN = strcat(groupAnalysisTablesFilePath, normalization, 'GroupAnalysis_ChannelCharacterization_Pre', selDetector, '.xls');
                    postGroupTableFN = strcat(groupAnalysisTablesFilePath, normalization, 'GroupAnalysis_ChannelCharacterization_Post', selDetector, '.xls');
                    groupTablePre = readtable(preGroupTableFN, 'Sheet', feature);
                    groupTablePost = readtable(postGroupTableFN, 'Sheet', feature);

                    %% Select channels to compare
                    tableDelIdxs = groupTablePre.rftcVals == 1;
                    groupTablePre(tableDelIdxs, :) = [];
                    tableDelIdxs = groupTablePost.rftcVals == 1;
                    groupTablePost(tableDelIdxs, :) = [];
                    
                    %% Select outcomes
                    if strcmp(patsSel, 'allPatients')
                        analysisSubType = strcat(analysisSubType, '_', patsSel);
                    elseif strcmp(patsSel, 'improvedPatients')
                        analysisSubType = strcat(analysisSubType, '_', patsSel);
                        tableDelIdxs = groupTablePre.outcome <= outcomeTh;
                        groupTablePre(tableDelIdxs, :) = [];
                        tableDelIdxs = groupTablePost.outcome <= outcomeTh;
                        groupTablePost(tableDelIdxs, :) = [];
                    elseif strcmp(patsSel, 'nonImprovedPatients')
                        analysisSubType = strcat(analysisSubType, '_', patsSel);
                        tableDelIdxs = groupTablePre.outcome > outcomeTh;
                        groupTablePre(tableDelIdxs, :) = [];
                        tableDelIdxs = groupTablePost.outcome > outcomeTh;
                        groupTablePost(tableDelIdxs, :) = [];
                    end
                    
                    analysisSubType = strcat(analysisSubType, '_', normalization);
                    %%

                    analysis = cat(1, analysis, analysisSubType);
                    threshold = cat(1, threshold, th);
                    nrPats = cat(1, nrPats, length(unique(groupTablePost.patName)));
                    nrChanns = cat(1, nrChanns, length(groupTablePost.patName));
                    
                                        patNameCol = groupTablePre.patName;
                    patNames = unique(patNameCol);
                    iaValsSum = 0;
                    allHFOValsSum = 0; iesHFOValsSum = 0; isolHFOValsSum = 0;
                    allRippleValsSum = 0; iesRippleValsSum = 0; isolRippleValsSum = 0;
                    allFR_ValsSum = 0; iesFR_ValsSum = 0; isolFR_ValsSum = 0;
                    nrPatsForAvg = length(patNames);
                    pValLim = 0.05;
                    for pi = 1:nrPatsForAvg
                        patName = patNames{pi};
                        patSelIdx = ismember(patNameCol, patName);
                        perPatTablePre = groupTablePre(patSelIdx, :);
                        perPatTablePost = groupTablePost(patSelIdx, :);
                            
                        iaValsSum = iaValsSum + (signrank(perPatTablePre.iaVals, perPatTablePost.iaVals, 'tail', 'right') < pValLim);
                        allHFOValsSum = allHFOValsSum + (signrank(perPatTablePre.allHFOVals, perPatTablePost.allHFOVals, 'tail', 'right') < pValLim);
                        iesHFOValsSum = iesHFOValsSum + (signrank(perPatTablePre.iesHFOVals, perPatTablePost.iesHFOVals, 'tail', 'right') < pValLim);
                        isolHFOValsSum = isolHFOValsSum + (signrank(perPatTablePre.isolHFOVals, perPatTablePost.isolHFOVals, 'tail', 'right') < pValLim);

                        allRippleValsSum = allRippleValsSum + (signrank(perPatTablePre.allRippleVals, perPatTablePost.allRippleVals, 'tail', 'right') < pValLim);
                        iesRippleValsSum = iesRippleValsSum + (signrank(perPatTablePre.iesRippleVals, perPatTablePost.iesRippleVals, 'tail', 'right') < pValLim);
                        isolRippleValsSum = isolRippleValsSum + (signrank(perPatTablePre.isolRippleVals, perPatTablePost.isolRippleVals, 'tail', 'right') < pValLim);

                        allFR_ValsSum = allFR_ValsSum + (signrank(perPatTablePre.allFR_Vals, perPatTablePost.allFR_Vals, 'tail', 'right') < pValLim);
                        iesFR_ValsSum = iesFR_ValsSum + (signrank(perPatTablePre.iesFR_Vals, perPatTablePost.iesFR_Vals, 'tail', 'right') < pValLim);
                        isolFR_ValsSum = isolFR_ValsSum + (signrank(perPatTablePre.isolFR_Vals, perPatTablePost.isolFR_Vals, 'tail', 'right') < pValLim);
                    end

                    iaVals = cat(1, iaVals, iaValsSum / nrPatsForAvg*100);
                    allHFOVals = cat(1, allHFOVals, allHFOValsSum / nrPatsForAvg*100);
                    iesHFOVals = cat(1, iesHFOVals, iesHFOValsSum / nrPatsForAvg*100);
                    isolHFOVals = cat(1, isolHFOVals, isolHFOValsSum / nrPatsForAvg*100);
                    allRippleVals = cat(1, allRippleVals, allRippleValsSum / nrPatsForAvg*100);
                    iesRippleVals = cat(1, iesRippleVals, iesRippleValsSum / nrPatsForAvg*100);
                    isolRippleVals = cat(1, isolRippleVals, isolRippleValsSum / nrPatsForAvg*100);
                    allFR_Vals = cat(1, allFR_Vals, allFR_ValsSum / nrPatsForAvg*100);
                    iesFR_Vals = cat(1, iesFR_Vals, iesFR_ValsSum / nrPatsForAvg*100);
                    isolFR_Vals = cat(1, isolFR_Vals, isolFR_ValsSum / nrPatsForAvg*100);
                end
            end
        end
        iaVals = round(iaVals);
        allHFOVals = round(allHFOVals); iesHFOVals = round(iesHFOVals); isolHFOVals = round(isolHFOVals);
        allRippleVals = round(allRippleVals); iesRippleVals = round(iesRippleVals); isolRippleVals = round(isolRippleVals);
        allFR_Vals = round(allFR_Vals); iesFR_Vals = round(iesFR_Vals); isolFR_Vals = round(isolFR_Vals);
        
        zonesAnalysisTable = table(analysis, threshold, nrPats, nrChanns, iaVals, allHFOVals, iesHFOVals, isolHFOVals, allRippleVals, iesRippleVals, isolRippleVals, allFR_Vals, iesFR_Vals, isolFR_Vals); 
        spaceT = zonesAnalysisTable(1,:);
        for ti = 1:length(spaceT.Properties.VariableNames)
            spaceT.(spaceT.Properties.VariableNames{ti}) = nan(1,1);
        end
        zonesAnalysisTableSpaces = [];
        for si = 0:(size(zonesAnalysisTable,1)/3)-1
            idx = 3*si;
            zonesAnalysisTableSpaces = cat(1, zonesAnalysisTableSpaces, zonesAnalysisTable(idx+1:idx+3,:), spaceT);
        end
        
        newSelDetector = selDetector;
        if th > 0
            analysisPerZonesTablePath(end) = [];
            analysisPerZonesTablePath = strcat(analysisPerZonesTablePath, '_Depurated\', 'Th', num2str(th), '\');mkdir(analysisPerZonesTablePath);
            selDetector = strcat(selDetector, '_Th', num2str(th));
        end
        analysisPerZonesTableFN = strcat(analysisPerZonesTablePath, 'AnalysisPerZones_', newSelDetector,'.xls');
        writetable(zonesAnalysisTableSpaces, analysisPerZonesTableFN, 'Sheet', analysisType{1});
    end
end

function get_HighEI_Channels_Analysis_PrePost(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh)
    analysisType = {'HighEI_Channels'};
    
    analysis = {};
    threshold = [];
    nrPats = [];
    nrChanns = [];
    iaVals = [];
    allHFOVals = []; iesHFOVals = []; isolHFOVals = [];
    allRippleVals = []; iesRippleVals = []; isolRippleVals = [];
    allFR_Vals = []; iesFR_Vals = []; isolFR_Vals = [];
            
    for th = thList
        for ni = 1:length(normOptions)
            for fi = 1:length(features)
                for psi = 1:length(patsSelList)
                    patsSel = patsSelList{psi};
                    feature = features{fi};
                    normalization = normOptions{ni};
                    analysisSubType = feature;

                    selDetectorFolder = selDetector;
                    if th > 0
                        selDetectorFolder = strcat('MOSSDET_Depurated\Th', num2str(th));
                    end
                    groupAnalysisTablesFilePath = strcat(paths.workspacePath, 'GroupAnalysis_ChannelCharacterizationTables\', selDetectorFolder, '\');
                    preGroupTableFN = strcat(groupAnalysisTablesFilePath, normalization, 'GroupAnalysis_ChannelCharacterization_Pre', selDetector, '.xls');
                    postGroupTableFN = strcat(groupAnalysisTablesFilePath, normalization, 'GroupAnalysis_ChannelCharacterization_Post', selDetector, '.xls');
                    groupTablePre = readtable(preGroupTableFN, 'Sheet', feature);
                    groupTablePost = readtable(postGroupTableFN, 'Sheet', feature);
                    
                    %% read normalized EI
                    preGroupTableFN = strcat(groupAnalysisTablesFilePath, '','GroupAnalysis_ChannelCharacterization_Pre', selDetector, '.xls');
                    postGroupTableFN = strcat(groupAnalysisTablesFilePath, '','GroupAnalysis_ChannelCharacterization_Post', selDetector, '.xls');
                    groupTablePreEI = readtable(preGroupTableFN, 'Sheet', feature);
                    groupTablePostEI = readtable(postGroupTableFN, 'Sheet', feature);
                    groupTablePre.eiVals = groupTablePreEI.eiVals;
                    groupTablePost.eiVals = groupTablePostEI.eiVals;

                    %% Select channels to compare
                    tableDelIdxs = not(perctThValsPerPat(groupTablePre.patName, groupTablePre.eiVals, eiTh));
                    %tableDelIdxs = groupTablePre.eiVals < (prctile(groupTablePre.eiVals, eiTh));
                    groupTablePre(tableDelIdxs, :) = [];
                    groupTablePost(tableDelIdxs, :) = [];
                    
                    tableDelIdxs = groupTablePre.rftcVals > 0;
                    groupTablePre(tableDelIdxs, :) = [];
                    groupTablePost(tableDelIdxs, :) = [];
                    
                    pairedChanns = sum(not(strcmp(groupTablePre.channelLabels , groupTablePost.channelLabels))) == 0;
                    if not(pairedChanns)
                        stopHere = 1;
                    end
                    
                    %% Select outcomes
                    if strcmp(patsSel, 'allPatients')
                        analysisSubType = strcat(analysisSubType, '_', patsSel);
                    elseif strcmp(patsSel, 'improvedPatients')
                        analysisSubType = strcat(analysisSubType, '_', patsSel);
                        tableDelIdxs = groupTablePre.outcome <= outcomeTh;
                        groupTablePre(tableDelIdxs, :) = [];
                        tableDelIdxs = groupTablePost.outcome <= outcomeTh;
                        groupTablePost(tableDelIdxs, :) = [];
                    elseif strcmp(patsSel, 'nonImprovedPatients')
                        analysisSubType = strcat(analysisSubType, '_', patsSel);
                        tableDelIdxs = groupTablePre.outcome > outcomeTh;
                        groupTablePre(tableDelIdxs, :) = [];
                        tableDelIdxs = groupTablePost.outcome > outcomeTh;
                        groupTablePost(tableDelIdxs, :) = [];
                    end
                    
                    analysisSubType = strcat(analysisSubType, '_', normalization);
                    %%

                    analysis = cat(1, analysis, analysisSubType);
                    threshold = cat(1, threshold, th);
                    nrPats = cat(1, nrPats, length(unique(groupTablePost.patName)));
                    nrChanns = cat(1, nrChanns, length(groupTablePost.patName));
                    
                                        patNameCol = groupTablePre.patName;
                    patNames = unique(patNameCol);
                    iaValsSum = 0;
                    allHFOValsSum = 0; iesHFOValsSum = 0; isolHFOValsSum = 0;
                    allRippleValsSum = 0; iesRippleValsSum = 0; isolRippleValsSum = 0;
                    allFR_ValsSum = 0; iesFR_ValsSum = 0; isolFR_ValsSum = 0;
                    nrPatsForAvg = length(patNames);
                    pValLim = 0.05;
                    for pi = 1:nrPatsForAvg
                        patName = patNames{pi};
                        patSelIdx = ismember(patNameCol, patName);
                        perPatTablePre = groupTablePre(patSelIdx, :);
                        perPatTablePost = groupTablePost(patSelIdx, :);
                            
                        iaValsSum = iaValsSum + (signrank(perPatTablePre.iaVals, perPatTablePost.iaVals, 'tail', 'right') < pValLim);
                        allHFOValsSum = allHFOValsSum + (signrank(perPatTablePre.allHFOVals, perPatTablePost.allHFOVals, 'tail', 'right') < pValLim);
                        iesHFOValsSum = iesHFOValsSum + (signrank(perPatTablePre.iesHFOVals, perPatTablePost.iesHFOVals, 'tail', 'right') < pValLim);
                        isolHFOValsSum = isolHFOValsSum + (signrank(perPatTablePre.isolHFOVals, perPatTablePost.isolHFOVals, 'tail', 'right') < pValLim);

                        allRippleValsSum = allRippleValsSum + (signrank(perPatTablePre.allRippleVals, perPatTablePost.allRippleVals, 'tail', 'right') < pValLim);
                        iesRippleValsSum = iesRippleValsSum + (signrank(perPatTablePre.iesRippleVals, perPatTablePost.iesRippleVals, 'tail', 'right') < pValLim);
                        isolRippleValsSum = isolRippleValsSum + (signrank(perPatTablePre.isolRippleVals, perPatTablePost.isolRippleVals, 'tail', 'right') < pValLim);

                        allFR_ValsSum = allFR_ValsSum + (signrank(perPatTablePre.allFR_Vals, perPatTablePost.allFR_Vals, 'tail', 'right') < pValLim);
                        iesFR_ValsSum = iesFR_ValsSum + (signrank(perPatTablePre.iesFR_Vals, perPatTablePost.iesFR_Vals, 'tail', 'right') < pValLim);
                        isolFR_ValsSum = isolFR_ValsSum + (signrank(perPatTablePre.isolFR_Vals, perPatTablePost.isolFR_Vals, 'tail', 'right') < pValLim);
                    end

                    iaVals = cat(1, iaVals, iaValsSum / nrPatsForAvg*100);
                    allHFOVals = cat(1, allHFOVals, allHFOValsSum / nrPatsForAvg*100);
                    iesHFOVals = cat(1, iesHFOVals, iesHFOValsSum / nrPatsForAvg*100);
                    isolHFOVals = cat(1, isolHFOVals, isolHFOValsSum / nrPatsForAvg*100);
                    allRippleVals = cat(1, allRippleVals, allRippleValsSum / nrPatsForAvg*100);
                    iesRippleVals = cat(1, iesRippleVals, iesRippleValsSum / nrPatsForAvg*100);
                    isolRippleVals = cat(1, isolRippleVals, isolRippleValsSum / nrPatsForAvg*100);
                    allFR_Vals = cat(1, allFR_Vals, allFR_ValsSum / nrPatsForAvg*100);
                    iesFR_Vals = cat(1, iesFR_Vals, iesFR_ValsSum / nrPatsForAvg*100);
                    isolFR_Vals = cat(1, isolFR_Vals, isolFR_ValsSum / nrPatsForAvg*100);
                end
            end
        end
        iaVals = round(iaVals);
        allHFOVals = round(allHFOVals); iesHFOVals = round(iesHFOVals); isolHFOVals = round(isolHFOVals);
        allRippleVals = round(allRippleVals); iesRippleVals = round(iesRippleVals); isolRippleVals = round(isolRippleVals);
        allFR_Vals = round(allFR_Vals); iesFR_Vals = round(iesFR_Vals); isolFR_Vals = round(isolFR_Vals);
        
        zonesAnalysisTable = table(analysis, threshold, nrPats, nrChanns, iaVals, allHFOVals, iesHFOVals, isolHFOVals, allRippleVals, iesRippleVals, isolRippleVals, allFR_Vals, iesFR_Vals, isolFR_Vals); 
        spaceT = zonesAnalysisTable(1,:);
        for ti = 1:length(spaceT.Properties.VariableNames)
            spaceT.(spaceT.Properties.VariableNames{ti}) = nan(1,1);
        end
        zonesAnalysisTableSpaces = [];
        for si = 0:(size(zonesAnalysisTable,1)/3)-1
            idx = 3*si;
            zonesAnalysisTableSpaces = cat(1, zonesAnalysisTableSpaces, zonesAnalysisTable(idx+1:idx+3,:), spaceT);
        end

        newSelDetector = selDetector;
        if th > 0
            analysisPerZonesTablePath(end) = [];
            analysisPerZonesTablePath = strcat(analysisPerZonesTablePath, '_Depurated\', 'Th', num2str(th), '\');mkdir(analysisPerZonesTablePath);
            selDetector = strcat(selDetector, '_Th', num2str(th));
        end
        analysisPerZonesTableFN = strcat(analysisPerZonesTablePath, 'AnalysisPerZones_', newSelDetector,'.xls');
        writetable(zonesAnalysisTableSpaces, analysisPerZonesTableFN, 'Sheet', analysisType{1});
    end
end

function get_RFTC_Connected_Channels_Analysis_PrePost(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh)

    for fbi = 1:length(freqBandConnList)    
        analysisType = {'RFTC_Connected_Channels'};

        analysis = {};
        threshold = [];
        nrPats = [];
        nrChanns = [];
        iaVals = [];
        allHFOVals = []; iesHFOVals = []; isolHFOVals = [];
        allRippleVals = []; iesRippleVals = []; isolRippleVals = [];
        allFR_Vals = []; iesFR_Vals = []; isolFR_Vals = [];
        for th = thList
            for ni = 1:length(normOptions)
                for fi = 1:length(features)
                    for psi = 1:length(patsSelList)
                        freqBandConn = freqBandConnList{fbi};
                        patsSel = patsSelList{psi};
                        feature = features{fi};
                        normalization = normOptions{ni};
                        analysisSubType = feature;

                        selDetectorFolder = selDetector;
                        if th > 0
                            selDetectorFolder = strcat('MOSSDET_Depurated\Th', num2str(th));
                        end
                        groupAnalysisTablesFilePath = strcat(paths.workspacePath, 'GroupAnalysis_ChannelCharacterizationTables\', selDetectorFolder, '\');
                        preGroupTableFN = strcat(groupAnalysisTablesFilePath, normalization, 'GroupAnalysis_ChannelCharacterization_Pre', selDetector,  '_', freqBandConn, '.xls');
                        postGroupTableFN = strcat(groupAnalysisTablesFilePath, normalization, 'GroupAnalysis_ChannelCharacterization_Post', selDetector, '_', freqBandConn,'.xls');
                        groupTablePre = readtable(preGroupTableFN, 'Sheet', feature);
                        groupTablePost = readtable(postGroupTableFN, 'Sheet', feature);

                        %% read rftcElectroPhysioConnect
                        preGroupTableNormFN = strcat(groupAnalysisTablesFilePath, 'Normalized', 'GroupAnalysis_ChannelCharacterization_Pre', selDetector,  '_', freqBandConn, '.xls');
                        postGroupTableNormFN = strcat(groupAnalysisTablesFilePath, 'Normalized', 'GroupAnalysis_ChannelCharacterization_Post', selDetector, '_', freqBandConn,'.xls');
                        groupTablePre_RFTCconn = readtable(preGroupTableNormFN, 'Sheet', feature);
                        groupTablePost_RFTCconn = readtable(postGroupTableNormFN, 'Sheet', feature);
                        groupTablePre.rftcElectroPhysioConnect = groupTablePre_RFTCconn.rftcElectroPhysioConnect;
                        groupTablePost.rftcElectroPhysioConnect = groupTablePost_RFTCconn.rftcElectroPhysioConnect;


                        %% Select channels to compare
                        tableDelIdxs = not(perctThValsPerPat(groupTablePre.patName, groupTablePre.rftcElectroPhysioConnect, connTh));
                        %tableDelIdxs = groupTablePre.rftcElectroPhysioConnect < (prctile(groupTablePre.rftcElectroPhysioConnect, connTh));
                        groupTablePre(tableDelIdxs, :) = [];
                        groupTablePost(tableDelIdxs, :) = [];
                        
                        tableDelIdxs = groupTablePre.rftcVals > 0;
                        groupTablePre(tableDelIdxs, :) = [];
                        groupTablePost(tableDelIdxs, :) = [];
                        
                        pairedChanns = sum(not(strcmp(groupTablePre.channelLabels , groupTablePost.channelLabels))) == 0;
                        if not(pairedChanns)
                            stopHere = 1;
                        end

                        %% Select outcomes
                        if strcmp(patsSel, 'allPatients')
                            analysisSubType = strcat(analysisSubType, '_', patsSel);
                        elseif strcmp(patsSel, 'improvedPatients')
                            analysisSubType = strcat(analysisSubType, '_', patsSel);
                            tableDelIdxs = groupTablePre.outcome <= outcomeTh;
                            groupTablePre(tableDelIdxs, :) = [];
                            tableDelIdxs = groupTablePost.outcome <= outcomeTh;
                            groupTablePost(tableDelIdxs, :) = [];
                        elseif strcmp(patsSel, 'nonImprovedPatients')
                            analysisSubType = strcat(analysisSubType, '_', patsSel);
                            tableDelIdxs = groupTablePre.outcome > outcomeTh;
                            groupTablePre(tableDelIdxs, :) = [];
                            tableDelIdxs = groupTablePost.outcome > outcomeTh;
                            groupTablePost(tableDelIdxs, :) = [];
                        end
                        
                        analysisSubType = strcat(analysisSubType, '_', normalization);
                        %%

                        analysis = cat(1, analysis, analysisSubType);
                        threshold = cat(1, threshold, th);
                        nrPats = cat(1, nrPats, length(unique(groupTablePost.patName)));
                        nrChanns = cat(1, nrChanns, length(groupTablePost.patName));

                        patNameCol = groupTablePre.patName;
                        patNames = unique(patNameCol);
                        iaValsSum = 0;
                        allHFOValsSum = 0; iesHFOValsSum = 0; isolHFOValsSum = 0;
                        allRippleValsSum = 0; iesRippleValsSum = 0; isolRippleValsSum = 0;
                        allFR_ValsSum = 0; iesFR_ValsSum = 0; isolFR_ValsSum = 0;
                        nrPatsForAvg = length(patNames);
                        pValLim = 0.05;
                        {freqBandConn patsSel feature normalization}
                        for pi = 1:nrPatsForAvg
                            patName = patNames{pi};
                            patSelIdx = ismember(patNameCol, patName);
                            perPatTablePre = groupTablePre(patSelIdx, :);
                            perPatTablePost = groupTablePost(patSelIdx, :);

                            iaValsSum = iaValsSum + (signrank(perPatTablePre.iaVals, perPatTablePost.iaVals, 'tail', 'right') < pValLim);
                            allHFOValsSum = allHFOValsSum + (signrank(perPatTablePre.allHFOVals, perPatTablePost.allHFOVals, 'tail', 'right') < pValLim);
                            iesHFOValsSum = iesHFOValsSum + (signrank(perPatTablePre.iesHFOVals, perPatTablePost.iesHFOVals, 'tail', 'right') < pValLim);
                            isolHFOValsSum = isolHFOValsSum + (signrank(perPatTablePre.isolHFOVals, perPatTablePost.isolHFOVals, 'tail', 'right') < pValLim);

                            allRippleValsSum = allRippleValsSum + (signrank(perPatTablePre.allRippleVals, perPatTablePost.allRippleVals, 'tail', 'right') < pValLim);
                            iesRippleValsSum = iesRippleValsSum + (signrank(perPatTablePre.iesRippleVals, perPatTablePost.iesRippleVals, 'tail', 'right') < pValLim);
                            isolRippleValsSum = isolRippleValsSum + (signrank(perPatTablePre.isolRippleVals, perPatTablePost.isolRippleVals, 'tail', 'right') < pValLim);

                            allFR_ValsSum = allFR_ValsSum + (signrank(perPatTablePre.allFR_Vals, perPatTablePost.allFR_Vals, 'tail', 'right') < pValLim);
                            iesFR_ValsSum = iesFR_ValsSum + (signrank(perPatTablePre.iesFR_Vals, perPatTablePost.iesFR_Vals, 'tail', 'right') < pValLim);
                            isolFR_ValsSum = isolFR_ValsSum + (signrank(perPatTablePre.isolFR_Vals, perPatTablePost.isolFR_Vals, 'tail', 'right') < pValLim);
                        end

                        iaVals = cat(1, iaVals, iaValsSum / nrPatsForAvg*100);
                        allHFOVals = cat(1, allHFOVals, allHFOValsSum / nrPatsForAvg*100);
                        iesHFOVals = cat(1, iesHFOVals, iesHFOValsSum / nrPatsForAvg*100);
                        isolHFOVals = cat(1, isolHFOVals, isolHFOValsSum / nrPatsForAvg*100);
                        allRippleVals = cat(1, allRippleVals, allRippleValsSum / nrPatsForAvg*100);
                        iesRippleVals = cat(1, iesRippleVals, iesRippleValsSum / nrPatsForAvg*100);
                        isolRippleVals = cat(1, isolRippleVals, isolRippleValsSum / nrPatsForAvg*100);
                        allFR_Vals = cat(1, allFR_Vals, allFR_ValsSum / nrPatsForAvg*100);
                        iesFR_Vals = cat(1, iesFR_Vals, iesFR_ValsSum / nrPatsForAvg*100);
                        isolFR_Vals = cat(1, isolFR_Vals, isolFR_ValsSum / nrPatsForAvg*100);
                    end                
                end
            end
            iaVals = round(iaVals);
            allHFOVals = round(allHFOVals); iesHFOVals = round(iesHFOVals); isolHFOVals = round(isolHFOVals);
            allRippleVals = round(allRippleVals); iesRippleVals = round(iesRippleVals); isolRippleVals = round(isolRippleVals);
            allFR_Vals = round(allFR_Vals); iesFR_Vals = round(iesFR_Vals); isolFR_Vals = round(isolFR_Vals);

            zonesAnalysisTable = table(analysis, threshold, nrPats, nrChanns, iaVals, allHFOVals, iesHFOVals, isolHFOVals, allRippleVals, iesRippleVals, isolRippleVals, allFR_Vals, iesFR_Vals, isolFR_Vals); 
            spaceT = zonesAnalysisTable(1,:);
            for ti = 1:length(spaceT.Properties.VariableNames)
                spaceT.(spaceT.Properties.VariableNames{ti}) = nan(1,1);
            end
            zonesAnalysisTableSpaces = [];
            for si = 0:(size(zonesAnalysisTable,1)/3)-1
                idx = 3*si;
                zonesAnalysisTableSpaces = cat(1, zonesAnalysisTableSpaces, zonesAnalysisTable(idx+1:idx+3,:), spaceT);
            end

            newSelDetector = selDetector;
            newAnalysisPerZonesTablePath = analysisPerZonesTablePath;
            if th > 0
                newAnalysisPerZonesTablePath(end) = [];
                newAnalysisPerZonesTablePath = strcat(newAnalysisPerZonesTablePath, '_Depurated\', 'Th', num2str(th), '\');mkdir(analysisPerZonesTablePath);
                newSelDetector = strcat(selDetector, '_Th', num2str(th));
            end
            analysisPerZonesTableFN = strcat(newAnalysisPerZonesTablePath, 'RFTC_ConnectedChannels_Analysis_', newSelDetector, '.xls');
            sheetName = strcat('RFTC-Conn-', freqBandConn);
            writetable(zonesAnalysisTableSpaces, analysisPerZonesTableFN, 'Sheet', sheetName);
        end
    end
end

function get_RFTC_NonConnected_Channels_Analysis_PrePost(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh)

    for fbi = 1:length(freqBandConnList)    
        analysisType = {'RFTC_Connected_Channels'};

        analysis = {};
        threshold = [];
        nrPats = [];
        nrChanns = [];
        iaVals = [];
        allHFOVals = []; iesHFOVals = []; isolHFOVals = [];
        allRippleVals = []; iesRippleVals = []; isolRippleVals = [];
        allFR_Vals = []; iesFR_Vals = []; isolFR_Vals = [];
        for th = thList
            for ni = 1:length(normOptions)
                for fi = 1:length(features)
                    for psi = 1:length(patsSelList)
                        freqBandConn = freqBandConnList{fbi};
                        patsSel = patsSelList{psi};
                        feature = features{fi};
                        normalization = normOptions{ni};
                        analysisSubType = feature;

                        selDetectorFolder = selDetector;
                        if th > 0
                            selDetectorFolder = strcat('MOSSDET_Depurated\Th', num2str(th));
                        end
                        groupAnalysisTablesFilePath = strcat(paths.workspacePath, 'GroupAnalysis_ChannelCharacterizationTables\', selDetectorFolder, '\');
                        preGroupTableFN = strcat(groupAnalysisTablesFilePath, normalization, 'GroupAnalysis_ChannelCharacterization_Pre', selDetector,  '_', freqBandConn, '.xls');
                        postGroupTableFN = strcat(groupAnalysisTablesFilePath, normalization, 'GroupAnalysis_ChannelCharacterization_Post', selDetector, '_', freqBandConn,'.xls');
                        groupTablePre = readtable(preGroupTableFN, 'Sheet', feature);
                        groupTablePost = readtable(postGroupTableFN, 'Sheet', feature);

                        %% read rftcElectroPhysioConnect
                        preGroupTableNormFN = strcat(groupAnalysisTablesFilePath, 'Normalized', 'GroupAnalysis_ChannelCharacterization_Pre', selDetector,  '_', freqBandConn, '.xls');
                        postGroupTableNormFN = strcat(groupAnalysisTablesFilePath, 'Normalized', 'GroupAnalysis_ChannelCharacterization_Post', selDetector, '_', freqBandConn,'.xls');
                        groupTablePre_RFTCconn = readtable(preGroupTableNormFN, 'Sheet', feature);
                        groupTablePost_RFTCconn = readtable(postGroupTableNormFN, 'Sheet', feature);
                        groupTablePre.rftcElectroPhysioConnect = groupTablePre_RFTCconn.rftcElectroPhysioConnect;
                        groupTablePost.rftcElectroPhysioConnect = groupTablePost_RFTCconn.rftcElectroPhysioConnect;


                        %% Select channels to compare                        
                        tableDelIdxs = perctThValsPerPat(groupTablePre.patName, groupTablePre.rftcElectroPhysioConnect, connTh);
                        %tableDelIdxs = groupTablePre.rftcElectroPhysioConnect < (prctile(groupTablePre.rftcElectroPhysioConnect, connTh));
                        groupTablePre(tableDelIdxs, :) = [];
                        groupTablePost(tableDelIdxs, :) = [];
                        
                        tableDelIdxs = groupTablePre.rftcVals > 0;
                        groupTablePre(tableDelIdxs, :) = [];
                        groupTablePost(tableDelIdxs, :) = [];
                        
                        pairedChanns = sum(not(strcmp(groupTablePre.channelLabels , groupTablePost.channelLabels))) == 0;
                        if not(pairedChanns)
                            stopHere = 1;
                        end
                        

                        %% Select outcomes
                        if strcmp(patsSel, 'allPatients')
                            analysisSubType = strcat(analysisSubType, '_', patsSel);
                        elseif strcmp(patsSel, 'improvedPatients')
                            analysisSubType = strcat(analysisSubType, '_', patsSel);
                            tableDelIdxs = groupTablePre.outcome <= outcomeTh;
                            groupTablePre(tableDelIdxs, :) = [];
                            tableDelIdxs = groupTablePost.outcome <= outcomeTh;
                            groupTablePost(tableDelIdxs, :) = [];
                        elseif strcmp(patsSel, 'nonImprovedPatients')
                            analysisSubType = strcat(analysisSubType, '_', patsSel);
                            tableDelIdxs = groupTablePre.outcome > outcomeTh;
                            groupTablePre(tableDelIdxs, :) = [];
                            tableDelIdxs = groupTablePost.outcome > outcomeTh;
                            groupTablePost(tableDelIdxs, :) = [];
                        end
                        
                        analysisSubType = strcat(analysisSubType, '_', normalization);
                        %%

                        analysis = cat(1, analysis, analysisSubType);
                        threshold = cat(1, threshold, th);
                        nrPats = cat(1, nrPats, length(unique(groupTablePost.patName)));
                        nrChanns = cat(1, nrChanns, length(groupTablePost.patName));

                        patNameCol = groupTablePre.patName;
                        patNames = unique(patNameCol);
                        iaValsSum = 0;
                        allHFOValsSum = 0; iesHFOValsSum = 0; isolHFOValsSum = 0;
                        allRippleValsSum = 0; iesRippleValsSum = 0; isolRippleValsSum = 0;
                        allFR_ValsSum = 0; iesFR_ValsSum = 0; isolFR_ValsSum = 0;
                        nrPatsForAvg = length(patNames);
                        pValLim = 0.05;
                        for pi = 1:nrPatsForAvg
                            patName = patNames{pi};
                            patSelIdx = ismember(patNameCol, patName);
                            perPatTablePre = groupTablePre(patSelIdx, :);
                            perPatTablePost = groupTablePost(patSelIdx, :);

                            iaValsSum = iaValsSum + (signrank(perPatTablePre.iaVals, perPatTablePost.iaVals, 'tail', 'right') < pValLim);
                            allHFOValsSum = allHFOValsSum + (signrank(perPatTablePre.allHFOVals, perPatTablePost.allHFOVals, 'tail', 'right') < pValLim);
                            iesHFOValsSum = iesHFOValsSum + (signrank(perPatTablePre.iesHFOVals, perPatTablePost.iesHFOVals, 'tail', 'right') < pValLim);
                            isolHFOValsSum = isolHFOValsSum + (signrank(perPatTablePre.isolHFOVals, perPatTablePost.isolHFOVals, 'tail', 'right') < pValLim);

                            allRippleValsSum = allRippleValsSum + (signrank(perPatTablePre.allRippleVals, perPatTablePost.allRippleVals, 'tail', 'right') < pValLim);
                            iesRippleValsSum = iesRippleValsSum + (signrank(perPatTablePre.iesRippleVals, perPatTablePost.iesRippleVals, 'tail', 'right') < pValLim);
                            isolRippleValsSum = isolRippleValsSum + (signrank(perPatTablePre.isolRippleVals, perPatTablePost.isolRippleVals, 'tail', 'right') < pValLim);

                            allFR_ValsSum = allFR_ValsSum + (signrank(perPatTablePre.allFR_Vals, perPatTablePost.allFR_Vals, 'tail', 'right') < pValLim);
                            iesFR_ValsSum = iesFR_ValsSum + (signrank(perPatTablePre.iesFR_Vals, perPatTablePost.iesFR_Vals, 'tail', 'right') < pValLim);
                            isolFR_ValsSum = isolFR_ValsSum + (signrank(perPatTablePre.isolFR_Vals, perPatTablePost.isolFR_Vals, 'tail', 'right') < pValLim);
                        end

                        iaVals = cat(1, iaVals, iaValsSum / nrPatsForAvg*100);
                        allHFOVals = cat(1, allHFOVals, allHFOValsSum / nrPatsForAvg*100);
                        iesHFOVals = cat(1, iesHFOVals, iesHFOValsSum / nrPatsForAvg*100);
                        isolHFOVals = cat(1, isolHFOVals, isolHFOValsSum / nrPatsForAvg*100);
                        allRippleVals = cat(1, allRippleVals, allRippleValsSum / nrPatsForAvg*100);
                        iesRippleVals = cat(1, iesRippleVals, iesRippleValsSum / nrPatsForAvg*100);
                        isolRippleVals = cat(1, isolRippleVals, isolRippleValsSum / nrPatsForAvg*100);
                        allFR_Vals = cat(1, allFR_Vals, allFR_ValsSum / nrPatsForAvg*100);
                        iesFR_Vals = cat(1, iesFR_Vals, iesFR_ValsSum / nrPatsForAvg*100);
                        isolFR_Vals = cat(1, isolFR_Vals, isolFR_ValsSum / nrPatsForAvg*100);
                    
                    end                
                end
            end
            iaVals = round(iaVals);
            allHFOVals = round(allHFOVals); iesHFOVals = round(iesHFOVals); isolHFOVals = round(isolHFOVals);
            allRippleVals = round(allRippleVals); iesRippleVals = round(iesRippleVals); isolRippleVals = round(isolRippleVals);
            allFR_Vals = round(allFR_Vals); iesFR_Vals = round(iesFR_Vals); isolFR_Vals = round(isolFR_Vals);
        
            zonesAnalysisTable = table(analysis, threshold, nrPats, nrChanns, iaVals, allHFOVals, iesHFOVals, isolHFOVals, allRippleVals, iesRippleVals, isolRippleVals, allFR_Vals, iesFR_Vals, isolFR_Vals); 
            spaceT = zonesAnalysisTable(1,:);
            for ti = 1:length(spaceT.Properties.VariableNames)
                spaceT.(spaceT.Properties.VariableNames{ti}) = nan(1,1);
            end
            zonesAnalysisTableSpaces = [];
            for si = 0:(size(zonesAnalysisTable,1)/3)-1
                idx = 3*si;
                zonesAnalysisTableSpaces = cat(1, zonesAnalysisTableSpaces, zonesAnalysisTable(idx+1:idx+3,:), spaceT);
            end
            
            newSelDetector = selDetector;
            newAnalysisPerZonesTablePath = analysisPerZonesTablePath;
            if th > 0
                newAnalysisPerZonesTablePath(end) = [];
                newAnalysisPerZonesTablePath = strcat(newAnalysisPerZonesTablePath, '_Depurated\', 'Th', num2str(th), '\');mkdir(analysisPerZonesTablePath);
                newSelDetector = strcat(selDetector, '_Th', num2str(th));
            end
            analysisPerZonesTableFN = strcat(newAnalysisPerZonesTablePath, 'RFTC_NonConnectedChannels_Analysis_', newSelDetector, '.xls');
            sheetName = strcat('RFTC-Conn-', freqBandConn);
            writetable(zonesAnalysisTableSpaces, analysisPerZonesTableFN, 'Sheet', sheetName);
        end
    end
end

function allPatsKeepIdxs = perctThValsPerPat(patNameCol, vals, th)
    patNames = unique(patNameCol);
    allPatsKeepIdxs = zeros(length(vals),1);
    for pi = 1:length(patNames)
        patName = patNames{pi};
        patSelIdx = ismember(patNameCol, patName);
        patVals = vals(patSelIdx);
        patKeepIdx = patVals >= prctile(patVals, th);
        allPatsKeepIdxs(patSelIdx) = patKeepIdx;        
    end
    allPatsKeepIdxs = logical(allPatsKeepIdxs);
end