clc; clear all; close all;
paths = getFilesPaths();
preFiles  = getPreFiles();%getAllFiles, getPreFiles, getPostFiles
postFiles  = getPostFiles();%getAllFiles, getPreFiles, getPostFiles
analysisPerZonesPath = 'F:\ForschungsProjekte\RFTC\MATLAB\DetectHFO\AnalysisPerZones\';

selDetector = 'MOSSDET'; 'Delphos'; 'MOSSDET'; 'MOSSDET_Depurated';
features = {'OccRate', 'Power', 'Frequency'};
normOptions = {'', 'Scaled', 'Normalized'};
normOptions = {'Normalized'};
freqBandConnList = {'delta', 'theta', 'alpha', 'beta', 'gamma', 'highGamma', 'ripple', 'fr', 'maxAllBands', 'meanAllBands'};
%freqBandConnList = {'maxAllBands'};
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

get_Biomarker_EI_Correlation(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh);

%get_Same_BrainLobe_Channels_Analysis_PrePost(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh);
%get_Opposite_BrainLobe_Channels_Analysis_PrePost(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh);

get_Same_BrainParcel_Channels_Analysis_PrePost(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh);
get_Opposite_BrainParcel_Channels_Analysis_PrePost(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh);

get_Same_Hemisphere_Channels_Analysis_PrePost(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh);
get_Opposite_Hemisphere_Channels_Analysis_PrePost(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh);

%All Pats
get_Biomarker_EI_Correlation(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh);
getAllChannelsPairedAnalysis_PrePost(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh);
getRFTC_Channels_PairedAnalysis_PrePost(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh);
get_AllMinusRFTC_Channels_PairedAnalysis_PrePost(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh);
get_HighEI_Channels_Analysis_PrePost(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh);
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
                    
                    %testNormality(groupTablePre);
                    
                    iaVals = cat(1, iaVals, getSigCorr(groupTablePre.iaVals, groupTablePre.eiVals, 'Type','Spearman'));

                    allHFOVals = cat(1, allHFOVals, getSigCorr(groupTablePre.allHFOVals, groupTablePre.eiVals, 'Type','Spearman'));
                    iesHFOVals = cat(1, iesHFOVals, getSigCorr(groupTablePre.iesHFOVals, groupTablePre.eiVals, 'Type','Spearman'));
                    isolHFOVals = cat(1, isolHFOVals, getSigCorr(groupTablePre.isolHFOVals, groupTablePre.eiVals, 'Type','Spearman'));

                    allRippleVals = cat(1, allRippleVals, getSigCorr(groupTablePre.allRippleVals, groupTablePre.eiVals, 'Type','Spearman'));
                    iesRippleVals = cat(1, iesRippleVals, getSigCorr(groupTablePre.iesRippleVals, groupTablePre.eiVals, 'Type','Spearman'));
                    isolRippleVals = cat(1, isolRippleVals, getSigCorr(groupTablePre.isolRippleVals, groupTablePre.eiVals, 'Type','Spearman'));

                    allFR_Vals = cat(1, allFR_Vals, getSigCorr(groupTablePre.allFR_Vals, groupTablePre.eiVals, 'Type','Spearman'));
                    iesFR_Vals = cat(1, iesFR_Vals, getSigCorr(groupTablePre.iesFR_Vals, groupTablePre.eiVals, 'Type','Spearman'));
                    isolFR_Vals = cat(1, isolFR_Vals, getSigCorr(groupTablePre.isolFR_Vals, groupTablePre.eiVals, 'Type','Spearman'));
                end
            end
        end    
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

                    iaVals = cat(1, iaVals, signrank(groupTablePre.iaVals, groupTablePost.iaVals, 'tail', 'right'));

                    allHFOVals = cat(1, allHFOVals, signrank(groupTablePre.allHFOVals, groupTablePost.allHFOVals, 'tail', 'right'));
                    iesHFOVals = cat(1, iesHFOVals, signrank(groupTablePre.iesHFOVals, groupTablePost.iesHFOVals, 'tail', 'right'));
                    isolHFOVals = cat(1, isolHFOVals, signrank(groupTablePre.isolHFOVals, groupTablePost.isolHFOVals, 'tail', 'right'));

                    allRippleVals = cat(1, allRippleVals, signrank(groupTablePre.allRippleVals, groupTablePost.allRippleVals, 'tail', 'right'));
                    iesRippleVals = cat(1, iesRippleVals, signrank(groupTablePre.iesRippleVals, groupTablePost.iesRippleVals, 'tail', 'right'));
                    isolRippleVals = cat(1, isolRippleVals, signrank(groupTablePre.isolRippleVals, groupTablePost.isolRippleVals, 'tail', 'right'));

                    allFR_Vals = cat(1, allFR_Vals, signrank(groupTablePre.allFR_Vals, groupTablePost.allFR_Vals, 'tail', 'right'));
                    iesFR_Vals = cat(1, iesFR_Vals, signrank(groupTablePre.iesFR_Vals, groupTablePost.iesFR_Vals, 'tail', 'right'));
                    isolFR_Vals = cat(1, isolFR_Vals, signrank(groupTablePre.isolFR_Vals, groupTablePost.isolFR_Vals, 'tail', 'right'));
                end
            end
        end    
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

                    iaVals = cat(1, iaVals, signrank(groupTablePre.iaVals, groupTablePost.iaVals, 'tail', 'right'));

                    allHFOVals = cat(1, allHFOVals, signrank(groupTablePre.allHFOVals, groupTablePost.allHFOVals, 'tail', 'right'));
                    iesHFOVals = cat(1, iesHFOVals, signrank(groupTablePre.iesHFOVals, groupTablePost.iesHFOVals, 'tail', 'right'));
                    isolHFOVals = cat(1, isolHFOVals, signrank(groupTablePre.isolHFOVals, groupTablePost.isolHFOVals, 'tail', 'right'));

                    allRippleVals = cat(1, allRippleVals, signrank(groupTablePre.allRippleVals, groupTablePost.allRippleVals, 'tail', 'right'));
                    iesRippleVals = cat(1, iesRippleVals, signrank(groupTablePre.iesRippleVals, groupTablePost.iesRippleVals, 'tail', 'right'));
                    isolRippleVals = cat(1, isolRippleVals, signrank(groupTablePre.isolRippleVals, groupTablePost.isolRippleVals, 'tail', 'right'));

                    allFR_Vals = cat(1, allFR_Vals, signrank(groupTablePre.allFR_Vals, groupTablePost.allFR_Vals, 'tail', 'right'));
                    iesFR_Vals = cat(1, iesFR_Vals, signrank(groupTablePre.iesFR_Vals, groupTablePost.iesFR_Vals, 'tail', 'right'));
                    isolFR_Vals = cat(1, isolFR_Vals, signrank(groupTablePre.isolFR_Vals, groupTablePost.isolFR_Vals, 'tail', 'right'));
                end
            end
        end    
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
                    
                    iaVals = cat(1, iaVals, signrank(groupTablePre.iaVals, groupTablePost.iaVals, 'tail', 'right'));

                    allHFOVals = cat(1, allHFOVals, signrank(groupTablePre.allHFOVals, groupTablePost.allHFOVals, 'tail', 'right'));
                    iesHFOVals = cat(1, iesHFOVals, signrank(groupTablePre.iesHFOVals, groupTablePost.iesHFOVals, 'tail', 'right'));
                    isolHFOVals = cat(1, isolHFOVals, signrank(groupTablePre.isolHFOVals, groupTablePost.isolHFOVals, 'tail', 'right'));

                    allRippleVals = cat(1, allRippleVals, signrank(groupTablePre.allRippleVals, groupTablePost.allRippleVals, 'tail', 'right'));
                    iesRippleVals = cat(1, iesRippleVals, signrank(groupTablePre.iesRippleVals, groupTablePost.iesRippleVals, 'tail', 'right'));
                    isolRippleVals = cat(1, isolRippleVals, signrank(groupTablePre.isolRippleVals, groupTablePost.isolRippleVals, 'tail', 'right'));

                    allFR_Vals = cat(1, allFR_Vals, signrank(groupTablePre.allFR_Vals, groupTablePost.allFR_Vals, 'tail', 'right'));
                    iesFR_Vals = cat(1, iesFR_Vals, signrank(groupTablePre.iesFR_Vals, groupTablePost.iesFR_Vals, 'tail', 'right'));
                    isolFR_Vals = cat(1, isolFR_Vals, signrank(groupTablePre.isolFR_Vals, groupTablePost.isolFR_Vals, 'tail', 'right'));
                end
            end
        end    
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
                    
                    iaVals = cat(1, iaVals, signrank(groupTablePre.iaVals, groupTablePost.iaVals, 'tail', 'right'));

                    allHFOVals = cat(1, allHFOVals, signrank(groupTablePre.allHFOVals, groupTablePost.allHFOVals, 'tail', 'right'));
                    iesHFOVals = cat(1, iesHFOVals, signrank(groupTablePre.iesHFOVals, groupTablePost.iesHFOVals, 'tail', 'right'));
                    isolHFOVals = cat(1, isolHFOVals, signrank(groupTablePre.isolHFOVals, groupTablePost.isolHFOVals, 'tail', 'right'));

                    allRippleVals = cat(1, allRippleVals, signrank(groupTablePre.allRippleVals, groupTablePost.allRippleVals, 'tail', 'right'));
                    iesRippleVals = cat(1, iesRippleVals, signrank(groupTablePre.iesRippleVals, groupTablePost.iesRippleVals, 'tail', 'right'));
                    isolRippleVals = cat(1, isolRippleVals, signrank(groupTablePre.isolRippleVals, groupTablePost.isolRippleVals, 'tail', 'right'));

                    allFR_Vals = cat(1, allFR_Vals, signrank(groupTablePre.allFR_Vals, groupTablePost.allFR_Vals, 'tail', 'right'));
                    iesFR_Vals = cat(1, iesFR_Vals, signrank(groupTablePre.iesFR_Vals, groupTablePost.iesFR_Vals, 'tail', 'right'));
                    isolFR_Vals = cat(1, isolFR_Vals, signrank(groupTablePre.isolFR_Vals, groupTablePost.isolFR_Vals, 'tail', 'right'));
                end
            end
        end    
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

                        iaVals = cat(1, iaVals, signrank(groupTablePre.iaVals, groupTablePost.iaVals, 'tail', 'right'));

                        allHFOVals = cat(1, allHFOVals, signrank(groupTablePre.allHFOVals, groupTablePost.allHFOVals, 'tail', 'right'));
                        iesHFOVals = cat(1, iesHFOVals, signrank(groupTablePre.iesHFOVals, groupTablePost.iesHFOVals, 'tail', 'right'));
                        isolHFOVals = cat(1, isolHFOVals, signrank(groupTablePre.isolHFOVals, groupTablePost.isolHFOVals, 'tail', 'right'));

                        allRippleVals = cat(1, allRippleVals, signrank(groupTablePre.allRippleVals, groupTablePost.allRippleVals, 'tail', 'right'));
                        iesRippleVals = cat(1, iesRippleVals, signrank(groupTablePre.iesRippleVals, groupTablePost.iesRippleVals, 'tail', 'right'));
                        isolRippleVals = cat(1, isolRippleVals, signrank(groupTablePre.isolRippleVals, groupTablePost.isolRippleVals, 'tail', 'right'));

                        allFR_Vals = cat(1, allFR_Vals, signrank(groupTablePre.allFR_Vals, groupTablePost.allFR_Vals, 'tail', 'right'));
                        iesFR_Vals = cat(1, iesFR_Vals, signrank(groupTablePre.iesFR_Vals, groupTablePost.iesFR_Vals, 'tail', 'right'));
                        isolFR_Vals = cat(1, isolFR_Vals, signrank(groupTablePre.isolFR_Vals, groupTablePost.isolFR_Vals, 'tail', 'right'));
                    end                
                end
            end 
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

                        iaVals = cat(1, iaVals, signrank(groupTablePre.iaVals, groupTablePost.iaVals, 'tail', 'right'));

                        allHFOVals = cat(1, allHFOVals, signrank(groupTablePre.allHFOVals, groupTablePost.allHFOVals, 'tail', 'right'));
                        iesHFOVals = cat(1, iesHFOVals, signrank(groupTablePre.iesHFOVals, groupTablePost.iesHFOVals, 'tail', 'right'));
                        isolHFOVals = cat(1, isolHFOVals, signrank(groupTablePre.isolHFOVals, groupTablePost.isolHFOVals, 'tail', 'right'));

                        allRippleVals = cat(1, allRippleVals, signrank(groupTablePre.allRippleVals, groupTablePost.allRippleVals, 'tail', 'right'));
                        iesRippleVals = cat(1, iesRippleVals, signrank(groupTablePre.iesRippleVals, groupTablePost.iesRippleVals, 'tail', 'right'));
                        isolRippleVals = cat(1, isolRippleVals, signrank(groupTablePre.isolRippleVals, groupTablePost.isolRippleVals, 'tail', 'right'));

                        allFR_Vals = cat(1, allFR_Vals, signrank(groupTablePre.allFR_Vals, groupTablePost.allFR_Vals, 'tail', 'right'));
                        iesFR_Vals = cat(1, iesFR_Vals, signrank(groupTablePre.iesFR_Vals, groupTablePost.iesFR_Vals, 'tail', 'right'));
                        isolFR_Vals = cat(1, isolFR_Vals, signrank(groupTablePre.isolFR_Vals, groupTablePost.isolFR_Vals, 'tail', 'right'));
                    end                
                end
            end    
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

