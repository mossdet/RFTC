clc; clear all; close all;
paths = getFilesPaths();
preFiles  = getPreFiles();%getAllFiles, getPreFiles, getPostFiles
postFiles  = getPostFiles();%getAllFiles, getPreFiles, getPostFiles
analysisPerZonesPath = 'F:\ForschungsProjekte\RFTC\MATLAB_ver2\DetectHFO\BiomarkersAndZones_Correlations\';

selDetector = 'MOSSDET'; 'Delphos'; 'MOSSDET'; 'MOSSDET_Depurated';
analysisPerZonesTablePath = strcat(analysisPerZonesPath, selDetector, '\');
features = {'OccRate', 'Power', 'Frequency'}; {'OccRate', 'Power', 'Frequency'};
normOptions = {'Normalized'}; {'Normalized', '', 'Scaled'}; %'Normalized'; % {'', 'Scaled', 'Normalized'};
freqBandConnList = { 'maxAllBands'};
thList = 0;
patsSelList = {'allPatients', 'improvedPatients', 'nonImprovedPatients'};
biomarkersList =  {'iaVals', 'allHFOVals', 'iesHFOVals', 'isolHFOVals', 'allRippleVals', 'iesRippleVals', 'isolRippleVals', 'allFR_Vals', 'iesFR_Vals', 'isolFR_Vals'};
biomarkersList =  {'iaVals', 'allHFOVals', 'iesHFOVals'};

params.analysisPerZonesTablePath = analysisPerZonesTablePath;
params.selDetector = selDetector;
params.connTh = 75;
params.eiTh = 75;
params.outcomeTh = 49;

%% Analysis Loop
for nmi = 1:length(normOptions)
    for psi = 1:length(patsSelList)
        for fi = 1:length(features)
            %% Set params
            params.th = 0;
            params.normalization = normOptions{nmi};
            params.feature = features{fi};
            params.patsSel = patsSelList{psi};
            params.analysisSubType = features{fi};
            params.freqBandConn = 'maxAllBands';
            params.biomarkersList = biomarkersList;

            [groupTablePre, groupTablePost] = readTables(paths, params);
            plotBiomarker_ZoneCorrelation(params, groupTablePre);
        end
    end
end

function biomarkerActvDiff = fillResultsStruct(params, groupTablePreIn, groupTablePostIn, groupTablePreOut, groupTablePostOut)

    biomarkerActvDiff.in.pre = groupTablePreIn.(params.biomarker);
    biomarkerActvDiff.in.post = groupTablePostIn.(params.biomarker);
    
    biomarkerActvDiff.out.pre = groupTablePreOut.(params.biomarker);
    biomarkerActvDiff.out.post = groupTablePostOut.(params.biomarker);
    
    biomarkerActvDiff.in.diff = groupTablePreIn.(params.biomarker) - groupTablePostIn.(params.biomarker);
    biomarkerActvDiff.out.diff = groupTablePreOut.(params.biomarker) - groupTablePostOut.(params.biomarker);
    
    biomarkerActvDiff.in.prePostP = signrank(groupTablePreIn.(params.biomarker), groupTablePostIn.(params.biomarker), 'tail', 'right');
    biomarkerActvDiff.out.prePostP = signrank(groupTablePreOut.(params.biomarker), groupTablePostOut.(params.biomarker), 'tail', 'right');
    
    biomarkerActvDiff.in.nrChanns = length(groupTablePreIn.patName);
    biomarkerActvDiff.out.nrChanns = length(groupTablePreOut.patName);
    
    biomarkerActvDiff.in.nrPats = length(unique(groupTablePreIn.patName));
    biomarkerActvDiff.out.nrPats = length(unique(groupTablePreOut.patName));
    
    [p, h] = ranksum(biomarkerActvDiff.in.diff, biomarkerActvDiff.out.diff);
    biomarkerActvDiff.inOutP = p;
end

function groupTable = removeRFTC_Channels(groupTable)
    tableDelChannsRFTC = groupTable.rftcVals > 0;
    groupTable(tableDelChannsRFTC, :) = [];
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

function [groupTablePre, groupTablePost] = readTables(paths, params)
    %% Read tables
    groupAnalysisTablesFilePath = strcat(paths.workspacePath, 'GroupAnalysis_ChannelCharacterizationTables\', params.selDetector, '\');
    preGroupTableFN = strcat(groupAnalysisTablesFilePath, params.normalization, 'GroupAnalysis_ChannelCharacterization_Pre', params.selDetector, '.xls');
    postGroupTableFN = strcat(groupAnalysisTablesFilePath, params.normalization, 'GroupAnalysis_ChannelCharacterization_Post', params.selDetector, '.xls');
    groupTablePre = readtable(preGroupTableFN, 'Sheet', params.feature);
    groupTablePost = readtable(postGroupTableFN, 'Sheet', params.feature);
    
    %% Truncate outliers
    %groupTablePre = truncateOutliers(groupTablePre);
    %groupTablePost = truncateOutliers(groupTablePost);
    
    %% read normalized EI
    preGroupTableFN = strcat(groupAnalysisTablesFilePath, 'Normalized','GroupAnalysis_ChannelCharacterization_Pre', params.selDetector, '.xls');%'Normalized' Scaled
    postGroupTableFN = strcat(groupAnalysisTablesFilePath, 'Normalized','GroupAnalysis_ChannelCharacterization_Post', params.selDetector, '.xls');
    groupTablePreEI_correct = readtable(preGroupTableFN, 'Sheet', 'OccRate');
    groupTablePostEI_correct = readtable(postGroupTableFN, 'Sheet', 'OccRate');
    groupTablePre.eiVals = groupTablePreEI_correct.eiVals;
    groupTablePost.eiVals = groupTablePostEI_correct.eiVals;
    
    %% read normalized rftcElectroPhysioConnect
    preGroupTableNormFN = strcat(groupAnalysisTablesFilePath, 'Normalized', 'GroupAnalysis_ChannelCharacterization_Pre', params.selDetector,  '_', params.freqBandConn, '.xls');
    postGroupTableNormFN = strcat(groupAnalysisTablesFilePath, 'Normalized', 'GroupAnalysis_ChannelCharacterization_Post', params.selDetector, '_', params.freqBandConn,'.xls');
    groupTablePre_RFTCconn = readtable(preGroupTableNormFN, 'Sheet', 'OccRate');
    groupTablePost_RFTCconn = readtable(postGroupTableNormFN, 'Sheet', 'OccRate');
    groupTablePre.rftcElectroPhysioConnect = groupTablePre_RFTCconn.rftcElectroPhysioConnect;
    groupTablePost.rftcElectroPhysioConnect = groupTablePost_RFTCconn.rftcElectroPhysioConnect;
    
    %% Select outcomes
    if strcmp(params.patsSel, 'improvedPatients')        
        tableDelIdxsIn = groupTablePre.outcome <= params.outcomeTh;
        groupTablePre(tableDelIdxsIn, :) = [];
        groupTablePost(tableDelIdxsIn, :) = [];        
    elseif strcmp(params.patsSel, 'nonImprovedPatients')
        tableDelIdxsIn = groupTablePre.outcome > params.outcomeTh;
        groupTablePre(tableDelIdxsIn, :) = [];
        groupTablePost(tableDelIdxsIn, :) = []; 
    end
end

function groupTable = truncateOutliers(groupTable)

    allPatNames = unique(groupTable.patName);
    nrPatients = length(allPatNames);
    columnNames = groupTable.Properties.VariableNames;

    for pi = 1:nrPatients
        patSelIdx = ismember(groupTable.patName, allPatNames{pi});
        sum(patSelIdx)
        for ci = 1:length(columnNames)
            column = columnNames{ci};
            columnPatVals = groupTable.(column)(patSelIdx);
            if ci >= 9 %strcmp(class(columnPatVals), 'double') % ci >= 9
                meanVal = mean(columnPatVals);
                medianVal = median(columnPatVals);
                stdVal = std(columnPatVals);
                interQR = iqr(columnPatVals);
                quart1 = prctile(columnPatVals, 25);
                quart3 = prctile(columnPatVals, 75);

                lowBound = quart1 - 10*interQR;
                highBound = quart3 + 10*interQR;
                selLowOutliers = (columnPatVals < lowBound);
                selHighOutliers = (columnPatVals > highBound);
                if(sum(selLowOutliers) > 0)
                    [quart1 quart3 interQR meanVal medianVal]
                    columnPatVals(selHighOutliers)'
                    columnPatVals(selLowOutliers) = lowBound;
                end
                if(sum(selHighOutliers) > 0)
                    [quart1 quart3 interQR meanVal medianVal]
                    columnPatVals(selHighOutliers)'
                    columnPatVals(selHighOutliers) = highBound;
                end
                groupTable.(column)(patSelIdx) = columnPatVals;
            end
        end
    end
end