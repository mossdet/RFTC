clc; clear all; close all;
paths = getFilesPaths();
preFiles  = getPreFiles();%getAllFiles, getPreFiles, getPostFiles
postFiles  = getPostFiles();%getAllFiles, getPreFiles, getPostFiles
analysisPerZonesPath = 'F:\ForschungsProjekte\RFTC\MATLAB\DetectHFO\AnalysisPerZones_Boxplots\';

selDetector = 'MOSSDET';
analysisPerZonesTablePath = strcat(analysisPerZonesPath, selDetector, '\');
normOptions = {'', 'Normalized', '', 'Scaled'}; %'Normalized'; % {'', 'Scaled', 'Normalized'};
freqBandConnList = { 'maxAllBands'};
features = {'OccRate', 'Power', 'Frequency'};
thList = 0;
patsSelList = {'allPatients'};
biomarkersList =  {'iesHFOVals', 'iaVals', 'allHFOVals', 'iesHFOVals', 'isolHFOVals', 'allRippleVals', 'iesRippleVals', 'isolRippleVals', 'allFR_Vals', 'iesFR_Vals', 'isolFR_Vals'};

params.analysisPerZonesTablePath = analysisPerZonesTablePath;
params.selDetector = selDetector;

params.connTh = 75;
params.eiTh = 75;
params.outcomeTh = 49.9;


%% Analysis Loop
%% Set params
params.patsSel = 'improvedPatients'; {'allPatients', 'improvedPatients', 'nonImprovedPatients'};
params.freqBandConn = 'maxAllBands';

params.patsSel = 'allPatients';

params.normalization = '';{'', 'Scaled', 'Normalized'};

params.feature = 'OccRate';
[patAvgPre, patAvgPost, patAvgDiffOccRate] = getPatientAverageParcelActivityDiff(paths, params);
params.feature = 'Power';
[patAvgPre, patAvgPost, patAvgDiffPower] = getPatientAverageParcelActivityDiff(paths, params);
params.feature = 'Frequency';
[patAvgPre, patAvgPost, patAvgDiffFreq] = getPatientAverageParcelActivityDiff(paths, params);

diffActComp = rescale(patAvgDiffOccRate.iesHFOVals,0,1) + rescale(patAvgDiffPower.iesHFOVals,0,1);
outcome = patAvgPre.outcome;
scatter(outcome, diffActComp, 100, outcome, 'filled')

P = polyfit(outcome, diffActComp,1);
yfit = P(1)*outcome+P(2);
hold on;
plot(outcome,yfit,'r-.');

ylabel('iesHFO-\Delta');
xlabel("Outcome");
a = colorbar;
a.Label.String = 'Improvement(%)';
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w');
set(gca,'FontSize',20)

[p, h] = signrank(patAvgPre.iesHFOVals, patAvgPost.iesHFOVals);
[p, h] = ranksum(patAvgPre.iesHFOVals, patAvgPost.iesHFOVals);

diffAct = patAvgPre.iesHFOVals-patAvgPost.iesHFOVals;
diffAct = patAvgDiff.iesHFOVals;
outcome = patAvgPre.outcome;
scatter(diffAct, outcome, 25, outcome, 'filled')

[patAvgPre, patAvgPost] = getPatientAverageTables(groupTablePre, groupTablePre);

improvedPats_iesHFO_OccRate_Delta = groupTablePre.iesHFOVals - groupTablePost.iesHFOVals;
nrImprovedPatients = length(unique(groupTablePre.patName));
outcome = groupTablePre.outcome;
rftcConnect = groupTablePre.rftcElectroPhysioConnect;

params.patsSel = 'nonImprovedPatients';
params.feature = 'Power';
[groupTablePre, groupTablePost] = readTables(paths, params);
tableDelIdxsStruct = not(sameRFTC_BrainParcelChannsPerPat(groupTablePre));
tableDelIdxsConn = not(perctThValsPerPat(groupTablePre.patName, groupTablePre.rftcElectroPhysioConnect, params.connTh));
tableDelIdxs = tableDelIdxsStruct | tableDelIdxsConn;
tableDelIdxs = tableDelIdxsConn;
%tableDelIdxs = false(length(tableDelIdxs),1);

groupTablePre(tableDelIdxs, :) = [];
groupTablePost(tableDelIdxs, :) = [];
groupTablePre = removeRFTC_Channels(groupTablePre);
groupTablePost = removeRFTC_Channels(groupTablePost);
nonImpPats_iesHFO_OccRate_Delta = groupTablePre.iesHFOVals - groupTablePost.iesHFOVals;
nrNonImprovedPatients = length(unique(groupTablePre.patName));


[p, h] = ranksum(improvedPats_iesHFO_OccRate_Delta, nonImpPats_iesHFO_OccRate_Delta);
nrImprovedPatients
nrNonImprovedPatients
p

function [patAvgPreTable, patAvgPostTable, patAvgDiffTable] = getPatientAverageParcelActivityDiff(paths, params)

    %% Get RFTC Parcel Channels minus the RFTC channels themselves
    [groupTablePre, groupTablePost] = readTables(paths, params);
    notRFTC_ParcelChIdx = not(sameRFTC_BrainParcelChannsPerPat(groupTablePre));
    groupTablePre(notRFTC_ParcelChIdx, :) = [];
    groupTablePost(notRFTC_ParcelChIdx, :) = [];
    groupTablePre = removeRFTC_Channels(groupTablePre);
    groupTablePost = removeRFTC_Channels(groupTablePost);

    %%
    allPatNames = unique(groupTablePre.patName);
    nrPatients = length(allPatNames);
    columnNames = groupTablePre.Properties.VariableNames;
    
    for ci = 1:length(columnNames)
        column = columnNames{ci};
        patAvgPre.(column) = [];
        patAvgPost.(column) = [];
        patAvgDiff.(column) = [];
    end
    
    for pi = 1:nrPatients
        patSelIdx = ismember(groupTablePre.patName, allPatNames{pi});
        sum(patSelIdx)
        for ci = 1:length(columnNames)
            column = columnNames{ci};
            if not(ismember(ci, [1 3 7 8]))
                meanValPre = getValMean(groupTablePre.(column)(patSelIdx));
                meanValPost = getValMean(groupTablePost.(column)(patSelIdx));
                meanValDiff = groupTablePre.(column)(patSelIdx) - groupTablePost.(column)(patSelIdx);
                meanValDiff = getValMean(meanValDiff);
                patAvgPre.(column) = cat(1, patAvgPre.(column), meanValPre);
                patAvgPost.(column) = cat(1, patAvgPost.(column), meanValPost);
                patAvgDiff.(column) = cat(1, patAvgDiff.(column), meanValDiff);
            else
                varValPre = groupTablePre.(column)(patSelIdx);
                varValPost = groupTablePost.(column)(patSelIdx);
                patAvgPre.(column) = cat(1, patAvgPre.(column), varValPre(1));
                patAvgPost.(column) = cat(1, patAvgPost.(column), varValPost(1));
                patAvgDiff.(column) = cat(1, patAvgDiff.(column), varValPost(1));
            end
        end
        a = mean(groupTablePre.iesHFOVals(patSelIdx) - groupTablePost.iesHFOVals(patSelIdx));
        b = getValMean(groupTablePre.iesHFOVals(patSelIdx) - groupTablePost.iesHFOVals(patSelIdx));
        [a b]
    end
    patAvgPreTable = struct2table(patAvgPre);
    patAvgPostTable = struct2table(patAvgPost);
    patAvgDiffTable = struct2table(patAvgDiff);
end

function meanVec = getValMean(vec)
    selNoPosOutliers = vec <= (median(vec) + 1.5*std(vec));
    selNoNegOutliers = vec >= (median(vec) - 1.5*std(vec));
    selNoOutliers = selNoPosOutliers & selNoNegOutliers;
    meanVec = mean(vec(selNoOutliers));
    if isnan(meanVec)
        meanVec = mean(vec);
    end
end
% %% Analysis Loop
% for nmi = 1:length(normOptions)
%     %% Set params
%     params.normalization = normOptions{nmi};
%     params.patsSel = 'allPatients'; {'allPatients', 'improvedPatients', 'nonImprovedPatients'};
%     params.freqBandConn = 'maxAllBands';
% 
%     params.feature = 'OccRate';
%     [groupTablePre, groupTablePost] = readTables(paths, params);
%     tableDelIdxsStruct = not(sameRFTC_BrainParcelChannsPerPat(groupTablePre));
%     tableDelIdxsConn = not(perctThValsPerPat(groupTablePre.patName, groupTablePre.rftcElectroPhysioConnect, params.connTh));
%     tableDelIdxs = tableDelIdxsStruct | tableDelIdxsConn;
%     tableDelIdxs = tableDelIdxsStruct;
%     %tableDelIdxs = false(length(tableDelIdxs),1);
%     
%     groupTablePre(tableDelIdxs, :) = [];
%     groupTablePost(tableDelIdxs, :) = [];
%     groupTablePre = removeRFTC_Channels(groupTablePre);
%     groupTablePost = removeRFTC_Channels(groupTablePost);
%     iesHFO_OccRate_Delta = groupTablePre.iesHFOVals - groupTablePost.iesHFOVals;
%     
%     params.feature = 'Power';
%     [groupTablePre, groupTablePost] = readTables(paths, params);
%     groupTablePre(tableDelIdxs, :) = [];
%     groupTablePost(tableDelIdxs, :) = [];
%     groupTablePre = removeRFTC_Channels(groupTablePre);
%     groupTablePost = removeRFTC_Channels(groupTablePost);
%     iesHFO_Power_Delta = groupTablePre.iesHFOVals - groupTablePost.iesHFOVals;
%     
%     
%     params.feature = 'Frequency';
%     [groupTablePre, groupTablePost] = readTables(paths, params);
%     groupTablePre(tableDelIdxs, :) = [];
%     groupTablePost(tableDelIdxs, :) = [];
%     groupTablePre = removeRFTC_Channels(groupTablePre);
%     groupTablePost = removeRFTC_Channels(groupTablePost);
%     iesHFO_Freq_Delta = groupTablePre.iesHFOVals - groupTablePost.iesHFOVals;
% 
%     
%     outcome = groupTablePre.outcome;
%     rftcConnect = groupTablePre.rftcElectroPhysioConnect;
%     
%     sz = 25;
%     %scatter(iesHFO_OccRate_Delta, outcome, sz, outcome, 'filled')
%     xVal = iesHFO_OccRate_Delta+iesHFO_Power_Delta;
%     yVal = outcome;
%     scatter(xVal, outcome, sz, outcome, 'filled')
%     colorbar
%     
%     P = polyfit(xVal,yVal,1);
%     yfit = P(1)*xVal+P(2);
%     hold on;
%     plot(xVal,yfit,'r-.');
%         
% end

function biomarkerActvDiff = getAllChannels_ActivityDiff(params, groupTablePre, groupTablePost)
    groupTablePreIn = groupTablePre;
    groupTablePostIn = groupTablePost;
    groupTablePreOut = groupTablePre;
    groupTablePostOut = groupTablePost;
    
    %%
    groupTablePreIn = removeRFTC_Channels(groupTablePreIn);
    groupTablePostIn = removeRFTC_Channels(groupTablePostIn);
    groupTablePreOut = removeRFTC_Channels(groupTablePreOut);
    groupTablePostOut = removeRFTC_Channels(groupTablePostOut);
    
    %%
    biomarkerActvDiff = fillResultsStruct(params, groupTablePreIn, groupTablePostIn, groupTablePreOut, groupTablePostOut);
end

function biomarkerActvDiff = getRFTC_Channels_ActivityDiff(params, groupTablePre, groupTablePost)
    groupTablePreIn = groupTablePre;
    groupTablePostIn = groupTablePost;
    groupTablePreOut = groupTablePre;
    groupTablePostOut = groupTablePost;
    
    %%
    tableDelIdxs = groupTablePreIn.rftcVals < 1;
    groupTablePreIn(tableDelIdxs, :) = [];
    groupTablePostIn(tableDelIdxs, :) = [];
    
    %%
    tableDelIdxs = groupTablePreOut.rftcVals > 0;
    groupTablePreOut(tableDelIdxs, :) = [];
    groupTablePostOut(tableDelIdxs, :) = [];
    
    %%
    biomarkerActvDiff = fillResultsStruct(params, groupTablePreIn, groupTablePostIn, groupTablePreOut, groupTablePostOut);
end

function biomarkerActvDiff = getHighEI_Channels_ActivityDiff(params, groupTablePre, groupTablePost)
    biomarkerActvDiff.in = [];
    biomarkerActvDiff.out = [];
    groupTablePreIn = groupTablePre;
    groupTablePostIn = groupTablePost;
    groupTablePreOut = groupTablePre;
    groupTablePostOut = groupTablePost;
    
    %%
    tableDelIdxs = not(perctThValsPerPat(groupTablePre.patName, groupTablePre.eiVals, params.eiTh));
    groupTablePreIn(tableDelIdxs, :) = [];
    groupTablePostIn(tableDelIdxs, :) = [];
    
    %%
    tableDelIdxs = not(tableDelIdxs);
    groupTablePreOut(tableDelIdxs, :) = [];
    groupTablePostOut(tableDelIdxs, :) = [];
        
    %% 
    groupTablePreIn = removeRFTC_Channels(groupTablePreIn);
    groupTablePostIn = removeRFTC_Channels(groupTablePostIn);
    groupTablePreOut = removeRFTC_Channels(groupTablePreOut);
    groupTablePostOut = removeRFTC_Channels(groupTablePostOut);
    
    %%
    biomarkerActvDiff = fillResultsStruct(params, groupTablePreIn, groupTablePostIn, groupTablePreOut, groupTablePostOut);
end

function biomarkerActvDiff = getRFTC_ElectroPhysioConnected_Channels_ActivityDiff(params, groupTablePre, groupTablePost)
    biomarkerActvDiff.in = [];
    biomarkerActvDiff.out = [];
    groupTablePreIn = groupTablePre;
    groupTablePostIn = groupTablePost;
    groupTablePreOut = groupTablePre;
    groupTablePostOut = groupTablePost;
    
    %%
    tableDelIdxs = not(perctThValsPerPat(groupTablePre.patName, groupTablePre.rftcElectroPhysioConnect, params.connTh));
    groupTablePreIn(tableDelIdxs, :) = [];
    groupTablePostIn(tableDelIdxs, :) = [];
    
    %%
    tableDelIdxs = not(tableDelIdxs);
    groupTablePreOut(tableDelIdxs, :) = [];
    groupTablePostOut(tableDelIdxs, :) = [];
    
    %% 
    groupTablePreIn = removeRFTC_Channels(groupTablePreIn);
    groupTablePostIn = removeRFTC_Channels(groupTablePostIn);
    groupTablePreOut = removeRFTC_Channels(groupTablePreOut);
    groupTablePostOut = removeRFTC_Channels(groupTablePostOut);
    
    %%
    biomarkerActvDiff = fillResultsStruct(params, groupTablePreIn, groupTablePostIn, groupTablePreOut, groupTablePostOut);
end

function biomarkerActvDiff = getSameRFTC_Parcel_ActivityDiff(params, groupTablePre, groupTablePost)
    biomarkerActvDiff.in = [];
    biomarkerActvDiff.out = []; 
    groupTablePreIn = groupTablePre;
    groupTablePostIn = groupTablePost;
    groupTablePreOut = groupTablePre;
    groupTablePostOut = groupTablePost;
    
    %%
    tableDelIdxs = not(sameRFTC_BrainParcelChannsPerPat(groupTablePre));
    groupTablePreIn(tableDelIdxs, :) = [];
    groupTablePostIn(tableDelIdxs, :) = [];
    
    %%
    tableDelIdxs = not(tableDelIdxs);
    groupTablePreOut(tableDelIdxs, :) = [];
    groupTablePostOut(tableDelIdxs, :) = [];
    
    %% 
    groupTablePreIn = removeRFTC_Channels(groupTablePreIn);
    groupTablePostIn = removeRFTC_Channels(groupTablePostIn);
    groupTablePreOut = removeRFTC_Channels(groupTablePreOut);
    groupTablePostOut = removeRFTC_Channels(groupTablePostOut);
    
    %%
    biomarkerActvDiff = fillResultsStruct(params, groupTablePreIn, groupTablePostIn, groupTablePreOut, groupTablePostOut);
end

function biomarkerActvDiff = getSameRFTC_Lobe_ActivityDiff(params, groupTablePre, groupTablePost)
    biomarkerActvDiff.in = [];
    biomarkerActvDiff.out = [];
    groupTablePreIn = groupTablePre;
    groupTablePostIn = groupTablePost;
    groupTablePreOut = groupTablePre;
    groupTablePostOut = groupTablePost;
    
    %%
    tableDelIdxs = not(sameRFTC_BrainLobeChannsPerPat(groupTablePre));
    groupTablePreIn(tableDelIdxs, :) = [];
    groupTablePostIn(tableDelIdxs, :) = [];
    
    %%
    tableDelIdxs = not(tableDelIdxs);
    groupTablePreOut(tableDelIdxs, :) = [];
    groupTablePostOut(tableDelIdxs, :) = [];
    
    %% 
    groupTablePreIn = removeRFTC_Channels(groupTablePreIn);
    groupTablePostIn = removeRFTC_Channels(groupTablePostIn);
    groupTablePreOut = removeRFTC_Channels(groupTablePreOut);
    groupTablePostOut = removeRFTC_Channels(groupTablePostOut);
    
    %%
    biomarkerActvDiff = fillResultsStruct(params, groupTablePreIn, groupTablePostIn, groupTablePreOut, groupTablePostOut);;
end

function biomarkerActvDiff = getSameRFTC_Hemisphere_ActivityDiff(params, groupTablePre, groupTablePost)
    biomarkerActvDiff.in = [];
    biomarkerActvDiff.out = [];
    groupTablePreIn = groupTablePre;
    groupTablePostIn = groupTablePost;
    groupTablePreOut = groupTablePre;
    groupTablePostOut = groupTablePost;
    
    %%
    tableDelIdxs = not(sameRFTC_HemisphereChannsPerPat(groupTablePre));
    groupTablePreIn(tableDelIdxs, :) = [];
    groupTablePostIn(tableDelIdxs, :) = [];
    
    %%
    tableDelIdxs = not(tableDelIdxs);
    groupTablePreOut(tableDelIdxs, :) = [];
    groupTablePostOut(tableDelIdxs, :) = [];
    
    %% 
    groupTablePreIn = removeRFTC_Channels(groupTablePreIn);
    groupTablePostIn = removeRFTC_Channels(groupTablePostIn);
    groupTablePreOut = removeRFTC_Channels(groupTablePreOut);
    groupTablePostOut = removeRFTC_Channels(groupTablePostOut);
    
    %%
    biomarkerActvDiff = fillResultsStruct(params, groupTablePreIn, groupTablePostIn, groupTablePreOut, groupTablePostOut);
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
    
    %% read normalized EI
    preGroupTableFN = strcat(groupAnalysisTablesFilePath, 'Normalized','GroupAnalysis_ChannelCharacterization_Pre', params.selDetector, '.xls');
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