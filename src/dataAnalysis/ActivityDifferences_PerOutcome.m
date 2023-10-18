clc; clear all; close all;
paths = getFilesPaths();
preFiles  = getPreFiles();%getAllFiles, getPreFiles, getPostFiles
postFiles  = getPostFiles();%getAllFiles, getPreFiles, getPostFiles
outcomePredictionAnalysisPath = 'F:\ForschungsProjekte\RFTC\MATLAB\DetectHFO\OutcomePrediction\';

selDetector = 'MOSSDET'; 'Delphos'; 'MOSSDET'; 'MOSSDET_Depurated';
outcomePredictionAnalysisPath = strcat(outcomePredictionAnalysisPath, selDetector, '\');
features = {'OccRate', 'Power', 'Frequency'};
normOptions = {'', 'Normalized', '', 'Scaled'}; %'Normalized'; % {'', 'Scaled', 'Normalized'};
freqBandConnList = { 'maxAllBands'};
thList = 0;
patsSelList = {'improvedPatients', 'nonImprovedPatients'};
biomarkersList =  {'allHFOVals', 'iesHFOVals', 'isolHFOVals', 'allRippleVals', 'iesRippleVals', 'isolRippleVals', 'allFR_Vals', 'iesFR_Vals', 'isolFR_Vals', 'iaVals'};

params.outcomePredictionAnalysisPath = outcomePredictionAnalysisPath;
params.selDetector = selDetector;
params.connTh = 75;
params.eiTh = 75;
params.outcomeTh = 49;

biomarkerName=biomarkersList';

%% Analysis Loop
for fi = 1:length(features)
    rftc = zeros(length(biomarkersList),1);
    highEI = zeros(length(biomarkersList),1);
    rftcConnected = zeros(length(biomarkersList),1);
    rftcStructure = zeros(length(biomarkersList),1);
    rftcLobe = zeros(length(biomarkersList),1);
    rftcHemisphere = zeros(length(biomarkersList),1);
    predictTablePVal = table(biomarkerName, rftc, highEI, rftcConnected, rftcStructure, rftcLobe, rftcHemisphere);
    for bmi = 1:length(biomarkersList)
        %% Set params
        params.th = 0;
        params.normalization = '';
        params.feature = features{fi};
        params.analysisSubType = features{fi};
        params.freqBandConn = 'maxAllBands';
        params.biomarker = biomarkersList{bmi};

        params.patsSel = 'improvedPatients';
        [groupTablePre, groupTablePost] = readTables(paths, params);
        %[groupTablePre, groupTablePost] = normalizePerPatient(groupTablePre, groupTablePost);
        
        
        %% Start analyses 
        biomarkerActvDiff_RFTC_Improved = getRFTC_Channels_ActivityDiff(params, groupTablePre, groupTablePost);
        biomarkerActvDiff_highEI_Improved = getHighEI_Channels_ActivityDiff(params, groupTablePre, groupTablePost);
        biomarkerActvDiff_rftcConnected_Improved = getRFTC_ElectroPhysioConnected_Channels_ActivityDiff(params, groupTablePre, groupTablePost);
        sameRFTC_Structure_ActvDiff_Improved = getSameRFTC_Structure_ActivityDiff(params, groupTablePre, groupTablePost);
        sameRFTC_Lobe_ActvDiff_Improved = getSameRFTC_Lobe_ActivityDiff(params, groupTablePre, groupTablePost);
        sameRFTC_Hemisphere_ActvDiff_Improved = getSameRFTC_Hemisphere_ActivityDiff(params, groupTablePre, groupTablePost);

        params.patsSel = 'nonImprovedPatients';
        [groupTablePre, groupTablePost] = readTables(paths, params);
        %% Start analyses 
        biomarkerActvDiff_RFTC_NonImproved = getRFTC_Channels_ActivityDiff(params, groupTablePre, groupTablePost);
        biomarkerActvDiff_highEI_NonImproved  = getHighEI_Channels_ActivityDiff(params, groupTablePre, groupTablePost);
        biomarkerActvDiff_rftcConnected_NonImproved  = getRFTC_ElectroPhysioConnected_Channels_ActivityDiff(params, groupTablePre, groupTablePost);
        sameRFTC_Structure_ActvDiff_NonImproved  = getSameRFTC_Structure_ActivityDiff(params, groupTablePre, groupTablePost);
        sameRFTC_Lobe_ActvDiff_NonImproved  = getSameRFTC_Lobe_ActivityDiff(params, groupTablePre, groupTablePost);
        sameRFTC_Hemisphere_ActvDiff_NonImproved  = getSameRFTC_Hemisphere_ActivityDiff(params, groupTablePre, groupTablePost);

        %%        
        predictTableDiff.(biomarkersList{bmi}).rftc.imp = biomarkerActvDiff_RFTC_Improved.in.diff;
        predictTableDiff.(biomarkersList{bmi}).rftc.nonImp= biomarkerActvDiff_RFTC_NonImproved.in.diff;
        predictTableDiff.(biomarkersList{bmi}).highEI.imp = biomarkerActvDiff_highEI_Improved.in.diff;
        predictTableDiff.(biomarkersList{bmi}).highEI.nonImp = biomarkerActvDiff_highEI_NonImproved.in.diff;
        predictTableDiff.(biomarkersList{bmi}).rftcConnected.imp = biomarkerActvDiff_rftcConnected_Improved.in.diff;
        predictTableDiff.(biomarkersList{bmi}).rftcConnected.nonImp = biomarkerActvDiff_rftcConnected_NonImproved.in.diff;
        predictTableDiff.(biomarkersList{bmi}).rftcStructure.imp = sameRFTC_Structure_ActvDiff_Improved.in.diff;
        predictTableDiff.(biomarkersList{bmi}).rftcStructure.nonImp = sameRFTC_Structure_ActvDiff_NonImproved.in.diff;
        predictTableDiff.(biomarkersList{bmi}).rftcLobe.imp = sameRFTC_Lobe_ActvDiff_Improved.in.diff;
        predictTableDiff.(biomarkersList{bmi}).rftcLobe.nonImp = sameRFTC_Lobe_ActvDiff_NonImproved.in.diff;
        predictTableDiff.(biomarkersList{bmi}).rftcHemisphere.imp = sameRFTC_Hemisphere_ActvDiff_Improved.in.diff;
        predictTableDiff.(biomarkersList{bmi}).rftcHemisphere.nonImp = sameRFTC_Hemisphere_ActvDiff_NonImproved.in.diff;

        %%
        [pRFTC, h] = ranksum(biomarkerActvDiff_RFTC_Improved.in.diff, biomarkerActvDiff_RFTC_NonImproved.in.diff);
        [pEI, h] = ranksum(biomarkerActvDiff_highEI_Improved.in.diff, biomarkerActvDiff_highEI_NonImproved.in.diff);
        [pRFTC_Connect, h] = ranksum(biomarkerActvDiff_rftcConnected_Improved.in.diff, biomarkerActvDiff_rftcConnected_NonImproved.in.diff);
        [pRFTC_Structure, h] = ranksum(sameRFTC_Structure_ActvDiff_Improved.in.diff, sameRFTC_Structure_ActvDiff_NonImproved.in.diff);
        [pRFTC_Lobe, h] = ranksum(sameRFTC_Lobe_ActvDiff_Improved.in.diff, sameRFTC_Lobe_ActvDiff_NonImproved.in.diff);
        [pRFTC_Hemisphere, h] = ranksum(sameRFTC_Hemisphere_ActvDiff_Improved.in.diff, sameRFTC_Hemisphere_ActvDiff_NonImproved.in.diff);

        predictTablePVal.biomarkerName{bmi} = biomarkersList{bmi};
        predictTablePVal.rftc(bmi) = pRFTC;
        predictTablePVal.highEI(bmi) = pEI;
        predictTablePVal.rftcConnected(bmi) = pRFTC_Connect;
        predictTablePVal.rftcStructure(bmi) = pRFTC_Structure;
        predictTablePVal.rftcLobe(bmi) = pRFTC_Lobe;
        predictTablePVal.rftcHemisphere(bmi) = pRFTC_Hemisphere;
    end
    %plotOutcomeDifferences_Boxplot(predictTableDiff, predictTablePVal, params);
    mkdir(params.outcomePredictionAnalysisPath);
    tableFN = strcat(params.outcomePredictionAnalysisPath, 'Outcome_Prediction.xls');
    writetable(predictTablePVal, tableFN, 'Sheet', features{fi})
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

function biomarkerActvDiff = getHighEI_Channels_ActivityDiff(params, groupTablePre, groupTablePost)
    biomarkerActvDiff.in = [];
    biomarkerActvDiff.out = [];
    groupTablePreIn = groupTablePre;
    groupTablePostIn = groupTablePost;
    groupTablePreOut = groupTablePre;
    groupTablePostOut = groupTablePost;
    
    %%
    %tableDelIdxs = not(groupTablePre.eiVals > median(groupTablePre.eiVals) + std(groupTablePre.eiVals));
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
    %tableDelIdxs = not(groupTablePre.rftcElectroPhysioConnect > median(groupTablePre.rftcElectroPhysioConnect) + std(groupTablePre.rftcElectroPhysioConnect));
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

function biomarkerActvDiff = getSameRFTC_Structure_ActivityDiff(params, groupTablePre, groupTablePost)
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
    rftcStructureIdx = sameRFTC_BrainParcelChannsPerPat(groupTablePre);
    tableDelIdxs = tableDelIdxs | rftcStructureIdx;
    
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
    rftcLobeIdxs = sameRFTC_BrainLobeChannsPerPat(groupTablePre);
    rftcStructureIdx = sameRFTC_BrainParcelChannsPerPat(groupTablePre);
    tableDelIdxs = tableDelIdxs | rftcLobeIdxs | rftcStructureIdx;
    
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

        patKeepIdx = patVals >= median(patVals)+0.75*std(patVals);
        
        allPatsKeepIdxs(patSelIdx) = patKeepIdx;        
    end
    allPatsKeepIdxs = logical(allPatsKeepIdxs);
end

function [groupTablePre, groupTablePost] = normalizePerPatient(groupTablePre, groupTablePost)

    groupTableAll = [];

    allPatNames = unique(groupTablePre.patName);
    nrPatients = length(allPatNames);
    columnNames = groupTablePre.Properties.VariableNames;

    for pi = 1:nrPatients
        patSelIdx = ismember(groupTablePre.patName, allPatNames{pi});
        patSelIdxPost = ismember(groupTablePost.patName, allPatNames{pi});
        if sum(patSelIdx & patSelIdxPost) ~= sum(patSelIdx)
            stop = 1;
        end
        for ci = 1:length(columnNames)
            column = columnNames{ci};
            columnPatValsPre = groupTablePre.(column)(patSelIdx);
            columnPatValsPost = groupTablePost.(column)(patSelIdx);
            columnPatValsAll = [columnPatValsPre columnPatValsPost];
            if length(columnPatValsPre) ~= length(columnPatValsPost)
                stop = 1;
            end

            if ci >= 9
                groupTablePre.(column)(patSelIdx) = (columnPatValsPre-mean(columnPatValsAll))/std(columnPatValsAll);
                groupTablePost.(column)(patSelIdx) = (columnPatValsPost-mean(columnPatValsAll))/std(columnPatValsAll);
            end
        end
    end

end

        
function [groupTablePre, groupTablePost] = readTables(paths, params)
    %% Read tables
    groupAnalysisTablesFilePath = strcat(paths.workspacePath, 'GroupAnalysis_ChannelCharacterizationTables\', params.selDetector, '\');
    preGroupTableFN = strcat(groupAnalysisTablesFilePath, params.normalization, 'GroupAnalysis_ChannelCharacterization_Pre', params.selDetector, '.xls');
    postGroupTableFN = strcat(groupAnalysisTablesFilePath, params.normalization, 'GroupAnalysis_ChannelCharacterization_Post', params.selDetector, '.xls');
    groupTablePre = readtable(preGroupTableFN, 'Sheet', params.feature);
    groupTablePost = readtable(postGroupTableFN, 'Sheet', params.feature);
    
    %groupTablePre = truncateOutliers(groupTablePre);
    %groupTablePost = truncateOutliers(groupTablePost);
    
    %% read normalized EI
    preGroupTableFN = strcat(groupAnalysisTablesFilePath, 'Normalized','GroupAnalysis_ChannelCharacterization_Pre', params.selDetector, '.xls');
    postGroupTableFN = strcat(groupAnalysisTablesFilePath, 'Normalized','GroupAnalysis_ChannelCharacterization_Post', params.selDetector, '.xls');
    groupTablePreEI_correct = readtable(preGroupTableFN, 'Sheet', params.feature);
    groupTablePostEI_correct = readtable(postGroupTableFN, 'Sheet', params.feature);
    groupTablePre.eiVals = groupTablePreEI_correct.eiVals;
    groupTablePost.eiVals = groupTablePostEI_correct.eiVals;
    
    %% read normalized rftcElectroPhysioConnect
    preGroupTableNormFN = strcat(groupAnalysisTablesFilePath, 'Normalized', 'GroupAnalysis_ChannelCharacterization_Pre', params.selDetector,  '_', params.freqBandConn, '.xls');
    postGroupTableNormFN = strcat(groupAnalysisTablesFilePath, 'Normalized', 'GroupAnalysis_ChannelCharacterization_Post', params.selDetector, '_', params.freqBandConn,'.xls');
    groupTablePre_RFTCconn = readtable(preGroupTableNormFN, 'Sheet', params.feature);
    groupTablePost_RFTCconn = readtable(postGroupTableNormFN, 'Sheet', params.feature);
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