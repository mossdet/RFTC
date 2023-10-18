clc; clear all; close all;
paths = getFilesPaths();
preFiles  = getPreFiles();%getAllFiles, getPreFiles, getPostFiles
postFiles  = getPostFiles();%getAllFiles, getPreFiles, getPostFiles
outcomePredictionAnalysisPath = 'F:\ForschungsProjekte\RFTC\MATLAB\DetectHFO\OutcomePrediction\';

selDetector = 'MOSSDET'; 'Delphos'; 'MOSSDET'; 'MOSSDET_Depurated';
outcomePredictionAnalysisPath = strcat(outcomePredictionAnalysisPath, selDetector, '\');
features = {'OccRate', 'Power', 'Frequency'};
features = {'OccRate', 'Power'};

normOptions = {''}; %'Normalized'; % {'', 'Scaled', 'Normalized'};
freqBandConnList = { 'maxAllBands'};
thList = 0;
patsSelList = {'allPatients', 'improvedPatients', 'nonImprovedPatients'};
biomarkersList =  {'allHFOVals', 'iesHFOVals', 'isolHFOVals', 'allRippleVals', 'iesRippleVals', 'isolRippleVals', 'allFR_Vals', 'iesFR_Vals', 'isolFR_Vals', 'iaVals'};
biomarkersList =  {'iesHFOVals', 'allHFOVals', 'iaVals'};

params.outcomePredictionAnalysisPath = outcomePredictionAnalysisPath;
params.selDetector = selDetector;
params.connTh = 80;
params.eiTh = 80;
params.outcomeTh = 90;

biomarkerName=biomarkersList';

predictionAllBiomarkers = cell(length(biomarkersList),6);
%% Analysis Loop
for bmi = 1:length(biomarkersList)
    zonesDiffData = [];
    for fi = 1:length(features)
        %% Set params
        params.th = 0;
        params.normalization = '';
        params.feature = features{fi};
        params.analysisSubType = features{fi};
        params.freqBandConn = 'maxAllBands';
        params.biomarker = biomarkersList{bmi};

        params.patsSel = 'allPatients';
        [groupTablePre, groupTablePost] = readTables(paths, params);
        [groupTablePre, groupTablePost] = normalizePerPatient(groupTablePre, groupTablePost);
        
        patNames = unique(groupTablePre.patName);
        nrPats = length(patNames);
        
        %% Outcome
        perPatient_Outcome = zeros(nrPats,1);
        for pi = 1:nrPats
            patName = patNames{pi};
            patSel = strcmp(groupTablePre.patName, patName);
            outcomes = groupTablePre.outcome(patSel);
            perPatient_Outcome(pi) = mean(outcomes);
        end
        
        %% 
        perPatient_RFTC_Delta = zeros(nrPats,1);
        rftcConn_Delta = zeros(nrPats,1);
        HighEI_Delta = zeros(nrPats,1);
        structure_Delta = zeros(nrPats,1);
        lobe_Delta = zeros(nrPats,1);
        hemisphere_Delta = zeros(nrPats,1);

        for pi = 1:nrPats
            pi
            patName = patNames{pi}
            patSel = strcmp(groupTablePre.patName, patName);
            patTablePre = groupTablePre(patSel, :);
            patTablePost = groupTablePost(patSel, :);
            
            perPatient_RFTC_Delta(pi) = zoneSummaryMetric(getRFTC_Channels_ActivityDiff(params, patTablePre, patTablePost).in.diff);
            rftcConn_Delta(pi) = zoneSummaryMetric(getRFTC_ElectroPhysioConnected_Channels_ActivityDiff(params, patTablePre, patTablePost).in.diff);
            HighEI_Delta(pi) = zoneSummaryMetric(getHighEI_Channels_ActivityDiff(params, patTablePre, patTablePost).in.diff);
            
            structure_Delta(pi) = zoneSummaryMetric(getSameRFTC_Structure_ActivityDiff(params, patTablePre, patTablePost).in.diff);
            lobe_Delta(pi) = zoneSummaryMetric(getSameRFTC_Lobe_ActivityDiff(params, patTablePre, patTablePost).in.diff);
            hemisphere_Delta(pi) = zoneSummaryMetric(getSameRFTC_Hemisphere_ActivityDiff(params, patTablePre, patTablePost).in.diff);
        end
        
        zonesDiffData = cat(2, zonesDiffData, perPatient_RFTC_Delta, rftcConn_Delta, HighEI_Delta, structure_Delta, lobe_Delta, hemisphere_Delta);
    end
    
    pwFeats = 6;
    if strcmp(biomarkersList{bmi}, 'iesHFOVals')
        zonesSelOccRate = [];
        zonesSelPow = [1 2 4]+pwFeats;
    elseif strcmp(biomarkersList{bmi}, 'allHFOVals')
        zonesSelOccRate = [];
        zonesSelPow = [1 2 4]+pwFeats;
    elseif strcmp(biomarkersList{bmi}, 'iaVals')
        zonesSelOccRate = [];
        zonesSelPow = [1 2 4]+pwFeats;
    end
    zonesSel = [zonesSelOccRate zonesSelPow];
    
    [maxSensitivity, maxSpecificity, maxAccuracy, maxFone, maxMCC, maxCentroids, distanceMt, zonesSel] = clusterPatients(patNames, perPatient_Outcome, zonesDiffData, zonesSel);
    [maxSensitivityPros, maxSpecificityPros, maxAccuracyPros, maxFonePros, maxMCCPros] = prospectivePrediction(patNames, perPatient_Outcome, zonesDiffData, zonesSel, maxCentroids, distanceMt);

    predictionAllBiomarkers{bmi,1} = biomarkersList{bmi};
    predictionAllBiomarkers{bmi,2} = maxSensitivityPros;
    predictionAllBiomarkers{bmi,3} = maxSpecificityPros;
    predictionAllBiomarkers{bmi,4} = maxAccuracyPros;
    predictionAllBiomarkers{bmi,5} = maxFonePros;
    predictionAllBiomarkers{bmi,6} = maxMCCPros;
    
    biomarkersList{bmi}
    maxMCCPros
end

plotOutcomePredictionResult_AllBiomarkers(predictionAllBiomarkers);

function summarizedMetric = zoneSummaryMetric(channelsVals)
%     summarizedMetric = prctile(channelsVals, 90);
%     summarizedMetric = max(channelsVals);
%     summarizedMetric = mean(channelsVals);
%     summarizedMetric = median(channelsVals);
    summarizedMetric = mean(channelsVals);
end

function mcc = getMCC(vecA, vecB)
    vecA = logical(vecA);
    vecB = logical(vecB);

    tp = sum(vecA & vecB);
    tn = sum(not(vecA) & not(vecB));
    fp = sum(not(vecA) & vecB);
    fn = sum(vecA & not(vecB));
    
    mccA = (tp * tn) - (fp * fn);
    mccB = ((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));
    mcc = mccA / sqrt(mccB);
end

function getPerPatientAvgDiff(biomarkerActvDiff)
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
    if sum(rftcLobeIdxs) == length(rftcLobeIdxs)
        rftcLobeIdxs = false(length(rftcLobeIdxs),1);
        stop = 1;
    end
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
%     biomarkerActvDiff.in.diff = groupTablePreIn.(params.biomarker);
%     biomarkerActvDiff.out.diff = groupTablePreOut.(params.biomarker);
    
    %biomarkerActvDiff.in.prePostP = signrank(groupTablePreIn.(params.biomarker), groupTablePostIn.(params.biomarker), 'tail', 'right');
    %biomarkerActvDiff.out.prePostP = signrank(groupTablePreOut.(params.biomarker), groupTablePostOut.(params.biomarker), 'tail', 'right');
    
    biomarkerActvDiff.in.nrChanns = length(groupTablePreIn.patName);
    biomarkerActvDiff.out.nrChanns = length(groupTablePreOut.patName);
    
    biomarkerActvDiff.in.nrPats = length(unique(groupTablePreIn.patName));
    biomarkerActvDiff.out.nrPats = length(unique(groupTablePreOut.patName));

    biomarkerActvDiff.in.outcome = groupTablePreIn.outcome;
    biomarkerActvDiff.out.outcome = groupTablePreOut.outcome;
    
    %[p, h] = ranksum(biomarkerActvDiff.in.diff, biomarkerActvDiff.out.diff);
    %biomarkerActvDiff.inOutP = p;
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
        patValsIdxs = 1:length(patVals);
        
        pi = randperm(length(patVals));
        patValsRand = patVals(pi);
        patValsIdxsRand = patValsIdxs(pi);
        [cid, centroids, sumd, D] = kmeans(patValsRand , 2, 'Start','cluster', 'Distance', 'sqeuclidean', 'MaxIter', 1000); %
        [maxVal, maxIdx] = max([mean(patValsRand(cid==1)), mean(patValsRand(cid==2)), mean(patValsRand(cid==3))]);
        patKeepIdxA = cid == maxIdx;        
        [sortedVec, sortingIdxs] = sort(patValsIdxsRand);
        patValsRand = patValsRand(sortingIdxs);
        patKeepIdxA = patKeepIdxA(sortingIdxs);
        
        interQR = iqr(patVals);
        quart1 = prctile(patVals, 25);
        patKeepIdxB = patVals >= prctile(patVals, th);
        patKeepIdxC = patVals >= median(patVals)+0.75*std(patVals);
        
        [min(patVals(patKeepIdxA)) prctile(patVals, th)  median(patVals)+0.75*std(patVals)]
        [sum(patKeepIdxA) sum(patKeepIdxB) sum(patKeepIdxC)]
        
        patKeepIdx = patKeepIdxC;
        
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
        tableDelIdxsIn = groupTablePre.outcome < params.outcomeTh;
        groupTablePre(tableDelIdxsIn, :) = [];
        groupTablePost(tableDelIdxsIn, :) = [];        
    elseif strcmp(params.patsSel, 'nonImprovedPatients')
        tableDelIdxsIn = groupTablePre.outcome >= params.outcomeTh;
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