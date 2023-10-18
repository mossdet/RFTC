clc; clear all; close all;

paths = getFilesPaths();

patTablesPath = strcat(paths.workspacePath, 'CharacterizationTables\FeatureCharacterizationTables\');
groupTablesPath = strcat(paths.workspacePath, 'CharacterizationTables\GroupCharacterizationTablesAvg\');
analysisTablesPath = strcat(paths.workspacePath, 'AnalysisResults\Avg\');mkdir(analysisTablesPath);

figuresOutputPath = 'F:\ForschungsProjekte\RFTC\Paper\Figures\TempOutput\';mkdir(figuresOutputPath);

biomarkersList =  {'HFO', 'iesHFO', 'IES'};
featuresList = {'rate', 'maxAmpl', 'variance', 'power'};

channSelStr = 'FlexK';
normStr = ''; %'', '_Normalized'};
zoneFormation = ''; % '', '_NonExclusive';
prctlTh = 75;
outcomeTh = 90;
corrLimitP = 0.001;

% Tables to read
spreadSheetNamePre = strcat(groupTablesPath, 'GroupCharacterization_', channSelStr, '_Pre', normStr,'.xls');
spreadSheetNamePost = strcat(groupTablesPath, 'GroupCharacterization_', channSelStr, '_Post', normStr,'.xls');
        
for bm = biomarkersList
    correlationAnalysis = {};
    biomarker = bm{1};
    groupTablePre = readtable(spreadSheetNamePre, 'Sheet', biomarker);
    groupTablePost = readtable(spreadSheetNamePost, 'Sheet', biomarker);
    featsNames = groupTablePre.Properties.VariableNames;

    if not(sum(strcmp(groupTablePre.chName, groupTablePost.chName)) == length(groupTablePost.chName))
        error(' newGroupAnalysis_Paired, Pre and Post EEG Channels don''t match')
    end
    
    plotZonesCorrelationCircles(groupTablePre, zoneFormation, figuresOutputPath);
    
end

function plotZonesCorrelationCircles(groupTablePre, zoneFormation, figuresOutputPath)
    close all;
    patNames = groupTablePre.patName; 
    if not(isempty(zoneFormation))
        % get Zones
        rftcChanns = logical(groupTablePre.rftcSite);
        highEI_Channs = logical(groupTablePre.highEI);
        highConn_Channs = logical(groupTablePre.rftcConn);
        rftcStructure_Channs = logical(groupTablePre.rftcStruct);
        rftcLobe_Channs = logical(groupTablePre.rftcLobe);
        rftcHemisphere_Channs = logical(groupTablePre.rftcHemis);
    else    
        rftcChanns = logical(groupTablePre.rftcSite);
        highEI_Channs = logical(groupTablePre.highEI) & not(rftcChanns);
        highConn_Channs = logical(groupTablePre.rftcConn) & not(rftcChanns);
        rftcStructure_Channs = logical(groupTablePre.rftcStruct) & not(rftcChanns);
        rftcLobe_Channs = logical(groupTablePre.rftcLobe) & not(groupTablePre.rftcSite) & not(groupTablePre.rftcStruct);
        rftcHemisphere_Channs = logical(groupTablePre.rftcHemis)  & not(groupTablePre.rftcSite) & not(groupTablePre.rftcStruct) & not(groupTablePre.rftcLobe);
    end
    
    %% Plot Channel Data
    [rftcChanns, newOrder] = sort(rftcChanns,'descend');
    patNames = patNames(newOrder);
    highEI_Channs = highEI_Channs(newOrder);
    highConn_Channs = highConn_Channs(newOrder);
    rftcStructure_Channs = rftcStructure_Channs(newOrder);
    rftcLobe_Channs = rftcLobe_Channs(newOrder);
    rftcHemisphere_Channs = rftcHemisphere_Channs(newOrder);

    allZones = [rftcChanns highConn_Channs highEI_Channs rftcStructure_Channs rftcLobe_Channs rftcHemisphere_Channs]';    
    labels = {'rftcSite', 'rftcConnected', 'highEI', 'rftcStructure', 'rftcLobe', 'rftcHemisphere'};
    mccMatrix = zeros(length(labels), length(labels));
    for fi = 1:size(allZones,1)
        for si = 1:size(allZones,1)
            [mcc, pVal] = getMCC_withSignificance(allZones(fi,:), allZones(si,:));
            mcc = round(mcc*100);
            mccMatrix(fi,si) = mcc;
        end
    end
    hm = heatmap(labels, labels, mccMatrix, 'FontSize', 18, 'Colormap',summer);
    hm = heatmap(mccMatrix);
    hm.YDisplayData = labels;
end

function [mcc,pVal] = getChannelsOverlap(vecA, vecB)

    mcc = sum(vecA & vecB)/sum(vecA);
    pVal = 0;

end
function [mcc, pVal] = getMCC_withSignificance(vecA, vecB)
    nrSamples = length(vecA);
    nrTests = 1000;
    mccRndDistr = zeros(nrTests,1);
    rng('default')
    for n = 1:nrTests
        vecATest = randi([0 1],nrSamples,1);
        vecBTest = randi([0 1],nrSamples,1);
        mccTest = getMCC(vecATest, vecBTest);
        mccRndDistr(n) = mccTest;
    end

    mcc = getMCC(vecA, vecB);
    pValHigh = (1-sum(mcc>=mccRndDistr)/nrTests);
    pValLow = (1-sum(mcc<=mccRndDistr)/nrTests);
    pVal = min([pValHigh pValLow]);

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
        [cid, centroids, sumd, D] = kmeans(patValsRand , 3, 'Start','cluster', 'Distance', 'sqeuclidean', 'MaxIter', 1000); %
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