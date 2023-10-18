clc; clear all; close all;

paths = getFilesPaths();
biomarkersList =  {'HFO', 'iesHFO', 'IES'};
biomarkersList =  {'HFO','iesHFO'};
featuresList = {'rate','maxAmpl','variance','power'};
niceFeaturesList = {'Occ.Rate','Max.Ampl.','Variance','Power'};
zonesNames = {'rftcSite', 'rftcStructure', 'rftcConnected', 'highEI', 'rftcLobe', 'rftcHemisphere'};

patTablesPath = strcat(paths.workspacePath, 'CharacterizationTables\FeatureCharacterizationTables\');
groupTablesPath = strcat(paths.workspacePath, 'CharacterizationTables\GroupCharacterizationTablesAvg\');
analysisTablesPath = strcat(paths.workspacePath, 'AnalysisResults\Avg\');mkdir(analysisTablesPath);

normStr = '';
channSelStr = 'FlexK';

prctlTh = 75;
outcomeTh = 90;
corrLimitP = 0.001;

% Tables to read
spreadSheetNamePre = strcat(groupTablesPath, 'GroupCharacterization_', channSelStr, '_Pre', normStr,'.xls');
spreadSheetNamePost = strcat(groupTablesPath, 'GroupCharacterization_', channSelStr, '_Post', normStr,'.xls');

% % Tables to write
% tableImproved_vs_NonImproved_Delta_SpreadSheetName = strcat(analysisTablesPath, 'Improved_vs_nonImproved_Delta_Analysis_', channSelStr, '_', normStr,'.xls');
% delete(tableImproved_vs_NonImproved_Delta_SpreadSheetName);

histEst = false;
nrHistPts = 25;
featsValsAllBiomarkers = [];
for bmi = 1:length(biomarkersList)
    biomarkerPatsCharacterization = [];
    nrZones = length(zonesNames);
    nrFeatures= length(featuresList);
    biomarker = biomarkersList{bmi};


    groupTablePre = readtable(spreadSheetNamePre, 'Sheet', biomarker);
    groupTablePost = readtable(spreadSheetNamePost, 'Sheet', biomarker);
    featsNames = groupTablePre.Properties.VariableNames;
    patNames = unique(groupTablePre.patName);

    if not(sum(strcmp(groupTablePre.chName, groupTablePost.chName)) == length(groupTablePost.chName))
        error(' newGroupAnalysis_Paired, Pre and Post EEG Channels don''t match')
    end

    for pi = 1:length(patNames)
        patName = patNames{pi};
        if strcmp(patName, 'GRE_2017_MESd_Inter_Sleep_2048')
            continue;
        end
        patSel = ismember(groupTablePre.patName, patName);
        patTablePre = groupTablePre(patSel, :);
        patTablePost = groupTablePost(patSel, :);
        singlePatCharacterization = {};

        posOutcome = sum(patTablePre.outcome >= outcomeTh) > 0;

        %% Get Feature Delta
        patDelta = patTablePre;
        cellA = table2cell(patTablePre(:,15:end));
        cellB = table2cell(patTablePost(:,15:end));
        patDelta(:,15:end) = cellfun(@minus,cellA,cellB,'Un',0);
        patDelta = normalizeFeatures(patDelta);

        for fi = 1:length(featuresList)
            feature = featuresList{fi};
            niceFeature = niceFeaturesList{fi};        

            % rftcSite
            channSel = patDelta.rftcSite == 1;
            siteDeltaVals = getZoneSummaryMetric(patDelta.(feature)(channSel));
            
            % rftcStructure
            channSel = (patDelta.rftcStruct == 1) & not(patDelta.rftcSite == 1);
            structDeltaVals = getZoneSummaryMetric(patDelta.(feature)(channSel));
            
            % rftcConnected
            channSel = (patDelta.rftcConn == 1) & not(patDelta.rftcSite == 1);
            connectedDeltaVals = getZoneSummaryMetric(patDelta.(feature)(channSel));
    
            % highEI
            channSel = (patDelta.highEI == 1) & not(patDelta.rftcSite == 1);
            eiDeltaVals = getZoneSummaryMetric(patDelta.(feature)(channSel));
    
            % rftcLobe
            channSel = (patDelta.rftcLobe == 1) & not(patDelta.rftcSite == 1) & not(patDelta.rftcStruct == 1);
            lobeDeltaVals = getZoneSummaryMetric(patDelta.(feature)(channSel));
    
            % rftcHemisphere
            channSel = (patDelta.rftcHemis == 1) & not(patDelta.rftcSite == 1) & not(patDelta.rftcStruct == 1) & not(patDelta.rftcLobe == 1);
            hemisDeltaVals = getZoneSummaryMetric(patDelta.(feature)(channSel));

            singlePatCharacterization = cat(1, singlePatCharacterization, {patName, biomarker, feature, posOutcome, siteDeltaVals, structDeltaVals, connectedDeltaVals, eiDeltaVals, lobeDeltaVals, hemisDeltaVals});
        end
        %relevantFeatVals = getBiomarkerRelevantFeatValsAll(singlePatCharacterization);
        relevantFeatVals = getBiomarkerRelevantFeatValsSelected(singlePatCharacterization);
        biomarkerPatsCharacterization = cat(1, biomarkerPatsCharacterization, [relevantFeatVals, posOutcome]);
    end
    outcomeVals = biomarkerPatsCharacterization(:,end)>0;
    featsVals = biomarkerPatsCharacterization(:,1:end-1);
    featsValsAllBiomarkers = [featsValsAllBiomarkers, featsVals];
    kappaMax = 0;
    mccMax = 0;
    cMax = [];
    bestPrediction = [];
    kappaCum = predictThermocoagulationOutcomeSVM(featsVals, outcomeVals);
    for ci = 1:100
        [idx, C] = kmeans(featsVals,2, 'Distance','sqeuclidean', 'MaxIter', 1000, 'Start','cluster'); % sqeuclidean, correlation
        idx = (idx-1) > 0;
        kappaScore = getKappa(outcomeVals, idx);
        mcc = getMCC(outcomeVals, idx);
        if kappaScore > kappaMax
            kappaMax  = kappaScore;
            cMax = C;
            bestPrediction = idx;
            mccMax = mcc;
        end
    end
    {biomarker, kappaMax, mccMax}

    plotFeatsOK = false;
    if plotFeatsOK
        if strcmp(biomarker, 'HFO')
            subplot(2,2,1)
            plot(featsVals(outcomeVals,1), 'og', 'MarkerFaceColor', 'g', 'MarkerSize',20); hold on;
            plot(featsVals(~outcomeVals,1), 'or', 'MarkerFaceColor', 'r', 'MarkerSize',20); hold on;
            legend('Positive Outcome', 'Negative Outcome')
            ylabel('rftcStructure Max Amplitude \Delta')
            fontsize(18,"points");
            set(gca,'XTick',[])
    
            subplot(2,2,2)
            plot(featsVals(outcomeVals,2), 'og', 'MarkerFaceColor', 'g', 'MarkerSize',20); hold on;
            plot(featsVals(~outcomeVals,2), 'or', 'MarkerFaceColor', 'r', 'MarkerSize',20); hold on;
            legend('Positive Outcome', 'Negative Outcome')
            ylabel('rftcConnected Max Amplitude \Delta')
            set(gca,'XTick',[])
            fontsize(18,"points");
    
            subplot(2,2,3)
            plot(featsVals(outcomeVals,3), 'og', 'MarkerFaceColor', 'g', 'MarkerSize',20); hold on;
            plot(featsVals(~outcomeVals,3), 'or', 'MarkerFaceColor', 'r', 'MarkerSize',20); hold on;
            legend('Positive Outcome', 'Negative Outcome')
            ylabel('rftcConnected Variance \Delta')
            set(gca,'XTick',[])
            fontsize(18,"points");
    
            subplot(2,2,4)
            plot(featsVals(outcomeVals,4), 'og', 'MarkerFaceColor', 'g', 'MarkerSize',20); hold on;
            plot(featsVals(~outcomeVals,4), 'or', 'MarkerFaceColor', 'r', 'MarkerSize',20); hold on;
            legend('Positive Outcome', 'Negative Outcome')
            ylabel('rftcConnected Power \Delta')
            set(gca,'XTick',[])
            fontsize(18,"points");
    
            set(gcf, 'Position', get(0, 'Screensize'), 'color','w');
            sgtitle(biomarker);
        elseif strcmp(biomarker, 'iesHFO')
            subplot(2,3,1)
            plot(featsVals(outcomeVals,1), 'og', 'MarkerFaceColor', 'g', 'MarkerSize',20); hold on;
            plot(featsVals(~outcomeVals,1), 'or', 'MarkerFaceColor', 'r', 'MarkerSize',20); hold on;
            legend('Positive Outcome', 'Negative Outcome')
            ylabel('rftcStructure Max Amplitude \Delta')
            set(gca,'XTick',[])
            fontsize(18,"points");
    
            subplot(2,3,2)
            plot(featsVals(outcomeVals,2), 'og', 'MarkerFaceColor', 'g', 'MarkerSize',20); hold on;
            plot(featsVals(~outcomeVals,2), 'or', 'MarkerFaceColor', 'r', 'MarkerSize',20); hold on;
            legend('Positive Outcome', 'Negative Outcome')
            ylabel('rftcStructure Variance \Delta')
            set(gca,'XTick',[])
            fontsize(18,"points");
    
            subplot(2,3,3)
            plot(featsVals(outcomeVals,3), 'og', 'MarkerFaceColor', 'g', 'MarkerSize',20); hold on;
            plot(featsVals(~outcomeVals,3), 'or', 'MarkerFaceColor', 'r', 'MarkerSize',20); hold on;
            legend('Positive Outcome', 'Negative Outcome')
            ylabel('rftcStructure Power \Delta')
            set(gca,'XTick',[])
            fontsize(18,"points");
    
            subplot(2,3,4)
            plot(featsVals(outcomeVals,4), 'og', 'MarkerFaceColor', 'g', 'MarkerSize',20); hold on;
            plot(featsVals(~outcomeVals,4), 'or', 'MarkerFaceColor', 'r', 'MarkerSize',20); hold on;
            legend('Positive Outcome', 'Negative Outcome')
            ylabel('rftcConnected Max Amplitude \Delta')
            set(gca,'XTick',[])
            fontsize(18,"points");
    
            subplot(2,3,5)
            plot(featsVals(outcomeVals,5), 'og', 'MarkerFaceColor', 'g', 'MarkerSize',20); hold on;
            plot(featsVals(~outcomeVals,5), 'or', 'MarkerFaceColor', 'r', 'MarkerSize',20); hold on;
            legend('Positive Outcome', 'Negative Outcome')
            ylabel('rftcConnected Variance \Delta')
            set(gca,'XTick',[])
            fontsize(18,"points");
    
            subplot(2,3,6)
            plot(featsVals(outcomeVals,6), 'og', 'MarkerFaceColor', 'g', 'MarkerSize',20); hold on;
            plot(featsVals(~outcomeVals,6), 'or', 'MarkerFaceColor', 'r', 'MarkerSize',20); hold on;
            legend('Positive Outcome', 'Negative Outcome')
            ylabel('rftcConnected Power \Delta')
            set(gca,'XTick',[])
            fontsize(18,"points");
    
            set(gcf, 'Position', get(0, 'Screensize'), 'color','w');
            sgtitle(biomarker); 
        end
    end 
end

kappaMax = 0;
mccMax = 0;
cMax = [];
bestPrediction = [];
kappaCum = predictThermocoagulationOutcomeSVM(featsValsAllBiomarkers, outcomeVals);
for ci = 1:100
    [idx, C] = kmeans(featsValsAllBiomarkers,2, 'Distance','sqeuclidean', 'MaxIter', 1000, 'Start','cluster'); % sqeuclidean, correlation
    idx = (idx-1) > 0;
    kappaScore = getKappa(outcomeVals, idx);
    mcc = getMCC(outcomeVals, idx);
    if kappaScore > kappaMax
        kappaMax  = kappaScore;
        cMax = C;
        bestPrediction = idx;
        mccMax = mcc;
    end
end
{biomarker, kappaMax, mccMax}

function normPatDelta = normalizeFeatures(patDelta)
    normPatDelta = patDelta;
    featVals = normPatDelta(:,15:end);
    featVals = table2cell(featVals);
    featVals = cell2mat(featVals);
    nrFeats = size(featVals,2);
    for fi = 1:nrFeats
        featVals(:,fi) = (featVals(:,fi)-mean(featVals(:,fi)))/std(featVals(:,fi));
    end
    normPatDelta(:,15:end) = num2cell(featVals);
end

function summaryMetric = getZoneSummaryMetric(featVector)
    summaryMetric = mean(featVector);
end

function relevantFeatVals = getBiomarkerRelevantFeatValsAll(singlePatCharacterization)
    biomarker = singlePatCharacterization{1,2};
    relevantFeatVals = [];
    for zi = 5:10
        relevantFeatVals = cat(2, relevantFeatVals,...
        [singlePatCharacterization{1, zi},...
        singlePatCharacterization{2, zi},...
        singlePatCharacterization{3, zi},...
        singlePatCharacterization{4, zi}]);
    end
end

function relevantFeatVals = getBiomarkerRelevantFeatValsSelected(singlePatCharacterization)
    biomarker = singlePatCharacterization{1,2};
    relevantFeatVals = [];
    % rftcSite, rftcStructure, rftcConnected, highEI, rftcLobe, rftcHemis
    if strcmp(biomarker, 'HFO')
        % relevantFeatVals = [...
        %     singlePatCharacterization{2, 6},...
        %     singlePatCharacterization{2, 7},...
        %     singlePatCharacterization{3, 7},...
        %     singlePatCharacterization{4, 7},...
        %     ];
        relevantFeatVals = [...
            singlePatCharacterization{2, 6},...
            singlePatCharacterization{2, 7},...
            singlePatCharacterization{4, 7},...
            ];
    elseif strcmp(biomarker, 'iesHFO')
        % relevantFeatVals = [...
        % singlePatCharacterization{2, 6},...
        % singlePatCharacterization{3, 6},...
        % singlePatCharacterization{4, 6},...
        % singlePatCharacterization{2, 7},...
        % singlePatCharacterization{3, 7},...
        % singlePatCharacterization{4, 7},...
        % ];
        relevantFeatVals = [...
        singlePatCharacterization{2, 6},...
        singlePatCharacterization{4, 6},...
        singlePatCharacterization{2, 7},...
        singlePatCharacterization{4, 7},...
        ];
    elseif  strcmp(biomarker, 'IES')
    end
end

function kappa = getKappa(vecA, vecB)
    tp = 0; fp = 0; tn = 0; fn = 0;

    tp = tp + sum(vecA & vecB);
    fp = fp + sum(not(vecA) & vecB);
    tn = tn + sum(not(vecA) & not(vecB));
    fn = fn + sum(vecA & not(vecB));
    kappaA = 100 * (2 * (tp*tn - fn*fp));
    kappaB = 100 * ((tp+fp) * (fp+tn) + (tp+fn) * (fn+tn));
    if kappaB <= 0.01
        kappa = 0;
    else
        kappa = kappaA/kappaB;
    end
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
