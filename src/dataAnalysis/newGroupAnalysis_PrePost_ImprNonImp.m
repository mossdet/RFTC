clc; clear all; close all;

paths = getFilesPaths();
biomarkersList =  {'HFO', 'iesHFO', 'IES'};
featuresList = {'rate','maxAmpl','variance','power'};
featuresList = {'rate','maxAmpl', 'power'};


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

% Tables to write
tableAllPatients_SpreadSheetName = strcat(analysisTablesPath, 'allPats_Pre_vs_Post_Analysis_', channSelStr, '_', normStr,'.xls');
tableImproved_vs_NonImproved_Pre_SpreadSheetName = strcat(analysisTablesPath, 'Improved_vs_nonImproved_Pre_Analysis_', channSelStr, '_', normStr,'.xls');
tableImproved_vs_NonImproved_Post_SpreadSheetName = strcat(analysisTablesPath, 'Improved_vs_nonImproved_Post_Analysis_', channSelStr, '_', normStr,'.xls');
tableImproved_vs_NonImproved_Delta_SpreadSheetName = strcat(analysisTablesPath, 'Improved_vs_nonImproved_Delta_Analysis_', channSelStr, '_', normStr,'.xls');
tableBiomCorrelAllPatients_SpreadSheetName = strcat(analysisTablesPath, 'allPatsBiomarkerZoneCorrelation', channSelStr, '_', normStr,'.xls');
delete(tableAllPatients_SpreadSheetName);
delete(tableImproved_vs_NonImproved_Pre_SpreadSheetName);
delete(tableImproved_vs_NonImproved_Delta_SpreadSheetName);
delete(tableImproved_vs_NonImproved_Post_SpreadSheetName);
delete(tableBiomCorrelAllPatients_SpreadSheetName);

for bm = biomarkersList
    biomarker = bm{1};
    allGroupAnalysisP = {};
    imprvd_vs_NonImprvd_Pre_P = {};
    imprvd_vs_NonImprvd_Post_P = {};
    imprvd_vs_NonImprvd_Delta_P = {};
    correlationAnalysis = {};

    groupTablePre = readtable(spreadSheetNamePre, 'Sheet', biomarker);
    groupTablePost = readtable(spreadSheetNamePost, 'Sheet', biomarker);
    featsNames = groupTablePre.Properties.VariableNames;
    
    if not(sum(strcmp(groupTablePre.chName, groupTablePost.chName)) == length(groupTablePost.chName))
        error(' newGroupAnalysis_Paired, Pre and Post EEG Channels don''t match')
    end

    % Improved vs non-Improved
    selImprvdPatGroup = groupTablePre.outcome >= outcomeTh;
    % Pre
    imprvdGroupTablePre = groupTablePre(selImprvdPatGroup,:);
    nonImprvdGroupTablePre = groupTablePre(not(selImprvdPatGroup),:);
    % Post
    imprvdGroupTablePost = groupTablePost(selImprvdPatGroup,:);
    nonImprvdGroupTablePost = groupTablePost(not(selImprvdPatGroup),:);
    % Delta
    imprvdDeltaPats = imprvdGroupTablePre;        
    cellA = table2cell(imprvdGroupTablePre(:,15:end));
    cellB = table2cell(imprvdGroupTablePost(:,15:end));
    imprvdDeltaPats(:,15:end) = cellfun(@minus,cellA,cellB,'Un',0);
    
    nonImprvsDeltaPats = nonImprvdGroupTablePre;
    cellA = table2cell(nonImprvdGroupTablePre(:,15:end));
    cellB = table2cell(nonImprvdGroupTablePost(:,15:end));
    nonImprvsDeltaPats(:,15:end) = cellfun(@minus,cellA,cellB,'Un',0);

    for ft = featuresList
        feature = ft{1};

        [analysisResP, analysisResNrChann] = zoneAnalysis_CompareDependentGroups(groupTablePre, 'groupTablePre', groupTablePost, 'groupTablePost', biomarker, feature);
        featCorrelAnal = biomarkerZonesCorrelationAnalysis(groupTablePre, 'groupTablePre', biomarker, feature, corrLimitP);

        % Tests                
        impVsNonImp_Pre_P_Res = zoneAnalysis_CompareIndependentGroups(imprvdGroupTablePre, 'imprvdGroupTablePre', nonImprvdGroupTablePre, 'nonImprvdGroupTablePre', biomarker, feature);
        impVsNonImp_Post_P_Res = zoneAnalysis_CompareIndependentGroups(imprvdGroupTablePost, 'imprvdGroupTablePost', nonImprvdGroupTablePost, 'nonImprvdGroupTablePost', biomarker, feature);
        impVsNonImp_Delta_P_Res = zoneAnalysis_CompareIndependentGroups(imprvdDeltaPats, 'imprvdDeltaPats', nonImprvsDeltaPats, 'nonImprvsDeltaPats', biomarker, feature);

        if isempty(allGroupAnalysisP)
            allGroupAnalysisP = cat(1, allGroupAnalysisP, analysisResP);
            imprvd_vs_NonImprvd_Pre_P = cat(1, imprvd_vs_NonImprvd_Pre_P, impVsNonImp_Pre_P_Res);
            imprvd_vs_NonImprvd_Post_P = cat(1, imprvd_vs_NonImprvd_Post_P, impVsNonImp_Post_P_Res);
            imprvd_vs_NonImprvd_Delta_P = cat(1, imprvd_vs_NonImprvd_Delta_P, impVsNonImp_Delta_P_Res);
            correlationAnalysis = cat(1, correlationAnalysis, featCorrelAnal);
        else
            allGroupAnalysisP = cat(1, allGroupAnalysisP, analysisResP(2:end,:));
            imprvd_vs_NonImprvd_Pre_P = cat(1, imprvd_vs_NonImprvd_Pre_P, impVsNonImp_Pre_P_Res(2:end,:));
            imprvd_vs_NonImprvd_Post_P = cat(1, imprvd_vs_NonImprvd_Post_P, impVsNonImp_Post_P_Res(2:end,:));
            imprvd_vs_NonImprvd_Delta_P = cat(1, imprvd_vs_NonImprvd_Delta_P, impVsNonImp_Delta_P_Res(2:end,:));
            correlationAnalysis = cat(1, correlationAnalysis, featCorrelAnal(2:end,:));
        end
    end
    % pcaMatrixPre = cell2mat(table2cell(groupTablePre(:,15:end)));
    % pcaMatrixPost = cell2mat(table2cell(groupTablePre(:,15:end)));
    % [coeffPre,scorePre,latentPre,tsquaredPre,explainedPre,muPre] = pca(pcaMatrixPre);
    % [coeffPost,scorePost,latentPost,tsquaredPost,explainedPost,muPost] = pca(pcaMatrixPost);

    plotSignificanceCircles(bm, imprvd_vs_NonImprvd_Delta_P);
    % % All Patients
    % tableAllPatients = cell2table(allGroupAnalysisP(2:end,:),'VariableNames', allGroupAnalysisP(1,:));
    % writetable(tableAllPatients, tableAllPatients_SpreadSheetName, 'Sheet', biomarker);   
    % 
    % % Improved Pre vs non-Improved Pre
    % tableImproved_vs_NonImproved_Pre = cell2table(imprvd_vs_NonImprvd_Pre_P(2:end,:),'VariableNames', imprvd_vs_NonImprvd_Pre_P(1,:));
    % writetable(tableImproved_vs_NonImproved_Pre, tableImproved_vs_NonImproved_Pre_SpreadSheetName, 'Sheet', biomarker); 
    % 
    % % Improved Post vs non-Improved Post
    % tableImproved_vs_NonImproved_Post = cell2table(imprvd_vs_NonImprvd_Post_P(2:end,:),'VariableNames', imprvd_vs_NonImprvd_Post_P(1,:));
    % writetable(tableImproved_vs_NonImproved_Post, tableImproved_vs_NonImproved_Post_SpreadSheetName, 'Sheet', biomarker);
    % 
    % % Improved Delta vs non-Improved Delta
    % tableImproved_vs_NonImproved_Delta = cell2table(imprvd_vs_NonImprvd_Delta_P(2:end,:),'VariableNames', imprvd_vs_NonImprvd_Delta_P(1,:));
    % writetable(tableImproved_vs_NonImproved_Delta, tableImproved_vs_NonImproved_Delta_SpreadSheetName, 'Sheet', biomarker);
    % 
    % % Biomarker Correlation
    % tableBiomCorrAllPats= cell2table(correlationAnalysis(2:end,:),'VariableNames', correlationAnalysis(1,:));
    % writetable(tableBiomCorrAllPats, tableBiomCorrelAllPatients_SpreadSheetName, 'Sheet', biomarker); 
end

function [analysisResP, analysisResNrChann] = zoneAnalysis_CompareDependentGroups(groupA, labelA, groupB, labelB, biomarker, feature)

    analysisResP = {'Biomarker', 'Feature', 'rftcSite', 'rftcStructure', 'rftcConnected', 'highEI', 'rftcLobe', 'rftcHemisphere'};
    analysisResNrChann = {'Biomarker', 'Feature', 'rftcSite', 'rftcStructure', 'rftcConnected', 'highEI', 'rftcLobe', 'rftcHemisphere'};
    pValEntry = {};
    channEntry = {};
    % rftcSite
    channSel = groupA.rftcSite == 1;
    pVal = signrank(groupA.(feature)(channSel), groupB.(feature)(channSel), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    channEntry = cat(2, channEntry, sum(channSel));
    % rftcStructure
    channSel = (groupA.rftcStruct == 1) & not(groupA.rftcSite == 1);
    channSelB = (groupB.rftcStruct == 1) & not(groupB.rftcSite == 1);
    if sum(channSel == channSelB) ~= length(groupA.rftcHemis)
        error('Difference in pre post data');
    end
    pVal = signrank(groupA.(feature)(channSel), groupB.(feature)(channSel), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    channEntry = cat(2, channEntry, sum(channSel));
    % rftcConnected
    channSel = (groupA.rftcConn == 1) & not(groupA.rftcSite == 1);
    channSelB = (groupB.rftcConn == 1) & not(groupB.rftcSite == 1);
    if sum(channSel == channSelB) ~= length(groupA.rftcHemis)
        error('Difference in pre post data');
    end
    pVal = signrank(groupA.(feature)(channSel), groupB.(feature)(channSel), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    channEntry = cat(2, channEntry, sum(channSel));
    % highEI
    channSel = (groupA.highEI == 1) & not(groupA.rftcSite == 1);
    channSelB = (groupB.highEI == 1) & not(groupB.rftcSite == 1);
    if sum(channSel == channSelB) ~= length(groupA.rftcHemis)
        error('Difference in pre post data');
    end
    pVal = signrank(groupA.(feature)(channSel), groupB.(feature)(channSel), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    channEntry = cat(2, channEntry, sum(channSel));
    % rftcLobe
    channSel = (groupA.rftcLobe == 1) & not(groupA.rftcSite == 1) & not(groupA.rftcStruct == 1);
    channSelB = (groupB.rftcLobe == 1) & not(groupB.rftcSite == 1) & not(groupB.rftcStruct == 1);
    if sum(channSel == channSelB) ~= length(groupA.rftcHemis)
        error('Difference in pre post data');
    end
    pVal = signrank(groupA.(feature)(channSel), groupB.(feature)(channSel), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    channEntry = cat(2, channEntry, sum(channSel));
    % rftcHemis
    channSel = (groupA.rftcHemis == 1) & not(groupA.rftcSite == 1) & not(groupA.rftcStruct == 1) & not(groupA.rftcLobe == 1);
    channSelB = (groupB.rftcHemis == 1) & not(groupB.rftcSite == 1) & not(groupB.rftcStruct == 1) & not(groupB.rftcLobe == 1);
    if sum(channSel == channSelB) ~= length(groupA.rftcHemis)
        error('Difference in pre post data');
    end
    pVal = signrank(groupA.(feature)(channSel), groupB.(feature)(channSel), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    channEntry = cat(2, channEntry, sum(channSel));
    
    analysisResP = cat(1, analysisResP, cat(2, {biomarker}, {feature}, pValEntry));
    analysisResNrChann = cat(1, analysisResNrChann, cat(2, {biomarker}, {feature}, channEntry));
end

function analysisResP = zoneAnalysis_CompareIndependentGroups(groupA, labelA, groupB, labelB, biomarker, feature)

    analysisResP = {'Biomarker', 'Feature', 'rftcSite', 'rftcStructure', 'rftcConnected', 'highEI', 'rftcLobe', 'rftcHemisphere'};
    pValEntry = {};
    % rftcSite
    channSelA = groupA.rftcSite == 1;
    channSelB = groupB.rftcSite == 1;
    pVal = ranksum(groupA.(feature)(channSelA), groupB.(feature)(channSelB), 'tail', 'right'); % 'both' | 'right' | 'left'
    pValEntry = cat(2, pValEntry, pVal);
    % rftcStruct
    channSelA = (groupA.rftcStruct == 1) & not(groupA.rftcSite == 1);
    channSelB = (groupB.rftcStruct == 1) & not(groupB.rftcSite == 1);
    pVal = ranksum(groupA.(feature)(channSelA), groupB.(feature)(channSelB), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    % rftcConn
    channSelA = (groupA.rftcConn == 1) & not(groupA.rftcSite == 1);
    channSelB = (groupB.rftcConn == 1) & not(groupB.rftcSite == 1);
    pVal = ranksum(groupA.(feature)(channSelA), groupB.(feature)(channSelB), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    % highEI
    channSelA = (groupA.highEI == 1) & not(groupA.rftcSite == 1);
    channSelB = (groupB.highEI == 1) & not(groupB.rftcSite == 1);
    pVal = ranksum(groupA.(feature)(channSelA), groupB.(feature)(channSelB), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    % rftcLobe
    channSelA = (groupA.rftcLobe == 1) & not(groupA.rftcSite == 1) & not(groupA.rftcStruct == 1);
    channSelB = (groupB.rftcLobe == 1) & not(groupB.rftcSite == 1) & not(groupB.rftcStruct == 1);
    pVal = ranksum(groupA.(feature)(channSelA), groupB.(feature)(channSelB), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    % rftcHemis
    channSelA = (groupA.rftcHemis == 1) & not(groupA.rftcSite == 1) & not(groupA.rftcStruct == 1) & not(groupA.rftcLobe == 1);
    channSelB = (groupB.rftcHemis == 1) & not(groupB.rftcSite == 1) & not(groupB.rftcStruct == 1) & not(groupB.rftcLobe == 1);
    pVal = ranksum(groupA.(feature)(channSelA), groupB.(feature)(channSelB), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    
    analysisResP = cat(1, analysisResP, cat(2, {biomarker}, {feature}, pValEntry));
end

function analysisResP = zoneAnalysis_CompareIndependentGroupsDeprecated(groupA, labelA, groupB, labelB, biomarker, feature)

    analysisResP = {'Biomarker', 'Feature', 'rftcSite', 'rftcStructure', 'rftcConnected', 'highEI', 'rftcLobe', 'rftcHemisphere'};
    pValEntry = {};
    % rftcSite
    channSelA = groupA.rftcSite == 1;
    channSelB = groupB.rftcSite == 1;
    pVal = ranksum(groupA.(feature)(channSelA), groupB.(feature)(channSelB), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    % rftcStructure
    channSelA = (groupA.rftcStruct == 1);
    channSelA = channSelA & not(groupA.rftcSite == 1);
    channSelB = (groupB.rftcStruct == 1);
    channSelB = channSelB & not(groupB.rftcSite == 1);
    pVal = ranksum(groupA.(feature)(channSelA), groupB.(feature)(channSelB), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    % rftcConnected
    channSelA = (groupA.rftcConn == 1);
    channSelA = channSelA & not(groupA.rftcSite == 1);
    channSelB = (groupB.rftcConn == 1);
    channSelB = channSelB & not(groupB.rftcSite == 1);    
    pVal = ranksum(groupA.(feature)(channSelA), groupB.(feature)(channSelB), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    % highEI
    channSelA = (groupA.highEI == 1);
    channSelA = channSelA & not(groupA.rftcSite == 1);
    channSelB = (groupB.highEI == 1);
    channSelB = channSelB & not(groupB.rftcSite == 1);
    pVal = ranksum(groupA.(feature)(channSelA), groupB.(feature)(channSelB), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    % rftcLobe
    channSelA = (groupA.rftcLobe == 1);
    channSelA = channSelA & not(groupA.rftcSite == 1) & not(groupA.rftcStruct == 1);
    channSelB = (groupB.rftcLobe == 1);
    channSelB = channSelB & not(groupB.rftcSite == 1) & not(groupB.rftcStruct == 1);
    pVal = ranksum(groupA.(feature)(channSelA), groupB.(feature)(channSelB), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    % rftcHemis
    channSelA = (groupA.rftcHemis == 1);
    channSelA = channSelA & not(groupA.rftcSite == 1) & not(groupA.rftcStruct == 1) & not(groupA.rftcLobe == 1);
    channSelB = (groupB.rftcHemis == 1);
    channSelB = channSelB & not(groupB.rftcSite == 1) & not(groupB.rftcStruct == 1) & not(groupB.rftcLobe == 1);
    pVal = ranksum(groupA.(feature)(channSelA), groupB.(feature)(channSelB), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    
    analysisResP = cat(1, analysisResP, cat(2, {biomarker}, {feature}, pValEntry));
end

function analysisCorr = biomarkerZonesCorrelationAnalysis(group, label, biomarker, feature, corrLimitP)

    analysisCorr = {'Biomarker', 'Feature', 'rftcSite', 'rftcStructure', 'rftcConnected', 'highEI', 'rftcLobe', 'rftcHemisphere'};
    corrValEntry = {};
    
    % rftcSite
    channSel = group.rftcSite == 1;
    [rho,pval] = corr(channSel, group.(feature), 'Type','Spearman');
    if pval > corrLimitP
        rho = 0;
    end
    corrValEntry = cat(2, corrValEntry, rho);    
    
    % rftcStructure
    channSel = (group.rftcStruct == 1);
    channSel = channSel & not(group.rftcSite == 1);
    [rho,pval] = corr(channSel, group.(feature), 'Type','Spearman');
    if pval > corrLimitP
        rho = 0;
    end
    corrValEntry = cat(2, corrValEntry, rho);
    
    % rftcConnected
    [rho,pval] = corr(group.rftcElectroPhysioConnect, group.(feature), 'Type','Spearman');
    if pval > corrLimitP
        rho = 0;
    end
    corrValEntry = cat(2, corrValEntry, rho);
    
    % highEI
    [rho,pval] = corr(group.eiVals, group.(feature), 'Type','Spearman');
    if pval > corrLimitP
        rho = 0;
    end
    corrValEntry = cat(2, corrValEntry, rho);
    
    % rftcLobe
    channSel = (group.rftcLobe == 1);
    channSel = channSel & not(group.rftcSite == 1);
    [rho,pval] = corr(channSel, group.(feature), 'Type','Spearman');
    if pval > corrLimitP
        rho = 0;
    end
    corrValEntry = cat(2, corrValEntry, rho);
    
    % rftcHemis
    channSel = (group.rftcHemis == 1);
    channSel = channSel & not(group.rftcSite == 1);
    [rho,pval] = corr(channSel, group.(feature), 'Type','Spearman');
    if pval > corrLimitP
        rho = 0;
    end
    corrValEntry = cat(2, corrValEntry, rho);
    
    analysisCorr = cat(1, analysisCorr, cat(2, {biomarker}, {feature}, corrValEntry));
end
