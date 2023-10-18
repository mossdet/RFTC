clc; clear all; close all;

paths = getFilesPaths();
biomarkersList =  {'HFO', 'iesHFO', 'IES'};
featuresList = {'rate', 'maxAmpl','variance','power'};

patTablesPath = strcat(paths.workspacePath, 'CharacterizationTables\FeatureCharacterizationTables\');
groupTablesPath = strcat(paths.workspacePath, 'CharacterizationTables\GroupCharacterizationTablesAvg\');
analysisTablesPath = strcat(paths.workspacePath, 'AnalysisResults\Avg\');mkdir(analysisTablesPath);

normStr = ''; % {'', '_Normalized'};
channSelStr = 'FlexK'; % {'Threshold', 'FlexK', 'Cluster2K'};

prctlTh = 75;
outcomeTh = 90;
corrLimitP = 0.001;

% Tables to read
spreadSheetNamePre = strcat(groupTablesPath, 'GroupCharacterization_', channSelStr, '_Pre', normStr,'.xls');
spreadSheetNamePost = strcat(groupTablesPath, 'GroupCharacterization_', channSelStr, '_Post', normStr,'.xls');

resultsTableP = {};
resultsTableNrChanns = {};
for bm = biomarkersList
    preData = [];
    postData= [];
    biomarker = bm{1}
    resultsTableBiomarkerNrChanns = {};
    resultsTableBiomarkerP = {};
    for ft = featuresList
        feature = ft{1};

        groupTablePre = readtable(spreadSheetNamePre, 'Sheet', biomarker);
        groupTablePost = readtable(spreadSheetNamePost, 'Sheet', biomarker);
        featsNames = groupTablePre.Properties.VariableNames;

        if not(sum(strcmp(groupTablePre.chName, groupTablePost.chName)) == length(groupTablePost.chName))
            error(' newGroupAnalysis_Paired, Pre and Post EEG Channels don''t match')
        end

        [analysisResP, analysisResNrChann] = zoneAnalysis_CompareDependentGroups(groupTablePre, 'Pre', groupTablePost, 'Post', biomarker, feature);
        
        if isempty(resultsTableBiomarkerP)
            resultsTableBiomarkerP = cat(1, resultsTableBiomarkerP, analysisResP);
            resultsTableBiomarkerNrChanns = cat(1, resultsTableBiomarkerNrChanns, analysisResNrChann);
        else
            resultsTableBiomarkerP = cat(1, resultsTableBiomarkerP, analysisResP(2,:));
            resultsTableBiomarkerNrChanns = cat(1, resultsTableBiomarkerNrChanns, analysisResNrChann(2,:));
        end
    end
    plotSignificanceCircles(bm, resultsTableBiomarkerP)

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
    channSel = (groupA.rftcStruct == 1);
    channSel = channSel & not(groupA.rftcSite == 1);
    pVal = signrank(groupA.(feature)(channSel), groupB.(feature)(channSel), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    channEntry = cat(2, channEntry, sum(channSel));
    % rftcConnected
    channSel = (groupA.rftcConn == 1);
    channSel = channSel & not(groupA.rftcSite == 1);
    pVal = signrank(groupA.(feature)(channSel), groupB.(feature)(channSel), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    channEntry = cat(2, channEntry, sum(channSel));
    % highEI
    channSel = (groupA.highEI == 1);
    channSel = channSel & not(groupA.rftcSite == 1);
    pVal = signrank(groupA.(feature)(channSel), groupB.(feature)(channSel), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    channEntry = cat(2, channEntry, sum(channSel));
    % rftcLobe
    channSel = (groupA.rftcLobe == 1);
    channSel = channSel & not(groupA.rftcSite == 1) & not(groupA.rftcStruct == 1);
    pVal = signrank(groupA.(feature)(channSel), groupB.(feature)(channSel), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    channEntry = cat(2, channEntry, sum(channSel));
    % rftcHemis
    channSel = (groupA.rftcHemis == 1);
    channSel = channSel & not(groupA.rftcSite == 1) & not(groupA.rftcStruct == 1) & not(groupA.rftcLobe == 1);
    pVal = signrank(groupA.(feature)(channSel), groupB.(feature)(channSel), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    channEntry = cat(2, channEntry, sum(channSel));
    
    analysisResP = cat(1, analysisResP, cat(2, {biomarker}, {feature}, pValEntry));
    analysisResNrChann = cat(1, analysisResNrChann, cat(2, {biomarker}, {feature}, channEntry));
end
