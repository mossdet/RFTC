clc; clear all; close all;

paths = getFilesPaths();
patTablesPath = strcat(paths.workspacePath, 'CharacterizationTables\FeatureCharacterizationTables\');
groupTablesPath = strcat(paths.workspacePath, 'CharacterizationTables\GroupCharacterizationTablesAvg\');
analysisTablesPath = strcat(paths.workspacePath, 'AnalysisResults\Avg\');mkdir(analysisTablesPath);

biomarkersList =  {'HFO', 'iesHFO', 'IES'};
featuresList = {'rate','duration','maxAmpl','sumAmpl','variance','BPLL','power','sumPower','mobility','complexity','PBRatio','wvltRipplePower','wvltFastRiplePower','peaks','zeroCrossPerEOI','spectCentroid','spectPeak','deltaOccRate','deltaAmpl','deltaPow','deltaVar','thetaOccRate','thetaAmpl','thetaPow','thetaVar','alphaOccRate','alphaAmpl','alphaPow','alphaVar','betaOccRate','betaAmpl','betaPow','betaVar','gammaOccRate','gammaAmpl','gammaPow','gammaVar','highGammaOccRate','highGammaAmpl','highGammaPow','highGammaVar','rippleOccRate','rippleAmpl','ripplePow','rippleVar','fRippleOccRate','fRippleAmpl','fRipplePow','fRippleVar'};
featuresList = {'rate','maxAmpl','variance','power'};
normStrList = {''}; % '', '_Normalized'
channSelStrList = {'FlexK'}; %{'Threshold', 'FlexK', 'Cluster2K'};

prctlTh = 75;
outcomeTh = 90;
corrLimitP = 0.001;


for normIdx = 1:length(normStrList)
    for chSelIdx = 1:length(channSelStrList)
        normStr = normStrList{normIdx};
        channSelStr = channSelStrList{chSelIdx};
        
        % Tables to read
        spreadSheetNamePre = strcat(groupTablesPath, 'GroupCharacterization_', channSelStr, '_Pre', normStr,'.xls');
        spreadSheetNamePost = strcat(groupTablesPath, 'GroupCharacterization_', channSelStr, '_Post', normStr,'.xls');

        % Tables to write
        tableAllPatients_NonExc_SpreadSheetName = strcat(analysisTablesPath, 'allPats_Pre_vs_Post_Analysis_', channSelStr, '_NonExc', normStr,'.xls');
        tableAllPatients_Exc_SpreadSheetName = strcat(analysisTablesPath, 'allPats_Pre_vs_Post_Analysis_', channSelStr, '_Exc', normStr,'.xls');
        delete(tableAllPatients_NonExc_SpreadSheetName);
        delete(tableAllPatients_Exc_SpreadSheetName);
        
        for bm = biomarkersList
            allGroupAnalysisP_NonExc = {};
            allGroupAnalysisP_Exc = {};

            biomarker = bm{1};
            groupTablePre = readtable(spreadSheetNamePre, 'Sheet', biomarker);
            groupTablePost = readtable(spreadSheetNamePost, 'Sheet', biomarker);
            featsNames = groupTablePre.Properties.VariableNames;
            featuresList = featsNames(find(ismember(featsNames, 'rate')):end);

            if not(sum(strcmp(groupTablePre.chName, groupTablePost.chName)) == length(groupTablePost.chName))
                error(' newGroupAnalysis_Paired, Pre and Post EEG Channels don''t match')
            end

            for ft = featuresList
                feature = ft{1};

                [analysisResP_NonExc, analysisResNrChann_NonExc] = zoneAnalysis_CompareDependentGroups_NonExclusive(groupTablePre, 'groupTablePre', groupTablePost, 'groupTablePost', biomarker, feature);
                [analysisResP_Exc, analysisResNrChann_Exc] = zoneAnalysis_CompareDependentGroups_Exclusive(groupTablePre, 'groupTablePre', groupTablePost, 'groupTablePost', biomarker, feature);

                if isempty(allGroupAnalysisP_NonExc)
                    allGroupAnalysisP_NonExc = cat(1, allGroupAnalysisP_NonExc, analysisResP_NonExc);
                    allGroupAnalysisP_Exc = cat(1, allGroupAnalysisP_Exc, analysisResP_Exc);
                else
                    allGroupAnalysisP_NonExc = cat(1, allGroupAnalysisP_NonExc, analysisResP_NonExc(2:end,:));
                    allGroupAnalysisP_Exc = cat(1, allGroupAnalysisP_Exc, analysisResP_Exc(2:end,:));
                end
            end
            % All Patients
            tableAllPatients_NonExc = cell2table(allGroupAnalysisP_NonExc(2:end,:),'VariableNames', allGroupAnalysisP_NonExc(1,:));
            writetable(tableAllPatients_NonExc, tableAllPatients_NonExc_SpreadSheetName, 'Sheet', biomarker);
            tableAllPatients_Exc = cell2table(allGroupAnalysisP_Exc(2:end,:),'VariableNames', allGroupAnalysisP_Exc(1,:));
            writetable(tableAllPatients_Exc, tableAllPatients_Exc_SpreadSheetName, 'Sheet', biomarker);

        end
    end
end

function [analysisResP, analysisResNrChann] = zoneAnalysis_CompareDependentGroups_NonExclusive(groupA, labelA, groupB, labelB, biomarker, feature)

    analysisResP = {'Biomarker', 'Feature', 'rftcSite', 'rftcStructure', 'rftcConnected', 'highEI', 'rftcLobe', 'rftcHemisphere'};
    analysisResNrChann = {'Biomarker', 'Feature', 'rftcSite', 'rftcStructure', 'rftcConnected', 'highEI', 'rftcLobe', 'rftcHemisphere'};
    pValEntry = {};
    channEntry = {};
    % rftcSite
    rfctSiteFlags = groupA.rftcSite > 0;
    channSel = rfctSiteFlags;
    pVal = signrank(groupA.(feature)(channSel), groupB.(feature)(channSel), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    channEntry = cat(2, channEntry, sum(channSel));
    % rftcStructure
    channSel = (groupA.rftcStruct > 0);
    pVal = signrank(groupA.(feature)(channSel), groupB.(feature)(channSel), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    channEntry = cat(2, channEntry, sum(channSel));
    % rftcConnected
    channSel = (groupA.rftcConn > 0);
    pVal = signrank(groupA.(feature)(channSel), groupB.(feature)(channSel), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    channEntry = cat(2, channEntry, sum(channSel));
    % highEI
    channSel = (groupA.highEI > 0);
    pVal = signrank(groupA.(feature)(channSel), groupB.(feature)(channSel), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    channEntry = cat(2, channEntry, sum(channSel));
    % rftcLobe
    channSel = (groupA.rftcLobe > 0);
    pVal = signrank(groupA.(feature)(channSel), groupB.(feature)(channSel), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    channEntry = cat(2, channEntry, sum(channSel));
    % rftcHemis
    channSel = (groupA.rftcHemis == 1);
    pVal = signrank(groupA.(feature)(channSel), groupB.(feature)(channSel), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    channEntry = cat(2, channEntry, sum(channSel));
    
    analysisResP = cat(1, analysisResP, cat(2, {biomarker}, {feature}, pValEntry));
    analysisResNrChann = cat(1, analysisResNrChann, cat(2, {biomarker}, {feature}, channEntry));
end

function [analysisResP, analysisResNrChann] = zoneAnalysis_CompareDependentGroups_Exclusive(groupA, labelA, groupB, labelB, biomarker, feature)

    analysisResP = {'Biomarker', 'Feature', 'rftcSite', 'rftcStructure', 'rftcConnected', 'highEI', 'rftcLobe', 'rftcHemisphere'};
    analysisResNrChann = {'Biomarker', 'Feature', 'rftcSite', 'rftcStructure', 'rftcConnected', 'highEI', 'rftcLobe', 'rftcHemisphere'};
    pValEntry = {};
    channEntry = {};
    % rftcSite
    rfctSiteFlags = groupA.rftcSite > 0;
    channSel = rfctSiteFlags;
    pVal = signrank(groupA.(feature)(channSel), groupB.(feature)(channSel), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    channEntry = cat(2, channEntry, sum(channSel));
    % rftcStructure
    channSel = (groupA.rftcStruct == 1);
    channSel = channSel & not(rfctSiteFlags);
    pVal = signrank(groupA.(feature)(channSel), groupB.(feature)(channSel), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    channEntry = cat(2, channEntry, sum(channSel));
    % rftcConnected
    channSel = (groupA.rftcConn == 1);
    channSel = channSel & not(rfctSiteFlags);
    pVal = signrank(groupA.(feature)(channSel), groupB.(feature)(channSel), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    channEntry = cat(2, channEntry, sum(channSel));
    % highEI
    channSel = (groupA.highEI == 1);
    channSel = channSel & not(rfctSiteFlags);
    pVal = signrank(groupA.(feature)(channSel), groupB.(feature)(channSel), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    channEntry = cat(2, channEntry, sum(channSel));
    % rftcLobe
    channSel = (groupA.rftcLobe > 0);
    channSel = channSel & not(rfctSiteFlags) & not(groupA.rftcStruct > 0);
    pVal = signrank(groupA.(feature)(channSel), groupB.(feature)(channSel), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    channEntry = cat(2, channEntry, sum(channSel));
    % rftcHemis
    channSel = (groupA.rftcHemis == 1);
    channSel = channSel & not(rfctSiteFlags) & not(groupA.rftcStruct > 0) & not(groupA.rftcLobe > 0);
    pVal = signrank(groupA.(feature)(channSel), groupB.(feature)(channSel), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    channEntry = cat(2, channEntry, sum(channSel));
    
    analysisResP = cat(1, analysisResP, cat(2, {biomarker}, {feature}, pValEntry));
    analysisResNrChann = cat(1, analysisResNrChann, cat(2, {biomarker}, {feature}, channEntry));
end