clc; clear all; close all;

paths = getFilesPaths();
patTablesPath = strcat(paths.workspacePath, 'CharacterizationTables\FeatureCharacterizationTables\');
groupTablesPath = strcat(paths.workspacePath, 'CharacterizationTables\GroupCharacterizationTablesAvg\');
analysisTablesPath = strcat(paths.workspacePath, 'AnalysisResults\Avg\Improved_vs_NonImproved\');mkdir(analysisTablesPath);

biomarkersList =  {'HFO', 'iesHFO', 'IES'};
featuresList = {'rate','duration','maxAmpl','sumAmpl','variance','BPLL','power','sumPower','mobility','complexity','PBRatio','wvltRipplePower','wvltFastRiplePower','peaks','zeroCrossPerEOI','spectCentroid','spectPeak','deltaOccRate','deltaAmpl','deltaPow','deltaVar','thetaOccRate','thetaAmpl','thetaPow','thetaVar','alphaOccRate','alphaAmpl','alphaPow','alphaVar','betaOccRate','betaAmpl','betaPow','betaVar','gammaOccRate','gammaAmpl','gammaPow','gammaVar','highGammaOccRate','highGammaAmpl','highGammaPow','highGammaVar','rippleOccRate','rippleAmpl','ripplePow','rippleVar','fRippleOccRate','fRippleAmpl','fRipplePow','fRippleVar'};
featuresList = {'rate','duration','maxAmpl','variance','power','PBRatio','peaks','zeroCrossPerEOI','spectCentroid','spectPeak'};
featuresList = {'maxAmpl','variance','power'};
normStrList = {''};% '', '_Normalized'
excList = {'Exc'}; % 'NonExc', 'Exc'

channSelStrList = {'FlexK'}; %{'Threshold', 'FlexK', 'Cluster2K'};

prctlTh = 75;
outcomeTh = 90;
corrLimitP = 0.001;


for normIdx = 1:length(normStrList)
    for exi = 1:length(excList)
        excStr = excList{exi};
        for chSelIdx = 1:length(channSelStrList)
            normStr = normStrList{normIdx};
            channSelStr = channSelStrList{chSelIdx};

            % Tables to read
            spreadSheetNamePre = strcat(groupTablesPath, 'GroupCharacterization_', channSelStr, '_Pre', normStr,'.xls');
            spreadSheetNamePost = strcat(groupTablesPath, 'GroupCharacterization_', channSelStr, '_Post', normStr,'.xls');

            % Tables to write

            tableImproved_vs_NonImproved_Pre_SpreadSheetName = strcat(analysisTablesPath, 'Improved_vs_nonImproved_Pre_Analysis_', channSelStr, '_', excStr, normStr,'.xls');
            tableImproved_vs_NonImproved_Post_SpreadSheetName = strcat(analysisTablesPath, 'Improved_vs_nonImproved_Post_Analysis_', channSelStr, '_', excStr, normStr,'.xls');
            tableImproved_vs_NonImproved_Delta_SpreadSheetName = strcat(analysisTablesPath, 'Improved_vs_nonImproved_Delta_Analysis_', channSelStr, '_', excStr, normStr,'.xls');
            delete(tableImproved_vs_NonImproved_Pre_SpreadSheetName);
            delete(tableImproved_vs_NonImproved_Delta_SpreadSheetName);
            delete(tableImproved_vs_NonImproved_Post_SpreadSheetName);

            for bm = biomarkersList
                imprvd_vs_NonImprvd_Pre_P = {};
                imprvd_vs_NonImprvd_Post_P = {};
                imprvd_vs_NonImprvd_Delta_P = {};

                biomarker = bm{1};
                groupTablePre = readtable(spreadSheetNamePre, 'Sheet', biomarker);
                groupTablePost = readtable(spreadSheetNamePost, 'Sheet', biomarker);
                featsNames = groupTablePre.Properties.VariableNames;
                %featuresList = featsNames(find(ismember(featsNames, 'rate')):end);

                if not(sum(strcmp(groupTablePre.chName, groupTablePost.chName)) == length(groupTablePost.chName))
                    error(' newGroupAnalysis_Paired, Pre and Post EEG Channels don''t match')
                end

                % Improved vs non-Improved
                % Pre
                selImprvdPatGroup = groupTablePre.outcome >= outcomeTh;
                imprvdGroupTablePre = groupTablePre(selImprvdPatGroup,:);
                nonImprvdGroupTablePre = groupTablePre(not(selImprvdPatGroup),:);
                % Post
                selImprvdPatGroup = groupTablePost.outcome >= outcomeTh;
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

                    % Tests
                    if strcmp(excStr, 'NonExc')
                    impVsNonImp_Pre_P_Res = zoneAnalysis_CompareIndependentGroups(imprvdGroupTablePre, 'imprvdGroupTablePre', nonImprvdGroupTablePre, 'nonImprvdGroupTablePre', biomarker, feature);
                    impVsNonImp_Post_P_Res = zoneAnalysis_CompareIndependentGroups(imprvdGroupTablePost, 'imprvdGroupTablePost', nonImprvdGroupTablePost, 'nonImprvdGroupTablePost', biomarker, feature);
                    impVsNonImp_Delta_P_Res = zoneAnalysis_CompareIndependentGroups(imprvdDeltaPats, 'imprvdDeltaPats', nonImprvsDeltaPats, 'nonImprvsDeltaPats', biomarker, feature);
                    elseif strcmp(excStr, 'Exc')
                        impVsNonImp_Pre_P_Res = zoneAnalysis_CompareIndependentGroupsExc(imprvdGroupTablePre, 'imprvdGroupTablePre', nonImprvdGroupTablePre, 'nonImprvdGroupTablePre', biomarker, feature, 'Pre');
                        impVsNonImp_Post_P_Res = zoneAnalysis_CompareIndependentGroupsExc(imprvdGroupTablePost, 'imprvdGroupTablePost', nonImprvdGroupTablePost, 'nonImprvdGroupTablePost', biomarker, feature, 'Post');
                        impVsNonImp_Delta_P_Res = zoneAnalysis_CompareIndependentGroupsExc(imprvdDeltaPats, 'imprvdDeltaPats', nonImprvsDeltaPats, 'nonImprvsDeltaPats', biomarker, feature, 'Delta');
                    end
                    
                    % Generate Table
                    if isempty(imprvd_vs_NonImprvd_Pre_P)
                        imprvd_vs_NonImprvd_Pre_P = cat(1, imprvd_vs_NonImprvd_Pre_P, impVsNonImp_Pre_P_Res);
                        imprvd_vs_NonImprvd_Post_P = cat(1, imprvd_vs_NonImprvd_Post_P, impVsNonImp_Post_P_Res);
                        imprvd_vs_NonImprvd_Delta_P = cat(1, imprvd_vs_NonImprvd_Delta_P, impVsNonImp_Delta_P_Res);
                    else
                        imprvd_vs_NonImprvd_Pre_P = cat(1, imprvd_vs_NonImprvd_Pre_P, impVsNonImp_Pre_P_Res(2:end,:));
                        imprvd_vs_NonImprvd_Post_P = cat(1, imprvd_vs_NonImprvd_Post_P, impVsNonImp_Post_P_Res(2:end,:));
                        imprvd_vs_NonImprvd_Delta_P = cat(1, imprvd_vs_NonImprvd_Delta_P, impVsNonImp_Delta_P_Res(2:end,:));
                    end
                end
                % Improved Pre vs non-Improved Pre
                tableImproved_vs_NonImproved_Pre = cell2table(imprvd_vs_NonImprvd_Pre_P(2:end,:),'VariableNames', imprvd_vs_NonImprvd_Pre_P(1,:));
                writetable(tableImproved_vs_NonImproved_Pre, tableImproved_vs_NonImproved_Pre_SpreadSheetName, 'Sheet', biomarker); 
                % Improved Post vs non-Improved Post
                tableImproved_vs_NonImproved_Post = cell2table(imprvd_vs_NonImprvd_Post_P(2:end,:),'VariableNames', imprvd_vs_NonImprvd_Post_P(1,:));
                writetable(tableImproved_vs_NonImproved_Post, tableImproved_vs_NonImproved_Post_SpreadSheetName, 'Sheet', biomarker);
                % Improved Delta vs non-Improved Delta
                tableImproved_vs_NonImproved_Delta = cell2table(imprvd_vs_NonImprvd_Delta_P(2:end,:),'VariableNames', imprvd_vs_NonImprvd_Delta_P(1,:));
                writetable(tableImproved_vs_NonImproved_Delta, tableImproved_vs_NonImproved_Delta_SpreadSheetName, 'Sheet', biomarker);  
            end
        end
    end
end

function analysisResP = zoneAnalysis_CompareIndependentGroups(groupA, labelA, groupB, labelB, biomarker, feature)

    analysisResP = {'Biomarker', 'Feature', 'rftcSite', 'rftcStructure', 'rftcConnected', 'highEI', 'rftcLobe', 'rftcHemisphere'};
    pValEntry = {};
    % rftcSite
    channSelA = groupA.rftcSite == 1;
    channSelB = groupB.rftcSite == 1;
    pVal = ranksum(groupA.(feature)(channSelA), groupB.(feature)(channSelB), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    % rftcStructure
    channSelA = (groupA.rftcStruct == 1);
    channSelB = (groupB.rftcStruct == 1);
    pVal = ranksum(groupA.(feature)(channSelA), groupB.(feature)(channSelB), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    % rftcConnected
    channSelA = (groupA.rftcConn == 1);
    channSelB = (groupB.rftcConn == 1);
    pVal = ranksum(groupA.(feature)(channSelA), groupB.(feature)(channSelB), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    % highEI
    channSelA = (groupA.highEI == 1);
    channSelB = (groupB.highEI == 1);
    pVal = ranksum(groupA.(feature)(channSelA), groupB.(feature)(channSelB), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    % rftcLobe
    channSelA = (groupA.rftcLobe == 1);
    channSelB = (groupB.rftcLobe == 1);
    pVal = ranksum(groupA.(feature)(channSelA), groupB.(feature)(channSelB), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    % rftcHemis
    channSelA = (groupA.rftcHemis == 1);
    channSelB = (groupB.rftcHemis == 1);
    pVal = ranksum(groupA.(feature)(channSelA), groupB.(feature)(channSelB), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    
    analysisResP = cat(1, analysisResP, cat(2, {biomarker}, {feature}, pValEntry));
end

function analysisResP = zoneAnalysis_CompareIndependentGroupsExc(groupA, labelA, groupB, labelB, biomarker, feature, descStr)

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
    
    close all
    subplot(1,2,1)
    nrBins = max([sum(channSelA), sum(channSelB)])
    histogram(groupB.(feature)(channSelB), nrBins, 'EdgeColor','none'); hold on;
    histogram(groupA.(feature)(channSelA), nrBins, 'EdgeColor','none'); hold on;
    title({strcat(biomarker, " from ", 'rftcStructure'); strcat(feature, " ", descStr)});
    xlabel(feature);
    ylabel('Nr. Channs');
    %plotPDF(groupA, channSelA, groupB, channSelB, biomarker, feature, 'rftcStructure', descStr);
    % rftcConnected
    channSelA = (groupA.rftcConn == 1);
    channSelA = channSelA & not(groupA.rftcSite == 1);
    channSelB = (groupB.rftcConn == 1);
    channSelB = channSelB & not(groupB.rftcSite == 1);    
    pVal = ranksum(groupA.(feature)(channSelA), groupB.(feature)(channSelB), 'tail', 'right');
    pValEntry = cat(2, pValEntry, pVal);
    
    subplot(1,2,2)
    nrBins = max([sum(channSelA), sum(channSelB)])
    histogram(groupB.(feature)(channSelB), nrBins, 'EdgeColor','none'); hold on;
    histogram(groupA.(feature)(channSelA), nrBins, 'EdgeColor','none'); hold on;
    title({strcat(biomarker, " from ", 'rftcConnected'); strcat(feature, " ", descStr)});
    xlabel(feature);
    ylabel('Nr. Channs');
    %plotPDF(groupA, channSelA, groupB, channSelB, biomarker, feature, 'rftcConnected', descStr);
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

function plotPDF(groupA, channSelA, groupB, channSelB, biomarker, feature, zone, descStr)
    close all
    subplot(1,2,1)
    nrBins = max([sum(channSelA), sum(channSelB)])
    histogram(groupB.(feature)(channSelB), nrBins, 'EdgeColor','none'); hold on;
    histogram(groupA.(feature)(channSelA), nrBins, 'EdgeColor','none'); hold on;
    title({strcat(biomarker, " from ", zone); strcat(feature, " ", descStr)});
    xlabel(feature);
    ylabel('Nr. Channs');
end