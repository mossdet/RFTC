clc; clear all; close all;

paths = getFilesPaths();
biomarkersList =  {'HFO', 'iesHFO', 'IES'};
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

histEst = false;
nrHistPts = 25;

for bm = biomarkersList
    imprvd_vs_NonImprvd_Delta_P = {};
    
    nrZones = length(zonesNames);
    nrFeatures= length(featuresList);
    t = tiledlayout(nrFeatures, nrZones,"TileSpacing","tight");

    for fi = 1:length(featuresList)
        feature = featuresList{fi};
        niceFeature = niceFeaturesList{fi};
        biomarker = bm{1};

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

        groupA = imprvdDeltaPats;
        groupB = nonImprvsDeltaPats;

        % rftcSite
        channSelA = groupA.rftcSite == 1;
        channSelB = groupB.rftcSite == 1;
        zoneName = 'rftcSite';
        ti = nexttile;
        if histEst
            [f,xi] = ksdensity(groupA.(feature)(channSelA)); plot(xi,f); hold on;
            [f,xi] = ksdensity(groupB.(feature)(channSelB)); plot(xi,f); hold on;
        else
            histogram(groupA.(feature)(channSelA), nrHistPts, 'Normalization','pdf'); hold on;
            histogram(groupB.(feature)(channSelB), nrHistPts, 'Normalization','pdf');
        end
        legend('Pos.Outcome', 'Neg.Outcome');
        xlabel(niceFeature); ylabel('PDF');
        title(zoneName);
        pVal = ranksum(groupA.(feature)(channSelA), groupB.(feature)(channSelB), 'tail', 'right'); % 'both' | 'right' | 'left'
        strTxt = strcat("p = ", num2str(pVal));
        colorStr = 'green';
        if pVal >= 0.0005
            colorStr = 'red';
        end
        text(min(xlim), max(ylim),strTxt, 'Horiz','left', 'Vert','top', 'BackgroundColor', 'w', 'Color', colorStr, 'FontWeight', 'bold');
        
        % rftcStructure
        channSelA = (groupA.rftcStruct == 1) & not(groupA.rftcSite == 1);
        channSelB = (groupB.rftcStruct == 1) & not(groupB.rftcSite == 1);
        zoneName = 'rftcStructure';
        ti = nexttile;
        if histEst
            [f,xi] = ksdensity(groupA.(feature)(channSelA)); plot(xi,f); hold on;
            [f,xi] = ksdensity(groupB.(feature)(channSelB)); plot(xi,f); hold on;
        else
            histogram(groupA.(feature)(channSelA), nrHistPts, 'Normalization','pdf'); hold on;
            histogram(groupB.(feature)(channSelB), nrHistPts, 'Normalization','pdf');
        end
        legend('Pos.Outcome', 'Neg.Outcome');
        xlabel(niceFeature);
        title(zoneName);
        pVal = ranksum(groupA.(feature)(channSelA), groupB.(feature)(channSelB), 'tail', 'right'); % 'both' | 'right' | 'left'
        strTxt = strcat("p = ", num2str(pVal));
        colorStr = 'green';
        if pVal >= 0.0005
            colorStr = 'red';
        end
        text(min(xlim), max(ylim),strTxt, 'Horiz','left', 'Vert','top', 'BackgroundColor', 'w', 'Color', colorStr, 'FontWeight', 'bold');
        
        % rftcConnected
        channSelA = (groupA.rftcConn == 1) & not(groupA.rftcSite == 1);
        channSelB = (groupB.rftcConn == 1) & not(groupB.rftcSite == 1);
        zoneName = 'rftcConnected';
        ti = nexttile;
        if histEst
            [f,xi] = ksdensity(groupA.(feature)(channSelA)); plot(xi,f); hold on;
            [f,xi] = ksdensity(groupB.(feature)(channSelB)); plot(xi,f); hold on;
        else
            histogram(groupA.(feature)(channSelA), nrHistPts, 'Normalization','pdf'); hold on;
            histogram(groupB.(feature)(channSelB), nrHistPts, 'Normalization','pdf');
        end
        legend('Pos.Outcome', 'Neg.Outcome');
        xlabel(niceFeature);
        title(zoneName);
        pVal = ranksum(groupA.(feature)(channSelA), groupB.(feature)(channSelB), 'tail', 'right'); % 'both' | 'right' | 'left'
        strTxt = strcat("p = ", num2str(pVal));
        colorStr = 'green';
        if pVal >= 0.0005
            colorStr = 'red';
        end
        text(min(xlim), max(ylim),strTxt, 'Horiz','left', 'Vert','top', 'BackgroundColor', 'w', 'Color', colorStr, 'FontWeight', 'bold');

        % highEI
        channSelA = (groupA.highEI == 1) & not(groupA.rftcSite == 1);
        channSelB = (groupB.highEI == 1) & not(groupB.rftcSite == 1);
        zoneName = 'highEI';
        ti = nexttile;
        if histEst
            [f,xi] = ksdensity(groupA.(feature)(channSelA)); plot(xi,f); hold on;
            [f,xi] = ksdensity(groupB.(feature)(channSelB)); plot(xi,f); hold on;
        else
            histogram(groupA.(feature)(channSelA), nrHistPts, 'Normalization','pdf'); hold on;
            histogram(groupB.(feature)(channSelB), nrHistPts, 'Normalization','pdf');
        end
        legend('Pos.Outcome', 'Neg.Outcome');
        xlabel(niceFeature);
        title(zoneName);
        pVal = ranksum(groupA.(feature)(channSelA), groupB.(feature)(channSelB), 'tail', 'right'); % 'both' | 'right' | 'left'
        strTxt = strcat("p = ", num2str(pVal));
        colorStr = 'green';
        if pVal >= 0.0005
            colorStr = 'red';
        end
        text(min(xlim), max(ylim),strTxt, 'Horiz','left', 'Vert','top', 'BackgroundColor', 'w', 'Color', colorStr, 'FontWeight', 'bold');

        % rftcLobe
        channSelA = (groupA.rftcLobe == 1) & not(groupA.rftcSite == 1) & not(groupA.rftcStruct == 1);
        channSelB = (groupB.rftcLobe == 1) & not(groupB.rftcSite == 1) & not(groupB.rftcStruct == 1);
        zoneName = 'rftcLobe';
        ti = nexttile;
        if histEst
            [f,xi] = ksdensity(groupA.(feature)(channSelA)); plot(xi,f); hold on;
            [f,xi] = ksdensity(groupB.(feature)(channSelB)); plot(xi,f); hold on;
        else
            histogram(groupA.(feature)(channSelA), nrHistPts, 'Normalization','pdf'); hold on;
            histogram(groupB.(feature)(channSelB), nrHistPts, 'Normalization','pdf');
        end
        legend('Pos.Outcome', 'Neg.Outcome');
        xlabel(niceFeature);
        title(zoneName);
        pVal = ranksum(groupA.(feature)(channSelA), groupB.(feature)(channSelB), 'tail', 'right'); % 'both' | 'right' | 'left'
        strTxt = strcat("p = ", num2str(pVal));
        colorStr = 'green';
        if pVal >= 0.0005
            colorStr = 'red';
        end
        text(min(xlim), max(ylim),strTxt, 'Horiz','left', 'Vert','top', 'BackgroundColor', 'w', 'Color', colorStr, 'FontWeight', 'bold');

        % rftcHemisphere
        channSelA = (groupA.rftcHemis == 1) & not(groupA.rftcSite == 1) & not(groupA.rftcStruct == 1) & not(groupA.rftcLobe == 1);
        channSelB = (groupB.rftcHemis == 1) & not(groupB.rftcSite == 1) & not(groupB.rftcStruct == 1) & not(groupB.rftcLobe == 1);
        zoneName = 'rftcHemisphere';
        ti = nexttile;
        if histEst
            [f,xi] = ksdensity(groupA.(feature)(channSelA)); plot(xi,f); hold on;
            [f,xi] = ksdensity(groupB.(feature)(channSelB)); plot(xi,f); hold on;
        else
            histogram(groupA.(feature)(channSelA), nrHistPts, 'Normalization','pdf'); hold on;
            histogram(groupB.(feature)(channSelB), nrHistPts, 'Normalization','pdf');
        end
        legend('Pos.Outcome', 'Neg.Outcome');
        xlabel(niceFeature);
        title(zoneName);
        pVal = ranksum(groupA.(feature)(channSelA), groupB.(feature)(channSelB), 'tail', 'right'); % 'both' | 'right' | 'left'
        strTxt = strcat("p = ", num2str(pVal));
        colorStr = 'green';
        if pVal >= 0.0005
            colorStr = 'red';
        end
        text(min(xlim), max(ylim),strTxt, 'Horiz','left', 'Vert','top', 'BackgroundColor', 'w', 'Color', colorStr, 'FontWeight', 'bold');
    end
    set(gcf, 'Position', get(0, 'Screensize'), 'color','w');
    title(t, bm, 'FontSize',24, 'FontWeight','Bold')
    close();
end
