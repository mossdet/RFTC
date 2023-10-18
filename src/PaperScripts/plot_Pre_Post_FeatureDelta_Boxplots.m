clc; clear all; close all;

paths = getFilesPaths();
biomarkersList =  {'HFO', 'iesHFO', 'IES'};
featuresList = {'rate','duration','maxAmpl','sumAmpl','variance','BPLL','power','sumPower','mobility','complexity','PBRatio','wvltRipplePower','wvltFastRiplePower','peaks','zeroCrossPerEOI','spectCentroid','spectPeak','deltaOccRate','deltaAmpl','deltaPow','deltaVar','thetaOccRate','thetaAmpl','thetaPow','thetaVar','alphaOccRate','alphaAmpl','alphaPow','alphaVar','betaOccRate','betaAmpl','betaPow','betaVar','gammaOccRate','gammaAmpl','gammaPow','gammaVar','highGammaOccRate','highGammaAmpl','highGammaPow','highGammaVar','rippleOccRate','rippleAmpl','ripplePow','rippleVar','fRippleOccRate','fRippleAmpl','fRipplePow','fRippleVar'};
featuresList = {'rate','duration','maxAmpl','variance','BPLL','power'};
featuresList = {'rate', 'maxAmpl','variance','power'};

%featuresList = {'rate','duration','maxAmpl','variance','power'};

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
    for ft = featuresList
        feature = ft{1};

        groupTablePre = readtable(spreadSheetNamePre, 'Sheet', biomarker);
        groupTablePost = readtable(spreadSheetNamePost, 'Sheet', biomarker);
        featsNames = groupTablePre.Properties.VariableNames;

        if not(sum(strcmp(groupTablePre.chName, groupTablePost.chName)) == length(groupTablePost.chName))
            error(' newGroupAnalysis_Paired, Pre and Post EEG Channels don''t match')
        end

        [analysisResP, analysisResNrChann] = zoneAnalysis_CompareDependentGroups(groupTablePre, 'Pre', groupTablePost, 'Post', biomarker, feature);
        
        if isempty(resultsTableP)
            resultsTableP = cat(1, resultsTableP, analysisResP);
            resultsTableNrChanns = cat(1, resultsTableNrChanns, analysisResNrChann);
        else
            resultsTableP = cat(1, resultsTableP, analysisResP(2,:));
            resultsTableNrChanns = cat(1, resultsTableNrChanns, analysisResNrChann(2,:));
        end
    end
end

for bm = biomarkersList
    preData = [];
    postData= [];
    for ft = featuresList
        feature = ft{1};
        biomarker = bm{1};

        groupTablePre = readtable(spreadSheetNamePre, 'Sheet', biomarker);
        groupTablePost = readtable(spreadSheetNamePost, 'Sheet', biomarker);
        featsNames = groupTablePre.Properties.VariableNames;

        if not(sum(strcmp(groupTablePre.chName, groupTablePost.chName)) == length(groupTablePost.chName))
            error(' newGroupAnalysis_Paired, Pre and Post EEG Channels don''t match')
        end

        % rftcSite
        chSel = (groupTablePre.rftcSite>0);
        rftcSitePre = groupTablePre.(feature)(groupTablePre.rftcSite>0);
        rftcSitePost = groupTablePost.(feature)(groupTablePost.rftcSite>0);

        % rftcConnect
        chSel = (groupTablePre.rftcConn>0) & not(groupTablePre.rftcSite>0);
        rftcConnPre = groupTablePre.(feature)(chSel);
        rftcConnPost = groupTablePost.(feature)(chSel);

        % highEI
        chSel = (groupTablePre.highEI>0) & not(groupTablePre.rftcSite>0);
        highEI_Pre = groupTablePre.(feature)(chSel);
        highEI_Post = groupTablePost.(feature)(chSel);

        % rftcStruct
        chSel = (groupTablePre.rftcStruct>0) & not(groupTablePre.rftcSite>0);
        rftcStructPre = groupTablePre.(feature)(chSel);
        rftcStructPost = groupTablePost.(feature)(chSel);

        % rftcLobe
        chSel = (groupTablePre.rftcLobe>0) & not(groupTablePre.rftcStruct>0) & not(groupTablePre.rftcSite>0);
        rftcLobePre = groupTablePre.(feature)(chSel);
        rftcLobePost = groupTablePost.(feature)(chSel);

        % rftcHemis
        chSel = (groupTablePre.rftcHemis>0) & not(groupTablePre.rftcLobe>0) & not(groupTablePre.rftcStruct>0) & not(groupTablePre.rftcSite>0);
        rftcHemisPre = groupTablePre.(feature)(chSel);
        rftcHemisPost = groupTablePost.(feature)(chSel);

        preData.rftcSite.(feature) = rftcSitePre;
        preData.rftcStruct.(feature) = rftcStructPre;
        preData.rftcConn.(feature) = rftcConnPre;
        preData.highEI.(feature) = highEI_Pre;
        preData.rftcLobe.(feature) = rftcLobePre;
        preData.rftcHemis.(feature) = rftcHemisPre;

        postData.rftcSite.(feature) = rftcSitePost;
        postData.rftcStruct.(feature) = rftcStructPost;
        postData.rftcConn.(feature) = rftcConnPost;
        postData.highEI.(feature) = highEI_Post;
        postData.rftcLobe.(feature) = rftcLobePost;
        postData.rftcHemis.(feature) = rftcHemisPost;
    end
    % Plot Boxplots
    ydata = cat()
    boxchart(ydata)
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
