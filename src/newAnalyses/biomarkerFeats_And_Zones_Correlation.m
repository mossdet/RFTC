clc; clear all; close all;

paths = getFilesPaths();
biomarkersList =  {'HFO', 'iesHFO', 'IES'};

%featuresToAnalyze = {'rate','duration','maxAmpl','sumAmpl','variance','BPLL','power','sumPower','mobility','complexity','wvltRipplePower','wvltFastRiplePower','PBRatio','peaks','zeroCrossPerEOI','spectCentroid','spectPeak'};
%featuresToAnalyze = {'rate','duration','maxAmpl','variance','power','wvltRipplePower','wvltFastRiplePower', 'PBRatio', 'peaks','zeroCrossPerEOI','spectCentroid','spectPeak'};
%featuresToAnalyze = {'rate','duration','maxAmpl','variance','power','spectCentroid','spectPeak'};
%featuresToAnalyze = {'rate', 'maxAmpl', 'variance','power','wvltRipplePower','wvltFastRiplePower','deltaAmpl','deltaPow','deltaVar','thetaAmpl','thetaPow','thetaVar','alphaAmpl','alphaPow','alphaVar','betaAmpl','betaPow','betaVar','gammaAmpl','gammaPow','gammaVar','highGammaAmpl','highGammaPow','highGammaVar','rippleAmpl','ripplePow','rippleVar','fRippleAmpl','fRipplePow','fRippleVar'};
featuresList = {'rate','duration','maxAmpl','sumAmpl','variance','BPLL','power','sumPower','mobility','complexity','PBRatio','wvltRipplePower','wvltFastRiplePower','peaks','zeroCrossPerEOI','spectCentroid','spectPeak','deltaOccRate','deltaAmpl','deltaPow','deltaVar','thetaOccRate','thetaAmpl','thetaPow','thetaVar','alphaOccRate','alphaAmpl','alphaPow','alphaVar','betaOccRate','betaAmpl','betaPow','betaVar','gammaOccRate','gammaAmpl','gammaPow','gammaVar','highGammaOccRate','highGammaAmpl','highGammaPow','highGammaVar','rippleOccRate','rippleAmpl','ripplePow','rippleVar','fRippleOccRate','fRippleAmpl','fRipplePow','fRippleVar'};

%featuresToAnalyze = {'rate','duration','maxAmpl','variance','power'};

patTablesPath = strcat(paths.workspacePath, 'CharacterizationTables\FeatureCharacterizationTables\');
groupTablesPath = strcat(paths.workspacePath, 'CharacterizationTables\GroupCharacterizationTablesAvg\');
analysisTablesPath = strcat(paths.workspacePath, 'AnalysisResults\Avg\');mkdir(analysisTablesPath);

normStrList = {''}; % '', '_Normalized'
channSelStrList = {'FlexK'}; % 'Threshold', 'FlexK', 'Cluster2K'

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
        tableBiomCorrelAllPatients_NonExc_SpreadSheetName = strcat(analysisTablesPath, 'allPatsBiomarkerZoneCorrelation', channSelStr, '_NonExclusive', normStr,'.xls');
        tableBiomCorrelAllPatients_Exc_SpreadSheetName = strcat(analysisTablesPath, 'allPatsBiomarkerZoneCorrelation', channSelStr, '_Exclusive', normStr,'.xls');
        delete(tableBiomCorrelAllPatients_NonExc_SpreadSheetName);
        delete(tableBiomCorrelAllPatients_Exc_SpreadSheetName);
        
        tableZoneZoneCorrelAllPatients_NonExc_SpreadSheetName = strcat(analysisTablesPath, 'allPatsZoneZoneCorrelation', channSelStr, '_NonExclusive', normStr,'.xls');
        tableZoneZoneCorrelAllPatients_Exc_SpreadSheetName = strcat(analysisTablesPath, 'allPatsZoneZoneCorrelation', channSelStr, '_Exclusive', normStr,'.xls');
        delete(tableZoneZoneCorrelAllPatients_NonExc_SpreadSheetName);
        delete(tableZoneZoneCorrelAllPatients_Exc_SpreadSheetName);
        
        for bm = biomarkersList
            allGroupAnalysisP = {};
            imprvd_vs_NonImprvd_Pre_P = {};
            imprvd_vs_NonImprvd_Post_P = {};
            imprvd_vs_NonImprvd_Delta_P = {};
            correlationAnalysis = {};
            correlationAnalysisExclusive = {};
            zoneZoneCorrAnalysis = {};
            zoneZoneCorrAnalysisExclusive = {};

            biomarker = bm{1};

            groupTablePre = readtable(spreadSheetNamePre, 'Sheet', biomarker);
            groupTablePost = readtable(spreadSheetNamePost, 'Sheet', biomarker);
            featsNames = groupTablePre.Properties.VariableNames;
            featuresList = featsNames(find(ismember(featsNames, 'rate')):end);
            %featuresList = featuresToAnalyze;
            
            for ft = featuresList
                feature = ft{1};
                
                featCorrelAnal = biomarkerZonesCorrelationAnalysis_NonExclusive(groupTablePre, 'groupTablePre', biomarker, feature, corrLimitP);
                featCorrelAnalExclusive = biomarkerZonesCorrelationAnalysis_Exclusive(groupTablePre, 'groupTablePre', biomarker, feature, corrLimitP);

                zoneZoneCorrelAnal = zonesZonesCorrelationAnalysis_NonExclusive(groupTablePre, 'groupTablePre', biomarker, feature, corrLimitP);
                zoneZoneCorrelAnalExcl = zonesZonesCorrelationAnalysis_Exclusive(groupTablePre, 'groupTablePre', biomarker, feature, corrLimitP);
                                
                if isempty(correlationAnalysis)
                    correlationAnalysis = cat(1, correlationAnalysis, featCorrelAnal);
                    correlationAnalysisExclusive = cat(1, correlationAnalysisExclusive, featCorrelAnalExclusive);
                    zoneZoneCorrAnalysis = cat(1, zoneZoneCorrAnalysis, zoneZoneCorrelAnal);
                    zoneZoneCorrAnalysisExclusive = cat(1, zoneZoneCorrAnalysisExclusive, zoneZoneCorrelAnalExcl);
                else
                    correlationAnalysis = cat(1, correlationAnalysis, featCorrelAnal(2:end,:));
                    correlationAnalysisExclusive = cat(1, correlationAnalysisExclusive, featCorrelAnalExclusive(2:end,:));
                    zoneZoneCorrAnalysis = cat(1, zoneZoneCorrAnalysis, zoneZoneCorrelAnal(2:end,:));
                    zoneZoneCorrAnalysisExclusive = cat(1, zoneZoneCorrAnalysisExclusive, zoneZoneCorrelAnalExcl(2:end,:));
                end
            end
            
            % Biomarker Correlation
            tableBiomCorrAllPatsNonExc= cell2table(correlationAnalysis(2:end,:),'VariableNames', correlationAnalysis(1,:));
            writetable(tableBiomCorrAllPatsNonExc, tableBiomCorrelAllPatients_NonExc_SpreadSheetName, 'Sheet', biomarker);
            tableBiomCorrAllPats_Exc= cell2table(correlationAnalysisExclusive(2:end,:),'VariableNames', correlationAnalysisExclusive(1,:));
            writetable(tableBiomCorrAllPats_Exc, tableBiomCorrelAllPatients_Exc_SpreadSheetName, 'Sheet', biomarker);
            
            tableZoneZoneCorrAllPats_NonExc = cell2table(zoneZoneCorrAnalysis(2:end,:),'VariableNames', zoneZoneCorrAnalysis(1,:));
            writetable(tableZoneZoneCorrAllPats_NonExc, tableZoneZoneCorrelAllPatients_NonExc_SpreadSheetName, 'Sheet', biomarker);
            tableZoneZoneCorrAllPats_Exc = cell2table(zoneZoneCorrAnalysisExclusive(2:end,:),'VariableNames', zoneZoneCorrAnalysisExclusive(1,:));
            writetable(tableZoneZoneCorrAllPats_Exc, tableZoneZoneCorrelAllPatients_Exc_SpreadSheetName, 'Sheet', biomarker);

        end
    end
end

function zoneZoneCorrelAnal = zonesZonesCorrelationAnalysis_NonExclusive(group, label, biomarker, feature, corrLimitP)
    zonesList = {'rftcSite', 'rftcStruct', 'rftcConn', 'highEI', 'rftcLobe', 'rftcHemis'};
    zoneZoneCorrelAnal = { 'Zone', 'rftcSite', 'rftcStructure', 'rftcConnected', 'highEI', 'rftcLobe', 'rftcHemisphere'};
    corrLimitP = 0.01;
    
    for zia = 1:length(zonesList)
        zoneA =  zonesList{zia};
        corrZoneA = {zoneA};
        for zib = 1:length(zonesList)
            zoneB =  zonesList{zib};
            
            zoneA_Flags = group.(zoneA) > 0;
            zoneB_Flags = group.(zoneB) > 0;

            [rho,pval] = corr(zoneA_Flags, zoneB_Flags, 'Type','Spearman');
            mcc = getMCC(zoneA_Flags, zoneB_Flags);
            if pval > corrLimitP
                rho = 0;
            end
            corrZoneA = cat(2, corrZoneA, mcc); 
        end
        zoneZoneCorrelAnal = cat(1, zoneZoneCorrelAnal, corrZoneA); 
    end    
end

function zoneZoneCorrelAnal = zonesZonesCorrelationAnalysis_Exclusive(group, label, biomarker, feature, corrLimitP)
    zonesList = {'rftcSite', 'rftcStruct', 'rftcConn', 'highEI', 'rftcLobe', 'rftcHemis'};
    zoneZoneCorrelAnal = { 'Zone', 'rftcSite', 'rftcStructure', 'rftcConnected', 'highEI', 'rftcLobe', 'rftcHemisphere'};
    corrLimitP = 0.01;
    
    for zia = 1:length(zonesList)
        zoneA =  zonesList{zia};
        zoneA_Flags = group.(zoneA) > 0;
        corrZoneA = {zoneA};
        if not(strcmp(zoneA, 'rftcSite'))
            zoneA_Flags = zoneA_Flags & not(group.rftcSite > 0);
        end
        
        for zib = 1:length(zonesList)
            zoneB =  zonesList{zib};            
            zoneB_Flags = group.(zoneB) > 0;
            if not(strcmp(zoneB, 'rftcSite'))
                zoneB_Flags = zoneB_Flags & (group.rftcSite < 1);
            end
            [rho,pval] = corr(zoneA_Flags, zoneB_Flags, 'Type','Spearman');       
            mcc = getMCC(zoneA_Flags, zoneB_Flags);
            
            if pval > corrLimitP
                rho = 0;
            end
            corrZoneA = cat(2, corrZoneA, mcc); 
        end
        zoneZoneCorrelAnal = cat(1, zoneZoneCorrelAnal, corrZoneA); 
    end    
end

function analysisCorr = biomarkerZonesCorrelationAnalysis_NonExclusive(group, label, biomarker, feature, corrLimitP)

    analysisCorr = {'Biomarker', 'Feature', 'rftcSite', 'rftcStructure', 'rftcConnected', 'highEI', 'rftcLobe', 'rftcHemisphere'};
    corrValEntry = {};
    
    % rftcSite
    rfctSiteFlags = group.rftcSite;
    [rho,pval] = corr(rfctSiteFlags, group.(feature), 'Type','Spearman');
    if pval > corrLimitP
        rho = 0;
    end
    corrValEntry = cat(2, corrValEntry, rho);    
    
    % rftcStructure
    channSel = (group.rftcStruct == 1);
    [rho,pval] = corr(channSel, group.(feature), 'Type','Spearman');
    if pval > corrLimitP
        rho = 0;
    end
    corrValEntry = cat(2, corrValEntry, rho);
    
    % rftcConnected
    connValsVec = group.rftcElectroPhysioConnect;
    featValsVec = group.(feature);
    [rho,pval] = corr(connValsVec, featValsVec, 'Type','Spearman');
    if pval > corrLimitP
        rho = 0;
    end
    corrValEntry = cat(2, corrValEntry, rho);
    
    % highEI
    eiValsVec = group.eiVals;
    featValsVec = group.(feature);
    [rho,pval] = corr(eiValsVec, featValsVec, 'Type','Spearman');
    if pval > corrLimitP
        rho = 0;
    end
    corrValEntry = cat(2, corrValEntry, rho);
    
    % rftcLobe
    channSel = (group.rftcLobe == 1);
    [rho,pval] = corr(channSel, group.(feature), 'Type','Spearman');
    if pval > corrLimitP
        rho = 0;
    end
    corrValEntry = cat(2, corrValEntry, rho);
    
    % rftcHemis
    channSel = (group.rftcHemis == 1);
    [rho,pval] = corr(channSel, group.(feature), 'Type','Spearman');
    if pval > corrLimitP
        rho = 0;
    end
    corrValEntry = cat(2, corrValEntry, rho);
    
    analysisCorr = cat(1, analysisCorr, cat(2, {biomarker}, {feature}, corrValEntry));
end

function analysisCorr = biomarkerZonesCorrelationAnalysis_Exclusive(group, label, biomarker, feature, corrLimitP)

    analysisCorr = {'Biomarker', 'Feature', 'rftcSite', 'rftcStructure', 'rftcConnected', 'highEI', 'rftcLobe', 'rftcHemisphere'};
    corrValEntry = {};
    
    % rftcSite
    rfctSiteFlags = group.rftcSite > 0;
    [rho,pval] = corr(rfctSiteFlags, group.(feature), 'Type','Spearman');
    if pval > corrLimitP
        rho = 0;
    end
    corrValEntry = cat(2, corrValEntry, rho);    
    
    % rftcStructure
    channSel = (group.rftcStruct == 1) & not(rfctSiteFlags);
    [rho,pval] = corr(channSel, group.(feature), 'Type','Spearman');
    if pval > corrLimitP
        rho = 0;
    end
    corrValEntry = cat(2, corrValEntry, rho);
    
    % rftcConnected
    connValsVec = group.rftcElectroPhysioConnect;
    featValsVec = group.(feature);
    connValsVec(rfctSiteFlags) = [];
    featValsVec(rfctSiteFlags) = [];
    [rho,pval] = corr(connValsVec, featValsVec, 'Type','Spearman');
    if pval > corrLimitP
        rho = 0;
    end
    corrValEntry = cat(2, corrValEntry, rho);
    
    % highEI
    eiValsVec = group.eiVals;
    featValsVec = group.(feature);
    eiValsVec(rfctSiteFlags) = [];
    featValsVec(rfctSiteFlags) = [];
    [rho,pval] = corr(eiValsVec, featValsVec, 'Type','Spearman');
    if pval > corrLimitP
        rho = 0;
    end
    corrValEntry = cat(2, corrValEntry, rho);
    
    % rftcLobe
    channSel = (group.rftcLobe == 1);
    channSel = channSel & not(rfctSiteFlags) & not(group.rftcStruct);
    [rho,pval] = corr(channSel, group.(feature), 'Type','Spearman');
    if pval > corrLimitP
        rho = 0;
    end
    corrValEntry = cat(2, corrValEntry, rho);
    
    % rftcHemis
    channSel = (group.rftcHemis == 1);
    channSel = channSel & not(rfctSiteFlags) & not(group.rftcStruct) & not(group.rftcLobe);
    [rho,pval] = corr(channSel, group.(feature), 'Type','Spearman');
    if pval > corrLimitP
        rho = 0;
    end
    corrValEntry = cat(2, corrValEntry, rho);
    
    analysisCorr = cat(1, analysisCorr, cat(2, {biomarker}, {feature}, corrValEntry));
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
