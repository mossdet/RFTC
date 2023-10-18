clc; clear all; close all;

paths = getFilesPaths();

patTablesPath = strcat(paths.workspacePath, 'CharacterizationTables\FeatureCharacterizationTables\');
groupTablesPath = strcat(paths.workspacePath, 'CharacterizationTables\GroupCharacterizationTablesAvg\');
analysisTablesPath = strcat(paths.workspacePath, 'AnalysisResults\Avg\');mkdir(analysisTablesPath);

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

% Tables to write
tableBiomCorrelAllPatients_SpreadSheetName = strcat(analysisTablesPath, 'allPatsBiomarkerZoneCorrelation', channSelStr, '_', normStr, zoneFormation,'.xls');
delete(tableBiomCorrelAllPatients_SpreadSheetName);
        
for bm = biomarkersList
    correlationAnalysis = {};
    for ft = featuresList
        feature = ft{1};
        biomarker = bm{1};
        groupTablePre = readtable(spreadSheetNamePre, 'Sheet', biomarker);
        groupTablePost = readtable(spreadSheetNamePost, 'Sheet', biomarker);
        featsNames = groupTablePre.Properties.VariableNames;

        if not(sum(strcmp(groupTablePre.chName, groupTablePost.chName)) == length(groupTablePost.chName))
            error(' newGroupAnalysis_Paired, Pre and Post EEG Channels don''t match')
        end

        if isempty(zoneFormation)
            featCorrelAnal = biomarkerZonesCorrelationAnalysis_Exc(groupTablePre, 'groupTablePre', biomarker, feature, corrLimitP);
        else
            featCorrelAnal = biomarkerZonesCorrelationAnalysis_NonExc(groupTablePre, 'groupTablePre', biomarker, feature, corrLimitP);
        end
        if isempty(correlationAnalysis)
            correlationAnalysis = cat(1, correlationAnalysis, featCorrelAnal);
        else
            correlationAnalysis = cat(1, correlationAnalysis, featCorrelAnal(2:end,:));
        end
    end
    % Biomarker Correlation
    tableBiomCorrAllPats= cell2table(correlationAnalysis(2:end,:),'VariableNames', correlationAnalysis(1,:));
    writetable(tableBiomCorrAllPats, tableBiomCorrelAllPatients_SpreadSheetName, 'Sheet', biomarker); 
end

function analysisCorr = biomarkerZonesCorrelationAnalysis_NonExc(group, label, biomarker, feature, corrLimitP)

    analysisCorr = {'Biomarker', 'Feature', 'rftcSite', 'rftcStructure', 'rftcConnected', 'highEI', 'rftcLobe', 'rftcHemisphere'};
    corrValEntry = {};

    clusteredFeatVec = clusterVector_flexK_Group(group,feature);

    % rftcSite
    mcc = getMCC(group.rftcSite, clusteredFeatVec);
    corrValEntry = cat(2, corrValEntry, mcc);    
    
    % rftcStructure
    channSel = group.rftcStruct;
    mcc = getMCC(channSel, clusteredFeatVec);
    corrValEntry = cat(2, corrValEntry, mcc);
    
    % rftcConnected
    channSel = group.rftcConn;
    mcc = getMCC(channSel, clusteredFeatVec);
    corrValEntry = cat(2, corrValEntry, mcc);
    
    % highEI
    channSel = group.highEI;
    mcc = getMCC(channSel, clusteredFeatVec);
    corrValEntry = cat(2, corrValEntry, mcc);
    
    % rftcLobe
    channSel = group.rftcLobe;
    mcc = getMCC(channSel, clusteredFeatVec);
    corrValEntry = cat(2, corrValEntry, mcc);
    
    % rftcHemis
    channSel = group.rftcHemis;
    [rho,pval] = corr(channSel, group.(feature), 'Type','Spearman');
    mcc = getMCC(channSel, clusteredFeatVec);
    corrValEntry = cat(2, corrValEntry, mcc);
    
    analysisCorr = cat(1, analysisCorr, cat(2, {biomarker}, {feature}, corrValEntry));
end

function analysisCorr = biomarkerZonesCorrelationAnalysis_Exc(group, label, biomarker, feature, corrLimitP)

    analysisCorr = {'Biomarker', 'Feature', 'rftcSite', 'rftcStructure', 'rftcConnected', 'highEI', 'rftcLobe', 'rftcHemisphere'};
    corrValEntry = {};

    clusteredFeatVec = clusterVector_flexK_Group(group,feature);

    % rftcSite
    mcc = getMCC(group.rftcSite, clusteredFeatVec);
    corrValEntry = cat(2, corrValEntry, mcc);    
    
    % rftcStructure
    channSel = group.rftcStruct & not(group.rftcSite);
    mcc = getMCC(channSel, clusteredFeatVec);
    corrValEntry = cat(2, corrValEntry, mcc);
    
    % rftcConnected
    channSel = group.rftcConn & not(group.rftcSite);
    mcc = getMCC(channSel, clusteredFeatVec);
    corrValEntry = cat(2, corrValEntry, mcc);
    
    % highEI
    channSel = group.highEI & not(group.rftcSite);
    mcc = getMCC(channSel, clusteredFeatVec);
    corrValEntry = cat(2, corrValEntry, mcc);
    
    % rftcLobe
    channSel = group.rftcLobe & not(group.rftcSite) & not(group.rftcStruct);
    mcc = getMCC(channSel, clusteredFeatVec);
    corrValEntry = cat(2, corrValEntry, mcc);
    
    % rftcHemis
    channSel = group.rftcHemis & not(group.rftcSite) & not(group.rftcStruct) & not(group.rftcLobe);
    [rho,pval] = corr(channSel, group.(feature), 'Type','Spearman');
    mcc = getMCC(channSel, clusteredFeatVec);
    corrValEntry = cat(2, corrValEntry, mcc);
    
    analysisCorr = cat(1, analysisCorr, cat(2, {biomarker}, {feature}, corrValEntry));
end

function clusteredFeatVec = clusterVector_flexK_Group(group, feature)
    patNrs = unique(group.patNr);
    groupFeatVec = group.(feature);
    clusteredFeatVec = false(length(group.patNr),1);
    for pni = 1:length(patNrs)
       patNr = patNrs(pni);
       patSel = group.patNr == patNr;
       patVec = groupFeatVec(patSel);
       clstrIdx = clusterVector_flexK(patVec);
       clusteredFeatVec(patSel) = clstrIdx;
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