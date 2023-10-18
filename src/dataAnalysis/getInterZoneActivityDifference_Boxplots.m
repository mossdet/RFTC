clc; clear all; close all;
paths = getFilesPaths();
preFiles  = getPreFiles();%getAllFiles, getPreFiles, getPostFiles
postFiles  = getPostFiles();%getAllFiles, getPreFiles, getPostFiles
analysisPerZonesPath = 'F:\ForschungsProjekte\RFTC\MATLAB\DetectHFO\AnalysisPerZones_Boxplots\';

selDetector = 'MOSSDET'; 'Delphos'; 'MOSSDET'; 'MOSSDET_Depurated';
analysisPerZonesTablePath = strcat(analysisPerZonesPath, selDetector, '\');
features = {'OccRate', 'Power', 'Frequency'};
normOptions = {'Normalized'};
freqBandConnList = { 'maxAllBands'};
thList = 0;
patsSelList = {'allPatients', 'improvedPatients', 'nonImprovedPatients'};
biomarkersList =  {'iaVals', 'allHFOVals', 'iesHFOVals', 'isolHFOVals', 'allRippleVals', 'iesRippleVals', 'isolRippleVals', 'allFR_Vals', 'iesFR_Vals', 'isolFR_Vals'};

params.analysisPerZonesTablePath = analysisPerZonesTablePath;
params.selDetector = selDetector;
params.connTh = 75;
params.eiTh = 75;
params.outcomeTh = 49;

%% Analysis Loop
for psi = 1:length(patsSelList) 
    stopere = 1;
    for fi = 1:length(features)
        for bmi = 1:length(biomarkersList)

            %% Set params
            params.th = 0;
            params.normalization = 'Normalized';
            params.feature = features{fi};
            params.patsSel = patsSelList{psi};
            params.analysisSubType = features{fi};
            params.freqBandConn = 'maxAllBands';
            params.biomarker = biomarkersList{bmi};

            [groupTablePre, groupTablePost] = readTables(paths, params);

            %% Start analyses
            biomarkerActvDiff_RFTC = getRFTC_Channels_ActivityDiff(params, groupTablePre, groupTablePost);
            biomarkerActvDiff_highEI = getHighEI_Channels_ActivityDiff(params, groupTablePre, groupTablePost);
            biomarkerActvDiff_rftcConnected = getRFTC_ElectroPhysioConnected_Channels_ActivityDiff(params, groupTablePre, groupTablePost);

            sameRFTC_Parcel_ActvDiff = getSameRFTC_Parcel_ActivityDiff(params, groupTablePre, groupTablePost);
            sameRFTC_Lobe_ActvDiff = getSameRFTC_Lobe_ActivityDiff(params, groupTablePre, groupTablePost);
            sameRFTC_Hemisphere_ActvDiff = getSameRFTC_Hemisphere_ActivityDiff(params, groupTablePre, groupTablePost);

            plotData = {biomarkerActvDiff_RFTC; biomarkerActvDiff_highEI; biomarkerActvDiff_rftcConnected;...
                        sameRFTC_Parcel_ActvDiff; sameRFTC_Lobe_ActvDiff; sameRFTC_Hemisphere_ActvDiff};
            plotDataLabels = {'RFTC',...
                              'highEI',...
                              'rftcConnected',...
                              'sameParcel',...
                              'sameLobe',...
                              'sameHemisphere'};
            plotDataBoxChart(bmi, params, plotData, plotDataLabels);
        end
    end 
end

function plotDataBoxChart(plotIdx, params, plotData, plotDataLabels)
    nrGroups = length(plotData);
    spi = 1;
    colors = {[0 0.4470 0.7410];...	
                [0.8500 0.3250 0.0980];...
                [0.9290 0.6940 0.1250];...
                [0.4940 0.1840 0.5560];...
                [0.4660 0.6740 0.1880];...
                [0.3010 0.7450 0.9330];...
                [0.6350 0.0780 0.1840]};

    absMin = 0;
    absMax = 0;
    for gi = 1:nrGroups
        biomarkerAct = plotData{gi,1};
        allDat = [ biomarkerAct.in; biomarkerAct.out];
        if absMin > min(allDat)
            absMin = min(allDat);
        end
        if absMax < max(allDat)
            absMax = max(allDat);
        end
    end
    absMin = floor(absMin);
    absMax = ceil(absMax);

    for gi = 1:nrGroups
        subplot(1, nrGroups, spi)
        
        xgroupdata = [];
        ydata = [];
        biomarkerAct = plotData{gi,1};
        ydata = cat(1, ydata, biomarkerAct.in, biomarkerAct.out);
        xgroupdata = cat(1, xgroupdata, zeros(length(biomarkerAct.in),1)+1, zeros(length(biomarkerAct.out),1)+2);
        boxchart(xgroupdata, ydata, 'BoxFaceColor',colors{gi},'MarkerColor',colors{gi});
        spi = spi+1;
        
        [p, h] = ranksum(biomarkerAct.in, biomarkerAct.out);
        legendStr = strcat('p = ', num2str(p, 3));
        legend(legendStr, 'FontSize', 12);
        %xlabel('Grouping according to the RFTC Channels')
        if gi == 1
            ylabel({'Activity \Delta', '(preRFTC - postRFTC)'}, 'FontSize', 16)
        else
            set(gca,'yticklabel',[]);
        end
        ylim([absMin absMax]);
        xticks([1 2]);
        xStrA = strcat("(", num2str(biomarkerAct.nrPatsIn), ", ", num2str(biomarkerAct.nrChannsIn), ") ", plotDataLabels{gi});
        xStrB = strcat("(", num2str(biomarkerAct.nrPatsOut),", ", num2str(biomarkerAct.nrChannsOut), ") ", 'Rest');
        set(gca,'XTickLabel',{xStrA, xStrB});
        ax = gca;
        ax.XAxis.FontSize = 10;
        xtickangle(20);
        set(findobj(gcf,'type','axes'),'FontWeight','Bold');
    end
    biomrkrStr = strrep(params.biomarker, 'Vals', '');
    biomrkrStr = strrep(biomrkrStr, 'all', '');
    sgtitle({biomrkrStr; params.feature},'FontSize',20);
    set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w');
        normStr = 'nonNormalized';
    if strcmp(params.normalization, 'Normalized')
        normStr = params.normalization;
    end
    boxchartPath = strcat(params.analysisPerZonesTablePath, 'Inter_Zone_Analysis\', normStr, '\', params.patsSel, '\', params.feature, '\'); mkdir(boxchartPath)
    boxchartFN = strcat(boxchartPath, num2str(plotIdx), '_', strcat(params.biomarker, params.feature));
    saveas(gcf,boxchartFN);
    boxchartFN = strcat(boxchartFN, '.jpg');
    saveas(gcf,boxchartFN);
    close();
end

function biomarkerActvDiff = getRFTC_Channels_ActivityDiff(params, groupTablePre, groupTablePost)
    biomarkerActvDiff.in = [];
    biomarkerActvDiff.out = [];
    groupTablePreIn = groupTablePre;
    groupTablePostIn = groupTablePost;
    groupTablePreOut = groupTablePre;
    groupTablePostOut = groupTablePost;
    
    %%
    tableDelIdxs = groupTablePreIn.rftcVals < 1;
    groupTablePreIn(tableDelIdxs, :) = [];
    groupTablePostIn(tableDelIdxs, :) = [];
    
    %%
    tableDelIdxs = groupTablePreOut.rftcVals > 0;
    groupTablePreOut(tableDelIdxs, :) = [];
    groupTablePostOut(tableDelIdxs, :) = [];
    
    %%
    biomarkerActvDiff.in = groupTablePreIn.(params.biomarker) - groupTablePostIn.(params.biomarker);
    biomarkerActvDiff.out = groupTablePreOut.(params.biomarker) - groupTablePostOut.(params.biomarker);
    biomarkerActvDiff.nrChannsIn = length(groupTablePreIn.patName);
    biomarkerActvDiff.nrChannsOut = length(groupTablePreOut.patName);
    biomarkerActvDiff.nrPatsIn = length(unique(groupTablePreIn.patName));
    biomarkerActvDiff.nrPatsOut = length(unique(groupTablePreOut.patName));
end

function biomarkerActvDiff = getHighEI_Channels_ActivityDiff(params, groupTablePre, groupTablePost)
    biomarkerActvDiff.in = [];
    biomarkerActvDiff.out = [];
    groupTablePreIn = groupTablePre;
    groupTablePostIn = groupTablePost;
    groupTablePreOut = groupTablePre;
    groupTablePostOut = groupTablePost;
    
    %%
    tableDelIdxs = not(perctThValsPerPat(groupTablePre.patName, groupTablePre.eiVals, params.eiTh));
    groupTablePreIn(tableDelIdxs, :) = [];
    groupTablePostIn(tableDelIdxs, :) = [];
    
    %%
    tableDelIdxs = not(tableDelIdxs);
    groupTablePreOut(tableDelIdxs, :) = [];
    groupTablePostOut(tableDelIdxs, :) = [];
        
    %% 
    groupTablePreIn = removeRFTC_Channels(groupTablePreIn);
    groupTablePostIn = removeRFTC_Channels(groupTablePostIn);
    groupTablePreOut = removeRFTC_Channels(groupTablePreOut);
    groupTablePostOut = removeRFTC_Channels(groupTablePostOut);
    
    %%
    biomarkerActvDiff.in = groupTablePreIn.(params.biomarker) - groupTablePostIn.(params.biomarker);
    biomarkerActvDiff.out = groupTablePreOut.(params.biomarker) - groupTablePostOut.(params.biomarker);
    biomarkerActvDiff.nrChannsIn = length(groupTablePreIn.patName);
    biomarkerActvDiff.nrChannsOut = length(groupTablePreOut.patName);
    biomarkerActvDiff.nrPatsIn = length(unique(groupTablePreIn.patName));
    biomarkerActvDiff.nrPatsOut = length(unique(groupTablePreOut.patName));
end

function biomarkerActvDiff = getRFTC_ElectroPhysioConnected_Channels_ActivityDiff(params, groupTablePre, groupTablePost)
    biomarkerActvDiff.in = [];
    biomarkerActvDiff.out = [];
    groupTablePreIn = groupTablePre;
    groupTablePostIn = groupTablePost;
    groupTablePreOut = groupTablePre;
    groupTablePostOut = groupTablePost;
    
    %%
    tableDelIdxs = not(perctThValsPerPat(groupTablePre.patName, groupTablePre.rftcElectroPhysioConnect, params.connTh));
    groupTablePreIn(tableDelIdxs, :) = [];
    groupTablePostIn(tableDelIdxs, :) = [];
    
    %%
    tableDelIdxs = not(tableDelIdxs);
    groupTablePreOut(tableDelIdxs, :) = [];
    groupTablePostOut(tableDelIdxs, :) = [];
    
    %% 
    groupTablePreIn = removeRFTC_Channels(groupTablePreIn);
    groupTablePostIn = removeRFTC_Channels(groupTablePostIn);
    groupTablePreOut = removeRFTC_Channels(groupTablePreOut);
    groupTablePostOut = removeRFTC_Channels(groupTablePostOut);
    
    %%
    biomarkerActvDiff.in = groupTablePreIn.(params.biomarker) - groupTablePostIn.(params.biomarker);
    biomarkerActvDiff.out = groupTablePreOut.(params.biomarker) - groupTablePostOut.(params.biomarker);
    biomarkerActvDiff.nrChannsIn = length(groupTablePreIn.patName);
    biomarkerActvDiff.nrChannsOut = length(groupTablePreOut.patName);
    biomarkerActvDiff.nrPatsIn = length(unique(groupTablePreIn.patName));
    biomarkerActvDiff.nrPatsOut = length(unique(groupTablePreOut.patName));
end

function biomarkerActvDiff = getSameRFTC_Parcel_ActivityDiff(params, groupTablePre, groupTablePost)
    biomarkerActvDiff.in = [];
    biomarkerActvDiff.out = []; 
    groupTablePreIn = groupTablePre;
    groupTablePostIn = groupTablePost;
    groupTablePreOut = groupTablePre;
    groupTablePostOut = groupTablePost;
    
    %%
    tableDelIdxs = not(sameRFTC_BrainParcelChannsPerPat(groupTablePre));
    groupTablePreIn(tableDelIdxs, :) = [];
    groupTablePostIn(tableDelIdxs, :) = [];
    
    %%
    tableDelIdxs = not(tableDelIdxs);
    groupTablePreOut(tableDelIdxs, :) = [];
    groupTablePostOut(tableDelIdxs, :) = [];
    
    %% 
    groupTablePreIn = removeRFTC_Channels(groupTablePreIn);
    groupTablePostIn = removeRFTC_Channels(groupTablePostIn);
    groupTablePreOut = removeRFTC_Channels(groupTablePreOut);
    groupTablePostOut = removeRFTC_Channels(groupTablePostOut);
    
    %%
    biomarkerActvDiff.in = groupTablePreIn.(params.biomarker) - groupTablePostIn.(params.biomarker);
    biomarkerActvDiff.out = groupTablePreOut.(params.biomarker) - groupTablePostOut.(params.biomarker);
    biomarkerActvDiff.nrChannsIn = length(groupTablePreIn.patName);
    biomarkerActvDiff.nrChannsOut = length(groupTablePreOut.patName);
    biomarkerActvDiff.nrPatsIn = length(unique(groupTablePreIn.patName));
    biomarkerActvDiff.nrPatsOut = length(unique(groupTablePreOut.patName));
end

function biomarkerActvDiff = getSameRFTC_Lobe_ActivityDiff(params, groupTablePre, groupTablePost)
    biomarkerActvDiff.in = [];
    biomarkerActvDiff.out = [];
    groupTablePreIn = groupTablePre;
    groupTablePostIn = groupTablePost;
    groupTablePreOut = groupTablePre;
    groupTablePostOut = groupTablePost;
    
    %%
    tableDelIdxs = not(sameRFTC_BrainLobeChannsPerPat(groupTablePre));
    groupTablePreIn(tableDelIdxs, :) = [];
    groupTablePostIn(tableDelIdxs, :) = [];
    
    %%
    tableDelIdxs = not(tableDelIdxs);
    groupTablePreOut(tableDelIdxs, :) = [];
    groupTablePostOut(tableDelIdxs, :) = [];
    
    %% 
    groupTablePreIn = removeRFTC_Channels(groupTablePreIn);
    groupTablePostIn = removeRFTC_Channels(groupTablePostIn);
    groupTablePreOut = removeRFTC_Channels(groupTablePreOut);
    groupTablePostOut = removeRFTC_Channels(groupTablePostOut);
    
    %%
    biomarkerActvDiff.in = groupTablePreIn.(params.biomarker) - groupTablePostIn.(params.biomarker);
    biomarkerActvDiff.out = groupTablePreOut.(params.biomarker) - groupTablePostOut.(params.biomarker);
    biomarkerActvDiff.nrChannsIn = length(groupTablePreIn.patName);
    biomarkerActvDiff.nrChannsOut = length(groupTablePreOut.patName);
    biomarkerActvDiff.nrPatsIn = length(unique(groupTablePreIn.patName));
    biomarkerActvDiff.nrPatsOut = length(unique(groupTablePreOut.patName));
end

function biomarkerActvDiff = getSameRFTC_Hemisphere_ActivityDiff(params, groupTablePre, groupTablePost)
    biomarkerActvDiff.in = [];
    biomarkerActvDiff.out = [];
    groupTablePreIn = groupTablePre;
    groupTablePostIn = groupTablePost;
    groupTablePreOut = groupTablePre;
    groupTablePostOut = groupTablePost;
    
    %%
    tableDelIdxs = not(sameRFTC_HemisphereChannsPerPat(groupTablePre));
    groupTablePreIn(tableDelIdxs, :) = [];
    groupTablePostIn(tableDelIdxs, :) = [];
    
    %%
    tableDelIdxs = not(tableDelIdxs);
    groupTablePreOut(tableDelIdxs, :) = [];
    groupTablePostOut(tableDelIdxs, :) = [];
    
    %% 
    groupTablePreIn = removeRFTC_Channels(groupTablePreIn);
    groupTablePostIn = removeRFTC_Channels(groupTablePostIn);
    groupTablePreOut = removeRFTC_Channels(groupTablePreOut);
    groupTablePostOut = removeRFTC_Channels(groupTablePostOut);
    
    %%
    biomarkerActvDiff.in = groupTablePreIn.(params.biomarker) - groupTablePostIn.(params.biomarker);
    biomarkerActvDiff.out = groupTablePreOut.(params.biomarker) - groupTablePostOut.(params.biomarker);
    biomarkerActvDiff.nrChannsIn = length(groupTablePreIn.patName);
    biomarkerActvDiff.nrChannsOut = length(groupTablePreOut.patName);
    biomarkerActvDiff.nrPatsIn = length(unique(groupTablePreIn.patName));
    biomarkerActvDiff.nrPatsOut = length(unique(groupTablePreOut.patName));
end

function groupTable = removeRFTC_Channels(groupTable)
    tableDelChannsRFTC = groupTable.rftcVals > 0;
    groupTable(tableDelChannsRFTC, :) = [];
end

function allPatsKeepIdxs = perctThValsPerPat(patNameCol, vals, th)
    patNames = unique(patNameCol);
    allPatsKeepIdxs = zeros(length(vals),1);
    for pi = 1:length(patNames)
        patName = patNames{pi};
        patSelIdx = ismember(patNameCol, patName);
        patVals = vals(patSelIdx);
        patKeepIdx = patVals >= prctile(patVals, th);
        allPatsKeepIdxs(patSelIdx) = patKeepIdx;        
    end
    allPatsKeepIdxs = logical(allPatsKeepIdxs);
end

function [groupTablePre, groupTablePost] = readTables(paths, params)
    %% Read tables
    groupAnalysisTablesFilePath = strcat(paths.workspacePath, 'GroupAnalysis_ChannelCharacterizationTables\', params.selDetector, '\');
    preGroupTableFN = strcat(groupAnalysisTablesFilePath, params.normalization, 'GroupAnalysis_ChannelCharacterization_Pre', params.selDetector, '.xls');
    postGroupTableFN = strcat(groupAnalysisTablesFilePath, params.normalization, 'GroupAnalysis_ChannelCharacterization_Post', params.selDetector, '.xls');
    groupTablePre = readtable(preGroupTableFN, 'Sheet', params.feature);
    groupTablePost = readtable(postGroupTableFN, 'Sheet', params.feature);
    
    %% read normalized EI
    preGroupTableFN = strcat(groupAnalysisTablesFilePath, 'Normalized','GroupAnalysis_ChannelCharacterization_Pre', params.selDetector, '.xls');
    postGroupTableFN = strcat(groupAnalysisTablesFilePath, 'Normalized','GroupAnalysis_ChannelCharacterization_Post', params.selDetector, '.xls');
    groupTablePreEI_correct = readtable(preGroupTableFN, 'Sheet', params.feature);
    groupTablePostEI_correct = readtable(postGroupTableFN, 'Sheet', params.feature);
    groupTablePre.eiVals = groupTablePreEI_correct.eiVals;
    groupTablePost.eiVals = groupTablePostEI_correct.eiVals;
    
    %% read normalized rftcElectroPhysioConnect
    preGroupTableNormFN = strcat(groupAnalysisTablesFilePath, 'Normalized', 'GroupAnalysis_ChannelCharacterization_Pre', params.selDetector,  '_', params.freqBandConn, '.xls');
    postGroupTableNormFN = strcat(groupAnalysisTablesFilePath, 'Normalized', 'GroupAnalysis_ChannelCharacterization_Post', params.selDetector, '_', params.freqBandConn,'.xls');
    groupTablePre_RFTCconn = readtable(preGroupTableNormFN, 'Sheet', params.feature);
    groupTablePost_RFTCconn = readtable(postGroupTableNormFN, 'Sheet', params.feature);
    groupTablePre.rftcElectroPhysioConnect = groupTablePre_RFTCconn.rftcElectroPhysioConnect;
    groupTablePost.rftcElectroPhysioConnect = groupTablePost_RFTCconn.rftcElectroPhysioConnect;
    
    %% Select outcomes
    if strcmp(params.patsSel, 'improvedPatients')        
        tableDelIdxsIn = groupTablePre.outcome <= params.outcomeTh;
        groupTablePre(tableDelIdxsIn, :) = [];
        groupTablePost(tableDelIdxsIn, :) = [];        
    elseif strcmp(params.patsSel, 'nonImprovedPatients')
        tableDelIdxsIn = groupTablePre.outcome > params.outcomeTh;
        groupTablePre(tableDelIdxsIn, :) = [];
        groupTablePost(tableDelIdxsIn, :) = []; 
    end
end