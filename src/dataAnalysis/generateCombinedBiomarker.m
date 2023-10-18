clc; clear all; close all;
paths = getFilesPaths();
preFiles  = getPreFiles();%getAllFiles, getPreFiles, getPostFiles
postFiles  = getPostFiles();%getAllFiles, getPreFiles, getPostFiles
analysisPerZonesPath = 'F:\ForschungsProjekte\RFTC\MATLAB\DetectHFO\AnalysisPerZones_Exclusive_Boxplots\';

selDetector = 'MOSSDET'; 'Delphos'; 'MOSSDET'; 'MOSSDET_Depurated';
analysisPerZonesTablePath = strcat(analysisPerZonesPath, selDetector, '\');
features = {'OccRate', 'Power', 'Frequency'};
%features = {'OccRate', 'Power'};
normOptions = {'', 'Normalized'}; %'Normalized'; % {'', 'Scaled', 'Normalized'};
freqBandConnList = { 'maxAllBands'};
thList = 0;
patsSelList = {'allPatients', 'improvedPatients', 'nonImprovedPatients'};%{'allPatients', 'improvedPatients', 'nonImprovedPatients'};
patsSelList = {'allPatients', 'improvedPatients', 'nonImprovedPatients'};
biomarkersList =  {'iesHFOVals', 'iaVals', 'allHFOVals', 'isolHFOVals', 'allRippleVals', 'iesRippleVals', 'isolRippleVals', 'allFR_Vals', 'iesFR_Vals', 'isolFR_Vals'};
biomarkersList =  {'iesHFOVals', 'allHFOVals', 'iaVals'};

params.analysisPerZonesTablePath = analysisPerZonesTablePath;
params.selDetector = selDetector;
params.connTh = 75;
params.eiTh = 75;
params.outcomeTh = 90;

load(strcat(paths.workspacePath,'AllZonesData\', 'AllZonesData_In.mat'))
load(strcat(paths.workspacePath,'AllZonesData\', 'AllZonesData_Out.mat'))

%% Analysis Loop
for bmi = 1:length(biomarkersList)
    for nmi = 1:length(normOptions)
        for psi = 1:length(patsSelList)
            for fi = 1:length(features)
                %% Set params
                params.th = 0;
                params.normalization = normOptions{nmi};
                params.feature = features{fi};
                params.patsSel = patsSelList{psi};
                params.analysisSubType = features{fi};
                params.freqBandConn = 'maxAllBands';
                params.biomarker = biomarkersList{bmi};
                if isempty(params.normalization)
                    params.normalization = 'noNorm';
                end

                %% rftcSite
                inRatePre = zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('OccRate').rftcPre;
                inPowPre = zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('Power').rftcPre;
                inRatePowPre = sqrt(inRatePre.^2 + inPowPre.^2);
                inRatePost = zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('OccRate').rftcPost;
                inPowPost = zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('Power').rftcPost;
                inRatePowPost = sqrt(inRatePost.^2 + inPowPost.^2);
                inRatePowDiff_rftcSite = inRatePowPre - inRatePowPost;
                zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').rftcPre = inRatePowPre;
                zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').rftcPost = inRatePowPost;
                zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').rftcDiff = inRatePowDiff_rftcSite;

                %% rftcConnected
                inRatePre = zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('OccRate').rftcConnPre;
                inPowPre = zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('Power').rftcConnPre;
                inRatePowPre = sqrt(inRatePre.^2 + inPowPre.^2);
                inRatePost = zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('OccRate').rftcConnPost;
                inPowPost = zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('Power').rftcConnPost;
                inRatePowPost = sqrt(inRatePost.^2 + inPowPost.^2);
                inRatePowDiff_rftcConn = inRatePowPre - inRatePowPost;
                zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').rftcConnPre = inRatePowPre;
                zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').rftcConnPost = inRatePowPost;
                zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').rftcConnDiff = inRatePowDiff_rftcConn;
 
                %% rftcStructure
                inRatePre = zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('OccRate').rftcStructPre;
                inPowPre = zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('Power').rftcStructPre;
                inRatePowPre = sqrt(inRatePre.^2 + inPowPre.^2);
                inRatePost = zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('OccRate').rftcStructPost;
                inPowPost = zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('Power').rftcStructPost;
                inRatePowPost = sqrt(inRatePost.^2 + inPowPost.^2);
                inRatePowDiff_rftcStruct = inRatePowPre - inRatePowPost;
                zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').rftcStructPre = inRatePowPre;
                zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').rftcStructPost = inRatePowPost;
                zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').rftcStructDiff = inRatePowDiff_rftcStruct;

                
                %% highEI
                inRatePre = zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('OccRate').highEI_Pre;
                inPowPre = zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('Power').highEI_Pre;
                inRatePowPre = sqrt(inRatePre.^2 + inPowPre.^2);
                inRatePost = zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('OccRate').highEI_Post;
                inPowPost = zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('Power').highEI_Post;
                inRatePowPost = sqrt(inRatePost.^2 + inPowPost.^2);
                inRatePowDiff_highEI = inRatePowPre - inRatePowPost;
                zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').highEI_Pre = inRatePowPre;
                zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').highEI_Post = inRatePowPost;
                zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').highEI_Diff = inRatePowDiff_highEI;


                %% highEI
                inRatePre = zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('OccRate').rftcLobePre;
                inPowPre = zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('Power').rftcLobePre;
                inRatePowPre = sqrt(inRatePre.^2 + inPowPre.^2);
                inRatePost = zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('OccRate').rftcLobePost;
                inPowPost = zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('Power').rftcLobePost;
                inRatePowPost = sqrt(inRatePost.^2 + inPowPost.^2);
                inRatePowDiff_rftcLobe = inRatePowPre - inRatePowPost;
                zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').rftcLobePre = inRatePowPre;
                zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').rftcLobePost = inRatePowPost;
                zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').rftcLobeDiff = inRatePowDiff_rftcLobe;


                %% highEI
                inRatePre = zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('OccRate').rftcHemisPre;
                inPowPre = zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('Power').rftcHemisPre;
                inRatePowPre = sqrt(inRatePre.^2 + inPowPre.^2);
                inRatePost = zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('OccRate').rftcHemisPost;
                inPowPost = zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('Power').rftcHemisPost;
                inRatePowPost = sqrt(inRatePost.^2 + inPowPost.^2);
                inRatePowDiff_rftcHemisphere = inRatePowPre - inRatePowPost;
                zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').rftcHemisPre = inRatePowPre;
                zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').rftcHemisPost = inRatePowPost;
                zonesDataIn.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').rftcHemisDiff = inRatePowDiff_rftcHemisphere;
            end
        end 
    end
end

%% Analysis Loop
for bmi = 1:length(biomarkersList)
    for nmi = 1:length(normOptions)
        for psi = 1:length(patsSelList)
            for fi = 1:length(features)
                %% Set params
                params.th = 0;
                params.normalization = normOptions{nmi};
                params.feature = features{fi};
                params.patsSel = patsSelList{psi};
                params.analysisSubType = features{fi};
                params.freqBandConn = 'maxAllBands';
                params.biomarker = biomarkersList{bmi};
                if isempty(params.normalization)
                    params.normalization = 'noNorm';
                end

                %% rftcSite
                inRatePre = zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('OccRate').rftcPre;
                inPowPre = zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('Power').rftcPre;
                inRatePowPre = sqrt(inRatePre.^2 + inPowPre.^2);
                inRatePost = zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('OccRate').rftcPost;
                inPowPost = zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('Power').rftcPost;
                inRatePowPost = sqrt(inRatePost.^2 + inPowPost.^2);
                inRatePowDiff_rftcSite = inRatePowPre - inRatePowPost;
                zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').rftcPre = inRatePowPre;
                zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').rftcPost = inRatePowPost;
                zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').rftcDiff = inRatePowDiff_rftcSite;

                %% rftcConnected
                inRatePre = zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('OccRate').rftcConnPre;
                inPowPre = zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('Power').rftcConnPre;
                inRatePowPre = sqrt(inRatePre.^2 + inPowPre.^2);
                inRatePost = zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('OccRate').rftcConnPost;
                inPowPost = zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('Power').rftcConnPost;
                inRatePowPost = sqrt(inRatePost.^2 + inPowPost.^2);
                inRatePowDiff_rftcConn = inRatePowPre - inRatePowPost;
                zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').rftcConnPre = inRatePowPre;
                zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').rftcConnPost = inRatePowPost;
                zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').rftcConnDiff = inRatePowDiff_rftcConn;
 
                %% rftcStructure
                inRatePre = zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('OccRate').rftcStructPre;
                inPowPre = zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('Power').rftcStructPre;
                inRatePowPre = sqrt(inRatePre.^2 + inPowPre.^2);
                inRatePost = zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('OccRate').rftcStructPost;
                inPowPost = zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('Power').rftcStructPost;
                inRatePowPost = sqrt(inRatePost.^2 + inPowPost.^2);
                inRatePowDiff_rftcStruct = inRatePowPre - inRatePowPost;
                zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').rftcStructPre = inRatePowPre;
                zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').rftcStructPost = inRatePowPost;
                zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').rftcStructDiff = inRatePowDiff_rftcStruct;

                
                %% highEI
                inRatePre = zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('OccRate').highEI_Pre;
                inPowPre = zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('Power').highEI_Pre;
                inRatePowPre = sqrt(inRatePre.^2 + inPowPre.^2);
                inRatePost = zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('OccRate').highEI_Post;
                inPowPost = zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('Power').highEI_Post;
                inRatePowPost = sqrt(inRatePost.^2 + inPowPost.^2);
                inRatePowDiff_highEI = inRatePowPre - inRatePowPost;
                zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').highEI_Pre = inRatePowPre;
                zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').highEI_Post = inRatePowPost;
                zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').highEI_Diff = inRatePowDiff_highEI;


                %% highEI
                inRatePre = zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('OccRate').rftcLobePre;
                inPowPre = zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('Power').rftcLobePre;
                inRatePowPre = sqrt(inRatePre.^2 + inPowPre.^2);
                inRatePost = zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('OccRate').rftcLobePost;
                inPowPost = zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('Power').rftcLobePost;
                inRatePowPost = sqrt(inRatePost.^2 + inPowPost.^2);
                inRatePowDiff_rftcLobe = inRatePowPre - inRatePowPost;
                zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').rftcLobePre = inRatePowPre;
                zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').rftcLobePost = inRatePowPost;
                zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').rftcLobeDiff = inRatePowDiff_rftcLobe;


                %% highEI
                inRatePre = zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('OccRate').rftcHemisPre;
                inPowPre = zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('Power').rftcHemisPre;
                inRatePowPre = sqrt(inRatePre.^2 + inPowPre.^2);
                inRatePost = zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('OccRate').rftcHemisPost;
                inPowPost = zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('Power').rftcHemisPost;
                inRatePowPost = sqrt(inRatePost.^2 + inPowPost.^2);
                inRatePowDiff_rftcHemisphere = inRatePowPre - inRatePowPost;
                zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').rftcHemisPre = inRatePowPre;
                zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').rftcHemisPost = inRatePowPost;
                zonesDataOut.(params.biomarker).(params.normalization).(params.patsSel).('RatePow').rftcHemisDiff = inRatePowDiff_rftcHemisphere;
            end
        end 
    end
end

function biomarkerActvDiff = getRFTC_Channels_ActivityDiff(params, groupTablePre, groupTablePost)
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
    biomarkerActvDiff = fillResultsStruct(params, groupTablePreIn, groupTablePostIn, groupTablePreOut, groupTablePostOut);
end

function biomarkerActvDiff = getAllChannels_ActivityDiff(params, groupTablePre, groupTablePost)
    groupTablePreIn = groupTablePre;
    groupTablePostIn = groupTablePost;
    groupTablePreOut = groupTablePre;
    groupTablePostOut = groupTablePost;
    
    %%
    groupTablePreIn = removeRFTC_Channels(groupTablePreIn);
    groupTablePostIn = removeRFTC_Channels(groupTablePostIn);
    groupTablePreOut = removeRFTC_Channels(groupTablePreOut);
    groupTablePostOut = removeRFTC_Channels(groupTablePostOut);
    
    %%
    biomarkerActvDiff = fillResultsStruct(params, groupTablePreIn, groupTablePostIn, groupTablePreOut, groupTablePostOut);
end

function biomarkerActvDiff = getHighEI_Channels_ActivityDiff(params, groupTablePre, groupTablePost)
    biomarkerActvDiff.in = [];
    biomarkerActvDiff.out = [];
    groupTablePreIn = groupTablePre;
    groupTablePostIn = groupTablePost;
    groupTablePreOut = groupTablePre;
    groupTablePostOut = groupTablePost;
    
    %%
    %tableDelIdxs = not(groupTablePre.eiVals > median(groupTablePre.eiVals) + std(groupTablePre.eiVals));
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
    biomarkerActvDiff = fillResultsStruct(params, groupTablePreIn, groupTablePostIn, groupTablePreOut, groupTablePostOut);
end

function biomarkerActvDiff = getRFTC_ElectroPhysioConnected_Channels_ActivityDiff(params, groupTablePre, groupTablePost)
    biomarkerActvDiff.in = [];
    biomarkerActvDiff.out = [];
    groupTablePreIn = groupTablePre;
    groupTablePostIn = groupTablePost;
    groupTablePreOut = groupTablePre;
    groupTablePostOut = groupTablePost;
    
    %%
    %tableDelIdxs = not(groupTablePre.rftcElectroPhysioConnect > median(groupTablePre.rftcElectroPhysioConnect) + std(groupTablePre.rftcElectroPhysioConnect));
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
    biomarkerActvDiff = fillResultsStruct(params, groupTablePreIn, groupTablePostIn, groupTablePreOut, groupTablePostOut);
end

function biomarkerActvDiff = getSameRFTC_Structure_ActivityDiff(params, groupTablePre, groupTablePost)
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
    biomarkerActvDiff = fillResultsStruct(params, groupTablePreIn, groupTablePostIn, groupTablePreOut, groupTablePostOut);
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
    rftcStructureIdx = sameRFTC_BrainParcelChannsPerPat(groupTablePre);
    tableDelIdxs = tableDelIdxs | rftcStructureIdx;
    
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
    biomarkerActvDiff = fillResultsStruct(params, groupTablePreIn, groupTablePostIn, groupTablePreOut, groupTablePostOut);;
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
    rftcLobeIdxs = sameRFTC_BrainLobeChannsPerPat(groupTablePre);
    rftcStructureIdx = sameRFTC_BrainParcelChannsPerPat(groupTablePre);
    tableDelIdxs = tableDelIdxs | rftcLobeIdxs | rftcStructureIdx;
    
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
    biomarkerActvDiff = fillResultsStruct(params, groupTablePreIn, groupTablePostIn, groupTablePreOut, groupTablePostOut);
end

function biomarkerActvDiff = fillResultsStruct(params, groupTablePreIn, groupTablePostIn, groupTablePreOut, groupTablePostOut)

    biomarkerActvDiff.in.pre = groupTablePreIn.(params.biomarker);
    biomarkerActvDiff.in.post = groupTablePostIn.(params.biomarker);
    
    biomarkerActvDiff.out.pre = groupTablePreOut.(params.biomarker);
    biomarkerActvDiff.out.post = groupTablePostOut.(params.biomarker);
    
    biomarkerActvDiff.in.diff = groupTablePreIn.(params.biomarker) - groupTablePostIn.(params.biomarker);
    biomarkerActvDiff.out.diff = groupTablePreOut.(params.biomarker) - groupTablePostOut.(params.biomarker);
    
    biomarkerActvDiff.in.prePostP = signrank(groupTablePreIn.(params.biomarker), groupTablePostIn.(params.biomarker), 'tail', 'right');
    biomarkerActvDiff.out.prePostP = signrank(groupTablePreOut.(params.biomarker), groupTablePostOut.(params.biomarker), 'tail', 'right');
    
    biomarkerActvDiff.in.nrChanns = length(groupTablePreIn.patName);
    biomarkerActvDiff.out.nrChanns = length(groupTablePreOut.patName);
    
    biomarkerActvDiff.in.nrPats = length(unique(groupTablePreIn.patName));
    biomarkerActvDiff.out.nrPats = length(unique(groupTablePreOut.patName));
    
    [p, h] = ranksum(biomarkerActvDiff.in.diff, biomarkerActvDiff.out.diff);
    biomarkerActvDiff.inOutP = p;
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

        patKeepIdxB = patVals >= prctile(patVals, th);
        patKeepIdxC = patVals >= median(patVals)+0.75*std(patVals);
        patKeepIdx = patKeepIdxC;
        
        allPatsKeepIdxs(patSelIdx) = patKeepIdx;        
    end
    allPatsKeepIdxs = logical(allPatsKeepIdxs);
end

function [groupTablePre, groupTablePost] = normalizePerPatient(groupTablePre, groupTablePost)

    groupTableAll = [];

    allPatNames = unique(groupTablePre.patName);
    nrPatients = length(allPatNames);
    columnNames = groupTablePre.Properties.VariableNames;

    for pi = 1:nrPatients
        patSelIdx = ismember(groupTablePre.patName, allPatNames{pi});
        patSelIdxPost = ismember(groupTablePost.patName, allPatNames{pi});
        if sum(patSelIdx & patSelIdxPost) ~= sum(patSelIdx)
            stop = 1;
        end
        for ci = 1:length(columnNames)
            column = columnNames{ci};
            columnPatValsPre = groupTablePre.(column)(patSelIdx);
            columnPatValsPost = groupTablePost.(column)(patSelIdx);
            columnPatValsAll = [columnPatValsPre columnPatValsPost];
            if length(columnPatValsPre) ~= length(columnPatValsPost)
                stop = 1;
            end

            if ci >= 9
                groupTablePre.(column)(patSelIdx) = (columnPatValsPre-mean(columnPatValsAll))/std(columnPatValsAll);
                groupTablePost.(column)(patSelIdx) = (columnPatValsPost-mean(columnPatValsAll))/std(columnPatValsAll);
            end
        end
    end

end

function [groupTablePre, groupTablePost] = readTables(paths, params)
    %% Read tables
    groupAnalysisTablesFilePath = strcat(paths.workspacePath, 'GroupAnalysis_ChannelCharacterizationTables\', params.selDetector, '\');
    preGroupTableFN = strcat(groupAnalysisTablesFilePath, 'GroupAnalysis_ChannelCharacterization_Pre', params.selDetector, '.xls');
    postGroupTableFN = strcat(groupAnalysisTablesFilePath, 'GroupAnalysis_ChannelCharacterization_Post', params.selDetector, '.xls');
    groupTablePre = readtable(preGroupTableFN, 'Sheet', params.feature);
    groupTablePost = readtable(postGroupTableFN, 'Sheet', params.feature);
    
    %groupTablePre = truncateOutliers(groupTablePre);
    %groupTablePost = truncateOutliers(groupTablePost);
    
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
        tableDelIdxsIn = groupTablePre.outcome < params.outcomeTh;
        groupTablePre(tableDelIdxsIn, :) = [];
        groupTablePost(tableDelIdxsIn, :) = [];        
    elseif strcmp(params.patsSel, 'nonImprovedPatients')
        tableDelIdxsIn = groupTablePre.outcome >= params.outcomeTh;
        groupTablePre(tableDelIdxsIn, :) = [];
        groupTablePost(tableDelIdxsIn, :) = []; 
    end
end

function groupTable = truncateOutliers(groupTable)

    allPatNames = unique(groupTable.patName);
    nrPatients = length(allPatNames);
    columnNames = groupTable.Properties.VariableNames;

    for pi = 1:nrPatients
        patSelIdx = ismember(groupTable.patName, allPatNames{pi});
        sum(patSelIdx)
        for ci = 1:length(columnNames)
            column = columnNames{ci};
            columnPatVals = groupTable.(column)(patSelIdx);
            if ci >= 9 %strcmp(class(columnPatVals), 'double') % ci >= 9
                meanVal = mean(columnPatVals);
                medianVal = median(columnPatVals);
                stdVal = std(columnPatVals);
                interQR = iqr(columnPatVals);
                quart1 = prctile(columnPatVals, 25);
                quart3 = prctile(columnPatVals, 75);

                lowBound = quart1 - 10*interQR;
                highBound = quart3 + 10*interQR;
                selLowOutliers = (columnPatVals < lowBound);
                selHighOutliers = (columnPatVals > highBound);
                if(sum(selLowOutliers) > 0)
                    [quart1 quart3 interQR meanVal medianVal]
                    columnPatVals(selHighOutliers)'
                    columnPatVals(selLowOutliers) = lowBound;
                end
                if(sum(selHighOutliers) > 0)
                    [quart1 quart3 interQR meanVal medianVal]
                    columnPatVals(selHighOutliers)'
                    columnPatVals(selHighOutliers) = highBound;
                end
                groupTable.(column)(patSelIdx) = columnPatVals;
            end
        end
    end
end