clc; clear all; close all;

paths = getFilesPaths();
preFiles  = getPreFiles();%getAllFiles, getPreFiles, getPostFiles
postFiles  = getPostFiles();%getAllFiles, getPreFiles, getPostFiles

biomarkersList =  {'HFO', 'iesHFO', 'IES'};
patTablesPath = strcat(paths.workspacePath, 'CharacterizationTables\FeatureCharacterizationTablesAvg\');
groupTablesPath = strcat(paths.workspacePath, 'CharacterizationTables\GroupCharacterizationTablesAvg\'); mkdir(groupTablesPath);

channSelStrList = {'Threshold', 'FlexK', 'Cluster2K'};
channSelStrList = {'FlexK'};

for channSelIdx = 1:length(channSelStrList)
    channSelStr = channSelStrList{channSelIdx};
    spreadSheetNamePre = strcat(groupTablesPath, 'GroupCharacterization_', channSelStr, '_Pre.xls');
    delete(spreadSheetNamePre); 
    spreadSheetNamePost = strcat(groupTablesPath, 'GroupCharacterization_', channSelStr, '_Post.xls');
    delete(spreadSheetNamePost);
    spreadSheetNamePreNormalized = strcat(groupTablesPath, 'GroupCharacterization_', channSelStr, '_Pre_Normalized.xls');
    delete(spreadSheetNamePreNormalized); 
    spreadSheetNamePostNormalized = strcat(groupTablesPath, 'GroupCharacterization_', channSelStr, '_Post_Normalized.xls');
    delete(spreadSheetNamePostNormalized);

    for bm = biomarkersList
        biomarker = bm{1};            
        groupTablePre = {};
        groupTablePost = {};
        groupTablePreNorm = {};
        groupTablePostNorm = {};
        
        highEI_NrChannsPerPat = {'Pat', 'Th_NrChs', 'Cluster_NrChs', 'FlexK_NrChs'};
        rftcConn_NrChannsPerPat = {'Pat', 'Th_NrChs', 'Cluster_NrChs', 'FlexK_NrChs'};

        usedAtlasesStruct = {};

        for fileIdx = 1:size(preFiles,1)
            eegFilenamePre = strcat(paths.eegFilesPath, preFiles{fileIdx})
            eegFilenamePost = strcat(paths.eegFilesPath, postFiles{fileIdx})
            [origFilepath, patNamePre, ext] = fileparts(eegFilenamePre);
            [origFilepath, patNamePost, ext] = fileparts(eegFilenamePost);

            % Pre Table
            spreadSheetName = strcat(patTablesPath, patNamePre, '.xls');
            table_Pre = readtable(spreadSheetName, 'Sheet', biomarker);
            % Post Table
            spreadSheetName = strcat(patTablesPath, patNamePost, '.xls');
            table_Post = readtable(spreadSheetName, 'Sheet', biomarker);

            patNr = cell(length(table_Pre.chNr),1); patNr(:,1) = {fileIdx};
            patNr = cell2table(patNr);
            table_Pre = [patNr table_Pre];
            patNr = cell(length(table_Post.chNr),1); patNr(:,1) = {fileIdx};
            patNr = cell2table(patNr);
            table_Post = [patNr table_Post];

            % Select only channels present in pre and post
            selPre = ismember(table_Pre.chName, table_Post.chName);
            selPost = ismember(table_Post.chName, table_Pre.chName);
            table_Pre = table_Pre(selPre, :);
            table_Post = table_Post(selPost, :);
            %Sort channels
            [vec, sortIdxPre] = sort(table_Pre.chName);
            table_Pre = table_Pre(sortIdxPre, :);
            [vec, sortIdxPost] = sort(table_Post.chName);
            table_Post = table_Post(sortIdxPost, :);

            if sum(strcmp(table_Pre.chName, table_Post.chName)) < length(table_Pre.chName)
                error('Pre and Post EEG Channels are different')
            end

            % Add zone
            zoneVal = getBrainZone(paths, patNamePre, table_Pre.chName);
            table_Pre.zone = zoneVal;
            table_Post.zone = zoneVal;
            % Add outcome
            outcomeVal = getOutcome(patNamePre);
            table_Pre.outcome(:) = outcomeVal;
            table_Post.outcome(:) = outcomeVal;
            % Add EI
            eiVals = getEI(eegFilenamePre, table_Pre.chName);
            table_Pre.eiVals = eiVals;
            table_Post.eiVals = eiVals;
            % Add highEI Flags
            highEI_FlagsTh = thresholdVector(eiVals);
            highEI_FlagsCl = clusterVector(eiVals);
            highEI_FlagsClFlexK = clusterVector_flexK(eiVals);
            highEI_Flags = highEI_FlagsTh;
            if strcmp(channSelStr, 'Cluster2K')
                highEI_Flags = highEI_FlagsCl;
            elseif strcmp(channSelStr, 'FlexK')
                highEI_Flags = highEI_FlagsClFlexK;
            end
            table_Pre.highEI = highEI_Flags;
            table_Post.highEI = highEI_Flags;
            % Add RFTC Connectivity
            rftcConn = getRFTC_Connectivity(paths, patNamePre, table_Pre.chName);
            table_Pre.rftcElectroPhysioConnect = rftcConn;
            table_Post.rftcElectroPhysioConnect = rftcConn;
            % Add highConn Flags
            rftcConnectFlags_Th = thresholdVector(rftcConn);
            rftcConnectFlags_Cl = clusterVector(rftcConn);
            rftcConnectFlagsClFlexK = clusterVector_flexK(rftcConn);
            rftcConnectFlags = rftcConnectFlags_Th;
            if strcmp(channSelStr, 'Cluster2K')
                rftcConnectFlags = rftcConnectFlags_Cl;
            elseif strcmp(channSelStr, 'FlexK')
                rftcConnectFlags = rftcConnectFlagsClFlexK;
            end
            table_Pre.rftcConn = rftcConnectFlags;
            table_Post.rftcConn = rftcConnectFlags;

            % Add RFTC Site Flag       
            rftcSiteLabels = getRFTC_Site(paths, patNamePre, table_Pre.chName);
            table_Pre.rftcSite = rftcSiteLabels;
            table_Post.rftcSite = rftcSiteLabels;
            % Add RFTC Hemisphere Flag
            rftcHemisphereFlags = getRFTC_Hemisphere(table_Pre.chName, table_Pre.rftcSite);
            rftcHemisphereFlagsCorr = rftcHemisphereFlags | rftcSiteLabels;
            table_Pre.rftcHemis = rftcHemisphereFlagsCorr;
            table_Post.rftcHemis = rftcHemisphereFlagsCorr;
            % Add RFTC Lobe Flag
            rftcLobeFlags = getRFTC_Lobe(table_Pre);
            rftcLobeFlagsCorr = rftcLobeFlags & rftcHemisphereFlags;
            table_Pre.rftcLobe= rftcLobeFlagsCorr;
            table_Post.rftcLobe = rftcLobeFlagsCorr;
            % Add RFTC Struct Flag
            [rftcStructFlags, usedAtlasName] = getRFTC_Struct(table_Pre);
            rftcStructFlagsCorr = rftcStructFlags & rftcHemisphereFlags;
            table_Pre.rftcStruct = rftcStructFlagsCorr;
            table_Post.rftcStruct = rftcStructFlagsCorr;

            usedAtlasesStruct = cat(1, usedAtlasesStruct, usedAtlasName);

            clc;
            patNamePre
            nrChs = length(highEI_FlagsTh);
            {'Zone', 'th_chNr', 'cluster_chNr', 'cluster4k_chNr';...
              'highEI', sum(highEI_FlagsTh)/nrChs*100, sum(highEI_FlagsCl)/nrChs*100, sum(highEI_FlagsClFlexK)/nrChs*100;...
              'highConn', sum(rftcConnectFlags_Th)/nrChs*100, sum(rftcConnectFlags_Cl)/nrChs*100, sum(rftcConnectFlagsClFlexK)/nrChs*100}
            %[sum(rftcHemisphereFlags) sum(rftcHemisphereFlagsCorr); sum(rftcLobeFlags) sum(rftcLobeFlagsCorr); sum(rftcStructFlags) sum(rftcStructFlagsCorr)]

            highEI_NrChannsPerPat = cat(1, highEI_NrChannsPerPat, {patNamePre, sum(highEI_FlagsTh)/nrChs*100, sum(highEI_FlagsCl)/nrChs*100, sum(highEI_FlagsClFlexK)/nrChs*100});
            rftcConn_NrChannsPerPat = cat(1, rftcConn_NrChannsPerPat, {patNamePre, sum(rftcConnectFlags_Th)/nrChs*100, sum(rftcConnectFlags_Cl)/nrChs*100, sum(rftcConnectFlagsClFlexK)/nrChs*100});

            groupTablePre = cat(1, groupTablePre, table_Pre);
            groupTablePost = cat(1, groupTablePost, table_Post);

            [normTable_Pre, normTable_Post] = normalizeFeatsInTables(table_Pre, table_Post);
            groupTablePreNorm = cat(1, groupTablePreNorm, normTable_Pre);
            groupTablePostNorm = cat(1, groupTablePostNorm, normTable_Post);
        end
        usedAtlasesStruct

        sum(strcmp(groupTablePre.chName, groupTablePost.chName)) == length(groupTablePost.chName)
        writetable(groupTablePre, spreadSheetNamePre, 'Sheet', biomarker);
        writetable(groupTablePost, spreadSheetNamePost, 'Sheet', biomarker);
        writetable(groupTablePreNorm, spreadSheetNamePreNormalized, 'Sheet', biomarker);
        writetable(groupTablePostNorm, spreadSheetNamePostNormalized, 'Sheet', biomarker);

    end
end

function zoneVals = getBrainZone(paths, patName, channelLabels)

    anatLocalizationPath = 'F:\ForschungsProjekte\RFTC\Project_Files\ChannelLocalizationAtlases\';
    anatLocalizationGrenoblePath = 'F:\ForschungsProjekte\RFTC\Project_Files\ChannelLocalization\';
    anatLocTable = readtable(strcat(anatLocalizationPath, patName, '.csv'));
    patNameSep = strfind(patName, '_');
    tabPatName = patName(1:patNameSep(3)-1); tabPatName(1:3) = 'Gre';
    anatLocTableGrenoble = readtable(strcat(anatLocalizationGrenoblePath, tabPatName, '.csv'));
    
    zoneVals = cell(length(channelLabels),1);
    
    for chi = 1:length(channelLabels)
        chName = channelLabels{chi};
        if contains(chName, '-')
            chName = chName(1:strfind(chName, '-')-1);
        end
        chIdx = find(ismember(anatLocTable.Channel, chName));
        if isempty(chIdx)
            parcelName = 'N/A';
        else
            parcelName = anatLocTable.Harvard_OxfordSubcorticalStructuralAtlas{chIdx};
        end
        zoneVals{chi,1} = parcelName;
    end
end

function outcomeVal = getOutcome(patName)
    patCodeSeparators = strfind(patName, '_');
    patCode = patName(patCodeSeparators(2):patCodeSeparators(3));patCode = strrep(patCode, '_', '');
    outcomeTableFilename = 'F:\ForschungsProjekte\RFTC\RFTC_HFO\OtherData\Lachner_DetectedFiles_List.xlsx';
    outcomeTable = readtable(outcomeTableFilename, 'Sheet', 'MicromedFiles(.TRC)');
    outcomeVal = outcomeTable.Post_RFTCImprovement___(find(ismember(outcomeTable.Code,patCode)));
end

function eiVals = getEI(eegFilename, channelLabels)
    [origFilepath, patName, ext] = fileparts(eegFilename);
    [eiChannels, eiVals] = readEI(strcat(patName, ext));
    
    % Select only channels present in pre and post
    selChann = ismember(eiChannels, channelLabels);
    eiChannels = eiChannels(selChann);
    eiVals = eiVals(selChann);
    %Sort channels
    [eiChannels, sortIdx] = sort(eiChannels);
    eiVals = eiVals(sortIdx);
    
    if sum(strcmp(eiChannels, channelLabels)) < length(channelLabels)
        error('EI Channels and EEG channels are different')
    end
    
    eiVals = cell2mat(eiVals);
end

function rftcConnect = getRFTC_Connectivity(paths, patName, channelLabels)
    rftcConnect = zeros(length(channelLabels),1);
    for chi = 1:length(channelLabels)
        channName = channelLabels{chi};
        rftcConnect(chi) = getMax_RFTC_ElectroPhysioConnect(paths, patName, channName, 'maxAllBands');
    end
end

function rftcVals = getRFTC_Site(paths, patName, channelLabels)
    rftcPatName = patName(1:strfind(patName, 'Inter')-2);
    rftcFilename = strcat(paths.rftcFlags, rftcPatName, '_RFTC_Channels.xls');
    readTable = readtable(rftcFilename);
    rftcChannelLabels = readTable.channelLabels;
    rftcVals = readTable.rftcVals;
    
    % Select only channels present in pre and post
    selChann = ismember(rftcChannelLabels, channelLabels);
    rftcChannelLabels = rftcChannelLabels(selChann);
    rftcVals = rftcVals(selChann);
    %Sort channels
    [rftcChannelLabels, sortIdx] = sort(rftcChannelLabels);
    rftcVals = rftcVals(sortIdx);
    rftcVals = rftcVals > 0;
    
    if sum(strcmp(rftcChannelLabels, channelLabels)) < length(channelLabels)
        error('RFTC Channels and EEG channels are different')
    end    
end

function rftcHemisphereFlags = getRFTC_Hemisphere(channelLabels, rftcSiteFlags)

    % get rftc hemisphere
    rftcHemisphereList = {};
    rftcChanns = channelLabels(rftcSiteFlags>0);
    for chir = 1:length(rftcChanns)
        rftcChann = rftcChanns{chir};
        hemisSide = 'right';
        isLeftSide = strfind(rftcChann, '''');
        if not(isempty(isLeftSide))
            hemisSide = 'left';
        end
        rftcHemisphereList =  cat(1, rftcHemisphereList, {hemisSide});
    end
    rftcHemisphereList = unique(rftcHemisphereList);
    
    % set flags
    rftcHemisphereFlags = false(length(channelLabels),1);
    for chi = 1:length(channelLabels)
        channName = channelLabels{chi};
        hemisSide = 'right';
        isLeftSide = strfind(channName, '''');
        if not(isempty(isLeftSide))
            hemisSide = 'left';
        end
        
        rftcHemisphereFlags(chi) = ismember(hemisSide, rftcHemisphereList);
    end
end

function [rftcStructFlags, usedAtlasName] = getRFTC_Struct(patTable)
    [rftcStructFlags, usedAtlasName] = sameRFTC_BrainStructChanns(patTable);
end

function rftcLobeFlags = getRFTC_Lobe(patTable)
    rftcLobeFlags = sameRFTC_BrainLobeChanns(patTable);
end

function clstrIdx = thresholdVector(vec)
    
    %th = median(vec) + std(vec);
    th = median(vec) + 0.5*std(vec);
    clstrIdx = vec > th;
    %clstrIdx = vec > prctile(vec, th);
    
    [length(vec(clstrIdx==0)) length(vec(clstrIdx==1))]
end

function clstrIdx = clusterVector(vec)
    
    [clstrIdx,C] = kmeans(vec, 2,'Distance','cityblock', 'Replicates', 1000, 'MaxIter',10000);
    clstrIdx = clstrIdx-1;
    clstrIdx = clstrIdx > 0;
    if mean(vec(clstrIdx==0)) > mean(vec(clstrIdx==1))
        clstrIdx = not(clstrIdx);
        CA = C(1);
        C(1) = C(2);
        C(2) = CA;
    end

    if mean(vec(clstrIdx==0)) > mean(vec(clstrIdx==1))
        error('cluster 1 has a lower mean')
    end

    [length(vec(clstrIdx==0)) length(vec(clstrIdx==1))]
    [mean(vec(clstrIdx==0)) mean(vec(clstrIdx==1))]
    
%     opts = statset('Display','final');
%     %[clstrIdx,C] = kmeans(vec, 2,'Distance','cityblock', 'Replicates', 10000, 'MaxIter',10000, 'Options',opts);
%     [clstrIdx,C] = kmeans(vec, 2,'Distance','cityblock', 'Replicates', 1000, 'MaxIter',10000);
%     [maxC, maxC_Idx] = max(C);
%     
%     clstrIdx(clstrIdx ~= maxC_Idx) = 0;
%     clstrIdx(clstrIdx == maxC_Idx) = 1;
%     
%     [length(vec(clstrIdx==0)) length(vec(clstrIdx==1))]
%     [mean(vec(clstrIdx==0)) mean(vec(clstrIdx==1))]
end

function [normTable_Pre, normTable_Post] = normalizeFeatsInTables(table_Pre, table_Post)

    columnNames = table_Pre.Properties.VariableNames;
    prePlusPost_Table = cat(1, table_Pre, table_Post);

    for ci = 15:length(columnNames)
        clName = columnNames{ci};
        vec = prePlusPost_Table.(clName);
        prePlusPost_Table.(clName) = (vec - mean(vec))/std(vec);
    end
    normTable_Pre = prePlusPost_Table(1:size(table_Pre,1),:);
    normTable_Post = prePlusPost_Table(size(table_Pre,1)+1:end,:);
end
