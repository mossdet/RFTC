clc; clear all; close all;

paths = getFilesPaths();
eegFileNames  = getAllFiles();
preFiles  = getPreFiles();
biomarkersList =  {'HFO'};
featuresList = {'rate'};

patTablesPath = strcat(paths.workspacePath, 'CharacterizationTables\FeatureCharacterizationTables\');
groupTablesPath = strcat(paths.workspacePath, 'CharacterizationTables\GroupCharacterizationTablesMedian\');
analysisTablesPath = strcat(paths.workspacePath, 'AnalysisResults\Median\');mkdir(analysisTablesPath);
normStrList = {''};
channSelStrList = {'FlexK'};

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
        tableAllPatients_SpreadSheetName = strcat(analysisTablesPath, 'allPats_Pre_vs_Post_Analysis_', channSelStr, '_', normStr,'.xls');

        for bm = biomarkersList
            allGroupAnalysisP = {};
            for ft = featuresList
                feature = ft{1};
                biomarker = bm{1};

                groupTablePre = readtable(spreadSheetNamePre, 'Sheet', biomarker);
                groupTablePost = readtable(spreadSheetNamePost, 'Sheet', biomarker);
                featsNames = groupTablePre.Properties.VariableNames;
                patNames = unique(groupTablePre.patName);
                patNamesPost = unique(groupTablePost.patName);
                
                patientsInfoTableHdr = {"PatName", "PatNr", "Sex", "SamplingRate", "SeegDuration", "BilateralImplantation", "LobeElectrodeLocation", "RFTC Electrode Lobe Location", "Post-RFTC Improvement"};
                patientsInfoTable = {};
                
                for patIdx = 1:length(patNames)
                    patName = patNames{patIdx};
                    patNamePost = patNamesPost{patIdx};
                    patNameDivs = strfind(patName, '_');
                    patNameShort = patName(patNameDivs(2)+1:patNameDivs(3)-1);
                    
                    patTableSel = strcmp(groupTablePre.patName, patName);
                    patTable = groupTablePre(patTableSel,:);
                    
                    if not(contains(patNamePost, patNameShort))
                       error("Wrong pre- post- file name allignment");
                    end
                    
                    % Pre
                    eegFileIdx = find(contains(eegFileNames, patName));
                    eegFilename = strcat(paths.eegFilesPath, eegFileNames{eegFileIdx});
                    rftcPatData = loadRFTC_Data(paths.workspacePath, eegFilename);
                    samplingRate = rftcPatData.fs;
                    durationMinPre = round(rftcPatData.nrSamples/rftcPatData.fs/60);
                    header = ft_read_header(eegFilename);
                    
                    % Post
                    eegFileIdxPost = find(contains(eegFileNames, patNamePost));
                    eegFilenamePost = strcat(paths.eegFilesPath, eegFileNames{eegFileIdxPost});
                    rftcPatDataPost = loadRFTC_Data(paths.workspacePath, eegFilenamePost);
                    samplingRatePost = rftcPatDataPost.fs;
                    durationMinPost = round(rftcPatDataPost.nrSamples/rftcPatDataPost.fs/60);
                                        
                    channelLobeLoc = getPatientBrainLobes(patTable);
                    
                    durationStr = strcat("Pre-RFTC:", num2str(durationMinPre),", Post-RFTC:", num2str(durationMinPost));
                    bilateralContcatsNrStr = strcat("Left (", num2str(sum(contains(patTable.chName, ''''))), "), Right (", num2str(sum(not(contains(patTable.chName, '''')))), ")");
                    
                    %RFTC Channels
                    patTableRFTC = patTable(patTable.rftcSite>0,:);
                    channelLobeLocRFTC = getPatientBrainLobes(patTableRFTC);
                    
                    improvementStr = patTable.outcome(1);
                    
                    patEntree = {patNameShort, patIdx, "Sex", samplingRate, durationStr, bilateralContcatsNrStr, channelLobeLoc, channelLobeLocRFTC, improvementStr};
                    
                    patientsInfoTable = cat(1, patientsInfoTable, patEntree);
                end
            end
            [sortedNames, sortIdx] = sort(patientsInfoTable(:,1));
            patientsInfoTable = patientsInfoTable(sortIdx, :);
            
            % All Patients
            cell2table(patientsInfoTable, 'VariableNames', [patientsInfoTableHdr{1,:}]);
        end
    end
end



