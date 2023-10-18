clc; clear all; close all;
paths = getFilesPaths();
preFiles  = getPreFiles();%getAllFiles, getPreFiles, getPostFiles
postFiles  = getPostFiles();%getAllFiles, getPreFiles, getPostFiles

selDetector = 'MOSSDET'; 'Delphos'; 'MOSSDET';
features = {'OccRate', 'Power', 'Frequency'};
normOptions = {'', 'Normalized', 'Scaled'};
freqBandConnList = {'delta', 'theta', 'alpha', 'beta', 'gamma', 'highGamma', 'ripple', 'fr', 'maxAllBands', 'meanAllBands'};
thList = 0:30;
thList = 0;

if strcmp(selDetector, 'Delphos')
    thList = 0;
end


for th = thList    
    for fi = 1:length(features)
        for ni = 1:length(normOptions)
            for fbi = 1:length(freqBandConnList)
                freqBandConn = freqBandConnList{fbi};
                feature = features{fi};
                normalization = normOptions{ni};

                allPatientsAllZonesTablePre = [];
                allPatientsAllZonesTablePost = [];
                for fileIdx = 1:size(preFiles,1)
                    selDetectorFolder = selDetector;
                    if th > 0
                        selDetectorFolder = 'MOSSDET_Depurated';
                    end

                    tablesFilePath = strcat(paths.workspacePath, 'ChannelCharacterizationTables\', selDetectorFolder, '\');

                    if th == 0
                        tablesFilePath = tablesFilePath;
                    else 
                        tablesFilePath = strcat(tablesFilePath, 'Th', num2str(th), '\');
                    end

                    prePatName = preFiles{fileIdx}; prePatName = prePatName(1:length(prePatName)-4);
                    postPatName = postFiles{fileIdx}; postPatName = postPatName(1:length(postPatName)-4);
                    preTableFN = strcat(tablesFilePath, prePatName, '_', 'ChannelCharacterization_', selDetector, '.xls');
                    postTableFN = strcat(tablesFilePath, postPatName, '_', 'ChannelCharacterization_', selDetector, '.xls');
                    preTable = readtable(preTableFN, 'Sheet', feature);
                    postTable = readtable(postTableFN, 'Sheet', feature);
                    %Complete RFTC values, not savd to Power and Freq sheets by mistake
                    preTableFixError = readtable(preTableFN, 'Sheet', 'OccRate');
                    postTableFixError = readtable(postTableFN, 'Sheet', 'OccRate');
                    preTable.rftcVals = preTableFixError.rftcVals;
                    postTable.rftcVals = postTableFixError.rftcVals;
                    
                    %%Anat Localization Table
                    anatLocalizationPath = 'F:\ForschungsProjekte\RFTC\Project_Files\ChannelLocalizationAtlases\';
                    anatLocalizationGrenoblePath = 'F:\ForschungsProjekte\RFTC\Project_Files\ChannelLocalization\';
                    anatLocTable = readtable(strcat(anatLocalizationPath, prePatName, '.csv'));
                    patNameSep = strfind(prePatName, '_');
                    tabPatName = prePatName(1:patNameSep(3)-1); tabPatName(1:3) = 'Gre';
                    anatLocTableGrenoble = readtable(strcat(anatLocalizationGrenoblePath, tabPatName, '.csv'));

                    patCodeSeparators = strfind(prePatName, '_');
                    patCode = prePatName(patCodeSeparators(2):patCodeSeparators(3));patCode = strrep(patCode, '_', '');
                    outcomeTableFilename = 'F:\ForschungsProjekte\RFTC\MATLAB\DetectHFO\OtherData\Lachner_DetectedFiles_List.xlsx';
                    outcomeTable = readtable(outcomeTableFilename, 'Sheet', 'MicromedFiles(.TRC)');
                    outcomeVal = outcomeTable.Post_RFTCImprovement___(find(ismember(outcomeTable.Code,patCode)));

                    %Delete channels which are in pre and not in post or the
                    %other way around, because we want to make a paired test
                    preChannCheck = not(ismember(preTable.channelLabels, postTable.channelLabels));
                    if sum(preChannCheck) > 0    
                        preTable(preChannCheck, :) = [];
                    end
                    postChannCheck = not(ismember(postTable.channelLabels, preTable.channelLabels));
                    if sum(postChannCheck) > 0    
                        postTable(postChannCheck, :) = [];
                    end
                    
                    %Rearrange channel order so that it matches in pre and
                    %post
                    newChannOrder = [];
                    for chi = 1:length(postTable.channelLabels)
                        tableChann = postTable.channelLabels{chi};
                        newChannOrder = cat(1, newChannOrder, find(ismember(preTable.channelLabels, tableChann)));
                    end
                    postTable = postTable(newChannOrder, :);

                    nrChannsPre = length(preTable.channelLabels);
                    nrChannsPost = length(postTable.channelLabels);
                    if not(nrChannsPre == nrChannsPost)
                        errChann = 1;
                    end
                    % Add variables to table
                    patNameVecPre = preTable.channelLabels;
                    patNameVecPost = postTable.channelLabels;
                    outcomePre = preTable.rftcVals;
                    outcomePost = postTable.rftcVals;
                    rftcElectroPhysioConnectPre = preTable.rftcVals;
                    rftcElectroPhysioConnectPost = postTable.rftcVals;
                    hemispherePre = preTable.channelLabels;
                    hemispherePost = postTable.channelLabels;
                    brainParcelPre = preTable.channelLabels;
                    brainParcelPost = postTable.channelLabels;

                    for chi = 1:nrChannsPre
                        if not(strcmp(preTable.channelLabels{chi}, postTable.channelLabels{chi}))
                            errChann = 1;
                        end
                        
                        chName = preTable.channelLabels{chi};
                        patNameVecPre{chi} = patCode;
                        patNameVecPost{chi} = patCode;
                        outcomePre(chi) = outcomeVal;
                        outcomePost(chi) = outcomeVal;

                        rftcConnect = getMax_RFTC_ElectroPhysioConnect(paths, prePatName, preTable.channelLabels{chi}, freqBandConn);
                        rftcElectroPhysioConnectPre(chi) = rftcConnect;
                        rftcElectroPhysioConnectPost(chi) = rftcConnect;
                        
                        isLeftSide = strfind(chName, '''');
                        hemisSide = 'right';
                        if isLeftSide
                            hemisSide = 'left';
                        end
                        hemispherePre{chi} = hemisSide;
                        hemispherePost{chi} = hemisSide;
                        brainParcelPre{chi} = getBrainParcel(anatLocTable, chName);                    
                        brainParcelPost{chi} = brainParcelPre{chi};
                    end

                    patName = patNameVecPre; preTable = addvars(preTable, patName, 'Before','channelLabels');
                    patName = patNameVecPost; postTable = addvars(postTable, patName, 'Before','channelLabels');
                    outcome = outcomePre; preTable = addvars(preTable, outcome, 'Before','channelLabels');
                    outcome = outcomePost; postTable = addvars(postTable, outcome, 'Before','channelLabels');
                    rftcElectroPhysioConnect = rftcElectroPhysioConnectPre; preTable = addvars(preTable, rftcElectroPhysioConnect, 'After','rftcVals');
                    rftcElectroPhysioConnect = rftcElectroPhysioConnectPost; postTable = addvars(postTable, rftcElectroPhysioConnect, 'After','rftcVals');
                    hemisphere = hemispherePre; preTable = addvars(preTable, hemisphere, 'After','rftcElectroPhysioConnect');
                    hemisphere = hemispherePost; postTable = addvars(postTable, hemisphere, 'After','rftcElectroPhysioConnect');
                    brainParcel = brainParcelPre; preTable = addvars(preTable, brainParcel, 'After','hemisphere');
                    brainParcel = brainParcelPost; postTable = addvars(postTable, brainParcel, 'After','hemisphere');

                    preTable = fillEmptySpaces(preTable);
                    postTable = fillEmptySpaces(postTable);

                    if strcmp(normalization, 'Normalized')
                        normPreTable = normalizePatientTable(preTable);
                        normPostTable = normalizePatientTable(postTable);
                    elseif strcmp(normalization, 'Scaled')
                        normPreTable = scalePatientTable(preTable);
                        normPostTable = scalePatientTable(postTable);
                    else
                        normPreTable = preTable;
                        normPostTable = postTable;
                    end

                    patCode
                    if isempty(allPatientsAllZonesTablePre)
                        allPatientsAllZonesTablePre = normPreTable;
                        allPatientsAllZonesTablePost = normPostTable;
                    else
                        allPatientsAllZonesTablePre = cat(1, allPatientsAllZonesTablePre, normPreTable);
                        allPatientsAllZonesTablePost = cat(1, allPatientsAllZonesTablePost, normPostTable);            
                    end
                    size(normPreTable, 1)
                    size(normPostTable, 1)
                end

                if th == 0
                    selDetectorFolder = strcat(selDetector, '\');
                else
                    selDetectorFolder = 'MOSSDET_Depurated\';
                    selDetectorFolder = strcat(selDetectorFolder, 'Th', num2str(th), '\');
                end


                groupAnalysisTablesFilePath = strcat(paths.workspacePath, 'GroupAnalysis_ChannelCharacterizationTables\', selDetectorFolder); mkdir(groupAnalysisTablesFilePath);
                preGroupTableFN = strcat(groupAnalysisTablesFilePath, normalization, 'GroupAnalysis_ChannelCharacterization_Pre', selDetector, '_', freqBandConn, '.xls');
                postGroupTableFN = strcat(groupAnalysisTablesFilePath, normalization, 'GroupAnalysis_ChannelCharacterization_Post', selDetector, '_', freqBandConn, '.xls');

                writetable(allPatientsAllZonesTablePre, preGroupTableFN, 'Sheet', feature);
                writetable(allPatientsAllZonesTablePost, postGroupTableFN, 'Sheet', feature);
                
                if fbi == 1
                    preGroupTableFN = strcat(groupAnalysisTablesFilePath, normalization, 'GroupAnalysis_ChannelCharacterization_Pre', selDetector, '.xls');
                    postGroupTableFN = strcat(groupAnalysisTablesFilePath, normalization, 'GroupAnalysis_ChannelCharacterization_Post', selDetector, '.xls');
                    writetable(allPatientsAllZonesTablePre, preGroupTableFN, 'Sheet', feature);
                    writetable(allPatientsAllZonesTablePost, postGroupTableFN, 'Sheet', feature);
                end
%                 normalization
%                 feature
%                 size(allPatientsAllZonesTablePre)
%                 size(allPatientsAllZonesTablePost)
            end
        end
    end
end

function table = normalizePatientTable(table)
    columns = table.Properties.VariableNames;
    for ci = 1:length(columns)
        cName = columns{ci};                    
        if strcmp(cName, 'patName') || strcmp(cName, 'outcome') || strcmp(cName, 'channelLabels') || strcmp(cName, 'rftcVals') || strcmp(cName, 'hemisphere') || strcmp(cName, 'brainParcel')
            continue;
        end
        table.(cName) = (table.(cName) - mean(table.(cName)))/std(table.(cName));
    end
end

function table = scalePatientTable(table)
    minPrctile = 5;
    maxPrctile = 95;
    columns = table.Properties.VariableNames;
    %%Truncate outliers
    for ci = 1:length(columns)
        cName = columns{ci};
        if strcmp(cName, 'patName') || strcmp(cName, 'outcome') || strcmp(cName, 'channelLabels') || strcmp(cName, 'rftcVals') || strcmp(cName, 'hemisphere') || strcmp(cName, 'brainParcel')
            continue;
        end        
        %Truncate outliers
        colVals = table.(cName);
        sel = colVals < prctile(table.(cName), minPrctile);
        colVals(sel) = prctile(table.(cName), minPrctile);

        sel = colVals > prctile(table.(cName), maxPrctile);
        colVals(sel) = prctile(table.(cName), maxPrctile);
        table.(cName) = colVals;       
    end
    
    for ci = 1:length(columns)
        cName = columns{ci};
        if strcmp(cName, 'patName') || strcmp(cName, 'outcome') || strcmp(cName, 'channelLabels') || strcmp(cName, 'rftcVals') || strcmp(cName, 'hemisphere') || strcmp(cName, 'brainParcel')
            continue;
        end
        %table.(cName) = (table.(cName) - prctile(table.(cName), minPrctile))/(prctile(table.(cName), maxPrctile) - prctile(table.(cName), minPrctile));
        table.(cName) = (table.(cName) - min(table.(cName)))/(max(table.(cName)) - min(table.(cName)));
        
        if sum(table.(cName) < 0) > 0 || sum(table.(cName) > 1) > 0
            stopHere = 1;
        end
    end
end

function table = fillEmptySpaces(table)
    tableVars = table.Properties.VariableNames;
    for ci = 1:length(tableVars)
        varName = tableVars{ci};
        vec = table(:, varName);
        table(find(ismissing(vec)), varName) = {0};
    end
end

function parcelName = getBrainParcel(anatLocTable, chName)
    if contains(chName, '-')
        chName = chName(1:strfind(chName, '-')-1);
    end
    
    chIdx = find(ismember(anatLocTable.Channel, chName));
    if isempty(chIdx)
        parcelName = 'N/A';
    else
        parcelName = anatLocTable.Harvard_OxfordSubcorticalStructuralAtlas{chIdx};
    end
end