clc; clear all; close all;

paths = getFilesPaths();
files  = getFaultyFiles();
files  = getAllFiles();

for fileIdx = 1:size(files,1)

    eegFilename = strcat(paths.eegFilesPath, files{fileIdx})
    
    rftcPatData = loadRFTC_Data(paths.workspacePath, eegFilename);
    
    rftcPatData = deleteArfeactualChannels(rftcPatData);

    [origFilepath, patName, ext] = fileparts(rftcPatData.origFilename);
    nrChanns = rftcPatData.nrChanns;
    

    channCharactHFO = {};
    channCharactIES_HFO = {};
    channCharactIES = {};
    
    for chi = 1:nrChanns
        chName = rftcPatData.channsLabels{chi};
        filePath = strcat(paths.workspacePath, 'HFO_Detections_MOSSDET\', patName, '\', chName, '\');        
        loadFN = strcat(filePath, chName, '_hfoCharacterizationShort');
        load(loadFN, 'hfoFeatsShort');
        hdrFeats = hfoFeatsShort(1,:); hdrFeats{7} = 'zone';
        hfoFeatsShort(1,:) = []; 
        
        if chi == 1
           channCharactHFO =  hdrFeats;
           channCharactIES_HFO =  hdrFeats;
           channCharactIES =  hdrFeats;
        end
        
        patFilename = hfoFeatsShort{1,1};
        chName = hfoFeatsShort{1,2};
        type = cell2mat(hfoFeatsShort(:,6));
        zone = cell2mat(hfoFeatsShort(:,7));
               
        nrMins = rftcPatData.nrSamples / rftcPatData.fs / 60;
        selHFO = type == 1 | type == 2;
        selIES_HFO = type == 4 | type == 5;
        selIES = type == 3;
        
        rateHFO = sum(selHFO) / nrMins;
        rateIES_HFO = sum(selIES_HFO) / nrMins;
        rateIES = sum(selIES) / nrMins;
        
        channCharactHFO(chi+1,1) = {patFilename};
        channCharactIES_HFO(chi+1,1) = {patFilename};
        channCharactIES(chi+1,1) = {patFilename};
        
        channCharactHFO(chi+1,2) = {chName};
        channCharactIES_HFO(chi+1,2) = {chName};
        channCharactIES(chi+1,2) = {chName};
        
        channCharactHFO(chi+1,3) = {chi};
        channCharactIES_HFO(chi+1,3) = {chi};
        channCharactIES(chi+1,3) = {chi};
        
        channCharactHFO(chi+1,6) = {1};
        channCharactIES_HFO(chi+1,6) = {2};
        channCharactIES(chi+1,6) = {3};
        
        channCharactHFO(chi+1,7) = {0};
        channCharactIES_HFO(chi+1,7) = {0};
        channCharactIES(chi+1,7) = {0};
        
        channCharactHFO(chi+1,8) = {rateHFO};
        channCharactIES_HFO(chi+1,8) = {rateIES_HFO};
        channCharactIES(chi+1,8) = {rateIES};
       
        % get average feature values across all detections
        for fi = 9:length(hdrFeats)
            featVals =  cell2mat(hfoFeatsShort(:,fi));

            metricFeatHFO = getVecMetric(featVals(selHFO));
            metricFeatIES_HFO = getVecMetric(featVals(selIES_HFO));
            metricFeatIES = getVecMetric(featVals(selIES));
                        
            channCharactHFO(chi+1,fi) = {metricFeatHFO};
            channCharactIES_HFO(chi+1,fi) = {metricFeatIES_HFO};
            channCharactIES(chi+1,fi) = {metricFeatIES};
        end                  
    end
    
    %% Delete artefactual and NaN channels
    channIdx = 2;
    varianceIdx = 12;
    delChannsList = [];
    
    varianceVals = cell2mat(channCharactHFO(2:end, varianceIdx));
    nanIdx = find(isnan(varianceVals))+1;
    delChannsList = [delChannsList channCharactHFO(nanIdx, channIdx)];
    
    varianceVals = cell2mat(channCharactIES_HFO(2:end, varianceIdx));
    nanIdx = find(isnan(varianceVals))+1;
    delChannsList = [delChannsList channCharactIES_HFO(nanIdx, channIdx)];
    
    varianceVals = cell2mat(channCharactIES(2:end, varianceIdx));
    nanIdx = find(isnan(varianceVals))+1;
    delChannsList = [delChannsList channCharactIES(nanIdx, channIdx)];
    
    delIdx_HFO = find(contains(channCharactHFO(:, channIdx), delChannsList));
    channCharactHFO(delIdx_HFO,:) = [];
    delIdx_iesHFO = find(contains(channCharactIES_HFO(:, channIdx), delChannsList));
    channCharactIES_HFO(delIdx_HFO,:) = [];
    delIdx_IES = find(contains(channCharactIES(:, channIdx), delChannsList));
    channCharactIES(delIdx_HFO,:) = [];

    if sum(strcmp(channCharactHFO(:, channIdx), channCharactIES_HFO(:, channIdx))) < size(channCharactHFO,1)
        stop = 1;
    elseif sum(strcmp(channCharactIES_HFO(:, channIdx), channCharactIES(:, channIdx))) < size(channCharactHFO,1)
        stop = 1;
    end

    % Delete Sample Entrees
    channCharactHFO(:, 4:6) = []; channCharactIES_HFO(:, 4:6) = []; channCharactIES(:, 4:6) = [];

    %% Append fields for new data
    newDataHdr = {'outcome', 'eiVals', 'rftcElectroPhysioConnect', 'rftcSite', 'rftcConn', 'rftcStruct', 'highEI', 'rftcLobe', 'rftcHemis'};
    newDataCell = cell(size(channCharactHFO,1), length(newDataHdr)); 
    newDataCell(1, :) = newDataHdr;
    newDataCell(2:end, :) = {0};
    channCharactHFO = [channCharactHFO(:,1:4), newDataCell, channCharactHFO(:,5:end)];
    channCharactIES_HFO = [channCharactIES_HFO(:,1:4), newDataCell, channCharactIES_HFO(:,5:end)];
    channCharactIES = [channCharactIES(:,1:4), newDataCell, channCharactIES(:,5:end)];
    
    tablesPath = strcat(paths.workspacePath, 'CharacterizationTables\FeatureCharacterizationTables\'); mkdir(tablesPath);
    
    table_HFO = cell2table(channCharactHFO(2:end,:), "VariableNames", channCharactHFO(1,:));
    table_iesHFO = cell2table(channCharactIES_HFO(2:end,:), "VariableNames", channCharactIES_HFO(1,:));
    table_IES = cell2table(channCharactIES(2:end,:), "VariableNames", channCharactIES(1,:));
    

    
    spreadSheetName = strcat(tablesPath, patName, '.xls');
    delete(spreadSheetName); 
    writetable(table_HFO, spreadSheetName, 'Sheet','HFO');
    writetable(table_iesHFO, spreadSheetName, 'Sheet','iesHFO');
    writetable(table_IES, spreadSheetName, 'Sheet','IES');
end

function metric = getVecMetric(vec)
    metric = mean(vec);
    %metric = median(vec);
end

function rftcPatData = deleteArfeactualChannels(rftcPatData)
    channsList = rftcPatData.channsLabels;
    time = (0:rftcPatData.nrSamples-1)*1/rftcPatData.fs;
    channsStat = {'ChannName', 'Mean Voltage', 'AmplitudeVariance', 'Power'};

    for si = 1:length(channsList)
        channName = channsList{si};
        signal = rftcPatData.signals{si};
        meanVolt = mean(signal);
        ammplVariance = var(signal);
        signalPower = mean(signal.*signal);
        channsStat = cat(1, channsStat, {channName meanVolt ammplVariance signalPower});
    end
    
    varVals = cell2mat(channsStat(2:end, 3));
    powVals = cell2mat(channsStat(2:end, 4));
    nanIdx = find(isnan(varVals))+1;
    
    minBound = prctile(varVals, 25) - 10*iqr(varVals);
    artefIdxLowBound = find(varVals < minBound)+1; 
        
    maxBound = prctile(varVals, 75) + 10*iqr(varVals);
    artefIdxHighBound = find(varVals > maxBound)+1;
    
    artefIdx = [artefIdxLowBound artefIdxHighBound]
    channsStat(artefIdx, 1)
    
    % delete channels of rftcData
    artefIdx = artefIdx-1;
    
    rftcPatData.nrChanns = rftcPatData.nrChanns - length(artefIdx);
    rftcPatData.channsLabels(artefIdx) = [];
    
    rftcPatData.ei(artefIdx) = [];
    rftcPatData.rftc(artefIdx) = [];
    rftcPatData.hfoZone(artefIdx) = [];
    rftcPatData.ia(artefIdx) = [];
    rftcPatData.signals(artefIdx) = [];

end