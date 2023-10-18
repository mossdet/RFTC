clc; clear all; close all;

paths = getFilesPaths();
files  = getFaultyFiles();
files  = getAllFiles();

allPatsDelChannsList = {};
for fileIdx = 1:size(files,1)

    eegFilename = strcat(paths.eegFilesPath, files{fileIdx})
    
    rftcPatData = loadRFTC_Data(paths.workspacePath, eegFilename);

    [origFilepath, patName, ext] = fileparts(rftcPatData.origFilename);
    
    [rftcPatData, delChannsList] = deleteArfeactualChannels(rftcPatData);
    
    newDelEntree = cell(length(delChannsList), 2);
    newDelEntree(:,1) = {patName};
    newDelEntree(:,2) = delChannsList;

    allPatsDelChannsList = cat(1, allPatsDelChannsList, newDelEntree);
    
    %continue;

    nrChanns = rftcPatData.nrChanns;
    

    channCharactHFO = {};
    channCharactIES_HFO = {};
    channCharactIES = {};
    channCharactBP = {};
    
    for chi = 1:nrChanns
        chName = rftcPatData.channsLabels{chi};
        signal = rftcPatData.signals{chi};
        fs = rftcPatData.fs;
        filePath = strcat(paths.workspacePath, 'HFO_Detections_MOSSDET\', patName, '\', chName, '\');        
        loadFN = strcat(filePath, chName, '_hfoCharacterizationShort');
        load(loadFN, 'hfoFeatsShort');
        hdrFeats = hfoFeatsShort(1,:); hdrFeats{7} = 'zone';
        hfoFeatsShort(1,:) = []; 
        
        if chi == 1
           channCharactHFO =  hdrFeats;
           channCharactIES_HFO =  hdrFeats;
           channCharactIES =  hdrFeats;
           channCharactBP =  hdrFeats;
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
        channCharactBP(chi+1,1) = {patFilename};
        
        channCharactHFO(chi+1,2) = {chName};
        channCharactIES_HFO(chi+1,2) = {chName};
        channCharactIES(chi+1,2) = {chName};
        channCharactBP(chi+1,2) = {chName};
        
        channCharactHFO(chi+1,3) = {chi};
        channCharactIES_HFO(chi+1,3) = {chi};
        channCharactIES(chi+1,3) = {chi};
        channCharactBP(chi+1,3) = {chi};
        
        channCharactHFO(chi+1,6) = {1};
        channCharactIES_HFO(chi+1,6) = {2};
        channCharactIES(chi+1,6) = {3};
        channCharactBP(chi+1,6) = {3};
        
        channCharactHFO(chi+1,7) = {0};
        channCharactIES_HFO(chi+1,7) = {0};
        channCharactIES(chi+1,7) = {0};
        channCharactBP(chi+1,7) = {0};
        
        channCharactHFO(chi+1,8) = {rateHFO};
        channCharactIES_HFO(chi+1,8) = {rateIES_HFO};
        channCharactIES(chi+1,8) = {rateIES};
        channCharactBP(chi+1,8) = {rateIES};
       
        % get average feature values across all detections
        for fi = 9:length(hdrFeats)
            featVals =  cell2mat(hfoFeatsShort(:,fi));

            metricFeatHFO = getVecMetric(featVals(selHFO));
            metricFeatIES_HFO = getVecMetric(featVals(selIES_HFO));
            metricFeatIES = getVecMetric(featVals(selIES));
                        
            channCharactHFO(chi+1,fi) = {metricFeatHFO};
            channCharactIES_HFO(chi+1,fi) = {metricFeatIES_HFO};
            channCharactIES(chi+1,fi) = {metricFeatIES};
            channCharactBP(chi+1,fi) = {metricFeatIES};
        end
        
        % Experimental Features
        nrOrigF = length(hdrFeats);
        [expFeatsHdr, expFeatVals] = getBP_Feats(fs, signal);
        nrNewF = length(expFeatsHdr);
        channCharactHFO(1,nrOrigF+1:nrOrigF+nrNewF) = expFeatsHdr(:,:);
        channCharactHFO(end,nrOrigF+1:nrOrigF+nrNewF) = expFeatVals(:,:);
        channCharactIES_HFO(1,nrOrigF+1:nrOrigF+nrNewF) = expFeatsHdr(:,:);
        channCharactIES_HFO(end,nrOrigF+1:nrOrigF+nrNewF) = expFeatVals(:,:);
        channCharactIES(1,nrOrigF+1:nrOrigF+nrNewF) = expFeatsHdr(:,:);
        channCharactIES(end,nrOrigF+1:nrOrigF+nrNewF) = expFeatVals(:,:);
    end
    
    %% Delete NaN channels
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
    
    delIdx_HFO = find(ismember(channCharactHFO(:, channIdx), delChannsList));
    channCharactHFO(delIdx_HFO,:) = [];
    delIdx_iesHFO = find(ismember(channCharactIES_HFO(:, channIdx), delChannsList));
    channCharactIES_HFO(delIdx_HFO,:) = [];
    delIdx_IES = find(ismember(channCharactIES(:, channIdx), delChannsList));
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
    
    tablesPath = strcat(paths.workspacePath, 'CharacterizationTables\FeatureCharacterizationTablesAvg\'); mkdir(tablesPath);
    
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

function [rftcPatData, delChannsList] = deleteArfeactualChannels(rftcPatData)
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
    
    iqrFact = 10;
    minBound = prctile(varVals, 25) - iqrFact*iqr(varVals);
    artefIdxLowBound = find(varVals < minBound)+1; 
        
    maxBound = prctile(varVals, 75) + iqrFact*iqr(varVals);
    artefIdxHighBound = find(varVals > maxBound)+1;
    
    artefIdx = [artefIdxLowBound artefIdxHighBound]
    delChannsList = channsStat(artefIdx, 1)
    
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

function [featsHdr, featVals] = getBP_Feats(fs, signal)
    
    lowFreq = 1; highFreq = 4;
    delta = getBandpassedSignal(fs, signal, lowFreq, highFreq);
    lowFreq = 4; highFreq = 8;
    theta = getBandpassedSignal(fs, signal, lowFreq, highFreq);
    lowFreq = 8; highFreq = 15;
    alpha = getBandpassedSignal(fs, signal, lowFreq, highFreq);
    lowFreq = 15; highFreq = 30;
    beta = getBandpassedSignal(fs, signal, lowFreq, highFreq);
    lowFreq = 30; highFreq = 45;
    gamma = getBandpassedSignal(fs, signal, lowFreq, highFreq);
    lowFreq = 55; highFreq = 90;
    highGamma = getBandpassedSignal(fs, signal, lowFreq, highFreq);
    lowFreq = 80; highFreq = 250;
    ripple = getBandpassedSignal(fs, signal, lowFreq, highFreq);
    lowFreq = 250; highFreq = 500;
    fRipple = getBandpassedSignal(fs, signal, lowFreq, highFreq);
    
    allFreqs = [delta; theta; alpha; beta; gamma; highGamma; ripple; fRipple];
    allFreqsFeat.occRate = zeros(size(allFreqs,1),1);
    allFreqsFeat.amplitude = zeros(size(allFreqs,1),1);
    allFreqsFeat.power = zeros(size(allFreqs,1),1);
    allFreqsFeat.variance = zeros(size(allFreqs,1),1);
    
    epochLength = round(fs*0.1);
    nrMins= length(signal)/fs/60;
    for fi = 1:size(allFreqs,1)
        fSig = allFreqs(fi,:);
        % overlaps or underlaps successive frames in the output matrix by p samples.
        [epochs, rest] = buffer(fSig, epochLength);
        epochAmpl= max(epochs, [], 1) - min(epochs, [], 1);
        epochAvgPow = mean(epochs.*epochs, 1);
        epochVar = var(epochs, 0, 1);

        allFreqsFeat.occRate(fi) = sum(epochAmpl > median(epochAmpl) + 2.5*std(epochAmpl))/nrMins;
        allFreqsFeat.amplitude(fi) = prctile(epochAmpl, 75); mean(epochAmpl);
        allFreqsFeat.power(fi) = prctile(epochAvgPow, 75); mean(epochAvgPow);
        allFreqsFeat.variance(fi) = prctile(epochVar, 75); mean(epochVar);
    end
        
        
    featsHdr = {'deltaOccRate', 'deltaAmpl', 'deltaPow', 'deltaVar', 'thetaOccRate', 'thetaAmpl', 'thetaPow', 'thetaVar', 'alphaOccRate', 'alphaAmpl', 'alphaPow', 'alphaVar',...
    'betaOccRate', 'betaAmpl', 'betaPow', 'betaVar', 'gammaOccRate', 'gammaAmpl', 'gammaPow', 'gammaVar', 'highGammaOccRate', 'highGammaAmpl', 'highGammaPow', 'highGammaVar',...
    'rippleOccRate', 'rippleAmpl', 'ripplePow', 'rippleVar', 'fRippleOccRate', 'fRippleAmpl', 'fRipplePow', 'fRippleVar'};

    featVals = {};
    for fi = 1:size(allFreqs,1)
        featVals = cat(2, featVals, allFreqsFeat.occRate(fi));
        featVals = cat(2, featVals, allFreqsFeat.amplitude(fi));
        featVals = cat(2, featVals, allFreqsFeat.power(fi));
        featVals = cat(2, featVals, allFreqsFeat.variance(fi));
    end
end

%%
function filteredSignal = getBandpassedSignal(fs, signal, lowFreq, highFreq)
    %Filter the signal to calculate features
    order = 128;
    filterDelay = order/2;
    h = fir1(order/2, [lowFreq/(fs/2) highFreq/(fs/2)], 'bandpass'); % 'low' | 'bandpass' | 'high' | 'stop' | 'DC-0' | 'DC-1'
    tempSignal = filter(h, 1, flip(signal));
    tempSignal = filter(h, 1, flip(tempSignal));
    tempSignal(1:filterDelay) = tempSignal(filterDelay+1);
    tempSignal(end-filterDelay:end) = tempSignal(end-filterDelay-1);
    filteredSignal = tempSignal;

end