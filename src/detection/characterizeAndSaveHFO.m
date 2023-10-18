function characterizeAndSaveHFO(rftcPatData, workspacePath)
    
    [origFilepath, patName, ext] = fileparts(rftcPatData.origFilename);
    nrChanns = rftcPatData.nrChanns;
    
    %Characterize HFO
    for chi = 1:nrChanns
        chName = rftcPatData.channsLabels{chi};

        filePath = strcat(workspacePath, 'HFO_Detections_MOSSDET\', patName, '\', chName, '\');
        saveFN = strcat(filePath, chName, '_characterizedHFO.mat');
        
        filePath = strcat(workspacePath, 'HFO_Detections_MOSSDET\', patName, '\', chName, '\'); mkdir(filePath);
        saveFN = strcat(filePath, chName, '_hfoCharacterizationShort.mat');

        if isfile(saveFN)
            continue;
        end
        
        %Read signal
        signal = rftcPatData.signals{chi};
        fs = rftcPatData.fs;

        %Load HFO
        hfoDetections = [];
        filePath = strcat(workspacePath, 'HFO_Detections_MOSSDET\', patName, '\', chName, '\');
        loadFN = strcat(filePath, chName, '_hfoDetections');
        load(loadFN, 'hfoDetections');
        nrHFO = length(hfoDetections.mark);
        
        if nrHFO < 2
            stop = 1;
        end
    
        filteredSignalR = getBandpassedSignal(fs, signal, 80, 250);
        filteredSignalFR = getBandpassedSignal(fs, signal, 250, 500);
        
        allHFO = cell(nrHFO,1);
       
%         header = {'chName', 'chNr', 'type', 'zone', 'SS', 'SE', 'BP_Power', 'wvltRipplePower', 'wvltFR_Power', 'Duration', 'SpctrlPeak', 'SpctrlCentroid'};
        hfoFeatsShort = cell(nrHFO, 24); % 21 features
        parfor i = 1:nrHFO %parfor
            %{patName, chName, i}
            hfo.patName = patName;
            hfo.chName = chName;
            hfo.chNr = chi;
            hfo.type = hfoDetections.mark(i);
            hfo.region = 0;%getChannArea(chName);
            hfo.ss = hfoDetections.startSample(i);
            hfo.se = hfoDetections.endSample(i);
            hfoDurMs = (((hfo.se - hfo.ss)+1) / rftcPatData.fs)*1000;
            if hfoDurMs < 20
                stop = 1;
            end
            
            waveletPowSpectrum = getPowerSpectrum(signal, fs, hfo);            
            hfo = characterizeDetectedEvent(fs, signal, filteredSignalR, filteredSignalFR, waveletPowSpectrum, hfo);

            %get spindle coincidence
%             hfo.spindle = false(1);%getSpindleCoincidence(spindleModels, hfo);
%             hfo.spindlePlus = false(1);
%             hfo.iSpindle = false(1);%getSpindleCoincidence(spindles.invasive, hfo);
%             hfo.iSpindlePlus = false(1);
            allHFO{i} = hfo;
                
            hdrB = fieldnames(hfo.features);
            entree = cell(1, length(hdrB)-2+5);
            entree(1,1:5) = {patName, chName, hfo.chNr, hfo.ss, hfo.se};            
            for fi = 1:length(hdrB)
                if fi > 2
                    entree{1,fi-2+5} = hfo.features.(hdrB{fi});
                end
            end
            hfoFeatsShort(i+1,:) = entree;
                
            %hfoFeatsShort(i+1,:) = {chName, hfo.chNr, hfo.type, hfo.ss, hfo.se, hfo.features.power, hfo.features.wvltRipplePower, hfo.features.wvltFastRiplePower, hfo.features.duration, hfo.features.spectPeak, hfo.features.spectCentroid};
        end
        
        hfoFeatsShort(1,1:5) = {'patName', 'chName', 'chNr', 'SS', 'SE'};
        hdrFeats = getFeaturesHeader(patName, chName, chi, fs, signal, filteredSignalR, filteredSignalFR, hfoDetections);
        hfoFeatsShort(1,6:end) = hdrFeats(3:end)';
                
%         filePath = strcat(workspacePath, 'HFO_Detections_MOSSDET\', patName, '\', chName, '\'); mkdir(filePath);
%         saveFN = strcat(filePath, chName, '_characterizedHFO');
%         save(saveFN, 'allHFO');      
        
        filePath = strcat(workspacePath, 'HFO_Detections_MOSSDET\', patName, '\', chName, '\'); mkdir(filePath);
        saveFN = strcat(filePath, chName, '_hfoCharacterizationShort.mat');
        save(saveFN, 'hfoFeatsShort');  
    end
    %END Characterize HFO

end


function hdrFeats = getFeaturesHeader(patName, chName, chi, fs, signal, filteredSignalR, filteredSignalFR, hfoDetections)
    hfo.patName = patName;
    hfo.chName = chName;
    hfo.chNr = chi;
    i = 1;
    hfo.type = hfoDetections.mark(i);
    hfo.region = 0;
    hfo.ss = hfoDetections.startSample(i);
    hfo.se = hfoDetections.endSample(i);
    waveletPowSpectrum = getPowerSpectrum(signal, fs, hfo);
    hfo = characterizeDetectedEvent(fs, signal, filteredSignalR, filteredSignalFR, waveletPowSpectrum, hfo);
    hdrFeats = fieldnames(hfo.features);
end

function hfo = characterizeDetectedEvent(fs, signal, filteredSignalR, filteredSignalFR, waveletPowSpectrum, hfo)
    type = hfo.type;
    startSample = hfo.ss;
    endSample = hfo.se;
    
    filteredSignal = filteredSignalR;
    if (type == 2) || (type == 5) || (type == 7)
        filteredSignal = filteredSignalFR;
    elseif (type == 3)
        filteredSignal = signal;
    end
    
    detectionSignal = filteredSignal(startSample : endSample);
    segmentLength = length(detectionSignal);
    signalDurationSec = double(segmentLength)/double(fs);

    features.patNr = 0;
    features.chNr = 0;
    features.type = hfo.type;
    features.area = 0;
    features.rate = 0;
    features.duration = signalDurationSec;
    features.maxAmpl = abs(max(detectionSignal)-min(detectionSignal));
    features.sumAmpl = sum(detectionSignal);
    eoiVariance = var(detectionSignal);
    features.variance = eoiVariance;
    features.BPLL = sum(abs(diff(detectionSignal,1)));
    features.power = sum(detectionSignal.*detectionSignal)/segmentLength;
    features.sumPower = sum(detectionSignal.*detectionSignal);

    nrPeaks = 0;
    peakBufferLength = 3;
    if segmentLength > peakBufferLength+1
        for i = peakBufferLength+1:segmentLength-peakBufferLength
            peakPreTestNeg = sum((detectionSignal(i-peakBufferLength:i-1)) < (detectionSignal(i)));
            peakPostTestNeg = sum((detectionSignal(i+1:i+peakBufferLength)) < (detectionSignal(i)));
            peakPreTestPos = sum((detectionSignal(i-peakBufferLength:i-1)) > (detectionSignal(i)));
            peakPostTestPos = sum((detectionSignal(i+1:i+peakBufferLength)) > (detectionSignal(i)));
            if(peakPreTestNeg >= peakBufferLength && peakPostTestNeg >= peakBufferLength) || (peakPreTestPos >= peakBufferLength && peakPostTestPos >= peakBufferLength)
                nrPeaks = nrPeaks+1;
            end
        end
    end
    peaks = (double(nrPeaks) / signalDurationSec);
    
    eoiMobility = (sqrt(var(diff(detectionSignal,1))/eoiVariance));
    eoiComplexity = sqrt( var(diff(detectionSignal,2)) / var(diff(detectionSignal,1)) ) / eoiMobility;
    features.mobility = eoiMobility;
    features.complexity = eoiComplexity;
    
    nrZC = 0;
    zeroRef = mean(detectionSignal);
    zcBuffer = 2;
    for i = 1:segmentLength
        preS = i-zcBuffer;
        preE = i-1;
        postS = i;
        postE = i+zcBuffer-1;
        if preS >= 1 && preE <= segmentLength && postS >= 1 && postE <= segmentLength
            zcPreTestDown = sum(detectionSignal(preS:preE) < zeroRef);
            zcPostTestDown = sum(detectionSignal(postS:postE) > zeroRef);
            zcPreTestUp = sum(detectionSignal(preS:preE) > zeroRef);
            zcPostTestUp = sum(detectionSignal(postS:postE) < zeroRef);
            if(zcPreTestDown >= zcBuffer && zcPostTestDown >= zcBuffer) || (zcPreTestUp >= zcBuffer && zcPostTestUp >= zcBuffer)
                nrZC = nrZC+1;
            end
        end
    end
    
    zeroCrossingsPerEOI = (double(nrZC) / signalDurationSec);
        
    selectIdx = (waveletPowSpectrum.freqBins >=80 & waveletPowSpectrum.freqBins <= 250); 
    rippleBandPower = max(mean(waveletPowSpectrum.spectrum(selectIdx, :), 1));
    selectIdx = (waveletPowSpectrum.freqBins > 250 & waveletPowSpectrum.freqBins <= 500); 
    frBandPower = max(mean(waveletPowSpectrum.spectrum(selectIdx, :), 1));
    features.PBRatio = frBandPower/rippleBandPower;
    features.wvltRipplePower = rippleBandPower;
    features.wvltFastRiplePower = frBandPower;

    features.peaks = peaks;
    features.zeroCrossPerEOI = zeroCrossingsPerEOI;
    features.spectCentroid = mean(sum(waveletPowSpectrum.spectrum(:, :).*waveletPowSpectrum.freqBins,1) ./ sum(waveletPowSpectrum.spectrum(:, :),1));
    
    avgMaxPowFreq = 0;
    for i = 1:size(waveletPowSpectrum.spectrum, 2)
        [power, freqIdx] = max(waveletPowSpectrum.spectrum(:, i));        
		avgMaxPowFreq = avgMaxPowFreq+waveletPowSpectrum.freqBins(freqIdx);
    end
	avgMaxPowFreq = avgMaxPowFreq/size(waveletPowSpectrum.spectrum, 2);
    features.spectPeak = avgMaxPowFreq;
    
    features.maxAmpl = abs(max(detectionSignal)-min(detectionSignal));    
    hfo.features = features;
end


%%
function waveletPowSpectrum = getPowerSpectrum(signal, fs, hfo)
        wdw = fs*1;
        ss = hfo.ss-wdw;
        se = hfo.se+wdw;
        sl = hfo.se-hfo.ss+1;
        if ss <= 0
           ss = 1;
           wdw = hfo.ss;
        end
        if se > length(signal)
           se = length(signal); 
        end

        hfoSignal = signal(ss:se);
        %Get wavelet transform
        [cfs,frq, coi] = cwt(hfoSignal,'amor', fs, 'FrequencyLimits',[80 500]);
        absCFS = abs(cfs);
                    absCFS = absCFS(:, wdw:wdw+sl-1);
        waveletPowSpectrum.freqBins = frq;
        waveletPowSpectrum.spectrum = absCFS;
        
        if size(absCFS,2) ~= sl
            stop = 1;
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