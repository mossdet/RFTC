clc; clear all; close all;

paths = getFilesPaths();
workspacePath = paths.workspacePath;
files  = getPreFiles();%getAllFiles();%getPreFiles, getPostFiles
groupTablesPath = strcat(workspacePath, 'CharacterizationTables\GroupCharacterizationTablesMedian\');
zonesNames = {'rftcSite', 'rftcStructure', 'rftcConnected', 'highEI', 'rftcLobe', 'rftcHemisphere'};
for fileIdx = 1:size(files,1)
    eegFilename = strcat(paths.eegFilesPath, files{fileIdx})
    [origFilepath, patName, ext] = fileparts(eegFilename);
  
    % Get Zones and Outcome
    spreadSheetNamePre = strcat(groupTablesPath, 'GroupCharacterization_FlexK_Pre.xls');
    groupTablePre = readtable(spreadSheetNamePre, 'Sheet', 'HFO');
    patSel = ismember(groupTablePre.patName, patName);
    patTable = groupTablePre(patSel, :);
    
    patOutcome = patTable.outcome(1);
    patChannels = table2cell(patTable(:,3));

    rftcPatData = loadRFTC_Data(paths.workspacePath, eegFilename);
    fs = rftcPatData.fs;

    durationMin = num2str(rftcPatData.nrSamples/rftcPatData.fs/60);
    nrChanns = rftcPatData.nrChanns;

    for zi = 1:length(zonesNames)
        zoneName = zonesNames{zi};

        if strcmp(zoneName, 'rftcSite')
            % rftcSite
            rftcSiteSel = patTable.rftcSite == 1;
            rftcSiteChannels = rftcPatData.channsLabels(ismember(rftcPatData.channsLabels, patChannels(rftcSiteSel)));
            analysisChannels = rftcSiteChannels;
        elseif strcmp(zoneName, 'rftcStructure')
            % rftcStructure
            rftcStructureSel = (patTable.rftcStruct == 1) & not(patTable.rftcSite == 1);
            rftcStructureChannels  = rftcPatData.channsLabels(ismember(rftcPatData.channsLabels, patChannels(rftcStructureSel)));
            analysisChannels = rftcStructureChannels;
        elseif strcmp(zoneName, 'rftcConnected')
            % rftcConnected
            rftcConnectedSel = (patTable.rftcConn == 1) & not(patTable.rftcSite == 1);
            rftcConnectedChannels  = rftcPatData.channsLabels(ismember(rftcPatData.channsLabels, patChannels(rftcConnectedSel)));
            analysisChannels = rftcConnectedChannels;
        elseif strcmp(zoneName, 'highEI')
            % highEI
            highEISel = (patTable.highEI == 1) & not(patTable.rftcSite == 1);
            highEIChannels  = rftcPatData.channsLabels(ismember(rftcPatData.channsLabels, patChannels(highEISel)));
            analysisChannels = highEIChannels;
        elseif strcmp(zoneName, 'rftcLobe')
            % rftcLobe
            rftcLobeSel = (patTable.rftcLobe == 1) & not(patTable.rftcSite == 1) & not(patTable.rftcStruct == 1);
            rftcLobeChannels  = rftcPatData.channsLabels(ismember(rftcPatData.channsLabels, patChannels(rftcLobeSel)));
            analysisChannels = rftcLobeChannels;
        elseif strcmp(zoneName, 'rftcHemisphere')
            % rftcHemisphere
            rftcHemisphereSel = (patTable.rftcHemis == 1) & not(patTable.rftcSite == 1) & not(patTable.rftcStruct == 1) & not(patTable.rftcLobe == 1);
            rftcHemishpereChannels = rftcPatData.channsLabels(ismember(rftcPatData.channsLabels, patChannels(rftcHemisphereSel)));
            analysisChannels = rftcHemishpereChannels;
        end

        zoneName
        analysisChannels
    
        %Characterize HFO
        for chi = 1:length(analysisChannels)
            chName = analysisChannels(chi); 
            rftcChIdx = find(ismember(rftcPatData.channsLabels, chName));
         
            %Read signal
            signal = rftcPatData.signals{rftcChIdx};
            time = (0:length(signal)-1)/fs;
    
            %Load HFO
            hfoDetections = [];
            filePath = strcat(workspacePath, 'HFO_Detections_MOSSDET\', patName, '\', chName{1}, '\');
            loadFN = strcat(filePath, chName{1}, '_hfoDetections');
            load(loadFN, 'hfoDetections');

            hfoDetsSel = (hfoDetections.mark == 1) || (hfoDetections.mark == 2);
            iesHfoDetsSel = (hfoDetections.mark == 4) || (hfoDetections.mark == 5);
            iesDetsSel = (hfoDetections.mark == 3);

            nrHFO = sum(hfoDetsSel);
            nrIesHFO = sum(iesHfoDetsSel);
            nrIES = sum(iesDetsSel);
            
            if nrHFO < 2
                stop = 1;
            end
        
            filteredSignalR = getBandpassedSignal(fs, signal, 80, 250);
            filteredSignalFR = getBandpassedSignal(fs, signal, 250, 500);
    
            allHFO = cell(nrHFO,1);       
            hfoFeatsShort = cell(nrHFO, 24); % 21 features
            for i = 1:nrHFO %parfor
                %{patName, chName, i}
                hfo.patName = patName;
                hfo.chName = chName;
                type = hfoDetections.mark(i);
                hfo.type = type;
                hfo.chNr = rftcChIdx;
    
                if type == 1% || type == 2
                    %     - All Ripples     (1) -> any Ripple
                    %     - All FR          (2) -> any FR
                    %     - All IES         (3) -> any IES
                    
                    %     - IES_Ripples     (4) -> any Ripple coinciding with a IES
                    %     - IES_FR          (5) -> any FR coinciding with a IES
                    %     - isolRipples     (6) -> any Ripple not coinciding with IES
                    %     - isolFR          (7) -> any FR not coinciding with IES
        
                    hfo.region = 0;%getChannArea(chName);
                    hfo.ss = hfoDetections.startSample(i);
                    hfo.se = hfoDetections.endSample(i);
                    hfoDurMs = (((hfo.se - hfo.ss)+1) / fs)*1000;
                    if hfoDurMs < 20
                        stop = 1;
                    end
    
                    %%
                    filteredSignal = filteredSignalR;
                    hfoName = 'Ripple';
                    if (type == 2) || (type == 5) || (type == 7)
                        filteredSignal = filteredSignalFR;
                        hfoName = 'FR';
                    end
                    hfoSignalRaw = signal(hfo.ss:hfo.se);
                    hfoSignalBP = filteredSignal(hfo.ss:hfo.se);
                    plotSS = hfo.ss - fs/2;
                    plotSE = hfo.ss + fs/2;
                    if plotSS < 1 || plotSE > length(signal)
                        continue;
                    end
    
                    waveletPowSpectrum = getPowerSpectrum(signal, fs, hfo);            
                    hfo = characterizeDetectedEvent(fs, signal, filteredSignalR, filteredSignalFR, waveletPowSpectrum, hfo);
    
                    strTxt = {strcat('Max. Ampl.: ', num2str(hfo.features.maxAmpl)),...
                    strcat('Power: ', num2str(hfo.features.power)),...
                    strcat('Variance: ', num2str(hfo.features.variance))};
        
                    subplot(2,1,1)
                    plot(time(plotSS:plotSE), signal(plotSS:plotSE)); hold on;
                    plot(time(hfo.ss:hfo.se), signal(hfo.ss:hfo.se)); hold on;
                    xlim([time(plotSS), time(plotSE)]);
                    title(chName);
                    ylabel('Voltage (\muV)');
                    xlabel('Time (s)');
                    legend('Raw Signal', hfoName);
                    text(max(xlim), min(ylim)+(max(ylim)-min(ylim))/3,strTxt, 'Horiz','right', 'Vert','top', 'BackgroundColor', 'w', 'Color', 'k', 'FontWeight', 'bold')
    
                    subplot(2,1,2)
                    plot(time(plotSS:plotSE), filteredSignal(plotSS:plotSE)); hold on;
                    plot(time(hfo.ss:hfo.se), filteredSignal(hfo.ss:hfo.se)); hold on;
                    xlim([time(plotSS), time(plotSE)]);
                    title(chName);
                    ylabel('Voltage (\muV)');
                    xlabel('Time (s)');
                    legend('Randpassed Signal', hfoName);
    
                    powerB = sum(filteredSignal(hfo.ss:hfo.se).*filteredSignal(hfo.ss:hfo.se))/length(filteredSignal(hfo.ss:hfo.se));
                    strTxt = {strcat('Max. Ampl.: ', num2str(max(filteredSignal(hfo.ss:hfo.se))-min(filteredSignal(hfo.ss:hfo.se)))),...
                    strcat('Power: ', num2str(powerB)),...
                    strcat('Variance: ', num2str(var(filteredSignal(hfo.ss:hfo.se))))};
                    text(max(xlim), min(ylim)+(max(ylim)-min(ylim))/3,strTxt, 'Horiz','right', 'Vert','top', 'BackgroundColor', 'w', 'Color', 'k', 'FontWeight', 'bold')
    
                    titlePatName = strrep(patName, '_', '-');
                    sgtitle({titlePatName; zoneName})

                    set(gcf, 'Position', get(0, 'Screensize'), 'color','w');
    
                    close();
                end
            end
            
            hfoFeatsShort(1,1:5) = {'patName', 'chName', 'chNr', 'SS', 'SE'};
            hdrFeats = getFeaturesHeader(patName, chName, chi, fs, signal, filteredSignalR, filteredSignalFR, hfoDetections);
            hfoFeatsShort(1,6:end) = hdrFeats(3:end)';
    
        end
    end
end

%%

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
    order = 512;
    filterDelay = order/2;
    h = fir1(order/2, [lowFreq/(fs/2) highFreq/(fs/2)], 'bandpass'); % 'low' | 'bandpass' | 'high' | 'stop' | 'DC-0' | 'DC-1'
    tempSignal = filter(h, 1, flip(signal));
    tempSignal = filter(h, 1, flip(tempSignal));
    tempSignal(1:filterDelay) = tempSignal(filterDelay+1);
    tempSignal(end-filterDelay:end) = tempSignal(end-filterDelay-1);
    filteredSignal = tempSignal;
end

