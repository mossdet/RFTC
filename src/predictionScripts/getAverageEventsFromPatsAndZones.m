clc; clear all; close all;

paths = getFilesPaths();
workspacePath = paths.workspacePath;
files  = getAllFiles(); % getPoster2023Files getAllFiles getAllFiles getPreFiles getPostFiles
groupTablesPath = strcat(workspacePath, 'CharacterizationTables\GroupCharacterizationTablesMedian\');
zonesNames = {'rftcSite', 'rftcStructure', 'rftcConnected', 'highEI', 'rftcLobe', 'rftcHemisphere'};
%zonesNames = {'rftcConnected'};

% Tables to read
spreadSheetNamePre = strcat(groupTablesPath, 'GroupCharacterization_FlexK_Pre.xls');
spreadSheetNamePost = strcat(groupTablesPath, 'GroupCharacterization_FlexK_Post.xls');
groupTablePre = readtable(spreadSheetNamePre, 'Sheet', 'HFO');
groupTablePost = readtable(spreadSheetNamePost, 'Sheet', 'HFO');

patOutcome = 0;
allPatAverageEOI = {};
ylimVec = [];   
for fileIdx = 1:size(files,1)
    eegFilename = strcat(paths.eegFilesPath, files{fileIdx})
    [origFilepath, patName, ext] = fileparts(eegFilename);

    if strcmp(patName, 'GRE_2017_MESd_Inter_Sleep_2048') || strcmp(patName, 'GRE_2017_MESd_Inter_Sleep_PostRFTC_2048')
        continue;
    end
  

    if sum(ismember(groupTablePre.patName, patName)) > 0
        patSel = ismember(groupTablePre.patName, patName);
        patTable = groupTablePre(patSel, :);
    elseif sum(ismember(groupTablePost.patName, patName)) > 0
        patSel = ismember(groupTablePost.patName, patName);
        patTable = groupTablePost(patSel, :);
    else
        error("Patient not found in Tables")
    end
    
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
        nrIncluded_HFO = 0;
        avg_HFO = [];
        nrIncluded_iesHFO = 0;
        avg_iesHFO = [];
        for chi = 1:length(analysisChannels)
            chName = analysisChannels(chi); 
            rftcChIdx = find(ismember(rftcPatData.channsLabels, chName));
         
            %Read signal
            signal = rftcPatData.signals{rftcChIdx};
            signalBP = getBandpassedSignal(fs, signal, 80, 500);
            time = (0:length(signal)-1)/fs;
    
            %Load HFO
            hfoDetections = [];
            filePath = strcat(workspacePath, 'HFO_Detections_MOSSDET\', patName, '\', chName{1}, '\');
            loadFN = strcat(filePath, chName{1}, '_hfoDetections.mat');
            load(loadFN, 'hfoDetections');

            hfoDetsSel = find(hfoDetections.mark == 1);
            iesHfoDetsSel = find(hfoDetections.mark == 4);
            iesDetsSel = find(hfoDetections.mark == 3);

            nrHFO = length(hfoDetsSel);
            nrIesHFO = length(iesHfoDetsSel);
            nrIES = length(iesDetsSel);
            
            for ii = 1:length(hfoDetsSel)
                hfoIdx = hfoDetsSel(ii);
                hfoStart = hfoDetections.startSample(hfoIdx);
                hfoEnd = hfoDetections.endSample(hfoIdx);
                [M,I] = max(signalBP(hfoStart:hfoEnd));

                centeredHFOStart = hfoStart+I-(fs/2);
                centeredHFOEnd = centeredHFOStart + fs - 1;

                if centeredHFOStart<1 || centeredHFOEnd>length(signal)
                    continue;
                end

                hfoSignalCtrd = signalBP(centeredHFOStart:centeredHFOEnd);
                if fs == 1024
                    hfoSignalCtrd = interp(hfoSignalCtrd,2);
                end
                if isempty(avg_HFO)
                    avg_HFO = hfoSignalCtrd;
                else
                    avg_HFO = avg_HFO + hfoSignalCtrd;
                end
                nrIncluded_HFO = nrIncluded_HFO+1;
            end

            for ii = 1:length(iesHfoDetsSel)
                hfoIdx = iesHfoDetsSel(ii);
                hfoStart = hfoDetections.startSample(hfoIdx);
                hfoEnd = hfoDetections.endSample(hfoIdx);
                [M,I] = max(signalBP(hfoStart:hfoEnd));

                centeredHFOStart = hfoStart+I-(fs/2);
                centeredHFOEnd = centeredHFOStart + fs - 1;

                if centeredHFOStart<1 || centeredHFOEnd>length(signal)
                    continue;
                end

                hfoSignalCtrd = signalBP(centeredHFOStart:centeredHFOEnd);
                if fs == 1024
                    hfoSignalCtrd = interp(hfoSignalCtrd,2);
                end
                if isempty(avg_iesHFO)
                    avg_iesHFO = hfoSignalCtrd;
                else
                    avg_iesHFO = avg_iesHFO + hfoSignalCtrd;
                end
                nrIncluded_iesHFO = nrIncluded_iesHFO+1;
            end 

        end
        dur_min = (rftcPatData.nrSamples/rftcPatData.fs)/60;
        nr_an_channs = length(analysisChannels)
        avg_HFO = avg_HFO/nrIncluded_HFO;
        avg_iesHFO = avg_iesHFO/nrIncluded_iesHFO;
        timeVec = (0:length(avg_HFO)-1)/fs;
        if mod(fileIdx, 2) ~=0
            ylimVec = [min(avg_HFO), max(avg_HFO); min(avg_iesHFO), max(avg_iesHFO)];
        end
        hfo_occ_rate = (nrIncluded_HFO/nr_an_channs)/dur_min;
        ieshfo_occ_rate = (nrIncluded_iesHFO/nr_an_channs)/dur_min;

        %newEntree = {patName, patOutcome, zoneName, 'HFO', avg_HFO, timeVec; patName, patOutcome, zoneName, 'iesHFO', avg_iesHFO, timeVec, hfo_occ_rate, ieshfo_occ_rate};
        newEntree = {patName, patOutcome, zoneName, 'HFO', avg_HFO, timeVec, hfo_occ_rate; patName, patOutcome, zoneName, 'iesHFO', avg_iesHFO, timeVec, ieshfo_occ_rate};

        allPatAverageEOI = cat(1, allPatAverageEOI, newEntree);

        continue;

        subplot(2,3,1)
        plot(timeVec, avg_HFO, 'k'); hold on;
        xlim([min(timeVec), max(timeVec)]);
        ylim(ylimVec(1,:));
        xlabel("Time (s)", 'FontSize', 26);
        ylabel("Amplitude (uV)", 'FontSize', 26);
        title("HFO", 'FontSize', 32);
        ax = gca;
        ax.XAxis.FontSize = 18;
        ax.XAxis.FontWeight = 'bold';
        ax.YAxis.FontSize = 18;
        ax.YAxis.FontWeight = 'bold';

        analSigSel = int32(length(avg_HFO)/2-fs*0.05:length(avg_HFO)/2+fs*0.05);
        analSig = avg_HFO(analSigSel);
        analSigTime = timeVec(analSigSel);
        plot(analSigTime, analSig, '-r')

        maxAmpl = abs(max(analSig)-min(analSig));
        power = sum(analSig.*analSig)/length(analSig);
        annotStr = strcat("Max.Amplitude = ", num2str(maxAmpl), "\muV");
        xPos = min(xlim)+(max(xlim)-min(xlim))*0.7;
        yPos = min(ylim)+(max(ylim)-min(ylim))*0.8;
        text(xPos, yPos, annotStr, 'Horiz','left', 'Vert','top', 'BackgroundColor', 'w', 'Color', 'black', 'FontWeight', 'bold', 'FontSize', 16);

        annotStr = strcat("Power = ", num2str(power), "\muV^2");
        xPos = min(xlim)+(max(xlim)-min(xlim))*0.7;
        yPos = min(ylim)+(max(ylim)-min(ylim))*0.65;
        text(xPos, yPos, annotStr, 'Horiz','left', 'Vert','top', 'BackgroundColor', 'w', 'Color', 'black', 'FontWeight', 'bold', 'FontSize', 16);


        subplot(2,3,4)
        plot(timeVec, avg_iesHFO, 'k'); hold on;
        xlim([min(timeVec), max(timeVec)]);
        ylim(ylimVec(2,:));
        xlabel("Time (s)", 'FontSize', 26);
        ylabel("Amplitude (uV)", 'FontSize', 26);
        title("iesHFO", 'FontSize', 32);
        ax = gca;
        ax.XAxis.FontSize = 18;
        ax.XAxis.FontWeight = 'bold';
        ax.YAxis.FontSize = 18;
        ax.YAxis.FontWeight = 'bold';

        analSigSel = int32(length(avg_iesHFO)/2-fs*0.05:length(avg_iesHFO)/2+fs*0.05);
        analSig = avg_iesHFO(analSigSel);
        analSigTime = timeVec(analSigSel);
        plot(analSigTime, analSig, '-r')

        maxAmpl = abs(max(analSig)-min(analSig));
        power = sum(analSig.*analSig)/length(analSig);
        annotStr = strcat("Max.Amplitude = ", num2str(maxAmpl), "\muV");
        xPos = min(xlim)+(max(xlim)-min(xlim))*0.7;
        yPos = min(ylim)+(max(ylim)-min(ylim))*0.8;
        text(xPos, yPos, annotStr, 'Horiz','left', 'Vert','top', 'BackgroundColor', 'w', 'Color', 'black', 'FontWeight', 'bold', 'FontSize', 16);

        annotStr = strcat("Power = ", num2str(power), "\muV^2");
        xPos = min(xlim)+(max(xlim)-min(xlim))*0.7;
        yPos = min(ylim)+(max(ylim)-min(ylim))*0.65;
        text(xPos, yPos, annotStr, 'Horiz','left', 'Vert','top', 'BackgroundColor', 'w', 'Color', 'black', 'FontWeight', 'bold', 'FontSize', 16);


        patNameFig = strrep(patName, '_', '-');
        titleStr = {patNameFig; strcat("Outcome: ", num2str(patOutcome)); zoneName};

        sgtitle(titleStr)
        set(gcf, 'Position', get(0, 'Screensize'), 'color','w');
        
        posterImagesPath = 'F:\Postdoc_Calgary\Conferences\AES2023\AES_Abstract_2023\Figures\MatlabFigures\';
        posterImagesPath = "F:\ForschungsProjekte\RFTC\RFTC_HFO_Python\Images\Avg_HFO\";
        images_Path = strcat(posterImagesPath, '\');mkdir(images_Path)
        imageFilename = strcat(images_Path, patNameFig, '_', zoneName, '.fig');
        saveas(gcf, imageFilename);
        
        close();
    end
end
save("F:\ForschungsProjekte\RFTC\RFTC_HFO_Python\Images\Avg_HFO\allPatAverageEOI", 'allPatAverageEOI')
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

