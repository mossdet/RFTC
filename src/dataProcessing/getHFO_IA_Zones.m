clc; clear all; close all;

paths = getFilesPaths();

files  = getAllFiles();%getPreFiles, getPostFiles

plotOk = 0;
thresholds = 0:30;

for thi = 1:length(thresholds)
    detTh = thresholds(thi);
    for fileIdx = 1:size(files,1)

        eegFilename = strcat(paths.eegFilesPath, files{fileIdx});
        rftcPatData = loadRFTC_Data(paths.workspacePath, eegFilename);

        [origFilepath, patName, ext] = fileparts(rftcPatData.origFilename);
        patName

        tableFilePath = strcat(paths.workspacePath, 'ChannelCharacterizationTables\MOSSDET\');mkdir(tableFilePath);
        if detTh > 0
            tableFilePath = strcat(paths.workspacePath, 'ChannelCharacterizationTables\MOSSDET_Depurated', '\Th', num2str(detTh),'\');mkdir(tableFilePath);
        end
        tableFN = strcat(tableFilePath, patName, '_', 'ChannelCharacterization_MOSSDET.xls');

        if (isfile(tableFN))% not(isfile(tableFN))
            channelLabels = rftcPatData.channsLabels;
            eiVals = rftcPatData.ei;
            rftcVals = rftcPatData.rftc;
            iaVals = rftcPatData.ia;

            allHFOVals = rftcPatData.hfoZone;
            iesHFOVals = rftcPatData.hfoZone;
            isolHFOVals = rftcPatData.hfoZone;

            allRippleVals = rftcPatData.hfoZone;
            iesRippleVals = rftcPatData.hfoZone;
            isolRippleVals = rftcPatData.hfoZone;

            allFR_Vals = rftcPatData.hfoZone;
            iesFR_Vals = rftcPatData.hfoZone;
            isolFR_Vals = rftcPatData.hfoZone;


            allZonesOccRateTable = table(channelLabels, eiVals, rftcVals, iaVals,...
                allHFOVals, iesHFOVals, isolHFOVals,...
                allRippleVals, iesRippleVals, isolRippleVals,...
                allFR_Vals, iesFR_Vals, isolFR_Vals);

            allZonesPowTable = allZonesOccRateTable;
            allZonesFreqTable = allZonesOccRateTable;

            nrMins = (rftcPatData.nrSamples/rftcPatData.fs)/60;

            rftcTable = getRFTC_Flags(paths, patName);

            for chIdx = 1:rftcPatData.nrChanns %parfor
                chName = allZonesOccRateTable.channelLabels{chIdx};
                done = HFO_DetectionsExist(patName, chName, paths.workspacePath);

                if (done)
                    hfoDetections = loadHFO_Detections(paths.workspacePath, patName, chName);
                    signal = rftcPatData.signals{chIdx};
                    filteredSignalR = getBandpassedSignal(rftcPatData.fs, signal, 80, 250);
                    filteredSignalFR = getBandpassedSignal(rftcPatData.fs, signal, 250, 500);
                    waveletPowSpectrumRipples = getPowerSpectrum(filteredSignalR, rftcPatData.fs);
                    waveletPowSpectrumFR = getPowerSpectrum(filteredSignalFR, rftcPatData.fs);
                    hfoDetections = characterizeBackground(rftcPatData.fs, signal, filteredSignalR, filteredSignalFR, hfoDetections);
                    hfoDetections = characterizeDetectedEvents(rftcPatData.fs, signal, filteredSignalR, filteredSignalFR, waveletPowSpectrumRipples, waveletPowSpectrumFR, hfoDetections);

                    if detTh > 0
                        plotsDir = strcat(paths.workspacePath, 'HFO_Plots_Mossdet_Depurated\', patName, '\', chName, '\Th', num2str(detTh), '\');mkdir(plotsDir);
                        hfoDetections = depurateMossdetDetections(hfoDetections, rftcPatData.signals{chIdx}, rftcPatData.fs, chName, detTh, plotsDir, plotOk);
                    end
                    allZonesOccRateTable = fillZonesTableOccRate(hfoDetections, allZonesOccRateTable, chIdx, nrMins);
                    allZonesPowTable = fillZonesTablePower(hfoDetections, allZonesPowTable, chIdx, nrMins);
                    allZonesFreqTable = fillZonesTableFrequency(hfoDetections, allZonesFreqTable, chIdx, nrMins);

                else
                    stop = 1;
                end

                channFoundIdx = find(ismember(rftcTable.channelLabels, chName));
                allZonesOccRateTable.rftcVals(chIdx) = 0;
                if not(isempty(channFoundIdx))
                    allZonesOccRateTable.rftcVals(chIdx) = rftcTable.rftcVals(channFoundIdx);
                    allZonesPowTable.rftcVals(chIdx) = rftcTable.rftcVals(channFoundIdx);
                    allZonesFreqTable.rftcVals(chIdx) = rftcTable.rftcVals(channFoundIdx);
                else
                    stop = 1;
                end

            end

            writetable(allZonesOccRateTable, tableFN, 'Sheet','OccRate');
            writetable(allZonesPowTable, tableFN, 'Sheet','Power');
            writetable(allZonesFreqTable, tableFN, 'Sheet','Frequency');
        end
    end
end
       
function rftcPatData = loadRFTC_Data(workspacePath, eegFilename)
    [origFilepath,filename,ext] = fileparts(eegFilename);
    filePath = strcat(workspacePath, 'PatientFiles\RFTC_Data\');
    loadFN = strcat(filePath, filename, '_RFTC_Data.mat');
    load(loadFN, 'rftcPatData');
end

function exists = HFO_DetectionsExist(patName, chName, workspacePath)
    filePath = strcat(workspacePath, '\HFO_Detections_MOSSDET\', patName, '\', chName, '\');
    loadFN = strcat(filePath, chName, '_hfoDetections.mat');
    exists = isfile(loadFN);
end

%%
function waveletPowSpectrum = getPowerSpectrum(signal, fs)
        %Get wavelet transform
        [cfs,frq, coi] = cwt(signal,'amor', fs, 'FrequencyLimits',[80 500]);
        absCFS = abs(cfs);
        waveletPowSpectrum.freqBins = frq;
        waveletPowSpectrum.spectrum = absCFS;        
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

%%
function hfoDetections = characterizeBackground(fs, signal, filteredSignalR, filteredSignalFR, hfoDetections)
    bkgrdHalfPeriod = fs*30;
    nrDets = length(hfoDetections.mark);
    eventSamples = [];
    for di = 1:nrDets
        eventSamples = cat(2, eventSamples, hfoDetections.startSample(di) : hfoDetections.endSample(di));
    end
    for di = 1:nrDets
        bss = hfoDetections.startSample(di) - bkgrdHalfPeriod;
        bes = hfoDetections.endSample(di) + bkgrdHalfPeriod;
        if bss < 1
            bss = hfoDetections.endSample(di) + 1;
            bes = hfoDetections.endSample(di) + 2*bkgrdHalfPeriod;
        end
        if bes > length(signal)
            bss = hfoDetections.startSample(di) - 2*bkgrdHalfPeriod;
            bes = hfoDetections.startSample(di) - 1;
        end
        
        bkgrdSamples = bss:bes;
        bkgrdSamples = bkgrdSamples(not(ismember(bkgrdSamples, eventSamples)));
        
        if hfoDetections.mark(di) == 3
            bkgrdSignal = signal(bkgrdSamples);
        elseif (hfoDetections.mark(di) == 1) || (hfoDetections.mark(di) == 4) || (hfoDetections.mark(di) == 6)
            bkgrdSignal = filteredSignalR(bkgrdSamples);
        elseif (hfoDetections.mark(di) == 2) || (hfoDetections.mark(di) == 5) || (hfoDetections.mark(di) == 7)
            bkgrdSignal = filteredSignalFR(bkgrdSamples);
        end
        hfoDetections.avgBackgroundAmplitude(di) = abs(max(bkgrdSignal) - min(bkgrdSignal));
        hfoDetections.avgBackgroundPower(di) = sum(bkgrdSignal.*bkgrdSignal)/length(bkgrdSignal);
        hfoDetections.avgBackgroundPowerSD(di) = std((bkgrdSignal.*bkgrdSignal));
    end
end

%%
function hfoDetections = characterizeDetectedEvents(fs, signal, filteredSignalR, filteredSignalFR, waveletPowSpectrumRipples, waveletPowSpectrumFR, hfoDetections)
    nrDets = length(hfoDetections.mark);
    
    for di = 1:nrDets
        type = hfoDetections.mark(di);
        startSample = hfoDetections.startSample(di);
        endSample = hfoDetections.endSample(di);

        filteredSignal = filteredSignalR;
        waveletPowSpectrum = waveletPowSpectrumRipples;
        if (type == 2) || (type == 5) || (type == 7)
            filteredSignal = filteredSignalFR;
            waveletPowSpectrum = waveletPowSpectrumFR;
        elseif (type == 3)
            filteredSignal = signal;
        end

        detectionSignal = filteredSignal(startSample : endSample);
        segmentLength = length(detectionSignal);
        signalDurationSec = double(segmentLength)/double(fs);

        features.duration = signalDurationSec;
        features.maxAmpl = abs(max(detectionSignal)-min(detectionSignal));
        features.sumAmpl = sum(detectionSignal);
        eoiVariance = var(detectionSignal);
        features.variance = eoiVariance;
        features.BPLL = sum(abs(diff(detectionSignal,1)));
        features.power = sum(detectionSignal.*detectionSignal)/segmentLength;
        features.sumPower = sum(detectionSignal.*detectionSignal);

        %selSampleIdx = zeros(1,length(filteredSignal));selSampleIdx(startSample:endSample) = 1;
        detectionSpectrum = waveletPowSpectrum.spectrum(:, startSample:endSample);
        detectionFreqBins = waveletPowSpectrum.freqBins;
        selectFreqIdx = (detectionFreqBins >=80 & detectionFreqBins <= 250);
        rippleBandPower = sum(detectionSpectrum(selectFreqIdx, :) .* detectionSpectrum(selectFreqIdx, :),'all');
        selectFreqIdx = (detectionFreqBins > 250 & detectionFreqBins <= 500); 
        frBandPower = sum(detectionSpectrum(selectFreqIdx, :) .* detectionSpectrum(selectFreqIdx, :),'all');
        features.PBRatio = frBandPower/rippleBandPower;

        features.spectCentroid = mean(sum(detectionSpectrum(:, :).*detectionFreqBins,1) ./ sum(detectionSpectrum(:, :),1));

        avgMaxPowFreq = 0;
        perSampleMaxFreq = zeros(1,size(detectionSpectrum, 2));
        for i = 1:size(detectionSpectrum, 2)
            [power, freqIdx] = max(detectionSpectrum(:, i));        
            perSampleMaxFreq(i) = detectionFreqBins(freqIdx);
        end
        avgMaxPowFreq = mean(perSampleMaxFreq);
        features.spectPeak = avgMaxPowFreq;
        
        hfoDetections.maxAmpl(di) = features.maxAmpl;
        hfoDetections.maxPower(di) = features.power;
        hfoDetections.maxSpectralPeak(di) = features.spectPeak;
    end
end