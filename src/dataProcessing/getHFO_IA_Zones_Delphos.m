clc; clear all; close all;

paths = getFilesPaths();

files  = getAllFiles();%getPreFiles, getPostFiles
plotOK = 1;

for fileIdx = 1:size(files,1)

    eegFilename = strcat(paths.eegFilesPath, files{fileIdx});
    [origFilepath, patName, ext] = fileparts(eegFilename);
    done = HFO_DetectionsExistDelphos(patName, paths.workspacePath);
    
    if not(done)
        continue;
    else
        perChannDets = loadHFO_DetectionsDelphos(paths.workspacePath, patName);
        plotsDir = strcat(paths.workspacePath, 'HFO_Plots_Delphos\'); mkdir(plotsDir)

        rftcPatData = loadRFTC_Data(paths.workspacePath, eegFilename);

        %plotDelphosDetectionsPerChannel(paths, eegFilename, rftcPatData, perChannDets, plotsDir, plotOK)
        %continue;

        patName
        
        tableFilePath = strcat(paths.workspacePath, 'ChannelCharacterizationTables\Delphos\');mkdir(tableFilePath);
        tableFN = strcat(tableFilePath, patName, '_', 'ChannelCharacterization_Delphos.xls');

        if not(isfile(tableFN)) %not

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

            detectChanns = perChannDets(:,1);

            for chIdx = 1:rftcPatData.nrChanns %parfor
                chName = allZonesOccRateTable.channelLabels{chIdx};
                detChIdx = find(ismember(detectChanns, chName));
                done = not(isempty(detChIdx));
                if (done)
                    hfoDetections = perChannDets{chIdx, 2};

                    signal = rftcPatData.signals{chIdx};
                    filteredSignalR = getBandpassedSignal(rftcPatData.fs, signal, 80, 250);
                    filteredSignalFR = getBandpassedSignal(rftcPatData.fs, signal, 250, 500);
                    waveletPowSpectrumRipples = getPowerSpectrum(filteredSignalR, rftcPatData.fs);
                    waveletPowSpectrumFR = getPowerSpectrum(filteredSignalFR, rftcPatData.fs);
                    hfoDetections = characterizeDetectedEvents(rftcPatData.fs, filteredSignalR, filteredSignalFR, waveletPowSpectrumRipples, waveletPowSpectrumFR, hfoDetections);

                    allZonesOccRateTable = fillZonesTableOccRate(hfoDetections, allZonesOccRateTable, chIdx, nrMins);
                    allZonesPowTable = fillZonesTablePower(hfoDetections, allZonesPowTable, chIdx, nrMins);
                    allZonesFreqTable = fillZonesTableFrequency(hfoDetections, allZonesFreqTable, chIdx, nrMins);
                else
                    hfoDetections.mark = [];
                    hfoDetections.startSample = [];
                    hfoDetections.endSample = []; 
                    hfoDetections.startTime = [];
                    hfoDetections.endTime = [];
                    hfoDetections.maxAmpl = [];
                    hfoDetections.maxPower = [];
                    hfoDetections.maxSpectralPeak = [];
                    hfoDetections.avgBackgroundAmplitude = [];
                    hfoDetections.avgBackgroundPower = [];
                    hfoDetections.avgBackgroundPowerSD = [];
                    allZonesOccRateTable = fillZonesTableOccRate(hfoDetections, allZonesOccRateTable, chIdx, nrMins);
                    allZonesPowTable = fillZonesTablePower(hfoDetections, allZonesPowTable, chIdx, nrMins);
                    allZonesFreqTable = fillZonesTableFrequency(hfoDetections, allZonesFreqTable, chIdx, nrMins);
                end

                channFoundIdx = find(ismember(rftcTable.channelLabels, chName));
                allZonesOccRateTable.rftcVals(chIdx) = 0;
                if not(isempty(channFoundIdx))
                    allZonesOccRateTable.rftcVals(chIdx) = rftcTable.rftcVals(channFoundIdx);
                else
                    stop = 1;
                end

            end
            allZonesPowTable.rftcVals = allZonesOccRateTable.rftcVals;
            allZonesFreqTable.rftcVals = allZonesOccRateTable.rftcVals;

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

function exists = HFO_DetectionsExistDelphos(patName, workspacePath)
    filePath = strcat(workspacePath, '\HFO_Detections_Delphos\');    
    loadFN = strcat(filePath, patName, '_delphos_detections.mat');
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
function hfoDetections = characterizeDetectedEvents(fs, filteredSignalR, filteredSignalFR, waveletPowSpectrumRipples, waveletPowSpectrumFR, hfoDetections)
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