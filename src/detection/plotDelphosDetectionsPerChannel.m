function plotDelphosDetectionsPerChannel(paths, eegFilename, rftcPatData, perChannDets, plotsDir, plotOK)

    detectChanns = perChannDets(:,1);
    nrDetectChanns = length(detectChanns); % length(channsList)
    [origFilepath, patName, ext] = fileparts(eegFilename);

    for chi = 1:nrDetectChanns
        montageName = detectChanns{chi};
        eegChIdx = find(ismember(rftcPatData.channsLabels, montageName));
        if isempty(eegChIdx)
            continue;
        end
        thisPlotsDir = strcat(plotsDir, patName, '\', montageName, '\'); mkdir(thisPlotsDir);

        hfoSignal = rftcPatData.signals{eegChIdx}; 
        samplingRate = rftcPatData.fs;
        hfoDetections = perChannDets{chi, 2};

        hfoPlotStruct.signal = hfoSignal;
        hfoPlotStruct.samplingRate = samplingRate;
        hfoPlotStruct.filterOrder = 128;
        hfoPlotStruct.periodToPlot = 1;   % in seconds
        hfoPlotStruct.hfoDetections = hfoDetections; 
        hfoPlotStruct.plotsDir = thisPlotsDir;
        hfoPlotStruct.montageName = montageName;
        hfoPlotStruct.lastSample = 600*samplingRate;

        eegLength = length(hfoSignal);
        maxDetectIdx = max(hfoDetections.endSample);
        if maxDetectIdx > eegLength
            stop = 1;
        end

        if plotOK > 0
            plotDetectionsHFO(hfoPlotStruct);
        end
    end
end