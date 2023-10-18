function perChannDets = loadHFO_DetectionsDelphos(workspacePath, patName)

    filePath = strcat(workspacePath, '\HFO_Detections_Delphos\');    
    loadFN = strcat(filePath, patName, '_delphos_detections.mat');    
    load(loadFN, 'results', 'detection_charac');
    
    [origFilepath, delphosPatName, ext] = fileparts(results.cfg.file);
    
    if not(strcmp(patName, delphosPatName))
        {patName, delphosPatName}
        error('Delphos Detected Filename and EEG Filename are different');
    end
    
    detLabels = {results.markers(:).label};
    detStartSec = {results.markers(:).position};
    detChann = {results.markers(:).channels};
    fs = results.cfg.Fs;
    

    detStartSec = cell2mat(detStartSec);
    detStartSample = detStartSec * fs;
    detChann = cellfun(@(x) x{1}, detChann, 'UniformOutput', false);
    genericDurSamples = ceil(0.05 * fs);
    
    channsList = results.labels;
    
    eventTypes = {'Ripple', 'Fast Ripple', 'Spike'};
    perChannDets = cell(length(channsList),2);
    for chi = 1:length(channsList)
        channDetsLabels = [];
        channDetsSS = [];
        channDetsES = [];

        for ei = 1:length(eventTypes)
            eventName = eventTypes{ei};
            montageName = channsList{chi};
            selChann = ismember(detChann, montageName);
            selEvent = ismember(detLabels, eventName);
            nrEvents = sum(selChann & selEvent);
            if ei == 3
                genericDurSamples = ceil(0.1 * fs);
            end
            
            labels = zeros(1, nrEvents)+ei;
            startSamples = detStartSample(selChann & selEvent);
            endSamples = startSamples + genericDurSamples;

            rate = nrEvents / (results.cfg.duration / 60);            

            channDetsLabels = cat(2, channDetsLabels, labels);
            channDetsSS = cat(2, channDetsSS, startSamples);
            channDetsES = cat(2, channDetsES, endSamples);
        end
        hfoDetections.mark = channDetsLabels;
        hfoDetections.startSample = channDetsSS;
        hfoDetections.endSample = channDetsES; 
        hfoDetections.startTime = channDetsSS / fs;
        hfoDetections.endTime = channDetsES / fs;
        hfoDetections.maxAmpl = channDetsSS * 0;
        hfoDetections.maxPower = channDetsSS * 0;
        hfoDetections.maxSpectralPeak = channDetsSS * 0;
        hfoDetections.avgBackgroundAmplitude = channDetsSS * 0;
        hfoDetections.avgBackgroundPower = channDetsSS * 0;
        hfoDetections.avgBackgroundPowerSD = channDetsSS * 0;

        hfoIES_Detections = getIES_CoincidentHFO(hfoDetections);
        sorted_HFO_IES_Detections = sortDetections(hfoIES_Detections);
    
        perChannDets{chi, 1} = montageName;
        perChannDets{chi, 2} = sorted_HFO_IES_Detections;
    end
end

function plotDetectionsPerChannel(workspacePath, patName, perChannDets, plotsDir)

    filePath = strcat(workspacePath, '\HFO_Detections_Delphos\');    
    loadFN = strcat(filePath, patName, '_delphos_detections.mat');    
    load(loadFN, 'results', 'detection_charac');
    
    cfgDetects =  results.cfg;
    fs = cfgDetects.Fs;    
    channsList = results.labels;
    nrChanns = length(perChannDets); % length(channsList)

    startup;
    filename = cfgDetects.file;
    header = ft_read_header(filename);
    signals = ft_read_data(filename);
    eegChannels = header.labels;
    for chi = 1:length(channsList)
        montageName = channsList{chi};        
        find(ismember(montageName))
        hfoSignal = []; 
        samplingRate = fs;
        hfoDetections = perChannDets(chi);

        hfoPlotStruct.signal = hfoSignal;
        hfoPlotStruct.samplingRate = samplingRate;
        hfoPlotStruct.filterOrder = 128;
        hfoPlotStruct.periodToPlot = 1;   % in seconds
        hfoPlotStruct.hfoDetections = hfoDetections; 
        hfoPlotStruct.plotsDir = plotsDir;
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
