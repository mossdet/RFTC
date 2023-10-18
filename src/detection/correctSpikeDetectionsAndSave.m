function correctSpikeDetectionsAndSave(rftcPatData, workspacePath)

    [origFilepath, patName, ext] = fileparts(rftcPatData.origFilename);
    nrChanns = rftcPatData.nrChanns;

    %Detect and save HFO
    for chi = 1:nrChanns %parfor
        
        channIdxDebug = find(ismember(rftcPatData.channsLabels, 'c1-c2'));
        chi = channIdxDebug;
        
        chName = rftcPatData.channsLabels{chi};
        {patName; chName}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%correct error with exclude channels
        detectsExist = HFO_DetectionsExist(patName, chName, workspacePath);
        if detectsExist
            %Read signal
            signal = detrend(rftcPatData.signals{chi});

            hfoDetections = loadHFO_Detections(workspacePath, patName, chName);
            
            ripplesAndFR = hfoDetections;
            ripplesAndFR.mark(not(hfoDetections.mark == 1 | hfoDetections.mark == 2)) = [];
            ripplesAndFR.startSample(not(hfoDetections.mark == 1 | hfoDetections.mark == 2)) = [];
            ripplesAndFR.endSample(not(hfoDetections.mark == 1 | hfoDetections.mark == 2)) = [];
            
            spikeDetections = hfoDetections;
            spikeDetections.mark(not(hfoDetections.mark == 3)) = [];
            spikeDetections.startSample(not(hfoDetections.mark == 3)) = [];
            spikeDetections.endSample(not(hfoDetections.mark == 3)) = [];

            correctedSpikeDetections =  thresholdSpikeDetections(signal, rftcPatData.fs, spikeDetections, workspacePath, patName, chName);

            MOSSDET_Detections = [correctedSpikeDetections.mark; correctedSpikeDetections.startSample; correctedSpikeDetections.endSample];
            MOSSDET_Detections = cat(2, MOSSDET_Detections, [ripplesAndFR.mark; ripplesAndFR.startSample; ripplesAndFR.endSample]);
            [B,I] = sort(MOSSDET_Detections(2,:));
            MOSSDET_DetectionsSorted = MOSSDET_Detections(:,I); 
            hfoDetectionsCorrected = getIES_CoincidentHFO(MOSSDET_DetectionsSorted);
            
            saveCorrectedHFO_Detections(workspacePath, patName, chName, hfoDetectionsCorrected);

            plotOK = 0;
            if plotOK > 0
                plotsDir = strcat(workspacePath, 'SpikePlots\Corrected\', patName, '\', chName, '\'); mkdir(plotsDir);
               
                hfoPlotStruct.signal = signal;
                hfoPlotStruct.samplingRate = rftcPatData.fs;
                hfoPlotStruct.filterOrder = 128;
                hfoPlotStruct.periodToPlot = 1;   % in seconds
                hfoPlotStruct.hfoDetections = hfoDetectionsCorrected; 
                hfoPlotStruct.plotsDir = plotsDir;
                hfoPlotStruct.montageName = chName;
                hfoPlotStruct.lastSample = 600*rftcPatData.fs;
            
                plotDetectionsHFO(hfoPlotStruct);
            end
            
        else
            stop;
        end
    end

end


function correctedSpikeDetections =  thresholdSpikeDetections(signal, fs, spikeDetections, workspacePath, patName, chName)
    correctedSpikeDetections = [];
    nrDetects = length(spikeDetections.mark);
    nrSamples = length(signal);
    background = [];
    spikesToKeep = false(nrDetects,1);
    
    for i = 1:nrDetects
        ds = spikeDetections.startSample(i);
        de = spikeDetections.endSample(i);
        spikeSignal = signal(ds:de);

        fbs = ds - fs/2;
        fbe = ds - 1;
        sbs = de + 1;
        sbe = de + fs/2;

        if fbs >= 1 && sbe <= nrSamples
            background = [signal(fbs:fbe) signal(sbs:sbe)];
        elseif fbs < 1
            background = signal(sbs:sbe+fs/2);
        else
            background;
        end

        meanAmplBG = mean(background);
        meanPowBG = mean(background.*background);
        meanLL_BG = mean(diff(background));
        
        meanAmplSpike = mean(spikeSignal);
        meanPowSpike = mean(spikeSignal.*spikeSignal);
        meanLL_Spike = mean(diff(spikeSignal));

        %keep = meanAmplSpike > meanAmplBG*2;
        keep = meanPowSpike > meanPowBG*2 && meanLL_Spike > meanLL_BG*2;
        
        spikesToKeep(i) = keep;
        
        plotSpikeToKeep = 0;
        if plotSpikeToKeep > 0
            if fbs >= 1 && sbe <= nrSamples
                plotSignal = signal(fbs:sbe);
                plotTime = double(fbs:sbe)/fs;
                spikePlotTime = double(ds:de)/fs;
                plot(plotTime,plotSignal, '-k'); hold on;
                if keep 
                    plot(spikePlotTime,spikeSignal, '-g'); hold on;
                else
                    plot(spikePlotTime,spikeSignal, '-r'); hold on;
                end
                close();
            end
        end  
    end
    correctedSpikeDetections.mark = spikeDetections.mark(spikesToKeep);    
    correctedSpikeDetections.startSample = spikeDetections.startSample(spikesToKeep);
    correctedSpikeDetections.endSample = spikeDetections.endSample(spikesToKeep);
end

function plotDetections(signal, fs, hfoDetections, plotsDir, chName, plotOK)
    %Plot Detections
    hfoPlotStruct.signal = signal;
    hfoPlotStruct.samplingRate = fs;
    hfoPlotStruct.filterOrder = 128;
    hfoPlotStruct.periodToPlot = 1;   % in seconds
    hfoPlotStruct.hfoDetections = hfoDetections; 
    hfoPlotStruct.plotsDir = plotsDir;
    hfoPlotStruct.montageName = chName;
    hfoPlotStruct.lastSample = 600*fs;

    if plotOK > 0
        plotDetectionsHFO(hfoPlotStruct);
    end
end

function exists = HFO_DetectionsExist(patName, chName, workspacePath)
    filePath = strcat(workspacePath, '\HFO_Detections_MOSSDET\', patName, '\', chName, '\');
    loadFN = strcat(filePath, chName, '_hfoDetections.mat');
    exists = isfile(loadFN);
end

function hfoDetections = getIES_CoincidentHFO(MOSSDET_Detections)
    nrDetections = size(MOSSDET_Detections, 2);
    
    for fdi = 1:nrDetections    %iterate through HFO
        iesCoincidence = 0;
        fdType = MOSSDET_Detections(1, fdi);
        fdStart = MOSSDET_Detections(2, fdi);
        fdEnd = MOSSDET_Detections(3, fdi);
        fdDuration = fdEnd - fdStart;
        
        if(fdType == 3)
            continue;
        end
        
        for sdi = 1:nrDetections    %iterate through IES
            if (fdi == sdi || MOSSDET_Detections(1, sdi) ~= 3)
                continue;
            end
            sdStart = MOSSDET_Detections(2, sdi);
            sdEnd = MOSSDET_Detections(3, sdi);
            sdDuration = sdEnd - sdStart;
            overlapTime = getEventsOverlap(fdStart, fdEnd, sdStart, sdEnd);
            overlapPerc = 100*(overlapTime / fdDuration);
            if (100 * (overlapTime / sdDuration) > overlapPerc)
                overlapPerc = 100 * (overlapTime / sdDuration);
            end
            if overlapPerc > 50.0
                iesCoincidence = 1;
                break;
            end
        end
        
        %     - All Ripples     (1) -> any Ripple
        %     - All FR          (2) -> any FR
        %     - All IES         (3) -> any IES

        %     - IES_Ripples     (4) -> any Ripple coinciding with a IES
        %     - IES_FR          (5) -> any FR coinciding with a IES
        %     - isolRipples     (6) -> any Ripple not coinciding with IES
        %     - isolFR          (7) -> any FR not coinciding with IES

        if fdType == 1
            if iesCoincidence > 0
                MOSSDET_Detections = cat(2, MOSSDET_Detections, [4; fdStart; fdEnd]);
            else
                MOSSDET_Detections = cat(2, MOSSDET_Detections, [6; fdStart; fdEnd]);
            end
        elseif fdType == 2 
            if iesCoincidence > 0
                MOSSDET_Detections = cat(2, MOSSDET_Detections, [5; fdStart; fdEnd]);
            else
                MOSSDET_Detections = cat(2, MOSSDET_Detections, [7; fdStart; fdEnd]);
            end
        end
    end
    
    if isempty(MOSSDET_Detections)
        hfoDetections.mark = [];
        hfoDetections.startSample = [];
        hfoDetections.endSample = []; 
    else
        [~,idx] = sort(MOSSDET_Detections(2,:)); % sort just the second row
        MOSSDET_Detections = MOSSDET_Detections(:,idx);   % sort the whole matrix using the sort indices

        hfoDetections.mark = int64(MOSSDET_Detections(1,:));
        hfoDetections.startSample = MOSSDET_Detections(2,:);
        hfoDetections.endSample = MOSSDET_Detections(3,:);
%         detectionStartTimesLocal = MOSSDET_Detections(2,:);
%         detectionEndTimesLocal = MOSSDET_Detections(3,:);
%         detectionStartSamplesLocal = int64(double(detectionStartTimesLocal).*double(samplingRate));
%         detectionEndSamplesLocal = int64(double(detectionEndTimesLocal).*double(samplingRate));
%         hfoDetections.startSample = detectionStartSamplesLocal;
%         hfoDetections.endSample = detectionEndSamplesLocal; 
    end
end

function overlapTime = getEventsOverlap(feStart, feEnd, seStart, seEnd)
    overlapTime = 0;
    feDuration = feEnd - feStart;
    seDuration = seEnd - seStart;
    
    if feStart <= seStart && feEnd >= seEnd % first fully encompassing second
        overlapTime = seDuration;        
    elseif seStart <= feStart && seEnd >= feEnd % second fully encompassing first
        overlapTime = feDuration;        
    elseif (feStart <= seStart && feEnd >= seStart && feEnd <= seEnd) %last part of first inside second
        overlapTime = feEnd - seStart;
    elseif (seStart <= feStart && seEnd >= feStart && seEnd <= feEnd) %last part of second inside first
        overlapTime = seEnd - feStart;
    else
        overlapTime = 0;
    end
end