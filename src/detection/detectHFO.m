%Example of MATLAB function used to detect HFO, the signal must be at least one minute long to allow the correct normalization of the features used for the detection
%System requirements: Windows 64 bit
function sorted_HFO_IES_Detections = detectHFO(hfoDetectorFolder, hfoSignal, samplingRate, montageName, plotsDir, plotOK)        
    %detect the HFO
    hfoDetections = [];
    hfoSignal = transpose(hfoSignal);
    savedSignalPath = strcat(hfoDetectorFolder, montageName, '.mat');
    save(savedSignalPath,'hfoSignal');
    outputFolder = strcat(hfoDetectorFolder, montageName, '\');
    outputFolder = strrep(outputFolder, '\', '\\'); % this string is passed to a c++ program so the backslash needs to be escaped by anoter backslash
    %outputFolder = strcat('D:\\MATLAB\\Projects\\CCEP_ver3\\MOSSDET_c\\', stimMontageName, '_', responseChannel,'\\');
    
    mossdetVariables.exePath = strcat(hfoDetectorFolder, 'MOSSDET_c.exe');
    mossdetVariables.signalFilePath = savedSignalPath;
    mossdetVariables.decFunctionsPath = hfoDetectorFolder;
    mossdetVariables.outputPath = outputFolder;
    mossdetVariables.startTime = 0;
    mossdetVariables.endTime = 60*60*243*65;    
    mossdetVariables.samplingRate = samplingRate;
    mossdetVariables.eoiType = 'HFO+IES'; %Options are 'HFO+IES' or 'SleepSpindles';
    mossdetVariables.verbose = 0;
    mossdetVariables.saveDetections = 1;

    command = strcat(mossdetVariables.exePath, {' '},...
                     mossdetVariables.signalFilePath, {' '},...
                     mossdetVariables.decFunctionsPath, {' '}, ...
                     mossdetVariables.outputPath, {' '},...
                     num2str(mossdetVariables.startTime), {' '},...
                     num2str(mossdetVariables.endTime), {' '},...
                     num2str(mossdetVariables.samplingRate), {' '},...
                     mossdetVariables.eoiType, {' '},...
                     num2str(mossdetVariables.verbose), {' '},...
                     num2str(mossdetVariables.saveDetections));

    system(command{1})
    
    %read detections from generated text files instead of generating a
    %matlab file, which fails often

    rippleDetections = readDetectionsTextfile(mossdetVariables.outputPath, montageName, 'Ripple');
    frDetections = readDetectionsTextfile(mossdetVariables.outputPath, montageName, 'FastRipple');
    iesDetections = readDetectionsTextfile(mossdetVariables.outputPath, montageName, 'Spike');
    
    hfoDetections = concatenateDetections(rippleDetections, frDetections);
    %cleanIES_Detections = cleanSpikeDetections(iesDetections, 2.5);
    hfoDetections = concatenateDetections(hfoDetections, iesDetections);

    delete(mossdetVariables.signalFilePath);
    rmdir(mossdetVariables.outputPath, 's');
    
    hfoIES_Detections = getIES_CoincidentHFO(hfoDetections);
    sorted_HFO_IES_Detections = sortDetections(hfoIES_Detections);
        
    hfoPlotStruct.signal = hfoSignal;
    hfoPlotStruct.samplingRate = samplingRate;
    hfoPlotStruct.filterOrder = 128;
    hfoPlotStruct.periodToPlot = 1;   % in seconds
    hfoPlotStruct.hfoDetections = sorted_HFO_IES_Detections; 
    hfoPlotStruct.plotsDir = plotsDir;
    hfoPlotStruct.montageName = montageName;
    hfoPlotStruct.lastSample = 600*samplingRate;
    
    eegLength = length(hfoSignal);
    maxDetectIdx = max(sorted_HFO_IES_Detections.endSample);
    if maxDetectIdx > eegLength
        stop = 1;
    end
    
    if plotOK > 0
        plotDetectionsHFO(hfoPlotStruct);
    end
end

