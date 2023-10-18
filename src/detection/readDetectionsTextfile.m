function hfoDetections = readDetectionsTextfile(outputFolder, channelsInfo, eventName)

    outputFolder = strrep(outputFolder, '\\', '\');
    detectionOutFilename = strcat(outputFolder, 'MOSSDET_Output\', channelsInfo, '\DetectionFiles\', channelsInfo, '_', eventName, 'DetectionsAndFeatures.txt');
    if not(isfile(detectionOutFilename))
        detectionOutFilename
        stop = 1;
        return;
    end
    %Description	ChannelName	StartTime(s)	EndTime(s)	MaxEventAmplitude	MaxEventPower	MaxEventSpectralPeak (Hz)	AvgBackgroundAmplitude	AvgBackgroundPower	BackgroundStdDev
    [Description, ChannelName, StartSample, EndSample, StartTime, EndTime, MaxEventAmplitude, MaxEventPower, MaxEventSpectralPeak, AvgBackgroundAmplitude, AvgBackgroundPower, AvgBackgroundStdDev] =...   
    textread(detectionOutFilename, '%s\t%s\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f', 'headerlines', 1);

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

    delete(detectionOutFilename);
    mark = 0;
    if strcmp(eventName, 'Ripple')
        mark = 1;
    elseif strcmp(eventName, 'FastRipple')
        mark = 2;
    elseif strcmp(eventName, 'Spike')
        mark = 3;    
    end
    hfoDetections.mark = zeros(1, length(Description)) + mark;
    hfoDetections.startSample = transpose(StartSample);
    hfoDetections.endSample = transpose(EndSample); 
    hfoDetections.startTime = transpose(StartTime);
    hfoDetections.endTime = transpose(EndTime);
    hfoDetections.maxAmpl = transpose(MaxEventAmplitude);
    hfoDetections.maxPower = transpose(MaxEventPower);
    hfoDetections.maxSpectralPeak = transpose(MaxEventSpectralPeak);
    hfoDetections.avgBackgroundAmplitude = transpose(AvgBackgroundAmplitude);
    hfoDetections.avgBackgroundPower = transpose(AvgBackgroundPower);
    hfoDetections.avgBackgroundPowerSD = transpose(AvgBackgroundStdDev);

end