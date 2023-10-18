function hfoDepuratedDetections = depurateMossdetDetections(hfoDetections, signal, fs, chName, th, plotsDir, plotOk)

selCleanDets = (hfoDetections.maxPower > hfoDetections.avgBackgroundPower + th * hfoDetections.avgBackgroundPowerSD); % threshold the detections

fields = fieldnames(hfoDetections);
nrFields = length(fields);
selRipples = hfoDetections.mark == 1 & selCleanDets;

rippleDetections = [];
for fi = 1:nrFields
    fieldname = fields{fi};
    fieldValues = getfield(hfoDetections, fieldname);
    fieldValues = fieldValues(selRipples);
    rippleDetections = setfield(rippleDetections, fieldname, fieldValues);
end

selFastRipples = hfoDetections.mark == 2 & selCleanDets;
frDetections = [];
for fi = 1:nrFields
    fieldname = fields{fi};
    fieldValues = getfield(hfoDetections, fieldname);
    fieldValues = fieldValues(selFastRipples);
    frDetections = setfield(frDetections, fieldname, fieldValues);
end

selIES = hfoDetections.mark == 3 & selCleanDets;
iesDetections = [];
for fi = 1:nrFields
    fieldname = fields{fi};
    fieldValues = getfield(hfoDetections, fieldname);
    fieldValues = fieldValues(selIES);
    iesDetections = setfield(iesDetections, fieldname, fieldValues);
end

hfoDetections = concatenateDetections(rippleDetections, frDetections);
hfoDetections = concatenateDetections(hfoDetections, iesDetections);

hfoIES_Detections = getIES_CoincidentHFO(hfoDetections);
hfoDepuratedDetections = sortDetections(hfoIES_Detections);

if plotOk
    hfoPlotStruct.signal = signal;
    hfoPlotStruct.samplingRate = fs;
    hfoPlotStruct.filterOrder = 128;
    hfoPlotStruct.periodToPlot = 1;   % in seconds
    hfoPlotStruct.hfoDetections = hfoDepuratedDetections; 
    hfoPlotStruct.plotsDir = plotsDir;
    hfoPlotStruct.montageName = chName;
    hfoPlotStruct.lastSample = 60*fs;

    eegLength = length(signal);
    maxDetectIdx = max(hfoDepuratedDetections.endSample);
    if maxDetectIdx > eegLength
        stop = 1;
    end

    plotDetectionsHFO(hfoPlotStruct);
end
        
end