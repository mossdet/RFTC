function cleanIES = cleanSpikeDetections(iesDetections, th)

selCleanIES = (iesDetections.maxPower > iesDetections.avgBackgroundPower + th * iesDetections.avgBackgroundPowerSD); % threshold the spikes
selCleanIES = (selCleanIES & iesDetections.mark == 3); % make sure we are selecting spikes
    
fields = fieldnames(iesDetections);
nrFields = length(fields);
cleanIES = [];

    for fi = 1:nrFields
        fieldname = fields{fi};
        fieldValues = getfield(iesDetections, fieldname);
        fieldValues = fieldValues(selCleanIES);
        cleanIES = setfield(cleanIES, fieldname, fieldValues);
    end

end