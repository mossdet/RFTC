function sortedDetections = sortDetections(detections)

    [~,sortIdx] = sort(detections.startSample); % sort just the second row

    fields = fieldnames(detections);
    nrFields = length(fields);
    sortedDetections = [];

    for fi = 1:nrFields
        fieldname = fields{fi};
        fieldValues = getfield(detections, fieldname);
        sortedFieldValues = fieldValues(sortIdx);
        sortedDetections = setfield(sortedDetections, fieldname, sortedFieldValues);
    end

end