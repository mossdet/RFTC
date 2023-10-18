function indexedDetections = getIndexedDetections(detections, idxs)

    fields = fieldnames(detections);
    nrFields = length(fields);
    indexedDetections = [];

    for fi = 1:nrFields
        fieldname = fields{fi};
        allFieldValues = getfield(detections, fieldname);
        indexedFieldValues = allFieldValues(idxs);
        indexedDetections = setfield(indexedDetections, fieldname, indexedFieldValues);
    end
end