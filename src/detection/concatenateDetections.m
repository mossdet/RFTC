function hfoDetections = concatenateDetections(detectionsA, detectionsB)

fields = fieldnames(detectionsA);
nrFields = length(fields);
hfoDetections = [];

    for fi = 1:nrFields
        fieldname = fields{fi};
        fieldValuesA = getfield(detectionsA, fieldname);
        fieldValuesB = getfield(detectionsB, fieldname);
        newFieldValues = cat(2, fieldValuesA, fieldValuesB);
        hfoDetections = setfield(hfoDetections, fieldname, newFieldValues);
    end

end