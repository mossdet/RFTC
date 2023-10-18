function hfoDetections = loadHFO_Detections(workspacePath, patName, chName)
    filePath = strcat(workspacePath, '\HFO_Detections_MOSSDET\', patName, '\', chName, '\');
    loadFN = strcat(filePath, chName, '_hfoDetections');
    load(loadFN, 'hfoDetections');
end
