function saveHFO_Detections(workspacePath, patName, chName, hfoDetections)
    filePath = strcat(workspacePath, '\HFO_Detections_MOSSDET\', patName, '\', chName, '\');
    mkdir(filePath);
    saveFN = strcat(filePath, chName, '_hfoDetections');
    save(saveFN, 'hfoDetections');
end
