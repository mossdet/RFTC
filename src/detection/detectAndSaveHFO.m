function detectAndSaveHFO(rftcPatData, paths)
    workspacePath = paths.workspacePath;
    hfoDetectorFolder = paths.hfoDetectorFolder;

    [origFilepath, patName, ext] = fileparts(rftcPatData.origFilename);
    nrChanns = rftcPatData.nrChanns;
    cleanDetectionWorkspace(paths.hfoDetectorFolder); % clean detection folder in case the previous attempt was interrupted and files remain
    patName
    
    %Detect and save HFO
    for chi = 1:nrChanns %parfor
        chName = rftcPatData.channsLabels{chi};
        {patName, chName}

        %%%%%%%%%%%%%%%%%%%%%%%%%%%correct error with exclude channels
        done = HFO_DetectionsExist(patName, chName, workspacePath)
        %if not(done)
            %Read signal
            signal = detrend(rftcPatData.signals{chi});

            %Detect HFO
            plotsDir = strcat(workspacePath, 'HFO_Plots\', patName, '\', chName, '\'); mkdir(plotsDir);

            plotOK = 1;
            hfoDetections = detectHFO(hfoDetectorFolder, signal, rftcPatData.fs, chName, plotsDir, plotOK);
            %saveHFO_Detections(workspacePath, patName, chName, hfoDetections);
        %end
    end
end

function exists = HFO_DetectionsExist(patName, chName, workspacePath)
    filePath = strcat(workspacePath, '\HFO_Detections_MOSSDET\', patName, '\', chName, '\');
    loadFN = strcat(filePath, chName, '_hfoDetections.mat');
    exists = isfile(loadFN);
end