function rftcPatData = loadRFTC_Data(workspacePath, eegFilename)
    [origFilepath,filename,ext] = fileparts(eegFilename);
    filePath = strcat(workspacePath, 'PatientFiles\RFTC_Data\');
    loadFN = strcat(filePath, filename, '_RFTC_Data.mat');
    load(loadFN, 'rftcPatData');
end