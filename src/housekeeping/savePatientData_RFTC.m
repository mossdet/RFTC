clc; clear all; close all;

paths = getFilesPaths();
workspacePath = paths.workspacePath;
eegFilesPath = paths.eegFilesPath;
eiFilesPath = paths.eiFilesPath;
filesMapPath = paths.filesMapPath;

files  = getAllFiles();%getPreFiles, getPostFiles
ignoreChannels  = [];

for fileIdx = 1:size(files,1)

    eegFilename = strcat(eegFilesPath, files{fileIdx});
    
    rftcExists = rftcDataExists(workspacePath, eegFilename);
    
    if not(rftcExists)
        [eiChannels, eiVals] = readEI(files{fileIdx});

        bipSignalsCell = [];
        for i = 1:length(eiChannels)
            montageName = eiChannels{i};
            eiVal = eiVals{i};
            bipolarSignal = loadBipolarSignal(workspacePath, eegFilename, montageName);
            if isempty(bipolarSignal)
                eegFilename
                montageName
                'Montage not found!'
                stop = 1
            else
                bipSignalsCell = cat(1, bipSignalsCell, {montageName, eiVal, bipolarSignal});
            end
        end    

        hdr = readBipolarHeader(workspacePath, eegFilename);

        rftcPatData.fs = hdr.fs;
        rftcPatData.nrChanns = length(bipSignalsCell);
        rftcPatData.nrSamples = hdr.nrSamples;
        rftcPatData.channsLabels = bipSignalsCell(:,1);
        rftcPatData.origFilename = hdr.origFilename;
        rftcPatData.origFilepath = hdr.origFilepath;
        rftcPatData.info = hdr.info;

        rftcPatData.ei = bipSignalsCell(:,2);   
        rftcPatData.rftc = zeros(length(rftcPatData.ei),1);
        rftcPatData.hfoZone = zeros(length(rftcPatData.ei),1);
        rftcPatData.ia = zeros(length(rftcPatData.ei),1);

        rftcPatData.signals = bipSignalsCell(:,3);

        saveRFTC_Data(workspacePath, eegFilename, rftcPatData);
    else
        ; % rftcPatData = loadRFTC_Data(workspacePath, eegFilename);
    end
end

function saveRFTC_Data(workspacePath, eegFilename, rftcPatData)
    [origFilepath,filename,ext] = fileparts(eegFilename);
    filePath = strcat(workspacePath, 'PatientFiles\RFTC_Data\');
    mkdir(filePath);
    saveFN = strcat(filePath, filename, '_RFTC_Data.mat');
    save(saveFN, 'rftcPatData', '-v7.3');
end

function rftcPatData = loadRFTC_Data(workspacePath, eegFilename)
    [origFilepath,filename,ext] = fileparts(eegFilename);
    filePath = strcat(workspacePath, 'PatientFiles\RFTC_Data\');
    loadFN = strcat(filePath, filename, '_RFTC_Data.mat');
    load(loadFN, 'rftcPatData');
end

function exists = rftcDataExists(workspacePath, eegFilename)
    [origFilepath,filename,ext] = fileparts(eegFilename);
    filePath = strcat(workspacePath, 'PatientFiles\RFTC_Data\');
    loadFN = strcat(filePath, filename, '_RFTC_Data.mat');
    exists = isfile(loadFN);
end

function bipolarSignal = loadBipolarSignal(workspacePath, eegFN, montageName)
    bipolarSignal = [];
    [origFilepath,filename,ext] = fileparts(eegFN);
    filePath = strcat(workspacePath, 'PatientFiles\BipolarSignals\', filename, '\');

    loadFN = strcat(filePath, filename, '_', montageName, '.mat');
    isfile(loadFN)

    if isfile(loadFN)
        load(loadFN, 'bipolarSignal');
        %bipolarSignal = 1;
    end
end

function hdr = readBipolarHeader(workspacePath, eegFilename)
    [origFilepath, filename,ext] = fileparts(eegFilename);
    filePath = strcat(workspacePath, 'PatientFiles\BipolarSignals\', filename, '\');
    loadFN = strcat(filePath, filename, '.mat');
    hdr = load(loadFN, 'hdr');
    hdr = hdr.hdr;
end