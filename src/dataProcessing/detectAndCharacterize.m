clc; clear all; close all;

paths = getFilesPaths();

files  = getAllFiles();%getPreFiles, getPostFiles
%files  = getFaultyFiles();
ignoreChannels  = [];

durationMin = {};

for fileIdx = 1:size(files,1)

    eegFilename = strcat(paths.eegFilesPath, files{fileIdx})
    
    %eegFilename = strcat(eegFilesPath, 'GRE_2017_MESd_Inter_Sleep_2048.TRC'); %Debugging
    
    rftcPatData = loadRFTC_Data(paths.workspacePath, eegFilename);
    rftcPatData.nrSamples/rftcPatData.fs/60
    durationMin = cat(1, durationMin, strcat(eegFilename, ", ", num2str(rftcPatData.nrSamples/rftcPatData.fs/60)));
    %saveAvailableMontagesFile(workspacePath, eegFilename, rftcPatData.channsLabels, rftcPatData.ei)
  
    %detectAndSaveHFO(rftcPatData, paths);
    %correctSpikeDetectionsAndSave(rftcPatData, workspacePath);    
    characterizeAndSaveHFO(rftcPatData, paths.workspacePath);
    
    
%     %use of parfor is possible here if enough RAM avaibalable
%     for bchi = 1:nrBipolarChanns
%         bipolarSignal = bipolarSignals(bchi, :);
%         montageName = bipolarChannelsList{bchi};
%         time = (0:length(bipolarSignal)-1)/header.Fs;
%         plot(bipolarSignal);
%         title(montageName)
%     end
end

%%
function saveAvailableMontagesFile(workspacePath, originalFN, bipolarMontages, eiVals)
    [origFilepath, filename,ext] = fileparts(originalFN);
    filePath = strcat(workspacePath, 'PatientFiles\', filename, '\');
    mkdir(filePath);
    
    EI = eiVals;
    RFTC = zeros(length(bipolarMontages), 1);
    IA = zeros(length(bipolarMontages), 1);
    HFOZ = zeros(length(bipolarMontages), 1);
    
    T = table(bipolarMontages, EI, RFTC, IA, HFOZ);
    T = sortrows(T);
    
    channsFilename = strcat(filename, '_AvailableMontages.txt');
    saveFN = strcat(filePath, channsFilename);
    writetable(T, saveFN,'Delimiter','\t') 
end




