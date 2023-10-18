clc; clear all; close all;

paths = getFilesPaths();
workspacePath = paths.workspacePath;
eegFilesPath = paths.eegFilesPath;
eiFilesPath = paths.eiFilesPath;
filesMapPath = paths.filesMapPath;
startup;

files  = getAllFiles();
ignoreChannels  = [];%getIgnoreChannels();

allPatsDate = {};
for fileIdx = 1:size(files,1)
    filename = strcat(eegFilesPath, files{fileIdx});
    [~, subjName, ~] = fileparts(filename)

    header = ft_read_header(filename);
    channelLabels = sort(header.label);
    dateStr = strcat(header.orig.year, " ", header.orig.month, " ", header.orig.day);
    allPatsDate = cat(1, allPatsDate, {subjName, dateStr});
    
    continue;
    %Delete data from channels to ignore and clean weird symbols from
    %channel names, call the next 7 lines exactly in this order!
    channsList = header.label;
    [channsList, delChannsIdxs] = deleteChannsToIgnore(channsList, ignoreChannels);
    [unipolarContacts, bipolarMontages] = getBipolarMontageLabels(channsList);

    %%
    saveBipolarHeader(workspacePath, filename, channsList, bipolarMontages, header);
    saveChannelsFile(workspacePath, filename, channsList, bipolarMontages, header);

    %%
    signals = readSignals(filename, delChannsIdxs);
    [bipolarSignals, bipolarChannelsList] = getBipolarMontageSignals(signals, bipolarMontages);
    nrBipolarChanns = length(bipolarMontages);
    %use of parfor is possible here if enough RAM avaibalable
    parfor bchi = 1:nrBipolarChanns
        bipolarSignal = bipolarSignals(bchi, :);
        montageName = bipolarChannelsList{bchi};
        saveBipolarSignal(workspacePath, filename, montageName, bipolarSignal);        
    end
end

function saveBipolarSignal(workspacePath, origFN, montageName, bipolarSignal)
    [origFilepath,filename,ext] = fileparts(origFN);
    filePath = strcat(workspacePath, '\PatientFiles\BipolarSignals\', filename, '\');
    mkdir(filePath);
    saveFN = strcat(filePath, filename, '_', montageName);
    save(saveFN, 'bipolarSignal');
end

function saveBipolarHeader(workspacePath, originalFN, unipChanns, bipolarMontages, header)
    [origFilepath,filename,ext] = fileparts(originalFN);
    hdr.fs = header.Fs;
    hdr.nrChanns = length(bipolarMontages);
    hdr.nrSamples = header.nSamples;
    hdr.unipolarLabels = unipChanns;
    hdr.bipolarLabels = bipolarMontages;
    hdr.origFilename = originalFN;
    hdr.origFilepath = origFilepath;
    hdr.info.surname = header.orig.surname;
    hdr.info.name = header.orig.name;
    hdr.info.day = header.orig.day;
    hdr.info.month = header.orig.month;
    hdr.info.year = header.orig.year;
    hdr.info.Multiplexer = header.orig.Multiplexer;
    hdr.info.unit = header.orig.elec(1).Unit;

    filePath = strcat(workspacePath, 'PatientFiles\BipolarSignals\', filename, '\');
    mkdir(filePath);
    saveFN = strcat(filePath, filename);
    save(saveFN, 'hdr');
end

function saveChannelsFile(workspacePath, originalFN, unipChanns, bipolarMontages, header)
    [origFilepath, filename,ext] = fileparts(originalFN);
    filePath = strcat(workspacePath, 'PatientFiles\BipolarSignals\', filename, '\');
    mkdir(filePath);
    
    EI = zeros(length(bipolarMontages), 1);
    RFTC = zeros(length(bipolarMontages), 1);
    IA = zeros(length(bipolarMontages), 1);
    HFOZ = zeros(length(bipolarMontages), 1);

    bipolarLabels = {bipolarMontages.montageName}'; 
%     for i=1:length(bipolarMontages)
%         bipolarLabels = cat(1, bipolarLabels, bipolarMontages.montageName(i));
%     end

    
    T = table(bipolarLabels, EI, RFTC, IA, HFOZ);
    T = sortrows(T);
    
    channsFilename = strcat(filename, '_BipolarChannels.txt');
    saveFN = strcat(filePath, channsFilename);
    writetable(T, saveFN,'Delimiter','\t') 
end

% make bipolar montages with same letter channles and only(!) consecutive chann.numbers
function [unipolarContacts, bipolarMontages] = getBipolarMontageLabels(allChannsList)
	unipolarContacts = [];
    bipolarMontages = [];
    nrUnipolarChanns = length(allChannsList);
    numbers = ['0','1','2','3','4','5','6','7','8','9'];
    %Get Unipolar contacts
    for uniChIdx = 1:nrUnipolarChanns
		chLabel = allChannsList{uniChIdx};
		contact.contactName = chLabel;
		contact.contactGlobalIdx = uniChIdx;
		foundNrIndices = [];
		for ni = 1:length(numbers)
			strIdx = strfind(chLabel, numbers(ni));
			foundNrIndices = cat(2, foundNrIndices, strIdx);
		end
		firstNumIdx = min(foundNrIndices);
		lastNumIdx = max(foundNrIndices);
		contact.contactNr = str2double(chLabel(firstNumIdx:lastNumIdx));
		contact.electrodeName = chLabel(1:firstNumIdx-1);
        unipolarContacts = cat(1, unipolarContacts, contact);        
    end
    
    %Get Bipolar Montages
    montageNr = 1;
    for upi = 1:size(unipolarContacts, 1)-1
        montage.firstElectrodeName = unipolarContacts(upi).electrodeName;
        montage.firstContactNr = unipolarContacts(upi).contactNr;
        montage.firstContactGlobalIdx = unipolarContacts(upi).contactGlobalIdx;
        
        for supi = 1:size(unipolarContacts, 1)-1
            if upi == supi
                continue;
            end
            montage.secondElectrodeName = unipolarContacts(supi).electrodeName;
            montage.secondContactNr = unipolarContacts(supi).contactNr;
            montage.secondContactGlobalIdx = unipolarContacts(supi).contactGlobalIdx;
            montage.montageName = strcat(unipolarContacts(upi).electrodeName, num2str(unipolarContacts(upi).contactNr),...
                                        '-',...
                                         unipolarContacts(supi).electrodeName, num2str(unipolarContacts(supi).contactNr));
            montage.montageMOSSDET_Nr = montageNr;

            
            if (strcmp(montage.firstElectrodeName, montage.secondElectrodeName) &&...
                (montage.secondContactNr - montage.firstContactNr)==1)
                bipolarMontages = cat(1, bipolarMontages, montage); 
                montageNr = montageNr+1;
            end
        end
    end
end

function [bipolarSignals, bipolarChannelsList] = getBipolarMontageSignals(allSignals, bipolarMontages)
    nrBipolarSignals = size(bipolarMontages,1);
    nrSamples = size(allSignals, 2);
    bipolarSignals = zeros(nrBipolarSignals, nrSamples);
    bipolarChannelsList = cell(nrBipolarSignals, 1);

    for bpi = 1:nrBipolarSignals
        signalA = allSignals(bipolarMontages(bpi).firstContactGlobalIdx, :);
        signalB = allSignals(bipolarMontages(bpi).secondContactGlobalIdx, :);

        bipolarSignals(bpi, :) = (signalA(:)-signalB(:))*-1;
        bipolarChannelsList{bpi,1} = bipolarMontages(bpi).montageName;
    end
end

function signals = readSignals(filename, delChannsIdxs)
    signals = ft_read_data(filename);
    signals(delChannsIdxs, :) = [];
end

function [cleanChannsList, delChannsIdxs] = deleteChannsToIgnore(channsList, ignoreChannels)
    cleanChannsList = {};
    delChannsIdxs = [];
    nrChanns = length(channsList);
    for chi = 1:nrChanns
        chName = channsList{chi};
        ignore = false(1);
        for ich = 1:length(ignoreChannels)
            if strcmp(chName, ignoreChannels{ich})
                ignore = true(1);
                delChannsIdxs = cat(1, delChannsIdxs, chi);
                break;
            end
        end
        if not(ignore)
            cleanChannsList = cat(1, cleanChannsList, chName);
        end   
    end
end

function ignoreChannels  = getIgnoreChannels()

ignoreChannels = {'C3', 'C4', 'Cz', 'F3', 'F4', 'F7', 'F8', 'Fp1', 'FP1', 'Fp2', 'FP2', 'Fz',...
'FZ', 'O1', 'O2', 'P3', 'P4', 'Pz', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'igger',...
'ekg', 'ECG1', 'ECG2', 'EKG', 'EMG', 'emg', 'EOG','EOG_o', 'EOG_u', 'EOG_li', 'EOG_re',...
'MKR1+', 'ecg1', 'ecg2', 'delg1', 'delg2', 'deld1', 'deld2', 'PULS+', 'BEAT+', 'SpO2+', 'MKR2+',...
'ECG1', 'ECG2', 'delG1', 'delG2', 'delD1', 'delD2', 'DelG1', 'DelG2', 'DelD1', 'DelD2', ...
'thor+', 'abdo+', 'xyz+', 'PULS+', 'BEAT+', 'SpO2+', 'MKR2+',...
'MyD1', 'MyD2', 'MyG1', 'MyG2',...
'Deld1', 'Deld2', 'Delg1', 'Delg2', 'my2', 'my1', 'jamG1', 'jamG2', 'jbeD1', 'jbeD2'};
end