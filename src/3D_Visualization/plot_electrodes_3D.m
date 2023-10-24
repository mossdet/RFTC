clc; clear all; close all;

[path,~,~] = fileparts(mfilename('fullpath'));
cutIdx = strfind(path, '\');
workspacePath = path(1:cutIdx(end));
cd(workspacePath)
addpath(genpath(workspacePath));
%startup

paths = getFilesPaths();
workspacePath = paths.workspacePath;
eegFilesPath = paths.eegFilesPath;
eiFilesPath = paths.eiFilesPath;
filesMapPath = paths.filesMapPath;

files  = getPreFiles();%getPostFiles();
ignoreChannels  = [];%getIgnoreChannels();
paths.new_eeg_path = 'F:\ForschungsProjekte\RFTC\RFTC_HFO_Python\Data\';

for fileIdx = 1:size(files,1)
    fileIdx
    filename = strcat(eegFilesPath, files{fileIdx})
    get_3d_coordinates(paths, filename)
    %convert_file(paths, filename)
end

function get_3d_coordinates(paths, filename)
    anatLocalizationGrenoblePath = paths.anatLocalizationGrenoblePath;

    hdr = ft_read_header(filename);
    [~, patName, ~] = fileparts(filename);

    listing = dir(anatLocalizationGrenoblePath);
    nrContents = length(listing);
    for ci = 1:nrContents
        if strfind(listing(ci).name, patName) > 0
            anatLocTableJJ = readtable(strcat(anatLocalizationGrenoblePath, listing(ci).name));
            break;
        end
    end

    for li = 1:length(anatLocTableJJ.contact)
        anatLocTableJJ.contact{li} = getCorrectChannelNames(anatLocTableJJ.contact{li});
    end
    %mni_coord_str = anatLocTableJJ.MNI;
    mni_coord_str = anatLocTableJJ.T1preScannerBased;
    mni_coord_str = strrep(mni_coord_str, '[','');
    mni_coord_str = strrep(mni_coord_str, ']','');
    mni_coord_str = strrep(mni_coord_str, ' ','');
    mni_coord = zeros(length(mni_coord_str),3);
    for chi = 1:length(mni_coord_str)
        coord_str = mni_coord_str{chi,:};
        sep_idx = strfind(coord_str, ',');
        mni_coord(chi,1) = str2double(coord_str(1:sep_idx(1)-1));
        mni_coord(chi,2) = str2double(coord_str(sep_idx(1)+1:sep_idx(2)-1));
        mni_coord(chi,3) = str2double(coord_str(sep_idx(2)+1:end));
    end 

    mni_contacts = anatLocTableJJ.contact;
    coord_cell = {};
    for ei = 1:length(hdr.orig.elec)
        contact_name = hdr.orig.elec(ei).Name;
        if ~strcmp(contact_name,'NaN')
            coord_idx = ismember(mni_contacts, contact_name);
            if sum(coord_idx)>0
                if sum(coord_idx)>1
                    error("More than one contact with this name was found")
                else
                    entree = {contact_name, mni_coord(coord_idx,1), mni_coord(coord_idx,2), mni_coord(coord_idx,3)};
                    coord_cell = cat(1, coord_cell, entree);
                end
            end
        end
    end
    coord_table = cell2table(coord_cell,"VariableNames", ["ChName" "XPos" "YPos" "ZPos"]);

    coord_table_fn = strcat(paths.new_eeg_path, patName, '.csv');
    writetable(coord_table,coord_table_fn,'Delimiter',',') 
end

function convert_file(paths, filename)

    anatLocalizationGrenoblePath = paths.anatLocalizationGrenoblePath;

    hdr = ft_read_header(filename);
    [~, patName, ~] = fileparts(filename);

    listing = dir(anatLocalizationGrenoblePath);
    nrContents = length(listing);
    for ci = 1:nrContents
        if strfind(listing(ci).name, patName) > 0
            anatLocTableJJ = readtable(strcat(anatLocalizationGrenoblePath, listing(ci).name));
            break;
        end
    end
    for li = 1:length(anatLocTableJJ.contact)
        anatLocTableJJ.contact{li} = getCorrectChannelNames(anatLocTableJJ.contact{li});
    end
    mni_coord_str = anatLocTableJJ.MNI;
    mni_coord_str = strrep(mni_coord_str, '[','');
    mni_coord_str = strrep(mni_coord_str, ']','');
    mni_coord_str = strrep(mni_coord_str, ' ','');
    mni_coord = zeros(length(mni_coord_str),3);
    for chi = 1:length(mni_coord_str)
        coord_str = mni_coord_str{chi,:};
        sep_idx = strfind(coord_str, ',');
        mni_coord(chi,1) = str2double(coord_str(1:sep_idx(1)-1));
        mni_coord(chi,2) = str2double(coord_str(sep_idx(1)+1:sep_idx(2)-1));
        mni_coord(chi,3) = str2double(coord_str(sep_idx(2)+1:end));
    end

    mni_contacts = anatLocTableJJ.contact;
    save_ch_vec = zeros(1, length(hdr.orig.elec));
    for ei = 1:length(hdr.orig.elec)
        contact_name = hdr.orig.elec(ei).Name;
        type_code = 0;
        if strcmp(contact_name,'NaN')
            contact_name = hdr.orig.elec(ei).Name;
            type_code = 1;
        end
        coord_idx = ismember(mni_contacts, contact_name);
        if sum(coord_idx)>0
            if sum(coord_idx)>1
                error("More than one contact with this name was found")
            else
                hdr.orig.elec(ei).XPos = mni_coord(coord_idx,1);
                hdr.orig.elec(ei).YPos = mni_coord(coord_idx,2);
                hdr.orig.elec(ei).ZPos = mni_coord(coord_idx,3);
            end
            hdr.orig.elec(ei).Type = type_code;
            hdr.chantype = type_code;
            save_ch_vec(ei) = true();
        else
            hdr.orig.elec(ei).Type = 0;
            hdr.chantype = 0;
        end
        hdr.label{ei} = contact_name;
        hdr.orig.elec(ei).Name = contact_name;
    end
    data = ft_read_data(filename);

    new_eeg_fn = strcat(paths.new_eeg_path, patName, '.edf');
    ft_write_data(new_eeg_fn, data, 'header', hdr,'dataformat','edf')

    new_eeg_fn = strcat(paths.new_eeg_path, patName, '.vhdr');
    ft_write_data(new_eeg_fn, data, 'header', hdr,'dataformat', 'brainvision_eeg');
end

function chNameCorr = getCorrectChannelNames(chName)
    chNameCorr = strrep(chName, ' ', '');
    chNameCorr = strrep(chNameCorr, 'p', '''');
    chNameCorr = lower(chNameCorr);
    if not(isempty(chNameCorr))
        if contains(chNameCorr, '-')
            chA = chNameCorr(1:strfind(chNameCorr, '-')-1);
            chB = chNameCorr(strfind(chNameCorr, '-')+1:end);
            chNameCorr = strcat(getCorrectChannelNumber(chA), '-', getCorrectChannelNumber(chB));
        else
            chNameCorr = getCorrectChannelNumber(chNameCorr);
        end
    end
end

function chNameCorr = getCorrectChannelNumber(chName)
    numbers = ['0','1','2','3','4','5','6','7','8','9'];
    foundNrIndices = [];
    for ni = 1:length(numbers)
        strIdx = strfind(chName, numbers(ni));
        foundNrIndices = cat(2, foundNrIndices, strIdx);
    end
    firstNumIdx = min(foundNrIndices);
    lastNumIdx = max(foundNrIndices);
    contactNr = str2double(chName(firstNumIdx:lastNumIdx));
    electrodeName = chName(1:firstNumIdx-1);
    chNameCorr = strcat(electrodeName, num2str(contactNr));
end

function ignoreChannels = getIgnoreChannels()
    ignoreChannels = {'C3', 'C4', 'Cz', 'F3', 'F4', 'F7', 'F8', 'Fp1', 'FP1', 'Fp2', 'FP2', 'Fz',...
    'FZ', 'O1', 'O2', 'P3', 'P4', 'Pz', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'igger',...
    'ekg', 'ECG1', 'ECG2', 'EKG', 'EMG', 'emg', 'EOG','EOG_o', 'EOG_u', 'EOG_li', 'EOG_re',...
    'MKR1+', 'ecg1', 'ecg2', 'delg1', 'delg2', 'deld1', 'deld2', 'PULS+', 'BEAT+', 'SpO2+', 'MKR2+',...
    'ECG1', 'ECG2', 'delG1', 'delG2', 'delD1', 'delD2', 'DelG1', 'DelG2', 'DelD1', 'DelD2', ...
    'thor+', 'abdo+', 'xyz+', 'PULS+', 'BEAT+', 'SpO2+', 'MKR2+',...
    'MyD1', 'MyD2', 'MyG1', 'MyG2',...
    'Deld1', 'Deld2', 'Delg1', 'Delg2', 'my2', 'my1', 'jamG1', 'jamG2', 'jbeD1', 'jbeD2'};
end
