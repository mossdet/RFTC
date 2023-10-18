clc; clear all; close all;

[path,~,~] = fileparts(mfilename('fullpath'));
cutIdx = strfind(path, '\');
workspacePath = path(1:cutIdx(end));
cd(workspacePath)
addpath(genpath(workspacePath));
startup

paths = getFilesPaths();
workspacePath = paths.workspacePath;
eegFilesPath = paths.eegFilesPath;
eiFilesPath = paths.eiFilesPath;
filesMapPath = paths.filesMapPath;

files  = getAllFiles();
ignoreChannels  = [];%getIgnoreChannels();

for fileIdx = 1:size(files,1)
    filename = strcat(eegFilesPath, files{fileIdx});
    anatLocalizationGrenoblePath = paths.anatLocalizationGrenoblePath;

    hdr = ft_read_header(filename);
    [~, patName, ~] = fileparts(filename);

    listing = dir(anatLocalizationGrenoblePath);
    nrContents = length(listing);
    for ci = 1:nrContents
        if strfind(listing(ci).name, patName) > 0
            anatLocTableJJ = readtable(strcat(anatLocalizationGrenoblePath, listing(ci).name));
        end
    end
    for li = 1:length(anatLocTableJJ.contact)
        anatLocTableJJ.contact{li} = getCorrectChannelNames(anatLocTableJJ.contact{li});
    end
    atlasParcels = anatLocTableJJ.BrainRegion_JJ;

    for chi = 1:length(channels)
    end


    %signals = ft_read_data(filename);
end