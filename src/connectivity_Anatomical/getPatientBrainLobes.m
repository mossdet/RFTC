function channelLobeLoc = getPatientBrainLobes(patTable)
    patName = patTable.patName{1};
    patChannels = unique(patTable.chName);
    patChannels_LobeLocation = getBrainLobes(patName, patChannels);
    patLobes = unique(patChannels_LobeLocation);
    patLobes(strcmp(patLobes, 'N/A')) = [];
    
    channelLobeLoc = "";
    for pli = 1:length(patLobes)
        lobeStr = patLobes{pli};
        nrChannsInLobe = sum(strcmp(patChannels_LobeLocation, lobeStr));
        channelLobeLoc = strcat(channelLobeLoc, strcat(lobeStr, " (", num2str(nrChannsInLobe), "), "));
    end
    channelLobeLoc{1}(end)=[];
    channelLobeLoc{1}(end)=[];
end

function channelLobes = getBrainLobes(patName, channels)

    channelLobes = getJuliaBrainLobes(patName, channels);

end

function channelParcelling = getJuliaBrainLobes(patName, channels)
    channelParcelling = {};
    
    %%Anat Localization Table
    anatLocalizationGrenoblePath = 'F:\ForschungsProjekte\RFTC\Project_Files\ChannelLocalization_JJ\';
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
        chName = channels{chi};
        lobeName = 'N/A';
        if contains(chName, '-')
            analChannName = chName(1:strfind(chName, '-')-1);
        end
        chIdx = find(ismember(anatLocTableJJ.contact, analChannName));
        if isempty(chIdx)
            if contains(chName, '-')
                analChannName = chName(strfind(chName, '-')+1:end);
                chIdx = find(ismember(anatLocTableJJ.contact, analChannName));
                if not(isempty(chIdx))
                    lobeName = atlasParcels{chIdx};
                end
            end
        else
            parcelCode = atlasParcels(chIdx);
            switch parcelCode
                case 1
                    lobeName = 'Mesiotemporal';
                case 2
                    lobeName = 'Mesiotemporal';
                case 3
                    lobeName = 'Mesiotemporal';
                case 4
                    lobeName = 'Temporal';
                case 5
                    lobeName = 'Frontal';
                case 6
                    lobeName = 'Parietal';
                case 7
                    lobeName = 'Occipital';
                case 8
                    lobeName = 'Insula';
                case 0
                    lobeName = 'N/A';
                otherwise
                    patName
                    stopHere = 1;
            end
        end
        if not(isempty(strfind(lobeName, 'hite') > 0)) || not(isempty(strfind(lobeName, 'not in a') > 0))
            lobeName = 'N/A';
        end
        channelParcelling = cat(1, channelParcelling, lobeName);
    end
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