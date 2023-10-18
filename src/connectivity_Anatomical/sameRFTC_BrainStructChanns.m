function [rftcStructFlag, usedAtlasName] = sameRFTC_BrainStructChanns(patTable)
    patName = patTable.patName{1};      
    patRFTC_Vals = patTable.rftcSite;
    patChannels = patTable.chName;
    [patParcelVals, usedAtlasName] = getBrainParcels(patName, patChannels); 

    rftcParcels = unique(patParcelVals(patRFTC_Vals));
    rftcParcels(strcmp(rftcParcels, 'N/A')) = [];

    rftcStructFlag = ismember(patParcelVals, rftcParcels);

end

function [channelParcelling, usedAtlasName] = getBrainParcels(patName, channels)
    
    usIdxs = strfind(patName, '_');
    patName = patName(1:usIdxs(3)-1);
    patName = strcat(patName, '.csv');

    %channelParcelling = getEstimatedBrainParcels(patName, channels);
    [channelParcelling, usedAtlasName]  = getGrenobleBrainParcels(patName, channels);
    %channelParcelling = getJuliaBrainParcels(patName, channels);

end

function channelParcelling = getJuliaBrainParcels(patName, channels)
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
        parcelName = 'N/A';
        if contains(chName, '-')
            analChannName = chName(1:strfind(chName, '-')-1);
        end
        chIdx = find(ismember(anatLocTableJJ.contact, analChannName));
        if isempty(chIdx)
            if contains(chName, '-')
                analChannName = chName(strfind(chName, '-')+1:end);
                chIdx = find(ismember(anatLocTableJJ.contact, analChannName));
                if not(isempty(chIdx))
                    parcelName = atlasParcels{chIdx};
                end
            end
        else
            parcelCode = atlasParcels(chIdx);
            switch parcelCode
                case 1
                    parcelName = 'hippocampus';
                case 2
                    parcelName = 'amygdala';
                case 3
                    parcelName = 'Parahippocampus';
                case 4
                    parcelName = 'Temporal neocortical';
                case 5
                    parcelName = 'Frontal';
                case 6
                    parcelName = 'Parietal';
                case 7
                    parcelName = 'Occipital';
                case 8
                    parcelName = 'Insula';
                case 0
                    parcelName = 'N/A';
                otherwise
                    patName
                    stopHere = 1;
            end
        end
        if not(isempty(strfind(parcelName, 'hite') > 0)) || not(isempty(strfind(parcelName, 'not in a') > 0))
            parcelName = 'N/A';
        end
        channelParcelling = cat(1, channelParcelling, parcelName);
    end
end

function [channelParcelling, usedAtlasName] = getGrenobleBrainParcels(patName, channels)
    channelParcelling = {};
    usedAtlasName = '';
    
    %%Anat Localization Table
    anatLocalizationGrenoblePath = 'F:\ForschungsProjekte\RFTC\Project_Files\ChannelLocalization\';
    listing = dir(anatLocalizationGrenoblePath);
    nrContents = length(listing);
    for ci = 1:nrContents
        listName = listing(ci).name(4:end);
        searchPatName = patName(4:end);
        if strcmp(listName, searchPatName)
            anatLocTableGrenoble = readtable(strcat(anatLocalizationGrenoblePath, listing(ci).name));
            break;
        end
    end
    for li = 1:length(anatLocTableGrenoble.contact)
        anatLocTableGrenoble.contact{li} = getCorrectChannelNames(anatLocTableGrenoble.contact{li});
    end

%     'Broadmann'
%     'BroadmannDilate'
%     'Hammers'
%     'AICHA'
%     if sum(ismember(anatLocTableGrenoble.Properties.VariableNames, 'Broadmann'))>0
%         'Broadmann'
%         atlasParcels = anatLocTableGrenoble.Broadmann;
%         atlasParcels = num2cell(atlasParcels);
%         atlasParcels = cellfun(@num2str, atlasParcels, 'UniformOutput', false);
% %     else
%     if sum(ismember(anatLocTableGrenoble.Properties.VariableNames, 'AAL'))>0
%         'AAL'
%         atlasParcels = anatLocTableGrenoble.AAL;
% %     else
%     if sum(ismember(anatLocTableGrenoble.Properties.VariableNames, 'AALDilate'))>0
%         'AALDilate'
%         atlasParcels = anatLocTableGrenoble.AALDilate;
%     if sum(ismember(anatLocTableGrenoble.Properties.VariableNames, 'BroadmannDilate'))>0
%         'BroadmannDilate'
%         atlasParcels = anatLocTableGrenoble.BroadmannDilate;
%         atlasParcels = num2cell(atlasParcels);
%         atlasParcels = cellfun(@num2str, atlasParcels, 'UniformOutput', false);
    if sum(ismember(anatLocTableGrenoble.Properties.VariableNames, 'AICHA'))>0
        usedAtlasName = 'AICHA';
        atlasParcels = anatLocTableGrenoble.AICHA;
    elseif sum(ismember(anatLocTableGrenoble.Properties.VariableNames, 'Freesurfer'))>0
        usedAtlasName = 'Freesurfer';
        atlasParcels = anatLocTableGrenoble.Freesurfer;
    elseif sum(ismember(anatLocTableGrenoble.Properties.VariableNames, 'MNI_DKT'))>0
        usedAtlasName = 'MNI_DKT';
        atlasParcels = anatLocTableGrenoble.MNI_DKT;
    elseif sum(ismember(anatLocTableGrenoble.Properties.VariableNames, 'MNI_Destrieux'))>0
        usedAtlasName = 'MNI_Destrieux';
        atlasParcels = anatLocTableGrenoble.MNI_Destrieux;
    elseif sum(ismember(anatLocTableGrenoble.Properties.VariableNames, 'Hammers'))>0
        usedAtlasName = 'Hammers';
        atlasParcels = anatLocTableGrenoble.Hammers;
    elseif sum(ismember(anatLocTableGrenoble.Properties.VariableNames, 'AAL'))>0
        usedAtlasName = 'AAL';
        atlasParcels = anatLocTableGrenoble.AAL;
    end

    for chi = 1:length(channels)
        chName = channels{chi};
        parcelName = 'N/A';
        if contains(chName, '-')
            analChannName = chName(1:strfind(chName, '-')-1);
        end
        chIdx = find(ismember(anatLocTableGrenoble.contact, analChannName));
        if isempty(chIdx) || not(isempty(strfind(atlasParcels{chIdx}, 'hite') > 0)) || not(isempty(strfind(atlasParcels{chIdx}, 'not in a') > 0))
            if contains(chName, '-')
                analChannName = chName(strfind(chName, '-')+1:end);
                chIdx = find(ismember(anatLocTableGrenoble.contact, analChannName));
                if not(isempty(chIdx))
                    if isa(atlasParcels,'double')
                        parcelName = atlasParcels(chIdx);
                    else
                        parcelName = atlasParcels{chIdx};
                    end
                end
            end
        else
            if isa(atlasParcels,'double')
                parcelName = atlasParcels(chIdx);
            else
                parcelName = atlasParcels{chIdx};
            end
        end
        if not(isempty(strfind(parcelName, 'hite') > 0)) || not(isempty(strfind(parcelName, 'not in a') > 0))
            parcelName = 'N/A';
        end
        channelParcelling = cat(1, channelParcelling, parcelName);
    end
end

function channelParcelling = getEstimatedBrainParcels(patName, channels)
    channelParcelling = {};
    %%Anat Localization Table
    anatLocalizationPath = 'F:\ForschungsProjekte\RFTC\Project_Files\ChannelLocalizationAtlases\';
    listing = dir(anatLocalizationPath);
    nrContents = length(listing);
    for ci = 1:nrContents
        if strfind(listing(ci).name, patName) > 0
            anatLocTable = readtable(strcat(anatLocalizationPath, listing(ci).name));
        end
    end

    anatLocTable.Channel;
    nrNAs = [sum(ismember(anatLocTable.JuelichHistologicalAtlas, 'N/A')), ...
        sum(ismember(anatLocTable.Harvard_OxfordCorticalStructuralAtlas, 'N/A')), ...
        sum(ismember(anatLocTable.Harvard_OxfordSubcorticalStructuralAtlas, 'N/A')), ...
        sum(ismember(anatLocTable.MNIStructuralAtlas, 'N/A'))];
    [minVal, minIdx] = min(nrNAs);
    minIdx = 3;
    switch minIdx
        case 1
            atlasParcels = anatLocTable.JuelichHistologicalAtlas;
        case 2
            atlasParcels = anatLocTable.Harvard_OxfordCorticalStructuralAtlas;
        case 3
            atlasParcels = anatLocTable.Harvard_OxfordSubcorticalStructuralAtlas;
        case 4
            atlasParcels = anatLocTable.MNIStructuralAtlas;
    end
    for chi = 1:length(channels)
        chName = channels{chi};
        parcelName = 'N/A';
        if contains(chName, '-')
            analChannName = chName(1:strfind(chName, '-')-1);
        end
        chIdx = find(ismember(anatLocTable.Channel, analChannName));
        if isempty(chIdx)
            if contains(chName, '-')
                analChannName = chName(strfind(chName, '-')+1:end);
                chIdx = find(ismember(anatLocTable.Channel, analChannName));
                if not(isempty(chIdx))
                    parcelName = atlasParcels{chIdx};
                end
            end
        else
            parcelName = atlasParcels{chIdx};
        end
        if not(isempty(strfind(parcelName, 'hite') > 0)) || not(isempty(strfind(parcelName, 'not in a') > 0))
            parcelName = 'N/A';
        end
        channelParcelling = cat(1, channelParcelling, parcelName);
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