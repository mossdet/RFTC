function [eiChannels, eiVals] = readEI(eegFilename)

paths = getFilesPaths();
workspacePath = paths.workspacePath;
eegFilesPath = paths.eegFilesPath;
eiFilesPath = paths.eiFilesPath;
filesMapPath = paths.filesMapPath;

T = readtable(filesMapPath);

idx = find(strcmp([T.EEG_Filename(:)], eegFilename)); % find line in Map where the EI Filename is at
eiFilename = T.EI_Filename{idx}

eiTable = [];

fid = fopen(strcat(eiFilesPath, eiFilename));
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
while ischar(tline)
    tline = strrep(tline, '\t', '    ');
    tline = strrep(tline, '	', '    ');
    
    spacesIdxs = strfind(tline, ' ');
    if isempty(spacesIdxs)
        eiFilename
        tline
        stop = 1; 
    end

    % Obtain and clean EI channel
    channName = tline(1:spacesIdxs(1));
    incorrectNrs = {'01','02','03','04','05','06','07','08','09'};
    corrctNrs = {'1','2','3','4','5','6','7','8','9'};
    newChannName = channName;
    
    for i = 1:length(incorrectNrs)
        newChannName = strrep(newChannName, incorrectNrs{i}, corrctNrs{i});
    end
    newChannName = strrep(newChannName, ' ', '');
    newChannName = strrep(newChannName, '	', '');

    % Obtain EI value
    diffSpacesIdxs = diff(spacesIdxs);
    spFidx = find(diffSpacesIdxs>1, 1, 'first');
    if isempty(spFidx)
        spFidx = length(diffSpacesIdxs)+1;
    end
    eiVal = tline(spacesIdxs(spFidx):end);
    eiVal = strrep(eiVal, ' ', '');
    eiVal = strrep(eiVal, '	', '');
    eiVal = str2num(eiVal);
    
    % Save ei data
    if ~isnan(eiVal)
        eiTable = cat(1, eiTable, {newChannName, eiVal});
    else
        eiFilename
        tline
        stop = 1;
    end
    
    tline = fgetl(fid);
end
fclose(fid);

eiChannels = eiTable(:,1);
eiVals = eiTable(:,2);


end