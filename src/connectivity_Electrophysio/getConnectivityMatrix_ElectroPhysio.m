clc; clear all; close all;
paths = getFilesPaths();
files  = getPreFiles();%getAllFiles, getPreFiles, getPostFiles
eegCorrelationPath = 'F:\ForschungsProjekte\RFTC\Project_Files\PatientFilesMicromed\Connectivity\';
tablesFilePath = strcat(paths.workspacePath, 'ChannelCharacterizationTables\MOSSDET\');

for fileIdx = 1:size(files,1)
    eegFilename = strcat(paths.eegFilesPath, files{fileIdx});
    [origFilepath, patName, ext] = fileparts(eegFilename);
    patName
    
    deltaConn = getCorrelationResults(tablesFilePath, eegCorrelationPath, patName, 'delta');
    thetaConn = getCorrelationResults(tablesFilePath, eegCorrelationPath, patName, 'theta');
    alphaConn = getCorrelationResults(tablesFilePath, eegCorrelationPath, patName, 'alpha');
    betaConn = getCorrelationResults(tablesFilePath, eegCorrelationPath, patName, 'beta');
    gammaConn = getCorrelationResults(tablesFilePath, eegCorrelationPath, patName, 'gamma');
    highGammaConn = getCorrelationResults(tablesFilePath, eegCorrelationPath, patName, 'highGamma');
    rippleConn = getCorrelationResults(tablesFilePath, eegCorrelationPath, patName, 'ripple');
    frConn = getCorrelationResults(tablesFilePath, eegCorrelationPath, patName, 'fr');
    
    %get max correlation all freq. bands    
    maxCorrelAllBands = deltaConn;
    maxCorrelAllBands.corr = max(maxCorrelAllBands.corr, thetaConn.corr);
    maxCorrelAllBands.corr = max(maxCorrelAllBands.corr, alphaConn.corr);
    maxCorrelAllBands.corr = max(maxCorrelAllBands.corr, betaConn.corr);
    maxCorrelAllBands.corr = max(maxCorrelAllBands.corr, gammaConn.corr);
    maxCorrelAllBands.corr = max(maxCorrelAllBands.corr, highGammaConn.corr);
    maxCorrelAllBands.corr = max(maxCorrelAllBands.corr, rippleConn.corr);
    maxCorrelAllBands.corr = max(maxCorrelAllBands.corr, frConn.corr);

    %get mean correlation all freq. bands    
    meanCorrelAllBands = deltaConn;
    meanCorrelAllBands.corr = (deltaConn.corr + thetaConn.corr + alphaConn.corr + betaConn.corr + gammaConn.corr + highGammaConn.corr + rippleConn.corr + frConn.corr)/8;
    
    electroPhysioConnectPath = strcat(paths.workspacePath, 'ConnectivityMatricesAllBands\');mkdir(electroPhysioConnectPath);
    electroPhysioConnectFN = strcat(electroPhysioConnectPath, patName);
    
    save(electroPhysioConnectFN, 'deltaConn', 'thetaConn', 'alphaConn', 'betaConn', 'gammaConn',...
                                'highGammaConn', 'rippleConn', 'frConn', 'maxCorrelAllBands', 'meanCorrelAllBands');    
end

function bandConn = getCorrelationResults(tablesFilePath, eegCorrelationPath, patName, freqBand)
    correlationFN = [];
    if strcmp(freqBand, 'delta')
        correlationFN = strcat(eegCorrelationPath, patName, '.TRC_algo-H2_hp-1Hz_lp-4Hz.mat');
    elseif strcmp(freqBand, 'theta')
        correlationFN = strcat(eegCorrelationPath, patName, '.TRC_algo-H2_hp-4Hz_lp-8Hz.mat');
    elseif strcmp(freqBand, 'alpha')
        correlationFN = strcat(eegCorrelationPath, patName, '.TRC_algo-H2_hp-8Hz_lp-15Hz.mat');
    elseif strcmp(freqBand, 'beta')
        correlationFN = strcat(eegCorrelationPath, patName, '.TRC_algo-H2_hp-15Hz_lp-30Hz.mat');
    elseif strcmp(freqBand, 'gamma')
        correlationFN = strcat(eegCorrelationPath, patName, '.TRC_algo-H2_hp-30Hz_lp-45Hz.mat');
    elseif strcmp(freqBand, 'highGamma')
        correlationFN = strcat(eegCorrelationPath, patName, '.TRC_algo-H2_hp-55Hz_lp-90Hz.mat');
    elseif strcmp(freqBand, 'ripple')
        correlationFN = strcat(eegCorrelationPath, patName, '.TRC_algo-H2_hp-80Hz_lp-250Hz.mat');
    elseif strcmp(freqBand, 'fr')
        correlationFN = strcat(eegCorrelationPath, patName, '.TRC_algo-H2_hp-250Hz_lp-500Hz.mat');
    else
        "Error"
    end

    if isempty(correlationFN)
        stopHere = 1;
    end
    avgCorr = getAvgCorrelation(correlationFN);

    %get channel names
    correlationResults = load(correlationFN);
    channelNames = correlationResults.electrode_names;
    channelNames = editChannelNames(channelNames)';
    channelNames = sort(channelNames);
    %Get RFTC information
    tableFN = strcat(tablesFilePath, patName, '_', 'ChannelCharacterization_MOSSDET.xls');
    table = readtable(tableFN, 'Sheet', 'OccRate');

    %Delete and sort channels 
    delChanns = find(not(ismember(channelNames, table.channelLabels)));
    channelNames(delChanns) = [];
    avgCorr(delChanns,:) = [];
    avgCorr(:,delChanns) = [];
    newChannOrder = [];
    for chi = 1:length(table.channelLabels)
        tableChann = table.channelLabels{chi};
        newChannOrder = cat(1, newChannOrder, find(ismember(channelNames, tableChann)));
    end
    avgCorr = avgCorr(newChannOrder, newChannOrder);
    channelNames = channelNames(newChannOrder);
    rftcVals = table.rftcVals;
    
    bandConn.corr = avgCorr;
    bandConn.channelNames = channelNames;
    bandConn.rftcVals = rftcVals;    
end

function avgCorrelResults = getAvgCorrelation(correlationFN)
    correlationResults = load(correlationFN);
    avgCorrelResults = zeros(size(correlationResults.aw_h2,1), size(correlationResults.aw_h2,2));
    for ri = 1:size(correlationResults.aw_h2,1)
        for ci = 1:size(correlationResults.aw_h2, 2)
            avgCorrelResults(ri, ci) = mean(correlationResults.aw_h2(ri,ci,:));
        end
    end    
end

function editedChannels = editChannelNames(channelNames)
    numbers = ['0','1','2','3','4','5','6','7','8','9'];
    editedChannels = channelNames;
    nrChannels = length(channelNames);
    for chi = 1:nrChannels
        chLabel = channelNames{chi};
        channA = chLabel(1:strfind(chLabel, '-')-1);
        channB = chLabel(strfind(chLabel, '-')+1:end);

        %ChannA
        foundNrIndices = [];
		for ni = 1:length(numbers)
			strIdx = strfind(channA, numbers(ni));
			foundNrIndices = cat(2, foundNrIndices, strIdx);
		end
		firstNumIdx = min(foundNrIndices);
		lastNumIdx = max(foundNrIndices);
        contactNrA = str2double(channA(firstNumIdx:lastNumIdx));
		electrodeNameA = channA(1:firstNumIdx-1);
        
        %ChannB
        foundNrIndices = [];
		for ni = 1:length(numbers)
			strIdx = strfind(channB, numbers(ni));
			foundNrIndices = cat(2, foundNrIndices, strIdx);
		end
		firstNumIdx = min(foundNrIndices);
		lastNumIdx = max(foundNrIndices);
        contactNrB = str2double(channB(firstNumIdx:lastNumIdx));
		electrodeNameB = channB(1:firstNumIdx-1);
        editedChannels{chi} = strcat(electrodeNameA, num2str(contactNrA), '-', electrodeNameB, num2str(contactNrB));
    end
end