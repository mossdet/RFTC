clc; clear all; close all;
paths = getFilesPaths();
preFiles  = getPreFiles();%getAllFiles, getPreFiles, getPostFiles
newTablePath = 'F:\ForschungsProjekte\RFTC\Project_Files\JuliaLocalizationTables\'; mkdir(newTablePath);

for fileIdx = 1:size(preFiles,1)
    prePatName = preFiles{fileIdx}; prePatName = prePatName(1:length(prePatName)-4);
    anatLocalizationGrenoblePath = 'F:\ForschungsProjekte\RFTC\Project_Files\ChannelLocalization\';
    patNameSep = strfind(prePatName, '_');
    tabPatName = prePatName(1:patNameSep(3)-1); tabPatName(1:3) = 'Gre';
    anatLocTableGrenoble = readtable(strcat(anatLocalizationGrenoblePath, tabPatName, '.csv'));
    
    channSizes = strfind(anatLocTableGrenoble.contact, '-');
    delIdx = cellfun(@length, channSizes) > 0;
    anatLocTableGrenoble(delIdx,:) = [];
    
    BrainRegion_JJ = repmat({''}, size(anatLocTableGrenoble.contact));
    anatLocTableGrenoble = addvars(anatLocTableGrenoble, BrainRegion_JJ);
    
    newTableFN = strcat(newTablePath, prePatName, '_JJ.xls');
    writetable(anatLocTableGrenoble, newTableFN);
end