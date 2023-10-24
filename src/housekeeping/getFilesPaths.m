function paths = getFilesPaths()

    mk = 1;
    dlp = 0;
    paths = [];
    
    computer = getenv('computername');

    if strcmp(computer, 'ZKJ-20-6DQR393')
        paths.workspacePath = 'C:\Users\lachner\Documents\MATLAB_ver2\';
        paths.eegFilesPath = 'C:\Users\lachner\Documents\RFTC\Project_Files\PatientFilesMicromed\AllPatients\';
        paths.eiFilesPath = 'C:\Users\lachner\Documents\MATLAB_ver2\EI\';
        paths.filesMapPath= 'C:\Users\lachner\Documents\RFTC\Project_Files\PatientLists\GrenobleFilesMap.xlsx';
        paths.hfoDetectorFolder = 'C:\Users\lachner\Documents\RFTC\MATLAB_ver2\MOSSDET_c\';
        paths.rftcFlags = 'C:\Users\lachner\Documents\RFTC\Project_Files\RFTC_Channels\Done_Only_RFTC_BipolarMontage_Labeled\';
        paths.anatLocalizationGrenoblePath = 'C:\Users\lachner\Documents\RFTC\Project_Files\ChannelLocalization_JJ\';
    elseif strcmp(computer, 'LAPTOP-TFQFNF6U')
        paths.workspacePath = 'F:\ForschungsProjekte\RFTC\RFTC_HFO\';
        paths.eegFilesPath = 'F:\ForschungsProjekte\RFTC\Project_Files\PatientFilesMicromed\AllPatients\';
        paths.eiFilesPath = 'F:\ForschungsProjekte\RFTC\RFTC_HFO\OtherData\EI\';
        paths.filesMapPath= 'F:\ForschungsProjekte\RFTC\Project_Files\PatientLists\GrenobleFilesMap.xlsx';
        paths.hfoDetectorFolder = 'F:\ForschungsProjekte\RFTC\RFTC_HFO\MOSSDET_c\';
        paths.rftcFlags = 'F:\ForschungsProjekte\RFTC\Project_Files\RFTC_Channels\Done_Only_RFTC_BipolarMontage_Labeled\';
        paths.mniCoordinates = 'F:\ForschungsProjekte\RFTC\Project_Files\ChannelLocalization\';
        paths.anatLocalizationGrenoblePath = 'F:\ForschungsProjekte\RFTC\Project_Files\ChannelLocalization_JJ\';
    end
    
end