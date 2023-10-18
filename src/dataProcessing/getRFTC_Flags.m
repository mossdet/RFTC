function rftcTable = getRFTC_Flags(paths, patName)

    rftcPatName = patName(1:strfind(patName, 'Inter')-2);
    rftcFilename = strcat(paths.rftcFlags, rftcPatName, '_RFTC_Channels.xls');
    readTable = readtable(rftcFilename);
    rftcTable.channelLabels = readTable.channelLabels;
    rftcTable.rftcVals = readTable.rftcVals;
end