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
allPatsDateDiff = [];
for fileIdx = 1:size(files,1)
    filename = strcat(eegFilesPath, files{fileIdx});
    hdr = in_fopen_micromed(filename);
    [~, subjName, ~] = fileparts(filename)

    patAcq = {subjName, hdr.acquisition.year,  hdr.acquisition.month, hdr.acquisition.day,...
    hdr.acquisition.hour, hdr.acquisition.min, hdr.acquisition.sec};
    allPatsDate = cat(1, allPatsDate, patAcq);

    if mod(fileIdx, 2) == 0        
        t1 = datetime(allPatsDate{fileIdx-1,2}, allPatsDate{fileIdx-1,3}, allPatsDate{fileIdx-1,4}, allPatsDate{fileIdx-1,5}, allPatsDate{fileIdx-1,6}, allPatsDate{fileIdx-1,7});
        t2 = datetime(allPatsDate{fileIdx,2}, allPatsDate{fileIdx,3}, allPatsDate{fileIdx,4}, allPatsDate{fileIdx,5}, allPatsDate{fileIdx,6}, allPatsDate{fileIdx,7});
        dateDiff = t2-t1;
        allPatsDateDiff = cat(1, allPatsDateDiff, [days(dateDiff), hours(dateDiff), minutes(dateDiff)]);
    end
end

nrFiles = size(allPatsDate,1);
postRFTC_DurHoursVec = [];  
for fi = 2:2:nrFiles
	allPatsDate{1,2}
	t1 = datetime(allPatsDate{fi-1,2}, allPatsDate{fi-1,3}, allPatsDate{fi-1,4}, allPatsDate{fi-1,5}, allPatsDate{fi-1,6}, allPatsDate{fi-1,7});
	t2 = datetime(allPatsDate{fi,2}, allPatsDate{fi,3}, allPatsDate{fi,4}, allPatsDate{fi,5}, allPatsDate{fi,6}, allPatsDate{fi,7});
	t11=datevec(datenum(t1));
	t22=datevec(datenum(t2));
	time_interval_in_seconds = etime(t22,t11)
	time_interval_in_hours = time_interval_in_seconds/3600
	postRFTC_DurHoursVec = cat(1, postRFTC_DurHoursVec, time_interval_in_hours);
end
