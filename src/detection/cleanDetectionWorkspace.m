function cleanDetectionWorkspace(hfoDetectorFolder)
filesToKeep = {'svm_spindle.dat', 'MOSSDET_c.exe', 'svm_fr.dat', 'svm_ies.dat', 'svm_ripple.dat', '.', '..'};
listing = dir(hfoDetectorFolder);
    for ci = 1:length(listing)
        contentName = listing(ci).name;
        keepContent = ismember(contentName, filesToKeep(:));
        if not(keepContent)
            if listing(ci).isdir
                rmdir(strcat(hfoDetectorFolder, contentName), 's');
            else
                delete(strcat(hfoDetectorFolder, contentName));
            end
        end
    end
end