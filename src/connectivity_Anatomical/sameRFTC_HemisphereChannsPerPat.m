function allPatsKeepIdxs = sameRFTC_HemisphereChannsPerPat(groupTablePre)
    patNameCol = groupTablePre.patName;
    rftcCols = groupTablePre.rftcVals;
    hemisphereCols = groupTablePre.hemisphere;

    patNames = unique(patNameCol);
    allPatsKeepIdxs = zeros(length(rftcCols),1);
    for pi = 1:length(patNames)
        patName = patNames{pi};
        patSelIdx = ismember(patNameCol, patName);
        
        patRFTC_Vals = rftcCols(patSelIdx);
        patHemVals = hemisphereCols(patSelIdx);
        patHemispheres = unique(patHemVals);
        rftcHemispheres = unique(patHemVals(logical(patRFTC_Vals)));
        patKeepIdx = ismember(patHemVals, rftcHemispheres);
        allPatsKeepIdxs(patSelIdx) = patKeepIdx;        
    end
    allPatsKeepIdxs = logical(allPatsKeepIdxs);
end