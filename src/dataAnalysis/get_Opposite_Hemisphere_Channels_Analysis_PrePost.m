function get_Opposite_Hemisphere_Channels_Analysis_PrePost(paths, selDetector, analysisPerZonesTablePath, thList, features, normOptions, freqBandConnList, patsSelList, outcomeTh, connTh, eiTh)
    analysisType = {'Opposite_Hemisphere_Channels'};    
    analysis = {};
    threshold = [];
    nrPats = [];
    nrChanns = [];
    iaVals = [];
    allHFOVals = []; iesHFOVals = []; isolHFOVals = [];
    allRippleVals = []; iesRippleVals = []; isolRippleVals = [];
    allFR_Vals = []; iesFR_Vals = []; isolFR_Vals = [];
            
    for th = thList
        for ni = 1:length(normOptions)
            for fi = 1:length(features)
                for psi = 1:length(patsSelList)
                    patsSel = patsSelList{psi};
                    feature = features{fi};
                    normalization = normOptions{ni};
                    analysisSubType = feature;

                    selDetectorFolder = selDetector;
                    if th > 0
                        selDetectorFolder = strcat('MOSSDET_Depurated\Th', num2str(th));
                    end
                    groupAnalysisTablesFilePath = strcat(paths.workspacePath, 'GroupAnalysis_ChannelCharacterizationTables\', selDetectorFolder, '\');
                    preGroupTableFN = strcat(groupAnalysisTablesFilePath, normalization, 'GroupAnalysis_ChannelCharacterization_Pre', selDetector, '.xls');
                    postGroupTableFN = strcat(groupAnalysisTablesFilePath, normalization, 'GroupAnalysis_ChannelCharacterization_Post', selDetector, '.xls');
                    groupTablePre = readtable(preGroupTableFN, 'Sheet', feature);
                    groupTablePost = readtable(postGroupTableFN, 'Sheet', feature);

                    %% Select channels to compare
                    tableDelIdxs = (sameRFTC_HemisphereChannsPerPat(groupTablePre.patName, groupTablePre.rftcVals, groupTablePre.hemisphere));
                    groupTablePre(tableDelIdxs, :) = [];
                    groupTablePost(tableDelIdxs, :) = [];
                    
                    tableDelIdxs = groupTablePre.rftcVals > 0;
                    groupTablePre(tableDelIdxs, :) = [];
                    groupTablePost(tableDelIdxs, :) = [];
                    
                    %% Select outcomes
                    if strcmp(patsSel, 'allPatients')
                        analysisSubType = strcat(analysisSubType, '_', patsSel);
                    elseif strcmp(patsSel, 'improvedPatients')
                        analysisSubType = strcat(analysisSubType, '_', patsSel);
                        tableDelIdxs = groupTablePre.outcome <= outcomeTh;
                        groupTablePre(tableDelIdxs, :) = [];
                        tableDelIdxs = groupTablePost.outcome <= outcomeTh;
                        groupTablePost(tableDelIdxs, :) = [];
                    elseif strcmp(patsSel, 'nonImprovedPatients')
                        analysisSubType = strcat(analysisSubType, '_', patsSel);
                        tableDelIdxs = groupTablePre.outcome > outcomeTh;
                        groupTablePre(tableDelIdxs, :) = [];
                        tableDelIdxs = groupTablePost.outcome > outcomeTh;
                        groupTablePost(tableDelIdxs, :) = [];
                    end
                    
                    analysisSubType = strcat(analysisSubType, '_', normalization);
                    %%

                    analysis = cat(1, analysis, analysisSubType);
                    threshold = cat(1, threshold, th);
                    nrPats = cat(1, nrPats, length(unique(groupTablePost.patName)));
                    nrChanns = cat(1, nrChanns, length(groupTablePost.patName));

                    iaVals = cat(1, iaVals, signrank(groupTablePre.iaVals, groupTablePost.iaVals, 'tail', 'right'));

                    allHFOVals = cat(1, allHFOVals, signrank(groupTablePre.allHFOVals, groupTablePost.allHFOVals, 'tail', 'right'));
                    iesHFOVals = cat(1, iesHFOVals, signrank(groupTablePre.iesHFOVals, groupTablePost.iesHFOVals, 'tail', 'right'));
                    isolHFOVals = cat(1, isolHFOVals, signrank(groupTablePre.isolHFOVals, groupTablePost.isolHFOVals, 'tail', 'right'));

                    allRippleVals = cat(1, allRippleVals, signrank(groupTablePre.allRippleVals, groupTablePost.allRippleVals, 'tail', 'right'));
                    iesRippleVals = cat(1, iesRippleVals, signrank(groupTablePre.iesRippleVals, groupTablePost.iesRippleVals, 'tail', 'right'));
                    isolRippleVals = cat(1, isolRippleVals, signrank(groupTablePre.isolRippleVals, groupTablePost.isolRippleVals, 'tail', 'right'));

                    allFR_Vals = cat(1, allFR_Vals, signrank(groupTablePre.allFR_Vals, groupTablePost.allFR_Vals, 'tail', 'right'));
                    iesFR_Vals = cat(1, iesFR_Vals, signrank(groupTablePre.iesFR_Vals, groupTablePost.iesFR_Vals, 'tail', 'right'));
                    isolFR_Vals = cat(1, isolFR_Vals, signrank(groupTablePre.isolFR_Vals, groupTablePost.isolFR_Vals, 'tail', 'right'));
                end
            end
        end    
        zonesAnalysisTable = table(analysis, threshold, nrPats, nrChanns, iaVals, allHFOVals, iesHFOVals, isolHFOVals, allRippleVals, iesRippleVals, isolRippleVals, allFR_Vals, iesFR_Vals, isolFR_Vals); 
        spaceT = zonesAnalysisTable(1,:);
        for ti = 1:length(spaceT.Properties.VariableNames)
            spaceT.(spaceT.Properties.VariableNames{ti}) = nan(1,1);
        end
        zonesAnalysisTableSpaces = [];
        for si = 0:(size(zonesAnalysisTable,1)/3)-1
            idx = 3*si;
            zonesAnalysisTableSpaces = cat(1, zonesAnalysisTableSpaces, zonesAnalysisTable(idx+1:idx+3,:), spaceT);
        end
        
        newSelDetector = selDetector;
        if th > 0
            analysisPerZonesTablePath(end) = [];
            analysisPerZonesTablePath = strcat(analysisPerZonesTablePath, '_Depurated\', 'Th', num2str(th), '\');mkdir(analysisPerZonesTablePath);
            newSelDetector = strcat(selDetector, '_Th', num2str(th));
        end
        analysisPerZonesTableFN = strcat(analysisPerZonesTablePath, 'AnalysisPerHemishpere_', newSelDetector,'.xls');
        writetable(zonesAnalysisTableSpaces, analysisPerZonesTableFN, 'Sheet', analysisType{1});
    end
end

function allPatsKeepIdxs = sameRFTC_HemisphereChannsPerPat(patNameCol, rftcCols, hemisphereCols)
    patNames = unique(patNameCol);
    allPatsKeepIdxs = zeros(length(rftcCols),1);
    for pi = 1:length(patNames)
        patName = patNames{pi};
        patSelIdx = ismember(patNameCol, patName);
        patRFTC_Vals = rftcCols(patSelIdx);
        patHemVals = hemisphereCols(patSelIdx);
        patHemispheres = unique(patHemVals)
        rftcHemispheres = unique(patHemVals(logical(patRFTC_Vals)))
        patKeepIdx = ismember(patHemVals, rftcHemispheres);
        allPatsKeepIdxs(patSelIdx) = patKeepIdx;        
    end
    allPatsKeepIdxs = logical(allPatsKeepIdxs);
end