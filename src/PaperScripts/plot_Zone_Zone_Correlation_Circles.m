clc; clear all; close all;

paths = getFilesPaths();

patTablesPath = strcat(paths.workspacePath, 'CharacterizationTables\FeatureCharacterizationTables\');
groupTablesPath = strcat(paths.workspacePath, 'CharacterizationTables\GroupCharacterizationTablesAvg\');
analysisTablesPath = strcat(paths.workspacePath, 'AnalysisResults\Avg\');mkdir(analysisTablesPath);

biomarkersList =  {'HFO', 'iesHFO', 'IES'};
featuresList = {'rate', 'maxAmpl', 'variance', 'power'};

channSelStr = 'FlexK';
normStr = ''; %'', '_Normalized'};
zoneFormation = ''; % '', '_NonExclusive';
prctlTh = 75;
outcomeTh = 90;
corrLimitP = 0.001;


% Tables to read
spreadSheetNamePre = strcat(groupTablesPath, 'GroupCharacterization_', channSelStr, '_Pre', normStr,'.xls');
spreadSheetNamePost = strcat(groupTablesPath, 'GroupCharacterization_', channSelStr, '_Post', normStr,'.xls');
        
for bm = biomarkersList
    correlationAnalysis = {};
    biomarker = bm{1};
    groupTablePre = readtable(spreadSheetNamePre, 'Sheet', biomarker);
    groupTablePost = readtable(spreadSheetNamePost, 'Sheet', biomarker);
    featsNames = groupTablePre.Properties.VariableNames;

    if not(sum(strcmp(groupTablePre.chName, groupTablePost.chName)) == length(groupTablePost.chName))
        error(' newGroupAnalysis_Paired, Pre and Post EEG Channels don''t match')
    end
    
    plotZonesCorrelationCircles(groupTablePre, zoneFormation);
    
end

function plotZonesCorrelationCircles(groupTablePre, zoneFormation)
    close all;

    colors = {[0 0.4470 0.7410];...	
                [0.8500 0.3250 0.0980];...
                [0.9290 0.6940 0.1250];...
                [0.4940 0.1840 0.5560];...
                [0.4660 0.6740 0.1880];...
                [0.6350 0.0780 0.1840];...
                [0.3010 0.7450 0.9330]};
    patNames = groupTablePre.patName; 

    if not(isempty(zoneFormation))
        % get Zones
        rftcChanns = logical(groupTablePre.rftcSite);
        highEI_Channs = logical(groupTablePre.highEI);
        highConn_Channs = logical(groupTablePre.rftcConn);
        rftcStructure_Channs = logical(groupTablePre.rftcStruct);
        rftcLobe_Channs = logical(groupTablePre.rftcLobe);
        rftcHemisphere_Channs = logical(groupTablePre.rftcHemis);
    else    
        rftcChanns = logical(groupTablePre.rftcSite);
        highEI_Channs = logical(groupTablePre.highEI) & not(rftcChanns);
        highConn_Channs = logical(groupTablePre.rftcConn) & not(rftcChanns);
        rftcStructure_Channs = logical(groupTablePre.rftcStruct) & not(rftcChanns);
        rftcLobe_Channs = logical(groupTablePre.rftcLobe) & not(groupTablePre.rftcSite) & not(groupTablePre.rftcStruct);
        rftcHemisphere_Channs = logical(groupTablePre.rftcHemis)  & not(groupTablePre.rftcSite) & not(groupTablePre.rftcStruct) & not(groupTablePre.rftcLobe);
    end

    figure(1)
    set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w');
    tiledlayout(2, 3, 'TileSpacing', 'loose', 'Padding', 'tight');  
    
    %% Plot Channel Data
    [rftcChanns, newOrder] = sort(rftcChanns,'descend');
    patNames = patNames(newOrder);
    highEI_Channs = highEI_Channs(newOrder);
    highConn_Channs = highConn_Channs(newOrder);
    rftcStructure_Channs = rftcStructure_Channs(newOrder);
    rftcLobe_Channs = rftcLobe_Channs(newOrder);
    rftcHemisphere_Channs = rftcHemisphere_Channs(newOrder);

    allZones = [rftcChanns highConn_Channs highEI_Channs rftcStructure_Channs rftcLobe_Channs rftcHemisphere_Channs]';    
    labels = {'rftcSite', 'rftcConnected', 'highEI', 'rftcStructure', 'rftcLobe', 'rftcHemisphere'};
    
    for si = 1:size(allZones,1)
        nexttile
        dataP = allZones(si,:);
        axesLimits = zeros(2, size(dataP,2));
        axesLimits(1,:) = min(dataP,[],'all');
        axesLimits(2,:) = max(dataP,[],'all');
        axesLabels = 1:size(dataP,2);
        axesLabels = strsplit(num2str(axesLabels));

        spider_plot_R2019b(dataP,...
            'Web', 'off',...
            'AxesLabels', 'none',...
            'AxesInterval', 1,...
            'AxesPrecision', 0, ...
            'AxesLimits', axesLimits,...                        
            'AxesDisplay', 'one',...
            'AxesLabelsEdge', 'none',...
            'AxesPrecision', 1,...
            'AxesFontSize', 12,...
            'LineStyle', '-',...
            'LineWidth', 0.1,...
            'LineTransparency', 1,...
            'Marker', 'none',...
            'MarkerSize', 20,...
            'MarkerTransparency', 0.3,...
            'FillOption', {'on'},...
            'FillTransparency', [0.6],...
            'Color', colors{si});

        nrChanns = sum(dataP);
        nrPats = length(unique(patNames(dataP')));
        subtitleStr = {labels{si}; strcat("(", num2str(nrPats), ", ", num2str(nrChanns), ")")};
        subtitleStr = {labels{si}; strcat("(", num2str(nrChanns), ")")};
        %title(subtitleStr, 'FontSize', 14, 'Units', 'normalized', 'Position', [0.5, 0.6, 0])
        title(subtitleStr, 'FontSize', 14);

    end
    

    boxchartPath = strcat(params.analysisPerZonesTablePath, 'Zone_Zone_Correlations\'); mkdir(boxchartPath)
    boxchartFN = strcat(boxchartPath, 'ZoneAnalysis_', normStr, '_', params.patsSel);
    saveas(gcf,boxchartFN);
    boxchartFN = strcat(boxchartFN, '.jpg');
    saveas(gcf,boxchartFN);
    close();
end

function mcc = getMCC(vecA, vecB)
    vecA = logical(vecA);
    vecB = logical(vecB);

    tp = sum(vecA & vecB);
    tn = sum(not(vecA) & not(vecB));
    fp = sum(not(vecA) & vecB);
    fn = sum(vecA & not(vecB));
    
    mccA = (tp * tn) - (fp * fn);
    mccB = ((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));
    mcc = mccA / sqrt(mccB);
end

function allPatsKeepIdxs = perctThValsPerPat(patNameCol, vals, th)
    patNames = unique(patNameCol);
    allPatsKeepIdxs = zeros(length(vals),1);
    for pi = 1:length(patNames)
        patName = patNames{pi};
        patSelIdx = ismember(patNameCol, patName);
        patVals = vals(patSelIdx);
        patValsIdxs = 1:length(patVals);
        
        pi = randperm(length(patVals));
        patValsRand = patVals(pi);
        patValsIdxsRand = patValsIdxs(pi);
        [cid, centroids, sumd, D] = kmeans(patValsRand , 3, 'Start','cluster', 'Distance', 'sqeuclidean', 'MaxIter', 1000); %
        [maxVal, maxIdx] = max([mean(patValsRand(cid==1)), mean(patValsRand(cid==2)), mean(patValsRand(cid==3))]);
        patKeepIdxA = cid == maxIdx;        
        [sortedVec, sortingIdxs] = sort(patValsIdxsRand);
        patValsRand = patValsRand(sortingIdxs);
        patKeepIdxA = patKeepIdxA(sortingIdxs);
        
        interQR = iqr(patVals);
        quart1 = prctile(patVals, 25);
        patKeepIdxB = patVals >= prctile(patVals, th);
        patKeepIdxC = patVals >= median(patVals)+0.75*std(patVals);
        
        [min(patVals(patKeepIdxA)) prctile(patVals, th)  median(patVals)+0.75*std(patVals)]
        [sum(patKeepIdxA) sum(patKeepIdxB) sum(patKeepIdxC)]
        
        patKeepIdx = patKeepIdxC;
        
        allPatsKeepIdxs(patSelIdx) = patKeepIdx;        
    end
    allPatsKeepIdxs = logical(allPatsKeepIdxs);
end