function plotZone_ZoneCorrelationCircles_Paper(params, groupTablePre)
    close all;

    colors = {[0 0.4470 0.7410];...	
                [0.8500 0.3250 0.0980];...
                [0.9290 0.6940 0.1250];...
                [0.4940 0.1840 0.5560];...
                [0.4660 0.6740 0.1880];...
                [0.6350 0.0780 0.1840];...
                [0.3010 0.7450 0.9330]};

    % get Zones
    patNames = groupTablePre.patName; 
    rftcChanns = logical(groupTablePre.rftcVals);
    highEI_Channs = logical(perctThValsPerPat(groupTablePre.patName, groupTablePre.eiVals, params.eiTh));
    highConn_Channs = logical(perctThValsPerPat(groupTablePre.patName, groupTablePre.rftcElectroPhysioConnect, params.connTh));
    rftcStructure_Channs = logical(sameRFTC_BrainParcelChannsPerPat(groupTablePre));
    rftcLobe_Channs = logical(sameRFTC_BrainLobeChannsPerPat(groupTablePre));
    rftcHemisphere_Channs = logical(sameRFTC_HemisphereChannsPerPat(groupTablePre));
    
    if(params.exlusive)
        highEI_Channs(rftcChanns) = false(1,1);
        highConn_Channs(rftcChanns) = false(1,1);
        rftcStructure_Channs(rftcChanns) = false(1,1);
        rftcLobe_Channs(rftcChanns | rftcStructure_Channs) = false(1,1);
        rftcHemisphere_Channs(rftcChanns|rftcStructure_Channs|rftcLobe_Channs) = false(1,1);
        %rftcHemisphere_Channs(rftcChanns|rftcStructure_Channs) = false(1,1);
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
    labels = {'rftc', 'rftcConnected', 'highEI', 'rftcStructure', 'rftcLobe', 'rftcHemisphere'};
    
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
        %title(subtitleStr, 'FontSize', 14, 'Units', 'normalized', 'Position', [0.5, 0.6, 0])
        title(subtitleStr, 'FontSize', 14);

    end
    
    normStr = 'nonNormalized';
    if not(strcmp(params.normalization, ''))
        normStr = params.normalization;
    end
    
  params.analysisPerZonesTablePath
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