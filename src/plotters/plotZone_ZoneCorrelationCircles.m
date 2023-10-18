function plotZone_ZoneCorrelationCircles(params, groupTablePre)
    close all;
    
    colors = {[0 0.4470 0.7410];...	
    [0.8500 0.3250 0.0980];...
    [0.9290 0.6940 0.1250];...
    [0.4940 0.1840 0.5560];...
    [0.4660 0.6740 0.1880];...
    [0.3010 0.7450 0.9330];...
    [0.6350 0.0780 0.1840]};

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
    
    length(unique(patNames(rftcHemisphere_Channs)))
    allZones = [rftcChanns highEI_Channs highConn_Channs rftcStructure_Channs rftcLobe_Channs rftcHemisphere_Channs]';
    
    % define Labels
    labels = {'RFTC', 'highEI', 'highRFTC-Conn.', 'rftcParcel', 'rftcLobe', 'rftcHemisphere'};
    
    % Make sure all patients are being considered
    nrPats = 0;
    nrPats = nrPats + length(unique(groupTablePre.patName(rftcChanns)));
    nrPats = nrPats + length(unique(groupTablePre.patName(highEI_Channs)));
    nrPats = nrPats + length(unique(groupTablePre.patName(highConn_Channs)));
    nrPats = nrPats + length(unique(groupTablePre.patName(rftcStructure_Channs)));
    nrPats = nrPats + length(unique(groupTablePre.patName(rftcLobe_Channs)));
    nrPats = nrPats + length(unique(groupTablePre.patName(rftcHemisphere_Channs)));
    nrPats  = nrPats / 6
    
    % build correlation matrix
    corrMatrix = zeros(size(allZones,1), size(allZones,1));
    for ri = 1:size(corrMatrix,1)
        for ci = 1:size(corrMatrix,2)
            mccVal = getMCC(allZones(ri,:), allZones(ci,:));
            if mccVal < 0
                mccVal  = 0;
            end
            if ri == ci
                mccVal = 0;
            end
            corrMatrix(ri,ci) = mccVal;
            strcat(labels{ri}, " - ", labels{ci})
        end
    end
    maxMCC = max(corrMatrix,[],'all');
    for ri = 1:size(corrMatrix,1)
        for ci = 1:size(corrMatrix,2)
            if ri == ci
                corrMatrix(ri,ci) = maxMCC;
            end
        end
    end

    figure(1)
    set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w');
    tiledlayout(3, 4, 'TileSpacing', 'loose', 'Padding', 'tight');  
    
    %% Plot Correlation Matrix
    nexttile([3 2]);
    %subplot(3,4,[1,2,5,6,9,10])

    imagesc(corrMatrix);
    set(gca, 'XTick', 1:size(corrMatrix,2)); % center x-axis ticks on bins
    set(gca, 'YTick', 1:size(corrMatrix,1)); % center y-axis ticks on bins
    set(gca, 'XTickLabel', labels); % set x-axis labels
    set(gca, 'YTickLabel', labels); % set y-axis labels
    xtickangle(30);
    ax = gca;
    ax.FontSize = 18; 
    ax.FontWeight = 'bold';
    title('Zones Correlation Matrix', 'FontSize', 20); % set title
    colormap('summer'); % Choose jet or any other color scheme
    c = colorbar; % 
    c.FontSize = 18;
    c.FontWeight = 'bold';
    for ri = 1:size(corrMatrix,1)
        for ci = 1:size(corrMatrix,2)
            label = num2str(corrMatrix(ri, ci),3);
            if ri == ci
                label = "1";
            end
            text(ri-0.25, ci, label,'FontSize', 16, 'FontWeight', 'bold');
        end
    end   
    
    %% Plot Channel Data
    %% Plot Channel Data
    spidxs = [3 4 7 8 11 12];  
    [rftcChanns, newOrder] = sort(rftcChanns,'descend');
    patNames = patNames(newOrder);
    highEI_Channs = highEI_Channs(newOrder);
    highConn_Channs = highConn_Channs(newOrder);
    rftcStructure_Channs = rftcStructure_Channs(newOrder);
    rftcLobe_Channs = rftcLobe_Channs(newOrder);
    rftcHemisphere_Channs = rftcHemisphere_Channs(newOrder);

    allZones = [rftcChanns highEI_Channs rftcLobe_Channs highConn_Channs  rftcHemisphere_Channs rftcStructure_Channs]';    
    labels = {'RFTC', 'highEI', 'rftcLobe', 'highRFTC-Conn.', 'rftcHemisphere', 'rftcStructure'};
    
    allZones = [rftcChanns rftcStructure_Channs highConn_Channs rftcLobe_Channs highEI_Channs rftcHemisphere_Channs]';
    labels = {'RFTC', 'rftcStructure', 'highRFTC-Conn.',  'rftcLobe', 'highEI', 'rftcHemisphere'};
    
    for si = 1:size(allZones,1)
        nexttile
        %subplot(3,4,spidxs(si))
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
        
        interQR = iqr(patVals);
        quart1 = prctile(patVals, 25);
        %patKeepIdxA = patVals >= quart1 + interQR;
        %patKeepIdxB = patVals >= prctile(patVals, th);
        patKeepIdxC = patVals >= median(patVals)+0.75*std(patVals);
%         if sum(patKeepIdxA) ~= sum(patKeepIdxB)
%             [sum(patKeepIdxA) sum(patKeepIdxB)]
%             stop = 1;
%         end
        patKeepIdx = patKeepIdxC;
        
        allPatsKeepIdxs(patSelIdx) = patKeepIdx;        
    end
    allPatsKeepIdxs = logical(allPatsKeepIdxs);
end