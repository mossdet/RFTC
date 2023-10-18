function plotInterZoneBoxChart(plotIdx, params, plotData, plotDataLabels)
    %for interZone analysis there is no rest for all channels data
    plotData(1) = [];
    plotDataLabels(1) = [];
    
    nrGroups = length(plotData);
    spi = 1;
    colors = {[0 0.4470 0.7410];...	
                [0.8500 0.3250 0.0980];...
                [0.9290 0.6940 0.1250];...
                [0.4940 0.1840 0.5560];...
                [0.4660 0.6740 0.1880];...
                [0.3010 0.7450 0.9330];...
                [0.6350 0.0780 0.1840]};

    absMin = 0;
    absMax = 0;
    for gi = 1:nrGroups
        biomarkerAct = plotData{gi,1};
        allDat = [ biomarkerAct.in.diff; biomarkerAct.out.diff];
        if absMin > min(allDat)
            absMin = min(allDat);
        end
        if absMax < max(allDat)
            absMax = max(allDat);
        end
    end
    absMin = floor(absMin);
    absMax = ceil(absMax);

    for gi = 1:nrGroups
        subplot(1, nrGroups, spi)
        
        xgroupdata = [];
        ydata = [];
        biomarkerAct = plotData{gi,1};
        ydata = cat(1, ydata, biomarkerAct.in.diff, biomarkerAct.out.diff);
        xgroupdata = cat(1, xgroupdata, zeros(length(biomarkerAct.in.diff),1)+1, zeros(length(biomarkerAct.out.diff),1)+2);
        boxchart(xgroupdata, ydata, 'BoxFaceColor',colors{gi},'MarkerColor',colors{gi});
        spi = spi+1;
        
        [p, h] = ranksum(biomarkerAct.in.diff, biomarkerAct.out.diff);
        legendStr = strcat('p = ', num2str(p, 3));
        legend(legendStr, 'FontSize', 12);
        %xlabel('Grouping according to the RFTC Channels')
        if gi == 1
            ylabel({'Activity \Delta', '(preRFTC - postRFTC)'}, 'FontSize', 16)
        else
            set(gca,'yticklabel',[]);
        end
        ylim([absMin absMax]);
        xticks([1 2]);
        xStrA = strcat("(", num2str(biomarkerAct.in.nrPats), ", ", num2str(biomarkerAct.in.nrChanns), ") ", plotDataLabels{gi});
        xStrB = strcat("(", num2str(biomarkerAct.out.nrPats),", ", num2str(biomarkerAct.out.nrChanns), ") ", 'Rest');
        set(gca,'XTickLabel',{xStrA, xStrB});
        ax = gca;
        ax.XAxis.FontSize = 10;
        xtickangle(20);
        set(findobj(gcf,'type','axes'),'FontWeight','Bold');
    end
    biomrkrStr = strrep(params.biomarker, 'Vals', '');
    biomrkrStr = strrep(biomrkrStr, 'all', '');
    normStr = 'nonNormalized';
    if not(strcmp(params.normalization, ''))
        normStr = params.normalization;
    end
    
    titleStr = {strcat(params.patsSel, ", ", normStr); strcat("Inter-Zone Analysis, ", biomrkrStr, " ", params.feature)};
    sgtitle(titleStr,'FontSize',20);
    set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w');

    boxchartPath = strcat(params.analysisPerZonesTablePath, 'Inter_Zone_Analysis\', normStr, '\', params.patsSel, '\', params.feature, '\'); mkdir(boxchartPath)
    boxchartFN = strcat(boxchartPath, num2str(plotIdx), '_', strcat(params.biomarker, params.feature));
    saveas(gcf,boxchartFN);
    boxchartFN = strcat(boxchartFN, '.jpg');
    saveas(gcf,boxchartFN);
    close();
end

function plotInterZoneBoxChartWithAll(plotIdx, params, plotData, plotDataLabels)
    %for interZone analysis there is no rest for all channels data
%     plotData(1) = [];
%     plotDataLabels(1) = [];
    
    nrGroups = length(plotData);
    spi = 1;
    colors = {[0 0.4470 0.7410];...	
                [0.8500 0.3250 0.0980];...
                [0.9290 0.6940 0.1250];...
                [0.4940 0.1840 0.5560];...
                [0.4660 0.6740 0.1880];...
                [0.3010 0.7450 0.9330];...
                [0.6350 0.0780 0.1840]};

    absMin = 0;
    absMax = 0;
    for gi = 1:nrGroups
        biomarkerAct = plotData{gi,1};
        allDat = [ biomarkerAct.in.diff; biomarkerAct.out.diff];
        if absMin > min(allDat)
            absMin = min(allDat);
        end
        if absMax < max(allDat)
            absMax = max(allDat);
        end
    end
    absMin = floor(absMin);
    absMax = ceil(absMax);

    for gi = 1:nrGroups
        subplot(1, nrGroups, spi)
        
        xgroupdata = [];
        ydata = [];
        biomarkerAct = plotData{gi,1};
        ydata = cat(1, biomarkerAct.in.diff, biomarkerAct.out.diff);
        if gi == 1
            ydata = cat(1, biomarkerAct.in.diff, zeros(length(biomarkerAct.out.diff),1));
        end
        xgroupdata = cat(1, zeros(length(biomarkerAct.in.diff),1)+1, zeros(length(biomarkerAct.out.diff),1)+2);
        boxchart(xgroupdata, ydata, 'BoxFaceColor',colors{gi},'MarkerColor',colors{gi});
        spi = spi+1;
        
        [p, h] = ranksum(biomarkerAct.in.diff, biomarkerAct.out.diff);
        legendStr = strcat('p = ', num2str(p, 3));
        legend(legendStr, 'FontSize', 12);
        %xlabel('Grouping according to the RFTC Channels')
        if gi == 1
            ylabel({'Activity \Delta', '(preRFTC - postRFTC)'}, 'FontSize', 16)
        else
            set(gca,'yticklabel',[]);
        end
        ylim([absMin absMax]);
        xticks([1 2]);
        xStrA = strcat("(", num2str(biomarkerAct.in.nrPats), ", ", num2str(biomarkerAct.in.nrChanns), ") ", plotDataLabels{gi});
        xStrB = strcat("(", num2str(biomarkerAct.out.nrPats),", ", num2str(biomarkerAct.out.nrChanns), ") ", 'Rest');
        if gi == 1
            xStrB = '';
        end
        set(gca,'XTickLabel',{xStrA, xStrB});
        ax = gca;
        ax.XAxis.FontSize = 10;
        xtickangle(20);
        set(findobj(gcf,'type','axes'),'FontWeight','Bold');
    end
    biomrkrStr = strrep(params.biomarker, 'Vals', '');
    biomrkrStr = strrep(biomrkrStr, 'all', '');
    titleStr = {'Inter-Zone Analysis'; strcat(biomrkrStr, " ", params.feature)};
    sgtitle(titleStr,'FontSize',20);
    set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w');
    normStr = 'nonNormalized';
    if not(strcmp(params.normalization, ''))
        normStr = params.normalization;
    end
    boxchartPath = strcat(params.analysisPerZonesTablePath, 'Inter_Zone_Analysis\', normStr, '\', params.patsSel, '\', params.feature, '\'); mkdir(boxchartPath)
    boxchartFN = strcat(boxchartPath, num2str(plotIdx), '_', strcat(params.biomarker, params.feature));
    saveas(gcf,boxchartFN);
    boxchartFN = strcat(boxchartFN, '.jpg');
    saveas(gcf,boxchartFN);
    close();
end