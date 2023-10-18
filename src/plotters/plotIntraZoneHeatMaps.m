function plotIntraZoneHeatMaps(plotIdx, params, plotData, plotDataLabels)   
    close all;
    pTh = 0.04;
    nrGroups = length(plotData);
    spi = 1;
    colors = {[0 0.4470 0.7410];...	
                [0.8500 0.3250 0.0980];...
                [0.9290 0.6940 0.1250];...
                [0.4940 0.1840 0.5560];...
                [0.4660 0.6740 0.1880];...
                [0.3010 0.7450 0.9330];...
                [0.6350 0.0780 0.1840]};
            

    ha = tight_subplot(1, nrGroups, [.01 .03],[.2 .1],[.05 .01]);
    set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w');

    for gi = 1:nrGroups
        axes(ha(gi));
                
        biomarkerAct = plotData{gi,1};
        activity = [rescale(biomarkerAct.in.pre,0,1); rescale(biomarkerAct.in.post,0,1); rescale(biomarkerAct.out.pre,0,1); rescale(biomarkerAct.out.post,0,1)];
        group = [ones(length(biomarkerAct.in.pre),1);ones(length(biomarkerAct.in.post),1)+1];
        group = [group; ones(length(biomarkerAct.out.pre),1)+2; ones(length(biomarkerAct.out.post),1)+3];
        channNrs = [1:length(biomarkerAct.in.pre) 1:length(biomarkerAct.in.post)];
        channNrs = [channNrs 1:length(biomarkerAct.out.pre) 1:length(biomarkerAct.out.post)];
        positions = [2 2.25 3 3.25];
        
        x = repmat(positions(1), [1 length(biomarkerAct.in.pre)]);
        x = [x repmat(positions(2), [1 length(biomarkerAct.in.post)])];
        x = [x repmat(positions(3), [1 length(biomarkerAct.out.pre)])];
        x = [x repmat(positions(4), [1 length(biomarkerAct.out.post)])];

        y = channNrs;
        sz = 5;
        c = activity;

        %bp = boxplot(ydata, group, 'positions', positions,'Symbol','','Color', 'k');
        colormap turbo;
        scatter(x,y,sz,c);
        colorbar;
        title(plotDataLabels{gi}, 'FontWeight', 'bold', 'FontSize', 14);
        xtickangle(30);
        set(gca,'xtick',[mean(positions(1:2)) mean(positions(3:4)) ])
        ylim([min(y) max(y)]);
        zlim([min(activity) max(activity)])
        set(gca,'color',[0 0 0]);
        
        spi = spi+1;
    end
    biomrkrStr = strrep(params.biomarker, 'Vals', '');
    biomrkrStr = strrep(biomrkrStr, 'all', '');
    normStr = 'nonNormalized';
    if not(strcmp(params.normalization, ''))
        normStr = params.normalization;
    end
    
    titleStr = {strcat(params.patsSel, ", ", normStr); strcat("Intra-Zone Analysis, ", biomrkrStr, " ", params.feature)};
    sgtitle(titleStr,'FontSize',20);
    set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w');

    boxchartPath = strcat(params.analysisPerZonesTablePath, 'Intra_Zone_Analysis\', normStr, '\', params.patsSel, '\', params.feature, '\'); mkdir(boxchartPath)
    boxchartFN = strcat(boxchartPath, num2str(plotIdx), '_', strcat(params.biomarker, params.feature));
    saveas(gcf,boxchartFN);
    boxchartFN = strcat(boxchartFN, '.jpg');
    saveas(gcf,boxchartFN);
    close();
end
