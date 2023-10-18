function plotIntraZoneWebPlots(plotIdx, params, plotData, plotDataLabels)   
    close all;
    pTh = 0.04;
    nrGroups = length(plotData);
    spi = 1;
    colors = [0 0.4470 0.7410 ;...	
                0.8500 0.3250 0.0980;...
                0.9290 0.6940 0.1250;...
                0.4940 0.1840 0.5560;...
                0.4660 0.6740 0.1880;...
                0.3010 0.7450 0.9330;...
                0.6350 0.0780 0.1840];            

    t = tiledlayout(2, nrGroups, 'TileSpacing', 'tight', 'Padding', 'tight');
    %t = tiledlayout('flow');

    set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w');

    for gi = 1:nrGroups*2
        nexttile
        ti = gi;
        if ti <= nrGroups
            biomarkerAct = plotData{ti,1};
            dataP = [rescale(biomarkerAct.in.pre,0,1)'; rescale(biomarkerAct.in.post,0,1)'];
            dataP = [biomarkerAct.in.pre'; biomarkerAct.in.post'];
            dataP = [sort(biomarkerAct.in.pre, 'descend')'; sort(biomarkerAct.in.post, 'descend')'];
        else
            ti = gi - nrGroups;
            biomarkerAct = plotData{ti,1};
            dataP = [rescale(biomarkerAct.out.pre,0,1)'; rescale(biomarkerAct.out.post,0,1)'];
            dataP = [biomarkerAct.out.pre'; biomarkerAct.out.post'];
            dataP = [sort(biomarkerAct.out.pre, 'descend')'; sort(biomarkerAct.out.post, 'descend')'];
        end
        dataP = circshift(dataP, floor(length(dataP)*0.75), 2);
        axesLimits = zeros(2, size(dataP,2));
        axesLimits(1,:) = min(dataP,[],'all');
        axesLimits(2,:) = max(dataP,[],'all');
        %axesLimits(1,:) = 0;
        %axesLimits(2,:) = 1;
        
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
            'LineWidth', 1,...
            'LineTransparency', 1,...
            'Marker', 'none',...
            'MarkerSize', 20,...
            'MarkerTransparency', 0.3,...
            'FillOption', {'off', 'off'},...
            'FillTransparency', [0.3 0.3],...
            'Color', [[1 0 0]; [0 0 1]]);
        
        zoneName = plotDataLabels{ti};
        title(zoneName, 'FontWeight', 'bold', 'FontSize', 14); % set title
        legend({'Pre', 'Post'})
        
        
        %%
%         ydata = [biomarkerAct.in.pre; biomarkerAct.in.post; biomarkerAct.out.pre; biomarkerAct.out.post;];
%         if gi == 1
%             ydata = [biomarkerAct.in.pre; biomarkerAct.in.post];
%             ydata = [ydata; zeros(length(biomarkerAct.out.pre),1); zeros(length(biomarkerAct.out.post),1)];
%         end
%         group = [ones(length(biomarkerAct.in.pre),1);ones(length(biomarkerAct.in.post),1)+1];
%         group = [group; ones(length(biomarkerAct.out.pre),1)+2; ones(length(biomarkerAct.out.post),1)+3];
%         positions = [1 1.25 2 2.25];
%         bp = boxplot(ydata, group, 'positions', positions,'Symbol','','Color', 'k');
%         title(plotDataLabels{gi}, 'FontWeight', 'bold', 'FontSize', 14);
% 
%         %% X Axis
%         xStrA = strcat("(", num2str(biomarkerAct.in.nrPats), ", ", num2str(biomarkerAct.in.nrChanns), ") ", plotDataLabels{gi});
%         xStrB = strcat("(", num2str(biomarkerAct.out.nrPats),", ", num2str(biomarkerAct.out.nrChanns), ") ", 'Rest');
%         set(gca,'xtick',[mean(positions(1:2)) mean(positions(3:4)) ])
%         set(gca,'xticklabel',{xStrA, xStrB});
        
    end
    biomrkrStr = strrep(params.biomarker, 'Vals', '');
    biomrkrStr = strrep(biomrkrStr, 'all', '');
    normStr = 'nonNormalized';
    if not(strcmp(params.normalization, ''))
        normStr = params.normalization;
    end
    
    titleStr = {strcat(params.patsSel, ", ", normStr); strcat("Intra-Zone Analysis, ", biomrkrStr, " ", params.feature)};
    title(t, titleStr,'FontSize',20);
    set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w');

    boxchartPath = strcat(params.analysisPerZonesTablePath, 'Intra_Zone_Analysis\', normStr, '\', params.patsSel, '\', params.feature, '\'); mkdir(boxchartPath)
    boxchartFN = strcat(boxchartPath, num2str(plotIdx), '_', strcat(params.biomarker, params.feature));
    %saveas(gcf,boxchartFN);
    boxchartFN = strcat(boxchartFN, '.jpg');
    %saveas(gcf,boxchartFN);
    close();
end
