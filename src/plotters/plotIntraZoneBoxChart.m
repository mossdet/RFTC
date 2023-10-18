function plotIntraZoneBoxChart(plotIdx, params, plotData, plotDataLabels)   
    pTh = 0.01;
    nrGroups = length(plotData);
    spi = 1;
    colors = {[0 0.4470 0.7410];...	
                [0.8500 0.3250 0.0980];...
                [0.9290 0.6940 0.1250];...
                [0.4940 0.1840 0.5560];...
                [0.4660 0.6740 0.1880];...
                [0.6350 0.0780 0.1840];...
                [0.3010 0.7450 0.9330]};
            

    %ha = tight_subplot(1, nrGroups, [.01 .03],[.2 .1],[.05 .01]);
    tiledlayout(2, nrGroups, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    useFS = 16;
    titleFS = useFS+4;

    for gi = 1:nrGroups
        
        %axes(ha(gi));
        nexttile
        
        %subplot(1, nrGroups, spi)
        
        xgroupdata = [];
        ydata = [];
        biomarkerAct = plotData{gi,1};
        
        ydata = [biomarkerAct.in.pre; biomarkerAct.in.post; biomarkerAct.out.pre; biomarkerAct.out.post; biomarkerAct.in.diff; biomarkerAct.out.diff];
%         if gi == 2
%             ydata = [biomarkerAct.in.pre; biomarkerAct.in.post];
%             ydata = [ydata; zeros(length(biomarkerAct.out.pre),1); zeros(length(biomarkerAct.out.post),1)];
%         end
        group = [ones(length(biomarkerAct.in.pre),1);ones(length(biomarkerAct.in.post),1)+1];
        group = [group; ones(length(biomarkerAct.out.pre),1)+2; ones(length(biomarkerAct.out.post),1)+3];
        group = [group; ones(length(biomarkerAct.in.diff),1)+4; ones(length(biomarkerAct.out.diff),1)+5];

        positions = [1 1.25 2 2.25 3 3.25];
        positions = [0.5 0.6 0.8 0.9 1.1 1.2];
        bp = boxplot(ydata, group, 'positions', positions,'Symbol','','Color', 'k');
        title(plotDataLabels{gi}, 'FontWeight', 'bold', 'FontSize', useFS);

        %% X Axis
        %xStrA = {plotDataLabels{gi}{1}; strcat("(", num2str(biomarkerAct.in.nrPats), ", ", num2str(biomarkerAct.in.nrChanns), ") ")};
        %xStrB = {'Rest'; strcat("(", num2str(biomarkerAct.out.nrPats),", ", num2str(biomarkerAct.out.nrChanns), ") ")};
        xStrA = plotDataLabels{gi}{1};
        xStrB = 'Rest';
        xStrC = strcat(plotDataLabels{gi}{1}, " ", "Δ");
        xStrD = "Rest Δ";
        
        xTickPositions= [mean(positions(1:2)) mean(positions(3:4)) positions(5) positions(6)];
        xticks(xTickPositions);
        set(gca,'xticklabel',{xStrA, xStrB, xStrC, xStrD});
        xtickangle(70);
        ax = gca;
        ax.XAxis.FontSize = useFS;
        
%         txt= {xStrA,'','',''};
%         xTT=text(xTickPositions,(min(ylim)+1)*ones(1,length(xTickPositions)),txt,...
%             'HorizontalAlignment','right','VerticalAlignment','bottom',...
%             'FontSize', useFS-2, 'FontWeight','bold');
%         set(xTT,'Rotation',70);
%         
%         txt= {'',xStrB,'',''};
%         xTT=text(xTickPositions,(min(ylim))*ones(1,length(xTickPositions)),txt,...
%             'HorizontalAlignment','right','VerticalAlignment','bottom',...
%             'FontSize', useFS-2, 'FontWeight','bold');
%         set(xTT,'Rotation',70);
%         
%         txt= {'','',xStrC,xStrD};
%         xTT=text(xTickPositions,(min(ylim))*ones(1,length(xTickPositions)),txt,...
%         'HorizontalAlignment','right','VerticalAlignment','bottom',...
%         'FontSize', useFS-2, 'FontWeight','bold');
%         set(xTT,'Rotation',70);
                
        %% Y Axis
        h = findobj(gcf, 'tag', 'Upper Adjacent Value'); 
        up_adj = cell2mat(get(h,'YData')); up_adj=up_adj(:,1);up_adj=unique(up_adj);
        h = findobj(gcf, 'tag', 'Lower Adjacent Value');
        low_adj = cell2mat(get(h,'YData')); low_adj=low_adj(:,1);low_adj=unique(low_adj);
        h = findobj(gcf, 'tag', 'Upper Whisker');
        up_whisker = cell2mat(get(h,'YData'));up_whisker=unique(up_whisker);
        
        ylim([min(low_adj) max(up_adj)])
        params.biomarker = strrep(params.biomarker, 'OccRate', 'Occ.Rate');
        if gi == 1
            ylabel(params.feature, 'FontSize', useFS)
        end
        
        %% Legend
        alphas = [0.15 0.65 0.15 0.65 0.15 0.65];
        h = findobj(gca,'Tag','Box');
        for j=1:length(h)
            selColor = colors{gi};
            if j <= 2
                selColor = 'k';
            end
           patch(get(h(j),'XData'),get(h(j),'YData'), selColor, 'FaceAlpha',alphas(j));
        end
        c = get(gca, 'Children');
        %hleg1 = legend(c(1:2), 'Pre-RFTC', 'Post-RFTC','Location','northeast');
        %hleg1.FontSize = useFS;
        
            schriftFarbe  = [0.4667    0.6745    0.1882]; 'yellow'; [0.4660 0.6740 0.1880];
            bgColor = [1, 1, 1];
            %bgColor = 'none';

        
            %% Annotations
            up_adj = up_adj(up_adj~=0);
            labels1 = {strcat("p=", num2str(biomarkerAct.in.prePostP,2))};
            xPos = positions(2);
            yPos = mean(up_adj) + mean(up_adj)/2; %mean(ydata(ydata~=0));
            ht = text(xPos, yPos, labels1,'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','bold','BackgroundColor',bgColor, 'Color', schriftFarbe,'FontSize', useFS);
            set(ht,'Rotation',45);
            if(biomarkerAct.in.prePostP > pTh)
                ht.Color = 'red';
            end

            labels2 = {strcat("p=", num2str(biomarkerAct.out.prePostP,2))};
            xPos = positions(4);
            yPos = mean(up_adj) + mean(up_adj)/2; %mean(ydata);
            ht = text(xPos, yPos, labels2,'HorizontalAlignment','center','VerticalAlignment','bottom', 'FontWeight', 'bold','BackgroundColor',bgColor, 'Color', schriftFarbe,'FontSize', useFS);
            set(ht,'Rotation',45);
            if(biomarkerAct.out.prePostP > pTh)
                ht.Color = 'red';
            end

            labels3 = {strcat("\Delta p=", num2str(biomarkerAct.inOutP,2))};
            xPos = mean([positions(5) positions(6)]);
            yPos = mean(up_adj) + mean(up_adj)/2;
            ht = text(xPos, yPos, labels3,'HorizontalAlignment','center','VerticalAlignment','bottom', 'FontWeight', 'bold','BackgroundColor',bgColor, 'Color', schriftFarbe,'FontSize', useFS);
            set(ht,'Rotation',45);
            if(biomarkerAct.inOutP > pTh)
                ht.Color = 'red';
            end
        
        set(findobj(gcf,'type','axes'),'FontWeight','Bold');
        spi = spi+1;
    end
    biomrkrStr = strrep(params.biomarker, 'iaVals', 'IES');
    biomrkrStr = strrep(biomrkrStr, 'Vals', '');
    biomrkrStr = strrep(biomrkrStr, 'all', '');
    normStr = 'nonNormalized';
    if not(strcmp(params.normalization, ''))
        normStr = params.normalization;
    end
    
    titleStr = {strcat(params.patsSel, ", ", normStr); strcat("Intra-Zone Analysis, ", biomrkrStr, " ", params.feature)};
    titleStr = {strcat(biomrkrStr, " ", params.feature)};
    sgtitle(titleStr,'FontSize',titleFS);
    set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w');

    normStr = params.normalization;
    if strcmp(normStr, '')
        normStr = 'nonNormalized';
    end
    boxchartPath = strcat(params.analysisPerZonesTablePath, 'Intra_Zone_Analysis\', normStr, '\', params.patsSel, '\', params.feature, '\'); mkdir(boxchartPath)
    boxchartPath = strcat(params.analysisPerZonesTablePath, 'WideNoLegend\'); mkdir(boxchartPath)
    boxchartFN = strcat(boxchartPath, num2str(plotIdx), '_', strcat(params.biomarker, params.feature), '_', normStr, '_', params.patsSel);
    %saveas(gcf,boxchartFN);
    boxchartFN = strcat(boxchartFN, '.jpg');
    saveas(gcf,boxchartFN);
    close();
end
