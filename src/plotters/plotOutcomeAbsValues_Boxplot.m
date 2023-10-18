function plotOutcomeDifferences_Boxplot(predictTableDiff, predictTablePVal, params)
    close all;
    pTh = 0.05;
    featName = params.feature;
    biomarkers = fieldnames(predictTableDiff);
    zones = fieldnames(predictTableDiff.(biomarkers{1}));
    nrBiomarkers = length(biomarkers);
    nrZones = length(zones);
    colors = {[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.9290 0.6940 0.1250]; [0.4940 0.1840 0.5560];...
              [0.4660 0.6740 0.1880]; [0.3010 0.7450 0.9330]; [0.6350 0.0780 0.1840]};
            
    tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w');
    
    %%
    for bi = 1
        biomarkerAct = predictTableDiff.(biomarkers{bi});
        biomarkerName = biomarkers{bi};
        biomarkerName = strrep(biomarkerName, 'ia', 'IES');
        biomarkerName = strrep(biomarkerName, 'Vals', '');
        biomarkerName = strrep(biomarkerName, 'all', '');
        biomarkerName = strrep(biomarkerName, '_', '');
        featName = strrep(featName, 'OccRate', 'Occ.Rate');
            
        group = [];
        ydata = [];
        outcomes = fieldnames(biomarkerAct.(zones{1}));
        
        groupIdx = 1;
        for zi = 1:nrZones
            for oi = 1:length(outcomes)
                ydata = cat(1, ydata, biomarkerAct.(zones{zi}).(outcomes{oi}));
                group = cat(1, group, zeros(length(biomarkerAct.(zones{zi}).(outcomes{oi})),1)+groupIdx);
                groupIdx = groupIdx+1;
            end
        end
        positions = sort([1:nrZones, (1:nrZones)+0.25], 'ascend');

        nexttile
        bp = boxplot(ydata, group, 'positions', positions,'Symbol','','Color', 'k');
        %title(biomarkerName, 'FontWeight', 'bold', 'FontSize', 14);

        %% X Axis
%         xStrA = strcat("(", num2str(biomarkerAct.in.nrPats), ", ", num2str(biomarkerAct.in.nrChanns), ") ", predictTableDiffLabels{zi}{1});
%         xStrB = strcat("(", num2str(biomarkerAct.out.nrPats),", ", num2str(biomarkerAct.out.nrChanns), ") ", 'Rest');
%         xStrC = strcat(predictTableDiffLabels{zi}{1}, " ", "Δ");
%         xStrD = strcat("Rest ","Δ");
%         xStrAll = {xStrA, xStrB, xStrC, xStrD};
        xTicksPositions = [];
        for xtpi = 1:2:length(positions)
            xTicksPositions = [xTicksPositions; mean(positions(xtpi:xtpi+1))];
        end
        %set(gca,'xtick',[mean(positions(1:2)) mean(positions(3:4)) mean(positions(5:6)) mean(positions(7:8)) mean(positions(9:10)) mean(positions(11:12))])
        set(gca,'xtick',xTicksPositions)
        set(gca,'xticklabel',zones);

        %% Y Axis
        h = findobj(gcf, 'tag', 'Upper Adjacent Value'); 
        up_adj = cell2mat(get(h,'YData')); up_adj=up_adj(:,1);up_adj=unique(up_adj);
        h = findobj(gcf, 'tag', 'Lower Adjacent Value');
        low_adj = cell2mat(get(h,'YData')); low_adj=low_adj(:,1);low_adj=unique(low_adj);
        h = findobj(gcf, 'tag', 'Upper Whisker');
        up_whisker = cell2mat(get(h,'YData'));up_whisker=unique(up_whisker);
        ylim([min(low_adj) max(up_adj)])
        ylabel(strcat('\Delta', featName), 'FontSize', 16)
        title(biomarkerName, 'FontSize', 16);

        %% Legend
        alphas = repmat([0.1 0.8], 1, 6);
        h = findobj(gca,'Tag','Box');
        selColorIdx = 0;
        for j=length(h):-1:1
            if mod(j,2) == 0
                selColorIdx = selColorIdx+1;
            end
            selColor = colors{selColorIdx};
            patch(get(h(j),'XData'),get(h(j),'YData'), selColor, 'FaceAlpha',alphas(j));
        end
        c = get(gca, 'Children');
        %hleg1 = legend(c(1:2), 'Positive PSO', 'Negative PSO','Location','northeast');


        %% Annotations
        groupsID = unique(group);
        groupsIdx = 1;
        for zi = 1:nrZones
            selectedGroup = [groupsID(groupsIdx) groupsID(groupsIdx+1)] 
            impData = ydata(group == groupsID(groupsIdx));
            nonImpData = ydata(group == groupsID(groupsIdx+1));
            
            if zi == 1 || zi == 3 || zi == 4
                [p, h] = ranksum(impData, nonImpData, 'tail','right');
            else
                [p, h] = ranksum(impData, nonImpData, 'tail','left');
            end
            [p, h] = ranksum(impData, nonImpData);

            pVal = num2str(p,3);
            labels1 = {strcat("p=", pVal)};

            up_adj = up_adj(up_adj~=0);
            xPos = positions(zi*2);
            yPos = mean(up_adj); %mean(ydata(ydata~=0));
            ht = text(xPos, yPos, labels1,'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','bold','BackgroundColor','black', 'Color', 'yellow','FontSize',12);
            set(ht,'Rotation',45);
            if(p >= pTh)
                ht.Color = 'red';
            end
            groupsIdx = groupsIdx+2;
        end

        ax = gca;
        ax.XAxis.FontSize = 12;
        xtickangle(45);
        set(findobj(gcf,'type','axes'),'FontWeight','Bold');
    end
    titleStr = {"Activity Differences Between PSO Groups"};
    sgtitle(titleStr,'FontSize',12);
    
    normStr = 'nonNormalized';
    if not(strcmp(params.normalization, ''))
        normStr = params.normalization;
    end
    
    boxchartPath = strcat(params.outcomePredictionAnalysisPath, 'BoxplotAnalysis\'); mkdir(boxchartPath)
    boxchartFN = strcat(boxchartPath, strcat('BoxplotAnalysis_', params.feature), '_', normStr);
    saveas(gcf,boxchartFN);
    boxchartFN = strcat(boxchartFN, '.jpg');
    saveas(gcf,boxchartFN);
    close();
    
end