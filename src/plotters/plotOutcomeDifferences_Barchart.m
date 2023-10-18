function plotOutcomeDifferences_Barchart(predictTableDiff, predictTablePVal, params)
    close all;
    pTh = 0.05;
    featName = params.feature;
    biomarkers = fieldnames(predictTableDiff);
    zones = fieldnames(predictTableDiff.(biomarkers{1}));
    nrBiomarkers = length(biomarkers);
    nrZones = length(zones);
    colors = {[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.9290 0.6940 0.1250]; [0.4940 0.1840 0.5560];...
              [0.4660 0.6740 0.1880]; [0.6350 0.0780 0.1840]; [0.3010 0.7450 0.9330]};
          
    centroidsOccRate = [0.726333400541575,0.0749739318202013,0.203307836060067,0.137016019650521,-0.00444910917728204,-0.0101146555969047;...
                        0.562825373958852,0.132123224609946,0.195901323987377,0.405743935617079,0.0279907854575977,-0.0139998862279889];
                    
    centroidsPower = [0.405827310122140,0.100803571678327,0.116258707370680,0.156700002745292,0.0281906086915328,-0.0540889213651772;...
                      0.105688770609661,-0.00140554562027272,-0.0347026508085116,-0.0415437957297722,0.0146090645097350,-0.0292212164319009];

    centroidsOccRate = [0.686991158438075, 0.0201907417737201, 0.0939234193956363, 0.0291835796321246,	-0.148996700977968,	-0.156572737228407;...
                        0.455596648929743, -0.134212314457678, 0.0998073268752962, 0.433588264401153, -0.108565733314542, -0.181303339229077];
    
    centroidsPower = [0.0971845951392287, -0.0824271142382000, -0.123429403130625, -0.133521071911289, -0.118074565532689, -0.164451901359607;...
                        0.150583013713516, -0.134534412196181, -0.102090581272958, -0.0395487174967840, -0.157211341400218, -0.282108814552270];
                    
    centroidsOccRate = [2.90635924702511	0.238622574620531	0.622962657941679	0.475286765604628	-0.0426005935095951	-0.0299115146795989;...
    1.65321216577806	0.509165568023412	0.754143939493825	1.24905726130928	0.0786441602961765	0.0373938229829019];
              
    centroidsPower = [0.674986035287437	-0.0112837456141711	-0.104439604442446	-0.0851733174292842	-0.0754210538382535	-0.0794376236112551;...
    0.915340766776227	0.246679019711350	0.321667641314846	0.376367662081187	0.240494348523058	-0.161689318484002];

    centroids = centroidsOccRate;
    if strcmp(featName, 'Power')
        centroids = centroidsPower;
    end

            
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
            
        ydata = [];
        stdDevVals = [];
        annotsCoord = [];
        outcomes = fieldnames(biomarkerAct.(zones{1}));
        
        positions = sort([1:nrZones, (1:nrZones)+0.25], 'ascend');
        bi = 1;
        for zi = 1:nrZones
            for oi = 1:length(outcomes)
                ydata = cat(2, ydata, mean(biomarkerAct.(zones{zi}).(outcomes{oi})));
                stdDevVals = cat(2, stdDevVals, std(biomarkerAct.(zones{zi}).(outcomes{oi})));
                b = bar(positions(bi), ydata(bi),0.22); hold on;
                b.FaceColor = 'flat';
                b.FaceAlpha = 0.7;
                if mod(oi,2) == 0
                    b.FaceAlpha = 0.2;
                    %annotsCoord = cat(1, annotsCoord, [mean(positions(bi-1:bi)),max(ydata(bi-1:bi))+max(stdDevVals(bi-1:bi))/2]);
                    annotsCoord = cat(1, annotsCoord, [mean(positions(bi-1:bi)), mean(ydata(bi-1:bi))]);

                end
                %errorbar(positions(bi), ydata(bi), stdDevVals(bi), stdDevVals(bi), '.k', 'LineWidth',0.5)
                
                b.CData(1,:) = colors{zi};
                xlim([positions(1)-1 positions(end)+0.5])
                
                xtips1 = b.XEndPoints;
                ytips1 = b.YEndPoints;
                labels1 = num2str(b.YData, '%4.3f');
                text(xtips1,ytips1,labels1,'HorizontalAlignment','center', 'VerticalAlignment','bottom')

                bi = bi+1;
            end
        end
        annotsCoord(:,2) = max(annotsCoord(:,2)) + 0.2;
        
        %% X Ticks
        xTicksPositions = [];
        for xtpi = 1:2:length(positions)
            xTicksPositions = [xTicksPositions; mean(positions(xtpi:xtpi+1))];
        end
        set(gca,'xtick',xTicksPositions)
        set(gca,'xticklabel',zones);
        
        %% Y Axis
        ylim([min(ydata-stdDevVals/10) max(ydata+stdDevVals/4)])
        %ylim([min(ydata) max(ydata)])
        ylabel(strcat("Mean \Delta", featName), 'FontSize', 16)

        %% Annotations
        for zi = 1:nrZones
            impData = biomarkerAct.(zones{zi}).imp;
            nonImpData = biomarkerAct.(zones{zi}).nonImp;            
            [p, h] = ranksum(impData, nonImpData);

            pVal = num2str(p,3);
            labels1 = {strcat("p=", pVal)};
            
            xPos = annotsCoord(zi,1);
            yPos = annotsCoord(zi,2);
            ht = text(xPos, yPos, labels1,'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','bold','BackgroundColor','black', 'Color', 'yellow','FontSize',12);
            set(ht,'Rotation',45);
%             if(p >= pTh)
%                 ht.Color = 'red';
%             end
        end

        ax = gca;
        ax.XAxis.FontSize = 12;
        xtickangle(45);
        set(findobj(gcf,'type','axes'),'FontWeight','Bold');
    end
    titleStr = {biomarkerName; "Activity Differences Between PSO Groups"};    
    sgtitle(titleStr,'FontSize',16);
    
    normStr = 'nonNormalized';
    if not(strcmp(params.normalization, ''))
        normStr = params.normalization;
    end
    
    %%
%     hold on;
%     for pi = 1:length(positions)
%         ci = round(pi/2);
%         oi =  2-mod(pi, 2);
%         xPosts = positions(pi)-0.1:0.025:positions(pi)+0.1;
%         yPosts = repmat(centroids(oi,ci), 1, length(xPosts));
%         symbol = '*';
%         if oi == 2
%             symbol = '*';
%         end
%         %p = plot(xPosts, yPosts, symbol,'Color', colors{ci}, 'LineWidth', 2); hold on;
%         p = plot(xPosts, yPosts, symbol,'Color', 'k', 'LineWidth', 2); hold on;
% 
%         p.MarkerFaceColor = colors{ci};
%         p.MarkerSize = 3;
%     end
    
    legend('Positive PSO', 'Negative PSO','Location','northeast');

    
    boxchartPath = strcat(params.outcomePredictionAnalysisPath, 'BoxplotAnalysis\'); mkdir(boxchartPath)
    boxchartFN = strcat(boxchartPath, strcat('BoxplotAnalysis_', params.feature), '_', normStr);
    saveas(gcf,boxchartFN);
    boxchartFN = strcat(boxchartFN, '.jpg');
    saveas(gcf,boxchartFN);
    close();
    
end