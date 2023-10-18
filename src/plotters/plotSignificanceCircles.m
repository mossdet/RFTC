function plotSignificanceCircles(bm, resultsTableBiomarkerP)
    featuresList = {'Occ.Rate', 'Max.Ampl.','Power'};
    zonesNames = {'rftcSite', 'rftcStructure', 'rftcConnected', 'highEI', 'rftcLobe', 'rftcHemisphere'};
    nrZones = length(zonesNames);
    nrFeatures= length(featuresList);
    p = 0.0005;%0.05 / ((nrZones*nrFeatures)*3);
    resultsTableBiomarkerP(1,:) = [];
    resultsTableBiomarkerP(:,1) = [];
    resultsTableBiomarkerP(:,1) = [];
    
    t = tiledlayout(nrFeatures, nrZones,"TileSpacing","tight");
    for fi = 1:nrFeatures
        for zi = 1:nrZones
            featureName = featuresList{fi};
            zoneName =  zonesNames{zi};
            nexttile
            pos = [0.5 0.5 2 2];
            if resultsTableBiomarkerP{fi,zi} >= p
                rectangle('Position',pos,'Curvature',[1 1], 'FaceColor', [0.6353 0.0784 0.1843])
            else
                rectangle('Position',pos,'Curvature',[1 1], 'FaceColor', [0.4667 0.6745 0.1882])
            end
            %axis equal
            set(gca,'XTick',[])
            set(gca,'YTick',[])
            set(gca,'FontSize',24, 'FontWeight','Bold')
            if fi == 1
                title(zoneName)
            end
            if zi == 1
                ylabel(featureName)
            end
            %axis off
        end
    end
    set(gcf, 'Position', get(0, 'Screensize'), 'color','w');
    title(t, bm, 'FontSize',60, 'FontWeight','Bold')
    
    dim = [.9 .68 .3 .3];
    str = strcat("p < ", num2str(p,'%.4f'));
    annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',18, 'FontWeight','Bold')

end