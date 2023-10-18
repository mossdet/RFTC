function [rftcCorr, eiCorr, rftcConnCorr] = plotBiomarker_ZoneCorrelation(params, groupTablePre)
    close all;
        
    %% Electrophysio Zones
    rftcCorr = [];
    eiCorr = [];
    rftcConnCorr = [];
    for bmi = 1:length(params.biomarkersList)
        params.biomarker = params.biomarkersList{bmi};
        corrVal = getSigCorr(groupTablePre.(params.biomarker), groupTablePre.rftcVals, 'Type','Spearman');
        rftcCorr = [rftcCorr corrVal];
        corrVal = getSigCorr(groupTablePre.(params.biomarker), groupTablePre.eiVals, 'Type','Spearman');
        eiCorr = [eiCorr corrVal];
        new_rftcElectroPhysioConnect = groupTablePre.rftcElectroPhysioConnect;
        new_rftcElectroPhysioConnect(groupTablePre.rftcVals>0) = max(new_rftcElectroPhysioConnect);
        corrVal = getSigCorr(groupTablePre.(params.biomarker), new_rftcElectroPhysioConnect, 'Type','Spearman');
        rftcConnCorr = [rftcConnCorr corrVal];
    end
    
    %% Anatomical Zones
    rftcParcel = [];
    rftcLobe = [];
    rftcHemisphere = [];
    rftcParcel_Channs = double(sameRFTC_BrainParcelChannsPerPat(groupTablePre));
    rftcLobe_Channs = double(sameRFTC_BrainLobeChannsPerPat(groupTablePre));
    rftcHemisphere_Channs = double(sameRFTC_HemisphereChannsPerPat(groupTablePre));
    for bmi = 1:length(params.biomarkersList)
        params.biomarker = params.biomarkersList{bmi};
        corrVal = getSigCorr(groupTablePre.(params.biomarker), rftcParcel_Channs, 'Type','Spearman');
        rftcParcel = [rftcParcel corrVal]; 
        corrVal = getSigCorr(groupTablePre.(params.biomarker), rftcLobe_Channs, 'Type','Spearman');
        rftcLobe = [rftcLobe corrVal];
        corrVal = getSigCorr(groupTablePre.(params.biomarker), rftcHemisphere_Channs, 'Type','Spearman');
        rftcHemisphere = [rftcHemisphere corrVal];
    end

         
    %% Prepare data for plot   
    colors = {[0 0.4470 0.7410];...	
    [0.8500 0.3250 0.0980];...
    [0.9290 0.6940 0.1250];...
    [0.4940 0.1840 0.5560];...
    [0.4660 0.6740 0.1880];...
    [0.3010 0.7450 0.9330];...
    [0.6350 0.0780 0.1840]};

    colors = {[0 0.4470 0.7410];...	
                [0.8500 0.3250 0.0980];...
                [0.9290 0.6940 0.1250];...
                [0.4940 0.1840 0.5560];...
                [0.4660 0.6740 0.1880];...
                [0.6350 0.0780 0.1840];...
                [0.3010 0.7450 0.9330]};

    labels = strrep(params.biomarkersList, 'Vals', '');
    labels = strrep(labels, '_', '');
    labels = strrep(labels, 'ia', 'IES');
    
    dataEP = [rftcCorr; rftcConnCorr; eiCorr];
    dataAN = [rftcParcel; rftcLobe; rftcHemisphere];
       
    
    figure(1)
    set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w');
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); 
    
    %% Plot Electrophysio
    nexttile
    axesLimits = zeros(2, length(params.biomarkersList));
    axesLimits(1,:) = min(dataEP,[],'all');
    axesLimits(2,:) = max(dataEP,[],'all');
    axesDirectionStr = 'normal';
    if abs(min(min(dataEP))) > abs(max(max(dataEP))) 
        axesDirectionStr = 'reverse';
    end
    
    spider_plot_R2019b(dataEP,...
                        'AxesLabels', labels,...
                        'AxesInterval', 5,...
                        'AxesPrecision', 2, ...
                        'AxesLimits', axesLimits,...
						'AxesLabelsOffset', 0.1,...
						'AxesDisplay', 'one',...
						'AxesLabelsEdge', 'none',...
						'AxesFontSize', 12,...
                        'AxesDirection', axesDirectionStr,...
                        'FillOption', {'on', 'on', 'on'},...
                        'FillTransparency', [0.4, 0.3, 0.2],...
                        'MarkerSize', 100,...
                        'AxesFontSize', 14,...
                        'LabelFontSize', 18);
    legend('rftcChannels', 'rftcConnected', 'HighEI', 'Location', 'southoutside', 'FontSize', 20);
    title('Electrophysiologic Channel Labels','FontSize',20);
    
    %% Plot Anatomical
    nexttile
    axesLimits = zeros(2, length(params.biomarkersList));
    axesLimits(1,:) = min(dataAN,[],'all');
    axesLimits(2,:) = max(dataAN,[],'all');
    axesDirectionStr = 'normal';
    if abs(min(min(dataAN))) > abs(max(max(dataAN))) 
        axesDirectionStr = 'reverse';
    end
    
    spider_plot_R2019b(dataAN,...
                        'AxesLabels', labels,...
                        'AxesInterval', 5,...
                        'AxesPrecision', 2, ...
                        'AxesLimits', axesLimits,...
						'AxesLabelsOffset', 0.1,...
						'AxesDisplay', 'one',...
						'AxesLabelsEdge', 'none',...
						'AxesFontSize', 12,...
                        'AxesDirection', axesDirectionStr,...
                        'FillOption', {'on', 'on', 'on'},...
                        'FillTransparency', [0.4, 0.3, 0.2],...
                        'MarkerSize', 100,...
                        'AxesFontSize', 14,...
                        'LabelFontSize', 18);
    legend('rftcStructure', 'rftcLobe', 'rftcHemisphere', 'Location', 'southoutside', 'FontSize', 20);
    title('Anatomic Channel Labels','FontSize',20);

    titleStr = strrep(params.feature, 'OccRate', 'Occurrence Rate');
    sgtitle(titleStr,'FontSize',24,'FontWeight','bold');
    
    %% Save Figure
    normStr = 'nonNormalized';
    if not(strcmp(params.normalization, ''))
        normStr = params.normalization;
    end

    epanData = [dataEP' dataAN'];
    allData = [epanData(3,:);epanData(2,:);epanData(1,:)];
    
    boxchartPath = strcat(params.analysisPerZonesTablePath, 'Biomarker_Zone_Correlations\'); mkdir(boxchartPath)    
    boxchartFN = strcat(boxchartPath, 'Biomarker_Zone_Correlations', '_', params.patsSel, '_', params.feature, '_', normStr);
    saveas(gcf,boxchartFN);
    boxchartFN = strcat(boxchartFN, '.jpg');
    saveas(gcf,boxchartFN);
    close();
end

function rho = getSigCorr(vecA, vecB, typeSet, corrType)
    [rho,pval] = corr(vecA, vecB, typeSet,corrType);
    if not(pval < 0.01)
        rho = 0;
    end
end