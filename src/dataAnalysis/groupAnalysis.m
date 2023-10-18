clc; clear all; close all;
paths = getFilesPaths();
%files  = getAllFiles();%getAllFiles, getPreFiles, getPostFiles
preFiles  = getPreFiles();%getAllFiles, getPreFiles, getPostFiles
postFiles  = getPostFiles();%getAllFiles, getPreFiles, getPostFiles

selDetector = 'MOSSDET'; 'Delphos'; 'MOSSDET'; 'MOSSDET_Depurated';
thList = 0:30;

if strcmp(selDetector, 'Delphos')
    thList = 0;
end

for th = thList
    for fileIdx = 1:size(preFiles,1)
        if th > 0
            selDetector = 'MOSSDET_Depurated';
        end

        tablesFilePath = strcat(paths.workspacePath, 'ChannelCharacterizationTables\', selDetector, '\');

        if th == 0
            tablesFilePath = tablesFilePath;
        else 
            tablesFilePath = strcat(tablesFilePath, 'Th', num2str(th), '\');
        end

        prePatName = preFiles{fileIdx}; prePatName = prePatName(1:length(prePatName)-4);
        postPatName = postFiles{fileIdx}; postPatName = postPatName(1:length(postPatName)-4);

        preTableFN = strcat(tablesFilePath, prePatName, '_', 'ChannelCharacterization_', selDetector, '.xls');
        postTableFN = strcat(tablesFilePath, postPatName, '_', 'ChannelCharacterization_', selDetector, '.xls');

        preTable = readtable(preTableFN);
        postTable = readtable(postTableFN);

        patCodeSeparators = strfind(prePatName, '_');
        patCode = prePatName(patCodeSeparators(2):patCodeSeparators(3));patCode = strrep(patCode, '_', '');
        outcomeTableFilename = 'F:\ForschungsProjekte\RFTC\MATLAB\DetectHFO\OtherData\Lachner_DetectedFiles_List.xlsx';
        outcomeTable = readtable(outcomeTableFilename, 'Sheet', 'MicromedFiles(.TRC)');
        outcomeStr = num2str(outcomeTable.Post_RFTCImprovement___(find(ismember(outcomeTable.Code,patCode))));

        pointSize = 100;
        % iaVals
        % allHFOVals, iesHFOVals, isolHFOVals
        % allRippleVals, iesRippleVals, isolRippleVals
        % allFR_Vals, iesFR_Vals, isolFR_Vals
        plotIdx = 1;
        h = figure;
        colormap parula;
        for evType = 1:3
            preSelRFTC = preTable.rftcVals == 1;
            xPre = (preTable.eiVals-min(preTable.eiVals));xPre= xPre/max(xPre);
            yPre = preTable.iaVals;

            postSelRFTC = postTable.rftcVals == 1;
            xPost = (postTable.eiVals-min(postTable.eiVals));xPost= xPost/max(xPost);
            yPost = postTable.iaVals;

            switch evType
                case 1
                    cPre = preTable.allHFOVals;
                    cPost = postTable.allHFOVals;
                    titleString = 'HFO Activity';
                case 2
                    cPre = preTable.allRippleVals;
                    cPost = postTable.allRippleVals;
                    titleString = 'Ripple Activity';
                case 3
                    cPre = preTable.allFR_Vals;
                    cPost = postTable.allFR_Vals;
                    titleString = 'FR Activity';
            end

            subplot(3, 2, plotIdx)
            scatter(xPre(preSelRFTC), yPre(preSelRFTC), pointSize+30, 'red', 'LineWidth', 5); hold on;
            scatter(xPre, yPre, pointSize, cPre, 'filled');
            title(strcat('Pre-RFTC-', titleString));
            %yticks(yPre); yticklabels(preTable.channelLabels)
            clblabel = colorbar; clblabel.Label.String = strcat(titleString);
            ylim([min([yPre; yPost]) max([yPre; yPost])])
            caxis([min([cPre; cPost]) max([cPre; cPost])])
            xlabel('EI');
            ylabel('IES Occ.Rate');
            plotIdx = plotIdx+1;

            subplot(3, 2, plotIdx)
            scatter(xPost(postSelRFTC), yPost(postSelRFTC), pointSize+30, 'red', 'LineWidth', 5); hold on;
            scatter(xPost, yPost, pointSize, cPost, 'filled');
            title(strcat('Post-RFTC-', titleString));
            %yticks(yPost); yticklabels(postTable.channelLabels)
            clblabel = colorbar; clblabel.Label.String = strcat(titleString);
            ylim([min([yPre; yPost]) max([yPre; yPost])])
            caxis([min([cPre; cPost]) max([cPre; cPost])])
            xlabel('EI')
            ylabel('IES Occ.Rate');
            plotIdx = plotIdx+1;

        end
        superTitle = {prePatName(1:strfind(prePatName, '_Inter')-1); strcat(outcomeStr, '% Improvement')};
        sgtitle(superTitle, 'Interpreter','none') 

        set(gcf, 'Position', get(0, 'Screensize'), 'color','w');
        figuresPath = strcat(paths.workspacePath, 'AnalysisPlots\', selDetector, '\');mkdir(figuresPath);
        figFileName =  strcat(figuresPath, patCode);
        %savefig(h, figFileName, 'compact');
        hgexport(gcf, figFileName, hgexport('factorystyle'), 'Format', 'jpeg');
        close()
        %plot3D_ScatterPlot();
    end
end

function [featureNames, allPats, allPatsNorm] = getPhaseAvgFeats(patientsInfo, dataLocation, excludeChannels, periPdgmMin)

    allEvntsAllChnnsAllPatsAvgFeatsPre = [];
    allEvntsAllChnnsAllPatsAvgFeatsIntra = [];
    allEvntsAllChnnsAllPatsAvgFeatsPost = [];
    allEvntsAllChnnsAllPatsAvgFeatsPreNorm = [];
    allEvntsAllChnnsAllPatsAvgFeatsIntraNorm = [];
    allEvntsAllChnnsAllPatsAvgFeatsPostNorm = [];
    featureNames = [];

    for patIdx = 1:size(patientsInfo,2)
        patInfo = patientsInfo{patIdx};
        patName = patInfo{1};
        hdr = dlpReadHeader(dataLocation, patName);    
        whiteMatterChannels = getWhiteMatterChannels(patName);
        validChannels = getValidBipolarChannels(hdr.labels, cat(2, excludeChannels,whiteMatterChannels));
        validChannels = sort(validChannels);

        nrChanns = length(validChannels);
        allEvntsAllChnnsAvgFeatsPre = [];
        allEvntsAllChnnsAvgFeatsIntra = [];
        allEvntsAllChnnsAvgFeatsPost = [];

        for chi = 1:nrChanns
            chName = validChannels{chi};
            hcChann = strcmp(getElectrodeName(chName), 'HAR') || strcmp(getElectrodeName(chName), 'HAL');
            hcParcel = getChannBrainParcel(hdr.patName, chName) == 1;
            %if hcChann && hcParcel
                [prePhaseHFO, intraPhaseHFO, postPhaseHFO, preDurM, intraDurM, postDurM] = readCharacterizedHFO(periPdgmMin, hdr,chName);
%             else
%                 prePhaseHFO = [];
%                 intraPhaseHFO = [];
%                 postPhaseHFO = [];
%                 preDurM = 0;
%                 intraDurM = 0;
%                 postDurM = 0;
%             end
            chName
            hdr.patName
            
            featureNames = fieldnames(prePhaseHFO{1}.features);
            
            allEventsAvgFeatsPre = getAllEventsAvgMatrix(hdr, patIdx, chi, chName, featureNames, prePhaseHFO, preDurM);
            allEventsAvgFeatsIntra = getAllEventsAvgMatrix(hdr, patIdx, chi, chName, featureNames, intraPhaseHFO, intraDurM);
            allEventsAvgFeatsPost = getAllEventsAvgMatrix(hdr, patIdx, chi, chName, featureNames, postPhaseHFO, postDurM);

            allEvntsAllChnnsAvgFeatsPre = cat(1, allEvntsAllChnnsAvgFeatsPre, allEventsAvgFeatsPre);
            allEvntsAllChnnsAvgFeatsIntra = cat(1, allEvntsAllChnnsAvgFeatsIntra, allEventsAvgFeatsIntra);
            allEvntsAllChnnsAvgFeatsPost = cat(1, allEvntsAllChnnsAvgFeatsPost, allEventsAvgFeatsPost);
        end
        
        if not(isempty(allEvntsAllChnnsAvgFeatsPre)) && not(isempty(allEvntsAllChnnsAvgFeatsIntra)) && not(isempty(allEvntsAllChnnsAvgFeatsPost))
            %normalize
            [allEvntsAllChnnsAvgFeatsPreNorm, allEvntsAllChnnsAvgFeatsIntraNorm, allEvntsAllChnnsAvgFeatsPostNorm] = normalizeFeatsAcrossChannsPerPatient(allEvntsAllChnnsAvgFeatsPre, allEvntsAllChnnsAvgFeatsIntra, allEvntsAllChnnsAvgFeatsPost);
    %         checkNormalization(allEvntsAllChnnsAvgFeatsPre, allEvntsAllChnnsAvgFeatsPreNorm, validChannels, featureNames, patName)
    %         checkNormalization(allEvntsAllChnnsAvgFeatsIntra, allEvntsAllChnnsAvgFeatsIntraNorm, validChannels, featureNames, patName)
    %         checkNormalization(allEvntsAllChnnsAvgFeatsPost, allEvntsAllChnnsAvgFeatsPostNorm, validChannels, featureNames, patName)

            %plotPhasesHeatmapAllChanns(allEvntsAllChnnsAvgFeatsPre, allEvntsAllChnnsAvgFeatsIntra, allEvntsAllChnnsAvgFeatsPost, validChannels, featureNames, patName, periPdgmMin);

            allEvntsAllChnnsAllPatsAvgFeatsPre = cat(1, allEvntsAllChnnsAllPatsAvgFeatsPre, allEvntsAllChnnsAvgFeatsPre);
            allEvntsAllChnnsAllPatsAvgFeatsIntra = cat(1, allEvntsAllChnnsAllPatsAvgFeatsIntra, allEvntsAllChnnsAvgFeatsIntra);
            allEvntsAllChnnsAllPatsAvgFeatsPost = cat(1, allEvntsAllChnnsAllPatsAvgFeatsPost, allEvntsAllChnnsAvgFeatsPost);

            allEvntsAllChnnsAllPatsAvgFeatsPreNorm = cat(1, allEvntsAllChnnsAllPatsAvgFeatsPreNorm, allEvntsAllChnnsAvgFeatsPreNorm);
            allEvntsAllChnnsAllPatsAvgFeatsIntraNorm = cat(1, allEvntsAllChnnsAllPatsAvgFeatsIntraNorm, allEvntsAllChnnsAvgFeatsIntraNorm);
            allEvntsAllChnnsAllPatsAvgFeatsPostNorm = cat(1, allEvntsAllChnnsAllPatsAvgFeatsPostNorm, allEvntsAllChnnsAvgFeatsPostNorm);
        end

    end
    
    allPats.pre = allEvntsAllChnnsAllPatsAvgFeatsPre;
    allPats.intra = allEvntsAllChnnsAllPatsAvgFeatsIntra;
    allPats.post = allEvntsAllChnnsAllPatsAvgFeatsPost;
    
    allPatsNorm.pre = allEvntsAllChnnsAllPatsAvgFeatsPreNorm;
    allPatsNorm.intra = allEvntsAllChnnsAllPatsAvgFeatsIntraNorm;
    allPatsNorm.post = allEvntsAllChnnsAllPatsAvgFeatsPostNorm;
end

function plotPhasesHeatmapAllChanns(allEvntsAllChnnsAvgFeatsPre, allEvntsAllChnnsAvgFeatsIntra, allEvntsAllChnnsAvgFeatsPost, channLabels, featureNames, patName, periPdgmMin)
    dotsize = 150;
    featuresToPlot = [5 6 7 19];
    zonesCompare = [1:6]; 
    featuresToPlot = 5;%:length(featureNames);
    zonesCompare = [1];
    
    featEventCodeIdx = 3;
    zoneCompareIdx = 4;
    for zi = 1:length(zonesCompare)
        zoneToCompare = zonesCompare(zi);
        nrFeats = length(featuresToPlot);
        for fi = 1:nrFeats
            featCode = featuresToPlot(fi);
            featName = featureNames{featCode};
            eventsToPlot = [1 4 6 3 2 5 7 3];
            eventsToPlot = [1 4 6 3];
            subplotIdx = 1;
            %ha = tight_subplot(2, 4,[.08 .04],[.1 .08],[.05 .01]);
            ha = tight_subplot(1, 4,[.08 .06],[.1 .08],[.07 .007]);
            for ei = 1:length(eventsToPlot)
                eventCode = eventsToPlot(ei);
                selPre = allEvntsAllChnnsAvgFeatsPre(:,featEventCodeIdx)==eventCode;
                sigToPlotPre = allEvntsAllChnnsAvgFeatsPre(selPre,featCode);
                nrChanns = length(sigToPlotPre);
                xPre = zeros(1,nrChanns)+1; yPre = 1:1:nrChanns; zPre = sigToPlotPre';
                chSelPre = allEvntsAllChnnsAvgFeatsPre(selPre,2);
                allMinSelParcel = selPre & allEvntsAllChnnsAvgFeatsPre(:,zoneCompareIdx) ~= zoneToCompare & allEvntsAllChnnsAvgFeatsPre(:,zoneCompareIdx) ~= 0;
                hcSel = selPre & allEvntsAllChnnsAvgFeatsPre(:,zoneCompareIdx) == zoneToCompare;
                avgAMHPre = mean(allEvntsAllChnnsAvgFeatsPre(allMinSelParcel,featCode));
                avgHCPre = mean(allEvntsAllChnnsAvgFeatsPre(hcSel,featCode));
                strPre = strcat('(',  num2str(mean(avgHCPre),'%.1f'), '\_vs\_',  num2str(mean(avgAMHPre),'%.1f'), ')', ' PRE');


                selIntra = allEvntsAllChnnsAvgFeatsIntra(:,featEventCodeIdx)==eventCode;
                sigToPlotIntra = allEvntsAllChnnsAvgFeatsIntra(selIntra,featCode);
                nrChanns = length(sigToPlotIntra);
                xIntra = zeros(1,nrChanns)+2; yIntra = 1:1:nrChanns; zIntra = sigToPlotIntra';
                chSelIntra = allEvntsAllChnnsAvgFeatsIntra(selIntra,2);
                allMinSelParcel = selIntra & allEvntsAllChnnsAvgFeatsIntra(:,zoneCompareIdx) ~= zoneToCompare & allEvntsAllChnnsAvgFeatsIntra(:,zoneCompareIdx) ~= 0;
                hcSel = selIntra & allEvntsAllChnnsAvgFeatsIntra(:,zoneCompareIdx) == zoneToCompare;
                avgAMHPIntra = mean(allEvntsAllChnnsAvgFeatsIntra(allMinSelParcel,featCode));
                avgHCIntra = mean(allEvntsAllChnnsAvgFeatsIntra(hcSel,featCode));
                strIntra = strcat('(',  num2str(mean(avgHCIntra),'%.1f'), '\_vs\_',  num2str(mean(avgAMHPIntra),'%.1f'), ')', ' INTRA');


                selPost = allEvntsAllChnnsAvgFeatsPost(:,featEventCodeIdx)==eventCode;
                sigToPlotPost = allEvntsAllChnnsAvgFeatsPost(selPost,featCode);
                nrChanns = length(sigToPlotPost);
                xPost = zeros(1,nrChanns)+3; yPost = 1:1:nrChanns; zPost = sigToPlotPost';
                chSelPost = allEvntsAllChnnsAvgFeatsPost(selPost,2);
                allMinSelParcel = selPost & allEvntsAllChnnsAvgFeatsPost(:,zoneCompareIdx) ~= zoneToCompare & allEvntsAllChnnsAvgFeatsPost(:,zoneCompareIdx) ~= 0;
                hcSel = selPost & allEvntsAllChnnsAvgFeatsPost(:,zoneCompareIdx) == zoneToCompare;
                avgAMHPost = mean(allEvntsAllChnnsAvgFeatsPost(allMinSelParcel,featCode));
                avgHCPost = mean(allEvntsAllChnnsAvgFeatsPost(hcSel,featCode));
                strPost = strcat('(',  num2str(mean(avgHCPost),'%.1f'), '\_vs\_',  num2str(mean(avgAMHPost),'%.1f'), ')', ' POST');
                [maxNrChanns idx]= max([length(xPre) length(xIntra) length(xPost)]);

                if length(xPre) <maxNrChanns
                    xPre = zeros(1, maxNrChanns)+1;
                    yPre = 1:1:maxNrChanns;
                    zPre = zeros(1, maxNrChanns);
                    strPre = [];
                elseif  length(xIntra) <maxNrChanns
                    xIntra = zeros(1, maxNrChanns)+2;
                    yIntra = 1:1:maxNrChanns;
                    zIntra = zeros(1, maxNrChanns);
                    strIntra = [];
                elseif  length(xPost) <maxNrChanns
                    xPost = zeros(1, maxNrChanns)+3;
                    yPost = 1:1:maxNrChanns;
                    zPost = zeros(1, maxNrChanns);
                    strPost = [];
                end

                xAllPhases = [xPre xIntra xPost];
                yAllPhases = [yPre yIntra yPost];
                zAllPhases = [zPre zIntra zPost];

                %subplot(2,4,subplotIdx)
                axes(ha(subplotIdx));
                s = scatter3(xAllPhases, yAllPhases, zAllPhases, dotsize, zAllPhases, 's', 'filled','HandleVisibility','off');
                view([0 90]);
                s.MarkerFaceAlpha = 0.9;
                set(gca,'Color',  'k')
                xticks(0:1:4);
                xlim([0 4]);
                set(gca,'TickDir','out');
                yticks(1:1:nrChanns);
                ylim([1 nrChanns]);
                xticklabels({''; strPre; strIntra; strPost; ''});
                xtickangle(20) 
                yticklabels(channLabels(chSelPost));
                title({patName; strcat(featName, '-', getEventName(eventCode))})
                grid on;
                
                ax = gca;
                channelsToColor = unique(allEvntsAllChnnsAvgFeatsPost(hcSel,2));
                for i = 1:length(channelsToColor)
                    yLabelIdx = channelsToColor(i);
                    ax.YTickLabel{yLabelIdx} = ['\color{red}' ax.YTickLabel{yLabelIdx}];     
                end

                set(gca,'GridColor', 'w', 'GridAlpha', 1)

                h  = colorbar;

                subplotIdx = subplotIdx+1;
            end

            
            compareType = strcat(getParcelName(zoneToCompare), ' vs Rest');
            set(gcf, 'Position', get(0, 'Screensize'));
            sgtitle(compareType); 
            channPlotDir = strcat('..\Heatmaps_', num2str(periPdgmMin),'m\');%, patName, '\');
            mkdir(channPlotDir)
            compareType = strcat(getParcelName(zoneToCompare), '_vs_Rest');
            figOneFileName = strcat(channPlotDir, patName, '_',compareType, '_',featName, '_HeatmapAllPhases_', num2str(periPdgmMin),'m');
            hgexport(gcf, figOneFileName, hgexport('factorystyle'), 'Format', 'jpeg');
            close();
        end
    end
end

function [allEvntsAllChnnsAvgFeatsPreNorm, allEvntsAllChnnsAvgFeatsIntraNorm, allEvntsAllChnnsAvgFeatsPostNorm] = normalizeFeatsAcrossChannsPerPatient(allEvntsAllChnnsAvgFeatsPre, allEvntsAllChnnsAvgFeatsIntra, allEvntsAllChnnsAvgFeatsPost)
    %normalize
    allEvntsAllChnnsAvgFeatsPreNorm = allEvntsAllChnnsAvgFeatsPre; allEvntsAllChnnsAvgFeatsIntraNorm = allEvntsAllChnnsAvgFeatsIntra; allEvntsAllChnnsAvgFeatsPostNorm = allEvntsAllChnnsAvgFeatsPost;
    allEvntsAllChnnsAvgFeatsAllPhases = cat(1, allEvntsAllChnnsAvgFeatsPreNorm, allEvntsAllChnnsAvgFeatsIntraNorm, allEvntsAllChnnsAvgFeatsPostNorm);
    
    featNormStartIdx = 5;
    featEventCodeIdx = 3;
    eventCodes = unique(allEvntsAllChnnsAvgFeatsAllPhases(:,featEventCodeIdx))';
    for et = 1:7
        selEv = allEvntsAllChnnsAvgFeatsAllPhases(:,featEventCodeIdx) == et;
        avgFeatVals = mean(allEvntsAllChnnsAvgFeatsAllPhases(selEv, :),1);
        stdDevFeatVals = std(allEvntsAllChnnsAvgFeatsAllPhases(selEv, :),1);

        %Pre
        selEvPre = allEvntsAllChnnsAvgFeatsPreNorm(:,featEventCodeIdx) == et;
        allEvntsAllChnnsAvgFeatsPreNorm(selEvPre,featNormStartIdx:end) = (allEvntsAllChnnsAvgFeatsPreNorm(selEvPre,featNormStartIdx:end)-avgFeatVals(featNormStartIdx:end))./stdDevFeatVals(featNormStartIdx:end);
        if isempty(selEvPre)
            stop =1;
        end


        %Intra
        selEvIntra = allEvntsAllChnnsAvgFeatsIntraNorm(:,featEventCodeIdx) == et;
        allEvntsAllChnnsAvgFeatsIntraNorm(selEvIntra,featNormStartIdx:end) = (allEvntsAllChnnsAvgFeatsIntraNorm(selEvIntra,featNormStartIdx:end)-avgFeatVals(featNormStartIdx:end))./stdDevFeatVals(featNormStartIdx:end);
        if isempty(selEvIntra)
            stop =1;
        end

        %Post
        selEvPost = allEvntsAllChnnsAvgFeatsPostNorm(:,featEventCodeIdx) == et;
        allEvntsAllChnnsAvgFeatsPostNorm(selEvPost,featNormStartIdx:end) = (allEvntsAllChnnsAvgFeatsPostNorm(selEvPost,featNormStartIdx:end)-avgFeatVals(featNormStartIdx:end))./stdDevFeatVals(featNormStartIdx:end);
        if isempty(selEvPost)
            stop =1;
        end

    end
end

function allEventsAvgFeats = getAllEventsAvgMatrix(hdr, patIdx, chi, chName, featureNames, phaseHFO, durM)
    nrFeats = length(featureNames);
    nrHFO = length(phaseHFO);
    chArea = getChannBrainParcel(hdr.patName, chName);
    %chArea = getChannArea(chName);
    allEventsMatrix = double(zeros(nrHFO,nrFeats));
    for ei = 1:nrHFO
        hfo = phaseHFO{ei};
        allEventsMatrix(ei,1) = patIdx;
        allEventsMatrix(ei,2) = chi;
        allEventsMatrix(ei,3) = hfo.type;
        allEventsMatrix(ei,4) = chArea;
        allEventsMatrix(ei,5) = 1;%count
        allEventsMatrix(ei,6) = double(hfo.se-hfo.ss)/hdr.fs;%duration
        for fi = 7:nrFeats
            factor = 1;
            if fi >= 7 && fi <= 12
                factor = 1000*1000;
            end
            allEventsMatrix(ei,fi) = hfo.features.(featureNames{fi})*factor; 
        end
    end
    eventCodes = 1:7;
    nrEventTypes = length(eventCodes);
    allEventsAvgFeats = zeros(nrEventTypes, nrFeats);
    for eci = 1:nrEventTypes
        eventCode = eventCodes(eci);
        sel = allEventsMatrix(:,3) == eventCode;
        if sum(sel) == 0
            stop = 1;
            allEventsAvgFeats(eci,:) = zeros(1, nrFeats);
            allEventsAvgFeats(eci,1) = patIdx;
            allEventsAvgFeats(eci,2) = chi;
            allEventsAvgFeats(eci,3) = eventCode;
            allEventsAvgFeats(eci,4) = chArea;
            allEventsAvgFeats(eci,5) = 0;
        else
            allEventsAvgFeats(eci,:) = mean(allEventsMatrix(sel,:),1);
            allEventsAvgFeats(eci,5) = sum(sel)/durM;            
        end
    end           
end

function checkNormalization(allEvntsAllChnnsAvgFeats, allEvntsAllChnnsAvgFeatsNorm, channLabels, featureNames, patName)
    nrFeats = length(featureNames);
    featEventCodeIdx = 3;
    for fi = 5:nrFeats
        featName = featureNames{fi};
        for et = 1:7
            subplot(2,1,1)
            sel = allEvntsAllChnnsAvgFeats(:,featEventCodeIdx)==et;
            chSel = allEvntsAllChnnsAvgFeats(sel,2);
            sigToPlot = allEvntsAllChnnsAvgFeats(sel,fi);
            nrChanns = length(sigToPlot);
            plot(sigToPlot)
            xticks(1:1:nrChanns)
            xlim([1 nrChanns])
            xticklabels(channLabels(chSel));
            xtickangle(60);
            title({patName; strcat(featName, '-', getEventName(et))})
            
            subplot(2,1,2)
            sel = allEvntsAllChnnsAvgFeatsNorm(:,featEventCodeIdx)==et;
            chSel = allEvntsAllChnnsAvgFeatsNorm(sel,2);
            sigToPlot = allEvntsAllChnnsAvgFeatsNorm(sel,fi);
            nrChanns = length(sigToPlot);
            plot(sigToPlot)
            xticks(1:1:nrChanns)
            xlim([1 nrChanns])
            xticklabels(channLabels(chSel));
            xtickangle(60);
            title({patName; strcat('Normalized-', featName, '-', getEventName(et))})
            close();
        end
    end
end

function eventName = getEventName(eventCode)
    if eventCode == 1
        eventName = 'Ripples';
    elseif eventCode == 2
        eventName = 'FR';
    elseif eventCode == 3
        eventName = 'IES';
    elseif eventCode == 4
        eventName = 'IES-Ripples';
    elseif eventCode == 5
        eventName = 'IES-FR';
    elseif eventCode == 6
        eventName = 'isolRipples';
    elseif eventCode == 7
        eventName = 'isolFR';
    end
end

        
function hdr = dlpReadHeader(dataLocation, patName)
    %Surname	GivenName	Day	Month	Year	Hour	Min	Sec	SR	NrSamples	DoubleSize
    headerFilename = strcat(dataLocation, patName, '.bin');
    [Surname, GivenName, Day, Month, Year, Hour, Min, Sec, samplingRate, totalNrSamples, DoubleSize] =...
        textread(headerFilename, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t', 'headerlines', 1);
    channInfoFileName =  strcat(dataLocation, patName, '_montageChannels.txt');
    [ch_name, ch_nr] = textread(channInfoFileName, '%s\t%d', 'headerlines', 1);
    
    hdr.patName = patName;
    hdr.location = dataLocation;
    hdr.fs = samplingRate;
    hdr.nrSamples = totalNrSamples;
    hdr.labels = ch_name;
end

function [prePhaseHFO, intraPhaseHFO, postPhaseHFO, preDurM, intraDurM, postDurM] = readCharacterizedHFO(periPdgmMin, hdr, chName)
    filePath = strcat('..\HFO_Detections\', hdr.patName, '\', chName, '\');
    
    %filePath = strcat('..\HFO_Detections_1h_PrePost\', hdr.patName, '\', chName, '\');
    loadFN = strcat(filePath, chName, '_prePhaseHFO', num2str(periPdgmMin),'m');
    load(loadFN, 'prePhaseHFO', 'preBound');
    loadFN = strcat(filePath, chName, '_intraPhaseHFO', num2str(periPdgmMin),'m');
    load(loadFN, 'intraPhaseHFO', 'intraBound');
    loadFN = strcat(filePath, chName, '_postPhaseHFO', num2str(periPdgmMin),'m');
    load(loadFN, 'postPhaseHFO', 'postBound');

    preDurM = (preBound(2)-preBound(1))/(hdr.fs*60);
    intraDurM = (intraBound(2)-intraBound(1))/(hdr.fs*60);
    postDurM = (postBound(2)-postBound(1))/(hdr.fs*60);
end

function patientsInfo = getPatientsInfoAll()
    pi = 1;
    patientsInfo{pi} = {'2', {''}}; pi = pi+1;
    patientsInfo{pi} = {'3', {''}}; pi = pi+1;
    patientsInfo{pi} = {'4', {''}}; pi = pi+1;
    patientsInfo{pi} = {'4_sleep', {''}}; pi = pi+1;
    patientsInfo{pi} = {'4b', {''}}; pi = pi+1;
    patientsInfo{pi} = {'7', {''}}; pi = pi+1;
    patientsInfo{pi} = {'7_sleep', {''}}; pi = pi+1;
    patientsInfo{pi} = {'7b', {''}}; pi = pi+1;
    patientsInfo{pi} = {'8', {''}}; pi = pi+1;
    patientsInfo{pi} = {'8_sleep', {''}}; pi = pi+1;
    patientsInfo{pi} = {'8b', {''}}; pi = pi+1;
    patientsInfo{pi} = {'11', {''}}; pi = pi+1;
    patientsInfo{pi} = {'12', {''}}; pi = pi+1;
    patientsInfo{pi} = {'12_sleep', {''}}; pi = pi+1;
    patientsInfo{pi} = {'12b', {''}}; pi = pi+1;
    patientsInfo{pi} = {'13', {''}}; pi = pi+1;
    patientsInfo{pi} = {'13_sleep', {''}}; pi = pi+1;
    patientsInfo{pi} = {'13b', {''}}; pi = pi+1;
    patientsInfo{pi} = {'14', {''}}; pi = pi+1;
    patientsInfo{pi} = {'14_sleep', {''}}; pi = pi+1;
    patientsInfo{pi} = {'14b', {''}}; pi = pi+1;
    patientsInfo{pi} = {'15', {''}}; pi = pi+1;
    patientsInfo{pi} = {'16', {''}}; pi = pi+1;
    patientsInfo{pi} = {'16_sleep', {''}}; pi = pi+1;
    patientsInfo{pi} = {'16b', {''}}; pi = pi+1;
    patientsInfo{pi} = {'18', {''}}; pi = pi+1;
    patientsInfo{pi} = {'19', {''}}; pi = pi+1;
    patientsInfo{pi} = {'19_sleep', {''}}; pi = pi+1;
    patientsInfo{pi} = {'19b', {''}}; pi = pi+1;
    patientsInfo{pi} = {'20', {''}}; pi = pi+1;
    patientsInfo{pi} = {'21', {''}}; pi = pi+1;
    patientsInfo{pi} = {'22', {''}}; pi = pi+1;
    patientsInfo{pi} = {'23', {''}}; pi = pi+1;
    patientsInfo{pi} = {'24', {''}}; pi = pi+1;
    patientsInfo{pi} = {'25', {''}}; pi = pi+1;
    patientsInfo{pi} = {'26', {''}}; pi = pi+1;
    patientsInfo{pi} = {'27', {''}}; pi = pi+1;
    patientsInfo{pi} = {'27_sleep', {''}}; pi = pi+1;
    patientsInfo{pi} = {'27b', {''}}; pi = pi+1;
    patientsInfo{pi} = {'28', {''}}; pi = pi+1;
    patientsInfo{pi} = {'28_sleep', {''}}; pi = pi+1;
    patientsInfo{pi} = {'28b', {''}}; pi = pi+1;
end

function patientsInfo = getPatientsInfoAwake()
    pi = 1;
    patientsInfo{pi} = {'2', {''}}; pi = pi+1;
    patientsInfo{pi} = {'3', {''}}; pi = pi+1;
    patientsInfo{pi} = {'4', {''}}; pi = pi+1;
    patientsInfo{pi} = {'4b', {''}}; pi = pi+1;
    patientsInfo{pi} = {'7', {''}}; pi = pi+1;
    patientsInfo{pi} = {'7b', {''}}; pi = pi+1;
    patientsInfo{pi} = {'8', {''}}; pi = pi+1;
    patientsInfo{pi} = {'8b', {''}}; pi = pi+1;
    patientsInfo{pi} = {'11', {''}}; pi = pi+1;
    patientsInfo{pi} = {'12', {''}}; pi = pi+1;
    patientsInfo{pi} = {'12b', {''}}; pi = pi+1;
    patientsInfo{pi} = {'13', {''}}; pi = pi+1;
    patientsInfo{pi} = {'13b', {''}}; pi = pi+1;
    patientsInfo{pi} = {'14', {''}}; pi = pi+1;
    patientsInfo{pi} = {'14b', {''}}; pi = pi+1;
    patientsInfo{pi} = {'15', {''}}; pi = pi+1;
    patientsInfo{pi} = {'16', {''}}; pi = pi+1;
    patientsInfo{pi} = {'16b', {''}}; pi = pi+1;
    patientsInfo{pi} = {'18', {''}}; pi = pi+1;
    patientsInfo{pi} = {'19', {''}}; pi = pi+1;
    patientsInfo{pi} = {'19b', {''}}; pi = pi+1;
    patientsInfo{pi} = {'20', {''}}; pi = pi+1;
    patientsInfo{pi} = {'21', {''}}; pi = pi+1;
    patientsInfo{pi} = {'22', {''}}; pi = pi+1;
    patientsInfo{pi} = {'23', {''}}; pi = pi+1;
    patientsInfo{pi} = {'24', {''}}; pi = pi+1;
    patientsInfo{pi} = {'25', {''}}; pi = pi+1;
    patientsInfo{pi} = {'26', {''}}; pi = pi+1;
    patientsInfo{pi} = {'27', {''}}; pi = pi+1;
    patientsInfo{pi} = {'27b', {''}}; pi = pi+1;
    patientsInfo{pi} = {'28', {''}}; pi = pi+1;
    patientsInfo{pi} = {'28b', {''}}; pi = pi+1;
end

%%
function patientsInfo = getAllSleepPatientsInfo()
    pi = 1;
    patientsInfo{pi} = {'4', {''}}; pi = pi+1;
    patientsInfo{pi} = {'4b', {''}}; pi = pi+1;
    patientsInfo{pi} = {'7', {''}}; pi = pi+1;
    patientsInfo{pi} = {'7b', {''}}; pi = pi+1;
    patientsInfo{pi} = {'8', {''}}; pi = pi+1;
%     patientsInfo{pi} = {'8b', {''}}; pi = pi+1;
%     patientsInfo{pi} = {'12', {''}}; pi = pi+1;
%     patientsInfo{pi} = {'12b', {''}}; pi = pi+1;
%     patientsInfo{pi} = {'13', {''}}; pi = pi+1;
%     patientsInfo{pi} = {'13b', {''}}; pi = pi+1;
%     patientsInfo{pi} = {'14', {''}}; pi = pi+1;
%     patientsInfo{pi} = {'14b', {''}}; pi = pi+1;
%     patientsInfo{pi} = {'16', {''}}; pi = pi+1;
%     patientsInfo{pi} = {'16b', {''}}; pi = pi+1;
%     patientsInfo{pi} = {'19', {''}}; pi = pi+1;
%     patientsInfo{pi} = {'19b', {''}}; pi = pi+1;
%     patientsInfo{pi} = {'27', {''}}; pi = pi+1;
%     patientsInfo{pi} = {'27b', {''}}; pi = pi+1;
%     patientsInfo{pi} = {'28', {''}}; pi = pi+1;
%     patientsInfo{pi} = {'28b', {''}}; pi = pi+1;
end

function patientsInfo = getScalpAndInvasiveSleepPatientsInfo()
    pi = 1;
    patientsInfo{pi} = {'4_sleep', {''}}; pi = pi+1;
    patientsInfo{pi} = {'7_sleep', {''}}; pi = pi+1;
    patientsInfo{pi} = {'8_sleep', {''}}; pi = pi+1;
    patientsInfo{pi} = {'13_sleep', {''}}; pi = pi+1;
    patientsInfo{pi} = {'16_sleep', {''}}; pi = pi+1;
    patientsInfo{pi} = {'19_sleep', {''}}; pi = pi+1;
end