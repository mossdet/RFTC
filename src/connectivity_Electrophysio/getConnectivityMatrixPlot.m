function getConnectivityMatrixPlot(paths, patName, avgCorrelResults, channelNames, rftcFlags, outcomeVal)
   
    pointSize = 50;
    plotIdx = 1;
    h = figure;
    colormap jet;"jet";
    
    xx = [];
    yy = [];
    zz = [];

    for ri = 1:size(avgCorrelResults,1)
        for ci = 1:size(avgCorrelResults, 2)
            xx = cat(1, xx, ri);
            yy = cat(1, yy, ci);
            zz = cat(1, zz, avgCorrelResults(ri,ci));
        end
    end
    
    scatter(xx, yy, pointSize, zz, 'filled', 's');
    ax=gca;
    ax.FontSize = 8;
    xticks(1:length(channelNames));
    xticklabels(channelNames);
    xtickangle(90)
    
    yticks(1:length(channelNames));
    yticklabels(channelNames);
    
    rftcChannels = find(rftcFlags);
    ax = gca;
    for chi = 1:length(rftcChannels)
        rftcChIdx = rftcChannels(chi);
        ax.XTickLabel{rftcChIdx} = ['\color{red}' ax.XTickLabel{rftcChIdx}];
        ax.YTickLabel{rftcChIdx} = ['\color{red}' ax.YTickLabel{rftcChIdx}];
    end

    clblabel = colorbar; clblabel.Label.String = strcat('Correlation');
    xlim([min(xx) max(xx)])
    ylim([min(yy) max(yy)])
    zlim([min(zz) max(zz)])
    
    set(gca,'Color','k');
    grid on;
    grid minor;
    sgtitle({patName; strcat('Improvement:', num2str(outcomeVal), '%')}, 'Interpreter','none') 
    set(gcf, 'Position', get(0, 'Screensize'), 'color','w');

    figuresPath = strcat(paths.workspacePath, 'ConnectivityPlots\');mkdir(figuresPath);
    figFileName =  strcat(figuresPath, patName);
    hgexport(gcf, figFileName, hgexport('factorystyle'), 'Format', 'jpeg');
    close()

end