function plotPooledPatientsAveragedHFO(allPatAverageEOI)
    close all;
    hfoSel = ismember(allPatAverageEOI(:,4), 'HFO');
    ieshfoSel = ismember(allPatAverageEOI(:,4), 'iesHFO');

    allHFO = allPatAverageEOI(hfoSel, :);
    iesHFO = allPatAverageEOI(ieshfoSel, :);
    
    goodOutcome_HFOSel = cell2mat(allHFO(:,2))>= 90;
    badOutcome_HFOSel = not(goodOutcome_HFOSel);
    goodOutcomeHFO = allHFO(goodOutcome_HFOSel,:);
    badOutcomeHFO = allHFO(badOutcome_HFOSel,:);
    
    goodOutcomePostSel = contains(lower(goodOutcomeHFO(:,1)), 'post');
    goodOutcomePreSel = not(goodOutcomePostSel);
    goodOutcomeHFOPre = goodOutcomeHFO(goodOutcomePreSel,:);
    goodOutcomeHFOPost = goodOutcomeHFO(goodOutcomePostSel,:);
    
    badOutcomePostSel = contains(lower(badOutcomeHFO(:,1)), 'post');
    badOutcomePreSel = not(badOutcomePostSel);
    badOutcomeHFOPre = badOutcomeHFO(badOutcomePreSel,:);
    badOutcomeHFOPost = badOutcomeHFO(badOutcomePostSel,:);

    %%
    goodOutcome_iesHFOSel = cell2mat(iesHFO(:,2))>= 90;
    badOutcome_iesHFOSel = not(goodOutcome_iesHFOSel);
    goodOutcomeiesHFO = iesHFO(goodOutcome_iesHFOSel,:);
    badOutcomeiesHFO = iesHFO(badOutcome_iesHFOSel,:);
    
    goodOutcomePostSel = contains(lower(goodOutcomeiesHFO(:,1)), 'post');
    goodOutcomePreSel = not(goodOutcomePostSel);
    goodOutcomeiesHFOPre = goodOutcomeiesHFO(goodOutcomePreSel,:);
    goodOutcomeiesHFOPost = goodOutcomeiesHFO(goodOutcomePostSel,:);

    badOutcomePostSel = contains(lower(badOutcomeiesHFO(:,1)), 'post');
    badOutcomePreSel = not(badOutcomePostSel);
    badOutcomeiesHFOPre = badOutcomeiesHFO(badOutcomePreSel,:);
    badOutcomeiesHFOPost = badOutcomeiesHFO(badOutcomePostSel,:);

    %%
    [timeVec, avgGoodOutcomeHFOPre] = getAverageEOI(goodOutcomeHFOPre);
    [timeVec, avgGoodOutcomeHFOPost] = getAverageEOI(goodOutcomeHFOPost);
    [timeVec, avgGoodOutcomeiesHFOPre] = getAverageEOI(goodOutcomeiesHFOPre);
    [timeVec, avgGoodOutcomeiesHFOPost] = getAverageEOI(goodOutcomeiesHFOPost);

    %%
    [timeVec, avgBadOutcomeHFOPre] = getAverageEOI(badOutcomeHFOPre);
    [timeVec, avgBadOutcomeHFOPost] = getAverageEOI(badOutcomeHFOPost);
    [timeVec, avgBadOutcomeiesHFOPre] = getAverageEOI(badOutcomeiesHFOPre);
    [timeVec, avgBadOutcomeiesHFOPost] = getAverageEOI(badOutcomeiesHFOPost);

    %%
    fs = 2048;
    subplot(2,4,1)
    plotThisAvgHFO(fs, timeVec, avgGoodOutcomeHFOPre, "Pre-RFTC HFO")
    subplot(2,4,2)
    plotThisAvgHFO(fs, timeVec, avgGoodOutcomeHFOPost, "Post-RFTC HFO")

    subplot(2,4,3)
    plotThisAvgHFO(fs, timeVec, avgBadOutcomeHFOPre, "Pre-RFTC HFO")
    subplot(2,4,4)
    plotThisAvgHFO(fs, timeVec, avgBadOutcomeHFOPost, "Post-RFTC HFO")

    subplot(2,4,5)
    plotThisAvgHFO(fs, timeVec, avgGoodOutcomeiesHFOPre, "Pre-RFTC iesHFO")
    subplot(2,4,6)
    plotThisAvgHFO(fs, timeVec, avgGoodOutcomeiesHFOPost, "Post-RFTC iesHFO")

    subplot(2,4,7)
    plotThisAvgHFO(fs, timeVec, avgBadOutcomeiesHFOPre, "Pre-RFTC iesHFO")
    subplot(2,4,8)
    plotThisAvgHFO(fs, timeVec, avgBadOutcomeiesHFOPost, "Post-RFTC iesHFO")
end

function [timeVec, avgVec] = getAverageEOI(groupVecs)
    avgVec = zeros(1,2049);
    for pi = 1:size(groupVecs,1)
        hfoVec = groupVecs{pi,5};
        if length(hfoVec)~=2049
            continue;
        end
        avgVec = avgVec + hfoVec;
        timeVec = groupVecs{pi,6};
    end
    avgVec = avgVec/size(groupVecs,1);
end

function plotThisAvgHFO(fs, timeVec, avg_HFO, eoiStr)
    ylimVec = [min(avg_HFO) max(avg_HFO)];
    plot(timeVec, avg_HFO, 'k'); hold on;
    xlim([min(timeVec), max(timeVec)]);
    ylim(ylimVec(1,:));
    xlabel("Time (s)", 'FontSize', 24);
    ylabel("Amplitude (uV)", 'FontSize', 24);
    title(eoiStr, 'FontSize', 28);
    ax = gca;
    ax.XAxis.FontSize = 18;
    ax.XAxis.FontWeight = 'bold';
    ax.YAxis.FontSize = 18;
    ax.YAxis.FontWeight = 'bold';

    analSigSel = int32(length(avg_HFO)/2-fs*0.05:length(avg_HFO)/2+fs*0.05);
    analSig = avg_HFO(analSigSel);
    analSigTime = timeVec(analSigSel);
    plot(analSigTime, analSig, '-k')

    maxAmpl = abs(max(analSig)-min(analSig));
    power = sum(analSig.*analSig)/length(analSig);
    annotStr = strcat("Max.Amplitude = ", num2str(maxAmpl,'%.2f'), "\muV");
    xPos = min(xlim)+(max(xlim)-min(xlim))*0.05;
    yPos = min(ylim)+(max(ylim)-min(ylim))*0.98;
    text(xPos, yPos, annotStr, 'Horiz','left', 'Vert','top', 'BackgroundColor', 'w', 'Color', 'black', 'FontWeight', 'bold', 'FontSize', 16);

    annotStr = strcat("Power = ", num2str(power,'%.2f'), "\muV^2");
    xPos = min(xlim)+(max(xlim)-min(xlim))*0.05;
    yPos = min(ylim)+(max(ylim)-min(ylim))*0.88;
    text(xPos, yPos, annotStr, 'Horiz','left', 'Vert','top', 'BackgroundColor', 'w', 'Color', 'black', 'FontWeight', 'bold', 'FontSize', 16);
end