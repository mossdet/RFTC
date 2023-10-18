function maxRftcConnect = getMax_RFTC_ElectroPhysioConnect(paths, patName, channelName, freqBandConn)

    electroPhysioConnectPath = strcat(paths.workspacePath, 'ConnectivityMatricesAllBands\');
    electroPhysioConnectFN = strcat(electroPhysioConnectPath, patName);
    load(electroPhysioConnectFN, 'deltaConn', 'thetaConn', 'alphaConn', 'betaConn', 'gammaConn',...
                            'highGammaConn', 'rippleConn', 'frConn', 'maxCorrelAllBands', 'meanCorrelAllBands');    
    
    bandConn = [];
    if strcmp(freqBandConn, 'delta')
        bandConn = deltaConn;
    elseif strcmp(freqBandConn, 'theta')
        bandConn = thetaConn;
    elseif strcmp(freqBandConn, 'alpha')
        bandConn = alphaConn;
    elseif strcmp(freqBandConn, 'beta')
        bandConn = betaConn;
    elseif strcmp(freqBandConn, 'gamma')
        bandConn = gammaConn;
    elseif strcmp(freqBandConn, 'highGamma')
        bandConn = highGammaConn;
    elseif strcmp(freqBandConn, 'ripple')
        bandConn = rippleConn;
    elseif strcmp(freqBandConn, 'fr')
        bandConn = frConn;
    elseif strcmp(freqBandConn, 'maxAllBands')
        bandConn = maxCorrelAllBands;
    elseif strcmp(freqBandConn, 'meanAllBands')
        bandConn = meanCorrelAllBands;
    else
        "Error"
    end
    
    avgCorr = bandConn.corr;
    channelNames = bandConn.channelNames;
    rftcVals = bandConn.rftcVals;  
    
    channSel = find(ismember(channelNames, channelName));
    if length(channSel) > 1
       error('getMax_RFTC_ElectroPhysioConnect, more than two channels found') 
    end
    rftcChannsSel = find(rftcVals > 0);
    maxRftcConnect = max(avgCorr(channSel, rftcChannsSel),[],'all');
    if isempty(maxRftcConnect)
        maxRftcConnect = 0;
    end
end