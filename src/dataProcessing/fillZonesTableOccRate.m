function allZonesTable = fillZonesTableOccRate(hfoDetections, allZonesTable, chIdx, nrMins)
    %%IES
    selDets = hfoDetections.mark == 3;
    avgOccRate = sum(selDets)/nrMins;
    avgPow = mean(hfoDetections.maxPower(selDets));
    avgSpectralPeak = mean(hfoDetections.maxSpectralPeak(selDets));
    allZonesTable.iaVals(chIdx) = avgOccRate;

    %% HFO
    selDets = hfoDetections.mark == 1 | hfoDetections.mark == 2;
    avgOccRate = sum(selDets)/nrMins;
    avgPow = mean(hfoDetections.maxPower(selDets));
    avgSpectralPeak = mean(hfoDetections.maxSpectralPeak(selDets));
    allZonesTable.allHFOVals(chIdx) = avgOccRate;

    %%iesHFO
    selDets = hfoDetections.mark == 4 | hfoDetections.mark == 5;
    avgOccRate = sum(selDets)/nrMins;
    avgPow = mean(hfoDetections.maxPower(selDets));
    avgSpectralPeak = mean(hfoDetections.maxSpectralPeak(selDets));
    allZonesTable.iesHFOVals(chIdx) = avgOccRate;

    %%isolHFO
    selDets = hfoDetections.mark == 6 | hfoDetections.mark == 7;
    avgOccRate = sum(selDets)/nrMins;
    avgPow = mean(hfoDetections.maxPower(selDets));
    avgSpectralPeak = mean(hfoDetections.maxSpectralPeak(selDets));
    allZonesTable.isolHFOVals(chIdx) = avgOccRate;

    %% Ripples
    %%All Ripples
    selDets = hfoDetections.mark == 1;
    avgOccRate = sum(selDets)/nrMins;
    avgPow = mean(hfoDetections.maxPower(selDets));
    avgSpectralPeak = mean(hfoDetections.maxSpectralPeak(selDets));
    allZonesTable.allRippleVals(chIdx) = avgOccRate;

    %%iesRipples
    selDets = hfoDetections.mark == 4;
    avgOccRate = sum(selDets)/nrMins;
    avgPow = mean(hfoDetections.maxPower(selDets));
    avgSpectralPeak = mean(hfoDetections.maxSpectralPeak(selDets));
    allZonesTable.iesRippleVals(chIdx) = avgOccRate;

    %%isoRipples
    selDets = hfoDetections.mark == 6;
    avgOccRate = sum(selDets)/nrMins;
    avgPow = mean(hfoDetections.maxPower(selDets));
    avgSpectralPeak = mean(hfoDetections.maxSpectralPeak(selDets));
    allZonesTable.isolRippleVals(chIdx) = avgOccRate;

    %% FastRipples
    selDets = hfoDetections.mark == 2;
    avgOccRate = sum(selDets)/nrMins;
    avgPow = mean(hfoDetections.maxPower(selDets));
    avgSpectralPeak = mean(hfoDetections.maxSpectralPeak(selDets));
    allZonesTable.allFR_Vals(chIdx) = avgOccRate;

    %%iesFR
    selDets = hfoDetections.mark == 5;
    avgOccRate = sum(selDets)/nrMins;
    avgPow = mean(hfoDetections.maxPower(selDets));
    avgSpectralPeak = mean(hfoDetections.maxSpectralPeak(selDets));
    allZonesTable.iesFR_Vals(chIdx) = avgOccRate;

    %%isolHFO
    selDets = hfoDetections.mark == 7;
    avgOccRate = sum(selDets)/nrMins;
    avgPow = mean(hfoDetections.maxPower(selDets));
    avgSpectralPeak = mean(hfoDetections.maxSpectralPeak(selDets));
    allZonesTable.isolFR_Vals(chIdx) = avgOccRate;

end