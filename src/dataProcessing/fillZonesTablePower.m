function allZonesTable = fillZonesTablePower(hfoDetections, allZonesTable, chIdx, nrMins)
    %%IES
    selDets = hfoDetections.mark == 3;
    avgOccRate = sum(selDets)/nrMins;
    avgPow = mean(hfoDetections.maxPower(selDets));
    avgSpectralPeak = mean(hfoDetections.maxSpectralPeak(selDets));
    allZonesTable.iaVals(chIdx) = avgPow;

    %% HFO
    selDets = hfoDetections.mark == 1 | hfoDetections.mark == 2;
    avgOccRate = sum(selDets)/nrMins;
    avgPow = mean(hfoDetections.maxPower(selDets));
    avgSpectralPeak = mean(hfoDetections.maxSpectralPeak(selDets));
    allZonesTable.allHFOVals(chIdx) = avgPow;

    %%iesHFO
    selDets = hfoDetections.mark == 4 | hfoDetections.mark == 5;
    avgOccRate = sum(selDets)/nrMins;
    avgPow = mean(hfoDetections.maxPower(selDets));
    avgSpectralPeak = mean(hfoDetections.maxSpectralPeak(selDets));
    allZonesTable.iesHFOVals(chIdx) = avgPow;

    %%isolHFO
    selDets = hfoDetections.mark == 6 | hfoDetections.mark == 7;
    avgOccRate = sum(selDets)/nrMins;
    avgPow = mean(hfoDetections.maxPower(selDets));
    avgSpectralPeak = mean(hfoDetections.maxSpectralPeak(selDets));
    allZonesTable.isolHFOVals(chIdx) = avgPow;

    %% Ripples
    %%All Ripples
    selDets = hfoDetections.mark == 1;
    avgOccRate = sum(selDets)/nrMins;
    avgPow = mean(hfoDetections.maxPower(selDets));
    avgSpectralPeak = mean(hfoDetections.maxSpectralPeak(selDets));
    allZonesTable.allRippleVals(chIdx) = avgPow;

    %%iesRipples
    selDets = hfoDetections.mark == 4;
    avgOccRate = sum(selDets)/nrMins;
    avgPow = mean(hfoDetections.maxPower(selDets));
    avgSpectralPeak = mean(hfoDetections.maxSpectralPeak(selDets));
    allZonesTable.iesRippleVals(chIdx) = avgPow;

    %%isoRipples
    selDets = hfoDetections.mark == 6;
    avgOccRate = sum(selDets)/nrMins;
    avgPow = mean(hfoDetections.maxPower(selDets));
    avgSpectralPeak = mean(hfoDetections.maxSpectralPeak(selDets));
    allZonesTable.isolRippleVals(chIdx) = avgPow;

    %% FastRipples
    selDets = hfoDetections.mark == 2;
    avgOccRate = sum(selDets)/nrMins;
    avgPow = mean(hfoDetections.maxPower(selDets));
    avgSpectralPeak = mean(hfoDetections.maxSpectralPeak(selDets));
    allZonesTable.allFR_Vals(chIdx) = avgPow;

    %%iesFR
    selDets = hfoDetections.mark == 5;
    avgOccRate = sum(selDets)/nrMins;
    avgPow = mean(hfoDetections.maxPower(selDets));
    avgSpectralPeak = mean(hfoDetections.maxSpectralPeak(selDets));
    allZonesTable.iesFR_Vals(chIdx) = avgPow;

    %%isolHFO
    selDets = hfoDetections.mark == 7;
    avgOccRate = sum(selDets)/nrMins;
    avgPow = mean(hfoDetections.maxPower(selDets));
    avgSpectralPeak = mean(hfoDetections.maxSpectralPeak(selDets));
    allZonesTable.isolFR_Vals(chIdx) = avgPow;

end