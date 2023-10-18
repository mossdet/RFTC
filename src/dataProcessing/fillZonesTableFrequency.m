function allZonesTable = fillZonesTableFrequency(hfoDetections, allZonesTable, chIdx, nrMins)
    %%IES
    selDets = hfoDetections.mark == 3;
    avgOccRate = sum(selDets)/nrMins;
    avgPow = mean(hfoDetections.maxPower(selDets));
    avgSpectralPeak = mean(hfoDetections.maxSpectralPeak(selDets));
    allZonesTable.iaVals(chIdx) = avgSpectralPeak;

    %% HFO
    selDets = hfoDetections.mark == 1 | hfoDetections.mark == 2;
    avgOccRate = sum(selDets)/nrMins;
    avgPow = mean(hfoDetections.maxPower(selDets));
    avgSpectralPeak = mean(hfoDetections.maxSpectralPeak(selDets));
    allZonesTable.allHFOVals(chIdx) = avgSpectralPeak;

    %%iesHFO
    selDets = hfoDetections.mark == 4 | hfoDetections.mark == 5;
    avgOccRate = sum(selDets)/nrMins;
    avgPow = mean(hfoDetections.maxPower(selDets));
    avgSpectralPeak = mean(hfoDetections.maxSpectralPeak(selDets));
    allZonesTable.iesHFOVals(chIdx) = avgSpectralPeak;

    %%isolHFO
    selDets = hfoDetections.mark == 6 | hfoDetections.mark == 7;
    avgOccRate = sum(selDets)/nrMins;
    avgPow = mean(hfoDetections.maxPower(selDets));
    avgSpectralPeak = mean(hfoDetections.maxSpectralPeak(selDets));
    allZonesTable.isolHFOVals(chIdx) = avgSpectralPeak;

    %% Ripples
    %%All Ripples
    selDets = hfoDetections.mark == 1;
    avgOccRate = sum(selDets)/nrMins;
    avgPow = mean(hfoDetections.maxPower(selDets));
    avgSpectralPeak = mean(hfoDetections.maxSpectralPeak(selDets));
    allZonesTable.allRippleVals(chIdx) = avgSpectralPeak;

    %%iesRipples
    selDets = hfoDetections.mark == 4;
    avgOccRate = sum(selDets)/nrMins;
    avgPow = mean(hfoDetections.maxPower(selDets));
    avgSpectralPeak = mean(hfoDetections.maxSpectralPeak(selDets));
    allZonesTable.iesRippleVals(chIdx) = avgSpectralPeak;

    %%isoRipples
    selDets = hfoDetections.mark == 6;
    avgOccRate = sum(selDets)/nrMins;
    avgPow = mean(hfoDetections.maxPower(selDets));
    avgSpectralPeak = mean(hfoDetections.maxSpectralPeak(selDets));
    allZonesTable.isolRippleVals(chIdx) = avgSpectralPeak;

    %% FastRipples
    selDets = hfoDetections.mark == 2;
    avgOccRate = sum(selDets)/nrMins;
    avgPow = mean(hfoDetections.maxPower(selDets));
    avgSpectralPeak = mean(hfoDetections.maxSpectralPeak(selDets));
    allZonesTable.allFR_Vals(chIdx) = avgSpectralPeak;

    %%iesFR
    selDets = hfoDetections.mark == 5;
    avgOccRate = sum(selDets)/nrMins;
    avgPow = mean(hfoDetections.maxPower(selDets));
    avgSpectralPeak = mean(hfoDetections.maxSpectralPeak(selDets));
    allZonesTable.iesFR_Vals(chIdx) = avgSpectralPeak;

    %%isolHFO
    selDets = hfoDetections.mark == 7;
    avgOccRate = sum(selDets)/nrMins;
    avgPow = mean(hfoDetections.maxPower(selDets));
    avgSpectralPeak = mean(hfoDetections.maxSpectralPeak(selDets));
    allZonesTable.isolFR_Vals(chIdx) = avgSpectralPeak;

end