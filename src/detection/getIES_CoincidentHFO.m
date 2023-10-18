function hfoIES_Detections = getIES_CoincidentHFO(hfoDetections)
    
    nrDetections = length(hfoDetections.mark);
    hfoIES_Detections = hfoDetections;
    
    for fdi = 1:nrDetections    %iterate through HFO
        iesCoincidence = 0;
        fdType = hfoDetections.mark(fdi);
        fdStart = hfoDetections.startTime(fdi);
        fdEnd = hfoDetections.endTime(fdi);
        fdDuration = fdEnd - fdStart;
        
        if(fdType == 3)
            continue;
        end
        
        for sdi = 1:nrDetections    %iterate through IES
            sdType = hfoDetections.mark(sdi);
            if (fdi == sdi || sdType ~= 3)
                continue;
            end
            sdStart = hfoDetections.startTime(sdi);
            sdEnd = hfoDetections.endTime(sdi);
            sdDuration = sdEnd - sdStart;
            overlapTime = getEventsOverlap(fdStart, fdEnd, sdStart, sdEnd);
            overlapPerc = 100*(overlapTime / fdDuration);
            if (100 * (overlapTime / sdDuration) > overlapPerc)
                overlapPerc = 100 * (overlapTime / sdDuration);
            end
            if overlapPerc > 50.0
                iesCoincidence = 1;
                break;
            end
        end
        
        newDetection = getIndexedDetections(hfoDetections, fdi);

        if fdType == 1
            if iesCoincidence > 0
                newDetection.mark = 4;
                hfoIES_Detections = concatenateDetections(hfoIES_Detections, newDetection);
            else
                newDetection.mark = 6;
                hfoIES_Detections = concatenateDetections(hfoIES_Detections, newDetection);
            end
        elseif fdType == 2 
            if iesCoincidence > 0
                newDetection.mark = 5;
                hfoIES_Detections = concatenateDetections(hfoIES_Detections, newDetection);
            else
                newDetection.mark = 7;
                hfoIES_Detections = concatenateDetections(hfoIES_Detections, newDetection);
            end
        end
        
        %     - All Ripples     (1) -> any Ripple
        %     - All FR          (2) -> any FR
        %     - All IES         (3) -> any IES

        %     - IES_Ripples     (4) -> any Ripple coinciding with a IES
        %     - IES_FR          (5) -> any FR coinciding with a IES
        %     - isolRipples     (6) -> any Ripple not coinciding with IES
        %     - isolFR          (7) -> any FR not coinciding with IES
    end
end

function overlapTime = getEventsOverlap(feStart, feEnd, seStart, seEnd)
    overlapTime = 0;
    feDuration = feEnd - feStart;
    seDuration = seEnd - seStart;
    
    if feStart <= seStart && feEnd >= seEnd % first fully encompassing second
        overlapTime = seDuration;        
    elseif seStart <= feStart && seEnd >= feEnd % second fully encompassing first
        overlapTime = feDuration;        
    elseif (feStart <= seStart && feEnd >= seStart && feEnd <= seEnd) %last part of first inside second
        overlapTime = feEnd - seStart;
    elseif (seStart <= feStart && seEnd >= feStart && seEnd <= feEnd) %last part of second inside first
        overlapTime = seEnd - feStart;
    else
        overlapTime = 0;
    end
end