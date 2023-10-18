function [sensitivityP, specificityP, accuracyP, foneP, mccP] = prospectivePrediction(patNames, perPatient_Outcome, zonesDiffData, zonesSel, maxCentroids, distanceMt)

    clusterData = zonesDiffData(:,zonesSel);
    [perPatient_Outcome_Sorted, newOrder] = sort(perPatient_Outcome, 'descend');
    clusterData = clusterData(newOrder, :);
    outcomeLabels = perPatient_Outcome_Sorted >= 50;
    [~,idx_test] = pdist2(maxCentroids, clusterData, distanceMt, 'Smallest', 1); % cosine squaredeuclidean 
    idx_test = idx_test - max(idx_test) + 1;
    clusterIdxT = idx_test';
    [sensitivityP, specificityP, accuracyP, foneP, mccP] = getMetrics(clusterIdxT, outcomeLabels); 
    
    
    [sensitivityP, specificityP, accuracyP, foneP, mccP] = manualPorspectivePrediction(patNames, perPatient_Outcome, zonesDiffData, zonesSel, maxCentroids)
end

function [sensitivity, specificity,accuracy, fone, mcc] = getMetrics(vecA, vecB)
    vecA = logical(vecA);
    vecB = logical(vecB);

    tp = sum(vecA & vecB);
    tn = sum(not(vecA) & not(vecB));
    fp = sum(vecA & not(vecB));
    fn = sum(not(vecA) & vecB);
    
    sensitivity = tp/(tp+fn);
    specificity = tn/(tn+fp);
    accuracy =  (tp+tn)/(tp+tn+fp+fn);
    
    fone = (2*tp) / (2*tp + fp + fn);
    
    mccA = (tp * tn) - (fp * fn);
    mccB = ((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));
    mcc = mccA / sqrt(mccB);
end

function [sensitivityP, specificityP, accuracyP, foneP, mccP] = manualPorspectivePrediction(patNames, perPatient_Outcome, zonesDiffData, zonesSel, maxCentroids)
    clusterData = zonesDiffData(:,zonesSel);
    outcomeLabels = perPatient_Outcome >= 50;

    posDiffcos = zeros(size(clusterData,1),1);
    negDiffcos = zeros(size(clusterData,1),1);
    for pi = 1:size(clusterData,1)
        patObs = clusterData(pi,:);
        posCent = maxCentroids(1,:);
        negCent = maxCentroids(2,:);

        posDiffcos(pi) = 1 - (patObs*posCent')/sqrt((patObs*patObs')*(posCent*posCent'));
        negDiffcos(pi) = 1 - (patObs*negCent')/sqrt((patObs*patObs')*(negCent*negCent'));
    end
    posNegCos = [posDiffcos';negDiffcos'];
    [posNegCosMin, posNegCosMinIdx] = min(posNegCos,[],1);
    posNegCosMinIdx = posNegCosMinIdx - max(posNegCosMinIdx) + 1;
    manClusterIdxT = posNegCosMinIdx';
    [sensitivityP, specificityP, accuracyP, foneP, mccP] = getMetrics(manClusterIdxT, outcomeLabels);
end