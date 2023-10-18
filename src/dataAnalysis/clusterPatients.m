function [maxSensitivity, maxSpecificity, maxAccuracy, maxFone, maxMCC, maxCentroids, distanceMt, zonesSel] = clusterPatients(patNames, perPatient_Outcome, zonesDiffData, zonesSel) 

    % 1. RFTC
    % 2. rftcConnc
    % 3. highEI
    % 4. rftcStruct
    % 5. rftcLobe
    % 6. rftcHemisphere   
    repeatCnt = 0;
    maxSensitivity = 0;
    maxSpecificity = 0;
    maxAccuracy = 0;
    maxFone = 0;
    maxMCC = 0;
    maxCentroids = [];
    distanceMt = 'cosine'; % sqeuclidean cityblock cosine correlation

    clusterData = zonesDiffData(:,zonesSel);
    clusterData(isnan(clusterData))=0;
    outcomeLabels = perPatient_Outcome >= 50;
   
    pi = randperm(length(perPatient_Outcome));
    reordClusterData = clusterData(pi, :);
    reordOutcomeLabels = outcomeLabels(pi);
    reordPatNames = patNames(pi);


    while repeatCnt < 1000
    
        [cidx, centroids, sumd, D] = kmeans(reordClusterData , 2, 'Start','cluster', 'Distance', distanceMt,'MaxIter',1000); %
        
        cidx = cidx - max(cidx) + 1;
        clusterIdx = cidx;
        invertClusters = false(1,1);
        if sum(clusterIdx == 0) < sum(clusterIdx == 1)
            sel = clusterIdx == 0;
            clusterIdx(sel) = true(1,1);
            clusterIdx(not(sel)) = false(1,1);
            invertClusters = true(1,1);
        end        
        
        [sensitivity, specificity,accuracy, fone, mcc] = getMetrics(clusterIdx, reordOutcomeLabels);
        if (fone > maxFone) && not(invertClusters)
            maxSensitivity = sensitivity;
            maxSpecificity = specificity;
            maxAccuracy = accuracy;
            maxFone = fone;
            maxMCC = mcc;
            maxCentroids = centroids;
        end

        repeatCnt = repeatCnt+1;
    end
    [maxSensitivity, maxSpecificity, maxAccuracy]    
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
