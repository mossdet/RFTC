function clstrIdx = clusterVector_flexK(vec)
    thirdQuartChannSum = sum(vec > prctile(vec, 75));

    [clstrIdx,C] = kmeans(vec, 2,'Distance','cityblock', 'Replicates', 1000, 'MaxIter',10000);
    [maxC, maxC_Idx] = max(C);
    clstrIdx(clstrIdx ~= maxC_Idx) = 0;
    clstrIdx(clstrIdx == maxC_Idx) = 1;
    
    if sum(clstrIdx) > thirdQuartChannSum
        [clstrIdx,C] = kmeans(vec, 3,'Distance','cityblock', 'Replicates', 1000, 'MaxIter',10000);
        [maxC, maxC_Idx] = max(C);
        clstrIdx(clstrIdx ~= maxC_Idx) = 0;
        clstrIdx(clstrIdx == maxC_Idx) = 1;
    end
    
%     if sum(clstrIdx) > thirdQuartChannSum
%         [clstrIdx,C] = kmeans(vec, 4,'Distance','cityblock', 'Replicates', 1000, 'MaxIter',10000);
%         [maxC, maxC_Idx] = max(C);
%         clstrIdx(clstrIdx ~= maxC_Idx) = 0;
%         clstrIdx(clstrIdx == maxC_Idx) = 1;
%     end
    
    if sum(clstrIdx) > thirdQuartChannSum
        clstrIdx = vec > prctile(vec, 75);
    end
    
end