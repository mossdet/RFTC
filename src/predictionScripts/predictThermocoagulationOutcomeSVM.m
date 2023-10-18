function kappaCum = predictThermocoagulationOutcomeSVM(featsVals, outcomeVals)
    nrObservations = length(outcomeVals);
    trainingPerc = 0.7;
    cvK = nrObservations;
    mccCum = 0;
    kappaCum = 0;

    rng(0,'twister');
    trainNrEvents = round(trainingPerc*nrObservations);
    trainingEvsIdxs = randi([1 nrObservations], cvK, trainNrEvents);
    allEvsIdxs = 1:nrObservations;
    for cvi = 1:cvK
        trainingEvs = featsVals(trainingEvsIdxs(cvi,:),:);
        trainingLabels = outcomeVals(trainingEvsIdxs(cvi,:));
        
        testEvsIdxs = allEvsIdxs;
        testEvsIdxs(trainingEvsIdxs(cvi,:)) = [];
        testEvs = featsVals(testEvsIdxs,:);
        testsLabels = outcomeVals(testEvsIdxs,:);

        % 'linear' | 'gaussian' | 'rbf' | 'polynomial'
        SVMModel = fitcsvm(trainingEvs, trainingLabels,'Standardize',true,'KernelFunction','linear',...
        'KernelScale','auto');
        [label,score] = predict(SVMModel, testEvs);
        mcc = getMCC(testsLabels, label);
        mccCum = mccCum + mcc;
        kappa = getKappa(testsLabels, label);
        kappaCum = kappaCum+kappa;
    end
    mccCum = mccCum/cvK;
    kappaCum = kappaCum/cvK;
end

function mcc = getMCC(vecA, vecB)
    vecA = logical(vecA);
    vecB = logical(vecB);

    tp = sum(vecA & vecB);
    tn = sum(not(vecA) & not(vecB));
    fp = sum(not(vecA) & vecB);
    fn = sum(vecA & not(vecB));
    
    mccA = (tp * tn) - (fp * fn);
    mccB = ((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));
    if mccB == 0
        mcc =  0;
    else
        mcc = mccA / sqrt(mccB);
    end
end

function kappa = getKappa(vecA, vecB)
    tp = 0; fp = 0; tn = 0; fn = 0;

    tp = tp + sum(vecA & vecB);
    fp = fp + sum(not(vecA) & vecB);
    tn = tn + sum(not(vecA) & not(vecB));
    fn = fn + sum(vecA & not(vecB));
    kappaA = 100 * (2 * (tp*tn - fn*fp));
    kappaB = 100 * ((tp+fp) * (fp+tn) + (tp+fn) * (fn+tn));
    if kappaB <= 0.01
        kappa = 0;
    else
        kappa = kappaA/kappaB;
    end
end