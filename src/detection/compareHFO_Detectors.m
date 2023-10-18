clc; clear all; close all;

paths = getFilesPaths();

files  = getAllFiles();%getPreFiles, getPostFiles
nrPats = size(files,1);

sz = [nrPats 13];
varTypes = {'string','double','double','double','double','double','double','double','double','double','double','double','double'};
varNames = {'Patient','eiVals','rftcVals','iaVals','allHFOVals','iesHFOVals','isolHFOVals','allRippleVals','iesRippleVals','isolRippleVals','allFR_Vals','iesFR_Vals','isolFR_Vals'};
rhoTable = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
pValTable = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);

thresholds = 0:17;
rhoTableSummary = table('Size',[length(thresholds) 13],'VariableTypes',varTypes,'VariableNames',varNames);
pValTableSummary = table('Size',[length(thresholds) 13],'VariableTypes',varTypes,'VariableNames',varNames);

for thi = 1:length(thresholds)

    detTh = thresholds(thi);
    
    for fi = 1:nrPats

        patName = files{fi}; 
        patName = patName(1:length(patName)-4);
        mossdetTablesFilePath = strcat(paths.workspacePath, 'ChannelCharacterizationTables\MOSSDET\');
        delphosTablesFilePath = strcat(paths.workspacePath, 'ChannelCharacterizationTables\Delphos\');

        if detTh > 0
            mossdetTablesFilePath = strcat(paths.workspacePath, 'ChannelCharacterizationTables\MOSSDET_Depurated', '\Th', num2str(detTh), '\');
        end

        mossdetTableFN = strcat(mossdetTablesFilePath, patName, '_', 'ChannelCharacterization_MOSSDET.xls');
        delphosTableFN = strcat(delphosTablesFilePath, patName, '_', 'ChannelCharacterization_Delphos.xls');

        mossdetTable = readtable(mossdetTableFN);
        delphosTable = readtable(delphosTableFN);

        feat = 'Patient';
        rhoTable.(feat)(fi) = patName;
        pValTable.(feat)(fi) = patName;

        feat = 'eiVals';
        [rho,pval] = corr(mossdetTable.(feat), delphosTable.(feat)); rhoTable.(feat)(fi) = rho; pValTable.(feat)(fi) = pval;
        feat = 'rftcVals';
        [rho,pval] = corr(mossdetTable.(feat), delphosTable.(feat)); rhoTable.(feat)(fi) = rho; pValTable.(feat)(fi) = pval;
        feat = 'iaVals';
        [rho,pval] = corr(mossdetTable.(feat), delphosTable.(feat)); rhoTable.(feat)(fi) = rho; pValTable.(feat)(fi) = pval;

        feat = 'allHFOVals';
        [rho,pval] = corr(mossdetTable.(feat), delphosTable.(feat)); rhoTable.(feat)(fi) = rho; pValTable.(feat)(fi) = pval;
        feat = 'iesHFOVals';
        [rho,pval] = corr(mossdetTable.(feat), delphosTable.(feat)); rhoTable.(feat)(fi) = rho; pValTable.(feat)(fi) = pval;
        feat = 'isolHFOVals';
        [rho,pval] = corr(mossdetTable.(feat), delphosTable.(feat)); rhoTable.(feat)(fi) = rho; pValTable.(feat)(fi) = pval;

        feat = 'allRippleVals';
        [rho,pval] = corr(mossdetTable.(feat), delphosTable.(feat)); rhoTable.(feat)(fi) = rho; pValTable.(feat)(fi) = pval;
        feat = 'iesRippleVals';
        [rho,pval] = corr(mossdetTable.(feat), delphosTable.(feat)); rhoTable.(feat)(fi) = rho; pValTable.(feat)(fi) = pval;
        feat = 'isolRippleVals';
        [rho,pval] = corr(mossdetTable.(feat), delphosTable.(feat)); rhoTable.(feat)(fi) = rho; pValTable.(feat)(fi) = pval;

        feat = 'allFR_Vals';
        [rho,pval] = corr(mossdetTable.(feat), delphosTable.(feat)); rhoTable.(feat)(fi) = rho; pValTable.(feat)(fi) = pval;
        feat = 'iesFR_Vals';
        [rho,pval] = corr(mossdetTable.(feat), delphosTable.(feat)); rhoTable.(feat)(fi) = rho; pValTable.(feat)(fi) = pval;
        feat = 'isolFR_Vals';
        [rho,pval] = corr(mossdetTable.(feat), delphosTable.(feat)); rhoTable.(feat)(fi) = rho; pValTable.(feat)(fi) = pval;

    end
    
    
    features = rhoTable.Properties.VariableNames;
    nrFeats = length(features);

    rhoTableSummary.('Patient')(thi) = num2str(detTh);
    pValTableSummary.('Patient')(thi) = num2str(detTh);

    for fi = 2:nrFeats
        featureName = features{fi};
        rhoFieldVals = rhoTable.(featureName)(:);
        pFieldVals = pValTable.(featureName)(:);
        rhoTableSummary.(featureName)(thi) = mean(rhoFieldVals);
        pValTableSummary.(featureName)(thi) = mean(pFieldVals);
    end

end