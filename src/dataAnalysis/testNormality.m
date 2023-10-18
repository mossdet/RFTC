function testNormality(groupTablePre)

    iaVals = [];
    allHFOVals = []; iesHFOVals = []; isolHFOVals = [];
    allRippleVals = []; iesRippleVals = []; isolRippleVals = [];
    allFR_Vals = []; iesFR_Vals = []; isolFR_Vals = [];
    
    plotDistro(groupTablePre.iaVals);
    plotDistro(groupTablePre.eiVals);
    iaVals = cat(1, iaVals, corr(groupTablePre.iaVals, groupTablePre.eiVals, 'Type','Spearman'));

    plotDistro(groupTablePre.iaVals);
    allHFOVals = cat(1, allHFOVals, corr(groupTablePre.allHFOVals, groupTablePre.eiVals, 'Type','Spearman'));
    iesHFOVals = cat(1, iesHFOVals, corr(groupTablePre.iesHFOVals, groupTablePre.eiVals, 'Type','Spearman'));
    isolHFOVals = cat(1, isolHFOVals, corr(groupTablePre.isolHFOVals, groupTablePre.eiVals, 'Type','Spearman'));

    allRippleVals = cat(1, allRippleVals, corr(groupTablePre.allRippleVals, groupTablePre.eiVals, 'Type','Spearman'));
    iesRippleVals = cat(1, iesRippleVals, corr(groupTablePre.iesRippleVals, groupTablePre.eiVals, 'Type','Spearman'));
    isolRippleVals = cat(1, isolRippleVals, corr(groupTablePre.isolRippleVals, groupTablePre.eiVals, 'Type','Spearman'));

    allFR_Vals = cat(1, allFR_Vals, corr(groupTablePre.allFR_Vals, groupTablePre.eiVals, 'Type','Spearman'));
    iesFR_Vals = cat(1, iesFR_Vals, corr(groupTablePre.iesFR_Vals, groupTablePre.eiVals, 'Type','Spearman'));
    isolFR_Vals = cat(1, isolFR_Vals, corr(groupTablePre.isolFR_Vals, groupTablePre.eiVals, 'Type','Spearman'));

                    
end

function plotDistro(x)
    cdfplot(x)
    hold on
    x_values = linspace(min(x),max(x));
    plot(x_values,normcdf(x_values,0,1),'r-')
    legend('Empirical CDF','Standard Normal CDF','Location','best')
    h = kstest(x)
end