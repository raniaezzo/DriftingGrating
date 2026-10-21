function runWithinAcrossComparisonPaired(bidsDir, projectName, paFolderName, roiName, runA, runB)
% runWithinAcrossComparisonPaired - runs runWithinAcrossComparison.m in
% both directions for a pair of runs (A within/B across, then B within/A
% across), forcing both plots to share the same y-axis limits so they're
% directly visually comparable.
%
% The shared range is computed from the union of both runs' observed and
% predicted traces (all four traces that appear across the two plots).

savePath = fullfile(bidsDir, 'derivatives', 'runtimeseries', projectName, paFolderName);

Sa = load(fullfile(savePath, sprintf('Run%d_ObservedandPredicted.mat', runA)));
Sb = load(fullfile(savePath, sprintf('Run%d_ObservedandPredicted.mat', runB)));

allVals = [Sa.meanObserved(:); Sa.predicted(:); Sb.meanObserved(:); Sb.predicted(:)];
allVals = allVals(isfinite(allVals));
pad = 0.1 * (max(allVals) - min(allVals) + eps);
ylimOverride = [min(allVals) - pad, max(allVals) + pad];

runWithinAcrossComparison(bidsDir, projectName, paFolderName, roiName, runA, runB, ylimOverride);
runWithinAcrossComparison(bidsDir, projectName, paFolderName, roiName, runB, runA, ylimOverride);

end
