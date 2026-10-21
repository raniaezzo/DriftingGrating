function runWithinAcrossComparison(bidsDir, projectName, paFolderName, roiName, runWithin, runAcross, ylimOverride)
% runWithinAcrossComparison - "within vs. across" run comparison for the
% cross-subject group-average run time series.
%
% Loads the already-computed Run<runWithin>_ObservedandPredicted.mat and
% Run<runAcross>_ObservedandPredicted.mat files (produced by
% run_groupAverageRunTimeseries.m / plot_groupAverageRun.m) from
% <bidsDir>/derivatives/runtimeseries/<projectName>/<paFolderName>/, and
% plots run <runWithin>'s observed trace against:
%   (within) its own predicted trace (same as the regular group-average plot)
%   (across) run <runAcross>'s predicted trace, reused against run
%            <runWithin>'s observed data
%
% This tests how much of the model's fit is run-specific vs. driven by the
% block-design structure shared across every run. Output is saved into the
% same <paFolderName> folder it read from.
%
% <ylimOverride> (optional) forces the y-axis to a fixed [ymin ymax] range
% instead of auto-scaling to this call's own data -- pass the same range to
% both directions of a pair (e.g. via runWithinAcrossComparisonPaired.m) so
% "run A within/B across" and "run B within/A across" are directly visually
% comparable.

if nargin < 7
    ylimOverride = [];
end

savePath = fullfile(bidsDir, 'derivatives', 'runtimeseries', projectName, paFolderName);

withinFile = fullfile(savePath, sprintf('Run%d_ObservedandPredicted.mat', runWithin));
acrossFile = fullfile(savePath, sprintf('Run%d_ObservedandPredicted.mat', runAcross));

if ~exist(withinFile, 'file')
    error('runWithinAcrossComparison:missingFile', 'Run %d group-average output not found: %s (run run_groupAverageRunTimeseries.m first)', runWithin, withinFile);
end
if ~exist(acrossFile, 'file')
    error('runWithinAcrossComparison:missingFile', 'Run %d group-average output not found: %s (run run_groupAverageRunTimeseries.m first)', runAcross, acrossFile);
end

Sw = load(withinFile);
Sa = load(acrossFile);

if numel(Sw.meanObserved) ~= numel(Sa.predicted)
    error('runWithinAcrossComparison:lengthMismatch', ...
        'Run %d and run %d have different trace lengths (%d vs %d TRs) -- cannot compare directly.', ...
        runWithin, runAcross, numel(Sw.meanObserved), numel(Sa.predicted));
end

% across R^2: run <runWithin>'s observed data vs. run <runAcross>'s
% predicted trace, same 1-SSres/SStot formula used everywhere else
observed = Sw.meanObserved(:);
predictedAcross = Sa.predicted(:);
ssRes = sum((observed - predictedAcross).^2);
ssTot = sum((observed - mean(observed)).^2);
R2across = (1 - ssRes / ssTot) * 100;

compData.meanObserved = Sw.meanObserved;
compData.predictedWithin = Sw.predicted;
compData.predictedAcross = Sa.predicted;
compData.R2within = Sw.R2;
compData.R2across = R2across;
compData.onsets = Sw.onsets;
compData.staticOnsetsByOrientation = Sw.staticOnsets;
compData.stimdur_s = Sw.stimdur_s;
compData.tr_s = Sw.tr_s;
compData.runWithin = runWithin;
compData.runAcross = runAcross;
compData.ylimOverride = ylimOverride;

outputInfo.projectName = projectName;
outputInfo.roiName = roiName;
outputInfo.savePath = savePath;

plot_withinAcrossComparison(compData, outputInfo);

fprintf('  %s/%s run %d vs %d: within R^2 = %.1f%%, across R^2 = %.1f%%\n', ...
    projectName, paFolderName, runWithin, runAcross, Sw.R2, R2across);

end
