function runTimeseries_compute(matrices, datafiles_origin, modelfiles_origin, surfaceROI, projectSettings, glmR2)
% <glmR2> is the full-surface GLMsingle full-model R^2 vector (modelOut{4}.R2,
% same vertex order as datafiles/modelfiles). Its mean is computed over the
% exact same vertex selection used for the ROI-average time series below,
% and passed through to the plot for display.
%
% runTimeseries_compute - ROI-average observed (%signal change) and
% predicted (GLMsingle) time series, per run, for the FULL run duration.
%
% This reuses the exact same PSC conversion and ROI-averaging logic as
% ttave_compute.m, but skips the trial-onset extraction/condition-grouping
% step entirely -- the point here is the continuous run time series, not a
% condition-triggered average.
%
% Computes and saves BOTH the all-polar-angles-combined trace (what's
% plotted) AND the 4 polar-angle-bin-pair traces (same grouping as
% ttave_compute.m: [0,180], [45,225], [90,270], [135,315]), so a later
% by-polar-angle analysis doesn't require rerunning the slow GLM
% reconstruction / raw data loading steps -- just this (fast) aggregation.

stimdur_s = projectSettings.stimdur_s;
tr_s = projectSettings.tr_s;

prfData = projectSettings.filteredPrfBins;   % vertex -> assigned polar angle bin, or NaN
polarAngleConditions = {[0,180], [45, 225], [90, 270], [135, 315]};

% which polar angle wedge(s) to include in the main (plotted) trace --
% defaults to all 8 canonical wedges (equivalent to the old ~isnan(prfData)
% behavior, since retriveRetData.m only ever assigns one of these 8 values
% or NaN) if projectSettings doesn't specify a restriction.
if isfield(projectSettings, 'polarAnglesToInclude') && ~isempty(projectSettings.polarAnglesToInclude)
    polarAnglesToInclude = projectSettings.polarAnglesToInclude;
else
    polarAnglesToInclude = [0, 45, 90, 135, 180, 225, 270, 315];
end

% NaN (not 0) outside the valid set, matching ttave_compute.m's convention --
% nanmedian below must be able to EXCLUDE non-qualifying vertices, not have
% them zeroed and drag the median toward 0.
surfacePA = nan(size(prfData))';
surfacePA(ismember(prfData', polarAnglesToInclude)) = 1;
surfaceSelection = surfaceROI .* surfacePA;

% per-polar-angle-bin selections (same NaN-outside convention)
surfaceSelection_byPA = cell(1, numel(polarAngleConditions));
for pa = 1:numel(polarAngleConditions)
    angleVals = polarAngleConditions{pa};
    surfacePA_bin = nan(size(prfData))';
    surfacePA_bin(ismember(prfData', angleVals)) = 1;
    surfaceSelection_byPA{pa} = surfaceROI .* surfacePA_bin;
end

% mean R2 (GLMsingle full model) over the exact vertex set used above --
% i.e. within the ROI, and meeting the eccentricity/vexpl/polar-angle
% criteria that produced surfacePA.
includedVerts = ~isnan(surfaceSelection);
meanR2 = nanmean(glmR2(includedVerts));

nRuns = numel(matrices);
dataRuns = cell(1, nRuns);
modelRuns = cell(1, nRuns);
dataRuns_byPA = cell(1, nRuns);   % each: 1 x 4 cell (one per polarAngleConditions bin)
modelRuns_byPA = cell(1, nRuns);
onsetsByType = cell(1, nRuns); % {motion_onsets, static_onsets, blank_onsets} per run, in TRs
staticOnsetsByOrientation = cell(1, nRuns); % struct .TR, .deg per run, static trials only

% static condition column order (9,10,11,12) -> orientation in degrees.
% 0=horizontal, 90=vertical, 45=right oblique, 135=left oblique (see
% AnalysisCode/01_process_singlesubjectGLM/format_desmats.m and
% createStimMap.m for the S_0/S_90/S_45/S_135 column order and the
% degree-to-orientation convention).
staticOrientationOrder = [0, 90, 45, 135];

for mi = 1:nRuns
    m = matrices{mi};

    % --- observed: raw BOLD -> %signal change ---
    df = datafiles_origin{mi};
    mean_overT = mean(abs(df), 2);
    timeseries_psc = ((df ./ mean_overT) - 1) * 100;
    dataRuns{mi} = nanmedian(timeseries_psc .* surfaceSelection);

    % --- predicted: GLMsingle fit, already in %signal change ---
    mf = modelfiles_origin{mi};
    modelRuns{mi} = nanmedian(mf .* surfaceSelection);

    % --- same, per polar-angle bin ---
    dataRuns_byPA{mi} = cell(1, numel(polarAngleConditions));
    modelRuns_byPA{mi} = cell(1, numel(polarAngleConditions));
    for pa = 1:numel(polarAngleConditions)
        dataRuns_byPA{mi}{pa} = nanmedian(timeseries_psc .* surfaceSelection_byPA{pa});
        modelRuns_byPA{mi}{pa} = nanmedian(mf .* surfaceSelection_byPA{pa});
    end

    % --- condition onsets, for background shading in the plot ---
    motionOnsets = []; staticOnsets = []; blankOnsets = [];
    [~, cond_n] = size(m);
    for ci = 1:cond_n
        onsets = find(m(:, ci) == 1);
        if isempty(onsets)
            continue
        end
        if ci <= 8
            motionOnsets = [motionOnsets; onsets]; %#ok<AGROW>
        elseif ci <= 12
            staticOnsets = [staticOnsets; onsets]; %#ok<AGROW>
        else
            blankOnsets = [blankOnsets; onsets]; %#ok<AGROW>
        end
    end
    onsetsByType{mi} = {sort(motionOnsets), sort(staticOnsets), sort(blankOnsets)};

    % --- per-orientation onsets for static (stationary) trials, for the
    % oriented tick marks drawn in the plot ---
    staticOnsetTRs = []; staticOnsetDegs = [];
    for ci = 9:12
        onsetsThisCond = find(m(:, ci) == 1);
        if isempty(onsetsThisCond)
            continue
        end
        staticOnsetTRs = [staticOnsetTRs; onsetsThisCond]; %#ok<AGROW>
        staticOnsetDegs = [staticOnsetDegs; repmat(staticOrientationOrder(ci - 8), numel(onsetsThisCond), 1)]; %#ok<AGROW>
    end
    [staticOnsetTRs, sortIdx] = sort(staticOnsetTRs);
    staticOnsetDegs = staticOnsetDegs(sortIdx);
    staticOnsetsByOrientation{mi} = struct('TR', staticOnsetTRs, 'deg', staticOnsetDegs);
end

runTimeseriesOutput.data = dataRuns;
runTimeseriesOutput.model = modelRuns;
runTimeseriesOutput.data_byPA = dataRuns_byPA;
runTimeseriesOutput.model_byPA = modelRuns_byPA;
runTimeseriesOutput.polarAngleConditions = polarAngleConditions;
runTimeseriesOutput.onsetsByType = onsetsByType;
runTimeseriesOutput.staticOnsetsByOrientation = staticOnsetsByOrientation;
runTimeseriesOutput.stimdur_s = stimdur_s;
runTimeseriesOutput.tr_s = tr_s;
runTimeseriesOutput.meanR2 = meanR2;
runTimeseriesOutput.nIncludedVerts = sum(includedVerts);

plot_runTimeseries(runTimeseriesOutput, projectSettings);

end