function subjTraces = computeSubjectTraces(bidsDir, projectName, subj, ses, projectSettingsBase, stimdur_s, tr_s, matrices_onset, nRunsCheck, roiName, polarAnglesToInclude)
% computeSubjectTraces - per-run ROI-average %signal-change trace for one
% subject, vertex-selected to <roiName> + eccentricity 4-8 deg + pRF R^2 >=
% 0.1, optionally restricted to a subset of polar angle wedges. Used as
% input to the cross-subject group-average run-timeseries refit.
%
% <roiName> defaults to 'V1' if omitted/empty.
% <polarAnglesToInclude> defaults to all 8 canonical 45-deg wedges (i.e. no
% polar angle restriction) if omitted/empty -- see runTimeseries_compute.m
% for the same convention.

if nargin < 10 || isempty(roiName)
    roiName = 'V1';
end
if nargin < 11 || isempty(polarAnglesToInclude)
    polarAnglesToInclude = [0, 45, 90, 135, 180, 225, 270, 315];
end

hSize = get_surfsize(subj);
nRuns_subj = numel(matrices_onset);
nRunsUse = min(nRuns_subj, nRunsCheck);
run_ = 1:nRunsUse;
datafiles = load_data(bidsDir, projectName, 'fsnative', '.mgh', subj, ses, run_);

projectSettings = projectSettingsBase;
projectSettings.projectName = projectName;
projectSettings.subject = subj;
projectSettings.bidsDir = bidsDir;
projectSettings.retFolder = 'prfvista_mov';
projectSettings.ses = ses;
projectSettings.minECC = 4;
projectSettings.maxECC = 8;
projectSettings.minVAREXP = .1;
projectSettings.stimdur_s = stimdur_s;
projectSettings.tr_s = tr_s;
projectSettings.roiName = roiName;
projectSettings.polarAngleBinWidth = 45;

primaryROIvertices = getROIidxs(subj, projectSettings.roiName, hSize);
surfaceROI = nan(sum(hSize), 1);
surfaceROI(primaryROIvertices) = 1;

filteredPrfBins = retriveRetData(projectSettings);
surfacePA = nan(size(filteredPrfBins))';
surfacePA(ismember(filteredPrfBins', polarAnglesToInclude)) = 1;
surfaceSelection = surfaceROI .* surfacePA;

subjTraces.nIncludedVerts = sum(~isnan(surfaceSelection));
subjTraces.traces = cell(1, nRunsUse);
for r = 1:nRunsUse
    df = datafiles{r};
    mean_overT = mean(abs(df), 2);
    psc = ((df ./ mean_overT) - 1) * 100;
    subjTraces.traces{r} = nanmedian(psc .* surfaceSelection);
end

end
