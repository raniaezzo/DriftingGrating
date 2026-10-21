function results = runRawTimeseriesPermutationTest(projectName, roiName, varargin)
% RUNRAWTIMESERIESPERMUTATIONTEST  Exact permutation test of GLMsingle's
% temporal specificity: does the group-average, GLMsingle-based predicted
% run time series match ITS OWN run's group-average observed data better
% than a MISMATCHED run's predicted trace would? Concatenates across ALL
% 8 canonical polar-angle locations AND all valid runs into one
% aggregate statistic (agreed 2026-09-01) -- a single location's own
% result is only used for the illustrative figure, not this stat.
%
%   results = runRawTimeseriesPermutationTest('dg', 'V1')
%   results = runRawTimeseriesPermutationTest('da', 'V1', 'subjectMode', 'matched')
%
% Name-value options:
%   'subjectMode'  - 'full' (default for dg: all 13 subjects) or 'matched'
%                    (the 7 subjects run in both dg and da, sub-0395
%                    excluded). da always uses the matched 7 regardless of
%                    this setting (it has no other cohort).
%   'nRunsCheck'   - cap on runs considered (default 8, matching
%                    run_groupAverageRunTimeseries.m's convention).
%   'bidsDir', 'githubDir' - as elsewhere in this pipeline.
%
% METHOD:
% For each subject, computeSubjectRunTraces_GLMsingle.m builds a
% continuous per-run ROI-median trace, separately for each of the 8
% canonical 45-deg wedges, for both the observed data (GLMsingle-meanvol
% %sc, THEN polynomial-denoised) and GLMsingle's own predicted
% reconstruction (modelmd betas x each vertex's own HRFindex, via
% GLMpredictresponses) -- no fresh/simplified refit anywhere, unlike
% refitGroupMeanTrace.m's canonical-HRF OLS approach. Per-(location,run)
% group averages (mean across subjects) use the SAME per-subject
% validity logic run_groupAverageRunTimeseries.m already uses (run
% exists, not a duplicate design per findDuplicateDesignRuns.m, not
% all-NaN, matching trace length) -- a run is kept in the aggregate ONLY
% if it is valid at EVERY one of the 8 locations, so there is one shared
% set of valid runs the derangements below are built over.
%
% TRUE statistic: for each of the 8 locations, concatenate every valid
% run i's group-average observed trace (obs_i) end to end across runs;
% separately concatenate that location's own group-average predicted
% trace (pred_i, SAME run i -- the correct match) in the same run order.
% Then concatenate all 8 locations' (observed) vectors together, and all
% 8 locations' (predicted) vectors together in the same order, and
% compute ONE R^2 = 100*(1 - sum((obsCat-predCat).^2)/sum((obsCat-mean(obsCat)).^2))
% on the whole thing (centered R^2, matching runWithinAcrossComparison.m's
% existing convention).
%
% NULL distribution: every DERANGEMENT of the valid run indices (a
% permutation with no fixed point). The SAME derangement is applied to
% EVERY location within one null draw (a run's trial-timing mismatch is
% a property of the run itself, shared across all locations measured in
% that run -- using a different random mismatch per location within one
% draw would conflate two different questions). With <=8 runs this is
% EXHAUSTIVE (all derangements enumerated, not Monte Carlo sampled -- at
% n=8 there are only 14,833, trivially fast to enumerate exactly).
%
% This is a TEMPORAL SPECIFICITY test, not a held-out generalization
% test: GLMsingle fits ONE set of betas per subject pooling across all
% runs (not separately per run), so pred_i and pred_j are the SAME
% underlying betas reconstructed through DIFFERENT runs' own trial-onset
% designs. The question this answers is whether the model's response
% genuinely tracks which run's specific trial sequence occurred, not
% whether betas trained on one run generalize to unseen data.
%
% Returns a struct: trueR2, nullR2 (1 x nDerangements), pValue (fraction
% of null >= true), validRuns, nSubjectsPerRun (nLocations x nValidRuns),
% wedgeAngles, roiName, projectName, subjectMode, nDerangements.

p = inputParser;
p.addParameter('subjectMode', 'full', @(x) ismember(x, {'full','matched'}));
p.addParameter('nRunsCheck', 8, @isnumeric);
p.addParameter('bidsDir', '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/', @ischar);
p.addParameter('githubDir', '~/Documents/GitHub', @ischar);
p.parse(varargin{:});
opt = p.Results;

bidsDir = opt.bidsDir;
githubDir = opt.githubDir;
nRunsCheck = opt.nRunsCheck;

addpath(genpath(fullfile(githubDir, 'DriftingGrating', 'AnalysisCode')));
addpath(genpath(fullfile(githubDir, 'atlasmgz')));
cd(fullfile(githubDir, 'DriftingGrating', 'AnalysisCode'));
setup_user('rania', bidsDir);

dg_subjects_13 = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
    'sub-0395', 'sub-0426', 'sub-0250', 'sub-0442', 'sub-wlsubj121', 'sub-wlsubj127', ...
    'sub-0397', 'sub-0427'};
matched_7 = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
    'sub-0426', 'sub-0250'}; % sub-0395 excluded (mismatched da pilot stimulus)

if strcmp(projectName, 'dg')
    if strcmp(opt.subjectMode, 'full')
        subjects = dg_subjects_13;
    else
        subjects = matched_7;
    end
elseif strcmp(projectName, 'da')
    subjects = matched_7; % da has no other cohort
else
    error('runRawTimeseriesPermutationTest:project', 'projectName must be ''dg'' or ''da''.');
end
nSubj = numel(subjects);

projectSettingsBase = loadConfig(githubDir);
jsonParams = jsondecode(fileread('setup.json'));
stimdur_s = jsonParams.stimdur_s.Val;
tr_s = jsonParams.tr_s.Val;

wedgeAngles = 0:45:315;
nWedges = numel(wedgeAngles);

%% per-subject pass (one call per subject computes all 8 locations at once)

subjData = struct('subj', {}, 'matrices_onset', {}, 'observed', {}, 'predicted', {}, 'nRuns', {}, 'excludedRunIdx', {});
for ss = 1:nSubj
    subj = subjects{ss};

    subjectDir = fullfile(bidsDir, 'derivatives', sprintf('%sGLM', projectName), 'hRF_glmsingle', subj);
    contents = dir(subjectDir);
    sesNames = {contents([contents.isdir] & startsWith({contents.name}, 'ses-')).name};
    % Keep only session folders that actually contain GLM output
    % (rawInfo.mat) -- some subjects (found 2026-09-01: dg/sub-wlsubj121)
    % have a stray EMPTY ses-* folder alongside the real one (here,
    % ses-yu3t02 vs. the genuine ses-nyu3t02, which also matches the
    % session name in that subject's actual fmriprep BOLD filenames) --
    % this is a leftover/incomplete folder, not a genuine ambiguity, so it
    % should be filtered out rather than treated as an error.
    sesNames = sesNames(cellfun(@(s) isfile(fullfile(subjectDir, s, 'rawInfo.mat')), sesNames));
    if numel(sesNames) > 1
        error('Multiple session folders WITH GLM output for %s/%s: %s.', projectName, subj, strjoin(sesNames, ', '))
    elseif isempty(sesNames)
        warning('No session folder with GLM output for %s/%s, skipping.', projectName, subj)
        continue
    end
    ses = sesNames{1};
    derivativesFolder = fullfile(bidsDir, 'derivatives', sprintf('%sGLM',projectName), 'hRF_glmsingle', subj, ses);

    rawS = load(fullfile(derivativesFolder, 'rawInfo.mat'), 'matrices_onset');

    fprintf('  computing GLMsingle-based run traces (all 8 locations): %s / %s ...\n', projectName, subj);
    st = computeSubjectRunTraces_GLMsingle(bidsDir, projectName, subj, ses, projectSettingsBase, ...
        stimdur_s, tr_s, rawS.matrices_onset, nRunsCheck, roiName);

    excludedRunIdx = findDuplicateDesignRuns(rawS.matrices_onset, nRunsCheck);
    if ~isempty(excludedRunIdx)
        fprintf('  %s/%s: run(s) %s have a design identical to an earlier run -- excluding from group average\n', ...
            projectName, subj, mat2str(excludedRunIdx));
    end

    idx = numel(subjData) + 1;
    subjData(idx).subj = subj;
    subjData(idx).matrices_onset = rawS.matrices_onset;
    subjData(idx).observed = st.observed;   % 1 x nWedges cell, each 1 x nRuns cell
    subjData(idx).predicted = st.predicted;
    subjData(idx).nRuns = numel(st.observed{1});
    subjData(idx).excludedRunIdx = excludedRunIdx;
end

%% per-(location, run) group average, same validity logic as run_groupAverageRunTimeseries.m

groupObserved = cell(nWedges, nRunsCheck);
groupPredicted = cell(nWedges, nRunsCheck);
nSubjectsPerRunLoc = zeros(nWedges, nRunsCheck);

for w = 1:nWedges
    for r = 1:nRunsCheck
        obsMat = []; predMat = [];
        for ss = 1:numel(subjData)
            hasRun = subjData(ss).nRuns >= r;
            isDuplicateDesign = ismember(r, subjData(ss).excludedRunIdx);
            obs_r = []; pred_r = [];
            if hasRun
                obs_r = subjData(ss).observed{w}{r};
                pred_r = subjData(ss).predicted{w}{r};
            end
            validRun = hasRun && ~isDuplicateDesign && ~all(isnan(obs_r)) && ~all(isnan(pred_r)) && ...
                (isempty(obsMat) || numel(obs_r) == size(obsMat, 2));
            if validRun
                obsMat = [obsMat; obs_r]; %#ok<AGROW>
                predMat = [predMat; pred_r]; %#ok<AGROW>
            end
        end
        if isempty(obsMat)
            continue
        end
        groupObserved{w,r} = nanmean(obsMat, 1);
        groupPredicted{w,r} = nanmean(predMat, 1);
        nSubjectsPerRunLoc(w,r) = size(obsMat, 1);
    end
end

% a run is usable for the aggregate stat only if it is valid at EVERY
% location, so all 8 locations share one common set of valid runs (and
% therefore one shared derangement, per the header note)
validPerRun = all(~cellfun(@isempty, groupObserved), 1);
validRuns = find(validPerRun);
nRuns = numel(validRuns);
if nRuns < 3
    error('runRawTimeseriesPermutationTest:tooFewRuns', ...
        'Only %d run(s) valid across all 8 locations for %s/%s -- need at least 3 for a meaningful derangement set.', ...
        nRuns, projectName, roiName);
end

%% true statistic: concatenate obs_i / pred_i (i=i) across all valid runs, THEN across all 8 locations

obsConcat = [];
predConcatTrue = [];
for w = 1:nWedges
    obsConcat = [obsConcat, cell2mat(groupObserved(w, validRuns))]; %#ok<AGROW>
    predConcatTrue = [predConcatTrue, cell2mat(groupPredicted(w, validRuns))]; %#ok<AGROW>
end
trueR2 = centeredR2(obsConcat, predConcatTrue);

%% null distribution: every derangement of the valid runs, SAME derangement applied to every location

allPerms = perms(1:nRuns); % nRuns! x nRuns
isDerangement = all(allPerms ~= (1:nRuns), 2);
derangements = allPerms(isDerangement, :);
nDerangements = size(derangements, 1);

nullR2 = nan(1, nDerangements);
for d = 1:nDerangements
    permutedRuns = validRuns(derangements(d,:));
    predConcatNull = [];
    for w = 1:nWedges
        predConcatNull = [predConcatNull, cell2mat(groupPredicted(w, permutedRuns))]; %#ok<AGROW>
    end
    nullR2(d) = centeredR2(obsConcat, predConcatNull);
end

pValue = mean(nullR2 >= trueR2);

results.trueR2 = trueR2;
results.nullR2 = nullR2;
results.pValue = pValue;
results.validRuns = validRuns;
results.nSubjectsPerRun = nSubjectsPerRunLoc(:, validRuns);
results.wedgeAngles = wedgeAngles;
results.roiName = roiName;
results.projectName = projectName;
results.subjectMode = opt.subjectMode;
results.nDerangements = nDerangements;

if pValue == 0
    pStr = sprintf('< %.5g (1/%d)', 1/(nDerangements+1), nDerangements+1);
else
    pStr = sprintf('= %.5f', pValue);
end
fprintf('\n%s (%s, n=%d) / %s / all 8 locations x %d runs: true R^2 = %.2f%%, null R^2 = %.2f%% [%.2f, %.2f] (median [2.5,97.5], %d derangements), p %s\n', ...
    projectName, opt.subjectMode, nSubj, roiName, nRuns, ...
    trueR2, median(nullR2), prctile(nullR2,2.5), prctile(nullR2,97.5), nDerangements, pStr);

end

function r2 = centeredR2(obs, pred)
    ssRes = sum((obs - pred).^2);
    ssTot = sum((obs - mean(obs)).^2);
    r2 = (1 - ssRes/ssTot) * 100;
end
