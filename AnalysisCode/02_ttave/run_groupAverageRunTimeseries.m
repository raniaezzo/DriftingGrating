% run_groupAverageRunTimeseries - cross-subject group-average observed vs.
% predicted run time series (group-average counterpart to
% run_runTimeseries.m). Individual-subject GLMsingle per-vertex HRF fits
% can't be averaged directly, so for each run, a fresh canonical-HRF GLM
% is refit to the cross-subject group-mean observed trace itself
% (refitGroupMeanTrace.m).
%
% Also includes a "within vs. across" run-comparison section: for two
% chosen runs (runWithin, runAcross) of a given project/ROI/polar-angle
% selection, overlay run <runWithin>'s observed trace against its own
% ("within") predicted trace and run <runAcross>'s ("across") predicted
% trace, and report both R^2 values. Uses the already-computed per-run
% group-average outputs, so it's fast to rerun with different run pairs.
%
% Output: <bidsDir>/derivatives/runtimeseries/<project>/<roiName>_pa_<sel>/
%         Run<N>_ObservedandPredicted.{pdf,fig,mat}
%         Run<W>within_Run<A>across_Comparison.{pdf,fig,mat}

clear all; close all; clc;

% must be in DriftingGrating directory to run

addpath(genpath(pwd));
bidsDir =  '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
githubDir = '~/Documents/GitHub';
addpath(genpath(fullfile(githubDir, 'atlasmgz')));
setup_user('rania', bidsDir)

projectSettingsBase = loadConfig(githubDir);

hRF_setting = 'glmsingle';
projects = {'dg', 'da'};

% ROI to analyze (matches a filename in general/jsons/ROIS_ALL.json).
roiName = 'V1';

% Polar angle wedges (deg) to include in the vertex selection, on top of
% the <roiName> + eccentricity + pRF R^2 criteria. Defaults to all 8
% canonical 45-deg-wide wedges (i.e. the whole visual field, no polar
% angle restriction) -- change this to a subset (e.g. [0] for the
% 0+/-22.5deg horizontal-right wedge) to restrict to specific polar angle
% location(s). Output is saved into a subfolder named after the ROI and
% selection: '<roiName>_pa_all' for the full canonical set, otherwise
% '<roiName>_pa_<angle>_<angle>...'.
polarAnglesToInclude = [0, 45, 90, 135, 180, 225, 270, 315];
polarAngleBinWidth = 45;

% cap on how many runs to include in the group average -- subjects with
% fewer runs are excluded (per-run) from the average once their run count
% is exceeded; subjects with more runs (e.g. a 9th run) are trimmed.
nRunsCheck = 8;

% subjects shared between the two projects
dgSubs = dir(fullfile(bidsDir, 'derivatives', 'dgGLM', 'hRF_glmsingle', 'sub-*'));
dgSubs = {dgSubs([dgSubs.isdir]).name};
daSubs = dir(fullfile(bidsDir, 'derivatives', 'daGLM', 'hRF_glmsingle', 'sub-*'));
daSubs = {daSubs([daSubs.isdir]).name};
subjects = intersect(dgSubs, daSubs);

fprintf('%d subject(s) shared across dg/da\n', numel(subjects));

jsonParams = jsondecode(fileread('setup.json'));
stimdur_s = jsonParams.stimdur_s.Val;
tr_s = jsonParams.tr_s.Val;

allCanonicalAngles = [0, 45, 90, 135, 180, 225, 270, 315];
if isequal(sort(polarAnglesToInclude), allCanonicalAngles)
    paFolderName = sprintf('%s_pa_all', roiName);
else
    angleStrs = arrayfun(@(a) num2str(a), sort(polarAnglesToInclude), 'UniformOutput', false);
    paFolderName = sprintf('%s_pa_%s', roiName, strjoin(angleStrs, '_'));
end

% static condition column order (9,10,11,12) -> orientation in degrees.
% 0=horizontal, 90=vertical, 45=right oblique, 135=left oblique (see
% AnalysisCode/01_process_singlesubjectGLM/format_desmats.m and
% createStimMap.m for the S_0/S_90/S_45/S_135 column order and the
% degree-to-orientation convention).
staticOrientationOrder = [0, 90, 45, 135];

%% ===== main pass: per-run cross-subject group average =====
for pp = 1:numel(projects)
    projectName = projects{pp};
    fprintf('=== %s ===\n', projectName);

    % --- per-subject pass: load design + compute traces once per subject
    % (one data load per subject, not per run) ---
    subjData = struct('subj', {}, 'matrices_onset', {}, 'traces', {}, 'nRuns', {}, 'excludedRunIdx', {});
    for ss = 1:numel(subjects)
        subj = subjects{ss};

        subjectDir = fullfile(bidsDir,'derivatives',strcat(projectName, 'GLM'), strcat('hRF_',hRF_setting), subj);
        contents = dir(subjectDir);
        sesNames = {};
        for i = 1:length(contents)
            if contents(i).isdir && startsWith(contents(i).name, 'ses-')
                sesNames{end+1} = contents(i).name;
            end
        end
        if length(sesNames) > 1
            error('Multiple session folders for %s/%s.', projectName, subj)
        elseif isempty(sesNames)
            warning('No session folder for %s/%s, skipping.', projectName, subj)
            continue
        end
        ses = sesNames{1};
        derivativesFolder = fullfile(bidsDir, 'derivatives', sprintf('%sGLM',projectName), ...
            sprintf('hRF_%s', hRF_setting), subj, ses);

        rawS = load(fullfile(derivativesFolder, 'rawInfo.mat'), 'matrices_onset');

        fprintf('  computing traces: %s / %s ...\n', projectName, subj);
        st = computeSubjectTraces(bidsDir, projectName, subj, ses, projectSettingsBase, ...
            stimdur_s, tr_s, rawS.matrices_onset, nRunsCheck, roiName, polarAnglesToInclude);

        % detect runs whose design is identical to an earlier run for this
        % subject (same trial sequence shown/used more than once) --
        % these aren't independent repeats and would double-count that
        % trial sequence in the group average. See findDuplicateDesignRuns.m
        excludedRunIdx = findDuplicateDesignRuns(rawS.matrices_onset, nRunsCheck);
        if ~isempty(excludedRunIdx)
            fprintf('  %s/%s: run(s) %s have a design identical to an earlier run -- excluding from group average\n', ...
                projectName, subj, mat2str(excludedRunIdx));
        end

        idx = numel(subjData) + 1;
        subjData(idx).subj = subj;
        subjData(idx).matrices_onset = rawS.matrices_onset;
        subjData(idx).traces = st.traces;
        subjData(idx).nRuns = numel(st.traces);
        subjData(idx).excludedRunIdx = excludedRunIdx;
    end

    savePath = fullfile(bidsDir, 'derivatives', 'runtimeseries', projectName, paFolderName);

    % --- per-run pass: cross-subject group average + refit ---
    for r = 1:nRunsCheck
        includedSubjects = {};
        excludedSubjects = {};
        traceMat = [];
        refSubjIdx = [];
        for ss = 1:numel(subjData)
            hasRun = subjData(ss).nRuns >= r;
            isDuplicateDesign = ismember(r, subjData(ss).excludedRunIdx);
            trace_r = [];
            if hasRun
                trace_r = subjData(ss).traces{r};
            end
            % exclude if the run doesn't exist, is all-NaN, has a design
            % identical to an earlier run for this subject (see
            % findDuplicateDesignRuns.m), or doesn't match the reference
            % trace length (defensive -- shouldn't normally happen, but a
            % partial/aborted run could differ)
            validRun = hasRun && ~isDuplicateDesign && ~all(isnan(trace_r)) && ...
                (isempty(traceMat) || numel(trace_r) == size(traceMat, 2));
            if validRun
                traceMat = [traceMat; trace_r]; %#ok<AGROW>
                includedSubjects{end+1} = subjData(ss).subj; %#ok<AGROW>
                if isempty(refSubjIdx), refSubjIdx = ss; end
            else
                excludedSubjects{end+1} = subjData(ss).subj; %#ok<AGROW>
            end
        end

        if isempty(traceMat)
            warning('No subjects have run %d for %s -- skipping.', r, projectName);
            continue
        end

        meanObserved = nanmean(traceMat, 1);

        designMatrix = subjData(refSubjIdx).matrices_onset{r};
        [predicted, R2, ~] = refitGroupMeanTrace(meanObserved, designMatrix, stimdur_s, tr_s);

        % condition onsets, for background shading + orientation ticks --
        % same extraction as runTimeseries_compute.m. The design matrix is
        % identical across subjects for a given run (shared block design),
        % so any included subject's matrices_onset{r} is a valid reference.
        motionOnsets = []; staticOnsetsTR = []; blankOnsets = [];
        [~, cond_n] = size(designMatrix);
        for ci = 1:cond_n
            onsetsHere = find(designMatrix(:, ci) == 1);
            if isempty(onsetsHere)
                continue
            end
            if ci <= 8
                motionOnsets = [motionOnsets; onsetsHere]; %#ok<AGROW>
            elseif ci <= 12
                staticOnsetsTR = [staticOnsetsTR; onsetsHere]; %#ok<AGROW>
            else
                blankOnsets = [blankOnsets; onsetsHere]; %#ok<AGROW>
            end
        end
        onsets = {sort(motionOnsets), sort(staticOnsetsTR), sort(blankOnsets)};

        staticOnsetTRs = []; staticOnsetDegs = [];
        for ci = 9:12
            onsetsHere = find(designMatrix(:, ci) == 1);
            if isempty(onsetsHere)
                continue
            end
            staticOnsetTRs = [staticOnsetTRs; onsetsHere]; %#ok<AGROW>
            staticOnsetDegs = [staticOnsetDegs; repmat(staticOrientationOrder(ci - 8), numel(onsetsHere), 1)]; %#ok<AGROW>
        end
        [staticOnsetTRs, sortIdx] = sort(staticOnsetTRs);
        staticOnsetDegs = staticOnsetDegs(sortIdx);
        staticOnsetsByOrientation = struct('TR', staticOnsetTRs, 'deg', staticOnsetDegs);

        groupData = struct();
        groupData.meanObserved = meanObserved;
        groupData.predicted = predicted;
        groupData.R2 = R2;
        groupData.onsets = onsets;
        groupData.staticOnsetsByOrientation = staticOnsetsByOrientation;
        groupData.stimdur_s = stimdur_s;
        groupData.tr_s = tr_s;
        groupData.includedSubjects = includedSubjects;
        groupData.excludedSubjects = excludedSubjects;

        outputInfo = struct();
        outputInfo.projectName = projectName;
        outputInfo.roiName = roiName;
        outputInfo.runIdx = r;
        outputInfo.savePath = savePath;
        outputInfo.polarAnglesToInclude = polarAnglesToInclude;
        outputInfo.polarAngleBinWidth = polarAngleBinWidth;

        plot_groupAverageRun(groupData, outputInfo);
        close all;

        fprintf('  run %d: n=%d included (%s), R^2 = %.1f%%\n', r, numel(includedSubjects), ...
            strjoin(includedSubjects, ','), R2);
    end
end

%% ===== within vs. across run comparison =====
% Uses the already-computed Run<N>_ObservedandPredicted.mat files from the
% pass above (fast -- no data reload). To compare a specific pair of runs
% for one project, call directly, e.g.:
%   runWithinAcrossComparison(bidsDir, 'dg', paFolderName, roiName, 1, 2)
% The loop below runs a default demo pair (run 1 within, run 2 across) for
% both projects.

runWithin = 1;
runAcross = 2;

for pp = 1:numel(projects)
    runWithinAcrossComparison(bidsDir, projects{pp}, paFolderName, roiName, runWithin, runAcross);
end
