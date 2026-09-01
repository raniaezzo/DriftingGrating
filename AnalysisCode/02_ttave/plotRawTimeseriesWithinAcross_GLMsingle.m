function plotRawTimeseriesWithinAcross_GLMsingle(projectName, roiName, targetWedgeAngles, runA, runB, varargin)
% PLOTRAWTIMESERIESWITHINACROSS_GLMSINGLE  For each requested polar-angle
% location, produces two figures: (1) run A's group-average observed data
% with run A's OWN ("within") GLMsingle-predicted trace and run B's
% ("across", mismatched) predicted trace overlaid, and (2) the same with
% A/B swapped. Each plot reports both R^2 values in the title.
%
%   plotRawTimeseriesWithinAcross_GLMsingle('dg', 'V1', [0, 90], 2, 5)
%
% Illustrative companion to runRawTimeseriesPermutationTest.m's aggregate
% stat -- this shows individual example run-pairs; the aggregate stat is
% the actual quantitative claim (all 8 locations x all valid runs,
% derangement-based null). Uses the SAME GLMsingle-based reconstruction
% (computeSubjectRunTraces_GLMsingle.m: GLMsingle's own modelmd betas x
% each vertex's own HRFindex, via GLMpredictresponses -- not
% refitGroupMeanTrace.m's fresh canonical-HRF OLS refit) and the SAME
% per-subject validity logic (run exists, not a duplicate design per
% findDuplicateDesignRuns.m, not all-NaN) as that script, so the R^2
% values reported here are directly comparable to (individual draws from)
% that aggregate test's null/true distinction -- run A's own fit vs run
% A's data is one "true" (i=i) comparison; run B's fit vs run A's data is
% one specific draw from what the null distribution is built from.
%
% Name-value options:
%   'subjectMode', 'nRunsCheck', 'bidsDir', 'githubDir' - as in
%   runRawTimeseriesPermutationTest.m.
%   'figureDir' - where to save PDFs (default:
%                 <bidsDir>/derivatives/runtimeseries_GLMsingle/<project>/<roi>/)

p = inputParser;
p.addParameter('subjectMode', 'full', @(x) ismember(x, {'full','matched'}));
p.addParameter('nRunsCheck', 8, @isnumeric);
p.addParameter('bidsDir', '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/', @ischar);
p.addParameter('githubDir', '~/Documents/GitHub', @ischar);
p.addParameter('figureDir', '', @ischar);
p.parse(varargin{:});
opt = p.Results;

bidsDir = opt.bidsDir;
githubDir = opt.githubDir;
nRunsCheck = opt.nRunsCheck;

addpath(genpath(fullfile(githubDir, 'DriftingGrating', 'AnalysisCode')));
addpath(genpath(fullfile(githubDir, 'atlasmgz')));
cd(fullfile(githubDir, 'DriftingGrating', 'AnalysisCode'));
setup_user('rania', bidsDir);

if isempty(opt.figureDir)
    figureDir = fullfile(bidsDir, 'derivatives', 'runtimeseries_GLMsingle', projectName, roiName);
else
    figureDir = opt.figureDir;
end
if ~isfolder(figureDir), mkdir(figureDir); end

dg_subjects_13 = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
    'sub-0395', 'sub-0426', 'sub-0250', 'sub-0442', 'sub-wlsubj121', 'sub-wlsubj127', ...
    'sub-0397', 'sub-0427'};
matched_7 = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
    'sub-0426', 'sub-0250'};

if strcmp(projectName, 'dg')
    if strcmp(opt.subjectMode, 'full')
        subjects = dg_subjects_13;
    else
        subjects = matched_7;
    end
elseif strcmp(projectName, 'da')
    subjects = matched_7;
else
    error('plotRawTimeseriesWithinAcross_GLMsingle:project', 'projectName must be ''dg'' or ''da''.');
end
nSubj = numel(subjects);

projectSettingsBase = loadConfig(githubDir);
jsonParams = jsondecode(fileread('setup.json'));
stimdur_s = jsonParams.stimdur_s.Val;
tr_s = jsonParams.tr_s.Val;

wedgeAngles = 0:45:315;

%% per-subject pass (all 8 locations computed in one pass; cheap to keep all, select later)

subjData = struct('subj', {}, 'observed', {}, 'predicted', {}, 'nRuns', {}, 'excludedRunIdx', {}, 'matrices_onset', {});
for ss = 1:nSubj
    subj = subjects{ss};

    subjectDir = fullfile(bidsDir, 'derivatives', sprintf('%sGLM', projectName), 'hRF_glmsingle', subj);
    contents = dir(subjectDir);
    sesNames = {contents([contents.isdir] & startsWith({contents.name}, 'ses-')).name};
    sesNames = sesNames(cellfun(@(s) isfile(fullfile(subjectDir, s, 'rawInfo.mat')), sesNames));
    if numel(sesNames) > 1
        error('Multiple session folders with GLM output for %s/%s: %s.', projectName, subj, strjoin(sesNames, ', '))
    elseif isempty(sesNames)
        warning('No session folder with GLM output for %s/%s, skipping.', projectName, subj)
        continue
    end
    ses = sesNames{1};
    derivativesFolder = fullfile(bidsDir, 'derivatives', sprintf('%sGLM',projectName), 'hRF_glmsingle', subj, ses);

    rawS = load(fullfile(derivativesFolder, 'rawInfo.mat'), 'matrices_onset');

    fprintf('  computing GLMsingle-based run traces: %s / %s ...\n', projectName, subj);
    st = computeSubjectRunTraces_GLMsingle(bidsDir, projectName, subj, ses, projectSettingsBase, ...
        stimdur_s, tr_s, rawS.matrices_onset, nRunsCheck, roiName);

    excludedRunIdx = findDuplicateDesignRuns(rawS.matrices_onset, nRunsCheck);

    idx = numel(subjData) + 1;
    subjData(idx).subj = subj;
    subjData(idx).observed = st.observed;
    subjData(idx).predicted = st.predicted;
    subjData(idx).nRuns = numel(st.observed{1});
    subjData(idx).excludedRunIdx = excludedRunIdx;
    subjData(idx).matrices_onset = rawS.matrices_onset;
end

%% group-average per (location, run), same validity logic as runRawTimeseriesPermutationTest.m

runsNeeded = unique([runA, runB]);
groupObserved = containers.Map('KeyType','double','ValueType','any'); % key: wedgeAngle*100+run
groupPredicted = containers.Map('KeyType','double','ValueType','any');

for wi = 1:numel(targetWedgeAngles)
    wedgeAngle = targetWedgeAngles(wi);
    w = find(wedgeAngles == wedgeAngle);
    if isempty(w)
        error('plotRawTimeseriesWithinAcross_GLMsingle:wedge', '%d is not one of the 8 canonical wedge angles.', wedgeAngle);
    end
    for r = runsNeeded
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
            error('plotRawTimeseriesWithinAcross_GLMsingle:noData', 'No valid subjects for %s/%s wedge=%d run=%d.', projectName, roiName, wedgeAngle, r);
        end
        key = wedgeAngle*100 + r;
        groupObserved(key) = nanmean(obsMat, 1);
        groupPredicted(key) = nanmean(predMat, 1);
    end
end

%% stimulus-onset categories per run, for the background shading, from the
% first processed subject (the block design is identical across subjects
% for a given run). Moving conditions (columns 1-8) are pooled into one
% grey-shaded group; each of the 4 stationary conditions (columns 9-12)
% gets its OWN color, by raw stimulus code -- code 0->cyan, 90->yellow,
% 45->dark green, 135->magenta. This mapping is project-invariant: for dg
% those codes are horizontal/vertical/right-leaning/left-leaning, for da
% they are annulus(tangential)/pinwheel(radial)/CW spiral/CCW spiral, but
% the color-to-raw-code assignment agreed 2026-09-01 is the same either
% way. Columns 9,10,11,12 are in raw-code order [0,90,45,135] (same
% convention as run_groupAverageRunTimeseries.m's staticOrientationOrder),
% so STATIC_COLORS below is indexed positionally in that same order.
% Column 13 (blank) is not shaded -- blends with the white background.

STATIC_COLORS = {[0 1 1], [1 1 0], [0 0.5 0], [1 0 1]}; % cyan, yellow, dark green, magenta, for raw codes [0,90,45,135] respectively
MOTION_COLOR = [0.6 0.6 0.6]; % grey, all 8 moving conditions

refSubjMatricesOnset = subjData(1).matrices_onset;
onsetsByRun = containers.Map('KeyType','double','ValueType','any');
for r = runsNeeded
    designMatrix = refSubjMatricesOnset{r};
    motionOnsets = [];
    staticOnsetsByCode = cell(1,4); % same order as STATIC_CODE_ORDER
    [~, cond_n] = size(designMatrix);
    for ci = 1:cond_n
        onsetsHere = find(designMatrix(:, ci) == 1);
        if isempty(onsetsHere)
            continue
        end
        if ci <= 8
            motionOnsets = [motionOnsets; onsetsHere]; %#ok<AGROW>
        elseif ci <= 12
            staticOnsetsByCode{ci-8} = onsetsHere;
        end
    end
    onsetsByRun(r) = [{sort(motionOnsets)}, staticOnsetsByCode];
end
shadeColorsAll = [{MOTION_COLOR}, STATIC_COLORS]; % same 5-group order as onsetsByRun's cells

%% shared y-axis limits across every plot this call produces (agreed 2026-09-01:
% one shared scale per project/call, not per-plot), computed from all 4
% (2 locations x 2 directions) traces up front.

allVals = [];
for wi = 1:numel(targetWedgeAngles)
    wedgeAngle = targetWedgeAngles(wi);
    vObsA = groupObserved(wedgeAngle*100+runA); vPredA = groupPredicted(wedgeAngle*100+runA);
    vObsB = groupObserved(wedgeAngle*100+runB); vPredB = groupPredicted(wedgeAngle*100+runB);
    allVals = [allVals; vObsA(:); vPredA(:); vObsB(:); vPredB(:)]; %#ok<AGROW>
end
allVals = allVals(isfinite(allVals));
padShared = 0.1 * (max(allVals) - min(allVals) + eps);
sharedYlim = [min(allVals)-padShared, max(allVals)+padShared];

%% plot: for each location, (A within/B across) and (B within/A across)

for wi = 1:numel(targetWedgeAngles)
    wedgeAngle = targetWedgeAngles(wi);
    obsA = groupObserved(wedgeAngle*100 + runA);
    predA = groupPredicted(wedgeAngle*100 + runA);
    obsB = groupObserved(wedgeAngle*100 + runB);
    predB = groupPredicted(wedgeAngle*100 + runB);

    plotOnePair(obsA, predA, predB, runA, runB, wedgeAngle, projectName, roiName, nSubj, opt.subjectMode, figureDir, onsetsByRun(runA), shadeColorsAll, stimdur_s, tr_s, sharedYlim);
    plotOnePair(obsB, predB, predA, runB, runA, wedgeAngle, projectName, roiName, nSubj, opt.subjectMode, figureDir, onsetsByRun(runB), shadeColorsAll, stimdur_s, tr_s, sharedYlim);
end

end

function plotOnePair(obsWithin, predWithin, predAcross, runWithinNum, runAcrossNum, wedgeAngle, projectName, roiName, nSubj, subjectMode, figureDir, onsets, shadeColorsAll, stimdur_s, tr_s, sharedYlim)
% observed = solid black, within fit = solid dark red, across fit =
% dashed light grey (agreed 2026-09-01); condition-colored, 30%-opacity
% shaded patches behind the lines (grey=moving, cyan/yellow/dark-green/
% magenta=the 4 stationary conditions by raw code, blank unshaded).

    r2within = centeredR2(obsWithin, predWithin);
    r2across = centeredR2(obsWithin, predAcross);

    dataColor = [0 0 0];
    withinColor = [0.85 0.1 0.1];
    acrossColor = [0.75 0.75 0.75]; % light grey, was blue

    TRs = 0:(numel(obsWithin)-1);
    stimdur_TR = stimdur_s / tr_s;

    fig = figure('Visible','off');
    hold on

    for k = 1:numel(onsets)
        theseOnsets = onsets{k};
        color = shadeColorsAll{k};
        for oi = 1:numel(theseOnsets)
            onsetT = theseOnsets(oi) - 1; % 0-indexed TR, matching the TRs axis above
            offsetT = min(onsetT + stimdur_TR, TRs(end));
            patch([onsetT onsetT offsetT offsetT], ...
                [-1000 1000 1000 -1000], color, ...
                'EdgeColor', 'none', 'FaceAlpha', 0.15, 'HandleVisibility', 'off');
        end
    end

    plot(TRs, obsWithin, '-', 'LineWidth', 1.5, 'Color', dataColor, 'DisplayName', sprintf('Run %d observed', runWithinNum));
    plot(TRs, predWithin, '-', 'LineWidth', 1.5, 'Color', withinColor, 'DisplayName', sprintf('Run %d fit (within), R^2=%.1f%%', runWithinNum, r2within));
    plot(TRs, predAcross, '--', 'LineWidth', 1.5, 'Color', acrossColor, 'DisplayName', sprintf('Run %d fit (across), R^2=%.1f%%', runAcrossNum, r2across));

    ylim(sharedYlim);
    xlim([TRs(1), TRs(end)]);
    box on;

    hold off
    xlabel('TR'); ylabel('% signal change');
    legend('Location','best');
    title(sprintf('%s (%s, n=%d) / %s / %d^\\circ: Run %d data vs. Run %d (within) and Run %d (across) fits', ...
        projectName, subjectMode, nSubj, roiName, wedgeAngle, runWithinNum, runWithinNum, runAcrossNum), 'Interpreter','tex');

    % Axes (plot area) sized to exactly 15.37 x 3.37 cm, EXCLUDING the
    % title -- the figure canvas is that plus fixed margins for the axis
    % labels/ticks (left/bottom) and the title's own space (top), so the
    % title sits above the 3.37cm-tall plot area rather than eating into it.
    axesWidth_cm = 15.37;
    axesHeight_cm = 3.37;
    leftMargin_cm = 1.5;
    rightMargin_cm = 0.5;
    bottomMargin_cm = 1.2;
    topMarginForTitle_cm = 1.2;
    figWidth_cm = axesWidth_cm + leftMargin_cm + rightMargin_cm;
    figHeight_cm = axesHeight_cm + bottomMargin_cm + topMarginForTitle_cm;

    ax = gca;
    fig.Units = 'centimeters';
    fig.Position = [1, 1, figWidth_cm, figHeight_cm];
    fig.PaperUnits = 'centimeters';
    fig.PaperPositionMode = 'manual';
    fig.PaperSize = [figWidth_cm, figHeight_cm];
    fig.PaperPosition = [0, 0, figWidth_cm, figHeight_cm];
    ax.Units = 'centimeters';
    ax.Position = [leftMargin_cm, bottomMargin_cm, axesWidth_cm, axesHeight_cm];

    filename = fullfile(figureDir, sprintf('%s_%s_%ddeg_run%d_within_run%d_across.pdf', projectName, roiName, wedgeAngle, runWithinNum, runAcrossNum));
    exportgraphics(fig, filename, 'ContentType', 'vector');
    fprintf('Saved: %s (R^2 within=%.2f%%, across=%.2f%%)\n', filename, r2within, r2across);
    close(fig);
end

function r2 = centeredR2(obs, pred)
    ssRes = sum((obs - pred).^2);
    ssTot = sum((obs - mean(obs)).^2);
    r2 = (1 - ssRes/ssTot) * 100;
end
