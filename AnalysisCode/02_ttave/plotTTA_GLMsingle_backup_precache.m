function plotTTA_GLMsingle(projectName, roiName, targetWedgeAngle, varargin)
% PLOTTTA_GLMSINGLE  Trial-triggered average (TTA), average-observer, for
% the 4 stationary/orientation conditions, blank-subtracted, both observed
% and GLMsingle-predicted. Prints two R^2 summaries and plots ONE location.
%
%   plotTTA_GLMsingle('dg', 'V1', 90)
%   plotTTA_GLMsingle('da', 'V1', 90)
%
% Prints (agreed 2026-09-01, replacing the derangement-test idea -- with
% only 4 conditions there are just !4=9 possible mismatches, capping any
% exact permutation p-value at 1/10 regardless of how much data feeds
% each draw, so a plain concatenated R^2 summary was used instead):
%   (1) ONE summary R^2 for the whole concatenated vector -- every one of
%       the 8 canonical polar-angle locations x all 4 conditions, group-
%       average observed and CORRECTLY-matched predicted traces
%       concatenated together, one R^2 per experiment. This needs all 8
%       locations' TTA computed even though only targetWedgeAngle is
%       plotted.
%   (2) The R^2 for each of the 4 conditions separately, AT
%       targetWedgeAngle only (matching the plotted location).
%
% METHOD (agreed 2026-09-01):
% Per subject, per trial, extract a peristimulus epoch (5 TRs before to 19
% TRs after each onset, 25 timepoints) for both the observed and
% GLMsingle-predicted continuous time series -- reusing the SAME
% corrected reconstruction as computeSubjectRunTraces_GLMsingle.m/
% runRawTimeseriesPermutationTest.m (NOT createTTaveTable.m's
% computeConditionTTA_rawdata/_modelfit, which have the same %sc-before-
% denoising order bug found and fixed 2026-09-01 -- see that file's
% header comment). Average across trials (pooled over that subject's
% runs) and across the included vertices (ROI + this wedge + ecc[4,8] +
% R^2>=0.1) to get one 25-timepoint trace per condition (4 stationary +
% blank) per subject, separately for observed and predicted. Subtract the
% blank trace from each of the 4 orientation traces (per subject, on the
% already-averaged traces -- same convention as
% plotOrientationConditionsAllVoxels.m). Average across subjects (13 for
% dg INCLUDING sub-0395, 7 for da excluding sub-0395 -- verified
% 2026-09-01 that dg's 13-subject literal here matches the one already
% confirmed correct in runRawTimeseriesPermutationTest.m/
% plotRawTimeseriesWithinAcross_GLMsingle.m) to get the "average observer"
% trace per condition.
%
% Colors match plotRawTimeseriesWithinAcross_GLMsingle.m exactly, by raw
% stimulus code (project-invariant): 0=cyan, 90=yellow, 45=dark green,
% 135=magenta. Observed = markers, predicted = solid line, both in that
% condition's color (matching plotOrientationConditionsAllVoxels.m's
% established convention). R^2 (centered, same formula as the raw-
% timeseries stat) reported per condition in the legend.
%
% Name-value options: 'subjectMode', 'nRunsCheck', 'bidsDir', 'githubDir',
% 'figureDir' -- as in plotRawTimeseriesWithinAcross_GLMsingle.m.

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
eventTRs_prior = 5;
eventTRs_after = 20;
nTimepoints = eventTRs_prior + eventTRs_after; % 25
TRvals = (-eventTRs_prior):(eventTRs_after-1);

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
    'sub-0397', 'sub-0427'}; % sub-0395 INCLUDED here -- verified 2026-09-01, this matches the confirmed-correct dg=13 literal elsewhere
matched_7 = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
    'sub-0426', 'sub-0250'}; % sub-0395 excluded (mismatched da pilot stimulus)

if strcmp(projectName, 'dg')
    if strcmp(opt.subjectMode, 'full')
        subjects = dg_subjects_13;
    else
        subjects = matched_7;
    end
elseif strcmp(projectName, 'da')
    subjects = matched_7;
else
    error('plotTTA_GLMsingle:project', 'projectName must be ''dg'' or ''da''.');
end
nSubj = numel(subjects);

projectSettingsBase = loadConfig(githubDir);
jsonParams = jsondecode(fileread('setup.json'));
stimdur_s = jsonParams.stimdur_s.Val;
tr_s = jsonParams.tr_s.Val;

STATIC_COLORS = {[0 1 1], [1 1 0], [0 0.5 0], [1 0 1]}; % cyan, yellow, dark green, magenta, for raw codes [0,90,45,135] -- same as plotRawTimeseriesWithinAcross_GLMsingle.m
condLabels_dg = {'Horizontal','Vertical','Right-leaning','Left-leaning'};
condLabels_da = {'Annulus (tangential)','Pinwheel (radial)','CW spiral','CCW spiral'};
if strcmp(projectName,'dg')
    condLabels = condLabels_dg;
else
    condLabels = condLabels_da;
end

%% per-subject: TTA epochs for the 4 orientation conditions + blank, observed & predicted, ALL 8 locations

wedgeAngles = 0:45:315;
nWedges = numel(wedgeAngles);
targetW = find(wedgeAngles == targetWedgeAngle);
if isempty(targetW)
    error('plotTTA_GLMsingle:wedge', '%d is not one of the 8 canonical wedge angles.', targetWedgeAngle);
end

subjTTA_obs = nan(nSubj, nWedges, 4, nTimepoints); % blank-subtracted, per subject
subjTTA_pred = nan(nSubj, nWedges, 4, nTimepoints);

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

    fprintf('  computing GLMsingle-based TTA (all 8 locations): %s / %s ...\n', projectName, subj);
    st = computeSubjectRunTraces_GLMsingle(bidsDir, projectName, subj, ses, projectSettingsBase, ...
        stimdur_s, tr_s, rawS.matrices_onset, nRunsCheck, roiName);

    excludedRunIdx = findDuplicateDesignRuns(rawS.matrices_onset, nRunsCheck);
    nRuns = numel(st.observed{1});

    for w = 1:nWedges
        % accumulate per-condition epoch sums/counts across all
        % (non-excluded) runs -- condition columns: 9-12 static (raw codes
        % [0,90,45,135] in that order), 13 blank
        sumObs = zeros(5, nTimepoints); countObs = zeros(5, nTimepoints);
        sumPred = zeros(5, nTimepoints); countPred = zeros(5, nTimepoints);

        for r = 1:nRuns
            if ismember(r, excludedRunIdx)
                continue
            end
            obsTrace = st.observed{w}{r};
            predTrace = st.predicted{w}{r};
            if all(isnan(obsTrace)) || all(isnan(predTrace))
                continue
            end
            designMatrix = rawS.matrices_onset{r};
            nTRsRun = numel(obsTrace);

            for cc = 1:5
                ci = cc + 8; % columns 9-13
                if ci > size(designMatrix,2)
                    continue
                end
                onsets = find(designMatrix(:,ci) == 1);
                for oi = 1:numel(onsets)
                    startIdx = onsets(oi) - eventTRs_prior;
                    endIdx = onsets(oi) + eventTRs_after - 1;
                    validRange = max(startIdx,1):min(endIdx,nTRsRun);
                    if isempty(validRange), continue; end
                    relIdx = (validRange(1)-startIdx+1):(validRange(1)-startIdx+numel(validRange));

                    sumObs(cc,relIdx) = sumObs(cc,relIdx) + obsTrace(validRange);
                    countObs(cc,relIdx) = countObs(cc,relIdx) + 1;
                    sumPred(cc,relIdx) = sumPred(cc,relIdx) + predTrace(validRange);
                    countPred(cc,relIdx) = countPred(cc,relIdx) + 1;
                end
            end
        end

        condTTA_obs = sumObs ./ countObs;   % 5 x nTimepoints (NaN where countObs==0)
        condTTA_pred = sumPred ./ countPred;

        % blank-subtract (row 5 = blank) and store the 4 orientation conditions
        subjTTA_obs(ss,w,:,:) = condTTA_obs(1:4,:) - condTTA_obs(5,:);
        subjTTA_pred(ss,w,:,:) = condTTA_pred(1:4,:) - condTTA_pred(5,:);
    end
end

%% group average across subjects

groupTTA_obs = squeeze(nanmean(subjTTA_obs, 1));   % nWedges x 4 x nTimepoints
groupTTA_pred = squeeze(nanmean(subjTTA_pred, 1));

%% (1) summary R^2: concatenate ALL 8 locations x 4 conditions, correctly matched

obsConcatAll = reshape(permute(groupTTA_obs, [3 2 1]), [], 1);   % nTimepoints*4*nWedges x 1
predConcatAll = reshape(permute(groupTTA_pred, [3 2 1]), [], 1);
summaryR2 = centeredR2(obsConcatAll, predConcatAll);

fprintf('\n%s (%s, n=%d) / %s: SUMMARY R^2 (all 8 locations x 4 conditions concatenated) = %.2f%%\n', ...
    projectName, opt.subjectMode, nSubj, roiName, summaryR2);

%% (2) per-condition R^2, AT targetWedgeAngle only

r2vals = nan(1,4);
fprintf('%s (%s, n=%d) / %s / %d^\\circ (plotted location): per-condition R^2:\n', ...
    projectName, opt.subjectMode, nSubj, roiName, targetWedgeAngle);
for cc = 1:4
    obsVals = squeeze(groupTTA_obs(targetW,cc,:))';
    predVals = squeeze(groupTTA_pred(targetW,cc,:))';
    r2vals(cc) = centeredR2(obsVals, predVals);
    fprintf('  %s: R^2 = %.2f%%\n', condLabels{cc}, r2vals(cc));
end

%% plot (targetWedgeAngle only)

fig = figure('Visible','off');
hold on
legendHandles = gobjects(4,1);
legendLabels = strings(4,1);

for cc = 1:4
    col = STATIC_COLORS{cc};
    obsVals = squeeze(groupTTA_obs(targetW,cc,:))';
    predVals = squeeze(groupTTA_pred(targetW,cc,:))';

    plot(TRvals, obsVals, 'o', 'LineWidth', 1.5, 'Color', col, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'w', 'HandleVisibility','off');
    h = plot(TRvals, predVals, '-', 'LineWidth', 1.5, 'Color', col);
    legendHandles(cc) = h;
    legendLabels(cc) = sprintf('%s, R^2=%.1f%%', condLabels{cc}, r2vals(cc));
end

yline(0, 'k:', 'HandleVisibility', 'off');
xline(0, 'k:', 'HandleVisibility', 'off');
box on;
hold off

xlabel('TR relative to stimulus onset'); ylabel('% signal change');
legend(legendHandles, legendLabels, 'Location','best','Interpreter','none');
title(sprintf('%s (%s, n=%d) / %s / %d^\\circ: trial-triggered average (blank-subtracted)', ...
    projectName, opt.subjectMode, nSubj, roiName, targetWedgeAngle), 'Interpreter','tex');

axesWidth_cm = 15.37/2; axesHeight_cm = 3.37; % width halved per request, height unchanged
leftMargin_cm = 1.5; rightMargin_cm = 0.5; bottomMargin_cm = 1.2; topMarginForTitle_cm = 1.2;
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

filename = fullfile(figureDir, sprintf('TTA_%s_%s_%ddeg.pdf', projectName, roiName, targetWedgeAngle));
exportgraphics(fig, filename, 'ContentType', 'vector');
fprintf('Saved: %s\n', filename);
close(fig);

end

function r2 = centeredR2(obs, pred)
    valid = isfinite(obs) & isfinite(pred);
    ssRes = sum((obs(valid) - pred(valid)).^2);
    ssTot = sum((obs(valid) - mean(obs(valid))).^2);
    r2 = (1 - ssRes/ssTot) * 100;
end
