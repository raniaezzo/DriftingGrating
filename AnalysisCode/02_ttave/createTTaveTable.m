clc; clear all; close all;

% createTTaveTable.m
%
% Builds two long-format, per-vertex tables of voxelwise trial-triggered
% averages (TTA) across all 13 conditions x 2 experiments (dg, da):
%   - ttaveTable_raw.mat    : TTA computed from the raw preprocessed BOLD data
%   - ttaveTable_model.mat  : TTA computed from the GLMsingle single-trial betas
%                             (HRF-convolved model reconstruction)
%
% Row = one vertex at one peristimulus TR. Columns mirror
% 01_process_singlesubjectGLM/createTables.m's allsubjectsTable.csv
% (subject, visual_area, pRF_angle_bin, pRF_angle, pRF_ecc, pRF_r2, pRF_sigma,
% included), plus cartexp_R2/polexp_R2 (the GLMsingle GLM fit R^2 per
% project, distinct from pRF_r2 which is the retinotopy fit), a TR column,
% and the 26 condition columns (13 dg + 13 da). Each condition column holds
% the scalar TTA value at that TR rather than a single GLM beta. Rows are
% restricted to vertices within one of the 10 named ROIs (unlike
% allsubjectsTable.csv, which keeps the full surface) to keep the ~25x
% row-count increase from adding TR tractable.

%% SET UP

bidsDir = '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
githubDir = '~/Documents/GitHub';
addpath(genpath(fullfile(githubDir, 'atlasmgz')));
setup_user('rania', bidsDir)

savedir = fullfile(bidsDir, 'derivatives', 'summaryTables');
if ~isfolder(savedir)
    mkdir(savedir)
end

subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
    'sub-0395', 'sub-0426', 'sub-0250'};

projects = {'dg', 'da'};
hRF_setting = 'glmsingle';

hemis = {'lh'; 'rh'};
roinames = {'V1', 'V2', 'V3', 'hV4', 'V3a', 'V3b', 'hMTcomplex', 'pMT', 'pMST', 'V2d'};

dg_stimNames = {'cartexp_vertical_grating_rightwards_motion', 'cartexp_horizontal_grating_upwards_motion', ...
    'cartexp_vertical_grating_leftwards_motion', 'cartexp_horizontal_grating_downwards_motion', ...
    'cartexp_leftleaning_grating_upperrightwards_motion', 'cartexp_rightleaning_grating_upperleftwards_motion', ...
    'cartexp_leftleaning_grating_lowerleftwards_motion', 'cartexp_rightleaning_grating_lowerrightwards_motion', ...
    'cartexp_horizontal_stationary', 'cartexp_vertical_stationary', 'cartexp_rightleaning_grating_stationary', ...
    'cartexp_leftleaning_grating_stationary', 'cartexp_blank'};

da_stimNames = {'polexp_pinwheel_grating_clockwise_motion', 'polexp_annulus_grating_outwards_motion', ...
    'polexp_pinwheel_grating_cclockwise_motion', 'polexp_annulus_grating_inwards_motion', ...
    'polexp_ccspiral_grating_clockoutwards_motion', 'polexp_cspiral_grating_cclockoutwards_motion', ...
    'polexp_ccspiral_grating_cclockinwards_motion', 'polexp_cspiral_grating_clockinwards_motion', ...
    'polexp_annulus_grating_stationary', 'polexp_pinwheel_grating_stationary', 'polexp_cspiral_grating_stationary', ...
    'polexp_ccspiral_grating_stationary', 'polexp_blank'};

nConditions = 13;
eventTRs_prior = 5;
eventTRs_after = 20;
nTimepoints = eventTRs_prior + eventTRs_after; % 25
TRvals = (-eventTRs_prior):(eventTRs_after-1);
tr_s = 1;
stimdur_s = 3;

allRawTables = cell(numel(subjects), 1);
allModelTables = cell(numel(subjects), 1);

% If true and the saved tables already exist on disk, load them instead of
% recomputing (recomputing reloads raw BOLD data and GLMsingle output for
% all 8 subjects x 2 experiments and is slow; loading is fast). Set to
% false to force recomputation even if the files exist.
loadExistingTables = true;
rawTablePath = fullfile(savedir, 'ttaveTable_raw.mat');
modelTablePath = fullfile(savedir, 'ttaveTable_model.mat');
canLoadExisting = loadExistingTables && isfile(rawTablePath) && isfile(modelTablePath);

if canLoadExisting

    %% LOAD EXISTING TABLES (instead of recomputing)

    fprintf('Loading existing tables from %s (set loadExistingTables=false above to force recomputation)...\n', savedir);

    loaded = load(rawTablePath, 'ttaveTable_raw');
    ttaveTable_raw = loaded.ttaveTable_raw;
    loaded = load(modelTablePath, 'ttaveTable_model');
    ttaveTable_model = loaded.ttaveTable_model;

    for si = 1:numel(subjects)
        allRawTables{si} = ttaveTable_raw(ttaveTable_raw.subject == subjects{si}, :);
        allModelTables{si} = ttaveTable_model(ttaveTable_model.subject == subjects{si}, :);
        if isempty(allRawTables{si})
            warning('%s not found in the loaded tables -- treated as unprocessed by the plotting functions below.', subjects{si});
        end
    end

    fprintf('Loaded %d rows (raw) / %d rows (model) for %d subjects.\n', ...
        height(ttaveTable_raw), height(ttaveTable_model), numel(subjects));

else

%% PER SUBJECT

for si = 1:numel(subjects)

    subjectname = subjects{si};
    fprintf('=== %s ===\n', subjectname);

    %% pRF metadata (same source/logic as createTables.m)

    retDir = dir(fullfile(bidsDir, 'derivatives', 'prfvista_mov', subjectname, '**/stimfiles.mat'));
    retDir = retDir.folder;

    for hi = 1:numel(hemis)
        hemi = hemis{hi};
        ret_moving.(sprintf('%s_pa', hemi)) = MRIread(fullfile(retDir, sprintf('%s.angle_adj.mgz', hemi)));
        ret_moving.(sprintf('%s_ecc', hemi)) = MRIread(fullfile(retDir, sprintf('%s.eccen.mgz', hemi)));
        ret_moving.(sprintf('%s_vexp', hemi)) = MRIread(fullfile(retDir, sprintf('%s.vexpl.mgz', hemi)));
        ret_moving.(sprintf('%s_sigma', hemi)) = MRIread(fullfile(retDir, sprintf('%s.sigma.mgz', hemi)));
    end

    hSize = [size(ret_moving.lh_pa.vol,2), size(ret_moving.rh_pa.vol,2)];
    N = sum(hSize);

    Tmeta = table( ...
        strings(N,1), ...    % subject
        strings(N,1), ...    % visual_area
        nan(N,1), ...        % pRF_angle_bin
        nan(N,1), ...        % pRF_angle
        nan(N,1), ...        % pRF_ecc
        nan(N,1), ...        % pRF_r2
        nan(N,1), ...        % pRF_sigma
        false(N,1), ...      % included
        'VariableNames', {'subject','visual_area', 'pRF_angle_bin', 'pRF_angle','pRF_ecc', ...
        'pRF_r2','pRF_sigma', 'included'});

    Tmeta.subject = repmat(string(subjectname), N, 1);
    Tmeta.pRF_angle = map_theta([ret_moving.lh_pa.vol, ret_moving.rh_pa.vol])';
    Tmeta.pRF_ecc = [ret_moving.lh_ecc.vol, ret_moving.rh_ecc.vol]';
    Tmeta.pRF_r2 = [ret_moving.lh_vexp.vol, ret_moving.rh_vexp.vol]';
    Tmeta.pRF_sigma = [ret_moving.lh_sigma.vol, ret_moving.rh_sigma.vol]';

    binCenters = 0:45:315;
    angles = Tmeta.pRF_angle;
    circDist = abs(mod(angles - binCenters + 180, 360) - 180);
    [minDist, idx] = min(circDist, [], 2);
    Tmeta.pRF_angle_bin = nan(height(Tmeta),1);
    valid = minDist <= 22.5;
    Tmeta.pRF_angle_bin(valid) = binCenters(idx(valid));

    Tmeta.included = ...
        ~isnan(Tmeta.pRF_angle_bin) & ...
        Tmeta.pRF_ecc >= 4 & Tmeta.pRF_ecc <= 8 & ...
        Tmeta.pRF_r2 >= 0.1;

    for ri = 1:length(roinames)
        roiname = roinames{ri};

        lh_label = read_label(subjectname, sprintf('retinotopy_RE/lh.%s_REmanual', roiname));
        rh_label = read_label(subjectname, sprintf('retinotopy_RE/rh.%s_REmanual', roiname));

        label_idx = [lh_label(:,1)+1 ; rh_label(:,1)+hSize(1)+1];

        if strcmp(roiname, 'V2d')
            label_idx = lh_label(:,1)+1;
            roiname = strcat('left_', roiname);
        end

        if strcmp(roiname, 'pMT')
            roiname = 'MT';
        elseif strcmp(roiname, 'pMST')
            roiname = 'MST';
        end

        Tmeta.visual_area(label_idx) = roiname;
    end

    % restrict all downstream processing to vertices inside a named ROI
    roiVertexIdx = find(Tmeta.visual_area ~= "");
    Tmeta = Tmeta(roiVertexIdx, :);
    nSel = numel(roiVertexIdx);
    fprintf('  %d / %d vertices fall within a named ROI\n', nSel, N);

    %% per project: raw-data TTA and model-fit TTA

    condTTA_raw = struct();
    condTTA_model = struct();
    glmR2 = struct();

    for pi = 1:numel(projects)
        projectName = projects{pi};

        subjectDir = fullfile(bidsDir, 'derivatives', sprintf('%sGLM', projectName), ...
            sprintf('hRF_%s', hRF_setting), subjectname);
        ses = findses(subjectDir);
        derivativesFolder = fullfile(subjectDir, ses);

        % load GLMsingle output once per project and reuse it for both the
        % model-fit reconstruction and the raw-data %-signal-change
        % normalization, so both use the same per-vertex mean-signal
        % reference (meanvol) instead of two independently-computed means
        % that can drift apart between the dg and da sessions
        modelData = load(fullfile(derivativesFolder, 'modelOutput.mat'));
        results = modelData.modelOut{1,4};
        designInfo = modelData.designSINGLE;

        meanvol = results.meanvol(roiVertexIdx,1);
        glmR2.(projectName) = results.R2(roiVertexIdx,1);

        fprintf('  [%s] computing raw-data TTA...\n', projectName);
        condTTA_raw.(projectName) = computeConditionTTA_rawdata( ...
            bidsDir, projectName, subjectname, ses, roiVertexIdx, meanvol, tr_s, eventTRs_prior, eventTRs_after);

        fprintf('  [%s] computing model-fit TTA...\n', projectName);
        condTTA_model.(projectName) = computeConditionTTA_modelfit( ...
            results, designInfo, tr_s, stimdur_s, roiVertexIdx, eventTRs_prior, eventTRs_after);
    end

    % GLMsingle's overall model R^2 per vertex per project (distinct from
    % pRF_r2, which is the retinotopy fit) -- kept alongside pRF_r2 for
    % later use, and used in place of it when picking/plotting example
    % voxels below, since it reflects how well the 13-condition GLM
    % explains this experiment's data specifically.
    Tmeta.cartexp_R2 = glmR2.dg;
    Tmeta.polexp_R2 = glmR2.da;

    %% assemble long-format (one row per vertex per TR) tables for this subject

    allRawTables{si} = buildLongTable(Tmeta, condTTA_raw.dg, condTTA_raw.da, ...
        dg_stimNames, da_stimNames, TRvals, nTimepoints);

    allModelTables{si} = buildLongTable(Tmeta, condTTA_model.dg, condTTA_model.da, ...
        dg_stimNames, da_stimNames, TRvals, nTimepoints);

end

%% SAVE

ttaveTable_raw = vertcat(allRawTables{:});
save(fullfile(savedir, 'ttaveTable_raw.mat'), 'ttaveTable_raw', '-v7.3');

ttaveTable_model = vertcat(allModelTables{:});
save(fullfile(savedir, 'ttaveTable_model.mat'), 'ttaveTable_model', '-v7.3');

fprintf('Done. Saved %d rows to ttaveTable_raw.mat and ttaveTable_model.mat in %s\n', ...
    height(ttaveTable_raw), savedir);

end

%% OPTIONAL: EXAMPLE PLOT
% For one subject/ROI (included==1 vertices only), plots the "stationary"
% (orientation, no-motion) response per polar angle bin -- raw (dashed) vs.
% model (solid), one color per bin -- for both experiments side by side.
% aggregationMode selects how each bin's curve is formed:
%   'peakR2'  - the single highest-GLM-R2 vertex in that bin (project's own
%               cartexp_R2/polexp_R2, not pRF_r2)
%   'meanROI' - the mean across all included vertices in that bin
% Uses allRawTables/allModelTables (the per-subject cell arrays built above)
% directly, no need to reload the saved .mat files.
%
% exampleSubject can be set to 'all' instead of a specific subject ID --
% both plotExampleVoxels and plotOrientationConditionsAllVoxels below will
% then run their exact same per-subject steps for every subject processed
% in this run and average the resulting curves across subjects before
% plotting.
%
% offsetRawToBaseline/offsetModelToBaseline (applied independently, one per
% curve type) subtract each curve's own pre-stimulus baseline mean (TR
% -5:-1) from the whole curve, so that baseline period is 0-centered. Set
% either to false to leave that curve type unshifted.
%
% savePlotsAsPDF saves each figure as a vector PDF into figSaveDir.

makeExamplePlot = true;
exampleSubject = 'sub-wlsubj124'; %'sub-wlsubj124'; % a specific subject ID, or 'all'
exampleROI = 'V1';
exampleAggregation = 'peakR2'; % 'peakR2' or 'meanROI'
offsetRawToBaseline = false;
offsetModelToBaseline = false;
savePlotsAsPDF = true;
figSaveDir = fullfile(savedir, 'figures');

if makeExamplePlot
    plotExampleVoxels(allRawTables, allModelTables, subjects, exampleSubject, exampleROI, exampleAggregation, ...
        offsetRawToBaseline, offsetModelToBaseline, savePlotsAsPDF, figSaveDir);
end
%%
% Per-condition version: instead of combining the 4 orientation conditions
% into one mean (blank-subtracted) curve per polar angle bin, this plots
% each of the 4 orientation conditions as its own blank-subtracted curve,
% averaged across all valid (included==1) voxels in the ROI -- no
% polar-angle grouping. *_blank itself is excluded from the 4 plotted
% conditions (it's the baseline being subtracted, not an orientation
% condition).
makeConditionPlot = true;

if makeConditionPlot
    plotOrientationConditionsAllVoxels(allRawTables, allModelTables, subjects, exampleSubject, exampleROI, ...
        offsetRawToBaseline, offsetModelToBaseline, savePlotsAsPDF, figSaveDir);
end
%%
% Across-subjects version: averages the 4 orientation conditions
% (blank-subtracted) across all valid (included==1) V1 voxels per subject,
% then plots each subject as its own line (rather than one subject across
% voxel/angle-bin groupings) -- again as two subplots, one per experiment.
% Only subjects actually processed in allRawTables/allModelTables above
% will appear.
makeSubjectComparisonPlot = true;

if makeSubjectComparisonPlot
    plotOrientationConditionsAllSubjects(allRawTables, allModelTables, subjects, exampleROI, ...
        offsetRawToBaseline, offsetModelToBaseline, savePlotsAsPDF, figSaveDir);
end
%%
% Single-bin, two-condition version: restricted to one specified
% pRF_angle_bin (instead of all 8) and two specified conditions per
% experiment (instead of averaging all 4) -- still one subplot per
% experiment (dg/cartexp_ left, da/polexp_ right). Within each subplot, one
% line per condition (not per subject); exampleSubject controls whether
% subjects are averaged together ('all') or restricted to one subject's
% data, same as plotOrientationConditionsAllVoxels.
makeOneBinPlot = true;
exampleSubject = 'all'; % a specific subject ID, or 'all'
exampleROI = 'V1';
targetBin = 90; % 0 deg = right horizontal meridian
cartConditions = {'cartexp_horizontal_stationary', 'cartexp_vertical_stationary', ...
    'cartexp_rightleaning_grating_stationary', 'cartexp_leftleaning_grating_stationary'}; % all 4 dg stationary conditions
polConditions = {'polexp_annulus_grating_stationary', 'polexp_pinwheel_grating_stationary', ...
    'polexp_cspiral_grating_stationary', 'polexp_ccspiral_grating_stationary'}; % all 4 da stationary conditions

if makeOneBinPlot
    plotOrientationConditionsOneBin(allRawTables, allModelTables, subjects, exampleSubject, exampleROI, targetBin, ...
        cartConditions, polConditions, ...
        offsetRawToBaseline, offsetModelToBaseline, savePlotsAsPDF, figSaveDir);
end

%% ===================== LOCAL FUNCTIONS =====================

function sesFolder = findses(subjdir)
    d = dir(fullfile(subjdir, 'ses-*'));
    d = d([d.isdir]);
    if isempty(d)
        error('No ses-* directory found in %s', subjdir);
    elseif numel(d) > 1
        error('Multiple ses-* directories found in %s', subjdir);
    end
    sesFolder = d.name;
end

function condTTA = computeConditionTTA_rawdata(bidsDir, projectName, subj, ses, roiVertexIdx, meanvol, tr_s, eventTRs_prior, eventTRs_after)
% Voxelwise trial-triggered average per condition, computed from the raw
% preprocessed BOLD time series (%-signal-change). Mirrors the fititer==1
% branch of ttave_compute.m, but keeps every vertex instead of collapsing
% the ROI with nanmedian. %-signal-change uses GLMsingle's own per-vertex
% mean-signal reference (meanvol, passed in by the caller from that
% project's modelOutput.mat) rather than a separately-computed per-run
% mean, so it is on the same scale as the model-fit reconstruction (which
% is built from GLMsingle's own %-BOLD betas, themselves scaled by this
% same meanvol).
%
% Before converting to %-signal-change, each run is polynomial-denoised
% the same way GLMestimatesingletrial.m denoises its own noise-pool data
% (see its line ~1089: projectionmatrix(constructpolynomialmatrix(n,
% 0:maxpolydeg))*data): a per-run polynomial nuisance basis is projected
% out to remove low-frequency drift, using the same default polynomial
% degree GLMsingle itself uses (round(L/2), L = run duration in minutes --
% GLMestimatesingletrial.m line ~525). This keeps the raw-data values that
% end up in allRawTables/ttaveTable_raw on the same denoising basis as the
% betas the model-fit table is built from. constructpolynomialmatrix and
% projectionmatrix are GLMsingle utilities (matlab/utilities/), assumed to
% already be on the path alongside GLMpredictresponses/
% getcanonicalhrflibrary.

    nConditions = 13;
    nTimepoints = eventTRs_prior + eventTRs_after;
    nSel = numel(roiVertexIdx);

    derivativesFolder = fullfile(bidsDir, 'derivatives', sprintf('%sGLM', projectName), ...
        'hRF_glmsingle', subj, ses);

    rawInfo = load(fullfile(derivativesFolder, 'rawInfo.mat'), 'matrices_onset');
    matrices_onset = rawInfo.matrices_onset;
    nRuns = numel(matrices_onset);

    datafiles = load_data(bidsDir, projectName, 'fsnative', '.mgh', subj, ses, 1:nRuns);

    sumTTA = zeros(nConditions, nSel, nTimepoints);
    countTTA = zeros(nConditions, nSel, nTimepoints);

    for r = 1:nRuns
        df = datafiles{r}(roiVertexIdx, :); % nSel x nTRsRun (voxels x time)

        % polynomial-denoise: project out a run-specific low-degree
        % polynomial basis (DC, linear, ...) along the time dimension
        nTRsRun = size(df, 2);
        maxpolydeg = round(((nTRsRun * tr_s) / 60) / 2);
        pmatrix = constructpolynomialmatrix(nTRsRun, 0:maxpolydeg);
        polymatrix = projectionmatrix(pmatrix);
        df = (polymatrix * df')'; % denoise in time x voxels, then back to voxels x time

        timeseries_psc = ((df ./ meanvol) - 1) * 100;

        onsetsByCondition = cell(1, nConditions);
        for ci = 1:nConditions
            onsetsByCondition{ci} = find(matrices_onset{r}(:,ci) == 1);
        end

        [sumTTA, countTTA] = accumulateConditionEpochs(sumTTA, countTTA, timeseries_psc, ...
            onsetsByCondition, eventTRs_prior, eventTRs_after, 1:nSel);
    end

    condTTA = sumTTA ./ countTTA;
end

function condTTA = computeConditionTTA_modelfit(results, designInfo, tr_s, stimdur_s, roiVertexIdx, eventTRs_prior, eventTRs_after)
% Voxelwise trial-triggered average per condition, reconstructed from the
% GLMsingle single-trial betas via HRF convolution (GLMpredictresponses),
% batched by HRF-library index so the convolution is shared across all
% vertices with the same HRF rather than looped per vertex (c.f.
% ttave_computeGLMsingle.m). Reconstruction is done and immediately reduced
% to per-condition epochs per HRF group -- a full nVertices x totalTRs
% reconstructed time series is never held in memory. results/designInfo are
% pre-loaded by the caller from that project's modelOutput.mat (results =
% modelOut{1,4}, designInfo = designSINGLE) so the file is read only once
% per project rather than once per helper function.

    nConditions = 13;
    nTimepoints = eventTRs_prior + eventTRs_after;
    nSel = numel(roiVertexIdx);

    modelmd = results.modelmd;             % nVertFull x 1 x 1 x nTrials, in %
    HRFindex = results.HRFindex;           % nVertFull x 1
    nTrials = size(modelmd, 4);

    hrflibrary = getcanonicalhrflibrary(stimdur_s, tr_s); % nHRFs x time

    designRuns = designInfo.designSINGLE;  % 1 x nRuns cell, TR x nTrials
    stimorder = designInfo.stimorder;      % 1 x nTrials
    nRuns = numel(designRuns);
    numtimepoints = cellfun(@(x) size(x,1), designRuns);

    % onset TRs per run per condition (trial-level design -> condition via stimorder)
    onsetsByRunCondition = cell(nRuns, nConditions);
    for r = 1:nRuns
        for ci = 1:nConditions
            trialCols = find(stimorder == ci);
            conDesign = designRuns{r}(:, trialCols);
            onsetsByRunCondition{r, ci} = find(any(conDesign, 2));
        end
    end

    sumTTA = zeros(nConditions, nSel, nTimepoints);
    countTTA = zeros(nConditions, nSel, nTimepoints);

    hrfii_sel = HRFindex(roiVertexIdx);
    uniqueHRFs = unique(hrfii_sel)';

    for hh = uniqueHRFs
        selPos = find(hrfii_sel == hh);      % rows (1..nSel) in the accumulator
        voxFull = roiVertexIdx(selPos);      % absolute vertex ids into modelmd

        betas = reshape(modelmd(voxFull,1,1,:), numel(voxFull), nTrials); % nVoxSubset x nTrials
        hrf = hrflibrary(hh,:)';             % time x 1

        mf = GLMpredictresponses({hrf, betas}, designRuns, tr_s, numtimepoints, 1); % 1 x nRuns cell

        for r = 1:nRuns
            signalMatrix = mf{r};
            % GLMpredictresponses is expected to return one row per input beta
            % row (i.e. per voxFull vertex). If it instead returns full-surface
            % rows (as in plottingGLMsingle.m's usage), fall back to indexing.
            if size(signalMatrix, 1) ~= numel(voxFull)
                signalMatrix = signalMatrix(voxFull, :);
            end

            [sumTTA, countTTA] = accumulateConditionEpochs(sumTTA, countTTA, signalMatrix, ...
                onsetsByRunCondition(r, :), eventTRs_prior, eventTRs_after, selPos);
        end

        clear mf
    end

    condTTA = sumTTA ./ countTTA;
end

function [sumTTA, countTTA] = accumulateConditionEpochs(sumTTA, countTTA, signalMatrix, onsetsByCondition, eventTRs_prior, eventTRs_after, rowIdx)
% Extracts the [-eventTRs_prior, +eventTRs_after) peristimulus window around
% each onset TR for each of the 13 conditions and accumulates sum/count
% (NaN-padded at run edges, exactly like ttave_compute.m's try/catch block)
% into sumTTA/countTTA at rowIdx. condTTA = sumTTA./countTTA afterward is an
% elementwise nanmean, robust to edge-padding NaNs.

    nConditions = numel(onsetsByCondition);
    nTimepoints = eventTRs_prior + eventTRs_after;
    nTRsRun = size(signalMatrix, 2);
    nRows = size(signalMatrix, 1);

    for ci = 1:nConditions
        onsets = onsetsByCondition{ci};

        for ti = 1:numel(onsets)
            startIdx = onsets(ti) - eventTRs_prior;
            endIdx = onsets(ti) + eventTRs_after - 1;

            epoch = nan(nRows, nTimepoints);
            validRange = max(startIdx,1):min(endIdx,nTRsRun);
            if ~isempty(validRange)
                epoch(:, (validRange(1)-startIdx+1):(validRange(1)-startIdx+numel(validRange))) = ...
                    signalMatrix(:, validRange);
            end

            mask = ~isnan(epoch);
            epoch(~mask) = 0;

            current = reshape(sumTTA(ci, rowIdx, :), nRows, nTimepoints) + epoch;
            sumTTA(ci, rowIdx, :) = reshape(current, 1, nRows, nTimepoints);

            currentC = reshape(countTTA(ci, rowIdx, :), nRows, nTimepoints) + mask;
            countTTA(ci, rowIdx, :) = reshape(currentC, 1, nRows, nTimepoints);
        end
    end
end

function T = buildLongTable(Tmeta, condTTA_dg, condTTA_da, dg_stimNames, da_stimNames, TRvals, nTimepoints)
% Expands the per-subject metadata table (one row per vertex) into a
% long-format table (one row per vertex per TR), attaching the 13 dg + 13 da
% condition values at each TR.

    nSel = height(Tmeta);
    blocks = cell(nTimepoints, 1);

    for t = 1:nTimepoints
        blockMeta = Tmeta;
        blockMeta.TR = repmat(TRvals(t), nSel, 1);

        for ci = 1:numel(dg_stimNames)
            blockMeta.(dg_stimNames{ci}) = reshape(condTTA_dg(ci,:,t), nSel, 1);
        end
        for ci = 1:numel(da_stimNames)
            blockMeta.(da_stimNames{ci}) = reshape(condTTA_da(ci,:,t), nSel, 1);
        end

        blocks{t} = blockMeta;
    end

    T = vertcat(blocks{:});
end

function plotExampleVoxels(allRawTables, allModelTables, subjects, exampleSubject, exampleROI, aggregationMode, ...
    offsetRaw, offsetModel, savePDFFlag, figSaveDir)
% For exampleSubject/exampleROI (included==1 vertices only), plots the
% "stationary" (orientation, no-motion) blank-subtracted TTA per polar
% angle bin -- raw data (dashed) vs. model fit (solid), one color per bin
% -- as two subplots side by side: dg (cartexp_) and da (polexp_).
%
% exampleSubject may be a specific subject ID, or 'all' -- in which case
% the exact same per-subject computation is run for every processed
% subject and the resulting per-bin curves are averaged across subjects
% before plotting.
%
% aggregationMode:
%   'peakR2'  - each bin's curve is the single vertex with the highest
%               combined GLM-R2 across BOTH experiments (mean of
%               cartexp_R2 and polexp_R2, not pRF_r2) -- the same vertex is
%               therefore used for both the dg and da subplots
%   'meanROI' - each bin's curve is the mean across all included vertices
%
% offsetRaw/offsetModel: if true, that curve type is shifted so its own
% pre-stimulus baseline (TR -5:-1) is 0-centered, applied independently per
% curve.
%
% savePDFFlag/figSaveDir: if savePDFFlag is true, saves the figure as a
% vector PDF into figSaveDir (created if needed).
%
% Reads directly from the allRawTables/allModelTables per-subject cell
% arrays (no need to reload the saved .mat files).

    if strcmpi(exampleSubject, 'all')
        subjectIndices = find(~cellfun(@isempty, allRawTables));
        subjectLabel = sprintf('all subjects (n=%d)', numel(subjectIndices));
    else
        si_example = find(strcmp(subjects, exampleSubject), 1);
        if isempty(si_example) || isempty(allRawTables{si_example})
            warning('%s not found/processed in this run -- skipping example plot.', exampleSubject);
            return
        end
        subjectIndices = si_example;
        subjectLabel = exampleSubject;
    end

    if isempty(subjectIndices)
        warning('No processed subjects found -- skipping example plot.');
        return
    end

    angleBins = 0:45:315;
    colors = cyclicPerceptualColors(angleBins, 71, 44);
    projectPrefixes = {'cartexp_', 'polexp_'};
    projectLabels = {'dg (Cartesian)', 'da (Polar)'};

    fig = figure('Name', sprintf('%s %s example voxels (%s)', subjectLabel, exampleROI, aggregationMode));
    fig.Position = [1 954 2383 383];

    for pj = 1:numel(projectPrefixes)
        projectPrefix = projectPrefixes{pj};
        blankCol = strcat(projectPrefix, 'blank');

        subplot(1, 2, pj)
        hold on
        legendHandles = gobjects(numel(angleBins), 1);
        legendLabels = strings(numel(angleBins), 1);

        for bi = 1:numel(angleBins)
            bin = angleBins(bi);
            col = colors(bi, :);

            rawCurves = {};
            modelCurves = {};
            rawTR = [];
            modelTR = [];

            % single subject: this loop runs once. exampleSubject=='all':
            % runs once per processed subject, and the resulting curves
            % are averaged across subjects below.
            for k = 1:numel(subjectIndices)
                Traw = allRawTables{subjectIndices(k)};
                Tmodel = allModelTables{subjectIndices(k)};

                Traw = Traw(Traw.visual_area == exampleROI & Traw.included == 1, :);
                Tmodel = Tmodel(Tmodel.visual_area == exampleROI & Tmodel.included == 1, :);

                if isempty(Traw)
                    continue
                end

                % "orientation conditions" = the stationary (no-motion)
                % condition columns for this project (excludes *_blank,
                % which has no "stationary" in its name).
                varNames = Traw.Properties.VariableNames;
                orientationCols = varNames(startsWith(varNames, projectPrefix) & contains(varNames, 'stationary'));

                if isempty(orientationCols) || ~ismember(blankCol, varNames) || ...
                        ~ismember('cartexp_R2', varNames) || ~ismember('polexp_R2', varNames)
                    continue
                end

                combinedRaw = mean(Traw{:, orientationCols}, 2) - Traw.(blankCol);
                combinedModel = mean(Tmodel{:, orientationCols}, 2) - Tmodel.(blankCol);

                % combine GLM R^2 across BOTH experiments (not just this
                % project's own) so 'peakR2' picks the same vertex per bin
                % regardless of which subplot (dg or da) is being drawn
                combinedR2Raw = mean([Traw.cartexp_R2, Traw.polexp_R2], 2);
                combinedR2Model = mean([Tmodel.cartexp_R2, Tmodel.polexp_R2], 2);

                [subjRawTR, subjRawVals, subjModelTR, subjModelVals] = computeBinCurve( ...
                    Traw, Tmodel, combinedRaw, combinedModel, bin, combinedR2Raw, combinedR2Model, aggregationMode);

                if isempty(subjRawVals)
                    continue
                end

                rawCurves{end+1} = subjRawVals(:)'; %#ok<AGROW>
                modelCurves{end+1} = subjModelVals(:)'; %#ok<AGROW>
                rawTR = subjRawTR;
                modelTR = subjModelTR;
            end

            if isempty(rawCurves)
                continue
            end

            rawVals = mean(cell2mat(rawCurves'), 1);
            modelVals = mean(cell2mat(modelCurves'), 1);

            rawVals = baselineCenter(rawTR, rawVals, offsetRaw);
            modelVals = baselineCenter(modelTR, modelVals, offsetModel);

            plot(rawTR, rawVals, 'o', 'LineWidth', 1.5, 'Color', col, ...
                'MarkerFaceColor', col, 'MarkerEdgeColor', 'w', 'HandleVisibility', 'off')
            h = plot(modelTR, modelVals, '-', 'LineWidth', 1.5, 'Color', col);
            legendHandles(bi) = h;
            legendLabels(bi) = sprintf('%d^\\circ', bin);
            

        end

        yline(0, 'k:', 'HandleVisibility', 'off')
        xline(0, 'k:', 'HandleVisibility', 'off')
        hold off

        xlabel('TR')
        ylabel('% signal change')
        title(sprintf('%s: %s', projectLabels{pj}, exampleROI))

        keep = isgraphics(legendHandles);
        legend(legendHandles(keep), legendLabels(keep), 'Location', 'best')

        if strcmpi(exampleSubject, 'all')
            if strcmpi(aggregationMode, 'peakR2')
                ylim([-0.5 3])
            elseif strcmpi(aggregationMode, 'meanROI')
                ylim([-.2 1])
            end
        end

    end

    sgtitle(sprintf('%s %s orientation response (stationary avg - blank), %s (dashed=raw, solid=model)', ...
        subjectLabel, exampleROI, aggregationMode))

    if savePDFFlag
        savePDF(fig, figSaveDir, sprintf('exampleVoxels_%s_%s_%s', subjectLabel, exampleROI, aggregationMode));
    end
end

function [rawTR, rawVals, modelTR, modelVals] = computeBinCurve(Traw, Tmodel, combinedRaw, combinedModel, bin, r2VecRaw, r2VecModel, aggregationMode)
% Returns the raw/model curve (across TR) for one polar angle bin of one
% subject's (already ROI/included-filtered) tables, using either the
% single highest-R^2 vertex ('peakR2') or the mean across all vertices in
% the bin ('meanROI'). Returns empty arrays if the bin has no vertices.
%
% r2VecRaw/r2VecModel are R^2 vectors aligned with Traw/Tmodel's rows (the
% caller decides what R^2 to use -- e.g. a single project's own R^2, or a
% combined score across experiments so the same vertex is picked for both).

    rawTR = []; rawVals = []; modelTR = []; modelVals = [];

    switch aggregationMode
        case 'peakR2'
            binRowsRaw = Traw.pRF_angle_bin == bin;
            if ~any(binRowsRaw)
                return
            end
            maxR2 = max(r2VecRaw(binRowsRaw));

            % identify the single highest-R^2 vertex's rows by exact match
            % (safe here since it's compared against a value taken from
            % this same vector)
            rowsRaw = binRowsRaw & r2VecRaw == maxR2;
            rowsModel = Tmodel.pRF_angle_bin == bin & r2VecModel == maxR2;

            [rawTR, sortIdx] = sort(Traw.TR(rowsRaw));
            rawTmp = combinedRaw(rowsRaw); rawVals = rawTmp(sortIdx);

            [modelTR, sortIdxM] = sort(Tmodel.TR(rowsModel));
            modelTmp = combinedModel(rowsModel); modelVals = modelTmp(sortIdxM);

        case 'meanROI'
            rowsRaw = Traw.pRF_angle_bin == bin;
            rowsModel = Tmodel.pRF_angle_bin == bin;
            if ~any(rowsRaw)
                return
            end

            [rawTR, ~, gRaw] = unique(Traw.TR(rowsRaw));
            rawVals = accumarray(gRaw, combinedRaw(rowsRaw), [], @mean);

            [modelTR, ~, gModel] = unique(Tmodel.TR(rowsModel));
            modelVals = accumarray(gModel, combinedModel(rowsModel), [], @mean);

        otherwise
            error('Unknown aggregationMode "%s" -- use "peakR2" or "meanROI".', aggregationMode);
    end
end

function plotOrientationConditionsAllVoxels(allRawTables, allModelTables, subjects, exampleSubject, exampleROI, ...
    offsetRaw, offsetModel, savePDFFlag, figSaveDir)
% For exampleSubject/exampleROI (included==1 vertices only), plots each of
% the 4 "stationary" (no-motion) orientation conditions as its own TTA
% curve, blank-subtracted per vertex per TR and then averaged across all
% valid vertices in the ROI -- raw data (markers) vs. model fit (solid
% line), one color per condition -- as two subplots side by side: dg
% (cartexp_) and da (polexp_). *_blank is excluded from the 4 plotted
% conditions (it's the baseline being subtracted, not an orientation
% condition itself).
%
% exampleSubject may be a specific subject ID, or 'all' -- in which case
% the exact same per-subject computation is run for every processed
% subject and the resulting per-condition curves are averaged across
% subjects before plotting.
%
% offsetRaw/offsetModel: if true, that curve type is shifted so its own
% pre-stimulus baseline (TR -5:-1) is 0-centered, applied independently per
% curve.
%
% savePDFFlag/figSaveDir: if savePDFFlag is true, saves the figure as a
% vector PDF into figSaveDir (created if needed).
%
% Reads directly from the allRawTables/allModelTables per-subject cell
% arrays (no need to reload the saved .mat files).

    if strcmpi(exampleSubject, 'all')
        subjectIndices = find(~cellfun(@isempty, allRawTables));
        subjectLabel = sprintf('all subjects (n=%d)', numel(subjectIndices));
    else
        si_example = find(strcmp(subjects, exampleSubject), 1);
        if isempty(si_example) || isempty(allRawTables{si_example})
            warning('%s not found/processed in this run -- skipping plot.', exampleSubject);
            return
        end
        subjectIndices = si_example;
        subjectLabel = exampleSubject;
    end

    if isempty(subjectIndices)
        warning('No processed subjects found -- skipping plot.');
        return
    end

    projectPrefixes = {'cartexp_', 'polexp_'};
    projectLabels = {'dg (Cartesian)', 'da (Polar)'};

    fig = figure('Name', sprintf('%s %s orientation conditions (all voxels)', subjectLabel, exampleROI));
    fig.Position = [1 954 2383 383];

    for pj = 1:numel(projectPrefixes)
        projectPrefix = projectPrefixes{pj};
        blankCol = strcat(projectPrefix, 'blank');

        % "motion_none" / orientation conditions = the *_stationary columns
        % for this project (excludes *_blank, which has no "stationary" in
        % its name). Taken from the first processed subject -- all subjects
        % share the same 13 condition columns per project.
        orientationCols = {};
        for k = 1:numel(subjectIndices)
            candidateVarNames = allRawTables{subjectIndices(k)}.Properties.VariableNames;
            orientationCols = candidateVarNames(startsWith(candidateVarNames, projectPrefix) & contains(candidateVarNames, 'stationary'));
            if ~isempty(orientationCols)
                break
            end
        end

        subplot(1, 2, pj)

        if isempty(orientationCols)
            warning('Could not find "%s*stationary" columns -- skipping %s panel.', projectPrefix, projectPrefix);
            continue
        end

        colors = lines(numel(orientationCols));

        hold on
        legendHandles = gobjects(numel(orientationCols), 1);
        legendLabels = strings(numel(orientationCols), 1);

        for ci = 1:numel(orientationCols)
            condCol = orientationCols{ci};
            col = colors(ci, :);

            rawCurves = {};
            modelCurves = {};
            rawTR = [];
            modelTR = [];

            % single subject: this loop runs once. exampleSubject=='all':
            % runs once per processed subject, and the resulting curves
            % are averaged across subjects below.
            for k = 1:numel(subjectIndices)
                Traw = allRawTables{subjectIndices(k)};
                Tmodel = allModelTables{subjectIndices(k)};

                Traw = Traw(Traw.visual_area == exampleROI & Traw.included == 1, :);
                Tmodel = Tmodel(Tmodel.visual_area == exampleROI & Tmodel.included == 1, :);

                if isempty(Traw) || ~ismember(condCol, Traw.Properties.VariableNames) || ...
                        ~ismember(blankCol, Traw.Properties.VariableNames)
                    continue
                end

                % blank-subtract per vertex per TR, then average within each
                % polar angle bin first and average those 8 bin-means
                % together (equal weight per bin, not per voxel)
                rawDiff = Traw.(condCol) - Traw.(blankCol);
                [subjRawTR, subjRawVals] = binWeightedMeanByTR(Traw.TR, Traw.pRF_angle_bin, rawDiff);

                modelDiff = Tmodel.(condCol) - Tmodel.(blankCol);
                [subjModelTR, subjModelVals] = binWeightedMeanByTR(Tmodel.TR, Tmodel.pRF_angle_bin, modelDiff);

                rawCurves{end+1} = subjRawVals(:)'; %#ok<AGROW>
                modelCurves{end+1} = subjModelVals(:)'; %#ok<AGROW>
                rawTR = subjRawTR;
                modelTR = subjModelTR;
            end

            if isempty(rawCurves)
                continue
            end

            rawVals = mean(cell2mat(rawCurves'), 1);
            modelVals = mean(cell2mat(modelCurves'), 1);

            rawVals = baselineCenter(rawTR, rawVals, offsetRaw);
            modelVals = baselineCenter(modelTR, modelVals, offsetModel);

            plot(rawTR, rawVals, 'o', 'LineWidth', 1.5, 'Color', col, ...
                'MarkerFaceColor', col, 'MarkerEdgeColor', 'w', 'HandleVisibility', 'off')
            h = plot(modelTR, modelVals, '-', 'LineWidth', 1.5, 'Color', col);
            legendHandles(ci) = h;
            legendLabels(ci) = strrep(erase(condCol, projectPrefix), '_', ' ');
        end

        yline(0, 'k:', 'HandleVisibility', 'off')
        xline(0, 'k:', 'HandleVisibility', 'off')
        hold off

        xlabel('TR')
        ylabel('% signal change')
        title(sprintf('%s: %s', projectLabels{pj}, exampleROI))

        keep = isgraphics(legendHandles);
        legend(legendHandles(keep), legendLabels(keep), 'Location', 'best', 'Interpreter', 'none')
    
        if strcmpi(exampleSubject, 'all')
            ylim([-.2 1])
        elseif strcmp(exampleSubject, 'sub-wlsubj124')
            ylim([-1 3])
        end
    
    end

    sgtitle(sprintf('%s %s: per-condition orientation response minus blank, all valid voxels (circles=raw, solid=model)', ...
        subjectLabel, exampleROI))

    if savePDFFlag
        savePDF(fig, figSaveDir, sprintf('orientationConditions_%s_%s', subjectLabel, exampleROI));
    end
end

function plotOrientationConditionsAllSubjects(allRawTables, allModelTables, subjects, exampleROI, ...
    offsetRaw, offsetModel, savePDFFlag, figSaveDir)
% For each subject, averages the 4 "stationary" (no-motion) orientation
% conditions (blank-subtracted, per vertex per TR) across all valid
% (included==1) vertices in the ROI, then plots each subject's resulting
% curve as its own line -- raw data (markers) vs. model fit (solid line),
% one color per subject -- as two subplots side by side: dg (cartexp_) and
% da (polexp_). Subjects not present/processed in allRawTables/
% allModelTables (empty cells) are skipped.
%
% offsetRaw/offsetModel: if true, that curve type is shifted so its own
% pre-stimulus baseline (TR -5:-1) is 0-centered, applied independently per
% curve.
%
% savePDFFlag/figSaveDir: if savePDFFlag is true, saves the figure as a
% vector PDF into figSaveDir (created if needed).
%
% Reads directly from the allRawTables/allModelTables per-subject cell
% arrays (no need to reload the saved .mat files).

    projectPrefixes = {'cartexp_', 'polexp_'};
    projectLabels = {'dg (Cartesian)', 'da (Polar)'};

    nSubjects = numel(subjects);
    bins = 0:45:315;
    colors = cyclicPerceptualColors(bins, 71, 44);
% 
%     colors = lines(nSubjects);

    fig = figure('Name', sprintf('%s orientation response, all subjects', exampleROI));
    fig.Position = [1 954 2383 383];

    for pj = 1:numel(projectPrefixes)
        projectPrefix = projectPrefixes{pj};

        subplot(1, 2, pj)
        hold on
        legendHandles = gobjects(nSubjects, 1);
        legendLabels = strings(nSubjects, 1);

        for si = 1:nSubjects
            if isempty(allRawTables{si})
                continue
            end

            Traw = allRawTables{si};
            Tmodel = allModelTables{si};

            Traw = Traw(Traw.visual_area == exampleROI & Traw.included == 1, :);
            Tmodel = Tmodel(Tmodel.visual_area == exampleROI & Tmodel.included == 1, :);

            if isempty(Traw)
                continue
            end

            varNames = Traw.Properties.VariableNames;
            orientationCols = varNames(startsWith(varNames, projectPrefix) & contains(varNames, 'stationary'));
            blankCol = strcat(projectPrefix, 'blank');

            if isempty(orientationCols) || ~ismember(blankCol, varNames)
                warning('Could not find "%s*stationary" / "%sblank" columns for %s -- skipping.', ...
                    projectPrefix, projectPrefix, subjects{si});
                continue
            end

            % blank-subtract per vertex per TR and average across the 4
            % orientation conditions, then average within each polar angle
            % bin first and average those 8 bin-means together (equal
            % weight per bin, not per voxel)
            combinedRaw = mean(Traw{:, orientationCols}, 2) - Traw.(blankCol);
            combinedModel = mean(Tmodel{:, orientationCols}, 2) - Tmodel.(blankCol);

            [rawTR, rawVals] = binWeightedMeanByTR(Traw.TR, Traw.pRF_angle_bin, combinedRaw);
            [modelTR, modelVals] = binWeightedMeanByTR(Tmodel.TR, Tmodel.pRF_angle_bin, combinedModel);

            rawVals = baselineCenter(rawTR, rawVals, offsetRaw);
            modelVals = baselineCenter(modelTR, modelVals, offsetModel);

            col = colors(si, :);
            plot(rawTR, rawVals, 'o', 'LineWidth', 1.5, 'Color', col, ...
                'MarkerFaceColor', col, 'MarkerEdgeColor', 'w', 'HandleVisibility', 'off')
            h = plot(modelTR, modelVals, '-', 'LineWidth', 1.5, 'Color', col);
            legendHandles(si) = h;
            legendLabels(si) = subjects{si};
        end

        yline(0, 'k:', 'HandleVisibility', 'off')
        xline(0, 'k:', 'HandleVisibility', 'off')
        hold off
        ylim([-0.7 2])

        xlabel('TR')
        ylabel('% signal change')
        title(sprintf('%s: %s', projectLabels{pj}, exampleROI))

        keep = isgraphics(legendHandles);
        legend(legendHandles(keep), legendLabels(keep), 'Location', 'best', 'Interpreter', 'none')

    end

    sgtitle(sprintf('%s orientation response (4-condition avg minus blank) per subject (circles=raw, solid=model)', ...
        exampleROI))

    if savePDFFlag
        savePDF(fig, figSaveDir, sprintf('orientationConditions_allSubjects_%s', exampleROI));
    end

    % model R^2 per subject, for the same included exampleROI voxels used
    % above, averaged separately for the dg (cartexp) and da (polexp)
    % experiments (cartexp_R2/polexp_R2 are identical in allRawTables and
    % allModelTables, so either table works here)
    fprintf('\nModel R^2 per subject (%s, included voxels):\n', exampleROI);
    for si = 1:nSubjects
        if isempty(allRawTables{si})
            continue
        end

        Tsi = allRawTables{si};
        Tsi = Tsi(Tsi.visual_area == exampleROI & Tsi.included == 1, :);

        if isempty(Tsi) || ~ismember('cartexp_R2', Tsi.Properties.VariableNames) || ...
                ~ismember('polexp_R2', Tsi.Properties.VariableNames)
            continue
        end

        fprintf('  %s: dg (cartexp) R^2 = %.3f, da (polexp) R^2 = %.3f\n', ...
            subjects{si}, mean(Tsi.cartexp_R2), mean(Tsi.polexp_R2));
    end
end

function plotOrientationConditionsOneBin(allRawTables, allModelTables, subjects, exampleSubject, exampleROI, targetBin, ...
    cartConditions, polConditions, offsetRaw, offsetModel, savePDFFlag, figSaveDir)
% Restricted variant of the orientation-condition plots above: instead of
% averaging/plotting across all 8 polar angle bins, uses only the single
% specified targetBin; instead of averaging all 4 orientation conditions
% per experiment, plots the specified conditions separately (each minus
% that experiment's own blank). Subplot structure matches
% plotOrientationConditionsAllSubjects: one subplot for dg/cartexp_ (left),
% one for da/polexp_ (right). Within each subplot, one line per condition
% (not per subject) -- color = condition, matching
% plotOrientationConditionsAllVoxels's pattern.
%
% cartConditions/polConditions: cell arrays of full column names (any
% length, e.g. 2 or all 4 stationary conditions) for the cartexp_/polexp_
% subplots respectively, e.g. cartConditions = {'cartexp_vertical_stationary',
% 'cartexp_horizontal_stationary'} (dg and da conditions are named
% completely differently, so these two lists are independent).
%
% exampleSubject may be a specific subject ID, or 'all' -- in which case
% the exact same per-subject computation is run for every processed
% subject and the resulting per-condition curves are averaged across
% subjects before plotting (same pattern as plotOrientationConditionsAllVoxels).
%
% offsetRaw/offsetModel: if true, that curve type is shifted so its own
% pre-stimulus baseline (TR -5:-1) is 0-centered, applied independently per
% curve.
%
% savePDFFlag/figSaveDir: if savePDFFlag is true, saves the figure as a
% vector PDF into figSaveDir (created if needed).
%
% After plotting, prints the R^2 of each condition's plotted model curve
% against its plotted raw curve (i.e. at this plot's own level of
% aggregation -- one bin, one condition, averaged across whichever
% subjects were selected), using GLMsingle's own R^2 definition (R^2 =
% 100*(1-sum((data-model)^2)/sum(data^2))). This is a different quantity
% from the whole-GLM cartexp_R2/polexp_R2 table columns -- those reflect
% the full 13-condition model fit at the single-vertex level, not this
% specific observer-averaged, single-bin, single-condition summary curve.
%
% Reads directly from the allRawTables/allModelTables per-subject cell
% arrays (no need to reload the saved .mat files).

    if strcmpi(exampleSubject, 'all')
        subjectIndices = find(~cellfun(@isempty, allRawTables));
        subjectLabel = sprintf('all subjects (n=%d)', numel(subjectIndices));
    else
        si_example = find(strcmp(subjects, exampleSubject), 1);
        if isempty(si_example) || isempty(allRawTables{si_example})
            warning('%s not found/processed in this run -- skipping plot.', exampleSubject);
            return
        end
        subjectIndices = si_example;
        subjectLabel = exampleSubject;
    end

    if isempty(subjectIndices)
        warning('No processed subjects found -- skipping plot.');
        return
    end

    projectConditions = {cartConditions, polConditions};
    projectBlanks = {'cartexp_blank', 'polexp_blank'};
    projectLabels = {'dg (Cartesian)', 'da (Polar)'};

    fig = figure('Name', sprintf('%s %s bin%d selected conditions', subjectLabel, exampleROI, targetBin));
    fig.Position = [1 954 2383 383];

    axHandles = gobjects(numel(projectConditions), 1);
    r2CondNames = {};
    r2Values = [];

    for pj = 1:numel(projectConditions)
        theseConditions = projectConditions{pj};
        blankCol = projectBlanks{pj};
        colors = lines(numel(theseConditions));

        axHandles(pj) = subplot(1, 2, pj);
        hold on
        legendHandles = gobjects(numel(theseConditions), 1);
        legendLabels = strings(numel(theseConditions), 1);

        for ci = 1:numel(theseConditions)
            condCol = theseConditions{ci};
            col = colors(ci, :);

            rawCurves = {};
            modelCurves = {};
            rawTR = [];
            modelTR = [];

            % single subject: this loop runs once. exampleSubject=='all':
            % runs once per processed subject, and the resulting curves
            % are averaged across subjects below.
            for k = 1:numel(subjectIndices)
                Traw = allRawTables{subjectIndices(k)};
                Tmodel = allModelTables{subjectIndices(k)};

                Traw = Traw(Traw.visual_area == exampleROI & Traw.included == 1 & Traw.pRF_angle_bin == targetBin, :);
                Tmodel = Tmodel(Tmodel.visual_area == exampleROI & Tmodel.included == 1 & Tmodel.pRF_angle_bin == targetBin, :);

                if isempty(Traw) || ~ismember(condCol, Traw.Properties.VariableNames) || ...
                        ~ismember(blankCol, Traw.Properties.VariableNames)
                    continue
                end

                rawDiff = Traw.(condCol) - Traw.(blankCol);
                [subjRawTR, ~, gRaw] = unique(Traw.TR);
                subjRawVals = accumarray(gRaw, rawDiff, [], @mean);

                modelDiff = Tmodel.(condCol) - Tmodel.(blankCol);
                [subjModelTR, ~, gModel] = unique(Tmodel.TR);
                subjModelVals = accumarray(gModel, modelDiff, [], @mean);

                rawCurves{end+1} = subjRawVals(:)'; %#ok<AGROW>
                modelCurves{end+1} = subjModelVals(:)'; %#ok<AGROW>
                rawTR = subjRawTR;
                modelTR = subjModelTR;
            end

            if isempty(rawCurves)
                continue
            end

            rawVals = mean(cell2mat(rawCurves'), 1);
            modelVals = mean(cell2mat(modelCurves'), 1);

            rawVals = baselineCenter(rawTR, rawVals, offsetRaw);
            modelVals = baselineCenter(modelTR, modelVals, offsetModel);

            % R^2 of this condition's observer-averaged model curve
            % against its observer-averaged raw curve, using GLMsingle's
            % own R^2 definition (R^2 = 100*(1-sum((data-model)^2)/sum(data^2)))
            % applied at this plot's level of aggregation (one bin, one
            % condition, averaged across the selected subjects) rather
            % than the whole-GLM cartexp_R2/polexp_R2 columns.
            r2CondNames{end+1} = condCol; %#ok<AGROW>
            r2Values(end+1) = 100 * (1 - sum((rawVals - modelVals).^2) / sum(rawVals.^2)); %#ok<AGROW>

            plot(rawTR, rawVals, 'o', 'LineWidth', 1.5, 'Color', col, ...
                'MarkerFaceColor', col, 'MarkerEdgeColor', 'w', 'HandleVisibility', 'off')
            h = plot(modelTR, modelVals, '-', 'LineWidth', 1.5, 'Color', col);
            legendHandles(ci) = h;
            legendLabels(ci) = strrep(condCol, '_', ' ');
        end

        yline(0, 'k:', 'HandleVisibility', 'off')
        xline(0, 'k:', 'HandleVisibility', 'off')
        hold off

        xlabel('TR')
        ylabel('% signal change')
        title(sprintf('%s: %s (%d^\\circ)', projectLabels{pj}, exampleROI, targetBin))

        keep = isgraphics(legendHandles);
        legend(legendHandles(keep), legendLabels(keep), 'Location', 'best', 'Interpreter', 'none')
    end

    % share the same y-axis range across both subplots, and give each a
    % 2:3 (y:x) plot box aspect ratio
    sharedYLim = [min(arrayfun(@(a) min(ylim(a)), axHandles)), max(arrayfun(@(a) max(ylim(a)), axHandles))];
    for pj = 1:numel(axHandles)
        ylim(axHandles(pj), sharedYLim);
        pbaspect(axHandles(pj), [3 2 1]);
    end

    sgtitle(sprintf('%s %s %d^\\circ: selected conditions minus blank (circles=raw, solid=model)', ...
        subjectLabel, exampleROI, targetBin))

    if savePDFFlag
        savePDF(fig, figSaveDir, sprintf('orientationConditionsOneBin_bin%d_%s_%s', targetBin, exampleROI, subjectLabel));
    end

    fprintf('\nR^2 (observer-averaged model vs. raw TTA, %s, %s, %d^\\circ):\n', subjectLabel, exampleROI, targetBin);
    for ri = 1:numel(r2CondNames)
        fprintf('  %s: R^2 = %.2f\n', r2CondNames{ri}, r2Values(ri));
    end
end

function [TRout, valsOut] = binWeightedMeanByTR(TRvec, angleBinVec, valueVec)
% Averages valueVec within each (TR, polar-angle-bin) group, then averages
% those per-bin means across the 8 polar angle bins (0:45:315) to get one
% value per TR. This weights each of the 8 locations equally, rather than
% pooling all vertices together (which would implicitly weight bins with
% more vertices more heavily). Bins with no data at a given TR are simply
% omitted from that TR's across-bin average rather than counted as zero.

    angleBins = 0:45:315;

    [TRout, ~, gTR] = unique(TRvec);
    nTR = numel(TRout);
    valsOut = nan(nTR, 1);

    for t = 1:nTR
        rowsT = gTR == t;
        binMeans = nan(numel(angleBins), 1);

        for b = 1:numel(angleBins)
            rowsB = rowsT & angleBinVec == angleBins(b);
            if any(rowsB)
                binMeans(b) = mean(valueVec(rowsB));
            end
        end

        valsOut(t) = mean(binMeans, 'omitnan');
    end
end

function vals = baselineCenter(TRvec, vals, doOffset)
% If doOffset, subtracts the mean of the pre-stimulus baseline (TR -5:-1)
% from the whole curve, so that baseline period is 0-centered. No-op if
% doOffset is false or if the curve has no timepoints in that range.

    if ~doOffset
        return
    end

    baselineIdx = TRvec >= -5 & TRvec <= -1;
    if any(baselineIdx)
        vals = vals - mean(vals(baselineIdx));
    end
end

function colors = cyclicPerceptualColors(hueAnglesDeg, L, C)
% Generates a circular/gradient color sequence by sweeping hue at fixed
% lightness L and chroma C in CIE L*a*b* (D65), rather than raw HSV hue.
% Raw HSV hue is not perceptually uniform -- e.g. its green region is much
% "wider" than red/blue, so evenly-spaced HSV hues can land two samples
% almost on top of each other there (this is why hsv(8)'s two green-ish
% entries look nearly identical). Sweeping hue in Lab instead keeps the
% same even-angle circular structure (0 degrees is still adjacent to 315,
% matching the polar angle bins 1:1) while spacing the 8 colors more
% evenly by how different they actually look.
%
% Defaults tuned (against the same six-check validator used for the
% project's categorical palettes) to clear the lightness-band and
% chroma-floor checks and to maximize worst-adjacent-pair separation for a
% true evenly-spaced 8-point circular sweep: L=71, C=44 -> worst-pair
% normal-vision CVD-simulated separation ~7.8 (vs ~4.8 for hsv(8)'s worst
% pair). This does not clear the stricter >=15 bar used for the arbitrary-
% identity 8-slot categorical palette elsewhere in this method -- that
% palette is free to place its 8 hues non-uniformly to maximize
% separation, which would break the circular/gradient property requested
% here.

    n = numel(hueAnglesDeg);
    colors = zeros(n, 3);
    for i = 1:n
        h = hueAnglesDeg(i);
        a = C * cosd(h);
        b = C * sind(h);
        colors(i, :) = labToSRGB(L, a, b);
    end
end

function rgb = labToSRGB(L, a, b)
% CIE L*a*b* (D65 white point) to sRGB (gamma-corrected, [0,1], clipped).

    fy = (L + 16) / 116;
    fx = fy + a / 500;
    fz = fy - b / 200;

    delta = 6/29;
    finv = @(t) (t > delta) .* (t.^3) + (t <= delta) .* (3*delta^2*(t - 4/29));

    Xn = 95.0489; Yn = 100.0; Zn = 108.8840;

    X = Xn * finv(fx) / 100;
    Y = Yn * finv(fy) / 100;
    Z = Zn * finv(fz) / 100;

    Rl =  3.2406*X - 1.5372*Y - 0.4986*Z;
    Gl = -0.9689*X + 1.8758*Y + 0.0415*Z;
    Bl =  0.0557*X - 0.2040*Y + 1.0570*Z;

    gammaCorrect = @(c) (c <= 0.0031308) .* (12.92*c) + (c > 0.0031308) .* (1.055*max(c,0).^(1/2.4) - 0.055);
    rgb = gammaCorrect([Rl, Gl, Bl]);
    rgb = min(max(rgb, 0), 1);
end

function savePDF(fig, saveDir, filenameBase)
% Saves fig as a vector PDF into saveDir (created if needed), sanitizing
% filenameBase (which may contain spaces/parens, e.g. "all subjects (n=8)")
% into a filesystem-safe file name.

    if ~isfolder(saveDir)
        mkdir(saveDir)
    end

    safeName = regexprep(filenameBase, '[^\w-]+', '_');
    outPath = fullfile(saveDir, [safeName '.pdf']);
    exportgraphics(fig, outPath, 'ContentType', 'vector');
    fprintf('Saved figure: %s\n', outPath);
end
