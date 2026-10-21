% computeObserverPrecisionWeights.m
%
% Computes one within-observer measurement-reliability estimate
% (sigma_i^2) per observer, per ROI, per experiment, following the same
% across-RUN resampling logic implemented by J. Winawer
% (github.com/JWinawer/DriftingGrating, branch 'winawer',
% Reproduction/cleanroom/diagnose_within_observer_error.m /
% collect_runwise_betas.m), generalized from V1-only to all 8 ROIs and
% collapsed from 4 asymmetry-specific reliabilities down to a single,
% asymmetry-agnostic number per (observer, ROI, experiment).
%
% METHOD (matches Jon's approach):
%   - Reliability is estimated by resampling RUNS (bootstrap, 500 draws),
%     not vertices -- vertex-resampling only characterizes which patch of
%     cortex happened to be sampled, not the reliability of the
%     measurement itself (this was tried and explicitly withdrawn in
%     Jon's repository for that reason).
%   - Restricted to vertices meeting eccentricity 4-8 deg and pRF R^2 >=
%     0.1. Note: Jon's version used strictly >0.1; this uses >=0.1 to
%     match meanWithinLabel.m, the pipeline that actually generates
%     meanBOLDpa for the main analysis.
%   - Computed on the orientation-minus-blank contrast -- the blank
%     condition plus the 4 orientation conditions (5 of the 13 total GLM
%     conditions) -- the same set of conditions the beta-asymmetry
%     analysis itself uses throughout this project
%     (comparisonName='orientation_minus_baseline'). The blank condition
%     is included deliberately: every orientation-vs-blank contrast
%     inherits the blank condition's own measurement noise through
%     subtraction, so omitting it would understate the true noise in the
%     quantity that is actually analyzed downstream.
%   - Each observer's betas are gain-rescaled (subjectScale = groupGain /
%     gain_i, geometric mean) before anything else, matching every other
%     script in this project that touches BOLD magnitudes.
%
% COLLAPSING TO ONE NUMBER PER (OBSERVER, ROI, EXPERIMENT):
%   Jon's own script computes a separate reliability for each of the 4
%   asymmetries. Here, at each bootstrap draw, the full
%   orientation-minus-blank, per-wedge matrix is computed instead (4
%   conditions x 8 polar-angle wedges = 32 cells) -- the same
%   intermediate quantity his 4 asymmetries are themselves built from --
%   and the reported sigma_i^2 is the MEAN, ACROSS ALL 32 CELLS, of each
%   cell's own bootstrap variance. This is asymmetry-agnostic (reusable
%   for any later analysis built from these same cells, in either
%   experiment) and reflects the typical noise level of the individual
%   (condition, location) values that actually enter a downstream
%   regression as rows -- the quantity a per-subject regression weight
%   should track.
%
% WHY THIS SAVES sigma_i^2, NOT A FINAL NORMALIZED WEIGHT w_i:
%   w_i = 1/(tau^2 + sigma_i^2) also needs tau^2, the between-observer
%   variance of whatever specific effect is being combined -- and tau^2
%   depends on the actual outcome values (y_i) of that specific analysis,
%   which differ by asymmetry/experiment/ROI in a way that isn't fixed in
%   advance. sigma_i^2 is the reusable, analysis-independent ingredient;
%   w_i should be computed at the point of use once tau^2 for that
%   specific analysis is known.
%
% OUTPUTS (saved to derivatives/summaryTables/):
%   observerPrecisionWeights.mat  - sigma2 (nSubj x nROI x nExp), subjects,
%                                   rois, experiments
%   observerPrecisionWeights.csv  - long format: subject, roi, experiment,
%                                   sigma2, nRun, nVertsUsed
%
% Set overwriteExtraction=true to force re-extraction of the per-run
% betas (slow: streams the full single-trial modelmd file per
% subject/experiment, ~2 min each per Jon's own cost note, ~13*2=26
% subject-experiment pairs here since not every subject has both dg and
% da). Extraction is skipped automatically if its output file already
% exists.

clear; clc

githubDir = '~/Documents/GitHub';
addpath(genpath(fullfile(githubDir, 'DriftingGrating', 'AnalysisCode')));
cd(fullfile(githubDir, 'DriftingGrating', 'AnalysisCode'));
bidsDir = '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
setup_user('rania', bidsDir);

projectSettings = loadConfig(githubDir);
rois = projectSettings.rois; % {'V1','V2','V3','V3a','V3b','hV4','pMT','pMST'}
nROIs = numel(rois);

subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
    'sub-0395', 'sub-0426', 'sub-0250', 'sub-0442', 'sub-wlsubj121', 'sub-wlsubj127', ...
    'sub-0397', 'sub-0427'};
nSubj = numel(subjects);

experiments = {'dg','da'};
nExp = numel(experiments);

minECC = 4; maxECC = 8; minVAREXP = 0.1; % >=0.1, unlike Jon's strict >0.1
polarAngleBinWidth = 45;
nBoot = 1000;
orientationCols = 9:12; % S_0, S_90, S_45, S_135 (see extraction note below)
blankCol = 13;

gainSummaryFile = fullfile(bidsDir, 'derivatives', 'summaryTables', 'gainSummary.mat');
gainWeights = retrieveObserverGainWeights(subjects, gainSummaryFile);
groupGain = exp(mean(log(gainWeights)));
subjectScale = groupGain ./ gainWeights;

extractDir = fullfile(bidsDir, 'derivatives', 'summaryTables', 'runwiseBetas');
if ~isfolder(extractDir), mkdir(extractDir); end
saveDir = fullfile(bidsDir, 'derivatives', 'summaryTables');

overwriteExtraction = false;

%% Step 1: extract per-run, per-vertex, per-condition betas for the union
% of all 8 ROIs' vertices (once per subject x experiment)

for si = 1:nSubj
    subjectname = subjects{si};
    for ei = 1:nExp
        projectName = experiments{ei};
        outFile = fullfile(extractDir, sprintf('runbetas_%s_%s.mat', subjectname, projectName));
        if isfile(outFile) && ~overwriteExtraction
            fprintf('%s %s: extraction exists, skipping\n', subjectname, projectName);
            continue
        end

        d = dir(fullfile(bidsDir, 'derivatives', strcat(projectName, 'GLM'), 'hRF_glmsingle', subjectname, 'ses-*'));
        d = d([d.isdir]);
        if isempty(d)
            fprintf('%s %s: no ses-* folder -- skipped (subject likely not run in this experiment)\n', subjectname, projectName);
            continue
        end
        srcFolder = fullfile(d(1).folder, d(1).name);

        designFile = fullfile(srcFolder, 'DESIGNINFO.mat');
        modelFile = fullfile(srcFolder, 'TYPED_FITHRF_GLMDENOISE_RR.mat');
        if ~isfile(designFile) || ~isfile(modelFile)
            warning('%s %s: missing DESIGNINFO.mat or TYPED_FITHRF_GLMDENOISE_RR.mat in %s -- skipped', ...
                subjectname, projectName, srcFolder);
            continue
        end

        t0 = tic;
        DI = load(designFile, 'stimorder', 'numtrialrun');
        runOfTrial = repelem(1:numel(DI.numtrialrun), DI.numtrialrun);
        cond = DI.stimorder(:).';
        nRun = numel(DI.numtrialrun);
        nCond = max(cond);

        M = load(modelFile, 'modelmd');
        B = squeeze(M.modelmd); % nVertex x nTrials
        clear M

        hSize = get_surfsize(subjectname);

        % union of vertex indices across all 8 ROIs, plus each ROI's own
        % indices (as row positions into that union) for later filtering
        roiVertexIdx = cell(nROIs,1);
        allIdx = [];
        for ri = 1:nROIs
            roiname = rois{ri};
            if strcmp(roiname,'V2') || strcmp(roiname,'V3')
                lh1 = read_label(subjectname, sprintf('retinotopy_RE/lh.%sv_REmanual', roiname));
                lh2 = read_label(subjectname, sprintf('retinotopy_RE/lh.%sd_REmanual', roiname));
                rh1 = read_label(subjectname, sprintf('retinotopy_RE/rh.%sv_REmanual', roiname));
                rh2 = read_label(subjectname, sprintf('retinotopy_RE/rh.%sd_REmanual', roiname));
                idx = [lh1(:,1)+1; lh2(:,1)+1; rh1(:,1)+hSize(1)+1; rh2(:,1)+hSize(1)+1];
            else
                lh = read_label(subjectname, sprintf('retinotopy_RE/lh.%s_REmanual', roiname));
                rh = read_label(subjectname, sprintf('retinotopy_RE/rh.%s_REmanual', roiname));
                idx = [lh(:,1)+1; rh(:,1)+hSize(1)+1];
            end
            roiVertexIdx{ri} = sort(unique(idx));
            allIdx = [allIdx; roiVertexIdx{ri}]; %#ok<AGROW>
        end
        allIdx = sort(unique(allIdx));

        if size(B,1) < max(allIdx)
            error('computeObserverPrecisionWeights:size', '%s %s: modelmd has %d vertices, ROI index max is %d.', ...
                subjectname, projectName, size(B,1), max(allIdx));
        end
        if size(B,2) ~= numel(cond)
            error('computeObserverPrecisionWeights:trials', '%s %s: modelmd has %d trials, DESIGNINFO has %d.', ...
                subjectname, projectName, size(B,2), numel(cond));
        end

        Ball = B(allIdx, :);
        clear B

        runBeta = nan(numel(allIdx), nCond, nRun, 'single');
        for c = 1:nCond
            for r = 1:nRun
                sel = (cond==c) & (runOfTrial==r);
                if any(sel), runBeta(:,c,r) = mean(Ball(:,sel), 2); end
            end
        end

        roiRows = cell(nROIs,1);
        for ri = 1:nROIs
            [~, roiRows{ri}] = ismember(roiVertexIdx{ri}, allIdx);
        end

        S = struct('subject', subjectname, 'project', projectName, 'sourceFolder', srcFolder, ...
            'allVertexIdx', allIdx, 'roiNames', {rois}, 'roiRows', {roiRows}, 'nLH', hSize(1), ...
            'runBeta', runBeta, 'nRun', nRun, 'nCond', nCond, ...
            'note', ['runBeta is nVert x nCond x nRun, mean single-trial beta per (vertex,condition,run). ' ...
                     'Conditions 1-8 motion, 9=S_0 10=S_90 11=S_45 12=S_135, 13=blank.'], ...
            'collected', datetime('now'));
        save(outFile, '-struct', 'S', '-v7.3');
        fprintf('%s %s: %d ROI-union vertices, %d cond x %d runs, %.0f s -> %s\n', ...
            subjectname, projectName, numel(allIdx), nCond, nRun, toc(t0), outFile);
    end
end

%% Step 2: pRF-based inclusion (ecc 4-8, R^2>=0.1) + wedge assignment,
% then bootstrap-over-runs reliability per (subject, ROI, experiment)

sigma2 = nan(nSubj, nROIs, nExp);
nRunUsed = nan(nSubj, nExp);
nVertsUsed = nan(nSubj, nROIs, nExp);
longRows = {};

for si = 1:nSubj
    subjectname = subjects{si};

    retDir = dir(fullfile(bidsDir, 'derivatives', 'prfvista_mov', subjectname, '**/stimfiles.mat'));
    if isempty(retDir)
        warning('%s: no prfvista_mov data -- skipped entirely', subjectname);
        continue
    end
    retDir = retDir(1).folder;
    hemis = {'lh','rh'};
    for hi = 1:numel(hemis)
        hemi = hemis{hi};
        ret.(sprintf('%s_pa', hemi)) = MRIread(fullfile(retDir, sprintf('%s.angle_adj.mgz', hemi)));
        ret.(sprintf('%s_ecc', hemi)) = MRIread(fullfile(retDir, sprintf('%s.eccen.mgz', hemi)));
        ret.(sprintf('%s_vexp', hemi)) = MRIread(fullfile(retDir, sprintf('%s.vexpl.mgz', hemi)));
    end
    angConv = map_theta([ret.lh_pa.vol, ret.rh_pa.vol]); % conventional degrees, 1 x nVertTotal
    eccAll = [ret.lh_ecc.vol, ret.rh_ecc.vol];
    r2All = [ret.lh_vexp.vol, ret.rh_vexp.vol];

    binCenters = 0:45:315;
    circDist = abs(mod(angConv(:) - binCenters + 180, 360) - 180);
    [minDist, binIdx] = min(circDist, [], 2);
    wedgeAll = nan(numel(angConv),1);
    validBin = minDist <= 22.5;
    wedgeAll(validBin) = binCenters(binIdx(validBin));

    includedAll = validBin & (eccAll(:) >= minECC) & (eccAll(:) <= maxECC) & (r2All(:) >= minVAREXP);

    for ei = 1:nExp
        projectName = experiments{ei};
        rbFile = fullfile(extractDir, sprintf('runbetas_%s_%s.mat', subjectname, projectName));
        if ~isfile(rbFile)
            continue
        end
        RB = load(rbFile);
        nRunUsed(si, ei) = RB.nRun;

        for ri = 1:nROIs
            rows = RB.roiRows{ri};
            rows = rows(rows > 0);
            if isempty(rows)
                continue
            end
            vIdx = RB.allVertexIdx(rows);
            incl = includedAll(vIdx);
            rows = rows(incl);
            wedge = wedgeAll(vIdx(incl));
            nVertsUsed(si, ri, ei) = numel(rows);
            if numel(rows) < 4 % need at least a few vertices for a meaningful bootstrap
                continue
            end

            rng(si, 'twister'); % reproducible per subject; explicit generator avoids
                                 % "current generator is the legacy generator" errors
                                 % from upstream toolbox state (e.g. rand('state',0))
            bootCells = nan(nBoot, 4, 8); % draw x condition x wedge
            for b = 1:nBoot
                runsSample = randi(RB.nRun, [1, RB.nRun]);
                Bmean = mean(RB.runBeta(rows, :, runsSample), 3, 'omitnan'); % nRows x nCond
                Bmean = double(Bmean) * subjectScale(si);
                C = Bmean(:, orientationCols) - Bmean(:, blankCol); % nRows x 4
                for w = 1:8
                    m = wedge == binCenters(w);
                    if any(m)
                        bootCells(b, :, w) = mean(C(m, :), 1);
                    end
                end
            end

            cellVar = squeeze(var(bootCells, 0, 1, 'omitnan')); % 4 x 8
            sigma2(si, ri, ei) = mean(cellVar(:), 'omitnan');

            longRows(end+1, :) = {subjectname, rois{ri}, projectName, sigma2(si,ri,ei), ...
                RB.nRun, nVertsUsed(si,ri,ei)}; %#ok<SAGROW>
        end
        fprintf('%s %s: reliability computed for %d/%d ROIs\n', subjectname, projectName, ...
            sum(~isnan(sigma2(si,:,ei))), nROIs);
    end
end

%% Save

save(fullfile(saveDir, 'observerPrecisionWeights.mat'), 'sigma2', 'subjects', 'rois', 'experiments', ...
    'nRunUsed', 'nVertsUsed', 'minECC', 'maxECC', 'minVAREXP', 'nBoot');

Tout = cell2table(longRows, 'VariableNames', {'subject','roi','experiment','sigma2','nRun','nVertsUsed'});
writetable(Tout, fullfile(saveDir, 'observerPrecisionWeights.csv'));

fprintf('\nSaved observerPrecisionWeights.mat and .csv to %s\n', saveDir);
