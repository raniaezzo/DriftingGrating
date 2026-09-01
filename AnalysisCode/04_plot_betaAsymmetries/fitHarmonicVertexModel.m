function fitHarmonicVertexModel(projectName, varargin)
% FITHARMONICVERTEXMODEL  Per-vertex harmonic-model alternative to
% fitAsymmetryRegression.m's wedge-based fit, following the SETTLED
% specification at github.com/JWinawer/DriftingGrating
% (Reproduction/SPECIFICATION.md), with one deliberate substitution: a
% 1000-draw percentile bootstrap CI in place of that spec's primary
% t-on-(n-1)-df interval (explicit choice, 2026-08-31 -- for consistency
% with every other CI in this pipeline).
%
%   fitHarmonicVertexModel('dg')
%   fitHarmonicVertexModel('da')
%   fitHarmonicVertexModel('dg', 'dgSubjectMode', 'matched')
%
% THE MODEL (see SUPPLEMENT_harmonic_model.md / SPECIFICATION.md sec 2 for
% the full derivation -- this is a standalone, from-scratch reproduction,
% not a port of that repo's code): for each vertex v, take its 4 raw
% stationary-orientation betas (s0/s90/s45/s135_v_b), subtract that
% vertex's own mean across the 4 (removes the blank baseline; no
% intercept needed since all 4 predictors below sum to zero across the 4
% conditions), and fit
%
%   y_vk = b1*cos(2*theta) + b2*cos(4*theta) + b3*cos(2*(theta-theta_V)) + b4*cos(4*(theta-theta_V))
%
% where theta is stimulus k's LOCAL Cartesian bar orientation at vertex v
% (theta = raw code {0,90,45,135} directly for dg, since dg's raw codes
% are already Cartesian everywhere -- REF_ORIENTATION.md sec 1; theta =
% mod(rawCode + theta_V - 90, 180) for da, the same UVM-rotation formula
% already verified elsewhere in this pipeline -- REF_ORIENTATION.md sec 2)
% and theta_V is that vertex's own CONTINUOUS pRF polar angle (not its
% wedge's nominal center -- that approximation is exactly what this model
% exists to avoid). b1=horizontal-vertical, b2=cardinal-oblique,
% b3=radial-tangential, b4=polar-cardinal-oblique -- REGARDLESS of
% project (unlike fitAsymmetryRegression.m's mainCardinal/mainSubset/etc,
% whose CONCEPT meaning swaps between dg/da, the harmonic b1..b4 are the
% same 4 physical concepts in both experiments). Coefficients are reported
% as 2*b (pro-minus-con scale), matching every other cache in this
% pipeline.
%
% DATA SOURCE: meanWithinLabel.m's allvoxelsBOLDpa.mat/allparamsBOLDpa.mat
% (contrast x wedge x roi x up-to-2000-vertex-slots x subject, params =
% [angle, eccen, R^2, wedge-distance]) -- already restricted to this
% project's standard ecc[4,8]/R^2>=0.1 inclusion criteria, and already
% storing each vertex's own CONTINUOUS pRF angle (verified directly
% against the .mat file, not assumed). No new vertex extraction needed.
% Since the per-wedge vertex-inclusion mask (built from pRF angle/ecc/R^2
% alone) does not depend on which of the 29 contrasts is being read, slot
% k of the vertex dimension refers to the SAME physical vertex for every
% contrast at a given (wedge, roi, subject) -- verified against
% meanWithinLabel.m's loop structure -- so the 4 raw betas needed for
% demeaning can be read at matching slot indices across contrasts.
%
% NOTE ON "BENSON ORDER": fitAsymmetryRegression.m must map wedge index ->
% nominal angle explicitly (anglevals=[90 45 0 315 270 225 180 135], NOT
% polarAngles=[0 45 90 ...] in sequence) because its design codes a
% vertex's condition from which WEDGE it's in. This script never does
% that: theta_V below is each vertex's own genuine measured angle, read
% directly from allparamsBOLDpa, never derived from its wedge index -- so
% the wedge-order convention is irrelevant here. The 8 wedges are only
% used as a vertex-partitioning device (they tile the full 360 degrees
% exactly, 8x45, with no gaps or overlap, so looping pa=1:8 and pooling
% every valid vertex covers the ROI's full angular range exactly once).
%
% VERTEX WEIGHTING: equal coverage over the 8 (45-degree) wedges --
% weight = 1 / (number of valid vertices in that subject's own wedge),
% i.e. the wedges each get equal total weight in the per-observer fit,
% matching the ROI route's own "equally weighs each polar-angle location"
% design (SPECIFICATION.md sec 5's settled choice; 15-degree bins were the
% supplement's superseded parameter, not used here).
%
% GAIN: applied per (observer, ROI) at the observer boundary -- i.e. to
% each vertex's demeaned y, before that observer's own regression is fit
% -- via the SAME retrieveObserverGainWeights2.m / groupGain-geometric-
% mean machinery already used throughout this pipeline (SPECIFICATION.md
% sec 4 calls for "gain per observer x map", which is exactly what this
% repo's gainSummaryByROI.mat already is).
%
% COMBINING OBSERVERS: fit per observer (weighted least squares, no
% intercept), THEN average across observers with EQUAL weight (settled
% spec's primary choice -- precision-weighting observers was found there
% to be ~4% worse and is not implemented here). This is a materially
% different combination than fitAsymmetryRegression.m's single pooled
% joint WLS across all subjects' rows at once.
%
% INTERVALS: 1000-draw percentile bootstrap, resampling which OBSERVERS
% contribute (not raw vertices) with replacement and averaging their
% already-fit per-observer coefficient vectors -- standard nonparametric
% bootstrap of a mean of independent per-unit estimates, and far cheaper
% than refitting the vertex-level WLS 1000 times per ROI. This replaces
% the settled spec's primary t-on-(n-1)-df interval by explicit request
% (2026-08-31); percentile bootstrap is what every other cache in this
% pipeline already uses.
%
% Saves, per ROI, under
% derivatives/summaryTables/regressionResultsHarmonic/<fitLabel>/<roi>.mat
% (a tree separate from both regressionResults/ and
% regressionResultsNative/ -- never collides with either):
%   harmonicEstimates   1x4, [b1,b2,b3,b4]*2, project-invariant order
%                        (horizontal-vertical, cardinal-oblique,
%                        radial-tangential, polar-cardinal-oblique)
%   harmonicCoeffs       4 x nBoot bootstrap draws, same order
%   harmonicTermNames    the 4 labels above, for the harmonic* fields
%   estimates, coeffs    the SAME numbers, remapped into
%                        {mainCardinal,derivedCardinal,mainSubset,
%                        derivedSubset} slot order (project-specific
%                        remap, matching fitAsymmetryRegression.m's
%                        termNames/colorKeys convention exactly) so this
%                        cache is a drop-in alternative source for
%                        plotROISummary.m/plotAsymmetryAcrossROIs.m
%   perObserverEstimates nSubj x 4, each observer's own [b1,b2,b3,b4]*2
%                        (equal-weighted mean of these rows = harmonicEstimates)
%   termNames, subjects, projectName, outputLabel, roiname, nSubj

p = inputParser;
p.addParameter('overwrite', false, @islogical);
p.addParameter('nBoot', 1000, @isnumeric);
p.addParameter('bidsDir', '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/', @ischar);
p.addParameter('githubDir', '~/Documents/GitHub', @ischar);
p.addParameter('dgSubjectMode', 'all', @(x) ismember(x, {'all','matched'}));
p.parse(varargin{:});
opt = p.Results;

githubDir = opt.githubDir;
bidsDir = opt.bidsDir;

addpath(genpath(fullfile(githubDir, 'DriftingGrating', 'AnalysisCode')));
cd(fullfile(githubDir, 'DriftingGrating', 'AnalysisCode'));
setup_user('rania', bidsDir);

projectSettings = loadConfig(githubDir);
rois = projectSettings.rois;
roi_idx = projectSettings.roi_idx;
nROIs = numel(rois);
contrasts_dict = projectSettings.contrasts_dict;
contrastnames = {contrasts_dict.contrasts.('dg_contrast_name')};

s0_idx = find(strcmp(contrastnames,'s0_v_b'));
s90_idx = find(strcmp(contrastnames,'s90_v_b'));
s45_idx = find(strcmp(contrastnames,'s45_v_b'));
s135_idx = find(strcmp(contrastnames,'s135_v_b'));
contrastIdx4 = [s0_idx, s90_idx, s45_idx, s135_idx];
id4 = [0, 90, 45, 135]; % raw stimulus code, same order as contrastIdx4 -- see REF_ORIENTATION.md sec 2 for da's rotation

harmonicTermNames = {'horizontalVsVertical','cardinalVsOblique','radialVsTangential','polarCardinalVsPolarOblique'};
harmonicDisplayLabels = {'Horizontal minus Vertical','Cardinal minus Oblique','Radial minus Tangential','Polar Cardinal minus Polar Oblique'}; % print-only, same order as harmonicTermNames
termNames = {'mainCardinal','derivedCardinal','mainSubset','derivedSubset'}; % fitAsymmetryRegression.m's slot convention, for drop-in compatibility

% Project-specific remap from harmonic b1..b4 (project-invariant concept
% order) into the mainCardinal/derivedCardinal/mainSubset/derivedSubset
% slots -- matches retrieveProConIdx.m / plotROISummary.m's colorKeys
% exactly: dg's native terms are mainCardinal=cardinal-oblique(b2),
% mainSubset=horizontal-vertical(b1); da's native terms are
% mainCardinal=polar-cardinal-oblique(b4), mainSubset=radial-tangential(b3).
if strcmp(projectName, 'dg')
    slotFromB = [2, 4, 1, 3]; % estimates(slot) = harmonicEstimates(slotFromB(slot))
elseif strcmp(projectName, 'da')
    slotFromB = [4, 2, 3, 1];
else
    error('fitHarmonicVertexModel:project', 'projectName must be ''dg'' or ''da''.');
end

if strcmp(projectName, 'dg')
    if strcmp(opt.dgSubjectMode, 'all')
        subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
            'sub-0395', 'sub-0426', 'sub-0250', 'sub-0442', 'sub-wlsubj121', 'sub-wlsubj127', ...
            'sub-0397', 'sub-0427'}; % all 13
        outputLabel = 'dg';
    else % 'matched'
        subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
            'sub-0426', 'sub-0250'}; % 7: same subjects also run in da (sub-0395 excluded)
        outputLabel = 'dgMatched7';
    end
    fullCohort = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
        'sub-0395', 'sub-0426', 'sub-0250', 'sub-0442', 'sub-wlsubj121', 'sub-wlsubj127', ...
        'sub-0397', 'sub-0427'};
elseif strcmp(projectName, 'da')
    subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
        'sub-0426', 'sub-0250'}; % 7: sub-0395 excluded (mismatched pilot stimulus)
    outputLabel = 'da';
    fullCohort = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
        'sub-0395', 'sub-0426', 'sub-0250'};
end
nSubj = numel(subjects);
subjIdx = cellfun(@(s) find(strcmp(fullCohort,s)), subjects);

saveDir = fullfile(bidsDir, 'derivatives', 'summaryTables', 'regressionResultsHarmonic', outputLabel);
if ~isfolder(saveDir), mkdir(saveDir); end

gainWeightsFile = fullfile(bidsDir, 'derivatives', 'summaryTables', 'gainSummaryByROI.mat');
Ggain = load(gainWeightsFile, 'gainTable');
gainWeightsSource = Ggain.gainTable;

glmResultsfolder = fullfile(bidsDir, 'derivatives', strcat(projectName, 'GLM'), 'hRF_glmsingle');
Sv = load(fullfile(glmResultsfolder, 'allvoxelsBOLDpa'));
Sp = load(fullfile(glmResultsfolder, 'allparamsBOLDpa'));
allvoxelsBOLDpa_full = Sv.allvoxelsBOLDpa; % contrast x wedge x roi x slot x subject
allparamsBOLDpa_full = Sp.allparamsBOLDpa; % contrast x wedge x roi x slot x [angle,ecc,R2,dist] x subject

for ri = 1:nROIs
    roiname = rois{ri};
    outFile = fullfile(saveDir, sprintf('%s.mat', roiname));
    if isfile(outFile) && ~opt.overwrite
        fprintf('%s / %s: cached harmonic fit exists, skipping (pass ''overwrite'',true to refit)\n', projectName, roiname);
        continue
    end
    roiCol = roi_idx{ri};

    gainWeights = retrieveObserverGainWeights2(subjects, roiname, gainWeightsSource);
    groupGain = exp(mean(log(gainWeights), 'omitnan'));
    subjectScale = groupGain ./ gainWeights; % 1 x nSubj

    perObserverEstimates = nan(nSubj, 4);

    for si = 1:nSubj
        fullSi = subjIdx(si);
        rowsY = []; rowsX = []; rowsW = [];

        for pa = 1:8
            % Canonical vertex-validity mask for this (wedge, roi, subject)
            % cell, read from the s0 slot -- identical across all 4
            % contrasts by construction (the inclusion mask is built from
            % pRF angle/ecc/R^2 alone, never from BOLD).
            s0vals = squeeze(allvoxelsBOLDpa_full(s0_idx, pa, roiCol, :, fullSi));
            validSlots = find(~isnan(s0vals));
            nValid = numel(validSlots);
            if nValid == 0
                continue % this observer has no usable vertices in this wedge -- contributes nothing (equal-coverage weighting has no missing-cell special case; see SPECIFICATION.md sec 5)
            end
            w = 1 / nValid; % equal-coverage: this wedge's total weight is the same as every other populated wedge

            thetaV = squeeze(allparamsBOLDpa_full(s0_idx, pa, roiCol, validSlots, 1, fullSi)); % vertex's own continuous pRF angle

            raw4 = nan(nValid, 4);
            for k = 1:4
                vals = squeeze(allvoxelsBOLDpa_full(contrastIdx4(k), pa, roiCol, validSlots, fullSi));
                raw4(:,k) = vals(:);
            end
            demeaned4 = (raw4 - mean(raw4, 2)) .* subjectScale(si); % per-vertex demeaning, then gain correction

            for k = 1:4
                if strcmp(projectName, 'dg')
                    theta = id4(k) * ones(nValid,1); % dg: raw code IS Cartesian angle everywhere
                else
                    theta = mod(id4(k) + thetaV - 90, 180); % da: UVM rotation, REF_ORIENTATION.md sec 2
                end
                rowsY = [rowsY; demeaned4(:,k)]; %#ok<AGROW>
                rowsX = [rowsX; cosd(2*theta), cosd(4*theta), cosd(2*(theta-thetaV)), cosd(4*(theta-thetaV))]; %#ok<AGROW>
                rowsW = [rowsW; w*ones(nValid,1)]; %#ok<AGROW>
            end
        end

        if isempty(rowsY)
            continue % no usable data anywhere in this ROI for this observer
        end
        lmS = fitlm(rowsX, rowsY, 'Weights', rowsW, 'Intercept', false);
        perObserverEstimates(si,:) = 2 * lmS.Coefficients.Estimate'; % 2*b = pro-minus-con scale
    end

    harmonicEstimates = mean(perObserverEstimates, 1, 'omitnan'); % equal-weighted average across observers (settled spec's primary choice)

    rng(1, 'twister');
    harmonicCoeffs = nan(4, opt.nBoot);
    validObsIdx = find(~any(isnan(perObserverEstimates), 2));
    for b = 1:opt.nBoot
        bIdx = validObsIdx(randi(numel(validObsIdx), numel(validObsIdx), 1));
        harmonicCoeffs(:,b) = mean(perObserverEstimates(bIdx,:), 1)';
    end

    estimates = nan(1,4);
    coeffs = nan(4, opt.nBoot);
    for slot = 1:4
        estimates(slot) = harmonicEstimates(slotFromB(slot));
        coeffs(slot,:) = harmonicCoeffs(slotFromB(slot),:);
    end

    % Print each asymmetry's paired-bootstrap difference (point estimate +
    % 68%/95% CI from the 1000-draw resample-across-observers bootstrap
    % above) as a table, in the project-invariant harmonicDisplayLabels
    % labeling (e.g. "Horizontal minus Vertical") so dg and da print the
    % same 4 concept labels regardless of which raw term slot each maps to.
    ci68 = prctile(harmonicCoeffs, [16 84], 2); % 4 x 2
    ci95 = prctile(harmonicCoeffs, [2.5 97.5], 2); % 4 x 2
    sig = cell(4,1);
    for hi = 1:4
        if ci95(hi,1) > 0 || ci95(hi,2) < 0
            sig{hi} = '**';
        elseif ci68(hi,1) > 0 || ci68(hi,2) < 0
            sig{hi} = '*';
        else
            sig{hi} = '';
        end
    end
    asymmetryTable = table(harmonicDisplayLabels', harmonicEstimates', ci68(:,1), ci68(:,2), ci95(:,1), ci95(:,2), sig, ...
        'VariableNames', {'Asymmetry','Estimate','CI68_lower','CI68_upper','CI95_lower','CI95_upper','Sig'});
    fprintf('  %s / %s / %s:\n', projectName, outputLabel, roiname);
    disp(asymmetryTable)

    save(outFile, 'harmonicEstimates', 'harmonicCoeffs', 'harmonicTermNames', ...
        'estimates', 'coeffs', 'termNames', 'perObserverEstimates', ...
        'subjects', 'projectName', 'outputLabel', 'roiname', 'nSubj');
    fprintf('%s (%s) / %s: harmonic fit and saved -> %s\n', projectName, outputLabel, roiname, outFile);
end

fprintf('fitHarmonicVertexModel(''%s'', ''%s''): done.\n', projectName, outputLabel);
end
