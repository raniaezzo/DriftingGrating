function fitAsymmetryRegression_dgVsDa(varargin)
% FITASYMMETRYREGRESSION_DGVSDA  Paired dg-vs-da comparison of the 4
% canonical asymmetries, restricted to the 7 subjects run in both
% experiments (sub-0395 excluded from da throughout this project). Each
% bootstrap draw resamples ONE set of subject indices and applies it to
% BOTH dg's and da's data in the SAME iteration, so the resulting
% dg-minus-da difference is a genuine paired comparison -- not two
% independently-bootstrapped estimates subtracted after the fact, which
% would inflate the apparent uncertainty with ordinary between-subject
% noise common to both experiments.
%
%   fitAsymmetryRegression_dgVsDa()
%   fitAsymmetryRegression_dgVsDa('overwrite', true)
%
% Independent of lme1_fit.m and of fitAsymmetryRegression.m's per-project
% caches: reads meanBOLDpa directly, writes to its own save location.
%
% dg and da each use 4 predictor columns (mainCardinal, derivedCardinal,
% mainSubset, derivedSubset), but which physical asymmetry each column
% means SWAPS between the two projects (established and verified
% elsewhere this session): dg's mainCardinal/derivedCardinal/mainSubset/
% derivedSubset = [cardinal-oblique, polar-cardinal-oblique,
% horizontal-vertical, radial-tangential]; da's same four raw slots =
% [polar-cardinal-oblique, cardinal-oblique, radial-tangential,
% horizontal-vertical]. Both are reindexed into one shared "conceptual"
% order (dgTermForConcept / daTermForConcept below) before differencing,
% so position k means the same physical asymmetry in both projects.
%
% Saves, per ROI: dg estimates/boot, da estimates/boot, and diff
% estimates/boot (dg-minus-da, paired), all in conceptual order, plus
% conceptLabels for reference. Same placeholder-precision-weight hook as
% fitAsymmetryRegression.m -- currently uniform (unweighted). Also saves,
% separately for dg and da: a fixed subject-intercept refit (sum-to-zero
% coding) -- dgGrandInterceptFE/daGrandInterceptFE (scalars),
% dgSubjectFixedEffects/daSubjectFixedEffects (nSubj x 1, each subject's
% deviation from the grand intercept, NaN if no data), and
% dgSubjectIntercepts/daSubjectIntercepts (nSubj x 1, raw per-subject
% intercepts). See fitAsymmetryRegression.m for the full derivation of why
% this is guaranteed to leave the asymmetry estimates unchanged.

p = inputParser;
p.addParameter('overwrite', false, @islogical);
p.addParameter('nBoot', 1000, @isnumeric);
p.addParameter('bidsDir', '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/', @ischar);
p.addParameter('githubDir', '~/Documents/GitHub', @ischar);
p.addParameter('precisionWeights', [], @(x) isempty(x) || istable(x)); % ROI-aware placeholder hook
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
mdirvals_dg = [0; 90; 45; 135];
anglevals = [90; 45; 0; 315; 270; 225; 180; 135];
maincardinalmDir = [0,90,180,270];
primaryMeridians = [90,0,270,180];
termNames = {'mainCardinal','derivedCardinal','mainSubset','derivedSubset'};

conceptLabels = {'Cardinal vs Oblique','Polar Cardinal vs Polar Oblique','Horizontal vs Vertical','Radial vs Tangential'};
dgTermForConcept = [1,2,3,4];
daTermForConcept = [2,1,4,3];

sharedSubjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', 'sub-0426', 'sub-0250'};
nSubj = numel(sharedSubjects);

dg_full13 = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
    'sub-0395', 'sub-0426', 'sub-0250', 'sub-0442', 'sub-wlsubj121', 'sub-wlsubj127', 'sub-0397', 'sub-0427'};
da_full8 = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', 'sub-0395', 'sub-0426', 'sub-0250'};

saveDir = fullfile(bidsDir, 'derivatives', 'summaryTables', 'regressionResults', 'dgVsDa7');
if ~isfolder(saveDir), mkdir(saveDir); end

% Gain is ROI-specific (see computeObserverGainWeightsByROI.m). The table
% is loaded once here; per-subject gain values are looked up per ROI
% inside the loop below via retrieveObserverGainWeights2.m, same pattern
% as precision -- see fitAsymmetryRegression.m for why computing this
% once outside the loop would be a bug.
gainWeightsFile = fullfile(bidsDir, 'derivatives', 'summaryTables', 'gainSummaryByROI.mat');
Ggain = load(gainWeightsFile, 'gainTable');
gainWeightsSource = Ggain.gainTable;

% Precision weights are ROI-specific -- looked up per ROI inside the loop
% below via retrieveObserverPrecisionWeights.m. opt.precisionWeights is
% currently [] (PLACEHOLDER: every (subject, ROI) gets weight 1).

Sdg = load(fullfile(bidsDir, 'derivatives', 'dgGLM', 'hRF_glmsingle', 'meanBOLDpa'));
dgSubjIdx = cellfun(@(s) find(strcmp(dg_full13,s)), sharedSubjects);
Sda = load(fullfile(bidsDir, 'derivatives', 'daGLM', 'hRF_glmsingle', 'meanBOLDpa'));
daSubjIdx = cellfun(@(s) find(strcmp(da_full8,s)), sharedSubjects);

for ri = 1:nROIs
    roiname = rois{ri};
    outFile = fullfile(saveDir, sprintf('%s.mat', roiname));
    if isfile(outFile) && ~opt.overwrite
        fprintf('dgVsDa7 / %s: cached fit exists, skipping (pass ''overwrite'',true to refit)\n', roiname);
        continue
    end

    precisionW = retrieveObserverPrecisionWeights(sharedSubjects, roiname, opt.precisionWeights)';

    gainWeights = retrieveObserverGainWeights2(sharedSubjects, roiname, gainWeightsSource);
    groupGain = exp(mean(log(gainWeights), 'omitnan')); % omitnan: see retrieveObserverGainWeights2.m
    subjectScale = groupGain ./ gainWeights; % same per-subject gain factor for both dg and da (gain is subject-level within an ROI, not experiment-specific)

    roiCol = roi_idx{ri};
    dataDg = squeeze(Sdg.meanBOLDpa([s0_idx,s90_idx,s45_idx,s135_idx], :, roiCol, dgSubjIdx));
    dataDa = squeeze(Sda.meanBOLDpa([s0_idx,s90_idx,s45_idx,s135_idx], :, roiCol, daSubjIdx));

    Tdg = buildRowTable(dataDg, 'dg', nSubj, mdirvals_dg, anglevals, maincardinalmDir, primaryMeridians, subjectScale, precisionW);
    Tda = buildRowTable(dataDa, 'da', nSubj, mdirvals_dg, anglevals, maincardinalmDir, primaryMeridians, subjectScale, precisionW);

    [dgEst, dgContrib, dgGrandInterceptFE, dgSubjectFixedEffects, dgSubjectIntercepts] = ...
        fitOneWithContributions(Tdg, termNames, nSubj, 'dg', roiname);
    [daEst, daContrib, daGrandInterceptFE, daSubjectFixedEffects, daSubjectIntercepts] = ...
        fitOneWithContributions(Tda, termNames, nSubj, 'da', roiname);
    dgConcept = dgEst(dgTermForConcept);
    daConcept = daEst(daTermForConcept);
    diffConcept = dgConcept - daConcept;

    % Per-subject contributions, reindexed to conceptual order (see
    % fitAsymmetryRegression.m for what a "contribution" is and why it
    % reduces to simple weighted averaging for a complete-data subject).
    % diffContribConcept is each subject's contribution to the PAIRED
    % (dg-minus-da) difference -- valid subject-by-subject because dg and
    % da are fit on the exact same 7 subjects here.
    dgContribConcept = dgContrib(:, dgTermForConcept);
    daContribConcept = daContrib(:, daTermForConcept);
    diffContribConcept = dgContribConcept - daContribConcept;

    rng(1, 'twister');
    nBoot = opt.nBoot;
    dgBoot = nan(4, nBoot); daBoot = nan(4, nBoot); diffBoot = nan(4, nBoot);
    rowSubjDg = repelem((1:nSubj)', 32);
    rowSubjDa = repelem((1:nSubj)', 32);
    for b = 1:nBoot
        bIdx = randi(nSubj, nSubj, 1); % SAME draw applied to both dg and da below
        TbDg = resampleRows(Tdg, rowSubjDg, bIdx);
        TbDa = resampleRows(Tda, rowSubjDa, bIdx);

        estDg = fitOne(TbDg, termNames);
        estDa = fitOne(TbDa, termNames);

        dgBoot(:,b) = estDg(dgTermForConcept);
        daBoot(:,b) = estDa(daTermForConcept);
        diffBoot(:,b) = dgBoot(:,b) - daBoot(:,b);
    end

    save(outFile, 'dgConcept', 'daConcept', 'diffConcept', 'dgBoot', 'daBoot', 'diffBoot', ...
        'dgContribConcept', 'daContribConcept', 'diffContribConcept', ...
        'dgGrandInterceptFE', 'daGrandInterceptFE', 'dgSubjectFixedEffects', 'daSubjectFixedEffects', ...
        'dgSubjectIntercepts', 'daSubjectIntercepts', ...
        'conceptLabels', 'sharedSubjects', 'roiname', 'nSubj');
    fprintf('dgVsDa7 / %s: fit and saved -> %s\n', roiname, outFile);
end

fprintf('fitAsymmetryRegression_dgVsDa: done.\n');
end

function T = buildRowTable(data, projectName, nSubj, mdirvals_dg, anglevals, maincardinalmDir, primaryMeridians, subjectScale, precisionW)
    rows = {};
    for si = 1:nSubj
        for li = 1:8
            pa = anglevals(li);
            for mi = 1:4
                md = mdirvals_dg(mi);
                val = data(mi, li, si) * subjectScale(si);

                mainCardinal = 2*ismember(md, maincardinalmDir) - 1;
                derivedCardinal = 2*((ismember(md,maincardinalmDir) & ismember(pa,primaryMeridians)) | ...
                    (~ismember(md,maincardinalmDir) & ~ismember(pa,primaryMeridians))) - 1;
                mainSubset = 0; derivedSubset = 0;
                if strcmp(projectName, 'dg')
                    if ismember(md,[0,180]); mainSubset=1; elseif ismember(md,[90,270]); mainSubset=-1; end
                    if (abs(md-pa)==0 || abs(md-pa)==180); derivedSubset=1;
                    elseif (abs(md-pa)==90 || abs(md-pa)==270); derivedSubset=-1; end
                else
                    if ismember(md,[90,270]); mainSubset=1; elseif ismember(md,[0,180]); mainSubset=-1; end
                    isCardMd = ismember(md,[0,90]); isOblMd = ismember(md,[45,135]);
                    diffv = abs(md-pa);
                    proH = (isCardMd && (diffv==90||diffv==270)) || (isOblMd && (diffv==0||diffv==180));
                    conV = (isCardMd && (diffv==0||diffv==180)) || (isOblMd && (diffv==90||diffv==270));
                    if proH; derivedSubset=1; elseif conV; derivedSubset=-1; end
                end

                rows(end+1,:) = {val, mainCardinal, derivedCardinal, mainSubset, derivedSubset, precisionW(si)}; %#ok<AGROW>
            end
        end
    end
    T = cell2table(rows, 'VariableNames', {'bold','mainCardinal','derivedCardinal','mainSubset','derivedSubset','w'});
end

function est = fitOne(T, termNames)
    lm = fitlm(T, 'bold ~ mainCardinal + derivedCardinal + mainSubset + derivedSubset', 'Weights', T.w);
    coefNames = lm.CoefficientNames;
    est = nan(4,1);
    for ti = 1:4
        idx = strcmp(coefNames, termNames{ti});
        if any(idx)
            est(ti) = 2 * lm.Coefficients.Estimate(idx);
        end
    end
end

function [est, subjectContributions, grandInterceptFE, subjectFixedEffects, subjectIntercepts] = ...
    fitOneWithContributions(T, termNames, nSubj, projectLabel, roiname)
% Same fit as fitOne, plus each subject's exact per-term contribution to
% est (nSubj x 4, beta*2 scale) -- see fitAsymmetryRegression.m for the
% full explanation of what this is, INCLUDING THE SCALE NOTE: summed
% (not averaged) across subjects this reproduces `est` exactly, so each
% individual contribution is already ~1/nSubj-scaled relative to that
% subject's own marginal difference -- multiply by nSubj before using it
% as if it were on that raw-difference scale. With that rescaling, it
% reduces to simple (precision-)weighted averaging for a complete-data
% subject. Also fits a
% fixed subject intercept (sum-to-zero/effect coding) as an alternative to
% omitting subject entirely -- mathematically guaranteed (and verified
% below) to leave `est` unchanged, since a per-subject-constant term is
% orthogonal to all 4 asymmetry predictors regardless of missing data. See
% fitAsymmetryRegression.m for the full derivation.
    est = fitOne(T, termNames);

    rowSubj = repelem((1:nSubj)', 32);
    validRows = ~isnan(T.bold);
    Xfull = [ones(height(T),1), T.mainCardinal, T.derivedCardinal, T.mainSubset, T.derivedSubset];
    Xv = Xfull(validRows,:);
    Wv = T.w(validRows);
    yv = T.bold(validRows);
    XtWX_inv = inv(Xv' * (Wv .* Xv)); %#ok<MINV>

    subjectContributions = nan(nSubj, 4);
    for si = 1:nSubj
        subjRows = validRows & (rowSubj == si);
        if ~any(subjRows)
            subjectContributions(si,:) = 0;
            continue
        end
        Xs = Xfull(subjRows,:);
        Ws = T.w(subjRows);
        ys = T.bold(subjRows);
        contrib = XtWX_inv * (Xs' * (Ws .* ys));
        subjectContributions(si,:) = 2 * contrib(2:5)';
    end

    subjHasData = false(nSubj,1);
    for si = 1:nSubj
        subjHasData(si) = any(validRows & (rowSubj == si));
    end
    subjIdxWithData = find(subjHasData);
    nSubjFE = numel(subjIdxWithData);
    rowSubjValid = rowSubj(validRows);
    dummyCols = zeros(numel(rowSubjValid), nSubjFE);
    for k = 1:nSubjFE
        dummyCols(:,k) = (rowSubjValid == subjIdxWithData(k));
    end
    Xfe = [dummyCols, Xv(:,2:5)];
    lmFE = fitlm(Xfe, yv, 'Weights', Wv, 'Intercept', false);
    feCoefs = lmFE.Coefficients.Estimate;

    subjectIntercepts = nan(nSubj,1);
    subjectIntercepts(subjIdxWithData) = feCoefs(1:nSubjFE);
    grandInterceptFE = mean(subjectIntercepts(subjIdxWithData));
    subjectFixedEffects = nan(nSubj,1);
    subjectFixedEffects(subjIdxWithData) = subjectIntercepts(subjIdxWithData) - grandInterceptFE;

    estFE = 2 * feCoefs(nSubjFE+1:nSubjFE+4)';
    feMismatch = max(abs(estFE(:) - est(:)));
    if feMismatch > 1e-8
        warning('fitAsymmetryRegression_dgVsDa:feMismatch', ...
            'dgVsDa7 / %s / %s: fixed-effect slopes differ from no-subject-term slopes by %.3e (expected ~0)', ...
            roiname, projectLabel, feMismatch);
    end
end

function Tb = resampleRows(T, rowSubj, bIdx)
    nSubj = numel(bIdx);
    Tb = table();
    for k = 1:nSubj
        Tb = [Tb; T(rowSubj == bIdx(k), :)]; %#ok<AGROW>
    end
end
