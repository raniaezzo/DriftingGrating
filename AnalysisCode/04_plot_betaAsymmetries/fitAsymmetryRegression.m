function fitAsymmetryRegression(projectName, varargin)
% FITASYMMETRYREGRESSION  Fits the 4 beta-asymmetries (cardinal-oblique,
% polar-cardinal-oblique, horizontal-vertical, radial-tangential) for
% every ROI via joint linear regression (no random effect), with observer
% gain-correction and observer precision-weighting applied, and caches
% the result. This is the SINGLE place the fit happens; every downstream
% plotting script (plot1/plot2, the ROI-summary figure, the dg-vs-da
% comparisons) should load this cached output rather than refitting.
%
%   fitAsymmetryRegression('dg')
%   fitAsymmetryRegression('da')
%   fitAsymmetryRegression('dg', 'overwrite', true)
%
% Fully independent of lme1_fit.m: does not read or write anything under
% LME_results/ -- lme1_fit.m and its saved data are untouched and remain
% independently runnable. Results are cached under
% derivatives/summaryTables/regressionResults/<projectName>/<roi>.mat
%
% Subjects: dg uses all 13; da uses the 7 with matched dg/da sessions
% (sub-0395 excluded -- mismatched pilot stimulus). This is the same
% exclusion already used throughout plot_NeuralAsymmetries.m and the
% validated scratch comparisons this session.
%
% Precision weighting: currently a PLACEHOLDER (uniform weights, i.e.
% mathematically an ordinary unweighted fit) via the 'precisionWeights'
% name-value pair, which defaults to a function returning all-ones. Once
% the precision-weighting method is finalized, pass a different weight
% source here and re-run with overwrite=true -- no other file needs to
% change, since everything downstream reads only the cached output.
%
% Saves, per ROI: estimates (1x4, beta*2 = pro-minus-con, same scale
% convention as lme1_fit.m's own printed stats), coeffs (4 x nBoot, the
% full bootstrap coefficient draws) so any CI level can be recomputed
% later without re-bootstrapping, and subjectContributions (nSubj x 4,
% each subject's own exact contribution to `estimates`, summing across
% subjects to reproduce it -- see the note where it's computed below for
% what this is and how it relates to simple averaging). Also saves a
% fixed subject-intercept refit (sum-to-zero/effect coding):
% grandInterceptFE (scalar), subjectFixedEffects (nSubj x 1, each
% subject's deviation from grandInterceptFE, NaN if no data in this ROI,
% sums to 0 across subjects with data), subjectIntercepts (nSubj x 1, raw
% per-subject intercepts, = grandInterceptFE + subjectFixedEffects), and
% estimatesFE (1x4, the 4 asymmetry slopes from this refit -- included
% only as a validation record, mathematically guaranteed and numerically
% verified to equal `estimates`; see the note where it's computed below).

p = inputParser;
p.addParameter('overwrite', false, @islogical);
p.addParameter('nBoot', 1000, @isnumeric);
p.addParameter('bidsDir', '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/', @ischar);
p.addParameter('githubDir', '~/Documents/GitHub', @ischar);
p.addParameter('precisionWeights', [], @(x) isempty(x) || istable(x)); % ROI-aware placeholder hook; see below
p.addParameter('dgSubjectMode', 'all', @(x) ismember(x, {'all','matched'})); % only used when projectName='dg'; 'matched' = the same 7 subjects also run in da (sub-0395 excluded), matching plot_NeuralAsymmetries.m's dg_subjectMode convention
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
elseif strcmp(projectName, 'da')
    subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
        'sub-0426', 'sub-0250'}; % 7: sub-0395 excluded (mismatched pilot stimulus)
    outputLabel = 'da';
else
    error('fitAsymmetryRegression:project', 'projectName must be ''dg'' or ''da''.');
end
nSubj = numel(subjects);

saveDir = fullfile(bidsDir, 'derivatives', 'summaryTables', 'regressionResults', outputLabel);
if ~isfolder(saveDir), mkdir(saveDir); end

% Gain is ROI-specific (see computeObserverGainWeightsByROI.m -- same
% eccentricity/R^2 vertex criteria as precision, generalized from
% dg_computeGain.m's V1-only original). The TABLE is loaded once here
% (cheap), but the actual per-subject gain VALUES are looked up per ROI
% inside the loop below via retrieveObserverGainWeights2.m -- same
% pattern as precision, and for the same reason: gain genuinely differs
% by cortical area, so computing subjectScale once outside the loop would
% silently apply one ROI's gain correction to every other ROI.
gainWeightsFile = fullfile(bidsDir, 'derivatives', 'summaryTables', 'gainSummaryByROI.mat');
Ggain = load(gainWeightsFile, 'gainTable');
gainWeightsSource = Ggain.gainTable;

% Precision weights are ROI-specific (reliability genuinely varies by
% cortical area -- see computeObserverPrecisionWeights.m), so they're
% looked up per ROI inside the loop below via
% retrieveObserverPrecisionWeights.m, not once here. opt.precisionWeights
% is currently [] everywhere this is called (PLACEHOLDER: every
% (subject, ROI) gets weight 1, making the WLS fit below mathematically
% identical to an ordinary unweighted fitlm) -- pass a table with columns
% subject/roi/weight once the method is finalized; nothing else in this
% function needs to change.

glmResultsfolder = fullfile(bidsDir, 'derivatives', strcat(projectName, 'GLM'), 'hRF_glmsingle');
S1 = load(fullfile(glmResultsfolder, 'meanBOLDpa'));
meanBOLDpa_full = S1.meanBOLDpa;

if strcmp(projectName,'dg')
    dg_full13 = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
        'sub-0395', 'sub-0426', 'sub-0250', 'sub-0442', 'sub-wlsubj121', 'sub-wlsubj127', ...
        'sub-0397', 'sub-0427'};
    subjIdx = cellfun(@(s) find(strcmp(dg_full13,s)), subjects); % identity mapping (1:13) when dgSubjectMode='all'; a 7-index subset when 'matched'
else
    da_full8 = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
        'sub-0395', 'sub-0426', 'sub-0250'};
    subjIdx = cellfun(@(s) find(strcmp(da_full8,s)), subjects);
end

for ri = 1:nROIs
    roiname = rois{ri};
    outFile = fullfile(saveDir, sprintf('%s.mat', roiname));
    if isfile(outFile) && ~opt.overwrite
        fprintf('%s / %s: cached fit exists, skipping (pass ''overwrite'',true to refit)\n', projectName, roiname);
        continue
    end

    precisionW = retrieveObserverPrecisionWeights(subjects, roiname, opt.precisionWeights)';

    gainWeights = retrieveObserverGainWeights2(subjects, roiname, gainWeightsSource);
    groupGain = exp(mean(log(gainWeights), 'omitnan')); % omitnan: a subject with no gain data for this ROI (see retrieveObserverGainWeights2.m) must not NaN out everyone else's correction factor
    subjectScale = groupGain ./ gainWeights;

    roiCol = roi_idx{ri};
    data = squeeze(meanBOLDpa_full([s0_idx,s90_idx,s45_idx,s135_idx], :, roiCol, subjIdx)); % 4 x 8 x nSubj

    % build long-format table: one row per (subject, direction, location),
    % gain-corrected; rows with NaN (no data) are dropped automatically by
    % fitlm -- this is what lets ROIs with missing locations still fit
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
    Tt = cell2table(rows, 'VariableNames', {'bold','mainCardinal','derivedCardinal','mainSubset','derivedSubset','w'});

    lm = fitlm(Tt, 'bold ~ mainCardinal + derivedCardinal + mainSubset + derivedSubset', 'Weights', Tt.w);
    coefNames = lm.CoefficientNames;
    estimates = nan(1,4);
    for ti = 1:4
        idx = strcmp(coefNames, termNames{ti});
        estimates(ti) = 2 * lm.Coefficients.Estimate(idx); % beta*2 = pro-minus-con, matches lme1_fit.m's convention
    end

    % Per-subject CONTRIBUTION to each term's coefficient (same beta*2
    % scale as `estimates`): contribution_s = (X'WX)^-1 * X_s' W_s y_s,
    % using the SAME (X'WX)^-1 from the pooled fit above (not refit per
    % subject). Summed -- not averaged -- across subjects, these exactly
    % reproduce `estimates`, for any ROI, balanced or not: this is the
    % model's own decomposition of its coefficient into per-subject
    % pieces, not each subject's independently-computed marginal
    % difference. It is well-defined even for a subject whose own data
    % alone would be too sparse to fit their own 4-term regression (their
    % individual X_s need not be full rank; only the pooled X'WX does).
    %
    % NOTE on scale: `estimates` is the group MEAN effect (a regression
    % coefficient), and subjectContributions is built so that SUMMING all
    % nSubj contributions reproduces it exactly -- so each individual
    % contribution is already scaled down by ~1/nSubj relative to that
    % subject's own marginal difference (this is correct and intentional:
    % it is literally "this subject's share of the group mean," not "this
    % subject's own independent estimate"). Anywhere a contribution is
    % used as if it were on the same scale as a raw per-subject difference
    % (e.g. plotting it as an individual line's slope, comparable to the
    % group dots), it must first be multiplied back by nSubj -- see
    % plot2_experimentalCond.m, where using this un-rescaled was a real
    % bug (found and fixed after the grey lines came out ~nSubj times too
    % shallow).
    %
    % With that rescaling (nSubj * subjectContributions(s,:)): when a
    % subject has complete data (all 8 locations, all 32 of their rows
    % valid), the four asymmetry predictors are mutually orthogonal within
    % that subject's own rows alone (a property of this design, verified
    % elsewhere in this project's analysis), so this reduces EXACTLY to
    % that subject's own (precision-)weighted marginal pro-minus-con
    % difference -- i.e. rescaled-contribution-weighting IS simple
    % precision-weighted averaging for a complete-data subject. It only
    % differs from that for subjects with missing locations, where it
    % additionally corrects for confounding between the four asymmetries
    % that the missingness can introduce.
    rowSubj = repelem((1:nSubj)', 32); % 32 rows per subject, same construction order as Tt
    validRows = ~isnan(Tt.bold);
    Xfull = [ones(height(Tt),1), Tt.mainCardinal, Tt.derivedCardinal, Tt.mainSubset, Tt.derivedSubset];
    Xv = Xfull(validRows,:);
    Wv = Tt.w(validRows);
    yv = Tt.bold(validRows);
    XtWX_inv = inv(Xv' * (Wv .* Xv)); %#ok<MINV>

    subjectContributions = nan(nSubj, 4);
    for si = 1:nSubj
        subjRows = validRows & (rowSubj == si);
        if ~any(subjRows)
            subjectContributions(si,:) = 0; % no data for this subject in this ROI -- contributes nothing
            continue
        end
        Xs = Xfull(subjRows,:);
        Ws = Tt.w(subjRows);
        ys = Tt.bold(subjRows);
        contrib = XtWX_inv * (Xs' * (Ws .* ys));
        subjectContributions(si,:) = 2 * contrib(2:5)'; % beta*2 scale, the 4 asymmetry terms only (intercept excluded -- never used downstream)
    end

    % Fixed subject intercept (effect / sum-to-zero coding): refit with a
    % separate intercept per subject instead of omitting subject entirely.
    % A per-subject-constant term is orthogonal to all 4 asymmetry
    % predictors regardless of missing data (every single location's own
    % 4 orientation rows already sum to zero for each predictor, and
    % missingness here always drops whole locations together), so this is
    % mathematically guaranteed to leave `estimates` unchanged -- verified
    % numerically below (feMismatch). B0fe is the grand intercept
    % (unweighted mean across subjects who have any data in this ROI,
    % i.e. the sum-to-zero/effect-coding convention), and
    % subjectFixedEffects(m) = C_m is subject m's deviation from B0fe,
    % with sum(C_m)=0 by construction. Subjects with zero data in this ROI
    % get NaN (their intercept is not estimable, not zero).
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

    estimatesFE = 2 * feCoefs(nSubjFE+1:nSubjFE+4)';
    feMismatch = max(abs(estimatesFE - estimates));
    if feMismatch > 1e-8
        warning('fitAsymmetryRegression:feMismatch', ...
            '%s / %s: fixed-effect slopes differ from no-subject-term slopes by %.3e (expected ~0)', ...
            projectName, roiname, feMismatch);
    end

    % bootstrap over subjects (resample WHICH subjects, refit each draw)
    rng(1, 'twister');
    coeffs = nan(4, opt.nBoot);
    for b = 1:opt.nBoot
        bIdx = randi(nSubj, nSubj, 1);
        rowSubj = repelem((1:nSubj)', 32); % 32 rows per subject (8 locations x 4 directions), same construction order as above
        % build the resampled table by concatenating each drawn subject's
        % 32 rows (with repeats, since bIdx can repeat a subject)
        Tb = table();
        for k = 1:nSubj
            Tb = [Tb; Tt(rowSubj == bIdx(k), :)]; %#ok<AGROW>
        end
        lmB = fitlm(Tb, 'bold ~ mainCardinal + derivedCardinal + mainSubset + derivedSubset', 'Weights', Tb.w);
        coefNamesB = lmB.CoefficientNames;
        for ti = 1:4
            idx = strcmp(coefNamesB, termNames{ti});
            if any(idx)
                coeffs(ti, b) = 2 * lmB.Coefficients.Estimate(idx);
            end
        end
    end

    save(outFile, 'estimates', 'coeffs', 'subjectContributions', 'termNames', 'subjects', 'projectName', 'outputLabel', 'roiname', 'nSubj', ...
        'grandInterceptFE', 'subjectFixedEffects', 'subjectIntercepts', 'estimatesFE');
    fprintf('%s (%s) / %s: fit and saved -> %s\n', projectName, outputLabel, roiname, outFile);
end

fprintf('fitAsymmetryRegression(''%s'', ''%s''): done.\n', projectName, outputLabel);
end
