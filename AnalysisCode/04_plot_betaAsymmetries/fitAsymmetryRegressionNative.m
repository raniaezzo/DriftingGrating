function fitAsymmetryRegressionNative(projectName, varargin)
% FITASYMMETRYREGRESSIONNATIVE  Native-only counterpart to
% fitAsymmetryRegression.m, for plotROISummary.m/
% plotAsymmetryAcrossROIs.m's nativeOnly=true mode.
%
%   fitAsymmetryRegressionNative('dg')
%   fitAsymmetryRegressionNative('da')
%   fitAsymmetryRegressionNative('dg', 'dgSubjectMode', 'matched')
%
% Reads computeMeanBOLDNative.m's meanBOLDnative.mat (per (contrast, roi,
% subject) vertex mean within ecc[0.5,12], no polar-angle restriction, no
% varexp filter) instead of meanBOLDpa.mat. No regression/design-matrix
% fit is needed here: mainCardinal and mainSubset -- the two "native"
% terms in fitAsymmetryRegression.m's 4-term model -- are defined purely
% from stimulus direction, never from location (see that file's
% mainCardinal/mainSubset assignment, which never references `pa`; the
% derivedCardinal/derivedSubset terms are the ones that need a location
% argument, which is exactly why only those two are "derived" --
% REF_ORIENTATION.md sec 4). So the "fit" here is just each subject's own
% (gain-corrected) pro-minus-con scalar, combined across subjects via a
% precision-weighted mean, with a paired subject-resample bootstrap for
% CI -- same bootstrap convention (resample subjects with replacement,
% nBoot=1000, rng(1,'twister')) as fitAsymmetryRegression.m, just without
% the WLS refit step since there is no design matrix to refit.
%
% mainCardinal (both projects): pro=mean(s0_v_b,s90_v_b), con=mean(s45_v_b,s135_v_b).
% mainSubset: dg pro=s0_v_b(horizontal), con=s90_v_b(vertical);
%             da pro=s90_v_b(radial), con=s0_v_b(tangential).
% (Matches retrieveProConIdx.m's orientation_minus_baseline native branch
% exactly.)
%
% Saves the SAME 4-slot cache format as fitAsymmetryRegression.m
% (estimates 1x4, coeffs 4 x nBoot, termNames = {mainCardinal,
% derivedCardinal, mainSubset, derivedSubset}) so
% plotROISummary.m/plotAsymmetryAcrossROIs.m's nativeOnly branch can
% reuse the exact same colorKeys/xLabels lookup code as the non-native
% path, just restricted to termIdx in {1,3} -- termIdx {2,4}
% (derivedCardinal, derivedSubset) are intentionally left NaN; nativeOnly
% never reads them. Cached under
% derivatives/summaryTables/regressionResultsNative/<fitLabel>/<roi>.mat
% -- a SEPARATE tree from fitAsymmetryRegression.m's
% regressionResults/<fitLabel>/, so neither can ever overwrite the other.

p = inputParser;
p.addParameter('overwrite', false, @islogical);
p.addParameter('nBoot', 1000, @isnumeric);
p.addParameter('bidsDir', '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/', @ischar);
p.addParameter('githubDir', '~/Documents/GitHub', @ischar);
p.addParameter('precisionWeights', [], @(x) isempty(x) || istable(x));
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
termNames = {'mainCardinal','derivedCardinal','mainSubset','derivedSubset'}; % same 4-slot layout as fitAsymmetryRegression.m; only 1 and 3 are ever filled here

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
    error('fitAsymmetryRegressionNative:project', 'projectName must be ''dg'' or ''da''.');
end
nSubj = numel(subjects);

saveDir = fullfile(bidsDir, 'derivatives', 'summaryTables', 'regressionResultsNative', outputLabel);
if ~isfolder(saveDir), mkdir(saveDir); end

gainWeightsFile = fullfile(bidsDir, 'derivatives', 'summaryTables', 'gainSummaryByROI.mat');
Ggain = load(gainWeightsFile, 'gainTable');
gainWeightsSource = Ggain.gainTable;

glmResultsfolder = fullfile(bidsDir, 'derivatives', strcat(projectName, 'GLM'), 'hRF_glmsingle');
S1 = load(fullfile(glmResultsfolder, 'meanBOLDnative'));
meanBOLDnative_full = S1.meanBOLDnative;

if strcmp(projectName,'dg')
    dg_full13 = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
        'sub-0395', 'sub-0426', 'sub-0250', 'sub-0442', 'sub-wlsubj121', 'sub-wlsubj127', ...
        'sub-0397', 'sub-0427'};
    subjIdx = cellfun(@(s) find(strcmp(dg_full13,s)), subjects);
else
    da_full8 = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
        'sub-0395', 'sub-0426', 'sub-0250'};
    subjIdx = cellfun(@(s) find(strcmp(da_full8,s)), subjects);
end

for ri = 1:nROIs
    roiname = rois{ri};
    outFile = fullfile(saveDir, sprintf('%s.mat', roiname));
    if isfile(outFile) && ~opt.overwrite
        fprintf('%s / %s: cached native fit exists, skipping (pass ''overwrite'',true to refit)\n', projectName, roiname);
        continue
    end

    % Gain and precision weights are ROI-specific, looked up fresh each
    % iteration -- same reasoning and same helper calls as
    % fitAsymmetryRegression.m (computing them once outside this loop
    % would silently apply one ROI's correction to every other ROI).
    precisionW = retrieveObserverPrecisionWeights(subjects, roiname, opt.precisionWeights); % 1 x nSubj

    gainWeights = retrieveObserverGainWeights2(subjects, roiname, gainWeightsSource); % 1 x nSubj
    groupGain = exp(mean(log(gainWeights), 'omitnan'));
    subjectScale = groupGain ./ gainWeights; % 1 x nSubj

    roiCol = roi_idx{ri};
    s0 = squeeze(meanBOLDnative_full(s0_idx, roiCol, subjIdx))' .* subjectScale;
    s90 = squeeze(meanBOLDnative_full(s90_idx, roiCol, subjIdx))' .* subjectScale;
    s45 = squeeze(meanBOLDnative_full(s45_idx, roiCol, subjIdx))' .* subjectScale;
    s135 = squeeze(meanBOLDnative_full(s135_idx, roiCol, subjIdx))' .* subjectScale;

    diffMainCardinal = mean([s0; s90], 1) - mean([s45; s135], 1); % 1 x nSubj, cardinal/polar-cardinal minus oblique/polar-oblique
    if strcmp(projectName, 'dg')
        diffMainSubset = s0 - s90; % horizontal minus vertical
    else
        diffMainSubset = s90 - s0; % radial minus tangential
    end

    estimates = nan(1,4);
    estimates(1) = weightedMeanRow(diffMainCardinal, precisionW);
    estimates(3) = weightedMeanRow(diffMainSubset, precisionW);

    rng(1, 'twister');
    coeffs = nan(4, opt.nBoot);
    for b = 1:opt.nBoot
        bIdx = randi(nSubj, nSubj, 1);
        coeffs(1,b) = weightedMeanRow(diffMainCardinal(bIdx), precisionW(bIdx));
        coeffs(3,b) = weightedMeanRow(diffMainSubset(bIdx), precisionW(bIdx));
    end

    save(outFile, 'estimates', 'coeffs', 'termNames', 'subjects', 'projectName', 'outputLabel', 'roiname', 'nSubj');
    fprintf('%s (%s) / %s: native fit and saved -> %s\n', projectName, outputLabel, roiname, outFile);
end

fprintf('fitAsymmetryRegressionNative(''%s'', ''%s''): done.\n', projectName, outputLabel);
end

function m = weightedMeanRow(x, w)
% 1xN weighted mean, dropping any entry where x or w is NaN (a subject
% missing gain/precision data for this ROI) from both the numerator and
% denominator.
    valid = ~isnan(x) & ~isnan(w);
    m = sum(x(valid) .* w(valid)) / sum(w(valid));
end
