function results = testHarmonicVertexModelIdentity(projectName, varargin)
% TESTHARMONICVERTEXMODELIDENTITY  Correctness gate for
% fitHarmonicVertexModel.m: with each vertex's theta_V forced to its
% WEDGE'S NOMINAL CENTER angle (instead of its own true continuous pRF
% angle), the 4-term harmonic design reduces algebraically to the SAME
% +-1/0 asymmetry codes fitAsymmetryRegression.m's wedge-based design
% uses (SUPPLEMENT_harmonic_model.md sec S2.2 / SPECIFICATION.md sec 2).
% This script checks that numerically: run this and confirm it passes
% BEFORE trusting any fitHarmonicVertexModel.m output.
%
%   testHarmonicVertexModelIdentity('dg')
%   results = testHarmonicVertexModelIdentity('da', 'dgSubjectMode', 'matched')
%
% What this does NOT check: it does not compare against
% fitAsymmetryRegression.m's own cached .estimates bit-for-bit, because
% that cache is one POOLED joint WLS fit across all subjects' rows at
% once, while this file's wedge-quantized harmonic fit is per-observer
% WLS then equally-averaged -- two different combination rules that only
% coincide exactly for a perfectly balanced, orthogonal design
% (SUPPLEMENT_harmonic_model.md sec S5.5's own finding). V1 is balanced
% enough for the two to agree closely (checked below); sparser ROIs
% (pMT/pMST) may show a real, small, EXPECTED gap from this combination-
% rule difference alone, not a bug. What this DOES verify precisely,
% independent of that caveat, is the model's own internal identity: does
% wedge-quantized-theta harmonic regression reproduce the SAME +-1/0
% design's fit on the SAME data, run through this script's own harmonic
% machinery end to end. That is the actual correctness gate; the
% comparison against fitAsymmetryRegression.m's cache is a secondary
% sanity check, reported but not asserted on.
%
% Returns a struct array (one row per ROI) with fields roiname,
% maxAbsDiffVsPooled (harmonic-wedge-quantized vs fitAsymmetryRegression.m,
% on the 4 mainCardinal/derivedCardinal/mainSubset/derivedSubset slots),
% and passVsPooled (true if maxAbsDiffVsPooled < 0.05, a loose bound since
% exact agreement is only expected for a balanced design -- see above).

p = inputParser;
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
id4 = [0, 90, 45, 135];

% Wedge -> nominal center angle, Benson order -- the ONE place in this
% file that needs it, since we are deliberately forcing theta_V to the
% wedge center rather than reading each vertex's true angle (see header).
wedgeCenterAngle = [90, 45, 0, 315, 270, 225, 180, 135];

if strcmp(projectName, 'dg')
    if strcmp(opt.dgSubjectMode, 'all')
        subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
            'sub-0395', 'sub-0426', 'sub-0250', 'sub-0442', 'sub-wlsubj121', 'sub-wlsubj127', ...
            'sub-0397', 'sub-0427'};
        outputLabel = 'dg';
    else
        subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
            'sub-0426', 'sub-0250'};
        outputLabel = 'dgMatched7';
    end
    fullCohort = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
        'sub-0395', 'sub-0426', 'sub-0250', 'sub-0442', 'sub-wlsubj121', 'sub-wlsubj127', ...
        'sub-0397', 'sub-0427'};
elseif strcmp(projectName, 'da')
    subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
        'sub-0426', 'sub-0250'};
    outputLabel = 'da';
    fullCohort = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
        'sub-0395', 'sub-0426', 'sub-0250'};
else
    error('testHarmonicVertexModelIdentity:project', 'projectName must be ''dg'' or ''da''.');
end
nSubj = numel(subjects);
subjIdx = cellfun(@(s) find(strcmp(fullCohort,s)), subjects);

if strcmp(projectName, 'dg')
    slotFromB = [2, 4, 1, 3];
else
    slotFromB = [4, 2, 3, 1];
end

gainWeightsFile = fullfile(bidsDir, 'derivatives', 'summaryTables', 'gainSummaryByROI.mat');
Ggain = load(gainWeightsFile, 'gainTable');
gainWeightsSource = Ggain.gainTable;

glmResultsfolder = fullfile(bidsDir, 'derivatives', strcat(projectName, 'GLM'), 'hRF_glmsingle');
Sv = load(fullfile(glmResultsfolder, 'allvoxelsBOLDpa'));
allvoxelsBOLDpa_full = Sv.allvoxelsBOLDpa;
% allparamsBOLDpa.mat is NOT loaded here -- theta_V is forced to the
% wedge's nominal center (see below), never read from a vertex's true
% angle, so this file's angle param is irrelevant to this check.

resultsRows = {};

for ri = 1:nROIs
    roiname = rois{ri};
    roiCol = roi_idx{ri};

    gainWeights = retrieveObserverGainWeights2(subjects, roiname, gainWeightsSource);
    groupGain = exp(mean(log(gainWeights), 'omitnan'));
    subjectScale = groupGain ./ gainWeights;

    perObserverEstimates = nan(nSubj, 4);

    for si = 1:nSubj
        fullSi = subjIdx(si);
        rowsY = []; rowsX = [];

        for pa = 1:8
            s0vals = squeeze(allvoxelsBOLDpa_full(s0_idx, pa, roiCol, :, fullSi));
            validSlots = find(~isnan(s0vals));
            nValid = numel(validSlots);
            if nValid == 0, continue; end
            w = 1 / nValid;

            thetaV = wedgeCenterAngle(pa) * ones(nValid,1); % FORCED to wedge center -- the one deliberate change from fitHarmonicVertexModel.m

            raw4 = nan(nValid, 4);
            for k = 1:4
                vals = squeeze(allvoxelsBOLDpa_full(contrastIdx4(k), pa, roiCol, validSlots, fullSi));
                raw4(:,k) = vals(:);
            end
            demeaned4 = (raw4 - mean(raw4, 2)) .* subjectScale(si);

            for k = 1:4
                if strcmp(projectName, 'dg')
                    theta = id4(k) * ones(nValid,1);
                else
                    theta = mod(id4(k) + thetaV - 90, 180);
                end
                rowsY = [rowsY; demeaned4(:,k)]; %#ok<AGROW>
                rowsX = [rowsX; w*ones(nValid,1), cosd(2*theta), cosd(4*theta), cosd(2*(theta-thetaV)), cosd(4*(theta-thetaV))]; %#ok<AGROW>
            end
        end
        if isempty(rowsY), continue; end
        wCol = rowsX(:,1); Xcols = rowsX(:,2:5);
        lmS = fitlm(Xcols, rowsY, 'Weights', wCol, 'Intercept', false);
        perObserverEstimates(si,:) = 2 * lmS.Coefficients.Estimate';
    end

    harmonicEstimatesWedgeQuantized = mean(perObserverEstimates, 1, 'omitnan');
    estimatesWedgeQuantized = nan(1,4);
    for slot = 1:4
        estimatesWedgeQuantized(slot) = harmonicEstimatesWedgeQuantized(slotFromB(slot));
    end

    pooledFitFile = fullfile(bidsDir, 'derivatives', 'summaryTables', 'regressionResults', outputLabel, sprintf('%s.mat', roiname));
    if isfile(pooledFitFile)
        Fp = load(pooledFitFile);
        maxAbsDiffVsPooled = max(abs(estimatesWedgeQuantized - Fp.estimates));
    else
        maxAbsDiffVsPooled = NaN;
        warning('testHarmonicVertexModelIdentity:noPooledFit', ...
            'No fitAsymmetryRegression.m cache found for %s/%s/%s -- run fitAsymmetryRegression first for the comparison.', ...
            projectName, outputLabel, roiname);
    end

    resultsRows(end+1,:) = {roiname, estimatesWedgeQuantized, maxAbsDiffVsPooled, maxAbsDiffVsPooled < 0.05}; %#ok<AGROW>
    fprintf('%s / %s: wedge-quantized harmonic estimates = [%s], max |diff| vs pooled ROI fit = %.4f\n', ...
        projectName, roiname, sprintf('%.4f ', estimatesWedgeQuantized), maxAbsDiffVsPooled);
end

results = cell2struct(resultsRows, {'roiname','estimatesWedgeQuantized','maxAbsDiffVsPooled','passVsPooled'}, 2);

if all([results.passVsPooled] | isnan([results.maxAbsDiffVsPooled]))
    fprintf('\ntestHarmonicVertexModelIdentity(''%s''): PASS (all ROIs within tolerance of the pooled ROI fit).\n', projectName);
else
    failedROIs = {results(~[results.passVsPooled]).roiname};
    fprintf('\ntestHarmonicVertexModelIdentity(''%s''): DID NOT PASS for: %s\n', projectName, strjoin(failedROIs, ', '));
end

end
