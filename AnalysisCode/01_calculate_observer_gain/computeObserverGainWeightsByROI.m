% computeObserverGainWeightsByROI.m
%
% Generalizes dg_computeGain.m from V1-only to all 8 ROIs: computes
% mean/median pRF gain per (observer, cortical area), using the SAME
% vertex inclusion criteria (eccentricity 4-8 deg, pRF R^2 >= 0.1, read
% from prfvista_mov's lh/rh .eccen.mgz/.vexpl.mgz maps, regardless of
% which protocol gain is computed from) as computeObserverPrecisionWeights.m
% -- so gain-weighting and precision-weighting for a given (subject, ROI)
% are always computed over the identical vertex set. Same ROI-label
% reading (V2/V3 split into dorsal+ventral REmanual labels, union'd) and
% same left-hemisphere-then-right-hemisphere vertex ordering as that
% script too.
%
% GAIN DEFINITION: unchanged from dg_computeGain.m -- rmModelGain's
% 'maxabs' metric (largest absolute excursion of the model's predicted
% time series from baseline, in percent BOLD), computed separately for the
% moving (prfvista_mov) and stationary (prfvista) carrier and then
% averaged (meanGain_avg/medianGain_avg), per (subject, ROI).
%
% SPEED: one rmModelGain call per subject per protocol (mov, stat), over
% the UNION of all 8 ROIs' vertices (unfiltered by pRF criteria at this
% stage) -- not one call per ROI -- then the pRF-criteria mask and
% per-ROI membership are applied afterward to aggregate into per-ROI
% mean/median. Same two-stage pattern as computeObserverPrecisionWeights.m
% (union-then-filter), which is what lets this stay fast despite covering
% 8 ROIs.
%
% OUTPUTS (saved to derivatives/summaryTables/):
%   gainSummaryByROI.mat - gainTable (long format: subject, roi, weight
%                           [=meanGain_avg, the column
%                           retrieveObserverGainWeights2.m reads],
%                           meanGain_mov, medianGain_mov, meanGain_stat,
%                           medianGain_stat, medianGain_avg, nVertices),
%                           subjects, rois, minECC, maxECC, minVAREXP, metric
%   gainSummaryByROI.csv - same table
%
% See also DG_COMPUTEGAIN (the original V1-only version),
% COMPUTEOBSERVERPRECISIONWEIGHTS (the ROI-vertex-selection pattern this
% mirrors), RETRIEVEOBSERVERGAINWEIGHTS2 (reads gainTable's 'weight' column)

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

minECC = 4; maxECC = 8; minVAREXP = 0.1; % matches computeObserverPrecisionWeights.m exactly
metric = 'maxabs'; % matches dg_computeGain.m's default

movRoot = fullfile(bidsDir, 'derivatives', 'prfvista_mov');
statRoot = fullfile(bidsDir, 'derivatives', 'prfvista');
saveDir = fullfile(bidsDir, 'derivatives', 'summaryTables');

longRows = {};

for si = 1:nSubj
    subjectname = subjects{si};

    retDir = dir(fullfile(movRoot, subjectname, '**/stimfiles.mat'));
    if isempty(retDir)
        warning('%s: no prfvista_mov data -- skipped entirely', subjectname);
        continue
    end
    retDir = retDir(1).folder;
    hemis = {'lh','rh'};
    ret = struct();
    for hi = 1:numel(hemis)
        hemi = hemis{hi};
        ret.(sprintf('%s_ecc', hemi)) = MRIread(fullfile(retDir, sprintf('%s.eccen.mgz', hemi)));
        ret.(sprintf('%s_vexp', hemi)) = MRIread(fullfile(retDir, sprintf('%s.vexpl.mgz', hemi)));
    end
    eccAll = [ret.lh_ecc.vol, ret.rh_ecc.vol];
    r2All = [ret.lh_vexp.vol, ret.rh_vexp.vol];
    nTotal = numel(eccAll);
    includedAll = (eccAll(:) >= minECC) & (eccAll(:) <= maxECC) & (r2All(:) >= minVAREXP);

    hSize = get_surfsize(subjectname);

    % union of ROI vertex indices, same V2/V3 dorsal+ventral special-casing
    % as computeObserverPrecisionWeights.m, plus each ROI's own row
    % positions into that union for later per-ROI aggregation
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

    if nTotal < max(allIdx)
        error('computeObserverGainWeightsByROI:size', '%s: eccen/vexpl maps have %d vertices, ROI index max is %d.', ...
            subjectname, nTotal, max(allIdx));
    end

    roiRows = cell(nROIs,1);
    for ri = 1:nROIs
        [~, roiRows{ri}] = ismember(roiVertexIdx{ri}, allIdx);
    end

    t0 = tic;
    gainMov = gainFromProtocol(movRoot, subjectname, metric, allIdx, nTotal);
    gainStat = gainFromProtocol(statRoot, subjectname, metric, allIdx, nTotal);
    fprintf('%s: gain computed over %d ROI-union vertices, %.0f s\n', subjectname, numel(allIdx), toc(t0));

    for ri = 1:nROIs
        rows = roiRows{ri};
        rows = rows(rows > 0);
        if isempty(rows)
            continue
        end
        vIdx = allIdx(rows);
        incl = includedAll(vIdx);
        rows = rows(incl);
        nVertsUsed = numel(rows);
        if nVertsUsed == 0
            continue
        end

        gMov = gainMov(rows);
        gStat = gainStat(rows);
        meanGain_mov = mean(gMov);
        medianGain_mov = median(gMov);
        meanGain_stat = mean(gStat);
        medianGain_stat = median(gStat);
        meanGain_avg = mean([meanGain_mov, meanGain_stat]);
        medianGain_avg = mean([medianGain_mov, medianGain_stat]);

        longRows(end+1, :) = {subjectname, rois{ri}, meanGain_avg, meanGain_mov, medianGain_mov, ...
            meanGain_stat, medianGain_stat, medianGain_avg, nVertsUsed}; %#ok<SAGROW>
    end
    fprintf('%s: gain summarized for %d/%d ROIs\n', subjectname, ...
        sum(strcmp(longRows(:,1), subjectname)), nROIs);
end

%% Save

gainTable = cell2table(longRows, 'VariableNames', {'subject','roi','weight','meanGain_mov','medianGain_mov', ...
    'meanGain_stat','medianGain_stat','medianGain_avg','nVertices'});

save(fullfile(saveDir, 'gainSummaryByROI.mat'), 'gainTable', 'subjects', 'rois', 'minECC', 'maxECC', 'minVAREXP', 'metric');
writetable(gainTable, fullfile(saveDir, 'gainSummaryByROI.csv'));

fprintf('\nSaved gainSummaryByROI.mat and .csv to %s\n', saveDir);


% ------------------------------------------------------------------------
function gain = gainFromProtocol(root, sub, metric, voxels, expectedNVerts)
% Same rm-struct-building and rmModelGain call as dg_computeGain.m's
% computeAndSaveGain, minus the save-to-disk step (this script only needs
% the in-memory gain vector, not a per-subject gain_<sub>.mat file).
hits = dir(fullfile(root, sub, 'ses-*', '*results.mat'));
if isempty(hits)
    error('no *results.mat under %s', fullfile(root, sub));
end
if numel(hits) > 1
    fprintf('    %d results files under %s, using %s\n', numel(hits), root, hits(1).name);
end
srcFile = fullfile(hits(1).folder, hits(1).name);

V = load(srcFile);
[model, params] = findModelParamsLocal(V);
if ~isfield(params, 'analysis') || ~isfield(params.analysis, 'allstimimages')
    error(['params.analysis.allstimimages is missing in %s, so the predicted ' ...
        'time series cannot be built.'], srcFile);
end

m = model; if iscell(m), m = m{1}; end
nV = numel(m.x0);
if nV ~= expectedNVerts
    error('%s has %d vertices, expected %d (from prfvista_mov''s mgz maps)', srcFile, nV, expectedNVerts);
end

rm.model = model;
rm.params = params;
gain = rmModelGain(rm, 'voxels', voxels, 'metric', metric);
gain = double(gain(:));
end


% ------------------------------------------------------------------------
function [model, params] = findModelParamsLocal(V)
% The results file may store model/params at the top level, or nested
% inside a variable such as 'results'. Handle both (same as dg_computeGain.m).
if isfield(V, 'model') && isfield(V, 'params')
    model = V.model; params = V.params; return
end
top = fieldnames(V);
for k = 1:numel(top)
    n = V.(top{k});
    if isstruct(n) && isscalar(n) && isfield(n, 'model') && isfield(n, 'params')
        model = n.model; params = n.params; return
    end
end
error(['could not find model and params in the results file. ' ...
    'Top-level variables were: %s'], strjoin(top', ', '));
end
