% countVerticesPerROI.m
%
% For each subject and each ROI, counts how many vertices meet the same
% pRF inclusion criteria used throughout the beta-asymmetry analysis
% (pRF R^2 >= 0.1, eccentricity 4-8 deg), broken down by which of the 8
% canonical 45-deg polar-angle wedges the vertex falls into. Polar angle
% is converted from the raw Benson convention (angle_adj.mgz) via
% map_theta.m and binned via retriveRetData.m -- the same helper used by
% the ttave scripts -- so location labels here (0:45:315) match the
% convention used elsewhere (e.g. lme1_fit.m's anglevals, retriveRetData.m).
%
% ROI vertex membership mirrors 03_process_groupBetas/meanWithinLabel.m
% exactly, including its V2/V3 special case (those two ROIs are stored on
% disk as separate dorsal/ventral labels and unioned here, since there is
% no single combined V2_REmanual/V3_REmanual label file).
%
% Prints, for each ROI, every subject's per-location vertex count plus the
% all-locations total, and saves the same numbers to
% derivatives/summaryTables/vertexCountsPerROI.mat (array) and
% vertexCountsPerROI.csv (long format: subject, roi, location, nVertices).

clear; clc

githubDir = '~/Documents/GitHub';
addpath(genpath(fullfile(githubDir, 'DriftingGrating', 'AnalysisCode')));
cd(fullfile(githubDir, 'DriftingGrating', 'AnalysisCode'));
bidsDir = '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
setup_user('rania', bidsDir);

projectSettings = loadConfig(githubDir);
rois = projectSettings.rois; % {'V1','V2','V3','V3a','V3b','hV4','pMT','pMST'}

% dg's 13 subjects are a superset of da's 8 -- vertex/pRF/ROI membership
% depends only on each subject's own retinotopy data, not on which
% experiment (dg/da) is being analyzed, so one pass over this list covers
% both.
subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
    'sub-0395', 'sub-0426', 'sub-0250', 'sub-0442', 'sub-wlsubj121', 'sub-wlsubj127', ...
    'sub-0397', 'sub-0427'};

polarAngles = 0:45:315; % conventional (map_theta'd) degrees, matches retriveRetData.m/lme1_fit.m
polarAngleBinWidth = 45;
minECC = 4;
maxECC = 8;
minVAREXP = 0.1;
retFolder = 'prfvista_mov';

nSubj = numel(subjects);
nROIs = numel(rois);
nBins = numel(polarAngles);

% (location, roi, subject); last row of dim 1 is the all-locations total
countTable = nan(nBins + 1, nROIs, nSubj);

% one pRF bin-assignment pass per subject (independent of ROI), reused
% across all 8 ROIs below
filteredPrfBinsBySubj = cell(nSubj, 1);
hSizeBySubj = cell(nSubj, 1);

for si = 1:nSubj
    subjectname = subjects{si};

    subjSettings = projectSettings;
    subjSettings.subject = subjectname;
    subjSettings.bidsDir = bidsDir;
    subjSettings.retFolder = retFolder;
    subjSettings.polarAngleBinWidth = polarAngleBinWidth;
    subjSettings.minECC = minECC;
    subjSettings.maxECC = maxECC;
    subjSettings.minVAREXP = minVAREXP;

    filteredPrfBinsBySubj{si} = retriveRetData(subjSettings); % 1 x nVertices, value = assigned bin or NaN
    hSizeBySubj{si} = get_surfsize(subjectname);
end

longRows = {};

for ri = 1:nROIs
    roiname = rois{ri};
    fprintf('\n================ ROI = %s ================\n', roiname);

    for si = 1:nSubj
        subjectname = subjects{si};
        hSize = hSizeBySubj{si};

        if strcmp(roiname, 'V2') || strcmp(roiname, 'V3') % combine dorsal and ventral, as in meanWithinLabel.m
            lh_label1 = read_label(subjectname, sprintf('retinotopy_RE/lh.%sv_REmanual', roiname));
            lh_label2 = read_label(subjectname, sprintf('retinotopy_RE/lh.%sd_REmanual', roiname));
            rh_label1 = read_label(subjectname, sprintf('retinotopy_RE/rh.%sv_REmanual', roiname));
            rh_label2 = read_label(subjectname, sprintf('retinotopy_RE/rh.%sd_REmanual', roiname));
            label_idx = [lh_label1(:,1)+1; lh_label2(:,1)+1; rh_label1(:,1)+hSize(1)+1; rh_label2(:,1)+hSize(1)+1];
        else
            lh_label = read_label(subjectname, sprintf('retinotopy_RE/lh.%s_REmanual', roiname));
            rh_label = read_label(subjectname, sprintf('retinotopy_RE/rh.%s_REmanual', roiname));
            label_idx = [lh_label(:,1)+1; rh_label(:,1)+hSize(1)+1];
        end

        roiPrfBins = filteredPrfBinsBySubj{si}(label_idx);

        fprintf('%s:\n', subjectname);
        for bi = 1:nBins
            n = sum(roiPrfBins == polarAngles(bi));
            countTable(bi, ri, si) = n;
            fprintf('  location %3d: %d vertices\n', polarAngles(bi), n);
            longRows(end+1, :) = {subjectname, roiname, polarAngles(bi), n}; %#ok<SAGROW>
        end

        nAll = sum(~isnan(roiPrfBins));
        countTable(nBins + 1, ri, si) = nAll;
        fprintf('  all locations: %d vertices\n', nAll);
        longRows(end+1, :) = {subjectname, roiname, "all", nAll}; %#ok<SAGROW>
    end
end

%% Save outputs

saveDir = fullfile(bidsDir, 'derivatives', 'summaryTables');
if ~isfolder(saveDir)
    mkdir(saveDir)
end

save(fullfile(saveDir, 'vertexCountsPerROI.mat'), 'countTable', 'subjects', 'rois', 'polarAngles');

Tout = cell2table(longRows, 'VariableNames', {'subject', 'roi', 'location', 'nVertices'});
writetable(Tout, fullfile(saveDir, 'vertexCountsPerROI.csv'));

fprintf('\nSaved vertexCountsPerROI.mat and vertexCountsPerROI.csv to %s\n', saveDir);

%% Reliability summary per ROI
%
% For each ROI, reports the primary "ROI size" statistic separately from
% its location-wise spread:
%   (1) average ROI total across subjects (each subject's count summed
%       across all 8 locations, then mean/sd/min/max taken across the 13
%       subjects) -- the headline reliability number.
%   (2) the weakest and strongest visual field region, i.e. which single
%       location has the lowest/highest mean vertex count across subjects
%       -- reported separately rather than folded into one pooled
%       mean/sd, since that pooled number mixes subject-to-subject and
%       location-to-location variability together.
%   (3) the full per-location breakdown, for completeness.

fprintf('\n================ Reliability summary (all ROIs) ================\n');

for ri = 1:nROIs
    roiname = rois{ri};

    subjTotals = squeeze(countTable(nBins + 1, ri, :));

    acrossSubjByLoc = squeeze(mean(countTable(1:nBins, ri, :), 3)); % 1 mean per location, across subjects
    [minLocMean, minLocIdx] = min(acrossSubjByLoc);
    [maxLocMean, maxLocIdx] = max(acrossSubjByLoc);
    minLocSD = std(squeeze(countTable(minLocIdx, ri, :)));
    maxLocSD = std(squeeze(countTable(maxLocIdx, ri, :)));

    fprintf('\n%s:\n', roiname);
    fprintf('  average total across subjects (n=%d): mean = %.1f, sd = %.1f, min = %d, max = %d\n', ...
        nSubj, mean(subjTotals), std(subjTotals), min(subjTotals), max(subjTotals));
    fprintf('  weakest location:  %3d deg -> mean = %.1f, sd = %.1f (across %d subjects)\n', ...
        polarAngles(minLocIdx), minLocMean, minLocSD, nSubj);
    fprintf('  strongest location: %3d deg -> mean = %.1f, sd = %.1f (across %d subjects)\n', ...
        polarAngles(maxLocIdx), maxLocMean, maxLocSD, nSubj);

    fprintf('  by visual field region (across %d subjects):\n', nSubj);
    for bi = 1:nBins
        acrossSubj = squeeze(countTable(bi, ri, :));
        fprintf('    location %3d: mean = %6.1f, sd = %6.1f, min = %3d, max = %3d\n', ...
            polarAngles(bi), mean(acrossSubj), std(acrossSubj), min(acrossSubj), max(acrossSubj));
    end
end

%% Table: number of subjects with zero vertices, per (ROI, location)

zeroCounts = squeeze(sum(countTable(1:nBins, :, :) == 0, 3)); % nBins x nROIs

fprintf('\n================ Subjects with 0 vertices, per (ROI, location) ================\n');
fprintf('%-10s', 'location');
for ri = 1:nROIs
    fprintf('%8s', rois{ri});
end
fprintf('\n');
for bi = 1:nBins
    fprintf('%-10d', polarAngles(bi));
    for ri = 1:nROIs
        fprintf('%8d', zeroCounts(bi, ri));
    end
    fprintf('\n');
end

zeroCountsTable = array2table(zeroCounts, 'VariableNames', rois);
zeroCountsTable = addvars(zeroCountsTable, polarAngles', 'Before', 1, 'NewVariableNames', 'location');
writetable(zeroCountsTable, fullfile(saveDir, 'zeroVertexSubjectCounts.csv'));
fprintf('\nSaved zeroVertexSubjectCounts.csv to %s\n', saveDir);
