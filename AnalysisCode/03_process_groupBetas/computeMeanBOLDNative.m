function computeMeanBOLDNative(projectName)
% COMPUTEMEANBOLDNATIVE  Per-(contrast,roi,subject) BOLD means for the
% "native" asymmetry pathway (nativeOnly=true in
% plotROISummary.m/plotAsymmetryAcrossROIs.m, wired via
% fitAsymmetryRegressionNative.m): averages ALL vertices in each ROI
% label within eccentricity [0.5,12] deg, with NO polar-angle restriction
% (every visual-field location contributes -- each native asymmetry is a
% single stimulus/pair of stimuli that is the same physical thing at
% every location, by construction; see REF_ORIENTATION.md secs 4-5) and
% NO variance-explained filter (explicit 2026-08-31 decision -- a
% deliberate departure from meanWithinLabel.m's meanBOLDpa, which uses
% ecc[4,8], varexp>=0.1, AND bins by polar angle; none of that machinery
% applies here).
%
% Only computes the 4 raw per-direction contrasts the native pathway ever
% reads (s0_v_b/s90_v_b/s45_v_b/s135_v_b -- see retrieveProConIdx.m's
% orientation_minus_baseline branch), not all 29 contrasts
% meanWithinLabel.m computes. Since eccentricity is the only inclusion
% criterion (no varexp/sigma/angle needed), this only loads each
% subject's eccen.mgz, not the full pRF stack -- much cheaper than
% meanWithinLabel.m's loop, and independent of it (does not read or
% overwrite meanBOLD.mat/meanBOLDpa.mat or the ROI masks that script
% writes).
%
% Output meanBOLDnative/medianBOLDnative are (contrast, roi, subject),
% indexed by the FULL contrastnames-length first dimension (all
% non-native entries left NaN) so downstream code can look up contrasts
% by name via the same `find(strcmp(contrastnames,...))` convention used
% throughout this pipeline, rather than a separate remapped index.
%
%   computeMeanBOLDNative('dg')
%   computeMeanBOLDNative('da')

bidsDir =  '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
githubDir = '~/Documents/GitHub';

addpath(genpath(fullfile(githubDir, 'DriftingGrating', 'AnalysisCode')));
addpath(genpath(fullfile(githubDir, 'atlasmgz')));
setup_user('rania', bidsDir);

hRF_setting = 'glmsingle';
hemis = {'lh'; 'rh'};

projectSettings = loadConfig(githubDir);
rois = projectSettings.rois;
contrasts_dict = projectSettings.contrasts_dict;
contrastnames = {contrasts_dict.contrasts.('dg_contrast_name')};

nativeContrastNames = {'s0_v_b','s90_v_b','s45_v_b','s135_v_b'};
nativeContrastIdx = cellfun(@(n) find(strcmp(contrastnames,n)), nativeContrastNames);

if strcmp(projectName, 'dg')
    subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', ...
        'sub-wlsubj124', 'sub-0395', 'sub-0426', 'sub-0250', ...
        'sub-0442', 'sub-wlsubj121', 'sub-wlsubj127',  'sub-0397', ...
        'sub-0427'};
elseif strcmp(projectName, 'da')
    subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
        'sub-0395', 'sub-0426', 'sub-0250'};
else
    error('computeMeanBOLDNative:project', 'projectName must be ''dg'' or ''da''.');
end

minECC_native = 0.5;
maxECC_native = 12;

meanBOLDnative = nan(length(contrastnames), length(rois), length(subjects));
medianBOLDnative = nan(length(contrastnames), length(rois), length(subjects));

for si = 1:numel(subjects)

    subjectname = subjects{si};
    fprintf('%s\n', subjectname);

    glmResultsfolder = fullfile(bidsDir, 'derivatives', strcat(projectName, 'GLM'), strcat('hRF_', hRF_setting), subjectname);
    glmFilelist = dir(fullfile(glmResultsfolder, '**', 'results.mat'));
    betaResults = load(fullfile(glmFilelist.folder, glmFilelist.name));
    betaResults = betaResults.results;

    movingRetDir = dir(fullfile(bidsDir, 'derivatives', 'prfvista_mov', subjectname, '**/stimfiles.mat'));
    movingRetDir = movingRetDir.folder;

    ret_ecc = struct();
    for hi = 1:numel(hemis)
        hemi = hemis{hi};
        ret_ecc.(hemi) = MRIread(fullfile(movingRetDir, sprintf('%s.eccen.mgz', hemi)));
    end
    eccData = [ret_ecc.lh.vol, ret_ecc.rh.vol]; % 1 x (lh+rh vertices)

    hSize = get_surfsize(subjectname);

    for ri = 1:numel(rois)

        roiname = rois{ri};

        if strcmp(roiname, 'V2') || strcmp(roiname, 'V3') % combine dorsal and ventral
            lh_label1 = read_label(subjectname, sprintf('retinotopy_RE/lh.%sv_REmanual', roiname));
            lh_label2 = read_label(subjectname, sprintf('retinotopy_RE/lh.%sd_REmanual', roiname));
            rh_label1 = read_label(subjectname, sprintf('retinotopy_RE/rh.%sv_REmanual', roiname));
            rh_label2 = read_label(subjectname, sprintf('retinotopy_RE/rh.%sd_REmanual', roiname));
            label_idx = [lh_label1(:,1)+1 ; lh_label2(:,1)+1; rh_label1(:,1)+hSize(1)+1; rh_label2(:,1)+hSize(1)+1];
        else
            lh_label = read_label(subjectname, sprintf('retinotopy_RE/lh.%s_REmanual', roiname));
            rh_label = read_label(subjectname, sprintf('retinotopy_RE/rh.%s_REmanual', roiname));
            label_idx = [lh_label(:,1)+1 ; rh_label(:,1)+hSize(1)+1];
        end

        eccBool = (eccData >= minECC_native) & (eccData <= maxECC_native);
        validIdx = find(eccBool);
        finalIdx = intersect(label_idx, validIdx);

        for ci = 1:numel(nativeContrastIdx)
            contrastGlobalIdx = nativeContrastIdx(ci);
            contrastname = contrastnames{contrastGlobalIdx};
            propername = strrep(contrastname, '_v_','V');
            currBold = betaResults.contrasts.(propername);

            meanBOLDnative(contrastGlobalIdx, ri, si) = mean(currBold(finalIdx));
            medianBOLDnative(contrastGlobalIdx, ri, si) = median(currBold(finalIdx));
        end
    end
end

saveFolder = fullfile(bidsDir, 'derivatives', strcat(projectName, 'GLM'), strcat('hRF_', hRF_setting));
save(fullfile(saveFolder, 'meanBOLDnative.mat'), 'meanBOLDnative');
save(fullfile(saveFolder, 'medianBOLDnative.mat'), 'medianBOLDnative');

fprintf('saved meanBOLDnative, medianBOLDnative to:\n  %s\n', saveFolder);

end
