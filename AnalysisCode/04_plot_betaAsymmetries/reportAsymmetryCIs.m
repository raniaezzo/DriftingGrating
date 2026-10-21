% reportAsymmetryCIs - compile a table of the pro-vs-con beta difference
% (beta*2) plus its 68% and 95% bootstrap CIs, for every asymmetry x ROI
% x project, from the saved lme1_fit.m outputs.
%
% Requires lme1_fit.m to have already been run (with overwrite=1) through
% completion -- including the bootstrap section -- for BOTH 'dg' and 'da',
% for every ROI, so that estimates.mat and boot.mat exist under
% <glmResultsfolder>/LME_results/<comparisonName>/<ROI>/ for each project.
%
% Reproduces the same CI logic as the per-roi diagnostic printout inside
% lme1_fit.m's plotting section (CIFcn = prctile of coeffs*2), just
% compiled across both projects/all ROIs/all asymmetries into one table
% instead of being printed once per ROI during that script's run.
%
% Writes results to a table `T` (also saved as asymmetryCIs.csv in this
% script's folder) with columns:
%   project, roi, asymmetry, estimateDiff, ci95_lower, ci95_upper, ci68_lower, ci68_upper

clear; clc

bidsDir = '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
githubDir = '~/Documents/GitHub';
hRF_setting = 'glmsingle';
comparisonName = 'orientation_minus_baseline';

projectSettings = loadConfig(githubDir);
rois = projectSettings.rois;
roi_plotnames = projectSettings.roi_plotnames;

% internal predictor name -> human-readable asymmetry name, per project
% (mainSubset/derivedSubset and mainCardinal/derivedCardinal swap meaning
% between dg and da -- same swap used in lme1_fit.m's plotting section)
asymLabels.dg = struct('mainCardinal', 'Cardinal vs Oblique', ...
                        'derivedCardinal', 'Polar Cardinal vs Polar Oblique', ...
                        'mainSubset', 'Horizontal vs Vertical', ...
                        'derivedSubset', 'Radial vs Tangential');
asymLabels.da = struct('mainCardinal', 'Polar Cardinal vs Polar Oblique', ...
                        'derivedCardinal', 'Cardinal vs Oblique', ...
                        'mainSubset', 'Radial vs Tangential', ...
                        'derivedSubset', 'Horizontal vs Vertical');

CIFcn = @(x,p) prctile(x, [100-p, p]);   % matches lme1_fit.m's CIFcn

termNames = {'mainCardinal','derivedCardinal','mainSubset','derivedSubset'};
termRow   = [2, 3, 4, 5];   % row in `coeffs`/index in `estimates` for each term

project   = {};
roiCol    = {};
asymmetry = {};
estimateDiff = [];
ci95_lower = []; ci95_upper = [];
ci68_lower = []; ci68_upper = [];

for projectName = {'dg', 'da'}
    projectName = projectName{1};
    glmResultsfolder = fullfile(bidsDir, 'derivatives', strcat(projectName, 'GLM'), strcat('hRF_', hRF_setting));

    for ri = 1:length(rois)
        saveDir = fullfile(glmResultsfolder, 'LME_results', comparisonName, rois{ri});

        estFile = fullfile(saveDir, 'LME_bold.mat');
        bootFile = fullfile(saveDir, 'boot.mat');
        if ~isfile(estFile) || ~isfile(bootFile)
            warning('Missing results for project=%s roi=%s -- skipping (did lme1_fit.m finish for this project/ROI?)', ...
                projectName, rois{ri});
            continue
        end

        load(estFile, 'estimates');    %#ok<LOAD>
        load(bootFile, 'coeffs');      %#ok<LOAD>

        for ti = 1:numel(termNames)
            term = termNames{ti};
            row = termRow(ti);

            diffVal = estimates(row) * 2;
            ci95 = CIFcn(coeffs(row,:) .* 2, 97.5);
            ci68 = CIFcn(coeffs(row,:) .* 2, 84);

            project{end+1,1}   = projectName;                          %#ok<SAGROW>
            roiCol{end+1,1}    = roi_plotnames{ri};                    %#ok<SAGROW>
            asymmetry{end+1,1} = asymLabels.(projectName).(term);      %#ok<SAGROW>
            estimateDiff(end+1,1) = diffVal;                           %#ok<SAGROW>
            ci95_lower(end+1,1) = ci95(1);                             %#ok<SAGROW>
            ci95_upper(end+1,1) = ci95(2);                             %#ok<SAGROW>
            ci68_lower(end+1,1) = ci68(1);                             %#ok<SAGROW>
            ci68_upper(end+1,1) = ci68(2);                             %#ok<SAGROW>
        end
    end
end

T = table(project, roiCol, asymmetry, estimateDiff, ci95_lower, ci95_upper, ci68_lower, ci68_upper, ...
    'VariableNames', {'project','roi','asymmetry','estimateDiff','ci95_lower','ci95_upper','ci68_lower','ci68_upper'});

disp(T)

outFile = fullfile(fileparts(mfilename('fullpath')), 'asymmetryCIs.csv');
writetable(T, outFile);
fprintf('\nSaved table to %s\n', outFile);
