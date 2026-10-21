clc; clear all; close all;

% lme1_fit_perSubjectContextDiagnostic.m
%
% Diagnostic for lme1_fit_compareContexts.m: for a given canonical
% asymmetry (default: radialVsTangential) and ROI (default: V1), shows
% each of the 8 subjects' OWN cross-context (dg - da) difference, so an
% outlier subject driving/suppressing the group-level effect can be
% spotted directly, rather than only seeing the group estimate.
%
% Uses each project's saved modeldata.mat (written by lme1_fit.m,
% row = subject x polarangle x motiondir, already gain-weighted per
% subject -- see lme1_fit.m's subjectScale). For each subject, computes
% the raw per-subject mean BOLD across "pro" rows (e.g. radial) minus
% "con" rows (e.g. tangential) for dg and for da separately, then takes
% that subject's own (dg - da) difference. This is a simpler quantity than
% the LME's beta (a plain per-subject mean contrast, not a mixed-model
% estimate), but that's the point here -- it's meant to be transparent
% enough to eyeball which subject is atypical, not to replace the LME.
%
% Requires lme1_fit.m already run (with overwrite=1) for both projects, so
% each ROI's modeldata.mat exists on disk.

%% SET UP

bidsDir = '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
githubDir = '~/Documents/GitHub';
addpath(genpath(fullfile(githubDir, 'atlasmgz')));
hRF_setting = 'glmsingle';
comparisonName = 'orientation_minus_baseline';

projectSettings = loadConfig(githubDir);
rois = projectSettings.rois;

figureDir = strrep(bidsDir, 'data_bids', 'figures');
compareFigureDir = fullfile(figureDir, 'dg_vs_da');
if ~isfolder(compareFigureDir)
    mkdir(compareFigureDir)
end

% same 8 subjects, same order, used by lme1_fit.m for both dg and da (see
% lme1_fit_compareContexts.m's header comment for how this was verified)
subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
    'sub-0395', 'sub-0426', 'sub-0250'};
nSubjs = numel(subjects);

% which canonical asymmetry to inspect, and in which ROI
asymmetryToCheck = 'verticalVsHorizontal'; % 'verticalVsHorizontal' | 'cardinalVsOblique' | 'radialVsTangential' | 'polarCardinalVsPolarOblique'
roiToCheck = 'V1';

% maps each canonical asymmetry to the raw modeldata column that encodes
% it in each project (same mapping as lme1_fit_compareContexts.m's
% dgPlotOrder/daPlotOrder, expressed directly in column-name form; +1 row
% = "pro" e.g. radial/horizontal/cardinal/polar-cardinal, -1 = "con")
colByProject = struct( ...
    'verticalVsHorizontal',        struct('dg', 'mainSubset',      'da', 'derivedSubset'), ...
    'cardinalVsOblique',           struct('dg', 'mainCardinal',    'da', 'derivedCardinal'), ...
    'radialVsTangential',          struct('dg', 'derivedSubset',   'da', 'mainSubset'), ...
    'polarCardinalVsPolarOblique', struct('dg', 'derivedCardinal', 'da', 'mainCardinal') ...
    );

if ~isfield(colByProject, asymmetryToCheck)
    error('Unknown asymmetryToCheck "%s".', asymmetryToCheck);
end
dgCol = colByProject.(asymmetryToCheck).dg;
daCol = colByProject.(asymmetryToCheck).da;

%% LOAD modeldata for both projects, this ROI only

projects = {'dg', 'da'};
modeldataByProject = struct();

for pi = 1:numel(projects)
    projectName = projects{pi};
    glmResultsfolder = fullfile(bidsDir, 'derivatives', strcat(projectName, 'GLM'), strcat('hRF_', hRF_setting));
    saveDir = fullfile(glmResultsfolder, 'LME_results', comparisonName, roiToCheck);

    modeldataFile = fullfile(saveDir, 'modeldata.mat');
    if ~isfile(modeldataFile)
        error(['No saved modeldata.mat found for project "%s", ROI "%s" (%s).\n' ...
            'Run lme1_fit.m with projectName=''%s'' and overwrite=1 first.'], ...
            projectName, roiToCheck, saveDir, projectName);
    end

    loaded = load(modeldataFile, 'modeldata');
    modeldataByProject.(projectName) = loaded.modeldata;
end

%% Per-subject pro-minus-con raw mean, per project, then per-subject (dg - da)

subjVal_dg = nan(nSubjs, 1);
subjVal_da = nan(nSubjs, 1);

for si = 1:nSubjs
    md = modeldataByProject.dg;
    proRows = md.sub == si & md.(dgCol) == 1;
    conRows = md.sub == si & md.(dgCol) == -1;
    subjVal_dg(si) = mean(md.bold(proRows)) - mean(md.bold(conRows));

    md = modeldataByProject.da;
    proRows = md.sub == si & md.(daCol) == 1;
    conRows = md.sub == si & md.(daCol) == -1;
    subjVal_da(si) = mean(md.bold(proRows)) - mean(md.bold(conRows));
end

subjDiff = subjVal_dg - subjVal_da; % nSubjs x 1: this subject's own (dg - da) cross-context difference
groupMeanDiff = mean(subjDiff);
groupMedianDiff = median(subjDiff);

%% Print a sorted-by-deviation table so the outlier is easy to spot

fprintf('\nPer-subject (dg - da) difference in %s pro-minus-con, ROI=%s\n', asymmetryToCheck, roiToCheck);
fprintf('(dg column: %s, da column: %s; group mean = %.5f, group median = %.5f)\n', dgCol, daCol, groupMeanDiff, groupMedianDiff);
fprintf('%-16s %12s %12s %12s %14s\n', 'subject', 'dg (pro-con)', 'da (pro-con)', 'dg - da', 'dev. from mean');

[~, sortIdx] = sort(abs(subjDiff - groupMeanDiff), 'descend');
for k = 1:nSubjs
    si = sortIdx(k);
    fprintf('%-16s %12.5f %12.5f %12.5f %14.5f\n', subjects{si}, subjVal_dg(si), subjVal_da(si), ...
        subjDiff(si), subjDiff(si) - groupMeanDiff);
end

%% Plot: one point/bar per subject, ordered as in `subjects` (not
% re-sorted, so it's easy to cross-reference against other figures), with
% the group mean/median shown for reference

figure; hold on

barColor = [0.3 0.3 0.3];
bar(1:nSubjs, subjDiff, 'FaceColor', barColor, 'EdgeColor', 'none', 'BarWidth', 0.6);

yline(0, '-', 'Color', [0 0 0], 'LineWidth', 1.5);
yline(groupMeanDiff, '--', 'Color', [0.85 0.1 0.1], 'LineWidth', 2);

xlim([0.5 nSubjs+0.5])
xticks(1:nSubjs)
xticklabels(subjects)
xtickangle(35)

box off
set(gca, 'linewidth', 2, 'YColor', [0 0 0], 'XColor', [0 0 0]);
set(gca, 'FontName', 'Arial', 'FontSize', 16);

ylim([-0.8 0.8])
ylabel('mean asymmetry (Cartesian minus polar)', 'FontSize', 18);
title(sprintf('%s, %s: per-subject cross-context difference', asymmetryToCheck, roiToCheck), ...
    'FontSize', 16, 'Interpreter', 'none');

%legend({'per-subject', 'zero', 'group mean'}, 'Location', 'best');

f = gcf;
f.Position = [298 843 750 494];

drawnow;
print(fullfile(compareFigureDir, sprintf('perSubject_crossContext_%s_%s', asymmetryToCheck, roiToCheck)), '-dpdf', '-bestfit');
