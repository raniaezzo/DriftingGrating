clc; clear all; close all;

% lme1_fit_compareContexts.m
%
% Compares the LME beta-asymmetry estimates between the dg (Cartesian) and
% da (Polar) experiments, using the PAIRED bootstrap structure already
% built into lme1_fit.m: both projects' bootstrap sections call
% rng('default'); rng(1) at the same point (immediately before the
% roi/bootstrap-iteration loop), with nothing project-dependent consuming
% random numbers in between, and datasample(1:8, 8) is called with
% identical inputs at each (roi, bootstrap-iteration) -- so bootstrap
% iteration bi resamples the SAME 8 subjects, in the SAME order, in both
% the dg run and the da run.
%
% This was verified against 03_process_groupBetas/meanWithinLabel.m (the
% script that builds each project's meanBOLDpa.mat/meanBOLD.mat, which
% lme1_fit.m's modeldata is ultimately derived from): its subjects list
% for 'da' -- {sub-0037, sub-0201, sub-0255, sub-wlsubj123, sub-wlsubj124,
% sub-0395, sub-0426, sub-0250} -- is exactly the first 8 (and in the same
% order) of its 13-subject 'dg' list, and the subject dimension of
% meanBOLDpa/meanBOLD is filled by iterating that list in order (si). So
% numeric subject index i corresponds to the same physical subject in
% both projects' data, not just the same index.
%
% For each of the 4 canonical asymmetries (vertical-vs-horizontal,
% cardinal-vs-oblique, radial-vs-tangential, polar-cardinal-vs-polar-
% oblique), the plotted quantity is the cross-context (dg - da) difference
% of each project's BETA ONLY (intercept excluded, same as lme1_fit.m's
% own meanRelative=1 convention) -- "pro" = beta_dg - beta_da (e.g.
% horizontal_dg - horizontal_da), "con" is exactly its negative (e.g.
% vertical_dg - vertical_da). Since con = -pro by construction (the
% intercept cancels, just as pro/con within a single project are always
% +-beta around the intercept), the two plotted series -- and their
% bootstrap CIs -- are always mirror images about 0, matching how the
% single-context plots are symmetric. Pro and con are still shown as two
% separate series (not collapsed), so they can cross zero and each other,
% they just can never be individually asymmetric about 0 as a pair.
%
% The separately-printed "stats" table instead reports the FULL asymmetry
% comparison per canonical asymmetry -- (horizontal-vertical) for dg vs.
% (horizontal-vertical) for da, i.e. 2*(beta_dg - beta_da) -- a single
% number per asymmetry, matching lme1_fit.m's own within-project stats
% convention (which reports pro-minus-con as beta*2).
%
% Requires that lme1_fit.m has already been run (with overwrite=1) for
% BOTH projectName='dg' and projectName='da', so each ROI's LME_bold.mat/
% boot.mat already exist on disk -- this script only loads those, it does
% not refit or re-bootstrap anything.

%% SET UP (same config as lme1_fit.m)

bidsDir = '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
githubDir = '~/Documents/GitHub';
addpath(genpath(fullfile(githubDir, 'atlasmgz')));
hRF_setting = 'glmsingle';

% can be 'motion_minus_orientation' ; 'motion_minus_baseline' ; 'orientation_minus_baseline'
comparisonName = 'orientation_minus_baseline';
metric = 'bold';

projectSettings = loadConfig(githubDir);
rois = projectSettings.rois; % rois{1} = 'V1' (see ROIS_ALL.json) -- included throughout below
colors_data = projectSettings.colors_data;

figureDir = strrep(bidsDir, 'data_bids', 'figures');
compareFigureDir = fullfile(figureDir, 'dg_vs_da');
if ~isfolder(compareFigureDir)
    mkdir(compareFigureDir)
end

ci_level = 84; % 68% CI -- same as the error bars lme1_fit.m plots (not its separate 95%/68% "stats" printout)

nROIs = numel(rois);
projects = {'dg', 'da'};

for pi = 1:numel(projects)
    checkFolder = fullfile(bidsDir, 'derivatives', strcat(projects{pi}, 'GLM'), strcat('hRF_', hRF_setting), ...
        'LME_results', comparisonName, rois{1});
    if ~isfile(fullfile(checkFolder, 'boot.mat'))
        error(['No saved boot.mat found for project "%s" (%s).\n' ...
            'Run lme1_fit.m with projectName=''%s'' and overwrite=1 first.'], ...
            projects{pi}, checkFolder, projects{pi});
    end
end

% Canonical asymmetry order + labels/colors, matching lme1_fit.m's master
% figure exactly (see its plotOrder): raw coeffs/estimates rows 2:5 are
% [mainCardinalVsMainOblique, derivedCardinalVsDerivedOblique, mainSubset,
% derivedSubset] for BOTH projects, but which raw slot means which
% semantic asymmetry swaps between dg and da -- these plotOrders undo that
% swap so canonical position k means the same physical asymmetry (and the
% same "pro" direction) in both projects. Colors are keyed by dg's own
% (canonical, unswapped) COLORS.json names -- verified identical to da's
% colors for the same semantic asymmetry, so it's safe to always read them
% from the dg namespace regardless of which project's data is being
% colored. canonicalProConLabels give the actual pro/con physical
% direction per asymmetry (verified against lme1_fit.m's own pro/con index
% definitions), used for plot legends.
dgPlotOrder = [3, 1, 4, 2];
daPlotOrder = [4, 2, 3, 1];
canonicalColorKeys    = {'verticalVsHorizontal', 'mainCardinalVsMainOblique', 'radialVsTangential', 'derivedCardinalVsDerivedOblique'};
canonicalDisplayLabels = {'verticalVsHorizontal', 'cardinalVsOblique', 'radialVsTangential', 'polarCardinalVsPolarOblique'};
canonicalProConLabels = {
    {'horizontal', 'vertical'}
    {'cardinal', 'oblique'}
    {'radial', 'tangential'}
    {'polar cardinal', 'polar oblique'}
    };
nA = numel(canonicalDisplayLabels);

%% LOAD both projects' saved LME fit + bootstrap results

rawEstimates = struct(); % per project: [5 x nROIs] (row 1 = intercept, rows 2:5 = betas, raw order)
rawCoeffs = struct();    % per project: {roi} -> [5 x nBoot]

for pi = 1:numel(projects)
    projectName = projects{pi};
    glmResultsfolder = fullfile(bidsDir, 'derivatives', strcat(projectName, 'GLM'), strcat('hRF_', hRF_setting));

    rawEstimates.(projectName) = nan(5, nROIs);
    rawCoeffs.(projectName) = cell(nROIs, 1);

    for roi = 1:nROIs
        saveDir = fullfile(glmResultsfolder, 'LME_results', comparisonName, rois{roi});

        load(fullfile(saveDir, strcat('LME_', metric)), 'estimates');
        load(fullfile(saveDir, 'boot'), 'coeffs');

        rawEstimates.(projectName)(:, roi) = estimates(:);
        rawCoeffs.(projectName){roi} = coeffs;
    end
end

nBoot = size(rawCoeffs.dg{1}, 2);
if size(rawCoeffs.da{1}, 2) ~= nBoot
    error('dg and da have different bootstrap sample counts (%d vs %d) -- the paired-bootstrap comparison below assumes they match.', ...
        nBoot, size(rawCoeffs.da{1}, 2));
end

%% For each project/ROI, pull out just the beta (NOT intercept-adjusted --
% same convention as lme1_fit.m's meanRelative=1 plots, where "pro"/"con"
% are plotted as +beta/-beta, not Gintercept+-beta), reindexed into
% canonical asymmetry order

plotOrderByProject = struct('dg', dgPlotOrder, 'da', daPlotOrder);

beta_point = struct(); % per project: [nA x nROIs]
beta_boot  = struct(); % per project: {roi} -> [nA x nBoot]

for pi = 1:numel(projects)
    projectName = projects{pi};
    thisPlotOrder = plotOrderByProject.(projectName);

    beta_point.(projectName) = rawEstimates.(projectName)(2:5, :);
    beta_point.(projectName) = beta_point.(projectName)(thisPlotOrder, :); % nA x nROIs, canonical order

    beta_boot.(projectName) = cell(nROIs, 1);
    for roi = 1:nROIs
        thisBoot = rawCoeffs.(projectName){roi}(2:5, :); % 4 x nBoot, raw order
        beta_boot.(projectName){roi} = thisBoot(thisPlotOrder, :); % nA x nBoot, canonical order
    end
end

%% Cross-context (dg - da) difference of the beta ALONE, per ROI (paired:
% same bootstrap iteration index = same resampled subjects in both
% projects). "pro" = this difference; "con" = its exact negative, so the
% two are always symmetric about 0.

proDiff_point = beta_point.dg - beta_point.da; % nA x nROIs
conDiff_point = -proDiff_point;

proDiff_boot = cell(nROIs, 1);
conDiff_boot = cell(nROIs, 1);

proDiff_errlow  = nan(nA, nROIs); proDiff_errhigh = nan(nA, nROIs);
conDiff_errlow  = nan(nA, nROIs); conDiff_errhigh = nan(nA, nROIs);

CIFcn = @(x, p) prctile(x, [100-p, p]);

for roi = 1:nROIs
    proDiff_boot{roi} = beta_boot.dg{roi} - beta_boot.da{roi}; % nA x nBoot
    conDiff_boot{roi} = -proDiff_boot{roi};

    for ai = 1:nA
        ciP = CIFcn(proDiff_boot{roi}(ai, :), ci_level);
        proDiff_errlow(ai, roi)  = proDiff_point(ai, roi) - ciP(1);
        proDiff_errhigh(ai, roi) = ciP(2) - proDiff_point(ai, roi);

        ciC = CIFcn(conDiff_boot{roi}(ai, :), ci_level);
        conDiff_errlow(ai, roi)  = conDiff_point(ai, roi) - ciC(1);
        conDiff_errhigh(ai, roi) = ciC(2) - conDiff_point(ai, roi);
    end
end

%% Stats printout: the FULL asymmetry comparison per canonical asymmetry --
% (pro-minus-con) for dg vs. (pro-minus-con) for da, i.e. 2*proDiff -- a
% single number per asymmetry, matching lme1_fit.m's own within-project
% "beta*2" stats convention. Reported for ROI 1 (V1).

fprintf('\nstats (roi=%s): CI of the cross-context (dg - da) difference in pro-minus-con effect (beta*2)\n', rois{1});
fprintf('%-28s %12s %24s %24s\n', 'asymmetry', 'estimate', '95% CI', '68% CI');
for ai = 1:nA
    ci95 = CIFcn(proDiff_boot{1}(ai, :) .* 2, 97.5);
    ci68 = CIFcn(proDiff_boot{1}(ai, :) .* 2, 84);
    fprintf('%-28s %12.5f   [%9.5f %9.5f]   [%9.5f %9.5f]\n', canonicalDisplayLabels{ai}, ...
        proDiff_point(ai, 1) * 2, ci95(1), ci95(2), ci68(1), ci68(2));
end

%% Master figure (ROI 1 = V1): x-axis = the 4 canonical asymmetries, pro
% and con cross-context differences (mirror images about 0) plotted side
% by side at each position -- same formatting as lme1_fit.m's own master
% figure

x = 1:nA;

figure; hold on

legendHandles = gobjects(2, 1);

for ai = 1:nA
    colP = colors_data.conditions.dg.(canonicalColorKeys{ai}).color_pro';
    colC = colors_data.conditions.dg.(canonicalColorKeys{ai}).color_con';

    boxchart(x(ai), proDiff_point(ai, 1), 'BoxFaceColor', colP, 'LineWidth', 4, 'BoxWidth', 0.6);
    hold on
    h1 = errorbar(x(ai), proDiff_point(ai, 1), proDiff_errlow(ai, 1), proDiff_errhigh(ai, 1), ...
        'LineStyle', 'none', 'LineWidth', 2, 'Color', colP);

    boxchart(x(ai), conDiff_point(ai, 1), 'BoxFaceColor', colC, 'LineWidth', 4, 'BoxWidth', 0.6);
    hold on
    h2 = errorbar(x(ai), conDiff_point(ai, 1), conDiff_errlow(ai, 1), conDiff_errhigh(ai, 1), ...
        'LineStyle', 'none', 'LineWidth', 2, 'Color', colC);

    if ai == 1
        legendHandles(1) = h1;
        legendHandles(2) = h2;
    end
end

yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2)

ylim([-0.4 0.4])
xlim([0.5 nA+0.5])

set(gca, 'XTick', [])
box off
set(gca, 'linewidth', 2, 'YColor', [0 0 0], 'XColor', [0 0 0]);
set(gca, 'FontName', 'Arial', 'FontSize', 20);

xticks(x)
xticklabels(canonicalDisplayLabels)
xtickangle(25)

ylabel('\Delta (dg - da) beta', 'FontSize', 20);
title(rois{1}, 'FontSize', 20, 'Interpreter', 'none');

legend(legendHandles, {'pro (e.g. horizontal/cardinal/radial/polar-cardinal)', 'con (e.g. vertical/oblique/tangential/polar-oblique)'}, ...
    'Location', 'best', 'FontSize', 10);

f = gcf;
f.Position = [298 843 651 494];

drawnow;
print(fullfile(compareFigureDir, sprintf('LME_crossContext_%s_%s', comparisonName, rois{1})), '-dpdf', '-bestfit');

%% Per-asymmetry figures: one figure per canonical asymmetry, x-axis = ALL
% ROIs (V1 included), pro/con cross-context differences (mirror images
% about 0) plotted side by side at each ROI position -- mirrors
% lme1_fit.m's per-asymmetryName (first) plot type, generalized across
% both projects

for ai = 1:nA
    colP = colors_data.conditions.dg.(canonicalColorKeys{ai}).color_pro';
    colC = colors_data.conditions.dg.(canonicalColorKeys{ai}).color_con';
    proLabel = canonicalProConLabels{ai}{1};
    conLabel = canonicalProConLabels{ai}{2};

    figure; hold on

    for roi = 1:nROIs
        boxchart(roi, proDiff_point(ai, roi), 'BoxFaceColor', colP, 'LineWidth', 4, 'BoxWidth', 0.75);
        hold on
        h1 = errorbar(roi, proDiff_point(ai, roi), proDiff_errlow(ai, roi), proDiff_errhigh(ai, roi), ...
            'LineStyle', 'none', 'LineWidth', 2, 'Color', colP);

        boxchart(roi, conDiff_point(ai, roi), 'BoxFaceColor', colC, 'LineWidth', 4, 'BoxWidth', 0.75);
        hold on
        h2 = errorbar(roi, conDiff_point(ai, roi), conDiff_errlow(ai, roi), conDiff_errhigh(ai, roi), ...
            'LineStyle', 'none', 'LineWidth', 2, 'Color', colC);
    end

    yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2)

    ylim([-0.4 0.4])
    xlim([0 nROIs+1])

    box off
    set(gca, 'linewidth', 2, 'YColor', [0 0 0], 'XColor', [0 0 0]);
    set(gca, 'FontName', 'Arial', 'FontSize', 20);

    xticks(1:nROIs)
    xticklabels(projectSettings.roi_plotnames)
    xtickangle(25)

    ylabel('\Delta (dg - da) beta', 'FontSize', 20);
    title(canonicalDisplayLabels{ai}, 'FontSize', 20, 'Interpreter', 'none');

    legend([h1, h2], {proLabel, conLabel}, 'Location', 'best');

    f = gcf;
    f.Position = [298 843 651 494];

    drawnow;
    print(fullfile(compareFigureDir, sprintf('LME_crossContext_%s_%s', comparisonName, canonicalDisplayLabels{ai})), '-dpdf', '-bestfit');

end

%% Cross-context difference of the WITHIN-project pro-minus-con (2*beta)
% value, ROI 1 (V1): dg's pro-minus-con MINUS da's pro-minus-con, plotted
% directly as a single combined number per canonical asymmetry (matching
% the printed stats table above, just visualized). Error bars are the CI
% of this same paired-difference bootstrap distribution (proDiff_boot is
% already beta_dg - beta_da per iteration), not derived independently
% from dg's and da's separate CIs.

dgMinusDa_point = proDiff_point(:, 1) .* 2; % ROI 1 = V1, nA x 1 (pro minus con = beta*2)
dgMinusDa_boot  = proDiff_boot{1} .* 2;     % nA x nBoot

dgMinusDa_errlow  = nan(nA, 1);
dgMinusDa_errhigh = nan(nA, 1);
for ai = 1:nA
    ci = CIFcn(dgMinusDa_boot(ai, :), ci_level);
    dgMinusDa_errlow(ai)  = dgMinusDa_point(ai) - ci(1);
    dgMinusDa_errhigh(ai) = ci(2) - dgMinusDa_point(ai);
end

figure; hold on

for ai = 1:nA
    col = colors_data.conditions.dg.(canonicalColorKeys{ai}).color_pro';

    boxchart(x(ai), dgMinusDa_point(ai), 'BoxFaceColor', col, 'LineWidth', 4, 'BoxWidth', 0.6);
    hold on
    errorbar(x(ai), dgMinusDa_point(ai), dgMinusDa_errlow(ai), dgMinusDa_errhigh(ai), ...
        'LineStyle', 'none', 'LineWidth', 2, 'Color', col);
end

yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2)

ylim([-0.4 0.4])
xlim([0.5 nA+0.5])

set(gca, 'XTick', [])
box off
set(gca, 'linewidth', 2, 'YColor', [0 0 0], 'XColor', [0 0 0]);
set(gca, 'FontName', 'Arial', 'FontSize', 20);

xticks(x)
xticklabels(canonicalDisplayLabels)
xtickangle(25)

ylabel('\Delta (dg - da) pro-minus-con (beta*2)', 'FontSize', 20);
title(rois{1}, 'FontSize', 20, 'Interpreter', 'none');

f = gcf;
f.Position = [298 843 651 494];

drawnow;
print(fullfile(compareFigureDir, sprintf('LME_proMinusConDiff_%s_%s', comparisonName, rois{1})), '-dpdf', '-bestfit');
