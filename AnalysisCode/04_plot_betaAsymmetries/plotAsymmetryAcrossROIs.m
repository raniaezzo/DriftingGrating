function plotAsymmetryAcrossROIs(projectSettings, varargin)
% PLOTASYMMETRYACROSSROIS  Successor to lme1_fit.m's "plot across ROIs
% per figure" visualization: one figure PER ASYMMETRY, with cortical
% areas along the x-axis (the inverse layout of plotROISummary.m, which
% is one figure per cortical area with asymmetries on the x-axis).
% Sourced from fitAsymmetryRegression.m's cached joint regression fit
% instead of the LME. Fully independent of lme1_fit.m -- does not read or
% write anything under LME_results/; lme1_fit.m and its output remain
% untouched and independently runnable.
%
%   plotAsymmetryAcrossROIs(projectSettings)
%   plotAsymmetryAcrossROIs(projectSettings, 'nativeOnly', true)
%
% projectSettings - same struct plot1/plot2_experimentalCond.m and
%                    plotROISummary.m use (must include .rois, .roi_idx,
%                    .comparisonName, .projectName, .colors_data,
%                    .figureDir, .gainSummaryFile). fitLabel: which
%                    regressionResults/<label>/ subfolder to read --
%                    defaults to projectName (dg's own full 13-subject
%                    cache, da's own 7 -- verified against the cached
%                    .mat files' nSubj field), but set
%                    projectSettings.fitLabel = 'dgMatched7' explicitly
%                    for dg's matched-to-da 7-subject variant, same
%                    convention plot1/plot2_experimentalCond.m use.
%
% Pro/con points are MODEL-DERIVED, zero-centered (pro = +beta,
% con = -beta) -- pro is always FILLED, con always UNFILLED (white face,
% colored edge), matching plot2_experimentalCond.m's 'pairwise' dot
% convention exactly (proFaceColor/proEdgeColor/conFaceColor/conEdgeColor,
% meanDotSize, proLineWidth, errorbarLineWidth -- all read from the same
% COLORS.json styleInfo fields, not re-derived or hand-tuned here). Error
% bars are the 68% CI of the cached fit's bootstrapped difference, full
% width, applied symmetrically to both -- unchanged from before (matches
% lme1_fit.m's ACTUAL plotted behavior given its meanRelative=1 setting,
% and the full-width CI approximates the conventional ~95% "non-overlap"
% significance heuristic). Significance asterisks at a fixed y=0.35 above
% each cortical area's point: ** if the 95% CI of the difference excludes
% 0, * if only the 68% CI does.
%
% Axes styling (line width/color, dashed zero line, physical cm sizing)
% also matches plot2_experimentalCond.m's pairwise convention: the SAME
% cm-per-x-unit density plot2 uses for its own pro/con spacing
% (~1.7625 cm/unit) is applied here across however many ROI x-positions
% this figure has, rather than a fixed on-screen pixel Position + bestfit
% export.
%
% nativeOnly (optional name-value, default false): when true, restricts
% to the 2 asymmetries that are "native" to this project (dg: Horizontal
% vs Vertical, Cardinal vs Oblique; da: Radial vs Tangential, Polar
% Cardinal vs Polar Oblique -- see REF_ORIENTATION.md sec 4), sourced from
% fitAsymmetryRegressionNative.m's cache (vertex mean within ecc[0.5,12],
% no polar-angle binning, no varexp filter -- a genuinely different
% analysis, not just a display filter) instead of
% fitAsymmetryRegression.m's. One figure is produced per native
% asymmetry (2 total) instead of per all 4.

p = inputParser;
p.addParameter('nativeOnly', false, @islogical);
p.parse(varargin{:});
nativeOnly = p.Results.nativeOnly;

projectName = projectSettings.projectName;
comparisonName = projectSettings.comparisonName;
colors_data = projectSettings.colors_data;
rois = projectSettings.rois;
figureDir = projectSettings.figureDir;

[summaryTablesDir,~,~] = fileparts(projectSettings.gainSummaryFile);
[derivativesDir,~,~] = fileparts(summaryTablesDir);
[bidsDir,~,~] = fileparts(derivativesDir);

% Concept label and COLORS.json style key per raw term (termIdx 1..4 =
% mainCardinal, derivedCardinal, mainSubset, derivedSubset) -- same
% mapping as plotROISummary.m's colorKeys, plus the human-readable
% concept label used there for xLabels but here indexed directly by
% termIdx (no cross-asymmetry x-axis reordering needed here, since each
% asymmetry gets its own figure).
if strcmp(projectName, 'dg')
    conceptLabels = {'Cardinal vs Oblique','Polar Cardinal vs Polar Oblique','Horizontal vs Vertical','Radial vs Tangential'};
    colorKeys = {'mainCardinalVsMainOblique','derivedCardinalVsDerivedOblique','verticalVsHorizontal','radialVsTangential'};
elseif strcmp(projectName, 'da')
    conceptLabels = {'Polar Cardinal vs Polar Oblique','Cardinal vs Oblique','Radial vs Tangential','Horizontal vs Vertical'};
    colorKeys = {'mainCardinalVsMainOblique','derivedCardinalVsDerivedOblique','radialVsTangential','verticalVsHorizontal'};
else
    error('plotAsymmetryAcrossROIs:project', 'projectName must be ''dg'' or ''da''.');
end

nROIs = length(rois);
statsRows = {};

% Load all cortical areas' cached fits once, up front (each figure below
% needs all of them, not just one). See header for the fitLabel default.
if isfield(projectSettings, 'fitLabel') && ~isempty(projectSettings.fitLabel)
    fitLabel = projectSettings.fitLabel;
else
    fitLabel = projectName;
end
if nativeOnly
    fitSubdir = 'regressionResultsNative';
else
    fitSubdir = 'regressionResults';
end
F = cell(nROIs,1);
for ri = 1:nROIs
    fitFile = fullfile(bidsDir,'derivatives','summaryTables',fitSubdir,fitLabel,sprintf('%s.mat',rois{ri}));
    if ~isfile(fitFile)
        if nativeOnly
            error('plotAsymmetryAcrossROIs:missingFit', ...
                'Cached native fit not found for %s / %s at %s -- run fitAsymmetryRegressionNative(''%s'') first.', ...
                fitLabel, rois{ri}, fitFile, projectName);
        else
            error('plotAsymmetryAcrossROIs:missingFit', ...
                'Cached fit not found for %s / %s at %s -- run fitAsymmetryRegression(''%s'') first.', ...
                fitLabel, rois{ri}, fitFile, projectName);
        end
    end
    F{ri} = load(fitFile);
end

% Styling constants, matching plot2_experimentalCond.m's pairwise
% convention exactly (same source fields, same formula).
axisLineWidth = 1;
axisLineColor = [0.25 0.25 0.25];
polarMarkerSize = 6 * 0.8;
meanDotSize = (pi/4) * polarMarkerSize^2;
dashLen = 0.12; gapLen = 0.06; % same NaN-segmented dash construction as plot2's subjectwiseDiff zero line

% Physical axis sizing: reuse plot2_experimentalCond.m's own
% cm-per-x-unit density exactly (derived the same way it is there), so
% ROI spacing here reads at the same visual scale as pro/con spacing in
% the pairwise plots, just extended across nROIs positions instead of 2.
previousXlim = [0.5, 2.5];
previousPairwisePlotWidth_cm = 3 * 1.25 * 0.94;
xScale_cm_per_unit = previousPairwisePlotWidth_cm / diff(previousXlim);
pairwisePlotHeight_cm = 2 * (3.6 * 0.94); % 2x plot2's own per-subplot height per feedback -- axis box itself (not just the figure canvas) is taller; ylim stays [-0.4,0.4] unchanged, so the same data range is stretched over more physical height
% padding_cm and extraLeftForYLabel_cm are ~2x their original values (2/4.5
% vs the prior 1/2.2) to accommodate the tick-label/asterisk font size
% doubling (8->16 / 10->20) requested 2026-08-31.
padding_cm = 2; % page margin on top/right/bottom (room for rotated ROI tick labels, title, and now-larger asterisks)
% Extra left margin -- the two-line ylabel ("BOLD response amplitude" /
% "(% signal change)") needs more room than the flat padding_cm gives it,
% same fix and same value as plotROISummary.m's extraLeftForXTickLabels_cm.
% Also needs to fit the (now larger) Y-tick labels themselves.
extraLeftForYLabel_cm = 4.5;

% nativeOnly restricts to termIdx {1,3} (mainCardinal, mainSubset) -- the
% two terms fitAsymmetryRegression.m's design never conditions on
% location (see that file's mainCardinal/mainSubset assignment) -- one
% figure per native asymmetry instead of per all 4.
if nativeOnly
    termIdxList = [1, 3];
else
    termIdxList = 1:4;
end

for termIdx = termIdxList
    figure; hold on

    styleInfo = colors_data.conditions.(projectName).(colorKeys{termIdx});
    proFaceColor = styleInfo.color_pro';
    conFaceColor = [1 1 1];
    proEdgeColor = styleInfo.color_pro';
    conEdgeColor = styleInfo.color_con'; % full-strength con color, no white blending
    proLineWidth = styleInfo.pro_lineWidth;
    errorbarLineWidth = styleInfo.errorbar_lineWidth / 2; % half width per feedback

    xlimRange = [0.5, nROIs+0.5];
    zeroLineX = buildDashedZeroLine(xlimRange, dashLen, gapLen);
    plot(zeroLineX, zeros(size(zeroLineX)), '-', 'LineWidth', axisLineWidth, 'Color', axisLineColor);
    hold on

    for ri = 1:nROIs
        beta = F{ri}.estimates(termIdx) / 2;
        dotPro = beta;
        dotCon = -beta;

        bootDraws = F{ri}.coeffs(termIdx, :)';
        estDiff = F{ri}.estimates(termIdx);
        ci68 = prctile(bootDraws, [16 84]);
        ci95 = prctile(bootDraws, [2.5 97.5]);
        ci68_halfwidth = (ci68(2) - ci68(1)) / 2;

        statsRows(end+1,:) = {rois{ri}, conceptLabels{termIdx}, dotPro, dotCon, estDiff, ci68(1), ci68(2), ci95(1), ci95(2)}; %#ok<AGROW>

        x = ri;
        % Pro always filled, con always unfilled -- matches
        % plot2_experimentalCond.m's pairwise dots exactly (fixed
        % convention, not swapped by which of pro/con is bigger).
        scatter(x, dotPro, meanDotSize, 'MarkerFaceColor', proFaceColor, 'MarkerEdgeColor', proEdgeColor, 'LineWidth', proLineWidth);
        hold on
        scatter(x, dotCon, meanDotSize, 'MarkerFaceColor', conFaceColor, 'MarkerEdgeColor', conEdgeColor, 'LineWidth', proLineWidth);
        hold on
        errorbar(x, dotPro, ci68_halfwidth, 'Color', proEdgeColor, 'LineWidth', errorbarLineWidth, 'CapSize', 0);
        hold on
        errorbar(x, dotCon, ci68_halfwidth, 'Color', conEdgeColor, 'LineWidth', errorbarLineWidth, 'CapSize', 0);
        hold on

        sig95 = ci95(1) > 0 || ci95(2) < 0;
        sig68 = ci68(1) > 0 || ci68(2) < 0;
        if sig95
            sigStr = '**';
        elseif sig68
            sigStr = '*';
        else
            sigStr = '';
        end
        if ~isempty(sigStr)
            text(x, 0.35, sigStr, 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'FontSize', 20, 'Color', axisLineColor);
        end
    end

    ylim([-0.4 0.4])
    xlim(xlimRange)
    set(gca, 'XTick', 1:nROIs, 'XTickLabel', rois, 'XTickLabelRotation', 25, 'FontSize', 16);
    ax = gca;
    ax.LineWidth = axisLineWidth;
    ax.XColor = axisLineColor;
    ax.YColor = axisLineColor;
    box off
    ylabel({'BOLD response amplitude', '(% signal change)'}, 'FontSize', 9);
    title(conceptLabels{termIdx}, 'FontSize', 10, 'Interpreter', 'none');

    % Physical cm sizing, matching plot2_experimentalCond.m's square-PDF
    % export technique (Units=centimeters, explicit PaperSize/
    % PaperPosition) instead of a pixel Position + '-bestfit'.
    axisWidth_cm = nROIs * xScale_cm_per_unit;
    figWidth_cm = axisWidth_cm + extraLeftForYLabel_cm + padding_cm;
    figHeight_cm = pairwisePlotHeight_cm + 2*padding_cm; % pairwisePlotHeight_cm above is already doubled, so this figure canvas grows to match the taller axis box directly (no separate extra-blank-margin doubling)

    gcf_edit = gcf;
    gcf_edit.Units = 'centimeters';
    gcf_edit.Position(3:4) = [figWidth_cm, figHeight_cm];
    gcf_edit.PaperUnits = 'centimeters';
    gcf_edit.PaperPositionMode = 'manual';
    gcf_edit.PaperSize = [figWidth_cm, figHeight_cm];
    gcf_edit.PaperPosition = [0, 0, figWidth_cm, figHeight_cm];
    ax.Units = 'centimeters';
    ax.Position = [extraLeftForYLabel_cm, padding_cm, axisWidth_cm, pairwisePlotHeight_cm];

    % nativeOnly tag always present (even when off) so current output is
    % never mistaken for the future nativeOnly-mode output once that
    % exists -- see header.
    if nativeOnly
        nativeOnlyTag = 'nativeonlyON';
    else
        nativeOnlyTag = 'nativeonlyOFF';
    end
    filename = fullfile(figureDir, sprintf('AsymmetryAcrossROIs_%s_%s_%s_%s', comparisonName, projectName, colorKeys{termIdx}, nativeOnlyTag));
    set(gcf_edit, 'Renderer', 'painters');
    print(gcf_edit, filename, '-dpdf', '-painters');
    close all;
end

if nativeOnly
    nativeOnlyTag = 'nativeonlyON';
else
    nativeOnlyTag = 'nativeonlyOFF';
end
statsTable = cell2table(statsRows, 'VariableNames', ...
    {'roi','asymmetry','pro_mean','con_mean','diff_estimate','ci68_lower','ci68_upper','ci95_lower','ci95_upper'});
writetable(statsTable, fullfile(figureDir, sprintf('AsymmetryAcrossROIs_%s_%s_%s_stats.csv', comparisonName, projectName, nativeOnlyTag)));

end

function zeroLineX = buildDashedZeroLine(xlimRange, dashLen, gapLen)
% Manually-constructed dash pattern (NaN-separated segments), same
% technique plot2_experimentalCond.m's subjectwiseDiff zero line uses.
    unitLen = dashLen + gapLen;
    zeroLineX = [];
    for u = 0:ceil(diff(xlimRange)/unitLen)-1
        segStart = xlimRange(1) + u*unitLen;
        if segStart >= xlimRange(2), break; end
        segEnd = min(segStart+dashLen, xlimRange(2));
        zeroLineX = [zeroLineX, segStart, segEnd, NaN]; %#ok<AGROW>
    end
end
