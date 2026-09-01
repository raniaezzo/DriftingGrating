function plotROISummary(projectSettings, varargin)
% PLOTROISUMMARY  Successor to lme1_fit.m's per-cortical-area "master
% figure" visualization: one figure per cortical area, all 4 asymmetries
% along the x-axis, pro vs con each shown as a point + error bar. Sourced
% from fitAsymmetryRegression.m's cached joint regression fit instead of
% the LME. Fully independent of lme1_fit.m -- does not read or write
% anything under LME_results/; lme1_fit.m and its output remain untouched
% and independently runnable.
%
%   plotROISummary(projectSettings)
%   plotROISummary(projectSettings, 'nativeOnly', true)
%
% nativeOnly (optional name-value, default false): when true, restricts
% to the 2 asymmetries that are "native" to this project (dg: Horizontal
% vs Vertical, Cardinal vs Oblique; da: Radial vs Tangential, Polar
% Cardinal vs Polar Oblique -- see REF_ORIENTATION.md sec 4), sourced from
% fitAsymmetryRegressionNative.m's cache (vertex mean within ecc[0.5,12],
% no polar-angle binning, no varexp filter -- a genuinely different
% analysis, not just a display filter) instead of
% fitAsymmetryRegression.m's.
%
% projectSettings - same struct plot1/plot2_experimentalCond.m use (must
%                    include .rois, .roi_idx, .comparisonName,
%                    .projectName, .colors_data, .figureDir,
%                    .gainSummaryFile -- the last only to locate bidsDir,
%                    same convention as plot2_experimentalCond.m).
%
% Pro/con points are MODEL-DERIVED, zero-centered: pro = +beta,
% con = -beta (beta = estimates(termIdx)/2, since estimates is already
% the beta*2 pro-minus-con scale) -- the direct readout of the +-1/0
% design coding, with no intercept added back in. This matches what
% lme1_fit.m's master figure ACTUALLY plots: it locally computes
% Gintercept+beta/Gintercept-beta, but then (since meanRelative=1 is its
% configured setting) subtracts baselineSub=Gintercept back out before
% storing/plotting (see its mean_pro(ai) = y1(x)-baselineSub, line ~543),
% so the intercept never actually appears in its figure -- only beta and
% -beta do. This is a DELIBERATE difference from
% plot1_experimentalCond.m/plot2_experimentalCond.m, whose dots are
% data-derived (absolute BOLD levels) rather than model-derived
% (asymmetry-only, zero-centered) -- this script's whole purpose is to
% reproduce lme1_fit.m's own design (all 4 asymmetries sharing a common
% center at exactly 0, representing the asymmetry effect alone, decoupled
% from the cortical area's overall response magnitude), not plot2's
% "dots stay data-derived" convention.
%
% Pro/con dots, error bars, and axes match plot2_experimentalCond.m's
% 'pairwise' convention exactly (same as plotAsymmetryAcrossROIs.m's
% equivalent rewrite): pro always FILLED, con always UNFILLED (white
% face, colored edge) -- a fixed convention, not swapped by which of
% pro/con is bigger; meanDotSize/proLineWidth/errorbarLineWidth read from
% the same COLORS.json styleInfo fields; axisLineColor/axisLineWidth
% instead of plain black; a NaN-segmented dashed zero line (not yline());
% and the SAME cm-per-x-unit physical density plot2 uses for its own
% pro/con spacing, applied across this figure's 4 asymmetry x-positions.
% Unlike plotAsymmetryAcrossROIs.m's rewrite, the plot/figure height here
% is NOT doubled -- kept at plot2's own per-subplot height.
%
% Error bars: 68% CI of the cached fit's bootstrapped difference
% (coeffs), FULL half-width applied symmetrically to both pro and con --
% the SAME convention plot1_experimentalCond.m/plot2_experimentalCond.m
% already use (confirmed deliberately, not lme1_fit.m's own convention of
% mirroring the raw, undoubled coefficient's own CI onto each side, which
% corresponds to a more liberal ~68%-of-the-difference significance
% boundary rather than the ~95%-ish one this convention approximates).

p = inputParser;
p.addParameter('nativeOnly', false, @islogical);
p.parse(varargin{:});
nativeOnly = p.Results.nativeOnly;

projectName = projectSettings.projectName;
comparisonName = projectSettings.comparisonName;
colors_data = projectSettings.colors_data;
rois = projectSettings.rois;
figureDir = projectSettings.figureDir;

% Locate fitAsymmetryRegression.m's cached output (same convention as
% plot2_experimentalCond.m).
[summaryTablesDir,~,~] = fileparts(projectSettings.gainSummaryFile);
[derivativesDir,~,~] = fileparts(summaryTablesDir);
[bidsDir,~,~] = fileparts(derivativesDir);

% x-axis order and display labels: matches lme1_fit.m's plotOrder/
% asymLabel exactly (desired order: vertical-vs-horizontal,
% cardinal-vs-oblique, radial-vs-tangential, polar-cardinal-vs-oblique,
% relabeled per project since which raw term is which concept swaps
% between dg/da -- see fitAsymmetryRegression_dgVsDa.m's
% dgTermForConcept/daTermForConcept for the same mapping). termIdx here
% indexes termNames = {mainCardinal, derivedCardinal, mainSubset,
% derivedSubset}, the order fitAsymmetryRegression.m saves estimates/
% coeffs in.
if strcmp(projectName, 'dg')
    xOrderTermIdx = [3, 1, 4, 2];
    xLabels = {'Horizontal vs Vertical','Cardinal vs Oblique','Radial vs Tangential','Polar Cardinal vs Polar Oblique'};
    colorKeys = {'mainCardinalVsMainOblique','derivedCardinalVsDerivedOblique','verticalVsHorizontal','radialVsTangential'}; % indexed by termIdx (1..4), not x-position
elseif strcmp(projectName, 'da')
    xOrderTermIdx = [4, 2, 3, 1];
    xLabels = {'Horizontal vs Vertical','Cardinal vs Oblique','Radial vs Tangential','Polar Cardinal vs Polar Oblique'};
    colorKeys = {'mainCardinalVsMainOblique','derivedCardinalVsDerivedOblique','radialVsTangential','verticalVsHorizontal'};
else
    error('plotROISummary:project', 'projectName must be ''dg'' or ''da''.');
end

% nativeOnly restricts to termIdx {3,1} (mainSubset, mainCardinal) -- the
% two terms fitAsymmetryRegression.m's design never conditions on
% location (see that file's mainCardinal/mainSubset assignment) -- and
% keeps only their entries from the xLabels/xOrderTermIdx built above, in
% the same relative order.
if nativeOnly
    nativeTermOrder = [3, 1];
    [~, keepPos] = ismember(nativeTermOrder, xOrderTermIdx);
    xOrderTermIdx = xOrderTermIdx(keepPos);
    xLabels = xLabels(keepPos);
end
nA = numel(xOrderTermIdx);

% fitLabel: which regressionResults/<label>/ (or regressionResultsNative/
% for nativeOnly) subfolder to read -- defaults to projectName, but for
% dg 'matched' mode (7 subjects, same ones run in da) set
% projectSettings.fitLabel = 'dgMatched7' explicitly.
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

% Styling constants, matching plot2_experimentalCond.m's pairwise
% convention exactly (same source fields, same formula).
axisLineWidth = 1;
axisLineColor = [0.25 0.25 0.25];
polarMarkerSize = 6 * 0.8;
meanDotSize = (pi/4) * polarMarkerSize^2;
dashLen = 0.12; gapLen = 0.06; % same NaN-segmented dash construction as plot2's subjectwiseDiff zero line

% Physical axis sizing: reuse plot2_experimentalCond.m's own
% cm-per-x-unit density exactly, so asymmetry spacing here reads at the
% same visual scale as pro/con spacing in the pairwise plots, just
% extended across nA=4 positions instead of 2.
previousXlim = [0.5, 2.5];
previousPairwisePlotWidth_cm = 3 * 1.25 * 0.94;
xScale_cm_per_unit = previousPairwisePlotWidth_cm / diff(previousXlim);
pairwisePlotHeight_cm = 3.6 * 0.94; % plot2's own per-subplot height, NOT doubled (unlike plotAsymmetryAcrossROIs.m)
padding_cm = 3; % page margin on top/left/right -- widened from 1 alongside the 2026-08-31 tick-label/asterisk font doubling (the rotated 2-line ylabel's far end otherwise clips the page top)
% Extra bottom margin -- these are full multi-word asymmetry names
% ("Horizontal vs Vertical"), much longer than plotAsymmetryAcrossROIs.m's
% short ROI labels ("V1"), so the standard padding_cm clips them at 25
% degrees rotation. Same fixed-margin approach plot2_experimentalCond.m's
% subjectwiseDiff mode uses for its own (shorter) rotated labels.
% Margins below are ~2x their original values (7/4.5/7 vs the prior
% 3.5/2.2/3.5) to accommodate the tick-label/asterisk font size doubling
% (8->16 / 10->20) requested 2026-08-31 -- text roughly doubled in
% physical extent, so the surrounding margin needed to as well.
extraBottomForXTickLabels_cm = 7;
% Extra left margin -- the leftmost label ("Horizontal vs Vertical")
% extends left of its own tick at 25 degrees rotation, past the axis's
% own left edge, same reason plot2_experimentalCond.m adds
% extraLeftForXTickLabels_cm for its own rotated labels. Also needs to
% fit the (now larger) Y-tick labels themselves.
extraLeftForXTickLabels_cm = 4.5;
% Extra right margin -- the rightmost label ("Polar Cardinal vs Polar
% Oblique") is the longest of the four and clips against the standard
% padding_cm on that side too.
extraRightForXTickLabels_cm = 7;
xlimRange = [0.5, nA+0.5];
axisWidth_cm = nA * xScale_cm_per_unit;
figWidth_cm = axisWidth_cm + extraLeftForXTickLabels_cm + extraRightForXTickLabels_cm;
figHeight_cm = pairwisePlotHeight_cm + padding_cm + extraBottomForXTickLabels_cm;

statsRows = {};

for ri = 1:length(rois)
    roiname = rois{ri};

    fitFile = fullfile(bidsDir,'derivatives','summaryTables',fitSubdir,fitLabel,sprintf('%s.mat',roiname));
    if ~isfile(fitFile)
        if nativeOnly
            error('plotROISummary:missingFit', ...
                'Cached native fit not found for %s / %s at %s -- run fitAsymmetryRegressionNative(''%s'') first.', ...
                fitLabel, roiname, fitFile, projectName);
        else
            error('plotROISummary:missingFit', ...
                'Cached fit not found for %s / %s at %s -- run fitAsymmetryRegression(''%s'') first.', ...
                fitLabel, roiname, fitFile, projectName);
        end
    end
    F = load(fitFile);

    figure; hold on

    zeroLineX = buildDashedZeroLine(xlimRange, dashLen, gapLen);
    plot(zeroLineX, zeros(size(zeroLineX)), '-', 'LineWidth', axisLineWidth, 'Color', axisLineColor);
    hold on

    for ai = 1:nA
        termIdx = xOrderTermIdx(ai);
        beta = F.estimates(termIdx) / 2;
        dotPro = beta;
        dotCon = -beta;

        bootDraws = F.coeffs(termIdx, :)';
        estDiff = F.estimates(termIdx);
        ci68 = prctile(bootDraws, [16 84]);
        ci95 = prctile(bootDraws, [2.5 97.5]);
        ci68_halfwidth = (ci68(2) - ci68(1)) / 2;

        statsRows(end+1,:) = {roiname, xLabels{ai}, dotPro, dotCon, estDiff, ci68(1), ci68(2), ci95(1), ci95(2)}; %#ok<AGROW>

        styleInfo = colors_data.conditions.(projectName).(colorKeys{termIdx});
        proFaceColor = styleInfo.color_pro';
        conFaceColor = [1 1 1];
        proEdgeColor = styleInfo.color_pro';
        conEdgeColor = styleInfo.color_con'; % full-strength con color, no white blending
        proLineWidth = styleInfo.pro_lineWidth;
        errorbarLineWidth = styleInfo.errorbar_lineWidth / 2; % half width per feedback

        x = ai;
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

        % Significance asterisk at a fixed y=0.35 for this asymmetry: two
        % asterisks if the 95% CI of the difference excludes 0, one if
        % only the 68% CI does (95% always implies 68%, since it's the
        % wider, nested interval) -- same convention as
        % plot2_experimentalCond.m's asterisks.
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
    ax = gca;
    ax.LineWidth = axisLineWidth;
    ax.XColor = axisLineColor;
    ax.YColor = axisLineColor;
    ax.Clipping = 'off'; % rotated XTickLabel text otherwise clips against the axes' own box, regardless of how much figure margin surrounds it
    box off
    ylabel({'BOLD response amplitude', '(% signal change)'}, 'FontSize', 9);
    title(roiname, 'FontSize', 10, 'Interpreter', 'none');

    % Physical cm sizing, matching plot2_experimentalCond.m's square-PDF
    % export technique (Units=centimeters, explicit PaperSize/
    % PaperPosition) instead of a pixel Position + '-bestfit'.
    gcf_edit = gcf;
    gcf_edit.Units = 'centimeters';
    gcf_edit.Position(3:4) = [figWidth_cm, figHeight_cm];
    gcf_edit.PaperUnits = 'centimeters';
    gcf_edit.PaperPositionMode = 'manual';
    gcf_edit.PaperSize = [figWidth_cm, figHeight_cm];
    gcf_edit.PaperPosition = [0, 0, figWidth_cm, figHeight_cm];
    ax.Units = 'centimeters';
    ax.Position = [extraLeftForXTickLabels_cm, extraBottomForXTickLabels_cm, axisWidth_cm, pairwisePlotHeight_cm];

    % XTick/XTickLabel set AFTER the figure/axes are resized to their
    % final physical cm dimensions -- setting them earlier (against the
    % default small figure size) left the rotated label text laid out/
    % cached for that smaller canvas, so it stayed clipped even once the
    % figure was later enlarged.
    set(ax, 'XTick', 1:nA, 'XTickLabel', xLabels, 'XTickLabelRotation', 25, 'FontSize', 16);

    if nativeOnly
        nativeOnlyTag = 'nativeonlyON';
    else
        nativeOnlyTag = 'nativeonlyOFF';
    end
    filename = fullfile(figureDir, sprintf('ROIsummary_%s_%s_%s_%s', comparisonName, projectName, roiname, nativeOnlyTag));
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
writetable(statsTable, fullfile(figureDir, sprintf('ROIsummary_%s_%s_%s_stats.csv', comparisonName, projectName, nativeOnlyTag)));

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
