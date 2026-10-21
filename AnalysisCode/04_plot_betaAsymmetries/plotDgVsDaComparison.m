function plotDgVsDaComparison(projectSettings)
% PLOTDGVSDACOMPARISON  dg-vs-da comparison, restricted to the 7 matched
% observers (sharedSubjects), sourced from
% fitAsymmetryRegression_dgVsDa.m's cached fit. One figure per cortical
% area, 4 columns (one per asymmetry). Each column has 3 rows:
%   1) dg (cartesian experiment): per-subject bar, all 7 observers
%   2) da (polar experiment): per-subject bar, all 7 observers
%   3) main plot: group-level DG minus DA point + 68% CI error bar
%      (full-width convention, same as plotROISummary.m/
%      plotAsymmetryAcrossROIs.m), spanning all 4 asymmetries on one axis
%
% Sign convention: Horizontal-vs-Vertical and Cardinal-vs-Oblique are
% multiplied by -1 throughout (bars AND the main point) so that all 4
% asymmetries' "expected" direction is positive, matching
% Radial-vs-Tangential and Polar-Cardinal-vs-Polar-Oblique's native sign
% -- per-subject bars and the main point are FILLED when positive (in
% this now-common "expected" direction) and UNFILLED (white fill,
% colored edge) when negative, so a reader can see at a glance which
% observers/asymmetries go against the expected direction even after the
% sign flip.
%
% Per-subject bar heights use the cached fit's contribution decomposition
% (dgContribConcept/daContribConcept), rescaled by nSubj (=7) -- these
% reduce to each subject's own raw marginal difference for a
% complete-data subject, and correctly handle missing-data confounds
% otherwise (see fitAsymmetryRegression.m's note on subjectContributions
% for the full derivation of why the nSubj rescale is needed).
%
%   plotDgVsDaComparison(projectSettings)
%
% projectSettings - must include .rois, .colors_data, .figureDir,
%                    .gainSummaryFile (used only to locate bidsDir, same
%                    convention as the other new scripts in this file).

comparisonName = projectSettings.comparisonName;
colors_data = projectSettings.colors_data;
rois = projectSettings.rois;
figureDir = projectSettings.figureDir;

[summaryTablesDir,~,~] = fileparts(projectSettings.gainSummaryFile);
[derivativesDir,~,~] = fileparts(summaryTablesDir);
[bidsDir,~,~] = fileparts(derivativesDir);

% Display order/labels/sign-flip, matching plotROISummary.m's x-axis
% convention. conceptIdx indexes dgConcept/daConcept/diffConcept/
% dgContribConcept/daContribConcept/diffBoot's own 1..4 concept order
% (Cardinal-Oblique, Polar Cardinal-Oblique, Horizontal-Vertical,
% Radial-Tangential -- fitAsymmetryRegression_dgVsDa.m's conceptLabels).
xOrderConceptIdx = [3, 1, 4, 2];
% Labels stay in their original (un-flipped) word order -- these are
% understood as (Horizontal minus Vertical)*-1 and (Cardinal minus
% Oblique)*-1, not relabeled to match the flipped sign. Fill state (see
% localBarPlot and the main-plot loop below) always reflects the
% ORIGINAL, un-flipped pro>con relationship (pro=filled, matching this
% pipeline's convention everywhere else) -- NOT the sign of the flipped
% displayed value -- so for these two columns, positive-displayed
% (=con>pro) is unfilled and negative-displayed (=pro>con) is filled,
% the reverse of the two un-flipped columns.
% "minus" rather than "vs" (more precise: this is literally pro-minus-con),
% and the two flipped asymmetries get an explicit "(x-1)" marker so the
% sign flip is visible directly in the label, not just inferable from the
% fill-state convention.
xLabels = {'Horizontal minus Vertical (x-1)','Cardinal minus Oblique (x-1)','Radial minus Tangential','Polar Cardinal minus Polar Oblique'};
signFlip = [-1, -1, 1, 1]; % matches xLabels order (and xOrderConceptIdx)

% Color per concept: dg's own COLORS.json key for that concept, used for
% BOTH the dg and da rows (dg and da's own per-concept colors are
% expected to already agree; using one canonical source per column keeps
% this figure internally consistent regardless).
dgColorKeysByConcept = {'mainCardinalVsMainOblique','derivedCardinalVsDerivedOblique','verticalVsHorizontal','radialVsTangential'};

nA = 4;

for ri = 1:length(rois)
    roiname = rois{ri};

    fitFile = fullfile(bidsDir,'derivatives','summaryTables','regressionResults','dgVsDa7',sprintf('%s.mat',roiname));
    if ~isfile(fitFile)
        error('plotDgVsDaComparison:missingFit', ...
            'Cached fit not found for dgVsDa7 / %s at %s -- run fitAsymmetryRegression_dgVsDa() first.', ...
            roiname, fitFile);
    end
    F = load(fitFile);
    nSubj = F.nSubj;

    % Shared, symmetric-about-0 ylim for all 8 per-subject bar panels
    % (sign flips/reordering don't change the max abs value, so this can
    % be computed directly from the raw contribution matrices).
    barMax = nSubj * max(abs([F.dgContribConcept(:); F.daContribConcept(:)])) * 1.1;
    barYlim = [-barMax, barMax];

    figure
    set(gcf, 'Position', [200 100 1100 900])

    for c = 1:nA
        conceptIdx = xOrderConceptIdx(c);
        sf = signFlip(c);
        colorKey = dgColorKeysByConcept{conceptIdx};
        styleInfo = colors_data.conditions.dg.(colorKey);
        barColor = styleInfo.color_pro';

        % Row 1: dg per-subject bars
        subplot(3, nA, c);
        dgVals = sf * nSubj * F.dgContribConcept(:, conceptIdx);
        localBarPlot(dgVals, barColor, sf, barYlim);
        % 2-line title (split at " minus ") so adjacent narrow columns
        % don't visually collide into each other; the "(x-1)" marker, if
        % present, rides along on the second line.
        titleParts = strsplit(xLabels{c}, ' minus ');
        title({titleParts{1}, ['minus ' titleParts{2}]}, 'FontSize', 8, 'Interpreter', 'none');
        if c == 1, ylabel('dg', 'FontSize', 10); end

        % Row 2: da per-subject bars
        subplot(3, nA, nA + c);
        daVals = sf * nSubj * F.daContribConcept(:, conceptIdx);
        localBarPlot(daVals, barColor, sf, barYlim);
        if c == 1, ylabel('da', 'FontSize', 10); end
    end

    % Row 3: main DG-minus-DA plot, one wide axis spanning all 4 columns
    mainAx = subplot(3, nA, [2*nA+1, 2*nA+2, 2*nA+3, 2*nA+4]);
    hold on
    for c = 1:nA
        conceptIdx = xOrderConceptIdx(c);
        sf = signFlip(c);
        colorKey = dgColorKeysByConcept{conceptIdx};
        styleInfo = colors_data.conditions.dg.(colorKey);
        barColor = styleInfo.color_pro';

        diffVal = sf * F.diffConcept(conceptIdx);
        flippedBoot = sf * F.diffBoot(conceptIdx, :);
        ci68 = prctile(flippedBoot, [16 84]);
        ci95 = prctile(flippedBoot, [2.5 97.5]);
        ci68_halfwidth = (ci68(2) - ci68(1)) / 2;

        if diffVal >= 0
            faceColor = barColor; edgeColor = barColor;
        else
            faceColor = [1 1 1]; edgeColor = barColor;
        end
        boxchart(c, diffVal, 'BoxFaceColor', faceColor, 'BoxEdgeColor', edgeColor, 'LineWidth', 3, 'BoxWidth', 0.5);
        errorbar(c, diffVal, ci68_halfwidth, ci68_halfwidth, 'LineStyle','none', 'LineWidth', 2, 'Color', barColor);

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
            text(c, 0.35, sigStr, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 14, 'Color', [0 0 0]);
        end
    end
    yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2)
    ylim([-0.4 0.4])
    xlim([0.5 nA+0.5])
    box off
    set(gca,'linewidth',2, 'YColor',[0 0 0], 'XColor',[0 0 0]);
    set(gca,'FontName','Arial','FontSize',12);
    xticks(1:nA)
    xticklabels(xLabels)
    xtickangle(25)
    ylabel('DG - DA', 'FontSize', 12);

    % Increase the main plot's y:x display ratio (a little more square,
    % not fully 1:1) -- a physical Position resize only, xlim/ylim/tick
    % values above are untouched. Grows upward from its current bottom
    % edge into the whitespace below row 2.
    mainAx.Units = 'normalized';
    pos = mainAx.Position; % [left bottom width height]
    heightScale = 1.6;
    newHeight = pos(4) * heightScale;
    mainAx.Position = [pos(1), pos(2) - (newHeight - pos(4)), pos(3), newHeight];

    sgtitle({sprintf('%s: dg vs da (n=%d matched), %s', roiname, nSubj, strrep(comparisonName,'_','-')), ...
        'DG-DA above 0 = cartesian dominant; below 0 = polar dominant'}, 'Interpreter', 'none', 'FontSize', 11)

    print(gcf, fullfile(figureDir, sprintf('DgVsDaComparison_%s_%s', comparisonName, roiname)), '-dpdf', '-bestfit');
    close all;
end

end

function localBarPlot(vals, barColor, sf, barYlim)
% One bar per subject (vals is nSubj x 1, already sign-flipped by sf for
% display height/position). Fill state reflects the ORIGINAL (un-flipped)
% pro>con relationship, not the displayed sign: pro=filled everywhere
% else in this pipeline, so multiplying vals(si) back by sf (undoing the
% flip, since sf is +-1) recovers that original sign for the fill
% decision -- filled when the original value was positive (pro>con),
% unfilled when negative (con>pro), regardless of which way sf flipped
% the bar's displayed height.
    hold on
    for si = 1:numel(vals)
        if vals(si)*sf >= 0
            faceColor = barColor;
        else
            faceColor = [1 1 1];
        end
        bar(si, vals(si), 'FaceColor', faceColor, 'EdgeColor', barColor, 'LineWidth', 1.2, 'BarWidth', 0.7);
    end
    yline(0, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    xlim([0.5 numel(vals)+0.5]);
    ylim(barYlim);
    set(gca, 'XTick', 1:numel(vals), 'FontSize', 8);
    box off
end
