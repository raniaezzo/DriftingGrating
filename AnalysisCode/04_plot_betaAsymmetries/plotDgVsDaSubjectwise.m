function plotDgVsDaSubjectwise(projectSettings)
% PLOTDGVSDASUBJECTWISE  dg-vs-da context comparison, in the same visual
% family as plot2_experimentalCond.m's 'subjectwiseDiff' mode (observer/
% group x-axis, jittered per-subject dots, dashed zero line, solid
% white-outlined dots) rather than plotDgVsDaComparison.m's bar/boxchart
% style. Restricted to the 7 matched observers, sourced from the same
% fitAsymmetryRegression_dgVsDa.m cached fit as plotDgVsDaComparison.m.
%
% One figure per cortical area, 4 columns (one per asymmetry, same
% display order/labels as plotDgVsDaComparison.m so the two figures'
% columns line up), 2 rows:
%
%   Row 1 (top): each subject's own value in BOTH contexts -- dg at
%   x~0.85, da at x~1.15 (both within the "observer" x=1 region), joined
%   by a short per-subject line so the dg-to-da change is visible
%   directly, same idea as plot2_experimentalCond.m's 'pairwise' mode's
%   per-subject connecting lines. "group" (x=2) mirrors this exactly at
%   the group level: mean(dg) at x~1.85, mean(da) at x~2.15, joined by
%   one line -- not a single collapsed difference point.
%
%   Row 2 (bottom): the per-subject (dg-minus-da) difference plotted
%   directly as a plain scatter (observer, x=1) with NO connecting lines
%   -- there is only one value per subject here, not a pair -- and the
%   group-level mean difference + 95% CI (paired bootstrap, from
%   F.diffBoot) at x=2, with significance asterisks (** = 95% CI excludes
%   0, * = only 68% CI does).
%
%   Row 1 and row 2 share fixed y-axis ranges across all 4 columns
%   ([-0.8, 0.4] and [-0.8, 0.1] respectively), rather than each column
%   auto-scaling to its own data.
%
% Natural (un-flipped) sign throughout -- unlike plotDgVsDaComparison.m,
% this does NOT multiply Horizontal-vs-Vertical/Cardinal-vs-Oblique by
% -1, so labels read literally (positive Horizontal-minus-Vertical means
% horizontal > vertical, not a fixed "expected direction").
%
%   plotDgVsDaSubjectwise(projectSettings)
%
% projectSettings - must include .rois, .colors_data, .figureDir,
%                    .gainSummaryFile (used only to locate bidsDir, same
%                    convention as plotDgVsDaComparison.m).

comparisonName = projectSettings.comparisonName;
colors_data = projectSettings.colors_data;
rois = projectSettings.rois;
figureDir = projectSettings.figureDir;

[summaryTablesDir,~,~] = fileparts(projectSettings.gainSummaryFile);
[derivativesDir,~,~] = fileparts(summaryTablesDir);
[bidsDir,~,~] = fileparts(derivativesDir);

% Same display order/labels as plotDgVsDaComparison.m (conceptIdx indexes
% dgContribConcept/daContribConcept/diffConcept/diffBoot's own 1..4
% concept order: Cardinal-Oblique, Polar Cardinal-Oblique,
% Horizontal-Vertical, Radial-Tangential), so the two figures' columns
% are directly comparable. No sign flip here (see header).
xOrderConceptIdx = [3, 1, 4, 2];
xLabels = {'Horizontal minus Vertical','Cardinal minus Oblique','Radial minus Tangential','Polar Cardinal minus Polar Oblique'};
dgColorKeysByConcept = {'mainCardinalVsMainOblique','derivedCardinalVsDerivedOblique','verticalVsHorizontal','radialVsTangential'};

nA = 4;

% Styling constants, matching plot2_experimentalCond.m's subjectwiseDiff
% conventions exactly (light grey axes/dashed line, tight NaN-segmented
% dashes, solid dots with a thin white outline sized to include the old
% full-stroke footprint) so this figure family reads consistently with
% the rest of the pipeline.
subjDiffAxisColor = [186, 188, 190] / 255;
subjDiffAxisLineWidth = 1;
subjectLineColor = [0.8196, 0.8275, 0.8314]; % per-subject dg-da connecting line, row 1 only
polarMarkerSize = 6 * 0.8;
meanDotSize = (pi/4) * polarMarkerSize^2;
xlimRange = [0.5 2.5];
dashLen = 0.12; gapLen = 0.06;

for ri = 1:length(rois)
    roiname = rois{ri};

    fitFile = fullfile(bidsDir,'derivatives','summaryTables','regressionResults','dgVsDa7',sprintf('%s.mat',roiname));
    if ~isfile(fitFile)
        error('plotDgVsDaSubjectwise:missingFit', ...
            'Cached fit not found for dgVsDa7 / %s at %s -- run fitAsymmetryRegression_dgVsDa() first.', ...
            roiname, fitFile);
    end
    F = load(fitFile);
    nSubj = F.nSubj;

    figure
    set(gcf, 'Position', [200 100 1300 700])

    for c = 1:nA
        conceptIdx = xOrderConceptIdx(c);
        colorKey = dgColorKeysByConcept{conceptIdx};
        styleInfo = colors_data.conditions.dg.(colorKey);
        proColor = styleInfo.color_pro';
        daColor = 0.5*proColor + 0.5*[1 1 1]; % lighter tint distinguishes da from dg within the paired display
        proLineWidth = styleInfo.pro_lineWidth;
        errorbarLineWidth = styleInfo.errorbar_lineWidth;

        dgVals = nSubj * F.dgContribConcept(:, conceptIdx);
        daVals = nSubj * F.daContribConcept(:, conceptIdx);
        diffVals = dgVals - daVals;
        groupDiff = F.diffConcept(conceptIdx);
        % F.diffBoot: each bootstrap draw resamples ONE set of subject
        % indices and applies it to BOTH dg's and da's data in the SAME
        % iteration (fitAsymmetryRegression_dgVsDa.m's own header), so
        % this is a genuine PAIRED bootstrap of the dg-minus-da
        % difference -- these percentile intervals describe the
        % difference itself, not two independently-bootstrapped context
        % means subtracted after the fact.
        bootDiff = F.diffBoot(conceptIdx, :);
        ci95 = prctile(bootDiff, [2.5 97.5]);
        ci95_halfwidth = (ci95(2) - ci95(1)) / 2;
        ci68 = prctile(bootDiff, [16 84]);
        sig95 = ci95(1) > 0 || ci95(2) < 0;
        sig68 = ci68(1) > 0 || ci68(2) < 0;
        if sig95
            sigStr = '**';
        elseif sig68
            sigStr = '*';
        else
            sigStr = '';
        end

        % Dot face grown to include the OLD full-stroke footprint, then a
        % thin (proLineWidth/4) white outline on top -- identical
        % construction to plot2_experimentalCond.m's subjectwiseDiff dots.
        circleSizeOld_ = meanDotSize*0.4;
        oldCircleDiameterPts_ = sqrt(4*circleSizeOld_/pi) + proLineWidth;
        circleFaceSize_ = (pi/4) * oldCircleDiameterPts_^2;
        meanDiameterPts_ = sqrt(4*meanDotSize/pi) + proLineWidth;
        meanFaceSize_ = (pi/4) * meanDiameterPts_^2;
        whiteOutlineLineWidth_ = proLineWidth/4;

        zeroLineX = buildDashedZeroLine(xlimRange, dashLen, gapLen);

        %% Row 1: paired dg/da per-subject dots + lines; group = paired
        %% dg-mean/da-mean dots + one connecting line (mirrors "observer")
        row1Ylim = [-0.8, 0.4];
        subplot(2, nA, c);
        hold on
        plot(zeroLineX, zeros(size(zeroLineX)), '-', 'LineWidth', subjDiffAxisLineWidth, 'Color', subjDiffAxisColor);
        dgX = 0.85; daX = 1.15;
        for si = 1:nSubj
            plot([dgX daX], [dgVals(si) daVals(si)], '-', 'Color', subjectLineColor, 'LineWidth', 1);
        end
        % DG and DA both circles here (color alone marks context: proColor
        % vs the lighter daColor tint), same at both observer and group level.
        scatter(repmat(dgX,nSubj,1), dgVals, circleFaceSize_, ...
            'Marker', 'o', 'MarkerFaceColor', proColor, 'MarkerEdgeColor', 'w', 'LineWidth', whiteOutlineLineWidth_);
        scatter(repmat(daX,nSubj,1), daVals, circleFaceSize_, ...
            'Marker', 'o', 'MarkerFaceColor', daColor, 'MarkerEdgeColor', 'w', 'LineWidth', whiteOutlineLineWidth_);
        dgGroupX = 1.85; daGroupX = 2.15;
        meanDg = mean(dgVals); meanDa = mean(daVals);
        plot([dgGroupX daGroupX], [meanDg meanDa], '-', 'Color', proColor, 'LineWidth', 2);
        scatter(dgGroupX, meanDg, meanFaceSize_, ...
            'Marker', 'o', 'MarkerFaceColor', proColor, 'MarkerEdgeColor', 'w', 'LineWidth', whiteOutlineLineWidth_);
        scatter(daGroupX, meanDa, meanFaceSize_, ...
            'Marker', 'o', 'MarkerFaceColor', daColor, 'MarkerEdgeColor', 'w', 'LineWidth', whiteOutlineLineWidth_);
        xlim(xlimRange)
        ylim(row1Ylim)
        set(gca, 'XTick', [1 2], 'XTickLabel', {'observer','group'}, 'XTickLabelRotation', 45, 'FontSize', 8);
        ax1 = gca;
        ax1.LineWidth = subjDiffAxisLineWidth;
        ax1.XColor = subjDiffAxisColor;
        ax1.YColor = subjDiffAxisColor;
        box off
        titleParts = strsplit(xLabels{c}, ' minus ');
        title({titleParts{1}, ['minus ' titleParts{2}]}, 'FontSize', 8, 'Interpreter', 'none');
        if c == 1, ylabel('dg (dark) vs da (light)', 'FontSize', 9); end

        %% Row 2: dg-minus-da delta, plain scatter, no lines, with the
        %% paired-bootstrap CI and significance asterisks
        row2Ylim = [-0.8, 0.1];
        subplot(2, nA, nA + c);
        hold on
        plot(zeroLineX, zeros(size(zeroLineX)), '-', 'LineWidth', subjDiffAxisLineWidth, 'Color', subjDiffAxisColor);
        jitterWidth = 0.45;
        xJitter = 1 + (rand(nSubj,1) - 0.5) * jitterWidth;
        scatter(xJitter, diffVals, circleFaceSize_, ...
            'Marker', '^', 'MarkerFaceColor', proColor, 'MarkerEdgeColor', 'w', 'LineWidth', whiteOutlineLineWidth_);
        errorbar(2, groupDiff, ci95_halfwidth, 'Color', proColor, 'LineWidth', errorbarLineWidth, 'CapSize', 0);
        scatter(2, groupDiff, meanFaceSize_, ...
            'Marker', '^', 'MarkerFaceColor', proColor, 'MarkerEdgeColor', 'w', 'LineWidth', whiteOutlineLineWidth_);
        if ~isempty(sigStr)
            text(1.5, row2Ylim(2) - 0.05, sigStr, 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'FontSize', 10, 'Color', [0 0 0]);
        end
        xlim(xlimRange)
        ylim(row2Ylim)
        set(gca, 'XTick', [1 2], 'XTickLabel', {'observer','group'}, 'XTickLabelRotation', 45, 'FontSize', 8);
        ax2 = gca;
        ax2.LineWidth = subjDiffAxisLineWidth;
        ax2.XColor = subjDiffAxisColor;
        ax2.YColor = subjDiffAxisColor;
        box off
        if c == 1, ylabel('dg - da', 'FontSize', 9); end
    end

    sgtitle(sprintf('%s: dg vs da subjectwise (n=%d matched), %s', roiname, nSubj, strrep(comparisonName,'_','-')), ...
        'Interpreter', 'none', 'FontSize', 11)

    print(gcf, fullfile(figureDir, sprintf('DgVsDaSubjectwise_%s_%s', comparisonName, roiname)), '-dpdf', '-bestfit');
    close all;
end

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
