function plot_withinAcrossComparison(compData, outputInfo)
% plot_withinAcrossComparison - overlay a group-average run's observed
% trace with (a) its own "within" predicted trace and (b) a DIFFERENT
% run's "across" predicted trace, to visualize/quantify how much of the
% model's explanatory power is specific to that run vs. driven by the
% block-design structure shared across runs.
%
% <compData> fields:
%   meanObserved     - 1 x nTRs, run <runWithin>'s group-mean observed trace
%   predictedWithin  - 1 x nTRs, run <runWithin>'s own predicted trace
%   predictedAcross  - 1 x nTRs, run <runAcross>'s predicted trace, compared
%                       against run <runWithin>'s observed data
%   R2within         - percent R^2 of meanObserved vs. predictedWithin
%   R2across         - percent R^2 of meanObserved vs. predictedAcross
%   onsets, staticOnsetsByOrientation, stimdur_s, tr_s - from run <runWithin>
%   runWithin, runAcross - the two run indices being compared
%
% Saves Run<W>within_Run<A>across_Comparison.{pdf,fig,mat} in
% outputInfo.savePath.

projectName = outputInfo.projectName;
roiName = outputInfo.roiName;
runWithin = compData.runWithin;
runAcross = compData.runAcross;

meanObserved = compData.meanObserved;
predictedWithin = compData.predictedWithin;
predictedAcross = compData.predictedAcross;
R2within = compData.R2within;
R2across = compData.R2across;
onsets = compData.onsets;
staticOnsets = compData.staticOnsetsByOrientation;
stimdur_s = compData.stimdur_s;
tr_s = compData.tr_s;

dataColor = [0 0 0];
withinColor = [0.85 0.1 0.1];
acrossColor = [0.1 0.4 0.85];
shadeColors = struct('motion', [0.75 0.85 1], 'static', [1 0.85 0.7]);
shadeNames = {'motion', 'static'};

nTRs = numel(meanObserved);
timeAxis = ((1:nTRs) - 1) * tr_s;

if ~isfolder(outputInfo.savePath)
    mkdir(outputInfo.savePath)
end

figure
hold on

for k = 1:numel(shadeNames)
    theseOnsets = onsets{k};
    color = shadeColors.(shadeNames{k});
    for oi = 1:numel(theseOnsets)
        onsetT = (theseOnsets(oi) - 1) * tr_s;
        offsetT = min(onsetT + stimdur_s, timeAxis(end));
        patch([onsetT onsetT offsetT offsetT], ...
            [-1000 1000 1000 -1000], color, ...
            'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
    end
end

hData = plot(timeAxis, meanObserved, '-', 'LineWidth', 1.5, 'Color', dataColor, ...
    'DisplayName', sprintf('run %d observed', runWithin));
hWithin = plot(timeAxis, predictedWithin, '-', 'LineWidth', 1.5, 'Color', withinColor, ...
    'DisplayName', sprintf('run %d predicted (within)', runWithin));
hAcross = plot(timeAxis, predictedAcross, '--', 'LineWidth', 1.5, 'Color', acrossColor, ...
    'DisplayName', sprintf('run %d predicted (across)', runAcross));

if isfield(compData, 'ylimOverride') && ~isempty(compData.ylimOverride)
    % shared y-axis limits, e.g. so a within/across pair and its reverse
    % (run A vs B, then run B vs A) are directly visually comparable --
    % see runWithinAcrossComparisonPaired.m
    yl = compData.ylimOverride;
else
    dataRange = [meanObserved(:); predictedWithin(:); predictedAcross(:)];
    dataRange = dataRange(isfinite(dataRange));
    if ~isempty(dataRange)
        pad = 0.1 * (max(dataRange) - min(dataRange) + eps);
        yl = [min(dataRange) - pad, max(dataRange) + pad];
    else
        yl = ylim;
    end
end
ylim(yl);
xlim([timeAxis(1), timeAxis(end)]);

title(sprintf('%s group average %s: run %d observed  (within R^2 = %.1f%%, across R^2 [run %d model] = %.1f%%)', ...
    projectName, roiName, runWithin, R2within, runAcross, R2across), 'Interpreter', 'none');
xlabel('time (s)');
ylabel('% signal change');
ax = gca;
ax.FontSize = 20;
box on;
legend([hData, hWithin, hAcross], 'Location', 'best');

f1 = gcf;
f1.Position = [72 712 1400 420]; % ~10:3 length-to-height ratio
set(f1, 'PaperPositionMode', 'auto');
set(f1, 'PaperUnits', 'points', 'PaperSize', f1.Position(3:4));

tickY = yl(2) - 0.08 * (yl(2) - yl(1));
tickPx = 18;
for si = 1:numel(staticOnsets.TR)
    onsetT = (staticOnsets.TR(si) - 1) * tr_s;
    centerT = onsetT + stimdur_s / 2;
    drawOrientationTick(ax, centerT, tickY, staticOnsets.deg(si), tickPx);
end

filename = sprintf('Run%dwithin_Run%dacross_Comparison', runWithin, runAcross);
saveOut = fullfile(outputInfo.savePath, filename);
saveas(gcf, sprintf('%s.pdf', saveOut));
savefig(gcf, sprintf('%s.fig', saveOut));

save(sprintf('%s.mat', saveOut), 'meanObserved', 'predictedWithin', 'predictedAcross', ...
    'R2within', 'R2across', 'runWithin', 'runAcross', 'onsets', 'staticOnsets', ...
    'projectName', 'roiName', 'stimdur_s', 'tr_s');

end


function drawOrientationTick(ax, xCenter, yPos, angleDeg, pixelLength)
% Same convention as plot_runTimeseries.m / plot_groupAverageRun.m:
% 0=horizontal, 90=vertical, 45=right oblique, 135=left oblique, drawn in
% pixel space so the visual angle on screen is correct regardless of the
% axes' data aspect ratio.

origUnits = ax.Units;
ax.Units = 'pixels';
axPos = ax.Position;
ax.Units = origUnits;

xl = xlim(ax);
yl = ylim(ax);
pxPerXData = axPos(3) / diff(xl);
pxPerYData = axPos(4) / diff(yl);

theta = deg2rad(angleDeg);
dxPx = pixelLength * cos(theta);
dyPx = pixelLength * sin(theta);

dxData = dxPx / pxPerXData;
dyData = dyPx / pxPerYData;

line(ax, [xCenter - dxData/2, xCenter + dxData/2], ...
         [yPos - dyData/2, yPos + dyData/2], ...
         'Color', [0.15 0.15 0.15], 'LineWidth', 2, 'HandleVisibility', 'off');

end
