function plot_groupAverageRun(groupData, outputInfo)
% plot_groupAverageRun - plot the cross-subject group-mean observed time
% series for one run, together with a freshly-refit GLM prediction (the
% per-vertex, per-subject GLMsingle HRF fits can't be averaged directly,
% so a new canonical-HRF GLM is fit to the group-mean trace itself -- see
% refitGroupMeanTrace.m). Same visual style as plot_runTimeseries.m.
%
% Saves Run<N>_ObservedandPredicted.{pdf,fig,mat} in outputInfo.savePath.

projectName = outputInfo.projectName;
roiName = outputInfo.roiName;
ri = outputInfo.runIdx;

meanObserved = groupData.meanObserved;
predicted = groupData.predicted;
R2 = groupData.R2;
onsets = groupData.onsets; % {motion, static, blank}, in TRs
staticOnsets = groupData.staticOnsetsByOrientation; % struct .TR, .deg
stimdur_s = groupData.stimdur_s;
tr_s = groupData.tr_s;
includedSubjects = groupData.includedSubjects;
excludedSubjects = groupData.excludedSubjects;

dataColor = [0 0 0];
modelColor = [0.85 0.1 0.1];
shadeColors = struct('motion', [0.75 0.85 1], 'static', [1 0.85 0.7]);
shadeNames = {'motion', 'static'};

nTRs = numel(meanObserved);
timeAxis = ((1:nTRs) - 1) * tr_s;

if isempty(excludedSubjects)
    exclLabel = 'none excluded';
else
    exclLabel = sprintf('excluded: %s', strjoin(excludedSubjects, ', '));
end

allCanonicalAngles = [0, 45, 90, 135, 180, 225, 270, 315];
if isfield(outputInfo, 'polarAnglesToInclude') && ~isempty(outputInfo.polarAnglesToInclude) ...
        && ~isequal(sort(outputInfo.polarAnglesToInclude), allCanonicalAngles)
    paLabel = sprintf(', pa %s%s\\pm%g%s', ...
        strjoin(arrayfun(@(a) num2str(a), sort(outputInfo.polarAnglesToInclude), 'UniformOutput', false), '/'), ...
        char(176), outputInfo.polarAngleBinWidth/2, char(176));
else
    paLabel = '';
end

infoLabel = sprintf('%s, group average, n=%d subjects (%s), vertex criteria: ecc 4-8%s, pRF R^2 \\geq 0.1%s', ...
    roiName, numel(includedSubjects), exclLabel, char(176), paLabel);

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

plot(timeAxis, meanObserved, '-', 'LineWidth', 1.5, 'Color', dataColor);
plot(timeAxis, predicted, '-', 'LineWidth', 1.5, 'Color', modelColor);

dataRange = [meanObserved(:); predicted(:)];
dataRange = dataRange(isfinite(dataRange));
if ~isempty(dataRange)
    pad = 0.1 * (max(dataRange) - min(dataRange) + eps);
    yl = [min(dataRange) - pad, max(dataRange) + pad];
else
    yl = ylim;
end
ylim(yl);
xlim([timeAxis(1), timeAxis(end)]);

title(sprintf('%s group average %s run %d  (R^2 = %.1f%%)', projectName, roiName, ri, R2), 'Interpreter', 'none');
subtitle(infoLabel, 'Interpreter', 'tex', 'FontSize', 11);
xlabel('time (s)');
ylabel('% signal change');
ax = gca;
ax.FontSize = 20;
box on;

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

filename = sprintf('Run%d_ObservedandPredicted', ri);
saveOut = fullfile(outputInfo.savePath, filename);
saveas(gcf, sprintf('%s.pdf', saveOut));
savefig(gcf, sprintf('%s.fig', saveOut));

save(sprintf('%s.mat', saveOut), 'meanObserved', 'predicted', 'R2', 'onsets', 'staticOnsets', ...
    'includedSubjects', 'excludedSubjects', 'projectName', 'roiName', 'ri', 'stimdur_s', 'tr_s');

end


function drawOrientationTick(ax, xCenter, yPos, angleDeg, pixelLength)
% Same convention as plot_runTimeseries.m: 0=horizontal, 90=vertical,
% 45=right oblique, 135=left oblique, drawn in pixel space so the visual
% angle on screen is correct regardless of the axes' data aspect ratio.

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
