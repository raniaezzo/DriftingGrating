function plot_runTimeseries(runTimeseriesOutput, projectSettings)
% plot_runTimeseries - plot observed vs. predicted %signal-change time
% series for each run, with light background shading marking motion and
% static stimulus periods for context, small oriented tick marks showing
% grating orientation for the stationary (static) trials, and the R^2
% between the plotted observed and predicted traces printed in the title.
%
% One figure per run, saved as PDF + FIG (Run<N>_ObservedandPredicted),
% in projectSettings.runTimeseriesSavePath. Also saves a .mat per run with
% the underlying numeric traces (combined AND per-polar-angle-bin) so later
% cross-subject or by-polar-angle analyses don't require rerunning the slow
% GLM reconstruction / raw data loading steps.

subj = projectSettings.subject;
roiName = projectSettings.roiName;
stimdur_s = runTimeseriesOutput.stimdur_s;
tr_s = runTimeseriesOutput.tr_s;
meanR2 = runTimeseriesOutput.meanR2;
nIncludedVerts = runTimeseriesOutput.nIncludedVerts;

dataColor = [0 0 0];            % observed: black
modelColor = [0.85 0.1 0.1];    % predicted: red
shadeColors = struct('motion', [0.75 0.85 1], 'static', [1 0.85 0.7]);
shadeNames = {'motion', 'static'};   % blank omitted -- blends with white background

allCanonicalAngles = [0, 45, 90, 135, 180, 225, 270, 315];
if isfield(projectSettings, 'polarAnglesToInclude') && ~isempty(projectSettings.polarAnglesToInclude) ...
        && ~isequal(sort(projectSettings.polarAnglesToInclude), allCanonicalAngles)
    paLabel = sprintf(', pa %s%s\\pm%g%s', ...
        strjoin(arrayfun(@(a) num2str(a), sort(projectSettings.polarAnglesToInclude), 'UniformOutput', false), '/'), ...
        char(176), projectSettings.polarAngleBinWidth/2, char(176));
else
    paLabel = '';
end

infoLabel = sprintf('%s, ecc %g-%g%s, pRF R^2 \\geq %g%s, n=%d vertices', ...
    roiName, projectSettings.minECC, projectSettings.maxECC, char(176), ...
    projectSettings.minVAREXP, paLabel, nIncludedVerts);

nRuns = numel(runTimeseriesOutput.data);

if ~isfolder(projectSettings.runTimeseriesSavePath)
    mkdir(projectSettings.runTimeseriesSavePath)
end

for ri = 1:nRuns

    dataTrace = runTimeseriesOutput.data{ri};
    modelTrace = runTimeseriesOutput.model{ri};
    dataTrace_byPA = runTimeseriesOutput.data_byPA{ri};   % 1x4 cell, one per polarAngleConditions bin
    modelTrace_byPA = runTimeseriesOutput.model_byPA{ri};
    polarAngleConditions = runTimeseriesOutput.polarAngleConditions;
    onsets = runTimeseriesOutput.onsetsByType{ri}; % {motion, static, blank}, in TRs
    staticOnsets = runTimeseriesOutput.staticOnsetsByOrientation{ri}; % struct .TR, .deg
    nTRs = numel(dataTrace);
    timeAxis = ((1:nTRs) - 1) * tr_s;

    % run-level R^2: predicted vs. observed, from the exact traces plotted
    valid = isfinite(dataTrace) & isfinite(modelTrace);
    ssRes = sum((dataTrace(valid) - modelTrace(valid)).^2);
    ssTot = sum((dataTrace(valid) - mean(dataTrace(valid))).^2);
    runR2 = (1 - ssRes / ssTot) * 100;

    figure
    hold on

    % background shading for stimulus periods (drawn first, behind traces)
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

    plot(timeAxis, dataTrace, '-', 'LineWidth', 1.5, 'Color', dataColor);
    plot(timeAxis, modelTrace, '-', 'LineWidth', 1.5, 'Color', modelColor);

    dataRange = [dataTrace, modelTrace];
    dataRange = dataRange(isfinite(dataRange));
    if ~isempty(dataRange)
        pad = 0.1 * (max(dataRange) - min(dataRange) + eps);
        yl = [min(dataRange) - pad, max(dataRange) + pad];
    else
        yl = ylim;
    end
    ylim(yl);
    xlim([timeAxis(1), timeAxis(end)]);

    title(sprintf('%s %s run %d  (run R^2 = %.1f%%, GLMsingle R^2 = %.1f%%)', subj, roiName, ri, runR2, meanR2), 'Interpreter', 'none');
    subtitle(infoLabel, 'Interpreter', 'tex', 'FontSize', 11);
    xlabel('time (s)');
    ylabel('% signal change');
    ax = gca;
    ax.FontSize = 20;
    box on;

    % figure size: ~10:3 length-to-height ratio
    f1 = gcf;
    f1.Position = [72 712 1400 420];
    % avoid "figure too large for the page" warning/clipping on PDF export
    set(f1, 'PaperPositionMode', 'auto');
    set(f1, 'PaperUnits', 'points', 'PaperSize', f1.Position(3:4));

    % small oriented tick marks for the static (stationary grating) trials,
    % centered on each stimulus period, near the top of the plot.
    % Orientation convention: 0=horizontal, 45=right oblique (/), 90=vertical,
    % 135=left oblique (\), i.e. degrees measured counter-clockwise from horizontal.
    tickY = yl(2) - 0.08 * (yl(2) - yl(1));
    tickPx = 18;
    for si = 1:numel(staticOnsets.TR)
        onsetT = (staticOnsets.TR(si) - 1) * tr_s;
        centerT = onsetT + stimdur_s / 2;
        drawOrientationTick(ax, centerT, tickY, staticOnsets.deg(si), tickPx);
    end

    filename = sprintf('Run%d_%s_%s_ObservedandPredicted', ri, subj, projectSettings.projectName);
    saveOut = fullfile(projectSettings.runTimeseriesSavePath, filename);
    saveas(gcf, sprintf('%s.pdf', saveOut));
    savefig(gcf, sprintf('%s.fig', saveOut));

    subj_s = subj; roiName_s = roiName; ri_s = ri; %#ok<NASGU> (naming clarity in the saved file)
    save(sprintf('%s.mat', saveOut), 'dataTrace', 'modelTrace', 'dataTrace_byPA', 'modelTrace_byPA', ...
        'polarAngleConditions', 'onsets', 'staticOnsets', 'runR2', 'meanR2', 'nIncludedVerts', 'projectSettings', ...
        'subj_s', 'roiName_s', 'ri_s');

end

end


function drawOrientationTick(ax, xCenter, yPos, angleDeg, pixelLength)
% Draws a short line centered at (xCenter, yPos), using the axes' current
% pixel-per-data-unit scale so the line's visual angle on screen is correct
% (e.g. a 45-degree oblique looks like 45 degrees) regardless of the axes'
% data aspect ratio. angleDeg is measured counter-clockwise from horizontal
% (standard convention): 0=horizontal, 90=vertical, 45=right oblique, 135=left oblique.

origUnits = ax.Units;
ax.Units = 'pixels';
axPos = ax.Position; % [x y width height] in pixels
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