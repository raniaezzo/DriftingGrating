function plotEachDirLocRegression(projectName, roiname, varargin)
% PLOTEACHDIRLOCREGRESSION  Offshoot of lme2_ploteachDirLoc.m: the same
% 8-location "compass" grid of polar subplots (model prediction overlaid
% on empirical data, one subplot per location, spatially arranged so each
% subplot sits in its own direction from the center), but sourced from
% fitAsymmetryRegression.m's cached joint regression fit instead of
% LME_bold.mat/modeldata.mat, and with the empirical data gain- and
% precision-weighted across subjects (the original script's "meta-subject"
% averaging was an unweighted mean across all subjects' rows -- this
% version applies the same per-subject gain correction and per-(subject,
% cortical area) precision weighting used everywhere else in this
% pipeline). Fully independent of lme1_fit.m/lme2_ploteachDirLoc.m and
% their LME_results/ output.
%
%   plotEachDirLocRegression('dg', 'V1')
%   plotEachDirLocRegression('da', 'V1', 'precisionWeights', T)
%
% Model prediction (red line/dots): M * [grandInterceptFE; beta1..beta4],
% where beta_k = estimates(k)/2 (estimates is already the beta*2
% pro-minus-con scale) and grandInterceptFE is
% fitAsymmetryRegression.m's fixed-subject-intercept refit -- the analog
% of lme1_fit.m's fitlme-derived intercept, reconstructing the SAME kind
% of design-matrix-times-coefficients prediction the original script
% computed from the LME's estimates.
%
% Empirical data (black dots): gain-corrected then precision-weighted
% mean across subjects at each (direction, location) cell -- same
% weightedNanMean logic used in plot1/plot2/plotROISummary -- replacing
% the original script's unweighted varfun(@mean,...) across subjects.
%
% Orientation stimuli are axial (a 0-degree grating looks the same
% rotated 180 degrees), so each of the 4 real stimulus directions
% (0/90/45/135) is duplicated at direction+180 for display across all 8
% compass positions, exactly matching the original script's "if
% orientation_minus_baseline, duplicate direction+180" step -- valid
% because every one of the 4 asymmetry predictors is invariant to a
% direction+180 shift by construction (verified: {0,90,180,270} and
% {45,135,225,315} are each closed under +180, so mainCardinal etc. don't
% change), so the duplicated rows correctly reuse the same term values.

p = inputParser;
p.addParameter('bidsDir', '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/', @ischar);
p.addParameter('githubDir', '~/Documents/GitHub', @ischar);
p.addParameter('precisionWeights', [], @(x) isempty(x) || istable(x));
p.addParameter('figureDir', '', @ischar); % '' = default production location (strrep(bidsDir,'data_bids','figures')/<projectName>)
p.addParameter('dgSubjectMode', 'all', @(x) ismember(x, {'all','matched'})); % only used when projectName='dg'; 'matched' = the same 7 subjects also run in da (sub-0395 excluded), matching fitAsymmetryRegression.m's dgSubjectMode
p.parse(varargin{:});
opt = p.Results;

githubDir = opt.githubDir;
bidsDir = opt.bidsDir;

addpath(genpath(fullfile(githubDir, 'DriftingGrating', 'AnalysisCode')));
cd(fullfile(githubDir, 'DriftingGrating', 'AnalysisCode'));
setup_user('rania', bidsDir);

comparisonName = 'orientation_minus_baseline';
projectSettings = loadConfig(githubDir);
roi_idx = projectSettings.roi_idx;
rois = projectSettings.rois;
riSel = find(strcmp(rois, roiname));
if isempty(riSel)
    error('plotEachDirLocRegression:roi', 'roiname ''%s'' not found in projectSettings.rois.', roiname);
end
regionIndex = roi_idx{riSel};

if isempty(opt.figureDir)
    figureDirSuffix = '';
    if strcmp(projectName, 'dg')
        figureDirSuffix = ['_', opt.dgSubjectMode]; % matches plot_NeuralAsymmetries.m's figureDir convention (e.g. figures/dg_all or figures/dg_matched)
    end
    figureDir = [strrep(bidsDir, 'data_bids', 'figures'), projectName, figureDirSuffix];
else
    figureDir = opt.figureDir;
end
if ~isfolder(figureDir), mkdir(figureDir); end

contrasts_dict = projectSettings.contrasts_dict;
contrastnames = {contrasts_dict.contrasts.('dg_contrast_name')};
s0_idx = find(strcmp(contrastnames,'s0_v_b'));
s90_idx = find(strcmp(contrastnames,'s90_v_b'));
s45_idx = find(strcmp(contrastnames,'s45_v_b'));
s135_idx = find(strcmp(contrastnames,'s135_v_b'));
mdirvals_dg = [0, 90, 45, 135];
anglevals = [90, 45, 0, 315, 270, 225, 180, 135]; % index order matching meanBOLDpa's 2nd dim
maincardinalmDir = [0,90,180,270];
primaryMeridians = [90,0,270,180];

if strcmp(projectName, 'dg')
    if strcmp(opt.dgSubjectMode, 'all')
        subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
            'sub-0395', 'sub-0426', 'sub-0250', 'sub-0442', 'sub-wlsubj121', 'sub-wlsubj127', ...
            'sub-0397', 'sub-0427'};
        fitLabel = 'dg';
    else % 'matched'
        subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
            'sub-0426', 'sub-0250'};
        fitLabel = 'dgMatched7';
    end
elseif strcmp(projectName, 'da')
    subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
        'sub-0426', 'sub-0250'};
    fitLabel = 'da';
else
    error('plotEachDirLocRegression:project', 'projectName must be ''dg'' or ''da''.');
end

gainWeightsFile = fullfile(bidsDir, 'derivatives', 'summaryTables', 'gainSummaryByROI.mat');
Ggain = load(gainWeightsFile, 'gainTable');
gainWeights = retrieveObserverGainWeights2(subjects, roiname, Ggain.gainTable);
groupGain = exp(mean(log(gainWeights), 'omitnan')); % omitnan: see retrieveObserverGainWeights2.m
subjectScale = groupGain ./ gainWeights;
% Row vector, matching vals' orientation below -- MATLAB vector indexing
% preserves the INDEXED array's own orientation (not the index's), so a
% column here would silently broadcast into a matrix instead of an
% elementwise product against vals(valid) (same bug found and fixed
% earlier in plot1_experimentalCond.m).
precisionW = retrieveObserverPrecisionWeights(subjects, roiname, opt.precisionWeights);

fitFile = fullfile(bidsDir,'derivatives','summaryTables','regressionResults',fitLabel,sprintf('%s.mat',roiname));
if ~isfile(fitFile)
    error('plotEachDirLocRegression:missingFit', ...
        'Cached fit not found for %s / %s at %s -- run fitAsymmetryRegression(''%s'') first.', ...
        fitLabel, roiname, fitFile, projectName);
end
F = load(fitFile);
betaVec = F.estimates(:) / 2; % raw (undoubled) coefficients, same scale as M below expects
estimatesVec = [F.grandInterceptFE; betaVec]; % [intercept; mainCardinal; derivedCardinal; mainSubset; derivedSubset]

glmResultsfolder = fullfile(bidsDir, 'derivatives', strcat(projectName, 'GLM'), 'hRF_glmsingle');
S1 = load(fullfile(glmResultsfolder, 'meanBOLDpa'));
meanBOLDpa_full = S1.meanBOLDpa;
if strcmp(projectName,'dg')
    dg_full13 = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
        'sub-0395', 'sub-0426', 'sub-0250', 'sub-0442', 'sub-wlsubj121', 'sub-wlsubj127', ...
        'sub-0397', 'sub-0427'};
    subjIdx = cellfun(@(s) find(strcmp(dg_full13,s)), subjects);
else
    da_full8 = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
        'sub-0395', 'sub-0426', 'sub-0250'};
    subjIdx = cellfun(@(s) find(strcmp(da_full8,s)), subjects);
end
data = squeeze(meanBOLDpa_full([s0_idx,s90_idx,s45_idx,s135_idx], :, regionIndex, subjIdx)); % 4(direction) x 8(location) x nSubj

% weighted empirical mean across subjects, and the 4 asymmetry term
% values, at every (direction, location) cell -- including the
% direction+180 duplicate, since orientation stimuli are axial.
N = [0 45 90 135 180 225 270 315]; % the 8 compass "direction" slots
weightedBold = nan(numel(N), 8); % rows=direction(N), cols=location(pa, in anglevals order)
mainCardinalGrid = nan(numel(N), 8);
derivedCardinalGrid = nan(numel(N), 8);
mainSubsetGrid = nan(numel(N), 8);
derivedSubsetGrid = nan(numel(N), 8);

for li = 1:8
    pa = anglevals(li);
    for mi = 1:4
        md = mdirvals_dg(mi);
        vals = squeeze(data(mi, li, :))' .* subjectScale; % 1 x nSubj, gain-corrected
        valid = ~isnan(vals);
        wBold = sum(precisionW(valid) .* vals(valid)) / sum(precisionW(valid));

        mainCardinal = 2*ismember(md, maincardinalmDir) - 1;
        derivedCardinal = 2*((ismember(md,maincardinalmDir) & ismember(pa,primaryMeridians)) | ...
            (~ismember(md,maincardinalmDir) & ~ismember(pa,primaryMeridians))) - 1;
        mainSubset = 0; derivedSubset = 0;
        if strcmp(projectName, 'dg')
            if ismember(md,[0,180]); mainSubset=1; elseif ismember(md,[90,270]); mainSubset=-1; end
            if (abs(md-pa)==0 || abs(md-pa)==180); derivedSubset=1;
            elseif (abs(md-pa)==90 || abs(md-pa)==270); derivedSubset=-1; end
        else
            if ismember(md,[90,270]); mainSubset=1; elseif ismember(md,[0,180]); mainSubset=-1; end
            isCardMd = ismember(md,[0,90]); isOblMd = ismember(md,[45,135]);
            diffv = abs(md-pa);
            proH = (isCardMd && (diffv==90||diffv==270)) || (isOblMd && (diffv==0||diffv==180));
            conV = (isCardMd && (diffv==0||diffv==180)) || (isOblMd && (diffv==90||diffv==270));
            if proH; derivedSubset=1; elseif conV; derivedSubset=-1; end
        end

        % same (direction, direction+180) duplication as the original
        % script -- both compass slots get the identical term values and
        % weighted BOLD value, since all 4 terms are +-180-invariant.
        for d = [md, md+180]
            rowIdx = find(N == d);
            weightedBold(rowIdx, li) = wBold;
            mainCardinalGrid(rowIdx, li) = mainCardinal;
            derivedCardinalGrid(rowIdx, li) = derivedCardinal;
            mainSubsetGrid(rowIdx, li) = mainSubset;
            derivedSubsetGrid(rowIdx, li) = derivedSubset;
        end
    end
end

%% Plot: same "compass" spatial layout as lme2_ploteachDirLoc.m

% Square figure (Position width=height), so the axes' own normalized
% width/height fractions (p(3), p(4) below) already represent the SAME
% physical scale once 'axis square' is applied -- this is what makes the
% trigonometric shifts below give genuinely equal on-page distances, and
% is also what's needed for the saved PDF (below) to come out square.
figSize = 1250;
figure
set(gcf, 'Position', [1118 87 figSize figSize])
plot(0,0,'+k','MarkerSize',12, 'linewidth',3)
axis square
hold on
xlim([-4 4])
ylim([-4 4])
xticks([])
yticks([])
xticklabels({})
yticklabels({})
p = get(gca, 'Position');

box on
set(gca,'linewidth',1, 'YColor', [0 0 0]);
set(gca,'linewidth',1, 'XColor', [0 0 0]);

% Subplot size/placement solved analytically to maximize how much of the
% square placeholder axes (p, drawn above with the box/+  border) the 8
% polar subplots cover, leaving only a small gap between adjacent ones --
% replacing the old fixed-fraction w/h (previously only ~16-20% of p) and
% empirically-averaged radius. Each subplot is a SQUARE bounding box of
% side s (w=h, since a polar plot's circular content is best served by a
% square box, and uniform squares at evenly-spaced compass angles give a
% uniform gap all the way around, unlike the old unequal w/h).
%
% Two constraints, both solved for s and the placement radius R:
%   1) Boundary: the two CARDINAL subplots (angle 0/90/180/270) reach
%      farthest along one axis, at R + s/2 from center -- this must not
%      exceed p's own half-width, p(3)/2 (p is square by construction, so
%      the same limit applies in x and y).
%   2) No-overlap: for two adjacent (45-degree-separated) axis-aligned
%      squares both at radius R, they clear each other once their center
%      separation along EITHER axis reaches s (they only need to clear on
%      one axis, not both); the smaller of the two axis separations
%      between 45-degree-adjacent points is R*(1-cosd(45)), so the
%      tightest non-overlap bound is s <= R*(1-cosd(45))... the LARGER
%      axis separation, R*cosd(45), is the one that actually matters
%      (only one axis needs to clear), giving s <= R*cosd(45).
% gapFactor<1 backs off from the exact touching point (s=R*cosd(45)) by a
% small margin, at the boundary constraint's full extent -- solving the
% two simultaneously for a given gapFactor:
%   R + gapFactor*cosd(45)*R/2 = p(3)/2  =>  R = (p(3)/2) / (1 + gapFactor*cosd(45)/2)
gapFactor = 0.92; % how close adjacent subplots get to touching (1 = touching, 0 = no size at all)
R = (p(3)/2) / (1 + gapFactor*cosd(45)/2);
s = gapFactor * cosd(45) * R;
w = s; h = s;
left = p(1)+p(3)/2-(w/2); bottom = p(2)+p(4)/2-(h/2);
origin = [left bottom w h];

radius = R;
shifts = cell(1,8);
for k = 1:8
    shifts{k} = [radius*cosd(anglevals(k)), radius*sind(anglevals(k)), 0, 0];
end

violetColor = [148, 0, 211] / 255; % DarkViolet, replaces the previous plain red (dots)
violetColorLight = 0.2*violetColor + 0.8*[1 1 1]; % lighter tint for the connecting curve, so the line reads distinctly from the black dots (lightened further per feedback: was 0.5, decreased 30 points)

% Dot face size and line width (both marker outlines and the connecting
% curve), scaled to match plot1_experimentalCond.m's own dot-size-to-
% axes-width and line-width-to-axes-width ratios, applied to THIS
% script's own (differently-sized, unchanged) polar axes -- so the two
% scripts' dots/lines read as visually consistent despite the different
% subplot sizes. plot1_experimentalCond.m's polar axes are a fixed 4x4cm
% square with MarkerSize=6*0.8=4.8pts and LineWidth=proLineWidth=1pt (both
% MATLAB points, independent of screen DPI). Every one of this script's 8
% polar subplots shares the same SIZE (only position shifts per location,
% per the trigonometric placement above), so this only needs computing
% once, via a single throwaway axes at 'origin'.
tmpAx = axes(gcf, 'Units', 'normalized', 'Position', origin, 'Visible', 'off');
tmpAx.Units = 'points';
polarAxesWidth_pts = min(tmpAx.Position(3:4)); % smaller dimension = the constraining one for a circular polar plot
delete(tmpAx);
plot1_polarPlotWidth_pts = 4 * 72/2.54; % plot1_experimentalCond.m's polarPlotWidth_cm=4, converted to points
plot1_markerSize_pts = 6 * 0.8; % plot1_experimentalCond.m's markerSize
plot1_lineWidth_pts = 1; % plot1_experimentalCond.m's proLineWidth (COLORS.json value for mainCardinalVsMainOblique/dg, used as the reference stroke width)
dotSize = 2 * polarAxesWidth_pts * (plot1_markerSize_pts / plot1_polarPlotWidth_pts); % 2x per feedback (SizeData, i.e. area, doubled -- not diameter)
lineW = polarAxesWidth_pts * (plot1_lineWidth_pts / plot1_polarPlotWidth_pts); % connecting-curve width; back to 1x per feedback (was briefly 2x, now 2x thinner than that -- markerOutlineW below follows proportionally)
markerOutlineW = lineW / 4; % plot1_experimentalCond.m's dots use proLineWidth/4 for their white outline, not the full proLineWidth the curve uses

% Redundant-marker fade: which of the 8 DISPLAYED compass positions get
% faded is fixed by DISPLAY angle (0-135 = full opacity, 180-315 =
% faded), not by which raw direction code originally landed there --
% da's plotShift rotation would otherwise rotate the faded half around
% with each subplot, making the pattern inconsistent location to
% location. Computed per-location below (plotShift varies by pa); for dg
% (plotShift=0 always) this reduces to the raw-code split (mdirvals_dg
% vs the +180 duplicates), same thing either way. Heavily faded (95%
% white blend, not a subtle one) so it's unambiguous at these marker
% sizes.
alphaBlend = @(c) 0.35*c + 0.65*[1 1 1]; % was 0.05 (too faint per feedback), increased 30 points

globalMin = -0.5;
globalMax = 1;

modelPlotVals = nan(8,8); % rows=location(anglevals order), cols=direction(N order)

for li = 1:8 % location, matches anglevals(li) -- same subplot spatial position convention as the original (positioned by location index, not by the location's own polar angle)
    shift = shifts{li};

    pax = polaraxes(gcf);
    set(pax, 'Position', origin+shift)

    M = [ones(8,1), mainCardinalGrid(:,li), derivedCardinalGrid(:,li), mainSubsetGrid(:,li), derivedSubsetGrid(:,li)];
    predicted = M * estimatesVec;
    modelPlotVals(li,:) = predicted';

    % Smooth harmonic curve connecting the model dots, replacing the
    % original straight-line interpolation between them. Each subplot
    % here fixes theta_V = pa (this location's own polar angle); theta
    % is the stimulus ORIENTATION, swept continuously -- the roles of
    % theta_stim/theta_V are swapped relative to plot1_experimentalCond.m's
    % polar plots (there theta was the polar-angle location, with
    % orientation implicit in which discrete asymmetry was plotted), but
    % the underlying model is the exact same 4-term harmonic:
    %   y(theta_stim) = grandIntercept
    %       + b_mainCardinal    * cos(4*theta_stim)
    %       + b_derivedCardinal * cos(4*(theta_stim-theta_V))
    %       + b_mainSubset      * mainSubsetFn(theta_stim)
    %       + b_derivedSubset   * derivedSubsetFn(theta_stim,theta_V)
    % mainCardinal/derivedCardinal are project-independent (cardinal-
    % oblique is always a pure 4-fold function of orientation alone, or
    % of orientation-relative-to-location, regardless of project).
    % mainSubset/derivedSubset are NOT project-independent -- da's raw
    % predictors use a different sign/structure than dg's (da's
    % derivedSubset is a PRODUCT with the 4-fold cardinal term, not a
    % plain cos(2*diff)). Both forms verified (separately, exhaustively)
    % to reproduce the discrete mainSubsetGrid/derivedSubsetGrid values
    % above exactly -- zero absolute error against the ismember/abs-
    % difference logic at all 8 locations x 8 directions, both projects
    % -- before being used here for the continuous sweep.
    pa = anglevals(li);

    % da's stimuli (annulus/pinwheel/spirals) are defined in a POLAR
    % reference frame, not Cartesian -- the same raw orientation code
    % (e.g. 90 = pinwheel) has a DIFFERENT local Cartesian orientation
    % depending on which location it's shown at (a pinwheel is always
    % locally RADIAL, so its Cartesian angle rotates with location). The
    % regression betas and the mainCardinal/derivedCardinal/mainSubset/
    % derivedSubset predictor labels above are already correct -- the
    % design matrix they were fit on is location-based and already
    % accounts for this. This is purely a PLOTTING-angle correction:
    % verified against 32 (orientation, location) pairs, the true local
    % Cartesian angle is (rawCode + pa - 90) mod 180, i.e. every raw
    % direction code needs an additional rotation of (pa-90) degrees to
    % land at its true local Cartesian angle within THIS location's own
    % subplot (no rotation at the UVM itself, pa=90, where the polar and
    % Cartesian frames coincide by design). dg's gratings are already
    % defined directly in Cartesian terms, so no rotation applies there.
    if strcmp(projectName, 'da')
        plotShift = pa - 90;
    else
        plotShift = 0;
    end

    % Fixed DISPLAY-angle fade split for this location's markers (see
    % comment above alphaBlend's definition) -- N positions that land at
    % displayed angle >= 180 after rotation get faded, < 180 stay full
    % opacity, regardless of which raw direction code they came from.
    isMirrored = mod(N+plotShift,360) >= 180;

    thetaFine = 0:1:360;
    mainCardinal_f = cosd(4*thetaFine);
    derivedCardinal_f = cosd(4*(thetaFine-pa));
    if strcmp(projectName, 'dg')
        mainSubset_f = cosd(2*thetaFine);
        derivedSubset_f = cosd(2*(thetaFine-pa));
    else
        mainSubset_f = -cosd(2*thetaFine);
        derivedSubset_f = -cosd(4*thetaFine) .* cosd(2*(thetaFine-pa));
    end
    predictedSmooth = estimatesVec(1) + estimatesVec(2)*mainCardinal_f + estimatesVec(3)*derivedCardinal_f ...
        + estimatesVec(4)*mainSubset_f + estimatesVec(5)*derivedSubset_f;

    % The curve is exactly 180-periodic in theta_stim (every harmonic
    % term above is built from cos(2*theta)/cos(4*theta)), so its
    % second half is mathematically redundant with its first, same as
    % the markers -- split it the same way, by DISPLAY angle, using NaN
    % gaps (not a boolean subset) so each color's polarplot call stays
    % correctly connected across the wrap rather than drawing a stray
    % chord between the two disjoint runs a rotated split can produce.
    dispTheta = mod(thetaFine+plotShift,360);
    curveIsMirrored = dispTheta >= 180;
    rhoOriginalHalf = predictedSmooth; rhoOriginalHalf(curveIsMirrored) = NaN;
    rhoMirroredHalf = predictedSmooth; rhoMirroredHalf(~curveIsMirrored) = NaN;
    polarplot(deg2rad(dispTheta), rhoOriginalHalf, '-', 'Color', violetColor, 'LineWidth', lineW)
    hold on
    polarplot(deg2rad(dispTheta), rhoMirroredHalf, '-', 'Color', violetColorLight, 'LineWidth', lineW)
    hold on

    % Model dots disabled per feedback -- curve only. Left commented
    % (not deleted) so they're easy to bring back.
    % q1 = polarplot(deg2rad(mod(N(~isMirrored)+plotShift,360)), predicted(~isMirrored), 'o', ...
    %     'LineWidth', markerOutlineW, 'MarkerFaceColor', violetColor, 'MarkerSize', dotSize);
    % q1.MarkerEdgeColor = 'w';
    % hold on
    % q2 = polarplot(deg2rad(mod(N(isMirrored)+plotShift,360)), predicted(isMirrored), 'o', ...
    %     'LineWidth', markerOutlineW, 'MarkerFaceColor', alphaBlend(violetColor), 'MarkerSize', dotSize);
    % q2.MarkerEdgeColor = 'w';
    % hold on

    dataDotSize = dotSize * 0.6 * 2; % doubled per feedback (was 0.6x the old model-dot size, now model dots are gone -- this is 1.2x that reference size)
    dataDot = weightedBold(:,li);
    pdot1 = polarplot(deg2rad(mod(N(~isMirrored)+plotShift,360)), dataDot(~isMirrored), 'o', ...
        'LineWidth', markerOutlineW, 'MarkerFaceColor', 'black', 'MarkerSize', dataDotSize);
    pdot1.MarkerEdgeColor = 'w';
    hold on
    pdot2 = polarplot(deg2rad(mod(N(isMirrored)+plotShift,360)), dataDot(isMirrored), 'o', ...
        'LineWidth', markerOutlineW, 'MarkerFaceColor', alphaBlend([0 0 0]), 'MarkerSize', dataDotSize);
    pdot2.MarkerEdgeColor = 'w';
    hold on
    pax.FontSize = 6;
    pax.RTickLabel = {''};
    pax.ThetaTickLabel = {''};
    thetaticks(0:45:315);
    rlim([globalMin globalMax])
    hold on
end

sgtitle(sprintf('Regression-based est %s %s %s (gain+precision-weighted)', projectName, strrep(comparisonName,'_','-'), roiname), 'Interpreter', 'none')

% Explicit square PaperSize/PaperPosition, matching the figure's own
% square Position -- NOT '-bestfit', which scales/fits into a standard
% (non-square) page size and would undo the squareness set up above.
gcf_edit = gcf;
gcf_edit.Units = 'inches';
figDims = gcf_edit.Position(3:4); % square by construction (figSize x figSize)
gcf_edit.PaperUnits = 'inches';
gcf_edit.PaperPositionMode = 'manual';
gcf_edit.PaperSize = figDims;
gcf_edit.PaperPosition = [0, 0, figDims(1), figDims(2)];

print(gcf_edit, fullfile(figureDir, sprintf('EachDirLocRegression_%s_%s_%s', projectName, comparisonName, roiname)), '-dpdf');
close all;

end
