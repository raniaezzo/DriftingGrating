function plot1_experimentalCond(medianBOLDpa, asymmetryName, projectSettings, varargin)

    %rng(0)
    rng('default')

    % this function will (polar plot) the averages of experimental conditions that
    % do not need to be derived from polar angle.
    % For dg experiment, this will include cardinal v oblique
    % For da experiment, this will include polar cardinal v polar oblique
    %                                   and radial v tangential

    plotModelToo = 0; %1;
    if plotModelToo
        % just for internal purposes (if I want to plot model fits ad hoc)
        modelfolder = fullfile(projectSettings.glmResultsfolder, 'LME_results', projectSettings.comparisonName);
    end

    % for now, ensure a third argument is always provided
    if nargin < 3 || isempty(varargin{1})
        error('A third input is required when projectName is "%s" and asymmetryName is "%s".', projectSettings.projectName, asymmetryName);
    else
        subset = varargin{1};
    end

    % Check project name & request
    if strcmp(projectSettings.projectName, 'da') || strcmp(projectSettings.projectName, 'dots')
        if strcmp(asymmetryName, 'mainCardinalVsMainOblique')
            derivedVals = 0; % this is needed because radialVsTang occurs for both conditions
            % Ensure a third argument is provided
            if subset == 1
                asymmetryName = 'radialVsTangential';
                %radialvstang = 1;
            end
        elseif strcmp(asymmetryName, 'derivedCardinalVsDerivedOblique')
             % Third argument is ignored
             derivedVals = 1;
             if subset == 1
                asymmetryName = 'verticalVsHorizontal';
                %radialvstang = 0;
             end
        end
    elseif strcmp(projectSettings.projectName, 'dg')
        if strcmp(asymmetryName, 'mainCardinalVsMainOblique')
            % Third argument is ignored
            derivedVals = 0; 
            if subset == 1
                asymmetryName = 'verticalVsHorizontal';
                %radialvstang = 0;
            end
        elseif strcmp(asymmetryName, 'derivedCardinalVsDerivedOblique')
            derivedVals = 1;
            if subset == 1
                asymmetryName = 'radialVsTangential';
                %radialvstang = 1;
            end
        end
    else
        error('Unknown projectName. Expected "da" or "dg".');
    end

    % Which of fitAsymmetryRegression.m's 4 cached terms (termNames =
    % {mainCardinal, derivedCardinal, mainSubset, derivedSubset}) this call
    % corresponds to -- same formula as plot2_experimentalCond.m, needed
    % below for the model-based harmonic overlay curve.
    termIdx = 1 + 2*subset + derivedVals;

    projectName = projectSettings.projectName;
    colors_data = projectSettings.colors_data;
    rois = projectSettings.rois;
    comparisonName = projectSettings.comparisonName;
    axes_limits = projectSettings.axes_limits;
    figureDir = projectSettings.figureDir;

    % Locate fitAsymmetryRegression.m's cached output -- same convention as
    % plot2_experimentalCond.m -- needed for the model-based harmonic
    % overlay curve below (termIdx==1, cardinal-vs-oblique, for now).
    [summaryTablesDir,~,~] = fileparts(projectSettings.gainSummaryFile);
    [derivativesDir,~,~] = fileparts(summaryTablesDir);
    [bidsDir,~,~] = fileparts(derivativesDir);

    styleInfo = colors_data.conditions.(projectName).(asymmetryName);
    colors = {styleInfo.color_pro', styleInfo.color_con'};

    % pro = filled marker, solid line; con = unfilled marker, solid line,
    % same full color as pro (no opacity/blend manipulation) -- see
    % COLORS.json's pro_filled/con_filled/pro_lineWidth/con_lineWidth/
    % errorbar_lineWidth fields for this asymmetry (color_pro == color_con
    % by design: pro/con are distinguished by fill state and line weight,
    % not by hue).
    proFaceColor = colors{1};
    conFaceColor = [1 1 1]; % white fill (not unfilled)
    proEdgeColor = colors{1};
    conEdgeColor = colors{2};
    proLineStyle = '-';
    conLineStyle = '-';
    proLineWidth = styleInfo.pro_lineWidth;
    conLineWidth = styleInfo.con_lineWidth;
    errorbarLineWidth = styleInfo.errorbar_lineWidth;
    markerSize = 6 * 0.8; % previous (half-of-original) size, now 20% smaller
    axisLineWidth = 1; % matches plot2_experimentalCond.m's pairwise-plot axis lines
    axisFontSize = 10 * 0.8; % 20% smaller than MATLAB's default axes FontSize (10) -- no explicit size was set before
    showTitleLegend = false;

    % retrieve the indices for specific asymmetries (e.g., motion -
    % orientation for main cardinal v main oblique)
    if ~derivedVals
        [proConditions, conConditions, ~] = retrieveProConIdx(projectName, comparisonName, subset);
    else
        proConditions = 1; conConditions = 2;
    end

    % Accumulates the exact numbers underlying this figure (one row per
    % ROI x polar-angle location) across the ROI loop below, so a full
    % top-to-bottom run leaves a persistent, inspectable record of the
    % plotted values -- not just console output that scrolls away, and
    % not something you have to re-run (with its own fresh bootstrap
    % draws) to see again.
    statsRows = {};

    figure
%     gap = [.04 .01]; % spacing between the subplots vertical gap - horizontal gap
%     marg_h = [.015 .13]; % margins of bottom - top of the figure
%     marg_w = [.01 .01]; % margina - left, right
% %    [ha, pos] = tight_subplot(2, 1, gap, marg_h, marg_w);
%      [ha, pos] = tight_subplot(2, 5, gap, marg_h, marg_w);
%     % Hide the 8th subplot
% %         emptyPlots = 10 - numel(rois);
% %     for ep=1:emptyPlots
% %         empty_idx = 10+1-ep;
% %         set(ha(empty_idx), 'Visible', 'off');
% %     end
    for ri=1:length(rois)

        % Specify the region index
        regionIndex = projectSettings.roi_idx{ri};
        
        % Extract the relevant conditions for the specified region
        conditions1 = medianBOLDpa(proConditions, :, regionIndex, :);
        conditions2 = medianBOLDpa(conConditions, :, regionIndex, :);
        conditions1 = squeeze(conditions1);
        conditions2 = squeeze(conditions2);
        
        % Average across the conditions within subjects
%         if (strcmp(projectName, 'da') || strcmp(projectName, 'dots')) && strcmp(comparisonName, 'orientation_minus_baseline') && radialvstang == 1 || ...
%                 derivedVals == 1 %strcmp(asymmetryName, 'derivedCardinalVsDerivedOblique')

        if strcmp(comparisonName, 'orientation_minus_baseline') && subset == 1 || derivedVals ==1    % updated to include vertical
            avgConditions1 = conditions1;  % for orientation DA - there is only 1 value
            avgConditions2 = conditions2; 
        else
            avgConditions1 = nanmean(conditions1, 1);
            avgConditions2 = nanmean(conditions2, 1);
            avgConditions1 = squeeze(avgConditions1);
            avgConditions2 = squeeze(avgConditions2);
        end

        % Gain-weight each observer: divide their value at every polar
        % angle by their own mean pRF gain (prfvista_mov/prfvista average),
        % BEFORE any averaging across observers. This down-weights
        % high-gain observers and up-weights low-gain observers, and all
        % downstream stats/bootstrapped error bars below inherit the
        % adjustment automatically since they are computed from these
        % arrays. See projectSettings.gainWeightsSource / retrieveObserverGainWeights2.m
        %
        % The across-observer average gain is then multiplied back in, so
        % the plotted scale/units resemble the original (unweighted) data.
        % Dividing by gain_i and then multiplying by groupGain is the same
        % as dividing by gain_i normalized to the group mean (gain_i /
        % groupGain), so the relative weighting -- and therefore the
        % relative pattern across data points/error bars -- is unchanged;
        % only the overall scale shifts.
        %
        % groupGain uses the geometric, not arithmetic, mean: the applied
        % factor (groupGain ./ gainWeights) is itself a multiplicative
        % scale factor, and geometric mean is the choice under which the
        % *geometric* mean of that factor across observers is exactly 1
        % (magnitude-neutral in the multiplicative sense, matching what
        % the factor actually is) -- arithmetic mean does not have this
        % property. Implemented via exp(mean(log(.))) to avoid a
        % dependency on the Statistics and Machine Learning Toolbox.
        gainWeights = retrieveObserverGainWeights2(projectSettings.subjects, rois{ri}, projectSettings.gainWeightsSource);
        if any(gainWeights <= 0)
            error('gainWeights must be strictly positive to take log() for the geometric mean (found %d non-positive value(s))', ...
                sum(gainWeights <= 0));
        end
        groupGain = exp(mean(log(gainWeights), 'omitnan')); % omitnan: see retrieveObserverGainWeights2.m
        avgConditions1 = avgConditions1 .* (groupGain ./ gainWeights);
        avgConditions2 = avgConditions2 .* (groupGain ./ gainWeights);

        % Precision-weight each observer for THIS cortical area --
        % ROI-specific (reliability genuinely varies by cortical area,
        % unlike gain), via retrieveObserverPrecisionWeights.m.
        % projectSettings.precisionWeightsSource is currently [] everywhere
        % this is called (PLACEHOLDER: every (subject, roi) gets weight 1,
        % a no-op) -- defaults gracefully if that field isn't set at all,
        % so existing call sites keep working unchanged.
        %
        % This is a fundamentally different operation from gain
        % correction above: gain rescales each observer's VALUES before
        % averaging; precision weighting changes how much each observer's
        % (already gain-corrected) value COUNTS in the average -- applied
        % below as a proper weighted mean (point estimates) and weighted
        % paired bootstrap (CIs), not another data rescale.
        if isfield(projectSettings, 'precisionWeightsSource')
            precisionSource = projectSettings.precisionWeightsSource;
        else
            precisionSource = [];
        end
        precisionWeights = retrieveObserverPrecisionWeights(projectSettings.subjects, rois{ri}, precisionSource);

        % Extract polar angles
        %anglevals = [90, 135, 180, 225, 270, 315, 0, 45];
        anglevals = [90, 45, 0, 315, 270, 225, 180, 135]; % <-- these were manually converted based on the order of polarAngles above (Noah's convention)
        
        nBoot = 1000;
        nLoc  = size(avgConditions1,1);

        % Precision-weighted mean across observers (identical to nanmean
        % when precisionWeights is uniform, as it currently is -- see
        % weightedNanMean at the end of this file). Transposed to match
        % nanmean(...,2)''s previous row-vector shape (nLoc columns).
        vals_1 = weightedNanMean(avgConditions1, precisionWeights)';
        %sem1 = nanstd(avgConditions1,0,2)' ./ sqrt(sum(~isnan(avgConditions1),2)');
        vals_2 = weightedNanMean(avgConditions2, precisionWeights)';
        %sem2 = nanstd(avgConditions2,0,2)' ./ sqrt(sum(~isnan(avgConditions2),2)');
 
        disp('Testing consistency with permutation')
        D_obs = vals_1 - vals_2;
        disp('number greater (asymmetry count):')
        D_obs>0
        observed_stat = mean(D_obs)

        % added to plot bootstrapped SED
        SED = nan(1,nLoc);
        CI_95_lower = nan(1,nLoc);
        CI_95_upper = nan(1,nLoc);

        for loc = 1:nLoc
        
            % Data for this location
            x1 = avgConditions1(loc,:);
            x2 = avgConditions2(loc,:);
        
            % Keep only subjects with data in BOTH conditions
            valid = ~isnan(x1) & ~isnan(x2);
            x1 = x1(valid);
            x2 = x2(valid);
            % Forced to a row (x1/x2 are rows): MATLAB vector indexing
            % preserves the INDEXED array's own orientation, not idx's, so
            % if w_loc stayed a column (precisionWeights' native shape)
            % while x1/x2 are rows, wb.*x1(idx) below would silently
            % broadcast into an nSub x nSub matrix instead of an
            % elementwise product.
            w_loc = precisionWeights(valid);
            w_loc = w_loc(:)';

            nSub = numel(x1);

            bootDiff = zeros(nBoot,1);

            for b = 1:nBoot

                % Bootstrap subjects WITH replacement
                idx = randi(nSub,nSub,1);

                % Precision-weighted difference of means for this
                % bootstrap sample (identical to the previous plain
                % mean(x1(idx))-mean(x2(idx)) when precisionWeights is
                % uniform, as it currently is). Weights renormalized
                % within each draw since resampling can duplicate/omit
                % subjects, changing which weights are in play.
                wb = w_loc(idx); wb = wb / sum(wb);
                bootDiff(b) = sum(wb.*x1(idx)) - sum(wb.*x2(idx));

            end
        
            % 95% percentile confidence interval
            CI = prctile(bootDiff,[2.5 97.5]);

            CI_95_lower(loc) = CI(1);
            CI_95_upper(loc) = CI(2);

            % 68% percentile confidence interval
            CI_68 = prctile(bootDiff,[16 84]);

            CI_68_lower(loc) = CI_68(1);
            CI_68_upper(loc) = CI_68(2);

        end

        CI_95_lower
        CI_95_upper

        %

%         n_perms = 10000;
%         null_stats = zeros(n_perms, 1);
%         for i = 1:n_perms
%             flip = randi([0, 1], 1, 8) * 2 - 1;  % randomly flip sign
%             permuted_diff = flip .* (D_obs);
%             %null_stats(i) = sum(permuted_diff > 0); % for 1 tailed
%             null_stats(i) = mean(permuted_diff); % for two-tailed
%         end
% 
%         %p = mean(null_stats >= observed_stat) % for 1 tailed
%         p = mean(abs(null_stats) >= abs(observed_stat));  % two-tailed
% 
%         if p==0
%             p = 1 / n_perms
%         else
%             p
%         end

        mean_diff = mean(D_obs)

%         std_effect = std(D_obs)
% 
%         cohens_d = mean_diff/std_effect
%         disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

        % using this method instead of subplot for TIGHT AXES
%         axes(ha(ri)); 
        
        % Half-width of the 68% percentile CI of the bootstrapped (pro-con)
        % difference, applied symmetrically to both markers -- same
        % simplification the previous std(bootDiff)-based error bar used
        % (a single per-location spread applied to both conditions), but
        % now driven by the actual percentile CI rather than the SD.
        ci68_halfwidth = (CI_68_upper - CI_68_lower) / 2;

        % record this ROI's per-location values for the stats export below
        for loc = 1:nLoc
            statsRows(end+1,:) = {rois{ri}, anglevals(loc), vals_1(loc), vals_2(loc), D_obs(loc), ...
                CI_95_lower(loc), CI_95_upper(loc), CI_68_lower(loc), CI_68_upper(loc), mean_diff}; %#ok<SAGROW>
        end

        % Draw order (bottom to top): model curves, then markers, then
        % error bars on top -- the model curves are the group-level
        % prediction, not the data itself, so they read best as a
        % background trend with the actual data (dots + error bars) drawn
        % on top of them, not obscuring them. Lines and markers are
        % deliberately drawn in SEPARATE calls (lines with no marker;
        % markers with no line) -- the previous combined 'o'+line calls
        % drew a marker at both the "closing segment" endpoints AND again
        % at those same points in the "all data" call, double-stroking the
        % outline at those 2 of 8 locations only (fill looked fine since a
        % flat opaque fill drawn twice looks identical to once, but the
        % outline stroke rendered visibly heavier there than at the other
        % 6 points). Drawing each marker exactly once fixes that.

        % Straight-line-segment connectors between the 8 dots used to be
        % drawn here (both the closing segment and the main 8-point
        % polyline, for pro and con) -- removed now that the model-based
        % continuous harmonic curve below (termIdx==1 so far) replaces
        % them with the regression's own predicted shape instead of a
        % naive straight-line interpolation between data points.

        % Model-based continuous overlay (replaces the old, now-superseded
        % plotModelToo dots -- see that block below, still present but
        % dead code with plotModelToo=0): the regression's own predicted
        % value as a SMOOTH function of location, not just the 8 discrete
        % dots connected by straight segments. Wired up for all 4 terms
        % (termIdx 1-4).
        %
        % termIdx==1 (cardinal-oblique): averaging the full model over
        % cardinal's 2 constituent directions (0/90 deg) exactly cancels
        % mainSubset/derivedSubset at every location, leaving
        %   f_pro(theta) = grandIntercept + beta(mainCardinal) + beta(derivedCardinal)*cos(4*theta)
        %   f_con(theta) = grandIntercept - beta(mainCardinal) - beta(derivedCardinal)*cos(4*theta)
        % cos(4*theta) is the exact continuous form of derivedCardinal's
        % +-1 coding at a cardinal reference direction (0 deg or 90 deg
        % give the IDENTICAL curve -- both are in maincardinalmDir, so
        % derivedCardinal's own formula collapses to the same thing either
        % way; same for oblique's 45/135 reference) -- verified to
        % reproduce +1/-1 at all 8 discrete wedge centers, matching
        % fitAsymmetryRegression.m's own predictor definition exactly.
        %
        % termIdx==2 (polar cardinal-oblique): the SAME 2nd-order harmonic
        % as cardinal-oblique (period 90 deg), just with the roles of the
        % two betas swapped -- it peaks/troughs at the opposite locations.
        % In termIdx==1's curve, beta(derivedCardinal) is already the
        % additive cos(4*theta)-varying term riding on top of
        % beta(mainCardinal)'s constant offset; polar-cardinal-oblique is
        % the mirror image of that same pair: beta(derivedCardinal)
        % becomes the constant term and beta(mainCardinal) becomes the
        % cos(4*theta)-varying one:
        %   f_pro(theta) = grandIntercept + beta(derivedCardinal) + beta(mainCardinal)*cos(4*theta)
        %   f_con(theta) = grandIntercept - beta(derivedCardinal) - beta(mainCardinal)*cos(4*theta)
        % Derived by averaging the full 4-beta model over the two real
        % directions theta_stim = theta_V and theta_stim = theta_V+90 (the
        % two directions that ARE "matching"/radial+tangential at location
        % theta_V, per compute_derivativeDirections.m's own definition of
        % polar cardinal) -- this cancels mainSubset/derivedSubset exactly,
        % the same cancellation mechanism as termIdx==1, just anchored to
        % theta_V instead of a fixed 0/90 reference. A previous version of
        % this formula used a piecewise sign(cos(4*theta))/abs(cos(4*theta))
        % construction that matched the 8 discrete points but introduced a
        % real discontinuity at the octant boundaries -- wrong, since this
        % is genuinely just another smooth 2nd-order harmonic, not a
        % step function. Verified to exactly reproduce the discrete pro/con
        % values AND vary smoothly (no jump) between them.
        %
        % termIdx==3 (horizontal-vertical, mainSubset) and termIdx==4
        % (radial-tangential, derivedSubset) are structurally different
        % from 1/2: pro/con here are each a SINGLE real direction (not an
        % averaged category), so there is no cancellation of the other 2
        % betas -- both curves genuinely contain all 4 betas, and pro/con
        % share (rather than mirror) the mainCardinal/derivedCardinal
        % contribution:
        %   termIdx==3: pro=horizontal (theta_stim=0 fixed), con=vertical (theta_stim=90 fixed)
        %     f_pro(theta) = grandIntercept + beta(mainCardinal) + beta(mainSubset) + beta(derivedCardinal)*cos(4*theta) + beta(derivedSubset)*cos(2*theta)
        %     f_con(theta) = grandIntercept + beta(mainCardinal) - beta(mainSubset) + beta(derivedCardinal)*cos(4*theta) - beta(derivedSubset)*cos(2*theta)
        %   termIdx==4: pro=radial (theta_stim=theta_V, i.e. motion direction
        %     pointing along the location itself), con=tangential (theta_stim=theta_V+90)
        %     f_pro(theta) = grandIntercept + beta(derivedCardinal) + beta(derivedSubset) + beta(mainCardinal)*cos(4*theta) + beta(mainSubset)*cos(2*theta)
        %     f_con(theta) = grandIntercept + beta(derivedCardinal) - beta(derivedSubset) + beta(mainCardinal)*cos(4*theta) - beta(mainSubset)*cos(2*theta)
        % Both derived the same way (evaluate the full 4-beta continuous
        % harmonic model at the fixed or location-tracking theta_stim that
        % defines pro/con for this term) and verified to exactly reproduce
        % the real discrete predictor values at all 8 wedge centers.
        if termIdx >= 1 && termIdx <= 4
            if isfield(projectSettings, 'fitLabel') && ~isempty(projectSettings.fitLabel)
                fitLabelModel = projectSettings.fitLabel;
            else
                fitLabelModel = projectName;
            end
            fitFileModel = fullfile(bidsDir,'derivatives','summaryTables','regressionResults',fitLabelModel,sprintf('%s.mat',rois{ri}));
            if ~isfile(fitFileModel)
                error('plot1_experimentalCond:missingFit', ...
                    'Cached fit not found for %s / %s at %s -- run fitAsymmetryRegression(''%s'') first.', ...
                    fitLabelModel, rois{ri}, fitFileModel, projectName);
            end
            Fmodel = load(fitFileModel);
            beta1 = Fmodel.estimates(1)/2; % mainCardinal
            beta2 = Fmodel.estimates(2)/2; % derivedCardinal
            beta3 = Fmodel.estimates(3)/2; % mainSubset
            beta4 = Fmodel.estimates(4)/2; % derivedSubset

            thetaFine = linspace(0, 360, 361); % 1-degree resolution, smooth curve
            switch termIdx
                case 1
                    modelPro = Fmodel.grandInterceptFE + beta1 + beta2*cosd(4*thetaFine);
                case 2
                    modelPro = Fmodel.grandInterceptFE + beta2 + beta1*cosd(4*thetaFine);
                case 3
                    modelPro = Fmodel.grandInterceptFE + beta1 + beta3 + beta2*cosd(4*thetaFine) + beta4*cosd(2*thetaFine);
                case 4
                    modelPro = Fmodel.grandInterceptFE + beta2 + beta4 + beta1*cosd(4*thetaFine) + beta3*cosd(2*thetaFine);
            end

            % Line widths match the original connecting-line widths exactly
            % (proLineWidth/conLineWidth genuinely differ per COLORS.json --
            % con is thinner than pro for this asymmetry, same as
            % everywhere else in this pipeline).
            polarplot(deg2rad(thetaFine), modelPro, '-', 'Color', proEdgeColor, 'LineWidth', proLineWidth);
            hold on

            % Con: a manually-constructed dash pattern (NaN-separated
            % segments), not MATLAB's native ':' -- native LineStyle
            % exposes no control over dash-length vs. gap independently.
            % dashDeg is the visible segment length, gapDeg the invisible
            % break between segments (both in degrees of the 360-degree
            % circle); gapDeg is half of the visual gap the native ':'
            % was producing, dashDeg left matching that same look.
            % Tunable if the ratio needs further adjustment.
            dashDeg = 3; gapDeg = 1.5;
            unitDeg = dashDeg + gapDeg;
            thetaCon = [];
            for u = 0:ceil(360/unitDeg)-1
                segStart = u*unitDeg;
                if segStart >= 360, break; end
                segEnd = min(segStart+dashDeg, 360);
                thetaCon = [thetaCon, segStart:0.5:segEnd, NaN]; %#ok<AGROW>
            end
            switch termIdx
                case 1
                    modelConDashed = Fmodel.grandInterceptFE - beta1 - beta2*cosd(4*thetaCon);
                case 2
                    modelConDashed = Fmodel.grandInterceptFE - beta2 - beta1*cosd(4*thetaCon);
                case 3
                    modelConDashed = Fmodel.grandInterceptFE + beta1 - beta3 + beta2*cosd(4*thetaCon) - beta4*cosd(2*thetaCon);
                case 4
                    modelConDashed = Fmodel.grandInterceptFE + beta2 - beta4 + beta1*cosd(4*thetaCon) - beta3*cosd(2*thetaCon);
            end
            polarplot(deg2rad(thetaCon), modelConDashed, '-', 'Color', conEdgeColor, 'LineWidth', conLineWidth);
            hold on
        end

        % Plot 68% CI per point -- drawn BEFORE the markers (below them in
        % z-order) so the pro marker's white outline (below) sits on top
        % of the whisker instead of the whisker cutting across it, same
        % fix used for the mean dot in plot2_experimentalCond.m's
        % subjectwiseDiff mode. Error bar color/alpha matches its own
        % condition (pro = full color, con = same 50%-white-blended color
        % used for con's marker/line).
        p1 = polarplot([deg2rad(anglevals); deg2rad(anglevals)], [vals_1 - ci68_halfwidth; vals_1 + ci68_halfwidth], '-', 'Color', proEdgeColor, 'LineWidth',errorbarLineWidth);
        hold on
        p2 = polarplot([deg2rad(anglevals); deg2rad(anglevals)], [vals_2 - ci68_halfwidth; vals_2 + ci68_halfwidth], '-', 'Color', conEdgeColor, 'LineWidth',errorbarLineWidth);
        hold on

        % Markers, exactly once per location, drawn on top of the model
        % curves and error bars above. Pro is solid-filled with a thin
        % white outline (half of plot2_experimentalCond.m's already-halved
        % subjectwiseDiff outline, i.e. proLineWidth/4) rather than its own
        % colored edge -- same "filled dot, white outline" convention used
        % there. Its MarkerSize is grown from markerSize to
        % markerSize+proLineWidth so the face alone now spans what the
        % OLD pro dot's (and the unchanged con/unfilled dot's) total
        % face+outline footprint occupied. Con stays exactly as before
        % (white face, colored edge, at the original markerSize/
        % proLineWidth) -- it's the "unfilled" dot whose current size is
        % the sizing reference for pro's new face.
        whiteOutlineLineWidth_ = proLineWidth/4;
        proMarkerSizeFilled_ = markerSize + proLineWidth;
        polarplot(deg2rad(anglevals),vals_1, 'o', 'LineStyle', 'none', 'MarkerSize', proMarkerSizeFilled_, 'MarkerFaceColor', proFaceColor, 'MarkerEdgeColor', 'w', 'LineWidth', whiteOutlineLineWidth_)
        hold on
        polarplot(deg2rad(anglevals),vals_2, 'o', 'LineStyle', 'none', 'MarkerSize', markerSize, 'MarkerFaceColor', conFaceColor, 'MarkerEdgeColor', conEdgeColor, 'LineWidth',proLineWidth)
        hold on

%        for subjectIndex = 1:size(medianBOLDpa, 4)
%             %p3  = polarplot(deg2rad(anglevals), avgConditions1(:, subjectIndex)', 'o', 'Color', [127/255, 191/255, 123/255]);
%             %hold on;
%             %p4 = polarplot(deg2rad(anglevals), avgConditions2(:, subjectIndex)', 'o', 'Color', [175/255, 141/255, 195/255]);
%             %hold on
%         %     hold on
% %             asymm = avgConditions1(:, subjectIndex) - avgConditions2(:, subjectIndex);
% %             if (sum(asymm<0)) ~= 0
% %                 sprintf('WaRNING: %s points are below 0', num2str(sum(asymm<0)))
% %             end
%             
%         end
        
        thetaticks(0:45:315);
    
    
        if ri==1 && showTitleLegend
            if strcmp(asymmetryName, 'radialVsTangential')
                hLegend = legend('Radial', 'Tangential', 'Location', 'northwest', 'Box', 'off', 'FontSize', 18);
            elseif strcmp(asymmetryName, 'verticalVsHorizontal')
                hLegend = legend('Horizontal', 'Vertical', 'Location', 'northwest', 'Box', 'off', 'FontSize', 18);
            else
                hLegend = legend('Cardinal', 'Oblique', 'Location', 'northwest', 'Box', 'off', 'FontSize', 18);
            end
            % Adjust the position of the legend
            newPosition = get(hLegend, 'Position'); % Get current position
            newPosition(1) = newPosition(1) + 0.1; % Shift the legend to the right by 0.1 normalized units
            newPosition = [0.7858 0.8972 0.1929 0.1131];
            nlegend = set(hLegend, 'Position', newPosition);
        end
        
        ax = gca;
            
        if strcmp(rois{ri}, 'hMTcomplex') || strcmp(rois{ri}, 'pMT') || strcmp(rois{ri}, 'pMST')
            ax.RLim = [0 2]; %[0 2];%[0 2]; %[1 3]; %
            ax.RLim = [axes_limits.(projectName).(comparisonName).ROIs_motion.min ...
                axes_limits.(projectName).(comparisonName).ROIs_motion.max];
    %     elseif ri==3  % remove later (only for orientation
    %         ax.RLim = [-0.25 2];
        else
            ax = gca;
            ax.RLim = [-.25 1.75]; %[-1 1]; %[-.25 1.75]; %[0 2]; %
            ax.RLim = [axes_limits.(projectName).(comparisonName).ROIs_early.min ...
                axes_limits.(projectName).(comparisonName).ROIs_early.max];
        end
    
        if strcmp(projectName, 'dots')
            ax.RLim = [0 0.5];
        end

        if plotModelToo
            load(fullfile(modelfolder, rois{ri}, 'modelPlotVals.mat')) % all already in absolute dir (dg and da)
            load(fullfile(modelfolder, 'modelPlotValsTags.mat')) % all already in absolute dir (dg and da)
            
            
            if strcmp(asymmetryName, 'verticalVsHorizontal')
                tagshow = tags(:,:,1);
                flipProCon = 1; % I need this until I reprocess the model with the correct Pro/Con order
            % cardinal vs oblique model fits
            elseif (strcmp(asymmetryName, 'mainCardinalVsMainOblique') && strcmp(projectSettings.projectName, 'dg')) || ...
                (strcmp(asymmetryName, 'derivedCardinalVsDerivedOblique') && strcmp(projectSettings.projectName, 'da'))
                tagshow = tags(:,:,2);
                flipProCon=0;
            elseif strcmp(asymmetryName, 'radialVsTangential')
                tagshow = tags(:,:,3);
                flipProCon=0;
            % polar cardinal vs polar oblique model fits
            elseif (strcmp(asymmetryName, 'mainCardinalVsMainOblique') && strcmp(projectSettings.projectName, 'da')) || ...
                (strcmp(asymmetryName, 'derivedCardinalVsDerivedOblique') && strcmp(projectSettings.projectName, 'dg'))
                tagshow = tags(:,:,4);
                flipProCon=0;
            end

            a = sum(modelPlotVals .* (tagshow == 1), 2) ./ sum(tagshow == 1, 2);
            b = sum(modelPlotVals .* (tagshow == -1), 2) ./ sum(tagshow == -1, 2);
            if flipProCon
                conMeans = a; proMeans=b;
            else
                proMeans=a; conMeans=b;
            end

            % plot pro and con for each location afer averaging:
            hold on % marker size was 6
            polarplot(deg2rad(0:45:315),proMeans, 'o', 'Color', colors{1},  'MarkerSize', 6, 'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'k','LineWidth',1.75)
            hold on
            polarplot(deg2rad(0:45:315),conMeans, 'o', 'Color', colors{2},  'MarkerSize', 6, 'MarkerFaceColor', colors{2}, 'MarkerEdgeColor', 'k','LineWidth',1.75)
        end

        ax = gca;
        ax.LineWidth = axisLineWidth;
        ax.FontSize = axisFontSize; % axis (tick) label size
        ax.GridColor = [0.25 0.25 0.25];
        ax.ThetaColor = [0.25 0.25 0.25]; % theta-axis line/tick/label color -- matched in plot2_experimentalCond.m's XColor/YColor
        ax.RColor = [0.25 0.25 0.25]; % r-axis line/tick/label color
        % ax.Layer left at its default ('bottom') -- 'top' drew grid/spoke
        % lines over solid data points. The rho-label overlap this was
        % originally meant to fix is instead solved by RAxisLocation
        % above, which moves the labels into an empty wedge with no data.
        ax.ThetaTickLabel = {};
        ax.RTick = [-0.5 0 0.5 1]; % explicit rho values, rather than relying on auto ticks matching these
        ax.RAxisLocation = 67.5; % midway between 45 and 90 deg -- an empty wedge, away from the 8 datapoint locations, so rho labels don't overlap data
        ax.Box = 0;
        %ax.RTickLabel = [];
        if showTitleLegend
            title(rois{ri}, 'FontSize', 18)
        end
    end

    fig1 = gcf;
    if showTitleLegend
        sgtitle(sprintf('%s: %s', projectName, strrep(comparisonName, "_", " ")), 'FontSize', 40)
    end
    fig1.Position = [34 228 1210 924]; %[152 569 2143 619];
    fig1.Color = 'w';
    hold off;
    
    % ROI name(s) plotted in this call, joined with '-' -- included in the
    % filename so re-running for a different cortical area (or set of
    % areas) doesn't silently overwrite another area's file. Today this is
    % always just the single ROI in projectSettings.rois (production is
    % currently scoped to rois(1), i.e. V1 only), but this works
    % unchanged if that scope is ever widened to multiple/other areas.
    roiNameStr = strjoin(rois, '-');
    if plotModelToo==1
        filename = fullfile(figureDir,sprintf('polarangle_%s_%s_%s_%s_wModel', comparisonName, projectName, asymmetryName, roiNameStr));
    else
        filename = fullfile(figureDir,sprintf('polarangle_%s_%s_%s_%s', comparisonName, projectName, asymmetryName, roiNameStr));
    end

    % Save the exact values underlying this figure alongside the PDF, so
    % a full top-to-bottom run leaves a persistent, reproducible record
    % (diffable across runs) rather than requiring a re-run -- with its
    % own fresh bootstrap draws -- to see the numbers again.
    statsTable = cell2table(statsRows, 'VariableNames', ...
        {'roi','polarAngle','pro_mean','con_mean','diff','ci95_lower','ci95_upper','ci68_lower','ci68_upper','meanDiff_acrossLocations'});
    writetable(statsTable, [filename, '_stats.csv']);

    % The axis (plot) box itself is the sized element (4 x 4 cm) -- the
    % PDF page is padded larger around it (padding_cm each side) so the
    % axis box is guaranteed exactly this size regardless of how much
    % room MATLAB's auto layout would otherwise reserve for tick labels.
    % NOT fitFig2Page, which scales the figure to fill a full
    % letter-landscape page and would override any custom size set here.
    polarPlotWidth_cm = 4;
    polarPlotHeight_cm = 4;
    basePadding_cm = 0.5; % previous padding, on each side
    areaScale = 1.25; % total figure area, not just padding, should grow by this factor
    % Solve for the (still-uniform, all 4 sides) padding p that gives
    % (axisW+2p)*(axisH+2p) = areaScale * (axisW+2*basePadding)*(axisH+2*basePadding),
    % i.e. a quadratic in p: 4p^2 + 2p(axisW+axisH) + axisW*axisH - targetArea = 0
    targetArea = areaScale * (polarPlotWidth_cm + 2*basePadding_cm) * (polarPlotHeight_cm + 2*basePadding_cm);
    qa = 4; qb = 2*(polarPlotWidth_cm + polarPlotHeight_cm); qc = polarPlotWidth_cm*polarPlotHeight_cm - targetArea;
    padding_cm = (-qb + sqrt(qb^2 - 4*qa*qc)) / (2*qa);
    figWidth_cm = polarPlotWidth_cm + 2*padding_cm;
    figHeight_cm = polarPlotHeight_cm + 2*padding_cm;

    gcf_edit = gcf;
    gcf_edit.Units = 'centimeters';
    gcf_edit.Position(3:4) = [figWidth_cm, figHeight_cm];
    gcf_edit.PaperUnits = 'centimeters';
    gcf_edit.PaperPositionMode = 'manual'; % otherwise 'auto' ignores PaperPosition below
    gcf_edit.PaperSize = [figWidth_cm, figHeight_cm];
    gcf_edit.PaperPosition = [0, 0, figWidth_cm, figHeight_cm];

    ax.Units = 'centimeters';
    ax.Position = [padding_cm, padding_cm, polarPlotWidth_cm, polarPlotHeight_cm];

    % Save as PDF
    set(gcf_edit,'Renderer','painters'); % new
    print(gcf_edit, filename, '-dpdf', '-vector'); %'-painters');
    close all;

end

function m = weightedNanMean(X, w)
% WEIGHTEDNANMEAN  Precision-weighted mean of X (nLoc x nSubj) across
% subjects (dim 2), ignoring NaNs, using per-subject weights w (nSubj x 1
% or 1 x nSubj). Returns an nLoc x 1 column. Identical to nanmean(X,2)
% when w is uniform (as it currently is, pending finalized precision
% weights) -- NaN entries get zero weight in both the numerator and the
% weight-sum denominator, so each location's mean is taken only over that
% location's own valid (non-NaN) subjects, same as nanmean.
    w = w(:)'; % row, so it broadcasts against X's columns (subjects)
    validMask = ~isnan(X);
    Wexpanded = repmat(w, size(X,1), 1);
    Wexpanded(~validMask) = 0;
    Xz = X;
    Xz(~validMask) = 0;
    m = sum(Xz .* Wexpanded, 2) ./ sum(Wexpanded, 2);
end