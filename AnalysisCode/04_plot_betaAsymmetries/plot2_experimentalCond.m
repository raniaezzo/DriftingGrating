function plot2_experimentalCond(medianBOLDpa, asymmetryName, projectSettings, varargin)

    rng('default')
    %rng(0)
    % Plot mean across polar angles
    % keep in mind that this equally weighs each PA, whereas there could be
    % differential # of voxels representing the PAs

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
    
    projectName = projectSettings.projectName;
    comparisonName = projectSettings.comparisonName;
    colors_data = projectSettings.colors_data;
    rois = projectSettings.rois;
    %pairaxes_limits = projectSettings.pairaxes_limits;
    pairaxes_PAew_limits = projectSettings.pairaxes_PAew_limits;
    figureDir = projectSettings.figureDir;

    % Which of fitAsymmetryRegression.m's 4 cached terms (termNames =
    % {mainCardinal, derivedCardinal, mainSubset, derivedSubset}) this
    % call corresponds to. derivedVals (0/1) selects the main-vs-derived
    % pair, subset (0/1) selects which of that pair; this is the SAME raw
    % term regardless of project (dg/da) -- verified numerically to match
    % this script's own (equal-weighted) point estimate exactly for V1,
    % both projects, both subset values (see AnalysisCode session notes /
    % validate_termMapping.m). No dg/da concept-reindexing is needed here
    % (that's only required in fitAsymmetryRegression_dgVsDa.m, which
    % aligns dg's and da's raw term slots to a shared concept order before
    % differencing across projects -- this function only ever looks at
    % one project's own cached fit).
    termIdx = 1 + 2*subset + derivedVals;

    % Display mode (5th positional arg, optional): toggles BOTH the dots
    % and the grey lines together between two internally-consistent
    % readings of the same figure --
    %   'model' (default): dots = grandInterceptFE +/- estimates/2 (the
    %     regression's own group-level asymmetry estimate); grey line
    %     slope = that subject's rescaled contribution to the regression
    %     (F.subjectContributions * F.nSubj); both centered on
    %     grandInterceptFE.
    %   'raw': dots = vals_1_overall/vals_2_overall (precision-weighted
    %     average of gain-corrected data); grey line slope = that
    %     subject's own raw (vals_1(s)-vals_2(s)) difference; both
    %     centered on mean([vals_1_overall,vals_2_overall]) -- the SAME
    %     precision-weighted anchor the dots use, not a separately
    %     unweighted grand mean.
    % Both computations are always performed regardless of mode (see
    % modelPro/modelCon and vals_1_overall/vals_2_overall below) so the
    % two are directly comparable -- verified to be numerically identical
    % for complete-data cortical areas (V1/V2, both projects, all 4
    % asymmetries) precisely because that's the condition under which the
    % pooled regression coefficient equals a precision-weighted average of
    % individual differences (see AnalysisCode session notes). Superedes
    % the older 'centeringMode' parameter (its 'raw'-per-subject-own-level
    % option is dropped -- no production caller ever used it, and the new
    % shared-precision-weighted-anchor convention subsumes its
    % 'recentered' option, which is now `displayMode='raw'`'s centering).
    if nargin >= 5 && ~isempty(varargin{2})
        displayMode = varargin{2};
    else
        displayMode = 'model';
    end
    if ~ismember(displayMode, {'model','raw'})
        error('plot2_experimentalCond:displayMode', 'displayMode must be ''model'' or ''raw'' (got ''%s'').', displayMode);
    end

    % plotType (6th positional arg, optional): toggles the WHOLE plot
    % style, independent of displayMode (which only matters for
    % 'pairwise').
    %   'pairwise' (default): the existing two-dot-plus-grey-lines figure,
    %     unchanged.
    %   'subjectwiseDiff': a different rendering of the SAME underlying
    %     asymmetry -- at x=1, each subject's own raw (pro-con) difference
    %     across all 8 locations (no intercept/model subtraction, no
    %     connecting lines), with x-jitter; at x=2, the group mean
    %     difference (same cached bootstrap fit as 'pairwise' uses) with a
    %     95% CI error bar (not 68% -- this is now a test of whether the
    %     difference itself is nonzero, not a pro-vs-con comparison, so
    %     the narrower 68% band used for visually separating two dots
    %     doesn't apply). Plus a thick dashed line at y=0.
    if nargin >= 6 && ~isempty(varargin{3})
        plotType = varargin{3};
    else
        plotType = 'pairwise';
    end
    if ~ismember(plotType, {'pairwise','subjectwiseDiff'})
        error('plot2_experimentalCond:plotType', 'plotType must be ''pairwise'' or ''subjectwiseDiff'' (got ''%s'').', plotType);
    end

    % Locate fitAsymmetryRegression.m's cached output (same bidsDir this
    % whole pipeline reads/writes; derived from gainSummaryFile rather
    % than re-hardcoded here so there's a single source of truth for the
    % path).
    [summaryTablesDir,~,~] = fileparts(projectSettings.gainSummaryFile);
    [derivativesDir,~,~] = fileparts(summaryTablesDir);
    [bidsDir,~,~] = fileparts(derivativesDir);

    styleInfo = colors_data.conditions.(projectName).(asymmetryName);
    colors = {styleInfo.color_pro', styleInfo.color_con'};

    % pro = filled marker, full color; con = unfilled marker, same
    % 50%-white-blended color used in plot1_experimentalCond.m's polar
    % plots (color_pro == color_con by design -- pro/con distinguished by
    % fill state and this alpha-equivalent blend, not by hue). Blending
    % toward white rather than true alpha for the same reason as the
    % polar plots: painters-rendered vector PDF export doesn't support
    % line/marker transparency.
    proFaceColor = colors{1};
    conFaceColor = [1 1 1]; % white fill (not unfilled)
    proEdgeColor = colors{1};
    conEdgeColor = 0.5*colors{2} + 0.5*[1 1 1];
    subjectLineColor = [0.8196, 0.8275, 0.8314]; % RGB(209,211,212)
    % subjectwiseDiff mode only: subjects present in dg but not da (only
    % populated when dg is run in its 'all' 13-subject mode) get an 'x'
    % marker instead of a dot in the per-subject scatter -- this is dg's
    % own 7-subject overlap with da, i.e. dg's 13 minus its 6 dg-only
    % subjects (5 dg-only subjects plus sub-0395, which da excludes
    % entirely as a pilot-mismatch subject and so never contributes to da).
    dgDaSharedSubjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', ...
        'sub-wlsubj124', 'sub-0426', 'sub-0250'};
    % Match plot1_experimentalCond.m's polar-plot marker size exactly.
    % plot()'s MarkerSize is a linear (roughly diameter, in points)
    % measure; scatter()'s size argument is marker AREA (points^2). To
    % render the same apparent diameter d, a circular scatter marker
    % needs area = pi/4 * d^2.
    polarMarkerSize = 6 * 0.8; % must match plot1_experimentalCond.m's markerSize
    meanDotSize = (pi/4) * polarMarkerSize^2;
    axisLineWidth = 1; % matches plot1_experimentalCond.m's polar-plot axis lines
    axisLineColor = [0.25 0.25 0.25]; % matches plot1_experimentalCond.m's polar-plot axis/grid color
    % subjectwiseDiff mode only: lighter/thinner axes + zero line than the
    % pairwise plot's axisLineColor/axisLineWidth above.
    subjDiffAxisColor = [186, 188, 190] / 255; % RGB(186,188,190)
    subjDiffAxisLineWidth = 1;
    errorbarLineWidth = styleInfo.errorbar_lineWidth; % same as plot1_experimentalCond.m's polar-plot error bars
    proLineWidth = styleInfo.pro_lineWidth; % marker outline thickness, matches polar plot's markers
    % NOT styleInfo.con_lineWidth: used to give the con dot a thinner
    % (0.5 vs 1) edge stroke than pro, which visually read as smaller
    % despite identical SizeData -- now using proLineWidth for both dots.
    showTitleLegend = false;

    % Narrow xlim (a shift, cropping in on the existing x-range) rather
    % than shrinking the axis Position width directly, so the physical
    % cm-per-data-unit scale stays exactly what it was -- i.e. the plotted
    % points/spacing are not compressed, only how much surrounding margin
    % is shown changes. previousXlim/previousPairwisePlotWidth_cm
    % describe the last (already-tuned) width=3.525cm-at-xlim=[0.5,2.5]
    % configuration; solving for the xlim that yields exactly 3cm at that
    % same scale, centered on the same midpoint.
    previousXlim = [0.5, 2.5];
    previousPairwisePlotWidth_cm = 3 * 1.25 * 0.94; % ~3.525cm
    targetPairwisePlotWidth_cm = 3;
    xScale_cm_per_unit = previousPairwisePlotWidth_cm / diff(previousXlim);
    newXlimRange = targetPairwisePlotWidth_cm / xScale_cm_per_unit;
    xlimCenter = mean(previousXlim);
    pairwiseXlim = xlimCenter + [-newXlimRange/2, newXlimRange/2];

    % retrieve the indices for specific asymmetries (e.g., motion -
    % orientation for main cardinal v main oblique)
    if ~derivedVals
        [proConditions, conConditions, ~] = retrieveProConIdx(projectName, comparisonName, subset);
    else
        proConditions = 1; conConditions = 2;
    end

    % Accumulates the exact numbers underlying this figure (one row per
    % ROI) across the loop below, so a full top-to-bottom run leaves a
    % persistent, inspectable record of the plotted values -- not just
    % console output that scrolls away, and not something you have to
    % re-run (with its own fresh bootstrap draws) to see again.
    statsRows = {};
    subjDiffRows = {}; % subjectwiseDiff mode only: one row per (roi, subject)

    % subjectwiseDiff mode only: one fixed y-range shared by every
    % asymmetry and both projects, so all subjectwiseDiff panels are
    % directly visually comparable to one another.
    diffYlim = [-0.8, 0.4];

    figure

    % Subplot grid dimensions: 1x1 for a single cortical area (matches
    % today's production V1-only configuration exactly), else the same
    % fixed 2x5 grid the ROI loop below has always used (reserved for up
    % to 10 cortical areas; unused cells simply stay empty for fewer).
    % Single source of truth for both the subplot() calls below and the
    % final figure/axes sizing, so they can't drift out of sync.
    if length(rois) == 1
        nCols = 1; nRows = 1;
    else
        nCols = 5; nRows = 2;
    end
    axHandles = gobjects(length(rois), 1);

    for ii=1:length(rois)
        rois{ii}
    
        % Specify the region index
        regionIndex = projectSettings.roi_idx{ii};
        roiname = rois{ii};

        % fitLabel: which regressionResults/<label>/ subfolder to read --
        % defaults to projectName, but for dg this must match whichever
        % dgSubjectMode (all=13 subj -> 'dg', matched=7 subj ->
        % 'dgMatched7') fitAsymmetryRegression.m was run with, since
        % projectSettings.subjects/observerGain (used for the dots above)
        % and this cached fit (used for the error bars) must describe the
        % SAME subject set. Set projectSettings.fitLabel explicitly
        % whenever using dg 'matched' mode; defaults to projectName
        % (i.e. dg 'all' mode) otherwise.
        if isfield(projectSettings, 'fitLabel') && ~isempty(projectSettings.fitLabel)
            fitLabel = projectSettings.fitLabel;
        else
            fitLabel = projectName;
        end
        fitFile = fullfile(bidsDir,'derivatives','summaryTables','regressionResults',fitLabel,sprintf('%s.mat',roiname));
        if ~isfile(fitFile)
            error('plot2_experimentalCond:missingFit', ...
                'Cached fit not found for %s / %s at %s -- run fitAsymmetryRegression(''%s'') first.', ...
                fitLabel, roiname, fitFile, projectName);
        end
        F = load(fitFile);

        % Model-derived dots (modelPro/modelCon): the group-level
        % regression's own predicted BOLD level for whichever direction(s)
        % this asymmetry's pro/con sides represent, i.e.
        %   grandInterceptFE + sum_k beta_k * predictor_k(direction,location),
        % averaged over locations (equal weight per location, matching
        % vals_1_overall/vals_2_overall's own convention below) and, for
        % the non-derived case, over multiple conditions when
        % proConditions/conConditions span more than one direction. This
        % is NOT the naive grandInterceptFE +/- estimates(termIdx)/2
        % shortcut -- that only works for mainCardinal/derivedCardinal;
        % mainSubset/derivedSubset need every predictor's own contribution
        % summed in, since their displayed direction(s) don't average that
        % contribution away to zero on their own (verified numerically
        % this session -- see AnalysisCode session notes).
        %
        % fittedByDirLoc(mi,li): the model's fitted value at direction
        % mdirvals_dg(mi), location anglevals_model(li) -- same predictor
        % formulas as fitAsymmetryRegression.m, same location order as
        % this file's own anglevals (defined below) so location index li
        % lines up row-for-row with medianBOLDpa's own location dimension.
        mdirvals_dg = [0, 90, 45, 135];
        anglevals_model = [90, 45, 0, 315, 270, 225, 180, 135];
        maincardinalmDir = [0,90,180,270];
        primaryMeridians = [90,0,270,180];
        beta1 = F.estimates(1)/2; beta2 = F.estimates(2)/2; beta3 = F.estimates(3)/2; beta4 = F.estimates(4)/2;
        fittedByDirLoc = nan(4,8);
        for mi_ = 1:4
            md_ = mdirvals_dg(mi_);
            for li_ = 1:8
                pa_ = anglevals_model(li_);
                mainCardinal_ = 2*ismember(md_, maincardinalmDir) - 1;
                derivedCardinal_ = 2*((ismember(md_,maincardinalmDir) & ismember(pa_,primaryMeridians)) | ...
                    (~ismember(md_,maincardinalmDir) & ~ismember(pa_,primaryMeridians))) - 1;
                mainSubset_ = 0; derivedSubset_ = 0;
                if strcmp(projectName, 'dg')
                    if ismember(md_,[0,180]); mainSubset_=1; elseif ismember(md_,[90,270]); mainSubset_=-1; end
                    if (abs(md_-pa_)==0 || abs(md_-pa_)==180); derivedSubset_=1;
                    elseif (abs(md_-pa_)==90 || abs(md_-pa_)==270); derivedSubset_=-1; end
                else
                    if ismember(md_,[90,270]); mainSubset_=1; elseif ismember(md_,[0,180]); mainSubset_=-1; end
                    isCardMd_ = ismember(md_,[0,90]); isOblMd_ = ismember(md_,[45,135]);
                    diffv_ = abs(md_-pa_);
                    proH_ = (isCardMd_ && (diffv_==90||diffv_==270)) || (isOblMd_ && (diffv_==0||diffv_==180));
                    conV_ = (isCardMd_ && (diffv_==0||diffv_==180)) || (isOblMd_ && (diffv_==90||diffv_==270));
                    if proH_; derivedSubset_=1; elseif conV_; derivedSubset_=-1; end
                end
                fittedByDirLoc(mi_,li_) = F.grandInterceptFE + beta1*mainCardinal_ + beta2*derivedCardinal_ + beta3*mainSubset_ + beta4*derivedSubset_;
            end
        end

        if ~derivedVals
            % Non-derived: proConditions/conConditions index directly into
            % contrasts_dict, each mapping to one physical direction (same
            % regex extraction compute_derivativeDirections.m uses) --
            % average fittedByDirLoc over locations for that direction,
            % then over conditions if more than one (e.g. mainCardinal's
            % pro = {s0,s90} averaged together).
            modelPro = mean(arrayfun(@(c) mean(fittedByDirLoc(mdirvals_dg==directionOfCondition(c, projectSettings), :)), proConditions));
            modelCon = mean(arrayfun(@(c) mean(fittedByDirLoc(mdirvals_dg==directionOfCondition(c, projectSettings), :)), conConditions));
        else
            % Derived: mirrors compute_derivativeDirections.m's own
            % per-location radial(0/180)/tangential(90/270)/oblique
            % classification -- applied here to the model's fitted values
            % instead of raw data -- then combined exactly the way
            % plot_NeuralAsymmetries.m's n_derivedConditions does:
            % subset=0 (derived cardinal) averages radial+tangential vs
            % oblique; subset=1 (derived subset/radial-vs-tangential or
            % vertical-vs-horizontal) compares radial vs tangential
            % directly. subset alone determines which recipe applies (see
            % is_mainsubset in plot_NeuralAsymmetries.m), so no extra
            % input is needed from the caller to reconstruct it here.
            radialByLoc = nan(1,8); tangByLoc = nan(1,8); obliqueByLoc = nan(1,8);
            for li_ = 1:8
                pa_ = anglevals_model(li_);
                rvals_ = []; tvals_ = []; ovals_ = [];
                for mi_ = 1:4
                    md_ = mdirvals_dg(mi_);
                    dist_ = abs(md_ - pa_);
                    if dist_==0 || dist_==180
                        rvals_ = [rvals_, fittedByDirLoc(mi_,li_)]; %#ok<AGROW>
                    elseif dist_==90 || dist_==270
                        tvals_ = [tvals_, fittedByDirLoc(mi_,li_)]; %#ok<AGROW>
                    else
                        ovals_ = [ovals_, fittedByDirLoc(mi_,li_)]; %#ok<AGROW>
                    end
                end
                radialByLoc(li_) = mean(rvals_); tangByLoc(li_) = mean(tvals_); obliqueByLoc(li_) = mean(ovals_);
            end
            if subset == 0
                modelPro = mean([radialByLoc, tangByLoc]);
                modelCon = mean(obliqueByLoc);
            else
                modelPro = mean(radialByLoc);
                modelCon = mean(tangByLoc);
            end
        end

        % Extract the relevant conditions for the specified region
        conditions1 = medianBOLDpa(proConditions, :, regionIndex, :);
        conditions2 = medianBOLDpa(conConditions, :, regionIndex, :);
        conditions1 = squeeze(conditions1);
        conditions2 = squeeze(conditions2);
        
        % Average across the conditions within subjects
        avgConditions1 = nanmean(conditions1, 1);
        avgConditions2 = nanmean(conditions2, 1);
        avgConditions1 = squeeze(avgConditions1);
        avgConditions2 = squeeze(avgConditions2);

        % Gain-weight each observer: divide their value (irrespective of
        % polar angle) by their own mean pRF gain (prfvista_mov/prfvista
        % average), BEFORE any averaging across observers. This down-weights
        % high-gain observers and up-weights low-gain observers; the
        % displayed statistics, plotted per-subject lines/scatter points,
        % and bootstrapped error bars below all inherit the adjustment
        % automatically since they are computed from these arrays. See
        % projectSettings.gainWeightsSource / retrieveObserverGainWeights2.m
        %
        % The across-observer average gain is then multiplied back in, so
        % the plotted scale/units resemble the original (unweighted) data.
        % Dividing by gain_i and then multiplying by groupGain is the same
        % as dividing by gain_i normalized to the group mean (gain_i /
        % groupGain), so the relative weighting -- and therefore the
        % relative pattern across per-subject lines, group markers, and
        % error bars -- is unchanged; only the overall scale shifts.
        %
        % groupGain uses the geometric, not arithmetic, mean: the applied
        % factor (groupGain ./ gainWeights) is itself a multiplicative
        % scale factor, and geometric mean is the choice under which the
        % *geometric* mean of that factor across observers is exactly 1
        % (magnitude-neutral in the multiplicative sense, matching what
        % the factor actually is) -- arithmetic mean does not have this
        % property. Implemented via exp(mean(log(.))) to avoid a
        % dependency on the Statistics and Machine Learning Toolbox.
        gainWeights = retrieveObserverGainWeights2(projectSettings.subjects, roiname, projectSettings.gainWeightsSource);
        if any(gainWeights <= 0)
            error('gainWeights must be strictly positive to take log() for the geometric mean (found %d non-positive value(s))', ...
                sum(gainWeights <= 0));
        end
        groupGain = exp(mean(log(gainWeights), 'omitnan')); % omitnan: see retrieveObserverGainWeights2.m
        avgConditions1 = avgConditions1 .* (groupGain ./ gainWeights);
        avgConditions2 = avgConditions2 .* (groupGain ./ gainWeights);

        % Precision-weight each observer for THIS cortical area --
        % ROI-specific, via retrieveObserverPrecisionWeights.m (see
        % plot1_experimentalCond.m for the full explanation of why this is
        % a different operation from gain correction above and must stay
        % separate from it). projectSettings.precisionWeightsSource is
        % currently [] everywhere this is called (PLACEHOLDER: every
        % (subject, roi) gets weight 1, a no-op).
        if isfield(projectSettings, 'precisionWeightsSource')
            precisionSource = projectSettings.precisionWeightsSource;
        else
            precisionSource = [];
        end
        precisionWeights = retrieveObserverPrecisionWeights(projectSettings.subjects, roiname, precisionSource);

        % Already averaged across the conditions within subjects < -- already did this in
        % the loop above

        % vals_1/vals_2 are each individual SUBJECT's own value (used for
        % the per-subject grey lines) -- not an aggregate across subjects,
        % so precision weighting (a statement about how much to trust one
        % subject RELATIVE to others) doesn't apply here; unweighted, as
        % before.
        vals_1 = nanmean(avgConditions1,1)';
        vals_2 = nanmean(avgConditions2,1)';

        % vals_1_overall/vals_2_overall (the group dots) DO combine across
        % subjects, so they get precision-weighted -- weighted mean across
        % subjects per location (weightedNanMean, at the end of this
        % file), then a plain mean across locations so each polar-angle
        % location still counts equally (this figure's established
        % "equally weighs each PA" design, unaffected by subject
        % precision-weighting). Identical to the previous
        % nanmean(...,'all') when precisionWeights is uniform, as it
        % currently is.
        vals_1_overall = mean(weightedNanMean(avgConditions1, precisionWeights), 'omitnan');
        vals_2_overall = mean(weightedNanMean(avgConditions2, precisionWeights), 'omitnan');
        %sem1 = nanstd(avgConditions1,0,2)' ./ sqrt(sum(~isnan(avgConditions1),2)');
        %sem2 = nanstd(avgConditions2,0,2)' ./ sqrt(sum(~isnan(avgConditions2),2)');
        
        % Plot the data on a polar plot
        subplot(nRows, nCols, ii)

        % Statistics come from the cached fitAsymmetryRegression.m output
        % (joint, gain- and precision-weighted regression, 1000
        % subject-level bootstrap draws already computed there) regardless
        % of displayMode/plotType -- there's no separate "raw" bootstrap
        % alternative in this pipeline.
        fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        bootDraws = F.coeffs(termIdx, :)';
        meanDiff = F.estimates(termIdx);
        ci_mean = prctile(bootDraws, [2.5 97.5]);

        % 68% CI: used as the pairwise plot's error-bar width (half-width
        % applied symmetrically to both vals_1 and vals_2, since this CI
        % describes the (pro-con) difference, not either condition's own
        % uncertainty separately, so there's no unique way to split it
        % between the two markers).
        ci_mean_68 = prctile(bootDraws, [16 84]);
        ci68_halfwidth = (ci_mean_68(2) - ci_mean_68(1)) / 2;

        % Significance, independent of plotType (and of the sign flip
        % below -- a CI excluding 0 still excludes 0 under negation): one
        % asterisk if the 68% CI excludes 0 (but the 95% CI does not), two
        % if the 95% CI excludes 0 (which always implies the narrower 68%
        % CI also does, since it's nested inside it) -- standard
        % convention, more asterisks = stronger evidence, mutually
        % exclusive tiers.
        sig95 = ci_mean(1) > 0 || ci_mean(2) < 0;
        sig68 = ci_mean_68(1) > 0 || ci_mean_68(2) < 0;
        if sig95
            sigStr = '**';
        elseif sig68
            sigStr = '*';
        else
            sigStr = '';
        end

        if strcmp(plotType, 'pairwise')
            % dotPro/dotCon: which of the two always-computed quantities
            % (modelPro/modelCon above, or vals_1_overall/vals_2_overall
            % below) actually gets drawn is the ONLY thing displayMode
            % controls for the dots -- both are computed unconditionally
            % every time so they stay directly comparable run to run.
            % centerVal is the single shared anchor both the dots AND the
            % grey lines use for this mode: mean([modelPro,modelCon]) for
            % 'model', or mean([vals_1_overall,vals_2_overall]) for 'raw'
            % -- the SAME precision-weighted group average the dots
            % themselves use, not a separately-computed unweighted grand
            % mean (that was the pre-toggle behavior; see displayMode's
            % own comment above for why it changed).
            %
            % NOTE: centerVal is NOT F.grandInterceptFE for 'model' --
            % that was a real bug caught by visual inspection
            % (mainSubset/derivedSubset's lines came out visibly offset
            % from the dots). mean([modelPro,modelCon]) only reduces to
            % grandInterceptFE for the two asymmetries whose naive
            % single-term formula is exact (mainCardinal, derivedCardinal);
            % mainSubset/derivedSubset's full formula carries an extra
            % shared term (e.g. both modelPro and modelCon include
            % +beta(mainCardinal)) that does NOT cancel out of the
            % average, so grandInterceptFE alone is the wrong anchor for
            % them. Using mean([modelPro,modelCon]) directly is correct
            % for all four asymmetries unconditionally, with no
            % special-casing needed.
            if strcmp(displayMode, 'model')
                dotPro = modelPro; dotCon = modelCon;
                centerVal = mean([modelPro, modelCon]);
            else
                dotPro = vals_1_overall; dotCon = vals_2_overall;
                centerVal = mean([vals_1_overall, vals_2_overall]);
            end

            % Grey line SLOPE, mode-dependent:
            %   'model': this subject's contribution to the group-level
            %   asymmetry coefficient (F.subjectContributions), RESCALED
            %   by F.nSubj. subjectContributions is built so that SUMMING
            %   all subjects reproduces F.estimates -- but F.estimates is
            %   a mean (a regression coefficient), not a sum, over
            %   subjects, so each individual raw contribution is already
            %   scaled down by ~1/nSubj relative to that subject's own
            %   marginal difference. Multiplying back by nSubj undoes that
            %   averaging and puts the slope back on the same scale as the
            %   group dots / raw data, while still reducing to exactly
            %   that subject's own raw marginal (vals_1-vals_2) difference
            %   for a subject with complete data (see the note in
            %   fitAsymmetryRegression.m) -- verified numerically for V1.
            %   For cortical areas/subjects with missing locations, the
            %   raw marginal difference can be confounded across the 4
            %   asymmetries in a way this model-based contribution
            %   corrects for.
            %   'raw': that subject's own raw (vals_1(s)-vals_2(s))
            %   difference -- no model involved at all.
            % Both slopes are centered on centerVal (this mode's shared
            % anchor) -- verified identical to each other for any
            % complete-data subject (e.g. every subject at V1), since
            % that's exactly the condition under which a subject's
            % rescaled contribution collapses to their own raw difference.
            %
            % A subject with genuinely zero data for this cortical area
            % has vals_1(subjectIndex)/vals_2(subjectIndex) = NaN
            % (propagated from the raw extraction above), so plot() below
            % silently draws nothing for them in 'raw' mode -- same
            % behavior as before this change. 'model' mode's slope only
            % relies on subjectContributions (0 for a subject with no data
            % in this ROI, not NaN -- see fitAsymmetryRegression.m), so it
            % draws a flat (zero-slope) line at centerVal for such a
            % subject instead of omitting them; this is a real, minor
            % behavioral difference between the two modes, worth knowing
            % if the two are ever compared side by side for an
            % incomplete-data cortical area.
            for subjectIndex = 1:size(medianBOLDpa, 4)
                if strcmp(displayMode, 'model')
                    slope = F.nSubj * F.subjectContributions(subjectIndex, termIdx);
                else
                    slope = vals_1(subjectIndex) - vals_2(subjectIndex);
                end
                linePos = centerVal + [slope/2, -slope/2];

                plot([1 2], linePos, 'Color', subjectLineColor);
                xlim(pairwiseXlim)
                hold on
            end
            hold on

            % record this ROI's aggregate stats for the export below --
            % both the raw and model dot values are always saved
            % (regardless of displayMode) so a diff against the cortical
            % area's own completeness (do rawPro/rawCon match
            % modelPro/modelCon?) is always available without re-running
            % under the other mode.
            statsRows(end+1,:) = {rois{ii}, meanDiff, ci_mean_68(1), ci_mean_68(2), ci_mean(1), ci_mean(2), ...
                vals_1_overall, vals_2_overall, modelPro, modelCon}; %#ok<SAGROW>

            %% Print to console
            fprintf('Mean difference: %.4f, 95%% CI: [%.4f, %.4f]\n', meanDiff, ci_mean(1), ci_mean(2));
            fprintf('rawPro=%.6f modelPro=%.6f (Delta=%.6f) | rawCon=%.6f modelCon=%.6f (Delta=%.6f)\n', ...
                vals_1_overall, modelPro, vals_1_overall-modelPro, vals_2_overall, modelCon, vals_2_overall-modelCon);
            fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

            % Mean dots drawn first. meanDotSize (SizeData, i.e. marker
            % area) is identical for both -- but COLORS.json sets
            % con_lineWidth (0.5) thinner than pro_lineWidth (1) for every
            % asymmetry, which makes the unfilled/con dot's edge stroke
            % visibly thinner and the dot read as smaller even though its
            % underlying size is the same. Using proLineWidth for both
            % here (pairwise dots only -- COLORS.json itself is shared
            % with other plots, e.g. plot1's polar markers, so left
            % unchanged) so the two dots read as the same size.
            scatter(1, dotPro, meanDotSize, 'MarkerFaceColor', proFaceColor, 'MarkerEdgeColor', proEdgeColor, 'LineWidth',proLineWidth);
            hold on
            scatter(2, dotCon,  meanDotSize, 'MarkerFaceColor', conFaceColor, 'MarkerEdgeColor', conEdgeColor, 'LineWidth',proLineWidth);
            hold on

            % Error bars -- top layer, drawn last -- matches
            % plot1_experimentalCond.m's polar-plot draw order. Centered
            % on dotPro/dotCon (whichever mode is displayed) rather than a
            % third, separately-computed quantity, so the error bar always
            % stays visually attached to the dot it belongs to regardless
            % of mode. Width (ci68_halfwidth) always comes from the
            % model's bootstrap either way -- that was already a
            % deliberate prior decision (see meanDiff/bootDraws above),
            % unrelated to which mode the dots/lines are in, and unchanged
            % by this toggle.
            errorbar(1, dotPro, ci68_halfwidth, 'Color', proEdgeColor, 'LineWidth', errorbarLineWidth, 'CapSize', 0);
            hold on
            errorbar(2, dotCon, ci68_halfwidth, 'Color', conEdgeColor, 'LineWidth', errorbarLineWidth, 'CapSize', 0);

            if ~isempty(sigStr)
                text(1.5, 0.6, sigStr, 'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', 'FontSize', 10, 'Color', axisLineColor);
            end

        else % 'subjectwiseDiff'
            % Per-subject: raw (pro-con) difference, no intercept/model
            % subtraction -- exactly this subject's own value in one
            % condition minus the other, already averaged (above, in
            % vals_1/vals_2) across all 8 locations equally. No sign
            % correction: pro-con follows this asymmetry's own convention
            % (see retrieveProConIdx.m/compute_derivativeDirections.m)
            % unchanged -- normalizing to a canonical "X minus Y" order is
            % only relevant for the cross-condition (dg vs da) comparison
            % script, not this plot.
            subjDiff = vals_1 - vals_2;

            % Group-level: same cached bootstrap the pairwise error bars
            % use, but the 95% CI (not 68%) -- this is now a test of
            % whether the difference itself is nonzero, not a pro-vs-con
            % comparison, so the narrower band used there for visually
            % separating two dots doesn't apply.
            ci95_halfwidth = (ci_mean(2) - ci_mean(1)) / 2;

            for si_ = 1:numel(subjDiff)
                subjDiffRows(end+1,:) = {rois{ii}, projectSettings.subjects{si_}, subjDiff(si_)}; %#ok<SAGROW>
            end
            statsRows(end+1,:) = {rois{ii}, meanDiff, ci_mean(1), ci_mean(2)}; %#ok<SAGROW>

            fprintf('Mean difference: %.4f, 95%% CI: [%.4f, %.4f]\n', meanDiff, ci_mean(1), ci_mean(2));
            fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

            % Dashed zero line, drawn under the data. A plain plot() line,
            % not yline() -- ConstantLine objects (what yline() creates)
            % are often rendered on top of other axes children regardless
            % of call order, so drawing this first would not reliably put
            % it behind the dots. A regular Line object obeys standard
            % MATLAB z-stacking (later calls draw on top), so drawing it
            % first here is sufficient.
            %
            % Manually-constructed dash pattern (NaN-separated segments),
            % same technique plot1_experimentalCond.m's polar con curves
            % use -- native '--' LineStyle exposes no control over
            % dash-length vs. gap independently. dashLen/gapLen are in
            % x-data units, at 2x the length/spacing of the polar plot's
            % own dashDeg=3/gapDeg=1.5 (same 2:1 ratio).
            dashLen = 0.12; gapLen = 0.06;
            unitLen = dashLen + gapLen;
            zeroLineX = [];
            for u = 0:ceil(diff(pairwiseXlim)/unitLen)-1
                segStart = pairwiseXlim(1) + u*unitLen;
                if segStart >= pairwiseXlim(2), break; end
                segEnd = min(segStart+dashLen, pairwiseXlim(2));
                zeroLineX = [zeroLineX, segStart, segEnd, NaN]; %#ok<AGROW>
            end
            plot(zeroLineX, zeros(size(zeroLineX)), '-', 'LineWidth', subjDiffAxisLineWidth, 'Color', subjDiffAxisColor);
            hold on

            % Per-subject dots at x=1, x-jittered, no connecting lines (a
            % single condition's worth of values per subject here, not a
            % pro/con pair, so there is nothing meaningful to connect).
            % All observer circles are rendered IDENTICALLY regardless of
            % sign now -- same face color (this asymmetry's own color),
            % same white outline, same size, same LineWidth. Sign is
            % conveyed only by y-position; there is no more filled-vs-
            % unfilled (or any other) visual distinction between a
            % positive and a negative dot. The white outline is invisible
            % against the page background and only becomes visible where
            % two dots overlap (a separating halo). SizeData is grown so
            % the FACE alone now spans what the OLD dot's total visual
            % footprint (face + stroke, at the OLD full LineWidth)
            % occupied. Subjects who only participated in the DG
            % experiment (absent from da's roster -- only possible here
            % for dg's 'all' 13-subject mode) are drawn as 'x' markers
            % instead of circles, unaffected by any of this -- sized 1.5x
            % the (now-grown) circle face size, with the asymmetry's own
            % colored edge at the full LineWidth, since the dg-vs-da
            % comparison this figure set supports doesn't apply to them.
            nSubjPlot_ = numel(subjDiff);
            jitterWidth = 0.45; % 1.5x the original 0.3 (doubling to 0.6 was too much)
            xJitter = 1 + (rand(nSubjPlot_,1) - 0.5) * jitterWidth;
            isDgOnly_ = reshape(strcmp(projectName, 'dg') & ...
                ~ismember(projectSettings.subjects, dgDaSharedSubjects), [], 1);
            isCircle_ = ~isDgOnly_;

            circleSizeOld_ = meanDotSize*0.4;
            oldCircleDiameterPts_ = sqrt(4*circleSizeOld_/pi) + proLineWidth;
            circleFaceSize_ = (pi/4) * oldCircleDiameterPts_^2;
            whiteOutlineLineWidth_ = proLineWidth/4; % half of the previous half-width outline
            xMarkerSize_ = circleFaceSize_*1.5; % 1.5x the CURRENT (grown) circle face size
            scatter(xJitter(isCircle_), subjDiff(isCircle_), circleFaceSize_, ...
                'Marker', 'o', 'MarkerFaceColor', proFaceColor, 'MarkerEdgeColor', 'w', 'LineWidth', whiteOutlineLineWidth_);
            hold on
            scatter(xJitter(isDgOnly_), subjDiff(isDgOnly_), xMarkerSize_, ...
                'Marker', 'x', 'MarkerEdgeColor', proEdgeColor, 'LineWidth', proLineWidth);
            hold on
            xlim(pairwiseXlim)

            % Group mean + 95% CI at x=2, same marker style the pairwise
            % plot's dots use -- always filled (solid) regardless of
            % asymmetry/sign or project, and carries the same half-width
            % white outline as the observer dots above. Its face SizeData
            % is grown from meanDotSize (the OLD face+full-stroke
            % footprint) the same way the observer dots' was, so it still
            % reads as the same overall size as before despite gaining an
            % outline. The error bar is drawn FIRST so the dot (face +
            % outline) sits on top of it, rather than the CI whisker
            % cutting across and interrupting the outline.
            errorbar(2, meanDiff, ci95_halfwidth, 'Color', proEdgeColor, 'LineWidth', errorbarLineWidth, 'CapSize', 0);
            hold on
            meanFaceColor = proFaceColor;
            meanDiameterPts_ = sqrt(4*meanDotSize/pi) + proLineWidth;
            meanFaceSize_ = (pi/4) * meanDiameterPts_^2;
            scatter(2, meanDiff, meanFaceSize_, 'MarkerFaceColor', meanFaceColor, 'MarkerEdgeColor', 'w', 'LineWidth', whiteOutlineLineWidth_);
            hold on

            if ~isempty(sigStr)
                text(1.5, diffYlim(2) - 0.1, sigStr, 'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', 'FontSize', 10, 'Color', axisLineColor);
            end
        end

        if showTitleLegend
            title(rois{ii});
        end
        %ylabel('zscored PSC')
        if strcmp(plotType, 'subjectwiseDiff')
            % Fixed rotation (not MATLAB's auto-choice, which varies with
            % the y-tick label widths -- those differ by asymmetry since
            % each has its own fixed diffYlim range, and an auto-rotation
            % that varies run to run made the required page margins
            % unpredictable and prone to clipping the label text).
            set(gca, 'XTick', [1 2], 'XTickLabel', {'observer','group'}, 'XTickLabelRotation', 45);
            % Tick marks every 0.2 across the fixed [-0.8, 0.4] diffYlim,
            % but text labels only every 0.4 (blank labels at the
            % in-between ticks) to keep the axis legible at this font size.
            set(gca, 'YTick', -0.8:0.2:0.4, ...
                'YTickLabel', {'-0.8','','-0.4','','0','','0.4'}, ...
                'FontSize', 8);
        else
            set(gca, 'XTick', []);
        end
        ax = gca;
        axHandles(ii) = ax;
        ax.LineWidth = axisLineWidth;
        ax.XColor = axisLineColor; % matches plot1_experimentalCond.m's polar-plot ThetaColor/RColor
        ax.YColor = axisLineColor;
        if strcmp(plotType, 'subjectwiseDiff')
            ax.LineWidth = subjDiffAxisLineWidth;
            ax.XColor = subjDiffAxisColor;
            ax.YColor = subjDiffAxisColor;
        end
        % ax.Layer left at its default ('bottom') -- matches
        % plot1_experimentalCond.m's polar plots, where 'top' drew grid
        % lines over solid data points.
        box off

        if ismember(rois{ii}, {'pMT', 'pMST', 'hMTcomplex'})
            ROI_category = 'ROIs_motion';
        else
            ROI_category = 'ROIs_early';
        end

        if strcmp(plotType, 'pairwise')
            ylim([pairaxes_PAew_limits.(projectName).(comparisonName).(ROI_category).min ...
                         pairaxes_PAew_limits.(projectName).(comparisonName).(ROI_category).max])
        end

        if ii==1 && showTitleLegend
            if strcmp(asymmetryName, 'radialVsTangential')
                lg1 = legend('Radial', 'Tangential', 'Location', 'northeast');
            elseif strcmp(asymmetryName, 'verticalVsHorizontal')
                lg1 = legend('Horizontal', 'Vertical', 'Location', 'northeast');
            else
                lg1 = legend('Card', 'Obl', 'Location', 'northeast');
            end
            lg1.Position = [0.7946 0.2614 0.1089 0.0631];
        end
    
        hold off;
    end
    
    if showTitleLegend
        sgtitle(sprintf('%s: %s', projectName, strrep(comparisonName, '_', ' ')), 'FontSize', 40)
    end

%     fa = gcf;
%     fa.Position = [1000 555 1514 782];

    % ROI name(s) plotted in this call, joined with '-' -- included in the
    % filename so re-running for a different cortical area (or set of
    % areas) doesn't silently overwrite another area's file. Today this is
    % always just the single ROI in projectSettings.rois (production is
    % currently scoped to rois(1), i.e. V1 only), but this works
    % unchanged if that scope is ever widened to multiple/other areas.
    roiNameStr = strjoin(rois, '-');
    if strcmp(plotType, 'pairwise')
        filename = fullfile(figureDir,sprintf('pairwise_PAequalweight_%s_%s_%s_%s_%s', comparisonName, projectName, asymmetryName, roiNameStr, displayMode));
    else
        filename = fullfile(figureDir,sprintf('subjectwiseDiff_PAequalweight_%s_%s_%s_%s', comparisonName, projectName, asymmetryName, roiNameStr));
    end

    % Save the exact values underlying this figure alongside the PDF, so
    % a full top-to-bottom run leaves a persistent, reproducible record
    % (diffable across runs) rather than requiring a re-run -- with its
    % own fresh bootstrap draws -- to see the numbers again. rawPro/rawCon
    % and modelPro/modelCon are both always saved regardless of
    % displayMode (see the statsRows note above). subjectwiseDiff mode
    % additionally saves every subject's own per-ROI difference (the x=1
    % dots), not just the group summary.
    if strcmp(plotType, 'pairwise')
        statsTable = cell2table(statsRows, 'VariableNames', ...
            {'roi','meanDiff','ci68_lower','ci68_upper','ci95_lower','ci95_upper','rawPro','rawCon','modelPro','modelCon'});
        writetable(statsTable, [filename, '_stats.csv']);
    else
        statsTable = cell2table(statsRows, 'VariableNames', ...
            {'roi','meanDiff','ci95_lower','ci95_upper'});
        writetable(statsTable, [filename, '_stats.csv']);
        subjDiffTable = cell2table(subjDiffRows, 'VariableNames', {'roi','subject','diff'});
        writetable(subjDiffTable, [filename, '_subjectDiffs.csv']);
    end

    % The axis (plot) box itself is the sized element (3.75 x 3.6 cm) --
    % the PDF page is padded larger around it (padding_cm each side) so
    % the axis box is guaranteed exactly this size regardless of how much
    % room MATLAB's auto layout would otherwise reserve for tick labels.
    % NOT fitFig2Page, which scales the figure to fill a full
    % letter-landscape page and would override any custom size set here.
    pairwisePlotWidth_cm = targetPairwisePlotWidth_cm; % 3cm, achieved via pairwiseXlim above rather than by scaling the axis box directly
    pairwisePlotHeight_cm = 3.6 * 0.94; % 6% smaller
    basePadding_cm = 0.5; % previous padding, on each side
    areaScale = 1.25; % total figure area, not just padding, should grow by this factor
    % Solve for the (still-uniform, all 4 sides) padding p that gives
    % (axisW+2p)*(axisH+2p) = areaScale * (axisW+2*basePadding)*(axisH+2*basePadding),
    % i.e. a quadratic in p: 4p^2 + 2p(axisW+axisH) + axisW*axisH - targetArea = 0
    targetArea = areaScale * (pairwisePlotWidth_cm + 2*basePadding_cm) * (pairwisePlotHeight_cm + 2*basePadding_cm);
    qa = 4; qb = 2*(pairwisePlotWidth_cm + pairwisePlotHeight_cm); qc = pairwisePlotWidth_cm*pairwisePlotHeight_cm - targetArea;
    padding_cm = (-qb + sqrt(qb^2 - 4*qa*qc)) / (2*qa);

    % Figure canvas scales with the actual subplot grid (nCols x nRows,
    % set once near the top of this function): each cortical area's axis
    % box keeps the same fixed size and per-side padding as the original
    % single-cortical-area design, just tiled across the grid instead of
    % only ever positioning one axes. For nCols=nRows=1 (today's
    % production V1-only case) this reduces EXACTLY to the original
    % formula/positions -- verified below by construction, not just by
    % eye: figWidth_cm/figHeight_cm and the single axes' Position are
    % algebraically identical to the pre-existing single-ROI-only code.
    % subjectwiseDiff mode's 'observer'/'group' XTickLabels are wider than
    % the narrow 3cm axis, so they're drawn at a fixed 45-degree rotation
    % (set above) -- that rotated text extends both below AND to the left
    % of the x=1 tick position, so extra room is needed on both the bottom
    % (or the text gets clipped by the page edge) and the left (or "observer"
    % gets clipped down to "erver"). Added only to those two sides (via the
    % axes Position loop below), not all 4, so top/right margins stay
    % exactly padding_cm as before.
    extraBottomForXTickLabels_cm = 0;
    extraLeftForXTickLabels_cm = 0;
    if strcmp(plotType, 'subjectwiseDiff')
        extraBottomForXTickLabels_cm = 1.2;
        extraLeftForXTickLabels_cm = 1.2;
    end

    figWidth_cm = nCols*pairwisePlotWidth_cm + (nCols+1)*padding_cm + extraLeftForXTickLabels_cm;
    figHeight_cm = nRows*pairwisePlotHeight_cm + (nRows+1)*padding_cm + extraBottomForXTickLabels_cm;

    gcf_edit = gcf;
    gcf_edit.Units = 'centimeters';
    gcf_edit.Position(3:4) = [figWidth_cm, figHeight_cm];
    gcf_edit.PaperUnits = 'centimeters';
    gcf_edit.PaperPositionMode = 'manual'; % otherwise 'auto' ignores PaperPosition below
    gcf_edit.PaperSize = [figWidth_cm, figHeight_cm];
    gcf_edit.PaperPosition = [0, 0, figWidth_cm, figHeight_cm];

    % Position EVERY cortical area's axes into its own grid cell (previously
    % only the last-created axes was ever repositioned here, which is
    % harmless for nCols=nRows=1 but left every other cortical area's axes at
    % its original, unrelated default-figure-sized subplot position once
    % the figure was shrunk -- causing them to visually overlap).
    for jj = 1:length(rois)
        row = floor((jj-1)/nCols) + 1; % 1 = top row, matches subplot()'s numbering
        col = mod(jj-1, nCols) + 1;
        left_cm = padding_cm + extraLeftForXTickLabels_cm + (col-1)*(pairwisePlotWidth_cm + padding_cm);
        bottom_cm = padding_cm + extraBottomForXTickLabels_cm + (nRows-row)*(pairwisePlotHeight_cm + padding_cm);
        axHandles(jj).Units = 'centimeters';
        axHandles(jj).Position = [left_cm, bottom_cm, pairwisePlotWidth_cm, pairwisePlotHeight_cm];
        if strcmp(plotType, 'pairwise')
            ylim(axHandles(jj), [-.25 0.75]) % if zero-meaning the data %%%% -- applied to every cortical area's axes, not just the last
        else
            ylim(axHandles(jj), diffYlim)
        end
    end
    % Save as PDF
    set(gcf_edit,'Renderer','painters'); % new
    print(gcf_edit, filename, '-dpdf', '-painters');
    close all;

end

function m = weightedNanMean(X, w)
% WEIGHTEDNANMEAN  Precision-weighted mean of X (nLoc x nSubj) across
% subjects (dim 2), ignoring NaNs, using per-subject weights w (nSubj x 1
% or 1 x nSubj). Returns an nLoc x 1 column. Identical to nanmean(X,2)
% when w is uniform (as it currently is, pending finalized precision
% weights). Duplicated from plot1_experimentalCond.m's identical helper
% rather than shared, matching this file's existing convention (e.g. the
% gain-correction block above is likewise duplicated, not shared).
    w = w(:)'; % row, so it broadcasts against X's columns (subjects)
    validMask = ~isnan(X);
    Wexpanded = repmat(w, size(X,1), 1);
    Wexpanded(~validMask) = 0;
    Xz = X;
    Xz(~validMask) = 0;
    m = sum(Xz .* Wexpanded, 2) ./ sum(Wexpanded, 2);
end

function d = directionOfCondition(condIdx, projectSettings)
% DIRECTIONOFCONDITION  Physical stimulus direction (0/90/45/135) encoded
% in a contrasts_dict condition's own name (e.g. 's0_v_b' -> 0), same
% regex extraction compute_derivativeDirections.m uses on
% dg_contrast_name. Used to look up that condition's row in
% fittedByDirLoc above -- only meaningful for orientation-vs-blank
% conditions (26-29), which is all this is ever called on.
    contrastValue = projectSettings.contrasts_dict.contrasts(condIdx).dg_contrast_name;
    before_v = regexp(contrastValue, '^(.*)_v_', 'tokens', 'once');
    if isempty(before_v)
        error('directionOfCondition:noMatch', 'No direction found in condition name ''%s''', contrastValue);
    end
    d = str2double(regexp(before_v{1}, '\d+', 'match', 'once'));
end
