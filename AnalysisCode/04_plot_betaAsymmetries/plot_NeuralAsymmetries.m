% plot BOLD asymmetries

clc; clear all; close all;

% set up
addpath(genpath(pwd));
projectName = 'da'; %'dots'; %'dg';
%bidsDir =  '/Volumes/server/Projects/Project_dg/data_bids/';
%bidsDir =  '/Volumes/EXTERNAL_US/Project_dg/data_bids/';
bidsDir =  '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
githubDir = '~/Documents/GitHub';
hRF_setting = 'glmsingle'; %'canonical' %'glmsingle';
fullfile(githubDir, 'DriftingGrating', 'AnalysisCode')
glmResultsfolder = fullfile(bidsDir, 'derivatives', strcat(projectName, 'GLM'), strcat('hRF_', hRF_setting));

% can be 'motion_minus_orientation' ; 'motion_minus_baseline' ; 'orientation_minus_baseline'
comparisonName = 'orientation_minus_baseline'; %'motion_minus_orientation';

projectSettings = loadConfig(githubDir);
projectSettings.glmResultsfolder = glmResultsfolder;

rois = projectSettings.rois;
axes_limits = projectSettings.axes_limits;
pairaxes_limits = projectSettings.pairaxes_limits;
pairaxes_PAew_limits = projectSettings.pairaxes_PAew_limits;
colors_data = projectSettings.colors_data;
contrasts_dict = projectSettings.contrasts_dict;

projectSettings.rois = rois(1); %:7); % remove once I process v3a / v3b

%% Load and remove subject with extreme motion

load(fullfile(glmResultsfolder, 'meanBOLDpa')) % contrasts x polarAngles x ROIs x subjects
load(fullfile(glmResultsfolder, 'meanBOLD')) % contrasts x ROIs x subjects

% dg_subjectMode is the SINGLE source of truth for which dg subject set
% this run uses -- 'all' (13 subjects) or 'matched' (the 8-subject subset
% matched to da). Everything downstream (subject list, meanBOLDpa/
% meanBOLD subsetting, projectSettings.subjects, projectSettings.
% observerGain, and figureDir) is derived from this one variable so they
% cannot drift out of sync with each other, including across partial/
% section-by-section re-runs of this script.
dg_subjectMode = 'matched'; %'all'; % 'all' or 'matched'

figureDirSuffix = '';
if strcmp(projectName, 'dg')
    figureDirSuffix = ['_', dg_subjectMode]; % e.g. figures/dg_all or figures/dg_matched
end
figureDir = [strrep(bidsDir, 'data_bids', 'figures'), projectName, figureDirSuffix];

if ~isfolder(figureDir)
    mkdir(figureDir)
end

% sub-0395 excluded from both dg's 'matched' mode and da: this observer's
% da session used a pilot stimulus not matched to the other observers, so
% their data isn't comparable in either the matched-subset dg analysis or
% the da analysis itself. (dg's 'all' 13-subject mode is unaffected --
% that comparison doesn't depend on da matching.) Excluding by name
% (rather than hand-editing the numeric index) so dg and da can't
% silently drift out of sync with each other.
excludedSubjects = {'sub-0395'};

if strcmp(projectName, 'dg')
    dg_subjects_8 = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
        'sub-0395', 'sub-0426', 'sub-0250'};
    dg_subjects_13 = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
        'sub-0395', 'sub-0426', 'sub-0250', 'sub-0442', 'sub-wlsubj121', 'sub-wlsubj127', ...
        'sub-0397', 'sub-0427'};

    if strcmp(dg_subjectMode, 'all')
        subjects = dg_subjects_13;
        % meanBOLDpa/meanBOLD already contain all 13 dg subjects, in this
        % exact order (see 03_process_groupBetas/meanWithinLabel.m) -- no
        % subsetting needed
    elseif strcmp(dg_subjectMode, 'matched')
        keepIdx = find(~ismember(dg_subjects_8, excludedSubjects));
        subjects = dg_subjects_8(keepIdx);
        meanBOLDpa = meanBOLDpa(:,:,:,keepIdx); % trying to match subjects across
        meanBOLD = meanBOLD(:,:,keepIdx);
    else
        error('dg_subjectMode must be ''all'' or ''matched''.');
    end
elseif strcmp(projectName, 'da')
    % load subjects -- da's raw meanBOLDpa/meanBOLD already contain
    % exactly these 8 subjects, in this order, so subsetting is needed
    % here (unlike dg's 'all' mode) to drop the excluded one
    da_subjects_8 = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
        'sub-0395', 'sub-0426', 'sub-0250'};
    keepIdx = find(~ismember(da_subjects_8, excludedSubjects));
    subjects = da_subjects_8(keepIdx);
    meanBOLDpa = meanBOLDpa(:,:,:,keepIdx);
    meanBOLD = meanBOLD(:,:,keepIdx);
elseif strcmp(projectName, 'dots')
    subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
        'sub-0426'};
end

% contrastnames = {
% 1- cardMsep
% 2- oblMsep
% 3 - allmValls
% 4 - allsVblank
% 5 - allmVblank
% 6 - cardmVblank
% 7 - oblmVblank
% 8 - m0_v_s90
% 9 - m90_v_s0
% 10 - m180_v_s90
% 11 - m270_v_s0
% 12 - m45_v_s135
% 13 - m135_v_s45
% 14 - m225_v_s135
% 15 - m315_v_s45
% 16 - cardsVblank
% 17 - oblsVblank
% 18 - m0_v_b
% 19 - m180_v_b 
% 20 - m90_v_b 
% 21 - m270_v_b 
% 22 - m45_v_b 
% 23 - m225_v_b 
% 24 - m135_v_b
% 25 - m315_v_b 
% 26 - s0_v_b 
% 27 - s90_v_b 
% 28 - s45_v_b 
% 29 - s135_v_b

% add to projectSettings
projectSettings.projectName = projectName;
projectSettings.comparisonName = comparisonName;
projectSettings.subjects = subjects;
projectSettings.figureDir = figureDir;
if strcmp(projectName, 'dg')
    projectSettings.dgSubjectMode = dg_subjectMode; % 'all' or 'matched' -- traceable alongside subjects/figureDir above
end

% per-observer, per-ROI pRF gain, used to gain-weight each observer's
% contribution in plot1_experimentalCond.m and plot2_experimentalCond.m --
% see retrieveObserverGainWeights2.m / computeObserverGainWeightsByROI.m.
% gainSummaryFile is kept (still the old V1-only file) purely because
% several downstream scripts derive bidsDir from its directory via
% fileparts(); it is no longer used to compute gain values directly.
% gainWeightsSource is the ROI-aware table (subject/roi/weight columns);
% unlike precisionWeightsSource there's no [] placeholder -- gain
% correction is not optional, so each script looks this table up per ROI
% inside its own loop, same pattern as precision.
projectSettings.gainSummaryFile = fullfile(bidsDir, 'derivatives', 'summaryTables', 'gainSummary.mat');
Ggain = load(fullfile(bidsDir, 'derivatives', 'summaryTables', 'gainSummaryByROI.mat'), 'gainTable');
projectSettings.gainWeightsSource = Ggain.gainTable;

% Precision weight source -- ROI-specific (reliability genuinely varies
% by cortical area, unlike gain), looked up per (subject, roi) via
% retrieveObserverPrecisionWeights.m in plot1_experimentalCond.m and
% plot2_experimentalCond.m. [] is the PLACEHOLDER (uniform/no-op, every
% (subject, roi) gets weight 1) until the precision-weighting method is
% finalized -- swap in a table with columns subject/roi/weight once it
% is; no other change needed in the scripts that read this.
projectSettings.precisionWeightsSource = [];


%% MAIN CONDITION: Plot pairwise plots (JUST FOR SANITY CHECK)
% e.g. condition indices 4 v 5: clearly shows MT as motion responsive
% these do not depend on polar angle / retinotopy (just contrasts)

% condIdx1 = [27]; %[26,27]; %[9,11];   
% condIdx2 = [26]; %[28,29]; %[8, 10];
% 
% plot0_experimentalCond(condIdx1, condIdx2, medianBOLD, projectSettings)

%% MAIN CONDITION: Plot polor plots

% % this will plot the main conditions for:
% if strcmp(projectName, 'da') || strcmp(projectName, 'dots') % experiment 2: da (polar cardinal vs oblique; and radial vs tangential)
%     isradial = [0, 1]; 
% elseif strcmp(projectName, 'dg') % experiment 1: dg (cartesian cardinal vs oblique)
%     isradial = 0; 
% end

% just controls the loop to run 2x for either experiment (dg, da)
is_mainsubset = [0, 1]; % if 0 do cardinal for dg or polar cardinal for da
                        % if 1 do vertical for dg or radial for da

counter = 0;

for ci=1:numel(is_mainsubset) %isradial)
    
    %radialvstang = isradial(ci);
    subset = is_mainsubset(ci);
    % this will plot the main conditions for experiment 1: dg (cardinal vs oblique) 
    % and for experiment 2: da (polar cardinal vs oblique; and radial vs
    % tangential)
    
    % this depends on retinotopy to isolate each direction with respect to its
    % location, but still only plots the main conditions
    plot1_experimentalCond(meanBOLDpa, 'mainCardinalVsMainOblique', projectSettings, subset) %radialvstang)

    % plot each absolute directions (only needed 1x)
    if strcmp(projectName, 'dg') && subset==0
        plotSepDirs(meanBOLDpa, 'mainCardinalVsMainOblique', 'Abs', projectSettings, subset) %radialvstang)
    elseif (strcmp(projectName, 'da') || strcmp(projectName, 'dots')) && subset==0
        plotSepDirs(meanBOLDpa, 'mainCardinalVsMainOblique', 'Rel', projectSettings, subset) %radialvstang)
    end

    % MAIN CONDITION: Plot pairwise plots (EQUALLY WEIGHING POLAR ANGLE)
    % this section will average these values to equally weigh each of the 8
    % visual field locations
    
    plot2_experimentalCond(meanBOLDpa, 'mainCardinalVsMainOblique', projectSettings, subset) %radialvstang)

    % same asymmetry, alternate rendering: per-subject raw differences
    % (x=1, jittered) + group mean/95% CI (x=2) instead of two dots/grey
    % lines -- see plot2_experimentalCond.m's plotType argument.
    plot2_experimentalCond(meanBOLDpa, 'mainCardinalVsMainOblique', projectSettings, subset, 'model', 'subjectwiseDiff')

    counter = counter+1;
end



%% DERIVED CONDITION
% (THESE DEPEND ON THE POLAR ANGLE WEDGE SIZE)

% % this will plot the derived conditions for:
% if strcmp(projectName, 'da') || strcmp(projectName, 'dots') % experiment 2: da (cartesian cardinal vs oblique)
%     n_derivedConditions = {{1:2, 3}};
%     %isradial = 0; 
% elseif strcmp(projectName, 'dg') % experiment 1: dg (polar cardinal vs oblique; and radial vs tangential)
%     n_derivedConditions = {{1:2, 3}, {1, 2}};
%     %isradial = [0, 1]; 
% end

% first pair are derived cardinals, second pair are derived subset
% note that indices are shared to derive:
% First set = derived cardinal, derived oblique; 
% Second set = derived subset (for dg= radial v tangential, for da= vertical v horizontal)
n_derivedConditions = {{1:2, 3}, {1, 2}};
is_mainsubset = [0, 1];

for ci=1:numel(n_derivedConditions)

    subset = is_mainsubset(ci);

    % NOTE: this is an approximate definition, since they are derived by
    % estimates with/without curvature (closer to the locations 0-315 in 45
    % degree increments are more exact).

    newMatrix = compute_derivativeDirections(meanBOLDpa, projectSettings, subset);

    % Average the values in the first two rows along the first dimension
    % (contrast)
    % for dg: 1:2,3 combines radial and tangential, whereas [1,2] compares radial/tang
    % for da: 1:2,3 combines vertical and horizontal
    averagedMatrix = mean(newMatrix(n_derivedConditions{ci}{1}, :, :, :), 1); 
        
    % Combine the averagedMatrix with the third row
    % for dg: this addition is "other" / polar oblique
    % for da: this addition is cartesian oblique
    proconMatrix = cat(1, averagedMatrix, newMatrix(n_derivedConditions{ci}{2}, :, :, :));

    % Plot polar plots 
    plot1_experimentalCond(proconMatrix, 'derivedCardinalVsDerivedOblique', projectSettings, subset)
    
    % plot each absolute directions
    if (strcmp(projectName, 'da') || strcmp(projectName, 'dots')) && subset==0 %~radialvstang
        newMatrix = compute_derivativeEachAbsDirection(meanBOLDpa, projectSettings);
        plotSepDirs(newMatrix, 'derivedCardinalVsDerivedOblique', 'Abs', projectSettings, subset)
    elseif strcmp(projectName, 'dg') && subset==0 %~radialvstang
        newMatrix = compute_derivativeEachRelDirection(meanBOLDpa, projectSettings);
        plotSepDirs(newMatrix, 'derivedCardinalVsDerivedOblique', 'Rel', projectSettings, subset)
    end

    % this section will average these values to equally weigh each of the 8
    % visual field locations
    
    plot2_experimentalCond(proconMatrix, 'derivedCardinalVsDerivedOblique', projectSettings, subset)

    % same asymmetry, alternate rendering -- see the MAIN CONDITION loop's
    % identical addition above for what this shows.
    plot2_experimentalCond(proconMatrix, 'derivedCardinalVsDerivedOblique', projectSettings, subset, 'model', 'subjectwiseDiff')

end

%% All 4 asymmetries, cortical areas along the x-axis (one figure per
% asymmetry) -- successor to lme1_fit.m's "plot across ROIs" section, not
% part of either loop above since it reads all 4 termIdx from the cached
% fit itself in one call. No projectSettings.fitLabel override -- the
% default (fitLabel=projectName) is what gives dg its own full 13-subject
% cache and da its own 7, per this project's cached fits.
%
% Uses the FULL 8-ROI list (rois, captured above before line 29 restricted
% projectSettings.rois to V1-only for the pairwise-style plots above) --
% this plot is deliberately across all cortical areas, unlike those.
projectSettingsAllROIs = projectSettings;
projectSettingsAllROIs.rois = rois;
plotAsymmetryAcrossROIs(projectSettingsAllROIs)

















