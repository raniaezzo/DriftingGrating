% plot BOLD asymmetries

clc; clear all; close all;

% set up
addpath(genpath(pwd));
projectName = 'dg';
bidsDir =  '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
githubDir = '~/Documents/GitHub';
hRF_setting = 'glmsingle';
fullfile(githubDir, 'DriftingGrating', 'AnalysisCode')
glmResultsfolder = fullfile(bidsDir, 'derivatives', strcat(projectName, 'GLM'), strcat('hRF_', hRF_setting));

% can be 'motion_minus_orientation' ; 'motion_minus_baseline' ; 'orientation_minus_baseline'
comparisonName = 'motion_minus_orientation';

projectSettings = loadConfig(githubDir);

rois = projectSettings.rois;
axes_limits = projectSettings.axes_limits;
pairaxes_limits = projectSettings.pairaxes_limits;
pairaxes_PAew_limits = projectSettings.pairaxes_PAew_limits;
colors_data = projectSettings.colors_data;
contrasts_dict = projectSettings.contrasts_dict;

%rois = rois(1:7); % remove once I process v3a / v3b

%% Load and remove subject with extreme motion

load(fullfile(glmResultsfolder, 'medianBOLDpa')) % contrasts x polarAngles x ROIs x subjects
load(fullfile(glmResultsfolder, 'medianBOLD')) % contrasts x ROIs x subjects

figureDir = [strrep(bidsDir, 'data_bids', 'figures'), projectName];

if ~isfolder(figureDir)
    mkdir(figureDir)
end

if strcmp(projectName, 'dg')
    medianBOLDpa = medianBOLDpa(:,:,:,[1,2,3,4,5,6,7,8,9,10,11,13]); % leave out 12 and include 13 instead
    medianBOLD = medianBOLD(:,:,[1,2,3,4,5,6,7,8,9,10,11,13]);
    % to leave out one subject
    subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-0397', ...
        'sub-0442', 'sub-wlsubj121', ...
        'sub-wlsubj123', 'sub-wlsubj124', 'sub-wlsubj127', 'sub-0395', 'sub-0426', ...
        'sub-0250'};
elseif strcmp(projectName, 'da')
    % load subjects
    subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
        'sub-0395', 'sub-0426'};
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


%% MAIN CONDITION: Plot pairwise plots (JUST FOR SANITY CHECK)
% e.g. condition indices 4 v 5: clearly shows MT as motion responsive
% these do not depend on polar angle / retinotopy (just contrasts)

condIdx1 = 9;   
condIdx2 = 8;

plot0_experimentalCond(condIdx1, condIdx2, medianBOLD, projectSettings)

%% MAIN CONDITION: Plot polor plots

% this will plot the derived conditions for:
if strcmp(projectName, 'da') % experiment 2: da (polar cardinal vs oblique; and radial vs tangential)
    isradial = [0, 1]; 
elseif strcmp(projectName, 'dg') % experiment 1: dg (cartesian cardinal vs oblique)
    isradial = 0; 
end

for ci=1:numel(isradial)

    radialvstang = isradial(ci);
    % this will plot the main conditions for experiment 1: dg (cardinal vs oblique) 
    % and for experiment 2: da (polar cardinal vs oblique; and radial vs
    % tangential)
    
    % this depends on retinotopy to isolate each direction with respect to its
    % location, but still only plots the main conditions
    plot1_experimentalCond(medianBOLDpa, 'mainCardinalVsMainOblique', projectSettings, radialvstang)

    % MAIN CONDITION: Plot pairwise plots (EQUALLY WEIGHING POLAR ANGLE)
    % this section will average these values to equally weigh each of the 8
    % visual field locations
    
    plot2_experimentalCond(medianBOLDpa, 'mainCardinalVsMainOblique', projectSettings, radialvstang)

end

%% DERIVED CONDITION
% (THESE DEPEND ON THE POLAR ANGLE WEDGE SIZE)

% this will plot the derived conditions for:
if strcmp(projectName, 'da') % experiment 2: da (cartesian cardinal vs oblique)
    n_derivedConditions = {{1:2, 3}};
    isradial = 0; 
elseif strcmp(projectName, 'dg') % experiment 1: dg (polar cardinal vs oblique; and radial vs tangential)
    n_derivedConditions = {{1:2, 3}, {1, 2}};
    isradial = [0, 1]; 
end

for ci=1:numel(n_derivedConditions)

    radialvstang = isradial(ci);

    % NOTE: this is an approximate definition, since they are derived by
    % estimates with/without curvature (closer to the locations 0-315 in 45
    % degree increments are more exact).

    newMatrix = compute_derivativeDirections(medianBOLDpa, projectSettings, radialvstang);

    % Average the values in the first two rows along the first dimension
    % (contrast)
    % for dg: 1:2 combines radial and tangential
    % for da: 1:2 combines vertical and horizontal
    averagedMatrix = mean(newMatrix(n_derivedConditions{ci}{1}, :, :, :), 1); 
        
    % Combine the averaged values with the third row
    % for dg: this addition is "other" / polar oblique
    % for da: this addition is cartesian oblique
    proconMatrix = cat(1, averagedMatrix, newMatrix(n_derivedConditions{ci}{2}, :, :, :));

    % Plot polar plots 
    plot1_experimentalCond(proconMatrix, 'derivedCardinalVsDerivedOblique', projectSettings, radialvstang)
    
    % this section will average these values to equally weigh each of the 8
    % visual field locations
    
    plot2_experimentalCond(proconMatrix, 'derivedCardinalVsDerivedOblique', projectSettings, radialvstang)

end


















