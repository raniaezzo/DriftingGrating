% call to plot trial triggered average

clear all; close all; clc;

% must be in DriftingGrating directory to run

% setup path
addpath(genpath(pwd));
projectName = 'da'; % 'dg', or 'dgl' or 'da'
bidsDir =  '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
%bidsDir = '/Volumes/server/Projects/Project_dg/data_bids/';
githubDir = '~/Documents/GitHub';
fsDir = '/Applications/freesurfer/7.2.0';
%fsDir = '/Applications/freesurfer';
designDir = fullfile('/Volumes/Vision/UsersShare/Rania/Project_dg/experimentalOutput', projectName);
%designDir = '/Volumes/server/Projects/Project_dg/experimentalOutput/dg/';
addpath(genpath(fullfile(githubDir, 'wpToolbox')));
addpath(genpath(fullfile(githubDir, 'atlasmgz')));
setup_user('rania', bidsDir)

projectSettings = loadConfig(githubDir);

hRF_setting = 'glmsingle'; %
subjects = {'sub-0426'}; %{'sub-0426', 'sub-wlsubj123', 'sub-wlsubj124', 'sub-0255', 'sub-0037'};
%{'sub-0201', 'sub-0250', 'sub-0427', 'sub-0426', 'sub-0397', 'sub-wlsubj123', 'sub-wlsubj124', 'sub-0255', 'sub-0037'};
%subj = 'sub-0442'; %'sub-0426';

%%
for ss=1:numel(subjects)
    subj = subjects{ss};

    % get ses name
    subjectDir = fullfile(bidsDir,'derivatives',strcat(projectName, 'GLM'), strcat('hRF_',hRF_setting), subj);
    contents = dir(subjectDir);
    sesNames = {};
    for i = 1:length(contents)
        % Check if the item is a directory and starts with 'ses-'
        if contents(i).isdir && startsWith(contents(i).name, 'ses-')
            % Add the subfolder name to the cell array
            sesNames{end+1} = contents(i).name;
        end
    end
    if length(sesNames)>1
        error('Multiple session folder.')
    else
        ses = sesNames{1};
    end

    retFolder = 'prfvista_mov';
    
    derivativesFolder = fullfile(bidsDir, 'derivatives', sprintf('%sGLM',projectName), ...
        sprintf('hRF_%s', hRF_setting), subj, ses);
    
    % load data and set up
    
    % get ROI surface
    hSize = get_surfsize(subj); 
    
    % load in design matrix
    load(fullfile(derivativesFolder, 'rawInfo.mat'));
    
    % load raw signal data (in wpToolbox) - convert to .mgh if not already
    run = 1:length(matrices_onset);
    datafiles = load_data(bidsDir,projectName,'fsnative','.mgh',subj,ses,run) ;
    
    %conditions = {'cardinaloblique', 'radialtangential'};
    conditions = {'mainCardinalVsMainOblique', 'derivedCardinalVsDerivedOblique'};
    
    projectSettings.projectName = projectName;
    projectSettings.subject= subj;
    projectSettings.bidsDir = bidsDir;
    projectSettings.retFolder = retFolder;
    
    % filter params
    projectSettings.polarAngleBinWidth = 45; %45; % degrees
    projectSettings.minECC = 1; %0.2; 
    projectSettings.maxECC = 15;
    projectSettings.minVAREXP = .1; %0.25;
    projectSettings.stimdur_s = 3;
    projectSettings.tr_s = 1;
    projectSettings.eventTRs_after = 20; % TRs to plot (ttave)
    projectSettings.eventTRs_prior = 5;
    projectSettings.padding = 10; % TRs
    
    %%
    
    % can be 'motion_minus_orientation' or 'motion_minus_baseline' (contains 'orientation_minus_baseline')
    comparisonNames = {'motion_minus_orientation','motion_minus_baseline'}; %%, };
    
    for cn=1:numel(comparisonNames)
    
        projectSettings.comparisonName = comparisonNames{cn};
    
        ttaveSavePath = fullfile(derivativesFolder, 'ttaveData', comparisonNames{cn});
        projectSettings.ttaveSavePath = ttaveSavePath;
    
        if ~isfolder(ttaveSavePath)
            mkdir(ttaveSavePath)
        end
    
        for ri=1:length(projectSettings.rois)
        
            % get ROI indices
            projectSettings.roiName = projectSettings.rois{ri};
            primaryROIvertices = getROIidxs(subj, projectSettings.roiName, hSize); % the name of Middle Temporal Cortex ROI
            surfaceROI = nan(sum(hSize),1);
            surfaceROI(primaryROIvertices) = 1;
        
            % get VF indices binned into 8 categories (after filtering for: ecc, r^2)
            projectSettings.filteredPrfBins = retriveRetData(projectSettings);
        
            % do MAIN CONDITION first (for dg, would only be card-v-oblique; for da, polCard-v-obl and radial-v-tang)
        
            % this will plot the derived conditions for:
            if strcmp(projectName, 'da') % experiment 2: da (polar cardinal vs oblique; and radial vs tangential)
                isradial = [0, 1]; 
            elseif strcmp(projectName, 'dg') % experiment 1: dg (cartesian cardinal vs oblique)
                isradial = 0; 
            end
            
            for ci=1:numel(isradial)
        
                radialvstang = isradial(ci);
        
                % compute and plot ttave main dir
                ttave_compute(matrices_onset,datafiles, 'mainCardinalVsMainOblique', surfaceROI, projectSettings, radialvstang);
        
            end
        
            % now do DERIVED CONDITION (for dg, would be polCard-v-oblique and
            % radial-v-tang; for da would only be cartCardinal-v-oblique
        
            % this will plot the derived conditions for:
            if strcmp(projectName, 'da') % experiment 2: da (cartesian cardinal vs oblique)
                %n_derivedConditions = {{1:2, 3}};
                isradial = 0; 
            elseif strcmp(projectName, 'dg') % experiment 1: dg (polar cardinal vs oblique; and radial vs tangential)
                %n_derivedConditions = {{1:2, 3}, {1, 2}};
                isradial = [0, 1]; 
            end
        
            for ci=1:numel(isradial)
        
                radialvstang = isradial(ci);
        
                % compute and plot ttave main dir
                ttave_compute(matrices_onset, datafiles, 'derivedCardinalVsDerivedOblique', surfaceROI, projectSettings, radialvstang);
        
            end
        
        end
        
        close all;
    end
end
