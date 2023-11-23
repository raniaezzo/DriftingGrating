clear all; close all; clc;
%%

stimdur_s = 3;
tr_s = 1;
hrf_opt = 'assume';
resampling = 0; % means fit fully (don't bootstrap or cross-validate)

% setup path
addpath(genpath(pwd));
projectName = 'dg';
bidsDir =  '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
githubDir = '~/Documents/GitHub';
fsDir = '/Applications/freesurfer/7.2.0';
designDir = '/Volumes/Vision/UsersShare/Rania/Project_dg/experimentalOutput/dg/';
addpath(genpath(fullfile(githubDir, 'wpToolbox')));
addpath(genpath(fullfile(githubDir, 'atlasmgz')));
setup_user(projectName,bidsDir,githubDir,fsDir);

% define params
subj = 'sub-wlsubj123';
ses = 'ses-01'; %'ses-nyu3t02'; %'ses-01';
space = 'fsnative';
rois = {'mt','mstL','mstR','fst'};

% Script Part A) Description
% This runs a basic GLM where each unique event (13) is a unique regressor
% (1) M_0, 
% (2) M_90, 
% (3) M_180, 
% (4) M_270, 
% (5) M_45, 
% (6) M_135, 
% (7) M_225, 
% (8) M_315, 
% (9) S_0, 
% (10) S_90, 
% (11) S_45, 
% (12) S_135, 
% (13) B
% Produces 13 beta maps (one for each regressor).

% check if design files exist, if not convert the matrices from experiment
% into the appropriate format
matrices_init = format_desmats(bidsDir, designDir, subj, ses, tr_s);

% add nuissance regressors here to matrices.

matrices = add_nuissancemats(bidsDir,subj, ses, matrices_init);

%%
% load data (in wpToolbox) - convert to .mgh if not already
run = 1:length(matrices);

datafiles = load_data(bidsDir,projectName,space,'.mgh',subj,ses,run) ;

%%

% run glm (from GLMdenoise repo)
results.(char('allevents')) = GLMestimatemodel(matrices, datafiles,...
    stimdur_s,tr_s,hrf_opt,[],resampling);

%%
% Script Part B) Description
% Combines beta maps per condition and produces contrast maps:
% average betas for cardinal motion - average betas for cardinal static
% average betas for oblique motion - average betas for oblique static

betamaps = results.allevents.modelmd{1,2};

derivativesFolder = fullfile(bidsDir, 'derivatives','dgGLM');

if ~isfolder(derivativesFolder)
    mkdir(derivativesFolder)
end

% average beta maps for cardinal moving - cardinal static
contrasts.cardmVcards = mean(betamaps(:,1:4),2) - mean(betamaps(:,9:10),2);
writeMGZfile(bidsDir, subj, ses, contrasts.cardmVcards, derivativesFolder, 'cardm_v_cards')

% average beta maps for oblique moving - oblique static
contrasts.oblmVobls = mean(betamaps(:,5:8),2) - mean(betamaps(:,11:12),2);
writeMGZfile(bidsDir, subj, ses, contrasts.oblmVobls, derivativesFolder, 'oblm_v_obls')

% average beta maps for all moving - all static
contrasts.allmValls = mean(betamaps(:,1:8),2) - mean(betamaps(:,9:12),2);
writeMGZfile(bidsDir, subj, ses, contrasts.allmValls, derivativesFolder, 'allm_v_alls')

% average beta maps for all static - blank
contrasts.allsVblank = mean(betamaps(:,9:12),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.allsVblank, derivativesFolder, 'alls_v_allb')

% average beta maps for all moving - blank
contrasts.allmVblank = mean(betamaps(:,1:8),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.allmVblank, derivativesFolder, 'allm_v_allb')

% average beta maps for all cardinal moving - blank
contrasts.cardmVblank = mean(betamaps(:,1:4),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.allsVblank, derivativesFolder, 'cardm_v_allb')

% average beta maps for all oblique moving - blank
contrasts.oblmVblank = mean(betamaps(:,5:8),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.allmVblank, derivativesFolder, 'oblm_v_allb')

results.contrasts = contrasts;

save(fullfile(derivativesFolder, subj, ses, 'results.mat'), 'results')

%% Visualize

view_fv(subj,bidsDir,'allm_v_alls')

%%








% mask with MT+ ROI ?
% do this later

%% visualize R2 for each roi
close all;

returndir = pwd;

% Generate wang and glasser atlases
glasser_dir = fullfile(bidsDir, 'derivatives', 'freesurfer', subj, 'label', 'Glasser2016');
if ~isfolder(glasser_dir)
    cd(fullfile(githubDir, 'atlasmgz'))
    system(['sh createAtlasLabels.sh ' extractAfter(subj,'-') ' ' bidsDir])
end

cd(returndir);

%'cardmVcards','oblmVobls','allmValls', 

visualize(results.contrasts, subj, {'allmValls'})
% visualize(results,subj,rois);


%view_fv(subj,bidsDir,'mt+2')