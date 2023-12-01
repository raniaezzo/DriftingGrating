%clear all; close all; clc;

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
subj = 'sub-wlsubj127';
ses = 'ses-nyu3t02'; %'ses-nyu3t02'; %'ses-01';
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

% average beta maps for cardinal static - baseline
contrasts.cardsVblank = mean(betamaps(:,9:10),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.cardsVblank, derivativesFolder, 'cards_v_allb')

% average beta maps for oblique static - baseline
contrasts.oblsVblank = mean(betamaps(:,11:12),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.oblsVblank, derivativesFolder, 'obls_v_allb')

% average beta maps for all cardinal moving - blank
contrasts.cardmVblank = mean(betamaps(:,1:4),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.allsVblank, derivativesFolder, 'cardm_v_allb')

% average beta maps for all oblique moving - blank
contrasts.oblmVblank = mean(betamaps(:,5:8),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.allmVblank, derivativesFolder, 'oblm_v_allb')

% below are contrast for each direction separately

% m_0 (horizontal motion) vs s_90 (vertical orientation)
contrasts.m0Vs90 = mean(betamaps(:,1),2) - mean(betamaps(:,10),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m0Vs90, derivativesFolder, 'm0_v_s90')

% m_0 (horizontal motion) vs blank
contrasts.m0Vb = mean(betamaps(:,1),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m0Vb, derivativesFolder, 'm0_v_b')

% m_180 (horizontal motion) vs s_90 (vertical orientation)
contrasts.m180Vs90 = mean(betamaps(:,3),2) - mean(betamaps(:,10),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m180Vs90, derivativesFolder, 'm180_v_s90')

% m_180 (horizontal motion) vs blank
contrasts.m180Vb = mean(betamaps(:,3),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m180Vb, derivativesFolder, 'm180_v_b')

% m_90 (vertical motion) vs s_0 (horizontal orientation)
contrasts.m90Vs0 = mean(betamaps(:,2),2) - mean(betamaps(:,9),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m90Vs0, derivativesFolder, 'm90_v_s0')

% m_90 (vertical motion) vs blank
contrasts.m90Vb = mean(betamaps(:,2),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m90Vb, derivativesFolder, 'm90_v_b')

% m_270 (vertical motion) vs s_0 (horizontal orientation)
contrasts.m270Vs0 = mean(betamaps(:,4),2) - mean(betamaps(:,9),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m270Vs0, derivativesFolder, 'm270_v_s0')

% m_270 (vertical motion) vs blank
contrasts.m270Vb = mean(betamaps(:,4),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m270Vb, derivativesFolder, 'm270_v_b')

% m_45 (upper right motion) vs s_135 (cc tilt orientation)
contrasts.m45Vs135 = mean(betamaps(:,5),2) - mean(betamaps(:,12),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m45Vs135, derivativesFolder, 'm45_v_s135')

% m_45 (upper right motion) vs blank
contrasts.m45Vb = mean(betamaps(:,5),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m45Vb, derivativesFolder, 'm45_v_b')

% m_225 (lower left motion) vs s_135 (cc tilt orientation)
contrasts.m225Vs135 = mean(betamaps(:,7),2) - mean(betamaps(:,12),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m225Vs135, derivativesFolder, 'm225_v_s135')

% m_225 (lower left motion) vs blank
contrasts.m225Vb = mean(betamaps(:,7),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m225Vb, derivativesFolder, 'm225_v_b')

% m_135 (upper left motion) vs s_45 (c tilt orientation)
contrasts.m135Vs45 = mean(betamaps(:,6),2) - mean(betamaps(:,11),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m135Vs45, derivativesFolder, 'm135_v_s45')

% m_135 (upper left motion) vs blank
contrasts.m135Vb = mean(betamaps(:,6),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m135Vb, derivativesFolder, 'm135_v_b')

% m_315 (lower right motion) vs s_45 (c tilt orientation)
contrasts.m315Vs45 = mean(betamaps(:,8),2) - mean(betamaps(:,11),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m315Vs45, derivativesFolder, 'm315_v_s45')

% m_315 (lower right motion) vs blank
contrasts.m315Vb = mean(betamaps(:,8),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m315Vb, derivativesFolder, 'm315_v_b')

% averaging after subtracting like condition

cardSep = [contrasts.m0Vs90, contrasts.m180Vs90, contrasts.m90Vs0, contrasts.m270Vs0];
contrasts.cardMsep = mean(cardSep,2);
writeMGZfile(bidsDir, subj, ses, contrasts.cardMsep, derivativesFolder, 'cardm_pairwise')

oblSep = [contrasts.m45Vs135, contrasts.m225Vs135, contrasts.m135Vs45, contrasts.m315Vs45];
contrasts.oblMsep = mean(oblSep,2);
writeMGZfile(bidsDir, subj, ses, contrasts.oblMsep, derivativesFolder, 'oblm_pairwise')

% orientations isolated

% s_0 vs blank
contrasts.s0Vb = mean(betamaps(:,9),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.s0Vb, derivativesFolder, 's0_v_b')

% s_90 vs blank
contrasts.s90Vb = mean(betamaps(:,10),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.s90Vb, derivativesFolder, 's90_v_b')

% s_45 vs blank
contrasts.s45Vb = mean(betamaps(:,11),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.s45Vb, derivativesFolder, 's45_v_b')

% s_135 vs blank
contrasts.s135Vb = mean(betamaps(:,12),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.s135Vb, derivativesFolder, 's135_v_b')

results.contrasts = contrasts;

save(fullfile(derivativesFolder, subj, ses, 'results.mat'), 'results')

%% Visualize

bidsDir = '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/'; %'/Volumes/Vision/MRI/recon-bank';
subj = 'sub-wlsubj123';

view_fv(subj,bidsDir,'cardm_v_cards', 'oblm_v_obls')




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