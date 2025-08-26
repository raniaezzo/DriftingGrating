clear all; close all; clc;

%% SET UP

% set fundamental vars to run a single subject / protocol
projectName = 'da'; % 'dg', or 'dgl' or 'da'
subj = 'sub-0250'; %'sub-0426';
ses = 'ses-03'; %'ses-01'; %'ses-nyu3t02'; %'ses-01'; <-- can I make this not manual?
space = 'fsnative';

% checks/sets working directory and adds dependencies
setup_user('rania', bidsDir)
jsonparamfile = 'setup.json';
jsonParams = jsondecode(fileread(jsonparamfile));

% load predefined paths
bidsDir = jsonParams.bidsdir.Path; expOutputDir = jsonParams.expoutputdir.Path; 
stimdur_s = jsonParams.stimdur_s.Val; tr_s = jsonParams.tr_s.Val;

% project specific paths
designDir = fullfile(expOutputDir, projectName);

% define params for hRF
hRF_setting = 'glmsingle'; % can be: 'canonical', 'glmdenoise', 'glmsingle';

if strcmp(hRF_setting, 'glmdenoise_canonical')
    hrf_opt = 'assume';
elseif strcmp(hRF_setting, 'glmdenoise_optimizeGlobal')
    hrf_opt = 'optimize';
elseif strcmp(hRF_setting, 'glmsingle') % GLMestimatesingletrial
    hrf_opt = [];
end

resampling = 0; % means fit fully (don't bootstrap or cross-validate)

derivativesFolder = fullfile(bidsDir, 'derivatives', sprintf('%sGLM',projectName), ...
    sprintf('hRF_%s', hRF_setting), subj, ses);

if ~isfolder(derivativesFolder)
    mkdir(derivativesFolder)
end

%% CONSTRUCT DESIGN MATRICES; LOAD DATA

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
[matrices_init, ~, ~, ~] = format_desmats(bidsDir, designDir, subj, tr_s);

% get design matrix with only trial onsets (rather than duration)
matrices_onset = filterOnset(matrices_init);

% add nuissance regressors here to matrices.
matrices = add_nuissancemats(bidsDir,subj, ses, matrices_onset); % for manual GLM

% for GLMsingle/GLMdenoise if including nuissance regressors (not needed)
C = cell(1, size(matrices_onset,2));
matrices_nuis = add_nuissancemats(bidsDir,subj, ses, C); 

% load data (in wpToolbox) - convert to .mgh if not already
run = 1:length(matrices_onset);

datafiles = load_data(bidsDir,projectName,space,'.mgh',subj,ses,run) ;

%% FIT GLM

% THREE ways to run GLM (nuissance regressors not used)-- all settings can be changed with opt struct:
%1- glmdenoise_canonical: runs GLM with canonical hRF, then provides denoising procedure
%2- glmdenoise_optimizeGlobal: optimizes hRF (finds 1 solution), then
%provides denoising procedure
%3- glmsingle: glm selection from a library of 20 options, GLMdenoise, then
% ridge regression

if strcmp(hRF_setting, 'glmdenoise_canonical') || strcmp(hRF_setting, 'glmdenoise_optimizeGlobal')

    opt = [];
    %opt = struct('extraregressors',{matrices_nuis},'hrfthresh',-Inf,'suppressoutput',1);

    results.(char('allevents')) = GLMestimatemodel(matrices_onset, datafiles,... % replaced matrices with matrices_onset
    stimdur_s,tr_s,hrf_opt,[],resampling,opt); % default says PSC is the result

    % % run glm (from GLMdenoise repo)
    % results.(char('allevents')) = GLMestimatemodel(matrices, datafiles,...
    %     stimdur_s,tr_s,hrf_opt,[],resampling); % default says PSC is the result

elseif strcmp(hRF_setting, 'glmsingle')

    opt = struct('wantmemoryoutputs',[1 1 1 1]);
    % opt.extraregressors = matrices_nuis;

    % default runs glm library selection, GLMdenoise, and ridge reg
    % for now, not including noise regressors - getting the function take
    % care of it
    [modelOut, designSINGLE] = GLMestimatesingletrial(matrices_onset,datafiles,...
        stimdur_s,tr_s,derivativesFolder,opt);
    results.(char('allevents')) = modelOut{1,4};
end

% save results  % default says PSC is the result
sprintf('SAVING VARIABLES..')
save(fullfile(derivativesFolder,'rawInfo.mat'), 'matrices_onset', 'matrices', ...
    'matrices_nuis', 'datafiles', 'tr_s', 'stimdur_s', '-v7.3');
save(fullfile(derivativesFolder,'modelOutput.mat'), 'modelOut', 'designSINGLE', '-v7.3');


%% NORMALIZE BETAS, SAVE CONTRASTS

% Script Part B) Description
% Combines beta maps per condition and produces contrast maps:
% average betas for cardinal motion - average betas for cardinal static
% average betas for oblique motion - average betas for oblique static

normalize = 1;

% average betas across trials if I ran glmSingle
if ~strcmp(hRF_setting, 'glmsingle')
    betamaps = results.allevents.modelmd{1,2};
else
    [nvert, ~, ~, ntrials] = size(results.allevents.modelmd);

    betamaps = nan(nvert, 13); % conditions
    % plot a condition
    for ci=1:13
        condSelect = designSINGLE.stimorder==ci;
        betas = results.allevents.modelmd(:,:,:,condSelect);
        newbetas = squeeze(betas);
        
        % trial mean
        betamaps(:,ci) = mean(newbetas,2);
    end
end

if normalize
    % new: convert to z-score (across betas)
    betamaps_condOnly = betamaps(:,1:13);
    betamaps = zscore(betamaps_condOnly, 0, 2);
end

mgzFolder = fullfile(bidsDir, 'derivatives', sprintf('%sGLM',projectName), ...
    sprintf('hRF_%s', hRF_setting));

% average beta maps for cardinal moving - cardinal static
contrasts.cardmVcards = mean(betamaps(:,1:4),2) - mean(betamaps(:,9:10),2);
writeMGZfile(bidsDir, subj, ses, contrasts.cardmVcards, mgzFolder, 'cardm_v_cards')

% average beta maps for oblique moving - oblique static
contrasts.oblmVobls = mean(betamaps(:,5:8),2) - mean(betamaps(:,11:12),2);
writeMGZfile(bidsDir, subj, ses, contrasts.oblmVobls, mgzFolder, 'oblm_v_obls')

% average beta maps for all moving - all static
contrasts.allmValls = mean(betamaps(:,1:8),2) - mean(betamaps(:,9:12),2);
writeMGZfile(bidsDir, subj, ses, contrasts.allmValls, mgzFolder, 'allm_v_alls')

% average beta maps for all static - blank
contrasts.allsVblank = mean(betamaps(:,9:12),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.allsVblank, mgzFolder, 'alls_v_allb')

% average beta maps for all moving - blank
contrasts.allmVblank = mean(betamaps(:,1:8),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.allmVblank, mgzFolder, 'allm_v_allb')

% average beta maps for cardinal static - baseline
contrasts.cardsVblank = mean(betamaps(:,9:10),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.cardsVblank, mgzFolder, 'cards_v_allb')

% average beta maps for oblique static - baseline
contrasts.oblsVblank = mean(betamaps(:,11:12),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.oblsVblank, mgzFolder, 'obls_v_allb')

% average beta maps for all cardinal moving - blank
contrasts.cardmVblank = mean(betamaps(:,1:4),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.allsVblank, mgzFolder, 'cardm_v_allb')

% average beta maps for all oblique moving - blank
contrasts.oblmVblank = mean(betamaps(:,5:8),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.allmVblank, mgzFolder, 'oblm_v_allb')

% below are contrast for each direction separately
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


% m_0 (horizontal motion) vs s_90 (vertical orientation)
contrasts.m0Vs90 = mean(betamaps(:,1),2) - mean(betamaps(:,10),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m0Vs90, mgzFolder, 'm0_v_s90')

% m_0 (horizontal motion) vs blank
contrasts.m0Vb = mean(betamaps(:,1),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m0Vb, mgzFolder, 'm0_v_b')

% m_180 (horizontal motion) vs s_90 (vertical orientation)
contrasts.m180Vs90 = mean(betamaps(:,3),2) - mean(betamaps(:,10),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m180Vs90, mgzFolder, 'm180_v_s90')

% m_180 (horizontal motion) vs blank
contrasts.m180Vb = mean(betamaps(:,3),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m180Vb, mgzFolder, 'm180_v_b')

% m_90 (vertical motion) vs s_0 (horizontal orientation)
contrasts.m90Vs0 = mean(betamaps(:,2),2) - mean(betamaps(:,9),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m90Vs0, mgzFolder, 'm90_v_s0')

% m_90 (vertical motion) vs blank
contrasts.m90Vb = mean(betamaps(:,2),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m90Vb, mgzFolder, 'm90_v_b')

% m_270 (vertical motion) vs s_0 (horizontal orientation)
contrasts.m270Vs0 = mean(betamaps(:,4),2) - mean(betamaps(:,9),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m270Vs0, mgzFolder, 'm270_v_s0')

% m_270 (vertical motion) vs blank
contrasts.m270Vb = mean(betamaps(:,4),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m270Vb, mgzFolder, 'm270_v_b')

% m_45 (upper right motion) vs s_135 (cc tilt orientation)
contrasts.m45Vs135 = mean(betamaps(:,5),2) - mean(betamaps(:,12),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m45Vs135, mgzFolder, 'm45_v_s135')

% m_45 (upper right motion) vs blank
contrasts.m45Vb = mean(betamaps(:,5),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m45Vb, mgzFolder, 'm45_v_b')

% m_225 (lower left motion) vs s_135 (cc tilt orientation)
contrasts.m225Vs135 = mean(betamaps(:,7),2) - mean(betamaps(:,12),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m225Vs135, mgzFolder, 'm225_v_s135')

% m_225 (lower left motion) vs blank
contrasts.m225Vb = mean(betamaps(:,7),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m225Vb, mgzFolder, 'm225_v_b')

% m_135 (upper left motion) vs s_45 (c tilt orientation)
contrasts.m135Vs45 = mean(betamaps(:,6),2) - mean(betamaps(:,11),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m135Vs45, mgzFolder, 'm135_v_s45')

% m_135 (upper left motion) vs blank
contrasts.m135Vb = mean(betamaps(:,6),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m135Vb, mgzFolder, 'm135_v_b')

% m_315 (lower right motion) vs s_45 (c tilt orientation)
contrasts.m315Vs45 = mean(betamaps(:,8),2) - mean(betamaps(:,11),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m315Vs45, mgzFolder, 'm315_v_s45')

% m_315 (lower right motion) vs blank
contrasts.m315Vb = mean(betamaps(:,8),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.m315Vb, mgzFolder, 'm315_v_b')

% averaging after subtracting like condition

cardSep = [contrasts.m0Vs90, contrasts.m180Vs90, contrasts.m90Vs0, contrasts.m270Vs0];
contrasts.cardMsep = mean(cardSep,2);
writeMGZfile(bidsDir, subj, ses, contrasts.cardMsep, mgzFolder, 'cardm_pairwise')

oblSep = [contrasts.m45Vs135, contrasts.m225Vs135, contrasts.m135Vs45, contrasts.m315Vs45];
contrasts.oblMsep = mean(oblSep,2);
writeMGZfile(bidsDir, subj, ses, contrasts.oblMsep, mgzFolder, 'oblm_pairwise')

% orientations isolated

% s_0 vs blank
contrasts.s0Vb = mean(betamaps(:,9),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.s0Vb, mgzFolder, 's0_v_b')

% s_90 vs blank
contrasts.s90Vb = mean(betamaps(:,10),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.s90Vb, mgzFolder, 's90_v_b')

% s_45 vs blank
contrasts.s45Vb = mean(betamaps(:,11),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.s45Vb, mgzFolder, 's45_v_b')

% s_135 vs blank
contrasts.s135Vb = mean(betamaps(:,12),2) - mean(betamaps(:,13),2);
writeMGZfile(bidsDir, subj, ses, contrasts.s135Vb, mgzFolder, 's135_v_b')

results.contrasts = contrasts;

save(fullfile(derivativesFolder, 'results.mat'), 'results')
