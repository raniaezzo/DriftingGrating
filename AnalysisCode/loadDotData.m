% make a file for dots -- with same format for easy analysis:
    %/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/derivatives/
    % dgGLM/hRF_glmsingle/sub-0037/ses-01/results.mat


clc; clear all; close all;

% set up
addpath(genpath(pwd));
projectName = 'dots';
bidsDir =  '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
githubDir = '~/Documents/GitHub';
fsDir = '/Applications/freesurfer'; %/7.2.0';
addpath(genpath(fullfile(githubDir, 'atlasmgz')));
setup_user('rania',bidsDir);

hRF_setting = 'canonical'; % can be: 'canonical', 'glmdenoise', 'glmsingle';

subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', 'sub-0426'};

% these two do not have separate mgz for in, out, cw, ccw -- need
% reprocessing 
% 'sub-0395',  %'sub-0250'};

dotConds = {'out.mgz','in'};

% for each subject, save a results.mat inside dots/canonical/ folder
for si=1:numel(subjects)

    subjectname = subjects{si}

    fspth = fullfile(bidsDir, 'derivatives', 'freesurfer', subjectname);
    glmResultsfolder = fullfile(bidsDir, 'derivatives', strcat(projectName, 'GLM'), strcat('hRF_', hRF_setting), subjectname);
    % load:  in.mgz, out.mgz, cw.mgz and ccw.mgz
    % assign to 270, 90, 0, 180
    
    % find out the length of volume:
    %mgz = MRIread(fullfile(fspth, 'mri', 'orig.mgz'));
    %lcurv = read_curv(fullfile(fspth, 'surf', 'lh.curv'));
    %rcurv = read_curv(fullfile(fspth, 'surf', 'rh.curv'));

    surfsize = get_surfsize(subjectname);
    surfconcat_length = surfsize(1)+surfsize(2);

    results.contrasts.cardmVcards = nan(surfconcat_length,1);
    results.contrasts.oblmVobls = nan(surfconcat_length,1);
    results.contrasts.allmValls = nan(surfconcat_length,1);
    results.contrasts.allsVblank = nan(surfconcat_length,1);
    results.contrasts.allmVblank = nan(surfconcat_length,1);
    results.contrasts.cardsVblank = nan(surfconcat_length,1);
    results.contrasts.oblsVblank = nan(surfconcat_length,1);
    results.contrasts.cardmVblank = nan(surfconcat_length,1);
    results.contrasts.oblmVblank = nan(surfconcat_length,1);
    
    results.contrasts.m0Vs90 = nan(surfconcat_length,1);

    lh_cw = MRIread(fullfile(glmResultsfolder, 'lh.cw.mgz'));
    rh_cw = MRIread(fullfile(glmResultsfolder, 'rh.cw.mgz'));
    results.contrasts.m0Vb = [lh_cw.vol; rh_cw.vol];     % clockwise dots

    results.contrasts.m180Vs90 = nan(surfconcat_length,1);

    lh_ccw = MRIread(fullfile(glmResultsfolder, 'lh.ccw.mgz'));
    rh_ccw = MRIread(fullfile(glmResultsfolder, 'rh.ccw.mgz'));
    results.contrasts.m180Vb = [lh_ccw.vol; rh_ccw.vol];     % counterclockwise dots

    results.contrasts.m90Vs0 = nan(surfconcat_length,1);

    lh_out = MRIread(fullfile(glmResultsfolder, 'lh.out.mgz'));
    rh_out = MRIread(fullfile(glmResultsfolder, 'rh.out.mgz'));
    results.contrasts.m90Vb = [lh_out.vol; rh_out.vol];     % outwards dots

    results.contrasts.m270Vs0 = nan(surfconcat_length,1);

    lh_in = MRIread(fullfile(glmResultsfolder, 'lh.in.mgz'));
    rh_in = MRIread(fullfile(glmResultsfolder, 'rh.in.mgz'));
    results.contrasts.m270Vb = [lh_in.vol; rh_in.vol];     % inwards dots

    results.contrasts.m45Vs135 = nan(surfconcat_length,1);
    results.contrasts.m45Vb = nan(surfconcat_length,1);

    results.contrasts.m225Vs135 = nan(surfconcat_length,1);
    results.contrasts.m225Vb = nan(surfconcat_length,1);

    results.contrasts.m135Vs45 = nan(surfconcat_length,1);
    results.contrasts.m135Vb = nan(surfconcat_length,1);

    results.contrasts.m315Vs45 = nan(surfconcat_length,1);
    results.contrasts.m315Vb = nan(surfconcat_length,1);
    
    results.contrasts.cardMsep = nan(surfconcat_length,1);
    results.contrasts.oblMsep = nan(surfconcat_length,1);
    
    results.contrasts.s0Vb = nan(surfconcat_length,1);
    results.contrasts.s90Vb = nan(surfconcat_length,1);
    results.contrasts.s45Vb = nan(surfconcat_length,1);
    results.contrasts.s135Vb = nan(surfconcat_length,1);

    save(fullfile(glmResultsfolder, 'results.mat'), 'results');

end



