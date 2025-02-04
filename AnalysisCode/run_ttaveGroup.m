% group ttave

clc; clear all; close all;

% set up
addpath(genpath(pwd));
projectName = 'dg';
bidsDir =  '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
%bidsDir = '/Volumes/server/Projects/Project_dg/data_bids/';
githubDir = '~/Documents/GitHub';
fsDir = '/Applications/freesurfer'; %/7.2.0';
%addpath(genpath(fullfile(githubDir, 'retinotopy-nyuad')));
addpath(genpath(fullfile(githubDir, 'wpToolbox')));
addpath(genpath(fullfile(githubDir, 'atlasmgz')));
%setup_user(projectName,bidsDir,githubDir,fsDir);
setup_user('rania',bidsDir);

hRF_setting = 'glmsingle'; % can be: 'canonical', 'glmdenoise', 'glmsingle';

conditions = {'cardinaloblique', 'radialtangential'};

roiNames = {'V1_REmanual', 'pMT_REmanual', 'hMTcomplex_REmanual', 'pMST_REmanual', 'V2_REmanual', 'V3_REmanual'};

subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-0397', ...
    'sub-wlsubj121', 'sub-0442', ...
    'sub-wlsubj124', 'sub-wlsubj127', 'sub-0395', 'sub-0426', ...
    'sub-0250','sub-wlsubj123'}; % last 3 were originally omitted  'sub-0427',

derivativesFolder = fullfile(bidsDir, 'derivatives', sprintf('%sGLM',projectName), ...
    sprintf('hRF_%s', hRF_setting), 'group');
ttaveSavePath = fullfile(derivativesFolder, 'ttaveData');

plotValue = 'subtractOri';

if ~isfolder(ttaveSavePath)
    mkdir(ttaveSavePath)
end

for ci=1:numel(conditions)

    condName = conditions{ci};

    for ri=1:numel(roiNames)

        roiName = roiNames{ri};

        % create empty matrices
        base = nan(length(subjects),10);
        blank = nan(length(subjects),25);
        advMotion = nan(length(subjects),25);
        disadvMotion = nan(length(subjects),25);
        advStatic = nan(length(subjects),25);
        disadvStatic = nan(length(subjects),25);

        for si=1:numel(subjects)
        
            subjectname = subjects{si}
            
            path2ses = fullfile(bidsDir, 'derivatives', 'dgGLM', sprintf('hRF_%s', hRF_setting), subjectname);
        
            % get session
            allFolders = dir(path2ses);
            sesFolders = {allFolders([allFolders.isdir]).name};
            sesFolders = sesFolders(startsWith(sesFolders, 'ses-'));
            if length(sesFolders) > 1
                warning('More than one session found. ensure correct session is retrieved for this subject.')
            end
            ses = sesFolders{1};
        
            path2Output = fullfile(path2ses, ses, 'ttaveData', sprintf('ttaveSignal_%s_%s_%s_%s.mat', subjectname, roiName, condName, plotValue));
        
            load(path2Output)
        
            % save subject mean here
            base(si,:) = nanmean(ttaveOutput{1},1);
            blank(si,:) = nanmean(ttaveOutput{2},1);
            advMotion(si,:) = nanmean(ttaveOutput{3},1);
            disadvMotion(si,:) = nanmean(ttaveOutput{4},1);
            advStatic(si,:) = nanmean(ttaveOutput{5},1);
            disadvStatic(si,:) = nanmean(ttaveOutput{6},1);
        end

        grandttaveOutput = {base, blank, advMotion, disadvMotion, advStatic, disadvStatic};
        % ttaveTime is constat across subjects
        
        subj = 'allsubjects';
        filename = sprintf('ttaveSignal_%s_%s_%s_%s', subj, roiName, condName, plotValue);
        ttaveSave = fullfile(ttaveSavePath,filename);
        [legendLabels] = plot_ttave(grandttaveOutput, ttaveTime, subj, roiName, ttaveSave, condName, plotValue);

        % save the variables
        save(ttaveSave, 'ttaveOutput', 'ttaveTime', 'legendLabels', 'roiName', 'subj');

    end
end

close all;