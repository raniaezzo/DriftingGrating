% call to plot observed vs. predicted (GLMsingle) run time series
%
% This is the run-timeseries counterpart to run_ttave.m: same subject/ROI
% setup and the same underlying data (raw BOLD -> %signal change, and the
% GLMsingle-predicted fit), but instead of extracting and averaging
% condition-triggered trial windows, it plots the full continuous run
% time series (observed vs. predicted) for each run, with the GLMsingle
% full-model R^2 (over the same vertex selection) printed on the plot.
%
% Runs for both projects (dg, da) and every subject shared between them.
% If derivedModelFit.mat doesn't exist yet for a given subject/project, it
% is generated automatically (ensureDerivedModelFit.m) -- this can take
% several minutes per subject/project the first time.
%
% Output: <bidsDir>/derivatives/<project>GLM/hRF_glmsingle/<subj>/<ses>/
%         ttaveData/ObservedandPredicted/Run<N>_ObservedandPredicted.{pdf,fig}

clear all; close all; clc;

% must be in DriftingGrating directory to run

% setup path
addpath(genpath(pwd));
bidsDir =  '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
githubDir = '~/Documents/GitHub';
fsDir = '/Applications/freesurfer/7.2.0';
addpath(genpath(fullfile(githubDir, 'atlasmgz')));
setup_user('rania', bidsDir)

projectSettingsBase = loadConfig(githubDir);

hRF_setting = 'glmsingle';
projects = {'dg', 'da'};

% ROI to analyze (matches a filename in general/jsons/ROIS_ALL.json).
roiName = 'V1';

% Polar angle wedges (deg) to include in the vertex selection, on top of
% the <roiName> + eccentricity + pRF R^2 criteria. Defaults to all 8
% canonical 45-deg-wide wedges (i.e. the whole visual field, no polar
% angle restriction) -- change this to a subset (e.g. [0] for the
% 0+/-22.5deg horizontal-right wedge) to restrict to specific polar angle
% location(s). Output is saved into a subfolder named after the ROI and
% selection: '<roiName>_pa_all' for the full canonical set, otherwise
% '<roiName>_pa_<angle>_<angle>...'.
polarAnglesToInclude = [0, 45, 90, 135, 180, 225, 270, 315];

% subjects shared between the two projects
dgSubs = dir(fullfile(bidsDir, 'derivatives', 'dgGLM', 'hRF_glmsingle', 'sub-*'));
dgSubs = {dgSubs([dgSubs.isdir]).name};
daSubs = dir(fullfile(bidsDir, 'derivatives', 'daGLM', 'hRF_glmsingle', 'sub-*'));
daSubs = {daSubs([daSubs.isdir]).name};
subjects = intersect(dgSubs, daSubs);

fprintf('%d subject(s) shared across dg/da\n', numel(subjects));

jsonParams = jsondecode(fileread('setup.json'));
stimdur_s = jsonParams.stimdur_s.Val;
tr_s = jsonParams.tr_s.Val;

%%
for pp = 1:numel(projects)
    projectName = projects{pp};

    for ss = 1:numel(subjects)
        subj = subjects{ss};
        fprintf('=== %s / %s ===\n', projectName, subj);

        % get ses name
        subjectDir = fullfile(bidsDir,'derivatives',strcat(projectName, 'GLM'), strcat('hRF_',hRF_setting), subj);
        contents = dir(subjectDir);
        sesNames = {};
        for i = 1:length(contents)
            if contents(i).isdir && startsWith(contents(i).name, 'ses-')
                sesNames{end+1} = contents(i).name;
            end
        end
        if length(sesNames) > 1
            error('Multiple session folders for %s/%s.', projectName, subj)
        elseif isempty(sesNames)
            warning('No session folder for %s/%s, skipping.', projectName, subj)
            continue
        end
        ses = sesNames{1};

        retFolder = 'prfvista_mov';
        derivativesFolder = fullfile(bidsDir, 'derivatives', sprintf('%sGLM',projectName), ...
            sprintf('hRF_%s', hRF_setting), subj, ses);

        % generate derivedModelFit.mat if it doesn't exist yet
        modelFitFile = ensureDerivedModelFit(bidsDir, projectName, subj, ses, hRF_setting, stimdur_s, tr_s);

        % load data and set up
        hSize = get_surfsize(subj);
        load(fullfile(derivativesFolder, 'rawInfo.mat'));

        run = 1:length(matrices_onset);
        datafiles = load_data(bidsDir,projectName,'fsnative','.mgh',subj,ses,run);
        modelfiles = loadDerivedModelFit(modelFitFile);

        % GLMsingle full-model R^2, same vertex order as datafiles/modelfiles
        modelOutS = load(fullfile(derivativesFolder, 'modelOutput.mat'), 'modelOut');
        glmR2 = double(modelOutS.modelOut{4}.R2(:));

        projectSettings = projectSettingsBase;
        projectSettings.projectName = projectName;
        projectSettings.subject = subj;
        projectSettings.bidsDir = bidsDir;
        projectSettings.retFolder = retFolder;
        projectSettings.ses = ses;

        % filter params (same criteria as ttave/dg_computeGain)
        projectSettings.polarAngleBinWidth = 45;
        projectSettings.minECC = 4;
        projectSettings.maxECC = 8;
        projectSettings.minVAREXP = .1;
        projectSettings.stimdur_s = stimdur_s;
        projectSettings.tr_s = tr_s;
        projectSettings.polarAnglesToInclude = polarAnglesToInclude;

        projectSettings.roiName = roiName;
        primaryROIvertices = getROIidxs(subj, projectSettings.roiName, hSize);
        surfaceROI = nan(sum(hSize),1);
        surfaceROI(primaryROIvertices) = 1;

        projectSettings.filteredPrfBins = retriveRetData(projectSettings);

        % subfolder name reflects the ROI + polar angle selection:
        % '<roiName>_pa_all' for the full canonical set (the default),
        % else '<roiName>_pa_<a>_<b>...'
        allCanonicalAngles = [0, 45, 90, 135, 180, 225, 270, 315];
        if isequal(sort(polarAnglesToInclude), allCanonicalAngles)
            paFolderName = sprintf('%s_pa_all', roiName);
        else
            angleStrs = arrayfun(@(a) num2str(a), sort(polarAnglesToInclude), 'UniformOutput', false);
            paFolderName = sprintf('%s_pa_%s', roiName, strjoin(angleStrs, '_'));
        end

        projectSettings.runTimeseriesSavePath = fullfile(derivativesFolder, 'ttaveData', 'ObservedandPredicted', paFolderName);

        runTimeseries_compute(matrices_onset, datafiles, modelfiles, surfaceROI, projectSettings, glmR2);

        close all;
    end
end
