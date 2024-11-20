% call to plot trial triggered average

clear all; close all; clc;

% must be in DriftingGrating directory to run

stimdur_s = 3;
tr_s = 1;
eventTRs_after = 20; % TRs to plot (ttave)
eventTRs_prior = 5;

% setup path
addpath(genpath(pwd));
projectName = 'dg'; % 'dg', or 'dgl' or 'da'
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

hRF_setting = 'glmsingle'; %
subj = 'sub-0427'; %'sub-0426';
ses = 'ses-01'; %'ses-01'; %'ses-nyu3t02'; %'ses-01';
retFolder = 'prfvista_mov';

derivativesFolder = fullfile(bidsDir, 'derivatives', sprintf('%sGLM',projectName), ...
    sprintf('hRF_%s', hRF_setting), subj, ses);

%%

% get ROI surface
% rois = {'V1', 'hMTcomplex'}; %, 'V2', 'V3', 'hV4', 'hMTcomplex', 'pMT', 'pMST'};
hSize = get_surfsize(subj); 

% for ri=1:numel(rois)
%         
%     roiname = rois{ri};
% 
%     if strcmp(roiname, 'V2') || strcmp(roiname, 'V3') % combine dorsal and ventral
%         lh_label1 = read_label(subj, sprintf('retinotopy_RE/lh.%sv_REmanual', roiname));
%         lh_label2 = read_label(subj, sprintf('retinotopy_RE/lh.%sd_REmanual', roiname));
%         rh_label1 = read_label(subj, sprintf('retinotopy_RE/rh.%sv_REmanual', roiname));
%         rh_label2 = read_label(subj, sprintf('retinotopy_RE/rh.%sd_REmanual', roiname));
% 
%         label_idx = [lh_label1(:,1)+1 ; lh_label2(:,1)+1; rh_label1(:,1)+hSize(1)+1; rh_label2(:,1)+hSize(1)+1];
%     else
%         % columns: vertex index, x, y, z, value
%         lh_label = read_label(subj, sprintf('retinotopy_RE/lh.%s_REmanual', roiname));
%         rh_label = read_label(subj, sprintf('retinotopy_RE/rh.%s_REmanual', roiname));
% 
%         % plus one because matlab starts from 1 not 0
%         label_idx = [lh_label(:,1)+1 ; rh_label(:,1)+hSize(1)+1];
%     end
% 
% end

ttaveSavePath = fullfile(derivativesFolder, 'ttaveData');

if ~isfolder(ttaveSavePath)
    mkdir(ttaveSavePath)
end

% load in design matrix
load(fullfile(derivativesFolder, 'rawInfo.mat'));

% load raw signal data (in wpToolbox) - convert to .mgh if not already
run = 1:length(matrices_onset);
datafiles = load_data(bidsDir,projectName,'fsnative','.mgh',subj,ses,run) ;

roiNames = {'V1_REmanual', 'V2_REmanual', 'V3_REmanual', 'pMST_REmanual', 'hMTcomplex_REmanual', 'pMT_REmanual'};

% filter params
polarAngleBinWidth = 45; %45; % degrees
minECC = 1; %0.2; 
maxECC = 15;
minVAREXP = .1; %0.25;

conditions = {'cardinaloblique', 'radialtangential'};
plotValue = 'subtractOri'; % make 'rawValues' (motion and orientation without subtraction)

%%
for ri=1:length(roiNames)

    % get ROI indices
    roiName = roiNames{ri};
    primaryROIvertices = getROIidxs(subj, roiName, hSize); % the name of Middle Temporal Cortex ROI
    surfaceROI = nan(sum(hSize),1);
    surfaceROI(primaryROIvertices) = 1;

    % get VF indices binned into 8 categories (after filtering for: ecc, r^2)
    filteredPrfBins = retriveRetData(bidsDir, retFolder, subj, polarAngleBinWidth, minECC, maxECC, minVAREXP);

    for ci=1:2

        condName = conditions{ci};

        filename = sprintf('ttaveSignal_%s_%s_%s_%s', subj, roiName, condName, plotValue);
        ttaveSave = fullfile(ttaveSavePath,filename);

        if strcmp(condName, 'cardinaloblique')
            % compute and plot ttave
            [ttaveOutput, ttaveTime] = ttave_compute(matrices_onset,datafiles, surfaceROI, tr_s, stimdur_s, subj, eventTRs_prior, eventTRs_after);
        elseif strcmp(condName, 'radialtangential')
            [ttaveOutput, ttaveTime] = ttave_compute(matrices_onset,datafiles, surfaceROI, tr_s, stimdur_s, subj, eventTRs_prior, eventTRs_after, filteredPrfBins);
        end
            
        [legendLabels] = plot_ttave(ttaveOutput, ttaveTime, subj, roiName, ttaveSave, condName, plotValue);
    
        % save the variables
        save(ttaveSave, 'ttaveOutput', 'ttaveTime', 'legendLabels', 'roiName', 'subj');
    end
end

close all;