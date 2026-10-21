% plotting GLMsingle model fits

clear all; close all; clc;

%% SET UP

% set fundamental vars to run a single subject / protocol
projectName = 'dg'; % 'dg', or 'dgl' or 'da'
subj = 'sub-0037'; %'sub-0426';
ses = 'ses-01'; %'ses-01'; %'ses-nyu3t02'; %'ses-01'; <-- can I make this not manual?
space = 'fsnative';
bidsDir =  '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';

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

%%

load(fullfile(derivativesFolder,'modelOutput.mat'))

%load(fullfile(derivativesFolder,'rawInfo.mat'))

results.(char('allevents')) = modelOut{1,4};

%% load V1 ROI, and Polar angle

hemis = {'lh'; 'rh'};
roiname = 'V1';

retDir = dir(fullfile(bidsDir, 'derivatives', 'prfvista_mov', subj, '**/stimfiles.mat'));
retDir = retDir.folder;

for hi=1:numel(hemis)
    hemi = hemis{hi};
    ret_moving.(sprintf('%s_pa', hemi)) = MRIread(fullfile(retDir, sprintf('%s.angle_adj.mgz', hemi)));
    ret_moving.(sprintf('%s_ecc', hemi)) = MRIread(fullfile(retDir, sprintf('%s.eccen.mgz', hemi)));
    ret_moving.(sprintf('%s_vexp', hemi)) = MRIread(fullfile(retDir, sprintf('%s.vexpl.mgz', hemi)));
end

retData = [ret_moving.lh_pa.vol, ret_moving.rh_pa.vol ; ...
        ret_moving.lh_ecc.vol, ret_moving.rh_ecc.vol ; ...
        ret_moving.lh_vexp.vol, ret_moving.rh_vexp.vol]';

hSize = [size(ret_moving.lh_pa.vol,2), size(ret_moving.rh_pa.vol,2)];
size_combinedSurfs = sum(hSize);

lh_label = read_label(subj, sprintf('retinotopy_RE/lh.%s_REmanual', roiname));
rh_label = read_label(subj, sprintf('retinotopy_RE/rh.%s_REmanual', roiname));

% plus one because matlab starts from 1 not 0
label_idx = [lh_label(:,1)+1 ; rh_label(:,1)+hSize(1)+1];

roi_mask = zeros(size_combinedSurfs, 1);
roi_mask(label_idx) = 1; 

% visual field bin
binCenters = 0:45:315;    % Bin labels (degrees)
angles = map_theta(retData(:,1));

% Circular distance from each angle to each bin center
circDist = abs(mod(angles - binCenters + 180, 360) - 180);

% Closest bin and its distance
[minDist, idx] = min(circDist, [], 2);

% Initialize as NaN
angle_bin = nan(height(angles),1);

% Assign only if within +/-22.5°
valid = minDist <= 22.5;
angle_bin(valid) = binCenters(idx(valid));

eccentricity = retData(:,2);
rsquared = retData(:,3);

% included in analysis? 
includeVoxels = ...
    ~isnan(angle_bin) & ...
    eccentricity >= 4 & eccentricity <= 8 & ...
    rsquared >= 0.1 & roi_mask;

%%

unique_hrfs = unique(modelOut{4}.HRFindex(includeVoxels));
% modelfit_roi = zeros(total_timepoints, includeVoxels); % pre-allocate
col_offset = 0;
numtrials = size(modelOut{4}.modelmd, 4);
hrflibrary = getcanonicalhrflibrary(stimdur_s, tr_s);

n_runs = numel(designSINGLE.designSINGLE);
n_vox = sum(includeVoxels);

numtimepoints = zeros(1, n_runs);

for n = 1:n_runs
    numtimepoints(n) = size(designSINGLE.designSINGLE{n}, 1);
end

modelfit = cell(1, n_runs);

n_vertices_total = numel(modelOut{4}.HRFindex);
for r = 1:n_runs
    modelfit{r} = nan(n_vertices_total, numtimepoints(r));
end

% figure
% plot(hrflibrary')
% xlabel('Sample')
% ylabel('Value')

for hh = unique_hrfs'

    % voxels in the ROI that use this HRF    
    vox_mask = (modelOut{4}.HRFindex == hh); %& includeVoxels;
    hrf = hrflibrary(hh,:)'; % hrf = hrflibrary(:, hh); % 
    betas_this_hrf = modelOut{4}.modelmd(repmat(vox_mask, [1 1 1 numtrials]));
    betas_this_hrf = reshape(betas_this_hrf, [], numtrials); % nvox x numtrials
    
    betas = reshape(modelOut{4}.modelmd .* vox_mask, [], numtrials);

    % call GLMpredictresponses for this HRF group
    mf = GLMpredictresponses({hrf, betas}, ...
        designSINGLE.designSINGLE, tr_s, numtimepoints, 1); % ... extract and accumulate

    %df = GLMestimatesingletrial(designSINGLE.designSINGLE,datafiles,stimdur_s,tr_s,outputdir,opt);

    for r = 1:numel(mf)
    
        valid = ~all(mf{r} == 0, 2);
    
        modelfit{r}(valid,:) = mf{r}(valid,:);
    
    end

end

save(fullfile(derivativesFolder, 'derivedModelFit.mat'), 'modelfit', '-v7.3');

disp('mean R^2 for valid voxels:')
nanmean(modelOut{4}.R2(includeVoxels))
disp('meidan R^2 for valid voxels:')
nanmedian(modelOut{4}.R2(includeVoxels))

%%

% tt = load(fullfile(derivativesFolder,'rawInfo.mat'));
% 
% %%
% run = 6;
% selectVoxels = modelfit_roi{run}(includeVoxels,:);
% meanData = mean(datafiles{run}(includeVoxels,:),2);
% selectVoxelsData = ((datafiles{run}(includeVoxels,:)-meanData)./(meanData)) *100;
% 
% figure
% plot(mean(selectVoxels), 'LineWidth',1)
% hold on
% plot(mean(selectVoxelsData), 'LineWidth',1)
% legend({'prediction', 'raw data'})