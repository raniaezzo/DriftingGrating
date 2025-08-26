clear all; close all; clc
addpath(genpath(pwd));
bidsDir =  '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
setup_user('rania', bidsDir)
%%

% betamaps = nVertices x 13
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

%%

% how to reconstruct the signal?


% To run this script we recommend running Example1 first to create all the
% necessary outpouts from GLMsingle that are going to be reused here.

% Code here shows how to create a predicted timeseries for an example 
% voxel. This predicted timecourse is used in the GLM and is fitted to the
% data. Code below will also show the user how to concatante data and 
% designs across runs. For the purpose of this experiment we will use 
% the data from example 1, pick the voxel with highest variance explained 
% from the ON-OFF model and calculate the predicted timecourse. The
% predicted timecourse can be used for various things, e.g ilustrating 
% in the paper how well your prediction matches the fMRI time series.

% load data
%load('./data/nsdflocexampledataset.mat')
projectName = 'da'; % 'dg', or 'dgl' or 'da'
bidsDir =  '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
space = 'fsnative';
subj = 'sub-wlsubj124';
ses = 'ses-05';
run = 1:8; % CHANGE BACK TO 8!!!!
datacell = load_data(bidsDir,projectName,space,'.mgh',subj,ses,run) ;

stimdur = 3;
tr = 1;


% % load in modelOut and designSingle
% load(sprintf('/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/derivatives/%sGLM/hRF_glmsingle/%s/%s/modelOutput.mat', projectName, subj, ses));
% 
% temp = modelOut{1,4}.modelmd; % size = nVertices x 1 x 1 x nTrials
% trials = designSINGLE.stimorder; % size = 1 x nTrials (contains trial ID 1-13)
% 
% for ci=1:13
%         condSelect = designSINGLE.stimorder==ci;
%         betas = temp(:,:,:,condSelect);
%         newbetas = squeeze(betas);
%         
%         % trial mean
%         betamaps(:,ci) = mean(newbetas,2);
% end

%%
% designSingle  - each cell holds a design marix for each run. Number of
% columns in each cell corresponds to the number of unique images showed
% during the fMRI acqusition.

% data          - is the raw time series.

% results       -  GLMsingle outpupts for the 4 models.

%%
%Find a voxel with highest variance explained of the ON-OFF model
% [vexpl,ind] = max(results.R2(:)); % results
% % locate it on the brain slice
% %[ri,ci] = ind2sub(size(results{1}.R2),ind); % took out the braces
% %everywhere
% [ri,ci] = ind2sub(size(results.R2),ind);

% load results of TYPED model
%results = load('./example2outputs/GLMsingle/TYPED_FITHRF_GLMDENOISE_RR.mat');
results = load(sprintf('/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/derivatives/%sGLM/hRF_glmsingle/%s/%s/TYPED_FITHRF_GLMDENOISE_RR.mat', projectName, subj, ses));


% pick indices ri:
hSize = get_surfsize(subj); 
ri = getROIidxs(subj, 'V1', hSize); % %V1 hMTcomplex

% pick betas from the TYPED_FITHRF_GLMDENOISE_RR model to build the
% prediction.
%betas = flatten(results.modelmd(ri,ci,1,:)); %results{4}
betas = squeeze(results.modelmd(ri, 1, 1, :)); % Average over ri, then flatten
% predict time series
hrflibrary = getcanonicalhrflibrary(stimdur,tr)';  % timepoinpts x HRFs
hrfii = results.HRFindex(ri,1);  % HRF index % results{4}
meansignal = results.meanvol(ri,1);  % mean signal results{4}

temp = load(sprintf('/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/derivatives/%sGLM/hRF_glmsingle/%s/%s/DESIGNINFO.mat', projectName, subj, ses));
%temp = load('./example2outputs/GLMsingle/DESIGNINFO.mat');
designSINGLE2 = temp.designSINGLE;

nVoxels = length(meansignal);
nRuns = length(designSINGLE2);
nTRs = size(designSINGLE2{1},1);
pts = nan(nTRs,nVoxels,nRuns);

% per run
for ii=1:length(designSINGLE2)
    
    for rr=1:length(hrfii) % for each voxel
        design0 = conv2(designSINGLE2{ii},hrflibrary(:,hrfii(rr)));  % convolve HRF into design matrix
        design0 = design0(1:size(designSINGLE2{ii},1),:);        % truncate

        % not sure about this: adding the PC nuissange regressors
        %design0 = horzcat(design0, results.pcregressors{ii}(:,1:results.pcnum));

        betatemp = betas(rr,:)/100 .* meansignal(rr,:);   % betas (nVoxels x nTrials); meansignal (nVoxels x 1)

        % design0 (TRs=280 x nTrials=416) ; betatemp' (nTrials=416 x nVoxels)
        pts(:,rr,ii) = zeromean(design0*betatemp');    % weight and sum
    end
end

% average the voxels in the ROI then flatten the array
pts2 = squeeze(mean(pts,2));
pts2 = pts2(:);

%%
% % get the demeaned raw data
% data = horzcat(datacell{:}); % added
% %dataVOX = cellfun(@(pts) zeromean(flatten(pts(ri,ci,1,:))),data,'UniformOutput',0);
% %dataVOX = cat(2,dataVOX{:});
% dataVOX = zeromean(data(ri, :));
% 
% % data is mean across all voxels
% data = mean(data, 1)';
% data = data-mean(data); % no correction
% 
% figure; hold on;
% plot(pts2,'.-');
% plot(dataVOX,'r-');
% hold on
% %plot(data, 'c-')
% ylabel('%BOLD')
% xlabel('TRs (concatanated across runs)')
% set(gca,'Fontsize',15)
% %xlim([0 size(data{1},4)*length(data)])
% %xticks([0:size(data{1},4):size(data{1},4)*length(data)]);
% legend({'predicted timeseries';'raw timeseries'})
% legend box off
% set(gcf,'Position',[1000        1047        1113         290])
% % note: the influence of the glmdenoise nuisance regressors is not included here.
% % also, the influence of the polynomial drift regressors is not included.

%%

% try plotting the mean for ci

data = horzcat(datacell{:});
%dataVOX = zeromean(data(ri, :));
dataVOX = data(ri, :)-mean(data(ri,:),2);
dataPlot = mean(dataVOX,1);

% get the indices (over 4480 total TRs) that meet ci

if strcmp(projectName, 'dg')
    blank = getMeanarray([13], designSINGLE, designSINGLE2, pts2);
    blankRaw = getMeanarray([13], designSINGLE, designSINGLE2, dataPlot);
    cardinalmotion = getMeanarray(1:4, designSINGLE, designSINGLE2, pts2);
    cardinalmotionRaw = getMeanarray(1:4, designSINGLE, designSINGLE2, dataPlot);
    obliquemotion = getMeanarray(5:8, designSINGLE, designSINGLE2, pts2);
    obliquemotionRaw = getMeanarray(5:8, designSINGLE, designSINGLE2, dataPlot);
    cardinalorientation = getMeanarray(9:10, designSINGLE, designSINGLE2, pts2);
    cardinalorientationRaw = getMeanarray(9:10, designSINGLE, designSINGLE2, dataPlot);
    obliqueorientation = getMeanarray(11:12, designSINGLE, designSINGLE2, pts2);
    obliqueorientationRaw = getMeanarray(11:12, designSINGLE, designSINGLE2, dataPlot);
    
    figure
    plot(cardinalmotion-blank, 'b')
    hold on
    scatter(1:length(cardinalmotionRaw), cardinalmotionRaw-blankRaw, 'bo', 'filled')
    hold on
    plot(cardinalorientation-blank, 'b--')
    hold on
    plot(cardinalorientationRaw-blankRaw, 'bo')
    hold on
    plot(obliquemotion-blank, 'r')
    hold on
    scatter(1:length(obliquemotionRaw), obliquemotionRaw-blankRaw, 'ro', 'filled')
    hold on
    plot(obliqueorientation-blank, 'r--')
    hold on
    plot(obliqueorientationRaw-blankRaw, 'ro')
    xlim([1, 26])
    xticks = get(gca, 'XTick');            % Get current x-tick positions
    xticklabels = str2double(get(gca, 'XTickLabel')); % Convert tick labels to numbers
    new_xticklabels = xticklabels - 5;     % Subtract 5 from each label
    set(gca, 'XTickLabel', new_xticklabels); % Update x-axis labels
elseif strcmp(projectName, 'da')
    blank = getMeanarray([13], designSINGLE, designSINGLE2, pts2);
    blankRaw = getMeanarray([13], designSINGLE, designSINGLE2, dataPlot);
    radialmotion = getMeanarray([2,4], designSINGLE, designSINGLE2, pts2);
    radialmotionRaw = getMeanarray([2,4], designSINGLE, designSINGLE2, dataPlot);
    tangmotion = getMeanarray([1,3], designSINGLE, designSINGLE2, pts2);
    tangmotionRaw = getMeanarray([1,3], designSINGLE, designSINGLE2, dataPlot);
    radialorientation = getMeanarray(10, designSINGLE, designSINGLE2, pts2);
    radialorientationRaw = getMeanarray(10, designSINGLE, designSINGLE2, dataPlot);
    tangorientation = getMeanarray(9, designSINGLE, designSINGLE2, pts2);
    tangorientationRaw = getMeanarray(9, designSINGLE, designSINGLE2, dataPlot);
    
    figure
    plot(radialmotion-blank, 'b')
    %hold on
    %scatter(1:length(radialmotionRaw), radialmotionRaw-blankRaw, 'bo', 'filled')
    hold on
    plot(radialorientation-blank, 'b--')
    hold on
    %plot(radialorientationRaw-blankRaw, 'bo')
    %hold on
    plot(tangmotion-blank, 'r')
    hold on
    %scatter(1:length(tangmotionRaw), tangmotionRaw-blankRaw, 'ro', 'filled')
    %hold on
    plot(tangorientation-blank, 'r--')
    hold on
    %plot(tangorientationRaw-blankRaw, 'ro')
    xlim([1, 26])
    xticks = get(gca, 'XTick');            % Get current x-tick positions
    xticklabels = str2double(get(gca, 'XTickLabel')); % Convert tick labels to numbers
    new_xticklabels = xticklabels - 5;     % Subtract 5 from each label
    set(gca, 'XTickLabel', new_xticklabels); % Update x-axis labels
end


function output = getMeanarray(ci, designSINGLE, designSINGLE2, pts2)
    condSelect = ismember(designSINGLE.stimorder, ci);
    
    designMatrixNEW = cat(1, designSINGLE2{:}); % 2240 (TRs) x 416 (trials)
    conDesign = designMatrixNEW(:,condSelect);
    selectedTimepoint = any(conDesign,2);
    
    % tta:
    
    epochLength = 26;  % Number of timepoints per epoch
    selectedIdx = find(selectedTimepoint);  % Get indices where selectedTimepoint == 1
    
    % Preallocate output matrix (each row is an epoch)
    numEpochs = numel(selectedIdx);
    epochs = NaN(numEpochs, epochLength); % Using NaN in case some indices exceed bounds
    
    for i = 1:numEpochs
        idx = selectedIdx(i) - 5; % Start epoch at idx - 5
        epoch_end = idx + epochLength - 1; % Define the end of the epoch
    
        % Initialize row with NaNs
        epoch_data = NaN(1, epochLength);  
    
        % Ensure valid indices within bounds
        valid_start = max(idx, 1); % Ensure starting index is at least 1
        valid_end = min(epoch_end, length(pts2)); % Ensure end index does not exceed bounds
    
        % Copy available data into the correct part of the epoch
        start_offset = valid_start - idx + 1; % Offset to place data correctly in epoch_data
        epoch_data(start_offset:(start_offset + (valid_end - valid_start))) = pts2(valid_start:valid_end);
    
        % Store in epochs matrix
        epochs(i, :) = epoch_data;
    end
    
    output = nanmean(epochs,1);

end

