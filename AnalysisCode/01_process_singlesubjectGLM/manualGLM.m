%% manual GLM

%load('hRF.mat')
% replace with fitted hRF
%hRF = results.allevents.modelmd{1,1}; % use the same hRF

nRuns = length(matrices);
[nTRs, nRegs] = size(matrices{1});
[nVoxels, ~] = size(datafiles{1});

drift = 1:nTRs;
const = ones(nTRs,1); 
nEvents = 13; % 13 events, the rest are noise regressors

%psc = nan(nVoxels, nRegs+2, nRuns); % +1 for global mean, and adding ISI and start/end clip

% should I concatenate all runs and then run GLM or fit per run and
% average? for now averaging each:
y = {}; yHat = {}; betas = {}; varexp = {};

for ri=1:nRuns
    %matrices_new = vertcat(matrices{:});
    %datafiles_new = horzcat(datafiles{:});
    design = [matrices{ri} const drift'];
    data = datafiles{ri};
    
    [y, yHat, betas, varexp] = calculateBeta(data, design, nEvents);

end

% % subtract out estimated noise:
% noise = betas(:,14:end)*[design(:,14:end)]';
% y2 = y-noise;
% 
% yHat2 = yHat-noise;
% 
% meanPerVoxel = betas(:,end-1)*[design(:,end-1)]';
% meanPerVoxel = mean(abs(meanPerVoxel),2);
% 
% y2 = ((y2-meanPerVoxel)./meanPerVoxel))*100;
% yHat2 = (yHat2./meanPerVoxel)*100;

% % abs (to prevent flipping sign) - done in GLMestimatemodel
% mean_overT = mean(abs(data),2); 
% 
% % doing PSC here - subtract 1 to remove baseline signal
% y = ((data./mean_overT)-1)*100;

%%
% to check a roi:

[~, row_index] = max(varexp(:, 1)); % row_index = 2846

% plot ROI
figure
roi_mean = nanmean(y.*surfaceROI);
roiHat_mean = nanmean(yHat.*surfaceROI);
roi_betas = nanmean(betas.*surfaceROI, 1);
subplot(3,1,1)
plot(1:280, roi_mean, 'k','LineWidth', 2)
hold on
plot(1:280, roiHat_mean, 'b','LineWidth', 2)
ylabel('PSC')
xlabel('Time (sec)')
xlim([0 280])
legend({'y','yHat'})
ax = gca;
ax.FontSize = 14;
eventMat = design(:,1:12);
load('hRF.mat')
subplot(3,1,2)
for i = 1:size(eventMat, 2)
    temp = conv(eventMat(:, i), hRF);
    % Truncate the convolution result to match the length of the original data
    pred_hRF(:, i) = temp(1:size(eventMat,1));
    hold on
    if i<9
        plot(pred_hRF(:,i)*roi_betas(i), 'b','LineWidth', 2)
    else
        plot(pred_hRF(:,i)*roi_betas(i), 'r','LineWidth', 2)
    end
end
xlim([0 280])
hold off
subplot(3,1,3)
mvs = [sum(design(:,1:8),2), sum(design(:,9:12),2)]';
imagesc(mvs)
xlim([0 280])
ax = gca;
%motionHat = (pred_hRF(1:8)' * betas(1:8)');

figure
histogram(betas(:,1:8).*surfaceROI)
hold on
histogram(betas(:,9:12).*surfaceROI)

motionBetamean = mean(betas(:,1:8),2);
staticBetamean = mean(betas(:,9:12),2);
roimotionBetas = motionBetamean.*surfaceROI;
roistaticBetas = staticBetamean.*surfaceROI;
figure
histogram(roimotionBetas, 'FaceColor', 'blue', 'FaceAlpha',0.3)
hold on
xline(nanmean(roimotionBetas), 'b', 'LineWidth', 2)
hold on
histogram(roistaticBetas, 'FaceColor', 'red', 'FaceAlpha',0.3)
hold on
xline(nanmean(roistaticBetas), 'r', 'LineWidth', 2)
title('Motion vs Static betas')


motionCBetamean = mean(betas(:,1:4),2);
motionOBetamean = mean(betas(:,5:8),2);
roimotionCBetas = motionCBetamean.*surfaceROI;
roimotionOBetas = motionOBetamean.*surfaceROI;
figure
histogram(roimotionCBetas, 'FaceColor', 'green', 'FaceAlpha',0.3)
hold on
xline(nanmean(roimotionCBetas), 'g', 'LineWidth', 2)
hold on
histogram(roimotionOBetas, 'FaceColor', [.5 .5 .5], 'FaceAlpha',0.3)
hold on
xline(nanmean(roimotionOBetas), 'k', 'LineWidth', 2)
title('Cardinal v Oblique Motion betas')

%%

% %% Visualize
% 
% bidsDir = '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/'; %'/Volumes/Vision/MRI/recon-bank';
% subj = 'sub-wlsubj123';
% 
% view_fv(subj,bidsDir,'cardm_v_cards', 'oblm_v_obls')


% %%
% 
% % mask with MT+ ROI ?
% % do this later
% 
% %% visualize R2 for each roi
% close all;
% 
% returndir = pwd;
% 
% % Generate wang and glasser atlases
% glasser_dir = fullfile(bidsDir, 'derivatives', 'freesurfer', subj, 'label', 'Glasser2016');
% if ~isfolder(glasser_dir)
%     cd(fullfile(githubDir, 'atlasmgz'))
%     system(['sh createAtlasLabels.sh ' extractAfter(subj,'-') ' ' bidsDir])
% end
% 
% cd(returndir);
% 
% %'cardmVcards','oblmVobls','allmValls', 
% 
% visualize(results.contrasts, subj, {'allmValls'})
% % visualize(results,subj,rois);
% 
% 
% %view_fv(subj,bidsDir,'mt+2')