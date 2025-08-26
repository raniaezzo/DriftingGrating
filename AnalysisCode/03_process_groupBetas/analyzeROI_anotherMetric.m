% compute myelin, etc.


clc; clear all; close all;

% set up
addpath(genpath(pwd));
projectName = 'dg';
%bidsDir =  '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
bidsDir = '/Volumes/server/Projects/Project_dg/data_bids/';
githubDir = '~/Documents/GitHub';
fsDir = '/Applications/freesurfer'; %/7.2.0';
%addpath(genpath(fullfile(githubDir, 'retinotopy-nyuad')));
addpath(genpath(fullfile(githubDir, 'wpToolbox')));
addpath(genpath(fullfile(githubDir, 'atlasmgz')));
%setup_user(projectName,bidsDir,githubDir,fsDir);
setup_user('rania',bidsDir);
masksFolder = fullfile(bidsDir, 'derivatives', 'masks');
retinotopyMode = 'moving'; % of 'flicker' to use for analysis

hemis = {'lh'; 'rh'};
%rois = {'V1', 'V2v', 'V3v', 'V2d', 'V3d', 'hV4', 'hMTcomplex', 'pMT', 'pMST'};
rois = {'V1', 'V2', 'V3', 'hV4', 'hMTcomplex', 'pMT', 'pMST'};

% replaced below to pairwise (sep): 'cardmVcards', 'oblmVobls'
subjects = {'sub-0037', 'sub-0201',  'sub-0255', 'sub-0397', ...
    'sub-wlsubj123', 'sub-wlsubj124', 'sub-0395', 'sub-0426', ...
    'sub-0250'}; % last 3 were originally omitted 

% 'sub-0442', 'sub-wlsubj121', ...'sub-wlsubj127', 'sub-0427', 

% update contrast names:
% contrastnames = [contrastnames, indiv_contrastnames];

polarAngles = [0, 45, 90, 135, 180, -135, -90, -45];
polarAngleBinWidth = 45; %45; % degrees
minECC = 1; %0.2; 
maxECC = 15;
minVAREXP = .1; %0.25;

stimulus = "MYELIN";
numConds = 1;
if strcmp(stimulus, "DG")
    contrastnames = {'allmVblank'}; %{'allmVallb'};
elseif strcmp(stimulus, "MYELIN")
elseif strcmp(stimulus, "TRANSPARENT")
end

% initialize a matrix for contrasts, rois, subjects
meanBOLD = nan(numConds, length(rois), length(subjects));
medianBOLD = nan(numConds, length(rois), length(subjects));

meanBOLDpa = nan(numConds, length(polarAngles), length(rois), length(subjects));
medianBOLDpa = nan(numConds,length(polarAngles), length(rois), length(subjects));
countBin = nan(length(polarAngles), length(rois), length(subjects));

nPreallVoxel = 2000;
allvoxelsBOLDpa = nan(numConds, length(polarAngles), length(rois), nPreallVoxel, length(subjects));
allparamsBOLDpa = nan(numConds, length(polarAngles), length(rois), nPreallVoxel, 4, length(subjects)); % pa, eccen, r^2
%%

for si=1:numel(subjects)

    subjectname = subjects{si}

    subj_FSfolder = fullfile(bidsDir, 'derivatives', 'freesurfer', subjectname);
    
    if strcmp(stimulus, "DG")
        glmResultsfolder = fullfile(bidsDir, 'derivatives', 'dgGLM', subjectname);
        glmFilelist = dir(fullfile(glmResultsfolder, '**', 'results.mat'));
        betaResults = load(fullfile(glmFilelist.folder, glmFilelist.name));
        betaResults = betaResults.results;
    elseif strcmp(stimulus, "MYELIN")
        betaResults = [];
        glmResultsfolder = fullfile(bidsDir, 'derivatives', 'T1MapMyelin', subjectname);
    elseif strcmp(stimulus, "TRANSPARENT")
        betaResults = [];
        glmResultsfolder = fullfile(bidsDir, 'derivatives', 'transparent', subjectname);
    end
    
    % load in pRF data for that subject (flicker)
%     flickerRetDir = dir(fullfile(bidsDir, 'derivatives', 'prfvista', subjectname, '**/stimfiles.mat'));
%     flickerRetDir = flickerRetDir.folder;
    movingRetDir = dir(fullfile(bidsDir, 'derivatives', 'prfvista_mov', subjectname, '**/stimfiles.mat'));
    movingRetDir = movingRetDir.folder;
    
    for hi=1:numel(hemis)
        hemi = hemis{hi};
        ret_moving.(sprintf('%s_pa', hemi)) = MRIread(fullfile(movingRetDir, sprintf('%s.angle_adj.mgz', hemi)));
        ret_moving.(sprintf('%s_ecc', hemi)) = MRIread(fullfile(movingRetDir, sprintf('%s.eccen.mgz', hemi)));
        ret_moving.(sprintf('%s_vexp', hemi)) = MRIread(fullfile(movingRetDir, sprintf('%s.vexpl.mgz', hemi)));
        ret_moving.(sprintf('%s_sigma', hemi)) = MRIread(fullfile(movingRetDir, sprintf('%s.sigma.mgz', hemi)));

%         for cc=1:numel(contrastnames)
%             cond = contrastnames{cc};
%             mot_v_stat.(sprintf('%s_%s', hemi,cond)) = MRIread(fullfile(glmFilelist.folder, sprintf('%s.%s.mgz', hemi, cond)));
%         end

    end
    
    if strcmp(stimulus, "MYELIN")
        betaResults_lh = MRIread(fullfile(glmResultsfolder, sprintf('%s.myelin0.5.mgz', 'lh')));
        betaResults_rh = MRIread(fullfile(glmResultsfolder, sprintf('%s.myelin0.5.mgz', 'rh')));
        betaResults = [betaResults_lh.vol, betaResults_rh.vol];
    elseif strcmp(stimulus, "TRANSPARENT")
        betaResults_lh = MRIread(fullfile(glmResultsfolder, sprintf('%s.oppo3.mgz', 'lh')));
        betaResults_rh = MRIread(fullfile(glmResultsfolder, sprintf('%s.oppo3.mgz', 'rh')));
        betaResults = [betaResults_lh.vol', betaResults_rh.vol'];
    end
    
    retMovingData = [ret_moving.lh_pa.vol, ret_moving.rh_pa.vol ; ...
        ret_moving.lh_ecc.vol, ret_moving.rh_ecc.vol ; ...
        ret_moving.lh_vexp.vol, ret_moving.rh_vexp.vol; ...
        ret_moving.lh_sigma.vol, ret_moving.rh_sigma.vol];


    if strcmp(retinotopyMode, 'moving')
        retData = retMovingData;
    elseif strcmp(retinotopyMode, 'flicker')
        retData = retFlickerData;
    end

    % get lh, rh sizes:
    hSize = get_surfsize(subjectname); 
    
    % create matrix of mean values per ROI (V1, V2, V3, MT+):
    
    for ci=1:numel(numConds)
    
        if strcmp(stimulus, "DG")
            contrastname = contrastnames{ci};
        end
    
        for ri=1:numel(rois)
        
            roiname = rois{ri};
        
            if strcmp(roiname, 'V2') || strcmp(roiname, 'V3') % combine dorsal and ventral
                lh_label1 = read_label(subjectname, sprintf('retinotopy_RE/lh.%sv_REmanual', roiname));
                lh_label2 = read_label(subjectname, sprintf('retinotopy_RE/lh.%sd_REmanual', roiname));
                rh_label1 = read_label(subjectname, sprintf('retinotopy_RE/rh.%sv_REmanual', roiname));
                rh_label2 = read_label(subjectname, sprintf('retinotopy_RE/rh.%sd_REmanual', roiname));

                label_idx = [lh_label1(:,1)+1 ; lh_label2(:,1)+1; rh_label1(:,1)+hSize(1)+1; rh_label2(:,1)+hSize(1)+1];
            else
                % columns: vertex index, x, y, z, value
                lh_label = read_label(subjectname, sprintf('retinotopy_RE/lh.%s_REmanual', roiname));
                rh_label = read_label(subjectname, sprintf('retinotopy_RE/rh.%s_REmanual', roiname));
        
                %rh_label = read_label(subjectname, 'Glasser2016/rh.MT');
        
                % plus one because matlab starts from 1 not 0
                label_idx = [lh_label(:,1)+1 ; rh_label(:,1)+hSize(1)+1];
            end
            
            
            if strcmp(stimulus, "DG")
                % load in results (units are in PSC)
                propername = strrep(contrastname, '_v_','V');
                currBold = betaResults.contrasts.(propername);
            elseif strcmp(stimulus, "MYELIN") || strcmp(stimulus, "TRANSPARENT")
                currBold = betaResults';
            end
                
            % get the currBOLD
            meanBOLD(ci, ri, si) = mean(currBold(label_idx));
            medianBOLD(ci, ri, si) = median(currBold(label_idx));

            % make mask (mgz) and save 
            mask = zeros(size(currBold));
            mask(label_idx) = 1;
            writeMGZfile(bidsDir, subjectname, '', mask, masksFolder, sprintf('%smask',roiname))
    
            % loop through pa bins:
            for pa=1:numel(polarAngles)
                pangle = polarAngles(pa);

                if pangle == -180 % cannot be -180
                    pangle = 180;
                end
                
                % Calculate the lower and upper bounds of the angular bin
                lower_bound = mod(pangle - polarAngleBinWidth/2 + 180, 360) - 180;
                upper_bound = mod(pangle + polarAngleBinWidth/2 + 180, 360) - 180;
                
                % Find values within the specified angular bin
                if lower_bound>upper_bound
                    currPA_bool = (retData(1,:) >= lower_bound) | (retData(1,:) <= upper_bound);
                else
                    currPA_bool = (retData(1,:) >= lower_bound) & (retData(1,:) <= upper_bound);
                end
                
                % log the magnitude away from center of the wedge bin:
                currPA_mag = abs(retData(1,:) - pangle);
                    
                % restrict eccentricities to be within a range
                currECC_bool = (retData(2,:) >= minECC) & (retData(2,:) <= maxECC);

                % restrict variance explained to be a certain minimum
                currVAREXP_bool = (retData(3,:) >= minVAREXP);

                % if it meets the 3 criteria
                filterRet = currPA_bool+currECC_bool+currVAREXP_bool == 3;
                validIdx = find(filterRet>0 );

                % find intersection between validIdx and label_idx
                finalIdx = intersect(label_idx, validIdx);

                % isolate ROI activity (with label)
                % mean per PA bin
                meanBOLDpa(ci,pa,ri,si) = mean(currBold(finalIdx));
                medianBOLDpa(ci,pa,ri,si) = median(currBold(finalIdx));
                
                % beta values (all voxels) - for histograms
                allvoxelsBOLDpa(ci,pa,ri,1:length(finalIdx),si) = currBold(finalIdx);
                allparamsBOLDpa(ci,pa,ri,1:length(finalIdx),1,si) = retData(1,finalIdx); % pa
                allparamsBOLDpa(ci,pa,ri,1:length(finalIdx),2,si) = retData(2,finalIdx); % eccen
                allparamsBOLDpa(ci,pa,ri,1:length(finalIdx),3,si) = retData(3,finalIdx); % r^2
                allparamsBOLDpa(ci,pa,ri,1:length(finalIdx),4,si) = currPA_mag(finalIdx); % pa magnitude
                

            end

        end
    end

end

%%

titlenames = {'v1', 'mt', 'mst'};
counter = 1;

n = length(subjects);
ones_vector = ones(n, 1);
jittered_vector = ones_vector + (rand(n, 1) - 0.5) * 0.1;

% figure
% for rr=[1,6,7]
%     titlename = titlenames{counter};
% 
%     %subplot(1,3, counter)
%     meanBOLDpanew = squeeze(meanBOLDpa);
% 
%     roi_meanBOLDpa = squeeze(meanBOLDpanew(:,rr,:));
%     
%     grandMeanPA = nanmean(roi_meanBOLDpa,1);
%     
%     scatter(counter+jittered_vector-1, grandMeanPA, 100, 'k', 'filled', 'MarkerFaceAlpha', 0.25);
%     
%     counter = counter+1;
%     hold on
% end

figure
x_values = [];  % To store x coordinates
y_values = [];  % To store y coordinates

for rr = [1, 6, 7]
    titlename = titlenames{counter};

    % Get the mean BOLD panel for the current ROI
    meanBOLDpanew = squeeze(meanBOLDpa);
    roi_meanBOLDpa = squeeze(meanBOLDpanew(:, rr, :));
    
    % Calculate the grand mean for the current ROI
    grandMeanPA = nanmean(roi_meanBOLDpa, 1);
    
    % Calculate x position for the current scatter
    x_pos = counter + jittered_vector - 1;
    
    % Scatter plot the points
    scatter(x_pos, grandMeanPA, 100, 'k', 'filled', 'MarkerFaceAlpha', 0.25);
    
    % Store the x and y values
    x_values = [x_values, x_pos(:)]; % Ensuring column vectors
    y_values = [y_values, grandMeanPA(:)];
    
    counter = counter + 1;
    hold on
end

% Plot the lines connecting the points across the ROIs
for i = 1:size(x_values, 1)
    plot(x_values(i, :), y_values(i, :), 'k-', 'LineWidth', 1);
end

hold on
mm = mean(y_values);
scatter(1, mm(1), 200, 'k', 'filled')
hold on
scatter(2, mm(2), 200, 'k', 'filled')
hold on
scatter(3, mm(3), 200, 'k', 'filled')
hold on



ylabel('R1 (1/T1)', 'FontSize', 24, 'FontName', 'Arial') % 'R1 (1/T1)' opponency (% bold)
xlim([0 4])
xticks([1 2 3]);
xticklabels({'V1', 'MT', 'MST'}); 
axis square
% Increase font size of axis values (tick labels)
set(gca, 'FontSize', 20, 'FontName', 'Arial');  % Here, 20 is the font size

ttest(y_values(:, 1), y_values(:, 2))