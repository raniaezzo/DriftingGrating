clc; clear all; close all;

% set up
addpath(genpath(pwd));
projectName = 'dg';
bidsDir =  '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
%bidsDir = '/Volumes/server/Projects/Project_dg/data_bids/';
githubDir = '~/Documents/GitHub';
fsDir = '/Applications/freesurfer'; %/7.2.0';
%addpath(genpath(fullfile(githubDir, 'retinotopy-nyuad')));
%addpath(genpath(fullfile(githubDir, 'wpToolbox')));
addpath(genpath(fullfile(githubDir, 'atlasmgz')));
%setup_user(projectName,bidsDir,githubDir,fsDir);
setup_user('rania',bidsDir);
masksFolder = fullfile(bidsDir, 'derivatives', 'masks');
retinotopyMode = 'moving'; % of 'flicker' to use for analysis

hRF_setting = 'glmsingle'; % can be: 'canonical', 'glmdenoise', 'glmsingle';

hemis = {'lh'; 'rh'};

projectSettings = loadConfig(githubDir);

rois = projectSettings.rois;
axes_limits = projectSettings.axes_limits;
pairaxes_limits = projectSettings.pairaxes_limits;
pairaxes_PAew_limits = projectSettings.pairaxes_PAew_limits;
colors_data = projectSettings.colors_data;
contrasts_dict = projectSettings.contrasts_dict;

% contrastnames = {contrasts_dict.contrasts.(strcat(projectName, '_contrast_name'))};
% for now, only use DG names -- they apply to both DG and DA b/c directions
% defined in absolute reference frame
%contrastnames = {contrasts_dict.contrasts.(strcat(projectName, '_contrast_name'))};
contrastnames = {contrasts_dict.contrasts.('dg_contrast_name')};

colors = colors_data.conditions.(projectName);

% % replaced below to pairwise (sep): 'cardmVcards', 'oblmVobls'
% contrastnames = {'cardMsep', 'oblMsep', 'allmValls', ...
%     'allsVblank', 'allmVblank', 'cardmVblank', 'oblmVblank', ...
%     'm0_v_s90','m90_v_s0','m180_v_s90','m270_v_s0', ...
%     'm45_v_s135','m135_v_s45','m225_v_s135','m315_v_s45', 'cardsVblank', 'oblsVblank', ...
%      'm0_v_b', 'm180_v_b', 'm90_v_b', 'm270_v_b', 'm45_v_b', 'm225_v_b', 'm135_v_b', ...
%      'm315_v_b', 's0_v_b', 's90_v_b', 's45_v_b', 's135_v_b'};

if strcmp(projectName, 'dg')
    subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', ...
        'sub-wlsubj124', 'sub-0395', 'sub-0426', 'sub-0250', ...
        'sub-0442', 'sub-wlsubj121', 'sub-wlsubj127',  'sub-0397', ...
        'sub-0427'}; % last 3 were originally omitted 
elseif strcmp(projectName, 'da')
    subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
        'sub-0395', 'sub-0426', 'sub-0250'};
end

indiv_contrastnames = {'m0_v_s90','m90_v_s0','m180_v_s90','m270_v_s0', ...
    'm45_v_s135','m135_v_s45','m225_v_s135','m315_v_s45'};
% Return contrast names as vector of motion dirs
pattern = 'm(\d+)_';
% Apply regular expression to each element in the cell array
matches = cellfun(@(str) regexp(str, pattern, 'tokens'), indiv_contrastnames, 'UniformOutput', false);
% Extract the numeric part and convert to a vector of integers
motiondirs = cellfun(@(match) str2double(match{1}), matches);

% update contrast names:
% contrastnames = [contrastnames, indiv_contrastnames];

polarAngles = [0, 45, 90, 135, 180, -135, -90, -45];
polarAngleBinWidth = 45; %45; % degrees
minECC = 4; %1; 
maxECC = 8; %15;
minVAREXP = .1; %0.25;

% initialize a matrix for contrasts, rois, subjects
meanBOLD = nan(length(contrastnames), length(rois), length(subjects));
medianBOLD = nan(length(contrastnames), length(rois), length(subjects));

meanBOLDpa = nan(length(contrastnames), length(polarAngles), length(rois), length(subjects));
medianBOLDpa = nan(length(contrastnames),length(polarAngles), length(rois), length(subjects));
countBin = nan(length(polarAngles), length(rois), length(subjects));

nPreallVoxel = 2000;
allvoxelsBOLDpa = nan(length(contrastnames), length(polarAngles), length(rois), nPreallVoxel, length(subjects));
allparamsBOLDpa = nan(length(contrastnames), length(polarAngles), length(rois), nPreallVoxel, 4, length(subjects)); % pa, eccen, r^2
%%

for si=1:numel(subjects)

    subjectname = subjects{si}

    subj_FSfolder = fullfile(bidsDir, 'derivatives', 'freesurfer', subjectname);
    %setenv('SUBJECT_NAME', subjectname)
    
    glmResultsfolder = fullfile(bidsDir, 'derivatives', strcat(projectName, 'GLM'), strcat('hRF_', hRF_setting), subjectname);
    glmFilelist = dir(fullfile(glmResultsfolder, '**', 'results.mat'));
    [~,ses_glm] = fileparts(glmFilelist.folder); % not needed
    glmFolder = fullfile(glmResultsfolder, ses_glm); % not needed
    betaResults = load(fullfile(glmFilelist.folder, glmFilelist.name));
    betaResults = betaResults.results;
    
    % load in pRF data for that subject (flicker)
%     flickerRetDir = dir(fullfile(bidsDir, 'derivatives', 'prfvista', subjectname, '**/stimfiles.mat'));
%     flickerRetDir = flickerRetDir.folder;
    movingRetDir = dir(fullfile(bidsDir, 'derivatives', 'prfvista_mov', subjectname, '**/stimfiles.mat'));
    movingRetDir = movingRetDir.folder;
    
    for hi=1:numel(hemis)
        hemi = hemis{hi};
%         ret_flicker.(sprintf('%s_pa', hemi)) = MRIread(fullfile(flickerRetDir, sprintf('%s.angle_adj.mgz', hemi)));
%         ret_flicker.(sprintf('%s_ecc', hemi)) = MRIread(fullfile(flickerRetDir, sprintf('%s.eccen.mgz', hemi)));
%         ret_flicker.(sprintf('%s_vexp', hemi)) = MRIread(fullfile(flickerRetDir, sprintf('%s.vexpl.mgz', hemi)));
%         ret_flicker.(sprintf('%s_sigma', hemi)) = MRIread(fullfile(flickerRetDir, sprintf('%s.sigma.mgz', hemi)));
        ret_moving.(sprintf('%s_pa', hemi)) = MRIread(fullfile(movingRetDir, sprintf('%s.angle_adj.mgz', hemi)));
        ret_moving.(sprintf('%s_ecc', hemi)) = MRIread(fullfile(movingRetDir, sprintf('%s.eccen.mgz', hemi)));
        ret_moving.(sprintf('%s_vexp', hemi)) = MRIread(fullfile(movingRetDir, sprintf('%s.vexpl.mgz', hemi)));
        ret_moving.(sprintf('%s_sigma', hemi)) = MRIread(fullfile(movingRetDir, sprintf('%s.sigma.mgz', hemi)));
        
        for cc=1:numel(indiv_contrastnames)
            cond = indiv_contrastnames{cc};
            mot_v_stat.(sprintf('%s_%s', hemi,cond)) = MRIread(fullfile(glmFilelist.folder, sprintf('%s.%s.mgz', hemi, cond)));
        end
    end
    
%     retFlickerData = [ret_flicker.lh_pa.vol, ret_flicker.rh_pa.vol ; ...
%         ret_flicker.lh_ecc.vol, ret_flicker.rh_ecc.vol ; ...
%         ret_flicker.lh_vexp.vol, ret_flicker.rh_vexp.vol; ...
%         ret_flicker.lh_sigma.vol, ret_flicker.rh_sigma.vol];
    
    retMovingData = [ret_moving.lh_pa.vol, ret_moving.rh_pa.vol ; ...
        ret_moving.lh_ecc.vol, ret_moving.rh_ecc.vol ; ...
        ret_moving.lh_vexp.vol, ret_moving.rh_vexp.vol; ...
        ret_moving.lh_sigma.vol, ret_moving.rh_sigma.vol];

    indvCondData = [mot_v_stat.lh_m0_v_s90.vol', mot_v_stat.rh_m0_v_s90.vol' ; ...
        mot_v_stat.lh_m90_v_s0.vol', mot_v_stat.rh_m90_v_s0.vol' ; ...
        mot_v_stat.lh_m180_v_s90.vol', mot_v_stat.rh_m180_v_s90.vol'; ...
        mot_v_stat.lh_m270_v_s0.vol', mot_v_stat.rh_m270_v_s0.vol'; ...
        mot_v_stat.lh_m45_v_s135.vol', mot_v_stat.rh_m45_v_s135.vol'; ...
        mot_v_stat.lh_m135_v_s45.vol', mot_v_stat.rh_m135_v_s45.vol'; ...
        mot_v_stat.lh_m225_v_s135.vol', mot_v_stat.rh_m225_v_s135.vol'; ...
        mot_v_stat.lh_m315_v_s45.vol', mot_v_stat.rh_m315_v_s45.vol'];


    if strcmp(retinotopyMode, 'moving')
        retData = retMovingData;
    elseif strcmp(retinotopyMode, 'flicker')
        retData = retFlickerData;
    end

    % get lh, rh sizes:
    hSize = get_surfsize(subjectname); 
    
    % create matrix of mean values per ROI (V1, V2, V3, MT+):
    
    for ci=1:numel(contrastnames)
    
        contrastname = contrastnames{ci};
    
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

            
    
            % load in results (units are in PSC)
            propername = strrep(contrastname, '_v_','V');
            currBold = betaResults.contrasts.(propername);
                
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
                
%                 if length(finalIdx)>nPreallVoxel
%                     disp('WARNING. NOT ENOUGH MEMORY ALLOCATED TO THE allvoxelsBOLDpa')
%                     disp('THIS WILL RESULT IN TRUCTATION.')
%                 else
%                     disp('ok')
%                 end

                if ci==1 % only count in the first interation to avoid redundancy
                    countBin(pa, ri, si) = length(finalIdx);
                end

            end

        end
    end

end


