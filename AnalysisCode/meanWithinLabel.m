clc; clear all; close all;

% set up
addpath(genpath(pwd));
projectName = 'dg';
bidsDir =  '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
githubDir = '~/Documents/GitHub';
fsDir = '/Applications/freesurfer/7.2.0';
addpath(genpath(fullfile(githubDir, 'retinotopy-nyuad')));
addpath(genpath(fullfile(githubDir, 'wpToolbox')));
addpath(genpath(fullfile(githubDir, 'atlasmgz')));
setup_user(projectName,bidsDir,githubDir,fsDir);
masksFolder = fullfile(bidsDir, 'derivatives', 'masks');
retinotopyMode = 'moving'; % of 'flicker' to use for analysis

hemis = {'lh'; 'rh'};
%rois = {'V1', 'V2v', 'V3v', 'V2d', 'V3d', 'hV4', 'hMTcomplex', 'pMT', 'pMST'};
rois = {'V1', 'V2', 'V3', 'hV4', 'hMTcomplex', 'pMT', 'pMST'};

% replaced below to pairwise (sep): 'cardmVcards', 'oblmVobls'
contrastnames = {'cardMsep', 'oblMsep', 'allmValls', ...
    'allsVblank', 'allmVblank', 'cardmVblank', 'oblmVblank', ...
    'm0_v_s90','m90_v_s0','m180_v_s90','m270_v_s0', ...
    'm45_v_s135','m135_v_s45','m225_v_s135','m315_v_s45', 'cardsVblank', 'oblsVblank', ...
     'm0_v_b', 'm180_v_b', 'm90_v_b', 'm270_v_b', 'm45_v_b', 'm225_v_b', 'm135_v_b', ...
     'm315_v_b', 's0_v_b', 's90_v_b', 's45_v_b', 's135_v_b'};
subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-0397', ...
    'sub-0442', 'sub-wlsubj121', ...
    'sub-wlsubj123', 'sub-wlsubj124', 'sub-wlsubj127', 'sub-0395'}; %'sub-0426', 'sub-0427', 

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
minECC = 0.2; maxECC = 15;
minVAREXP = 0.25;

% initialize a matrix for contrasts, rois, subjects
meanBOLD = nan(length(contrastnames), length(rois), length(subjects));
medianBOLD = nan(length(contrastnames), length(rois), length(subjects));

meanBOLDpa = nan(length(contrastnames), length(polarAngles), length(rois), length(subjects));
medianBOLDpa = nan(length(contrastnames),length(polarAngles), length(rois), length(subjects));
countBin = nan(length(polarAngles), length(rois), length(subjects));

%%

for si=1:numel(subjects)

    subjectname = subjects{si}

    subj_FSfolder = fullfile(bidsDir, 'derivatives', 'freesurfer', subjectname);
    
    glmResultsfolder = fullfile(bidsDir, 'derivatives', 'dgGLM', subjectname);
    glmFilelist = dir(fullfile(glmResultsfolder, '**', 'results.mat'));
    [~,ses_glm] = fileparts(glmFilelist.folder); % not needed
    glmFolder = fullfile(glmResultsfolder, ses_glm); % not needed
    betaResults = load(fullfile(glmFilelist.folder, glmFilelist.name));
    betaResults = betaResults.results;
    
    % load in pRF data for that subject (flicker)
    flickerRetDir = dir(fullfile(bidsDir, 'derivatives', 'prfvista', subjectname, '**/stimfiles.mat'));
    flickerRetDir = flickerRetDir.folder;
    movingRetDir = dir(fullfile(bidsDir, 'derivatives', 'prfvista_mov', subjectname, '**/stimfiles.mat'));
    movingRetDir = movingRetDir.folder;
    
    for hi=1:numel(hemis)
        hemi = hemis{hi};
        ret_flicker.(sprintf('%s_pa', hemi)) = MRIread(fullfile(flickerRetDir, sprintf('%s.angle_adj.mgz', hemi)));
        ret_flicker.(sprintf('%s_ecc', hemi)) = MRIread(fullfile(flickerRetDir, sprintf('%s.eccen.mgz', hemi)));
        ret_flicker.(sprintf('%s_vexp', hemi)) = MRIread(fullfile(flickerRetDir, sprintf('%s.vexpl.mgz', hemi)));
        ret_flicker.(sprintf('%s_sigma', hemi)) = MRIread(fullfile(flickerRetDir, sprintf('%s.sigma.mgz', hemi)));
        ret_moving.(sprintf('%s_pa', hemi)) = MRIread(fullfile(movingRetDir, sprintf('%s.angle_adj.mgz', hemi)));
        ret_moving.(sprintf('%s_ecc', hemi)) = MRIread(fullfile(movingRetDir, sprintf('%s.eccen.mgz', hemi)));
        ret_moving.(sprintf('%s_vexp', hemi)) = MRIread(fullfile(movingRetDir, sprintf('%s.vexpl.mgz', hemi)));
        ret_moving.(sprintf('%s_sigma', hemi)) = MRIread(fullfile(movingRetDir, sprintf('%s.sigma.mgz', hemi)));
        
        for cc=1:numel(indiv_contrastnames)
            cond = indiv_contrastnames{cc};
            mot_v_stat.(sprintf('%s_%s', hemi,cond)) = MRIread(fullfile(glmFilelist.folder, sprintf('%s.%s.mgz', hemi, cond)));
        end
    end
    
    retFlickerData = [ret_flicker.lh_pa.vol, ret_flicker.rh_pa.vol ; ...
        ret_flicker.lh_ecc.vol, ret_flicker.rh_ecc.vol ; ...
        ret_flicker.lh_vexp.vol, ret_flicker.rh_vexp.vol; ...
        ret_flicker.lh_sigma.vol, ret_flicker.rh_sigma.vol];
    
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

                if ci==1 % only count in the first interation to avoid redundancy
                    countBin(pa, ri, si) = length(finalIdx);
                end

            end

        end
    end

end

%% Color values
% removing subject with extreme motion
meanBOLDpa = meanBOLDpa(:,:,:,[1,2,3,4,5,6,7,8,9,10]); %,11,12]);
meanBOLD = meanBOLD(:,:,[1,2,3,4,5,6,7,8,9,10]); %,11,12]);
% to leave out one subject
subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-0397', ...
    'sub-0442', 'sub-wlsubj121', ... 
    'sub-wlsubj123', 'sub-wlsubj124', 'sub-wlsubj127', 'sub-0395'}; % 'sub-0426',

hue_values = linspace(0, 1, numel(subjects)+1);
hue_values = hue_values(1:end-1);  % Exclude the last point to avoid duplication

% Set saturation and value for all colors
saturation = 1;
value = .8; %1; % lower brightness

% Convert HSV to RGB
colors_hsv = [hue_values' repmat(saturation, numel(subjects), 1) repmat(value, numel(subjects), 1)];
colors_rgb = hsv2rgb(colors_hsv);

%%
% contrastnames = {
% 1- cardMsep
% 2- oblMsep
% 3 - allmValls
% 4 - allsVblank
% 5 - allmVblank
% 6 - cardmVblank
% 7 - oblmVblank
% 8 - m0_v_s90
% 9 - m90_v_s0
% 10 - m180_v_s90
% 11 - m270_v_s0
% 12 - m45_v_s135
% 13 - m135_v_s45
% 14 - m225_v_s135
% 15 - m315_v_s45
% 16 - cardsVblank
% 17 - oblsVblank
% 18 - m0_v_b
% 19 - m180_v_b 
% 20 - m90_v_b 
% 21 - m270_v_b 
% 22 - m45_v_b 
% 23 - m225_v_b 
% 24 - m135_v_b
% 25 - m315_v_b 
% 26 - s0_v_b 
% 27 - s90_v_b 
% 28 - s45_v_b 
% 29 - s135_v_b

% 4 v 5: clearly shows MT as motion responsive

figure

condIdx1 = 26; %16; 
condIdx2 = 29; %17;

% Select relevant rows and calculate the ratio
condition1 = squeeze(meanBOLD(condIdx1, :, :));  % First condition
condition2 = squeeze(meanBOLD(condIdx2, :, :));  % Second condition

meanMAT = nan(numel(rois), 2);
semMAT = nan(numel(rois), 2);
diffMAT = nan(numel(rois), 1);
diffSEM = nan(numel(rois), 1);

for region=1:numel(rois)
    meanMAT(region,1) = nanmean(condition1(region,:));
    semMAT(region,1) = nanstd(condition1(region, :)) / sqrt(sum(~isnan(condition1(region, :))));
    meanMAT(region,2) = nanmean(condition2(region,:));
    semMAT(region,2) = nanstd(condition2(region, :)) / sqrt(sum(~isnan(condition2(region, :))));
    diffMAT(region) = nanmean(condition1(region,:) - condition2(region,:));
    diffSEM(region) = nanstd(condition1(region,:) - condition2(region,:)) / sqrt(sum(~isnan(condition1(region,:) - condition2(region,:))));
end

colors = colors_rgb; %rand(length(subjects), 3);

for region = 1:size(meanMAT, 1)
    subplot(2, 5, region);  % Adjust subplot layout as needed

    % Plot points with different colors for conditions
    %scatter(1:2, meanMAT(region, :), 'o', 'DisplayName', 'Mean', 'LineWidth', 1.5);
    %hold on;
    scatter(ones(1, length(condition1(region,:))), condition1(region,:), 55, colors, 'filled','MarkerEdgeColor', 'white','LineWidth', 1.5)
    hold on
    scatter(2*ones(1, length(condition2(region,:))), condition2(region,:), 55, colors, 'filled','MarkerEdgeColor', 'white','LineWidth', 1.5)
    hold on

    for si=1:numel(subjects)
        plot([1,2], [condition1(region,si), condition2(region,si)], 'color', [colors(si, :), 0.25], 'LineWidth',1)
        hold on
    end

    % Error bars with different colors for conditions
    errorbar(1:2, meanMAT(region, :), semMAT(region, :), '.', 'Color', 'k', 'DisplayName', 'SEM', 'LineWidth',2);
    hold on

    errorbar(1.5, diffMAT(region), semMAT(region), '.', 'Color', [.5 .5 .5], 'DisplayName', 'SEM', 'LineWidth',2);

    xlim([0.5,2.5])
    if condIdx2==7 || condIdx1==6
        ylim([-.1, .38])
    else
        ylim([-.1, .25])
    end
    title([rois(region)]);
    xlabel('Contrast');
    ylabel('Mean PSC (Motion - Static)');
    f1 = gcf;
    f1.Position = [-3 137 1508 660];
    ax = gca;
    ax.XTick = [1,2];
    condName1 = strrep(contrastnames{condIdx1}, 'sep', '');
    condName2 = strrep(contrastnames{condIdx2}, 'sep', '');
    ax.XTickLabel = {condName1,condName2}; %{contrastnames{1},'',contrastnames{2},''};
    %legend('show');

    hold off;
end

sgtitle(sprintf('Mean BOLD per contrast: %s and %s', condName1,condName2));

%%

roi=1;
for si=1:numel(subjects)
    subjects(si)
    rois(roi)
    cardinaldir_cardloc = meanBOLDpa(1,[1,3,5,7],roi,si);
    obliquedir_cardloc = meanBOLDpa(2,[1,3,5,7],roi,si);
    cardinaldir_obliqloc = meanBOLDpa(1,[2,4,6,8],roi,si);
    obliquedir_obliqloc = meanBOLDpa(2,[2,4,6,8],roi,si);
    
    disp('Cardinal advantage at cardinal locs:')
    %cardinaldir_cardloc-obliquedir_cardloc
    nanmean(cardinaldir_cardloc-obliquedir_cardloc)
    disp('Cardinal advantage at oblique locs:')
    %cardinaldir_obliqloc-obliquedir_obliqloc
    nanmean(cardinaldir_obliqloc-obliquedir_obliqloc)

%     disp('Cardinal effect at cardinal locs vs oblique locs:')
%     cardinaldir_cardloc-cardinaldir_obliqloc

%     disp('Oblique effect at Oblique locs vs Cardinal locs:')
%     obliquedir_obliqloc-obliquedir_cardloc
end

%% Polar plot

figure

for ri=1:length(rois)
    % Specify the region index
    regionIndex = ri;
    
    % Extract the relevant conditions for the specified region
    conditions1 = meanBOLDpa(8:11, :, regionIndex, :);
    conditions2 = meanBOLDpa(12:15, :, regionIndex, :);
    conditions1 = squeeze(conditions1);
    conditions2 = squeeze(conditions2);
    
    % Average across the conditions within subjects
    avgConditions1 = nanmean(conditions1, 1);
    avgConditions2 = nanmean(conditions2, 1);
    avgConditions1 = squeeze(avgConditions1);
    avgConditions2 = squeeze(avgConditions2);
    
    min_axis1 = min(min(squeeze(min(avgConditions1))));
    min_axis2 = min(min(squeeze(min(avgConditions2))));
    axMIN = customRound(min([min_axis1, min_axis2]));
    max_axis1 = max(max(squeeze(max(avgConditions1))));
    max_axis2 = max(max(squeeze(max(avgConditions2))));
    axMAX = customRound(max([max_axis1, max_axis2]));
    axMIN = axMIN - abs(axMAX-axMIN).*0.05;
    axMAX = axMAX + abs(axMAX-axMIN).*0.05;
    
    % Extract polar angles
    anglevals = [90, 45, 0, 315, 270, 225, 180, 135]; % <-- these were manually converted based on the order of polarAngles above (Noah's convention)
    
    vals_1 = nanmean(avgConditions1,2)';
    sem1 = nanstd(avgConditions1,0,2)' ./ sqrt(sum(~isnan(avgConditions1),2)');
    vals_2 = nanmean(avgConditions2,2)';
    sem2 = nanstd(avgConditions2,0,2)' ./ sqrt(sum(~isnan(avgConditions2),2)');
    
    
    % Plot the data on a polar plot
    subplot(2,4, ri)
    
    for subjectIndex = 1:size(meanBOLDpa, 4)
        polarplot(deg2rad(anglevals),vals_1, 'o-', 'Color', [127/255, 191/255, 123/255], 'MarkerFaceColor', [127/255, 191/255, 123/255], 'LineWidth',1.5)
        hold on
        polarplot(deg2rad(anglevals),vals_2, 'o-', 'Color', [175/255, 141/255, 195/255], 'MarkerFaceColor', [175/255, 141/255, 195/255], 'LineWidth',1.5)
        hold on
        polarplot([deg2rad(anglevals(end)), deg2rad(anglevals(1))],[vals_1(end), vals_1(1)], 'o-', 'Color', [127/255, 191/255, 123/255], 'LineWidth',1.5,  'MarkerFaceColor', [127/255, 191/255, 123/255])
        hold on
        polarplot([deg2rad(anglevals(end)), deg2rad(anglevals(1))],[vals_2(end), vals_2(1)], 'o-', 'Color', [175/255, 141/255, 195/255], 'LineWidth',1.5,  'MarkerFaceColor', [175/255, 141/255, 195/255])
        hold on
        p1 = polarplot([deg2rad(anglevals); deg2rad(anglevals)], [vals_1 - sem1; vals_1 + sem1], '-', 'Color', [127/255, 191/255, 123/255], 'LineWidth',1.5);
        hold on
        p2 = polarplot([deg2rad(anglevals); deg2rad(anglevals)], [vals_2 - sem2; vals_2 + sem2], '-', 'Color', [175/255, 141/255, 195/255], 'LineWidth',1.5);
        hold on
        p3  = polarplot(deg2rad(anglevals), avgConditions1(:, subjectIndex)', 'o', 'Color', [127/255, 191/255, 123/255]);
        hold on;
        p4 = polarplot(deg2rad(anglevals), avgConditions2(:, subjectIndex)', 'o', 'Color', [175/255, 141/255, 195/255]);
        hold on
    %     hold on
        asymm = avgConditions1(:, subjectIndex) - avgConditions2(:, subjectIndex);
        if (sum(asymm<0)) ~= 0
            sprintf('WaRNING: %s points are below 0', num2str(sum(asymm<0)))
    %     else
    %         polarplot(deg2rad(anglevals), asymm', 'b.');
        end
        
    end
    
    rlim([axMIN, axMAX])
    %rticks(linspace(axMIN, axMAX, 4));
    thetaticks(0:45:315);
    title(sprintf('BOLD per polar angle in %s (PSC)',rois{regionIndex}));


    if ri==1
        legend('Cardinal', 'Oblique', 'Location', 'northwest');
    end
end


fig1 = gcf;
fig1.Position = [152 247 1702 941];
hold off;

figureDir = strrep(bidsDir, 'data_bids', 'figures');
filename = fullfile(figureDir,sprintf('PAcardvoblmotionVSstatic%s', rois{regionIndex}));
saveas(gcf, filename, 'pdf');

%% CARTESIAN CARDINAL vs OBLIQUE - Combine Polar Angle
% Plot mean across polar angles 
% keep in mind that this equally weighs each PA, whereas there could be
% differential # of voxels representing the PAs

figure

for ii=1:length(rois)

    % Specify the region index
    regionIndex = ii;

    % Extract the relevant conditions for the specified region
    conditions1 = meanBOLDpa(18:21, :, regionIndex, :);
    conditions2 = meanBOLDpa(22:25, :, regionIndex, :);
    conditions1 = squeeze(conditions1);
    conditions2 = squeeze(conditions2);
    
    % Average across the conditions within subjects
    avgConditions1 = nanmean(conditions1, 1);
    avgConditions2 = nanmean(conditions2, 1);
    avgConditions1 = squeeze(avgConditions1);
    avgConditions2 = squeeze(avgConditions2);
    
    min_axis1 = min(min(squeeze(min(avgConditions1))));
    min_axis2 = min(min(squeeze(min(avgConditions2))));
    axMIN = customRound(min([min_axis1, min_axis2]));
    max_axis1 = max(max(squeeze(max(avgConditions1))));
    max_axis2 = max(max(squeeze(max(avgConditions2))));
    axMAX = customRound(max([max_axis1, max_axis2]));
    axMIN = axMIN - abs(axMAX-axMIN).*0.05;
    axMAX = axMAX + abs(axMAX-axMIN).*0.05;
    
    % Already averaged across the conditions within subjects < -- already did this in
    % the loop above
    
    % Extract polar angles
    anglevals = [90, 45, 0, 315, 270, 225, 180, 135]; % <-- these were manually converted based on the order of polarAngles above (Noah's convention)
    
    vals_1 = nanmean(avgConditions1,1)';
    vals_2 = nanmean(avgConditions2,1)';

    vals_1_overall = nanmean(avgConditions1,'all')';
    vals_2_overall = nanmean(avgConditions2,'all')';
    %sem1 = nanstd(avgConditions1,0,2)' ./ sqrt(sum(~isnan(avgConditions1),2)');
    %sem2 = nanstd(avgConditions2,0,2)' ./ sqrt(sum(~isnan(avgConditions2),2)');
    
    % Plot the data on a polar plot
    subplot(2,5,ii)
    
    for subjectIndex = 1:size(meanBOLDpa, 4)
        scatter(1, vals_1(subjectIndex),  30, 'MarkerFaceColor', colors_rgb(subjectIndex,:), 'MarkerEdgeColor', 'none'); % [127/255, 191/255, 123/255]
        hold on
        scatter(2, vals_2(subjectIndex), 30, 'MarkerFaceColor', colors_rgb(subjectIndex,:), 'MarkerEdgeColor', 'none'); %  [175/255, 141/255, 195/255]
        plot([1 2], [vals_1(subjectIndex) vals_2(subjectIndex)], 'Color', colors_rgb(subjectIndex,:))
        xlim([0 3])
        ylim([-0.15 0.25])
    end
    
    scatter(1, vals_1_overall, 70, 'MarkerFaceColor', [127/255, 191/255, 123/255], 'MarkerEdgeColor', 'black');
    hold on
    scatter(2, vals_2_overall,  70, 'MarkerFaceColor', [175/255, 141/255, 195/255], 'MarkerEdgeColor', 'black');
    plot([1 2], [vals_1_overall vals_2_overall], 'k', 'LineWidth', 3)
    
    title(sprintf('BOLD per polar angle in %s',rois{regionIndex}));
    ylabel('PSC')
    set(gca, 'XTick', []);

    if ii==1
        legend('CartCard', 'CartObl', 'Location', 'northeast');
    end

    hold off;
end

fa = gcf;
fa.Position = [1000 555 1514 782];

%%

% Extract relevant motion direction indices
motionDirectionIndices = 8:15;

motiondirs = [nan,nan,nan,nan,nan,nan,nan,0,90,180,270,45,135,225,315];

% Polar angles in the index order
polarAngles = [90, 45, 0, 315, 270, 225, 180, 135];

% Initialize the new matrix for radial, tangential, and other
newMatrix = zeros(3, size(meanBOLDpa, 2), size(meanBOLDpa, 3), size(meanBOLDpa, 4));

for si = 1:size(meanBOLDpa, 4) % subjects
    for ri = 1:size(meanBOLDpa, 3) % region
        for polarIndex = 1:length(polarAngles)
            % Extract data for the current polar angle
            currentPolarData = meanBOLDpa(motionDirectionIndices, polarIndex, ri, si);
            
            radialvals = [];
            tangvals = [];
            othervals = [];
        
            for mi=8:15
        
                motionValue = motiondirs(mi);

                currentVal = meanBOLDpa(mi, polarIndex, ri, si);
        
                % Compute the distances between motion directions and polar angles
                distances = abs(polarAngles(polarIndex) - motionValue);
                
                % Identify indices for radial, tangential, and other
                if distances == 0 || distances == 180 % radial
                    radialvals = [radialvals, currentVal];
                elseif distances == 90 || distances == 270 % tangential
                    tangvals = [tangvals, currentVal];
                else % other
                    othervals = [othervals, currentVal];
                end
                
            end
        
            % Compute the averages and store in the new matrix
            newMatrix(1, polarIndex, ri, si) = mean(radialvals);
            newMatrix(2, polarIndex, ri, si) = mean(tangvals);
            newMatrix(3, polarIndex, ri, si) = mean(othervals);
        end
    end
end


%% Polar plot (RADIAL / TANGENTIAL)

% Plot the data on a polar plot
figure;

for ri=1:length(rois)

    % Specify the region index
    regionIndex = ri;
    
    % Extract the relevant conditions for the specified region
    conditions1 = newMatrix(1, :, regionIndex, :);
    conditions2 = newMatrix(2, :, regionIndex, :);
    avgConditions1 = squeeze(conditions1);
    avgConditions2 = squeeze(conditions2);
    
    min_axis1 = min(min(squeeze(min(avgConditions1))));
    min_axis2 = min(min(squeeze(min(avgConditions2))));
    axMIN = customRound(min([min_axis1, min_axis2]));
    max_axis1 = max(max(squeeze(max(avgConditions1))));
    max_axis2 = max(max(squeeze(max(avgConditions2))));
    axMAX = customRound(max([max_axis1, max_axis2]));
    axMIN = axMIN - abs(axMAX-axMIN).*0.05;
    axMAX = axMAX + abs(axMAX-axMIN).*0.05;
    
    % Already averaged across the conditions within subjects < -- already did this in
    % the loop above
    
    % Extract polar angles
    anglevals = [90, 45, 0, 315, 270, 225, 180, 135]; % <-- these were manually converted based on the order of polarAngles above (Noah's convention)
    
    vals_1 = nanmean(avgConditions1,2)';
    sem1 = nanstd(avgConditions1,0,2)' ./ sqrt(sum(~isnan(avgConditions1),2)');
    vals_2 = nanmean(avgConditions2,2)';
    sem2 = nanstd(avgConditions2,0,2)' ./ sqrt(sum(~isnan(avgConditions2),2)');
    

    subplot(2,4,ri)
    for subjectIndex = 1:size(meanBOLDpa, 4)
        p1  = polarplot(deg2rad(anglevals), avgConditions1(:, subjectIndex)', 'o', 'Color', [176/255, 224/255, 230/255]);
        hold on;
        p2 = polarplot(deg2rad(anglevals), avgConditions2(:, subjectIndex)', 'o', 'Color', 'r');
        hold on
    
    end
    
    hold on
    polarplot(deg2rad(anglevals),vals_1, 'o-', 'Color', [176/255, 224/255, 230/255], 'LineWidth',1.5,  'MarkerFaceColor', [176/255, 224/255, 230/255])
    hold on
    polarplot(deg2rad(anglevals),vals_2, 'o-', 'Color', 'r', 'LineWidth',1.5,  'MarkerFaceColor', 'r')
    hold on
    polarplot([deg2rad(anglevals(end)), deg2rad(anglevals(1))],[vals_1(end), vals_1(1)], 'o-', 'Color', [176/255, 224/255, 230/255], 'LineWidth',1.5,  'MarkerFaceColor', [176/255, 224/255, 230/255])
    hold on
    polarplot([deg2rad(anglevals(end)), deg2rad(anglevals(1))],[vals_2(end), vals_2(1)], 'o-', 'Color', 'r', 'LineWidth',1.5,  'MarkerFaceColor', 'r')
    hold on
    p3 = polarplot([deg2rad(anglevals); deg2rad(anglevals)], [vals_1 - sem1; vals_1 + sem1], '-', 'Color', [176/255, 224/255, 230/255], 'LineWidth',1.5);
    hold on
    p4 = polarplot([deg2rad(anglevals); deg2rad(anglevals)], [vals_2 - sem2; vals_2 + sem2], '-', 'Color', 'r', 'LineWidth',1.5);
    
    
    rlim([axMIN, axMAX])
    thetaticks(0:45:315);
    
    title(sprintf('BOLD per polar angle in %s (PSC)',rois{regionIndex}));

    if ri==1
        legend('Radial', 'Tangential', 'Location', 'northwest');
    end
end


fig1 = gcf;
fig1.Position = [152 247 1702 941];
hold off;

filename = fullfile(figureDir,sprintf('PAradialvtangmotionVSstatic%s', rois{regionIndex}));
saveas(gcf, filename, 'pdf');


%% RADIAL / TANGENTIAL - Combine Polar Angle
% Plot mean across polar angles 
% keep in mind that this equally weighs each PA, whereas there could be
% differential # of voxels representing the PAs


figure

for ii=1:length(rois)

    % Specify the region index
    regionIndex = ii;

    % Extract the relevant conditions for the specified region
    conditions1 = newMatrix(1, :, regionIndex, :);
    conditions2 = newMatrix(2, :, regionIndex, :);
    avgConditions1 = squeeze(conditions1);
    avgConditions2 = squeeze(conditions2);
    
    min_axis1 = min(min(squeeze(min(avgConditions1))));
    min_axis2 = min(min(squeeze(min(avgConditions2))));
    axMIN = customRound(min([min_axis1, min_axis2]));
    max_axis1 = max(max(squeeze(max(avgConditions1))));
    max_axis2 = max(max(squeeze(max(avgConditions2))));
    axMAX = customRound(max([max_axis1, max_axis2]));
    axMIN = axMIN - abs(axMAX-axMIN).*0.05;
    axMAX = axMAX + abs(axMAX-axMIN).*0.05;
    
    % Already averaged across the conditions within subjects < -- already did this in
    % the loop above
    
    % Extract polar angles
    anglevals = [90, 45, 0, 315, 270, 225, 180, 135]; % <-- these were manually converted based on the order of polarAngles above (Noah's convention)
    
    vals_1 = nanmean(avgConditions1,1)';
    vals_2 = nanmean(avgConditions2,1)';

    vals_1_overall = nanmean(avgConditions1,'all')';
    vals_2_overall = nanmean(avgConditions2,'all')';
    %sem1 = nanstd(avgConditions1,0,2)' ./ sqrt(sum(~isnan(avgConditions1),2)');
    %sem2 = nanstd(avgConditions2,0,2)' ./ sqrt(sum(~isnan(avgConditions2),2)');
    
    % Plot the data on a polar plot
    subplot(2,5,ii)
    
    for subjectIndex = 1:size(meanBOLDpa, 4)
        scatter(1, vals_1(subjectIndex),  30, 'MarkerFaceColor', colors_rgb(subjectIndex,:), 'MarkerEdgeColor', 'none'); %, [176/255, 224/255, 230/255]
        hold on
        scatter(2, vals_2(subjectIndex), 30, 'MarkerFaceColor', colors_rgb(subjectIndex,:), 'MarkerEdgeColor', 'none'); %, 'r'
        plot([1 2], [vals_1(subjectIndex) vals_2(subjectIndex)], 'Color', colors_rgb(subjectIndex,:));   %'k')
        xlim([0 3])
        ylim([-0.15 0.26])
    end
    
    scatter(1, vals_1_overall, 70, 'MarkerFaceColor', [176/255, 224/255, 230/255], 'MarkerEdgeColor', 'black');
    hold on
    scatter(2, vals_2_overall,  70, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'black');
    plot([1 2], [vals_1_overall vals_2_overall], 'k', 'LineWidth', 3)
    
    title(sprintf('BOLD per polar angle in %s',rois{regionIndex}));
    ylabel('PSC')
    set(gca, 'XTick', []);

    if ii==1
        legend('Radial', 'Tangential', 'Location', 'northeast');
    end

    hold off;
end

fa = gcf;
fa.Position = [1000 555 1514 782];

% filename = fullfile(figureDir,sprintf('PAradialvtangmotionVSstatic%s', rois{regionIndex}));
% saveas(gcf, filename, 'pdf');

%% Polar plot (POLAR CARDINAL vs POLAR OBLIQUE)

figure
for ri=1:length(rois)

    % Specify the region index
    regionIndex = ri;
    
    % Average the values in the first two rows along the first dimension
    averagedMatrix = mean(newMatrix(1:2, :, :, :), 1); % radial and tangential
    
    % Combine the averaged values with the third row
    polarcardMatrix = cat(1, averagedMatrix, newMatrix(3, :, :, :));
    
    % Extract the relevant conditions for the specified region
    conditions1 = polarcardMatrix(1, :, regionIndex, :);
    conditions2 = polarcardMatrix(2, :, regionIndex, :);
    avgConditions1 = squeeze(conditions1);
    avgConditions2 = squeeze(conditions2);
    
    min_axis1 = min(min(squeeze(min(avgConditions1))));
    min_axis2 = min(min(squeeze(min(avgConditions2))));
    axMIN = customRound(min([min_axis1, min_axis2]));
    max_axis1 = max(max(squeeze(max(avgConditions1))));
    max_axis2 = max(max(squeeze(max(avgConditions2))));
    axMAX = customRound(max([max_axis1, max_axis2]));
    axMIN = axMIN - abs(axMAX-axMIN).*0.05;
    axMAX = axMAX + abs(axMAX-axMIN).*0.05;
    
    % Already averaged across the conditions within subjects < -- already did this in
    % the loop above
    
    % Extract polar angles
    anglevals = [90, 45, 0, 315, 270, 225, 180, 135]; % <-- these were manually converted based on the order of polarAngles above (Noah's convention)
    
    vals_1 = nanmean(avgConditions1,2)';
    sem1 = nanstd(avgConditions1,0,2)' ./ sqrt(sum(~isnan(avgConditions1),2)');
    vals_2 = nanmean(avgConditions2,2)';
    sem2 = nanstd(avgConditions2,0,2)' ./ sqrt(sum(~isnan(avgConditions2),2)');
    
    % Plot the data on a polar plot
    subplot(2,4,ri)
    
    for subjectIndex = 1:size(meanBOLDpa, 4)
        polarplot(deg2rad(anglevals),vals_1, 'o-', 'Color', [139/255, 69/255, 19/255], 'LineWidth',1.5,  'MarkerFaceColor', [139/255, 69/255, 19/255])
        hold on
        polarplot(deg2rad(anglevals),vals_2, 'o-', 'Color', [0/255, 139/255, 139/255], 'LineWidth',1.5,  'MarkerFaceColor', [0/255, 139/255, 139/255])
        hold on
        polarplot([deg2rad(anglevals(end)), deg2rad(anglevals(1))],[vals_1(end), vals_1(1)], 'o-', 'Color', [139/255, 69/255, 19/255], 'LineWidth',1.5,  'MarkerFaceColor', [139/255, 69/255, 19/255])
        hold on
        polarplot([deg2rad(anglevals(end)), deg2rad(anglevals(1))],[vals_2(end), vals_2(1)], 'o-', 'Color', [0/255, 139/255, 139/255], 'LineWidth',1.5,  'MarkerFaceColor', [0/255, 139/255, 139/255])
        hold on
        p1 = polarplot([deg2rad(anglevals); deg2rad(anglevals)], [vals_1 - sem1; vals_1 + sem1], '-', 'Color', [139/255, 69/255, 19/255], 'LineWidth',1.5);
        hold on
        p2 = polarplot([deg2rad(anglevals); deg2rad(anglevals)], [vals_2 - sem2; vals_2 + sem2], '-', 'Color', [0/255, 139/255, 139/255], 'LineWidth',1.5);
        hold on
        p3  = polarplot(deg2rad(anglevals), avgConditions1(:, subjectIndex)', 'o', 'Color', [139/255, 69/255, 19/255]);
        hold on;
        p4 = polarplot(deg2rad(anglevals), avgConditions2(:, subjectIndex)', 'o', 'Color', [0/255, 139/255, 139/255]);
        hold on
    
    end
    
    rlim([axMIN, axMAX])
    thetaticks(0:45:315);
    
    title(sprintf('BOLD per polar angle in %s (PSC)',rois{regionIndex}));


    if ri==1
        legend('Polar cardinal', 'Polar oblique', 'Location', 'northwest');
    end
end


fig1 = gcf;
fig1.Position = [152 247 1702 941];
hold off;

filename = fullfile(figureDir,sprintf('PApolarcardvpolaroblmotionVSstatic%s', rois{regionIndex}));
saveas(gcf, filename, 'pdf');


%% POLAR CARDINAL vs POLAR OBLIQUE - Combine Polar Angle
% Plot mean across polar angles 
% keep in mind that this equally weighs each PA, whereas there could be
% differential # of voxels representing the PAs


figure

for ii=1:length(rois)

    % Specify the region index
    regionIndex = ii;

    % Average the values in the first two rows along the first dimension
    averagedMatrix = mean(newMatrix(1:2, :, :, :), 1); % radial and tangential
    
    % Combine the averaged values with the third row
    polarcardMatrix = cat(1, averagedMatrix, newMatrix(3, :, :, :));
    
    % Extract the relevant conditions for the specified region
    conditions1 = polarcardMatrix(1, :, regionIndex, :);
    conditions2 = polarcardMatrix(2, :, regionIndex, :);
    avgConditions1 = squeeze(conditions1);
    avgConditions2 = squeeze(conditions2);
    
    min_axis1 = min(min(squeeze(min(avgConditions1))));
    min_axis2 = min(min(squeeze(min(avgConditions2))));
    axMIN = customRound(min([min_axis1, min_axis2]));
    max_axis1 = max(max(squeeze(max(avgConditions1))));
    max_axis2 = max(max(squeeze(max(avgConditions2))));
    axMAX = customRound(max([max_axis1, max_axis2]));
    axMIN = axMIN - abs(axMAX-axMIN).*0.05;
    axMAX = axMAX + abs(axMAX-axMIN).*0.05;
    
    % Already averaged across the conditions within subjects < -- already did this in
    % the loop above
    
    % Extract polar angles
    anglevals = [90, 45, 0, 315, 270, 225, 180, 135]; % <-- these were manually converted based on the order of polarAngles above (Noah's convention)
    
    vals_1 = nanmean(avgConditions1,1)';
    vals_2 = nanmean(avgConditions2,1)';

    vals_1_overall = nanmean(avgConditions1,'all')';
    vals_2_overall = nanmean(avgConditions2,'all')';
    %sem1 = nanstd(avgConditions1,0,2)' ./ sqrt(sum(~isnan(avgConditions1),2)');
    %sem2 = nanstd(avgConditions2,0,2)' ./ sqrt(sum(~isnan(avgConditions2),2)');
    
    % Plot the data on a polar plot
    subplot(2,5,ii)
    
    for subjectIndex = 1:size(meanBOLDpa, 4)
        scatter(1, vals_1(subjectIndex),  30, 'MarkerFaceColor', colors_rgb(subjectIndex,:), 'MarkerEdgeColor', 'none'); %, [139/255, 69/255, 19/255]
        hold on
        scatter(2, vals_2(subjectIndex), 30, 'MarkerFaceColor', colors_rgb(subjectIndex,:), 'MarkerEdgeColor', 'none'); %, [0/255, 139/255, 139/255]
        plot([1 2], [vals_1(subjectIndex) vals_2(subjectIndex)], 'Color', colors_rgb(subjectIndex,:));  %; %'k')
        xlim([0 3])
        ylim([-0.15 0.25])
    end
    
    scatter(1, vals_1_overall, 70, 'MarkerFaceColor', [139/255, 69/255, 19/255], 'MarkerEdgeColor', 'black');
    hold on
    scatter(2, vals_2_overall,  70, 'MarkerFaceColor', [0/255, 139/255, 139/255], 'MarkerEdgeColor', 'black');
    plot([1 2], [vals_1_overall vals_2_overall], 'k', 'LineWidth', 3)
    
    title(sprintf('BOLD per polar angle in %s',rois{regionIndex}));
    ylabel('PSC')
    set(gca, 'XTick', []);

    if ii==1
        legend('PolarCard', 'PolarObl', 'Location', 'northeast');
    end

    hold off;
end

fa = gcf;
fa.Position = [1000 555 1514 782];

%%

function roundedValue = customRound(value)
    % Check if the value is positive or negative
    if value >= 0
        % Round up for positive values
        roundedValue = ceil(value * 1000) / 1000;
    else
        % Round down for negative values
        roundedValue = floor(value * 1000) / 1000;
    end
end
