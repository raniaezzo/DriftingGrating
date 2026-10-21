clc; clear all; close all

% set overwrite=1 to refit the LME and rerun the bootstrap (slow); set
% overwrite=0 to skip both and just load the previously saved modeldata/
% estimates/boot files for plotting.
overwrite = 0;

% set up
addpath(genpath(pwd));
projectName = 'da';
bidsDir =  '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
%bidsDir =  '/Volumes/server/Projects/Project_dg/data_bids/';
githubDir = '~/Documents/GitHub';
hRF_setting = 'glmsingle';
fullfile(githubDir, 'DriftingGrating', 'AnalysisCode')
glmResultsfolder = fullfile(bidsDir, 'derivatives', strcat(projectName, 'GLM'), strcat('hRF_', hRF_setting));

% can be 'motion_minus_orientation' ; 'motion_minus_baseline' ; 'orientation_minus_baseline'
comparisonName = 'orientation_minus_baseline';

projectSettings = loadConfig(githubDir);

rois = projectSettings.rois;
roi_idx = projectSettings.roi_idx; % rois{ri} is the ri-th ROI in THIS list, but
    % roi_idx{ri} is its actual column in meanBOLDpa/medianBOLDpa's ROI
    % dimension -- these only coincide if rois lists every ROI in the data
    % matrix in the matrix's own order, which is not guaranteed in general.
axes_limits = projectSettings.axes_limits;
pairaxes_limits = projectSettings.pairaxes_limits;
pairaxes_PAew_limits = projectSettings.pairaxes_PAew_limits;
colors_data = projectSettings.colors_data;
contrasts_dict = projectSettings.contrasts_dict;

%rois = rois(1:7); % remove once I process v3a / v3b
metric = 'bold';

figureDir = [strrep(bidsDir, 'data_bids', 'figures'), projectName];

%% Load and remove subject with extreme motion

load(fullfile(glmResultsfolder, 'meanBOLDpa')) % contrasts x polarAngles x ROIs x subjects
load(fullfile(glmResultsfolder, 'meanBOLD')) % contrasts x ROIs x subjects

load(fullfile(glmResultsfolder, 'medianBOLDpa')) % contrasts x polarAngles x ROIs x subjects
load(fullfile(glmResultsfolder, 'medianBOLD')) % contrasts x ROIs x subjects


figureDir = [strrep(bidsDir, 'data_bids', 'figures'), projectName];

if ~isfolder(figureDir)
    mkdir(figureDir)
end

if strcmp(projectName, 'dg')
    meanBOLDpa = meanBOLDpa(:,:,:,[1,2,3,4,5,6,7,8]); % to ensure the same subjects as DA
    meanBOLD = meanBOLD(:,:,[1,2,3,4,5,6,7,8]);
    medianBOLDpa = medianBOLDpa(:,:,:,[1,2,3,4,5,6,7,8]); % to ensure the same subjects as DA
    medianBOLD = medianBOLD(:,:,[1,2,3,4,5,6,7,8]);
    %to include only the 8 repeat subjects
    subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
       'sub-0395', 'sub-0426', 'sub-0250'};

    % medianBOLDpa = medianBOLDpa(:,:,:,[1,2,3,4,5,6,7,8,9,10,11,12,13]);
    % medianBOLD = medianBOLD(:,:,[1,2,3,4,5,6,7,8,9,10,11,12,13]);
    % subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', ...
    %     'sub-wlsubj124', 'sub-0395', 'sub-0426', 'sub-0250', ...
    %     'sub-0442', 'sub-wlsubj121', 'sub-wlsubj127',  'sub-0397', ...
    %     'sub-0427'};

elseif strcmp(projectName, 'da')
    % load subjects
    subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
        'sub-0395', 'sub-0426', 'sub-0250'};
%     subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
%         'sub-0395', 'sub-0426'};
end

% per-observer mean pRF gain (prfvista_mov/prfvista average), used to
% gain-weight each observer's contribution to the LME below -- same
% correction as in plot1_experimentalCond.m/plot2_experimentalCond.m, see
% retrieveObserverGainWeights.m. subjectScale(i) = groupGain / gain_i, so
% multiplying a subject's data by subjectScale(i) both re-weights them
% relative to the group (low gain -> up-weighted) AND restores the overall
% scale to be comparable to the original (unweighted) BOLD units.
%
% groupGain uses the geometric, not arithmetic, mean of gainWeights:
% subjectScale is itself a multiplicative scale factor (a ratio), so
% geometric mean is the choice under which the *geometric* mean of
% subjectScale across observers is exactly 1 (magnitude-neutral on
% average, in the multiplicative sense that matches what subjectScale
% actually is) -- the arithmetic mean of gainWeights does not have this
% property and systematically inflates the restored units. Implemented
% via exp(mean(log(.))) rather than geomean() to avoid a dependency on
% the Statistics and Machine Learning Toolbox.
gainSummaryFile = fullfile(bidsDir, 'derivatives', 'summaryTables', 'gainSummary.mat');
gainWeights = retrieveObserverGainWeights(subjects, gainSummaryFile);
if any(gainWeights <= 0)
    error('gainWeights must be strictly positive to take log() for the geometric mean (found %d non-positive value(s))', ...
        sum(gainWeights <= 0));
end
groupGain = exp(mean(log(gainWeights)));
subjectScale = groupGain ./ gainWeights;

radialvstang = 0;
[proConditions, conConditions, allConditions] = retrieveProConIdx(projectName, comparisonName, radialvstang);

filtered_meanBOLDpa = meanBOLDpa(allConditions, :, :, :);

% for now, only use DG names -- they apply to both DG and DA b/c directions
% defined in absolute reference frame
%contrastnames = {contrasts_dict.contrasts.(strcat(projectName, '_contrast_name'))};
contrastnames = {contrasts_dict.contrasts.('dg_contrast_name')};

% Extract integers before '_v_'
mdirvals = cellfun(@(x) str2double(regexp(x, '\d+(?=_v_)', 'match', 'once')), {contrastnames{allConditions}});
mdirvals = mdirvals';

[nMotDirs, nPAs, nROIs, nSubjs] = size(filtered_meanBOLDpa);

%%

% polar angle
anglevals = [90; 45; 0; 315; 270; 225; 180; 135];

maincardinalmDir = [0; 90; 180; 270]; % this is up/down/left/right for DG
                                      % and in/out/cc/c for DA
primaryMeridians = [90; 0; 270; 180];


if overwrite
rng(0)
for roi=1:length(rois)  % just 1 ROI at a time (makes interprettability easier)

    saveDir = fullfile(glmResultsfolder,'LME_results', comparisonName, rois{roi});
    
    if ~isfolder(saveDir)
        mkdir(saveDir)
    end
    
    % Loop through each element in the new matrix
    index = 1;

    % Preallocate the new matrix (1 long column)
    reshaped_mat = zeros(nMotDirs*nPAs*nSubjs, 1);

    dirCol = repmat(mdirvals, nPAs*nSubjs, 1); 
    paRep = repelem(anglevals, nPAs);
    paCol = repmat(paRep, nSubjs, 1);
    subCol = repelem(1:nSubjs, 64)';

    for subject = 1:nSubjs
        for pa = 1:nPAs
            for md = 1:nMotDirs
               % fill the new reshaped matrix with the foiled values
               % (condition, polar angle, roi, subject) -- index into the
               % data matrix via roi_idx{roi}, NOT the raw loop counter
               % roi, since rois{roi} is not guaranteed to be the roi-th
               % column of the underlying matrix
               reshaped_mat(index, 1) = filtered_meanBOLDpa(md, pa, roi_idx{roi}, subject);
               index = index + 1;
            end
        end
    end

    %%
    % now determine which are main cardinal, derived cardinal, and radial:

    if strcmp(comparisonName, 'orientation_minus_baseline')
        % just for orientation, because horizontal is both 0 and 180 deg
        % this repeats each value consecutively
        reshaped_mat = repelem(reshaped_mat, 2);
        dirCol = repelem(dirCol, 2);
    end

    finalMat = [reshaped_mat, dirCol, paCol, subCol];

    if strcmp(comparisonName, 'orientation_minus_baseline')
        finalMat = finalMat(1:2:end, :);
    end

    % initialize value for each asymmetry
    mainCardinal = zeros(length(finalMat),1);
    derivedCardinal = zeros(length(finalMat),1);
    mainSubset = zeros(length(finalMat),1);
    derivedSubset = zeros(length(finalMat),1);

    % this will select which directions are main cardinal (dg: up, down,
    % left, right = 1 vs NOT = -1 ; and da: in, out, cc, c = 1 vs NOT = -1)
    mainCardinal_idx = ismember(finalMat(:,2), maincardinalmDir);
    mainCardinal(mainCardinal_idx) = 1;
    mainCardinal(~mainCardinal_idx) = -1;

    % this will select which directions are derived cardinal (dg: in, out, 
    % cc, c = 1 vs NOT = -1 ; and da: up, down, left, right = 1 vs NOT = -1 )
    derivedCardinal_idx = (ismember(finalMat(:,2), maincardinalmDir) & ismember(finalMat(:,3), primaryMeridians)) ...
        | ((~ismember(finalMat(:,2), maincardinalmDir) & (~ismember(finalMat(:,3), primaryMeridians))));
    derivedCardinal(derivedCardinal_idx) = 1;
    derivedCardinal(~derivedCardinal_idx) = -1;

    if strcmp(projectName, 'dg')
        derived_subset_pro_idx = (abs(finalMat(:,2)-finalMat(:,3)) == 0 | abs(finalMat(:,2)-finalMat(:,3)) == 180); % in/out, radial
        derived_subset_con_idx = abs(finalMat(:,2)-finalMat(:,3)) == 90 | abs(finalMat(:,2)-finalMat(:,3)) == 270; % c/cc, tangential

        % added
        main_subset_pro_idx = ismember(finalMat(:,2), [0, 180]); % for dg, this is right/left, horizontal
        main_subset_con_idx = ismember(finalMat(:,2), [90, 270]); % for dg, this is up/down, vertical
    elseif strcmp(projectName, 'da')
        main_subset_pro_idx = ismember(finalMat(:,2), [90, 270]); % for da, this is in/out, radial
        main_subset_con_idx = ismember(finalMat(:,2), [0, 180]); % for da, this is c/cc, tangential

        % added
%         derived_subset_pro_idx = abs(finalMat(:,2)-finalMat(:,3)) == 90 | abs(finalMat(:,2)-finalMat(:,3)) == 270; % horizontal
%         derived_subset_con_idx = (abs(finalMat(:,2)-finalMat(:,3)) == 0 | abs(finalMat(:,2)-finalMat(:,3)) == 180); % vertical
        
        derived_subset_pro_idx = (ismember(finalMat(:,2), [0, 90]) & (abs(finalMat(:,2)-finalMat(:,3)) == 90 | abs(finalMat(:,2)-finalMat(:,3)) == 270)) | ...
            (ismember(finalMat(:,2), [45, 135]) & (abs(finalMat(:,2)-finalMat(:,3)) == 0 | abs(finalMat(:,2)-finalMat(:,3)) == 180)); % horizontal
        
        derived_subset_con_idx = (ismember(finalMat(:,2), [0, 90]) & (abs(finalMat(:,2)-finalMat(:,3)) == 0 | abs(finalMat(:,2)-finalMat(:,3)) == 180)) | ...
            (ismember(finalMat(:,2), [45, 135]) & (abs(finalMat(:,2)-finalMat(:,3)) == 90 | abs(finalMat(:,2)-finalMat(:,3)) == 270)); % vertical

    end

    mainSubset(main_subset_pro_idx) = 1;
    mainSubset(main_subset_con_idx) = -1;

    derivedSubset(derived_subset_pro_idx) = 1;
    derivedSubset(derived_subset_con_idx) = -1;

    finalMat = [finalMat, mainCardinal, derivedCardinal, mainSubset, derivedSubset];

    %%

    variable_names = {'bold', 'motiondir', 'polarangle', 'sub', 'mainCardinal', 'derivedCardinal', 'mainSubset', 'derivedSubset'};
    modeldata = array2table(finalMat, 'VariableNames', variable_names);

    % Gain-weight each observer's BOLD values (see subjectScale above):
    % divide by their own mean pRF gain, then multiply the across-observer
    % average gain back in. Applied per row via that row's subject index
    % (modeldata.sub), so every row for a given subject -- across all polar
    % angles/motion directions -- gets the SAME scale factor. This is
    % computed BEFORE fitting and BEFORE saving modeldata to disk, so the
    % LME fit/anova/estimates below, and the bootstrap section further down
    % (which reloads this saved modeldata), all inherit the gain-weighted,
    % rescaled values automatically.
    modeldata.bold = modeldata.bold .* subjectScale(modeldata.sub)';

    modeldata.subject = categorical(modeldata.sub);

    lme = fitlme(modeldata, 'bold ~ mainCardinal + derivedCardinal + mainSubset + derivedSubset + (1|sub)');

    anova(lme, 'dfmethod', 'satterthwaite')
    
    global_idx = find(contains(cellstr(lme.Coefficients(:,1)), 'Intercept'));
    maincard_idx = find(contains(cellstr(lme.Coefficients(:,1)), 'mainCardinal'));
    derivedcard_idx = find(contains(cellstr(lme.Coefficients(:,1)), 'derivedCardinal'));
    mainsub_idx = find(contains(cellstr(lme.Coefficients(:,1)), 'mainSubset'));
    derivedsub_idx = find(contains(cellstr(lme.Coefficients(:,1)), 'derivedSubset'));

    global_est = double(lme.Coefficients(global_idx,2));
    maincard_est = double(lme.Coefficients(maincard_idx,2));
    derivedcard_est = double(lme.Coefficients(derivedcard_idx,2));
    mainsub_est = double(lme.Coefficients(mainsub_idx,2));
    derivedsub_est = double(lme.Coefficients(derivedsub_idx,2));

    estimates = [global_est, maincard_est, derivedcard_est, mainsub_est, derivedsub_est];

    save(fullfile(saveDir,strcat('LME_',metric)), 'estimates');

    disp(lme)
    save(strcat(saveDir, '/modeldata'), 'modeldata');
end
end % if overwrite

%% bootstrap data

if overwrite
subjects = {'ALL'};

rng('default'); rng(1);

subjID = 1:nSubjs; %10;
bootN = 1000;

for roi=1:length(rois)

    roi

    saveDir = fullfile(glmResultsfolder,'LME_results', comparisonName, rois{roi});
    
    saveboot = {};
    coeffs = nan(5,bootN);
    load(strcat(saveDir, '/modeldata'), 'modeldata');

    for bi=1:bootN

        % changed this to nSubjects == #random samples (second unit)-- used
        % to be 10
        y = datasample(subjID,nSubjs); %10); % take a random sample of subjects

        %modeldata = readtable(strcat('/Users/rania/Desktop/RadialBias_pilot1/Data/ALLSUBSwFULLDATA/Analysis_1_PF/datatable.csv'));
        
        randomsample = [];
        for i=1:length(y)
            selectsub = y(i);
            temp = modeldata(modeldata.sub == selectsub,:);
            randomsample = [randomsample ; temp];
        end

        lme = fitlme(randomsample,'bold ~ mainCardinal + derivedCardinal + mainSubset + derivedSubset + (1|sub)');

        saveboot{bi} = lme;

        global_idx = find(contains(cellstr(lme.Coefficients(:,1)), 'Intercept'));
        maincard_idx = find(contains(cellstr(lme.Coefficients(:,1)), 'mainCardinal'));
        derivedcard_idx = find(contains(cellstr(lme.Coefficients(:,1)), 'derivedCardinal'));
        mainsub_idx = find(contains(cellstr(lme.Coefficients(:,1)), 'mainSubset'));
        derivedsub_idx = find(contains(cellstr(lme.Coefficients(:,1)), 'derivedSubset'));
    
        global_est = double(lme.Coefficients(global_idx,2));
        maincard_est = double(lme.Coefficients(maincard_idx,2));
        derivedcard_est = double(lme.Coefficients(derivedcard_idx,2));
        mainsub_est = double(lme.Coefficients(mainsub_idx,2));
        derivedsub_est = double(lme.Coefficients(derivedsub_idx,2));

        estimates = [global_est, maincard_est, derivedcard_est, mainsub_est, derivedsub_est];
        coeffs(:,bi) = estimates;
        clear lme
        disp(bi)



    end

    save(strcat(saveDir, '/boot'), 'saveboot','coeffs');
end
end % if overwrite



%% plot

% plot across rois per figure
asymmetryNames = {'mainCardinalVsMainOblique', 'derivedCardinalVsDerivedOblique', 'mainSubset', 'derivedSubset'};

% initialize arrays to store data for 1 ROI : for a later plot
nA = numel(asymmetryNames);
nROIs_plot = numel(rois);
box_pro = cell(nA,1);
box_con = cell(nA,1);
mean_pro = nan(nA,1);
mean_con = nan(nA,1);
errlow_pro = nan(nA,1);
errhigh_pro = nan(nA,1);
errlow_con = nan(nA,1);
errhigh_con = nan(nA,1);
color_pro = nan(nA,3);
color_con = nan(nA,3);
asymLabel = strings(nA,1);

% same as above, but for EVERY roi (not just roi==1), so a per-ROI version
% of the master figure can be made below -- indexed (ai, roi)
mean_pro_allroi = nan(nA, nROIs_plot);
mean_con_allroi = nan(nA, nROIs_plot);
errlow_pro_allroi = nan(nA, nROIs_plot);
errhigh_pro_allroi = nan(nA, nROIs_plot);
errlow_con_allroi = nan(nA, nROIs_plot);
errhigh_con_allroi = nan(nA, nROIs_plot);

%all_labels = {{'Main Cardinal', 'Main Oblique'}, {'Derived Cardinal', 'Derived Oblique'}, {'Radial', 'Tangential'}};

ci_level = 84; %68; %68;
meanRelative=1;
%colors = [[127 191 123]/255; [166 97 26]/255; [146 197 222]/255];
%colors2 = [[175 141 195]/255; [64 176 166]/255; [202 0 32]/255];

for ai=1:numel(asymmetryNames)

    x = ai; % this is condition 1, 2, or 3

    asymmetryName = asymmetryNames{ai};

    % specify the derived subset for dg and da (or later rename items in COLORS.json)
    if strcmp(projectName, 'dg')
        if strcmp(asymmetryName, 'mainSubset')
            asymmetryName = 'verticalVsHorizontal';
        elseif strcmp(asymmetryName, 'derivedSubset')
            asymmetryName = 'radialVsTangential';
        end
    elseif strcmp(projectName, 'da')
        if strcmp(asymmetryName, 'mainSubset')
            asymmetryName = 'radialVsTangential';
        elseif strcmp(asymmetryName, 'derivedSubset')
            asymmetryName = 'verticalVsHorizontal';
        end
    end

    % human-readable label for the master figure's x-tick (kept separate
    % from asymmetryName, which is also used as a COLORS.json lookup key
    % that only has entries for the literal 'mainCardinalVsMainOblique'/
    % 'derivedCardinalVsDerivedOblique' names, not these renamed ones)
    displayLabel = asymmetryName;
    if strcmp(projectName, 'dg')
        if strcmp(asymmetryName, 'mainCardinalVsMainOblique')
            displayLabel = 'cardinalVsOblique';
        elseif strcmp(asymmetryName, 'derivedCardinalVsDerivedOblique')
            displayLabel = 'polarCardinalVsPolarOblique';
        end
    elseif strcmp(projectName, 'da')
        if strcmp(asymmetryName, 'mainCardinalVsMainOblique')
            displayLabel = 'polarCardinalVsPolarOblique';
        elseif strcmp(asymmetryName, 'derivedCardinalVsDerivedOblique')
            displayLabel = 'cardinalVsOblique';
        end
    end

    colors = colors_data.conditions.(projectName).(asymmetryName).color_pro';
    colors2 = colors_data.conditions.(projectName).(asymmetryName).color_con';

    %labelnames = all_labels(x);
    labelnames = lower(strsplit(asymmetryNames{ai}, 'Vs'));

    % save for later plot
    color_pro(ai,:) = colors;
    color_con(ai,:) = colors2;
    asymLabel(ai) = string(displayLabel);

    figure
    xlim([0 8]);
    ylim([0 8]);
    for roi=1:length(rois)
        
        % load('bootLME_sensitivity.mat')
        % load('LME_sensitivity.mat')
    
        saveDir = fullfile(glmResultsfolder,'LME_results', comparisonName, rois{roi});
    
        load(strcat(saveDir, '/modeldata'), 'modeldata');
        load(fullfile(saveDir,strcat('LME_',metric)), 'estimates');
        load(strcat(saveDir, '/boot'), 'saveboot','coeffs');
    
        Gintercept = estimates(1);
        main_cardinal_est = Gintercept + estimates(2);
        derived_cardinal_est = Gintercept + estimates(3);
        main_subset_est = Gintercept + estimates(4);
        derived_subset_est = Gintercept + estimates(5);
    
        CIFcn = @(x,p)prctile(x, [100-p, p]); 
        
        p = ci_level;
        % coefficient order is: intercept, betas
        CI_maincardinality = CIFcn(coeffs(2,:)+Gintercept,p);
        CI_derivedcardinality = CIFcn(coeffs(3,:)+Gintercept,p);
        CI_mainsubset = CIFcn(coeffs(4,:)+Gintercept,p);
        CI_derivedsubset = CIFcn(coeffs(5,:)+Gintercept,p);

        % this is for the CI of the DIFFERENCE (stats, not the plot)
        % requires multiplying the beta by 2 (.*2) to compute the difference
        % (they have identical magnitude)
        % do not need the mean
        % coefficient order is: intercept, betas
        CI_maincardinality_stats = CIFcn(coeffs(2,:).*2,97.5);
        CI_derivedcardinality_stats = CIFcn(coeffs(3,:).*2,97.5);
        CI_mainsubset_stats = CIFcn(coeffs(4,:).*2,97.5);
        CI_derivedsubset_stats = CIFcn(coeffs(5,:).*2,97.5);

        CI_maincardinality_stats68 = CIFcn(coeffs(2,:).*2,84);
        CI_derivedcardinality_stats68 = CIFcn(coeffs(3,:).*2,84);
        CI_mainsubset_stats68 = CIFcn(coeffs(4,:).*2,84);
        CI_derivedsubset_stats68 = CIFcn(coeffs(5,:).*2,84);

        % labeled term-by-term printout, so each estimate sits next to its
        % own CI unambiguously. Guarded by ai==1 so this prints once per
        % roi instead of once per asymmetryName (these values only depend
        % on roi -- modeldata/estimates/coeffs are reloaded identically
        % regardless of asymmetryName, so printing on every ai iteration
        % just reprinted the same numbers numel(asymmetryNames) times).
        if ai == 1
            diffNames = {'mainCardinal','derivedCardinal','mainSubset','derivedSubset'};
            beta_estimate_difference = [estimates(2) estimates(3) estimates(4) estimates(5)] .* 2;
            diffCIs95 = [CI_maincardinality_stats; CI_derivedcardinality_stats; ...
                CI_mainsubset_stats; CI_derivedsubset_stats];
            diffCIs68 = [CI_maincardinality_stats68; CI_derivedcardinality_stats68; ...
                CI_mainsubset_stats68; CI_derivedsubset_stats68];

            fprintf('\nstats (roi=%s): CI of the pro-vs-con difference (beta*2)\n', rois{roi});
            fprintf('%-16s %12s %24s %24s\n', 'term', 'estimate*2', '95% CI', '68% CI');
            for ni = 1:numel(diffNames)
                fprintf('%-16s %12.5f   [%9.5f %9.5f]   [%9.5f %9.5f]\n', diffNames{ni}, ...
                    beta_estimate_difference(ni), diffCIs95(ni,1), diffCIs95(ni,2), ...
                    diffCIs68(ni,1), diffCIs68(ni,2));
            end
        end
    
        y = [estimates(2) estimates(3) estimates(4) estimates(5)]; %[0.11446, 0.033442, 0.015738];
    
        errlow = [CI_maincardinality(1) CI_derivedcardinality(1) CI_mainsubset(1) CI_derivedsubset(1)];
        errhigh = [CI_maincardinality(2) CI_derivedcardinality(2) CI_mainsubset(2) CI_derivedsubset(2)];
    
        y1 = Gintercept + y;
        y2 = Gintercept - y;

        if ai == 1
            disp('estimates:')
            y1-y2
        end
    
        errlow1 = y1-errlow;
        errhigh1 = errhigh -y1;
    
        errlow2 = errhigh1; %errlow1;
        errhigh2 = errlow1; %errhigh1;
    
        %ax = axes();
        %hold(ax);
    
        if meanRelative
            baselineSub = Gintercept;
        else
            baselineSub = 0;
        end
        %for i=1:length(x)
            boxchart(roi*1, y1(:,x)-baselineSub, 'BoxFaceColor', colors, 'LineWidth', 4, 'BoxWidth', .75)
            hold on
            %plot(roi*1, y1(:,x)-baselineSub, 'Color', colors, 'Marker', '.', 'MarkerSize', 10, 'LineStyle','none')
            %hold on
            errorbar(roi*1,y1(x)-baselineSub,errlow1(x), errhigh1(x), 'LineStyle','none', 'LineWidth', 2, 'Color', colors);
            hold on
            boxchart(roi*1, y2(:,x)-baselineSub, 'BoxFaceColor', colors2, 'LineWidth', 4, 'BoxWidth', .75)
            hold on
            errorbar(roi*1,y2(x)-baselineSub,errlow2(x), errhigh2(x), 'LineStyle','none', 'LineWidth', 2, 'Color', colors2);
            hold on
        %end
        
        hold on

        if roi==1
            % save for later
            k = ai;  % 2..5
        
            box_pro{ai} = y1(:,x) - baselineSub;
            box_con{ai} = y2(:,x) - baselineSub;
        
            mean_pro(ai) = y1(x) - baselineSub;
            mean_con(ai) = y2(x) - baselineSub;
        
            errlow_pro(ai)  = errlow1(x);
            errhigh_pro(ai) = errhigh1(x);
        
            errlow_con(ai)  = errlow2(x);
            errhigh_con(ai) = errhigh2(x);
        end

        % same values as above, stashed for every roi (used by the
        % per-ROI figures below)
        mean_pro_allroi(ai,roi)    = y1(x) - baselineSub;
        mean_con_allroi(ai,roi)    = y2(x) - baselineSub;
        errlow_pro_allroi(ai,roi)  = errlow1(x);
        errhigh_pro_allroi(ai,roi) = errhigh1(x);
        errlow_con_allroi(ai,roi)  = errlow2(x);
        errhigh_con_allroi(ai,roi) = errhigh2(x);

    end
    
    hold on
    % Create dummy plot objects for legend
    h1 = plot(nan, nan, 'Color', colors, 'LineWidth', 3); % Dummy red line
    hold on;
    h2 = plot(nan, nan, 'Color', colors2, 'LineWidth', 3); % Dummy blue line
    
    if meanRelative
        yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2)
        ylim([-0.05 0.05])
    else
        yline(Gintercept, '--', 'Color', [0 0 0], 'LineWidth', 2)
        ylim([-0.05+Gintercept 0.05+Gintercept])
    end
    xlim([0 10])
    ylim([-0.03 0.03])
    set(gca,'XTick',[])
    box off
    set(gca,'linewidth',2, 'YColor', [0 0 0]);
    set(gca,'linewidth',2, 'XColor', [0 0 0]);
    %title(condNames{x})
    set(gca, 'FontName', 'Arial', 'FontSize', 20);
    ax1 = gca;
    %ax1.YTick = [-0.03, -0.015, 0, 0.015, 0.03];
    
    % yticks = get(gca, 'ytick'); % Get current y-axis tick positions
    % text(0.5, mean(yticks), '\Delta', 'VerticalAlignment', 'baseline', 'FontSize', 20); % Add triangle symbol at (0.5, y) position
    
    % Get current y-axis label position
    ylabelHandle = ylabel('temp', 'FontSize', 20);
    ylabelPosition = get(ylabelHandle, 'Position');
    % Remove y-axis label
    delete(ylabelHandle);
    
    % Add triangle symbol in place of y-axis label
    y_pos = ylim;
    %text(ylabelPosition(1), mean(y_pos), '\Delta BOLD signal (%)', 'FontSize', 20, 'Rotation', 90, 'HorizontalAlignment', 'center');
    
    ylim([-0.4 0.4])
    ax1 = gca;
    ax1.YTick = -0.4:0.2:0.4;
    text(ylabelPosition(1), mean(y_pos), '\Delta standardized BOLD response', 'FontSize', 20, 'Rotation', 90, 'HorizontalAlignment', 'center');
    
    xticks(1:length(rois))
    
    % use the config's display names (e.g. 'V3A' instead of the on-disk
    % 'V3a', 'hMT+' instead of 'hMTcomplex') rather than hardcoding which
    % position needs relabeling, since that depended on a specific,
    % now-changed ROI list/order
    roinamesEdit = projectSettings.roi_plotnames;
    
    % Create a custom legend with the dummy plot objects
    legend([h1, h2], labelnames, 'Location', 'best');
    
    xticklabels(roinamesEdit); % Set x-axis tick labels
    
    f1 = gcf;
    f1.Position = [298 843 651 494];

    % Force all pending graphics updates (legend, xticklabels, position)
    % to render before capturing -- without this, print() can grab a
    % stale frame left over from the previous asymmetryName's figure
    % (observed: mainSubset/derivedSubset pair rendered identically).
    drawnow;

    % Save the figure as a TIFF file with specific options
    print(fullfile(figureDir, sprintf('LME_%s_%s_%s', comparisonName, projectName, asymmetryName)), '-dpdf', '-bestfit');

end


%% Make master figure

% Desired left-to-right order: vertical vs horizontal, cardinal vs oblique,
% radial vs tangential, polar cardinal vs polar oblique. asymmetryNames was
% collected in the fixed order {mainCardinalVsMainOblique,
% derivedCardinalVsDerivedOblique, mainSubset, derivedSubset}, but which of
% those is "cardinal vs oblique" vs "polar cardinal vs polar oblique" (and
% "radial vs tangential" vs "vertical vs horizontal") swaps between dg and
% da, so the permutation into the desired order is project-dependent.
if strcmp(projectName, 'dg')
    % desired order <- [mainSubset, mainCardinalVsMainOblique, derivedSubset, derivedCardinalVsDerivedOblique]
    plotOrder = [3, 1, 4, 2];
elseif strcmp(projectName, 'da')
    % desired order <- [derivedSubset, derivedCardinalVsDerivedOblique, mainSubset, mainCardinalVsMainOblique]
    plotOrder = [4, 2, 3, 1];
end

figure; hold on

x = 1:nA;
dx = 0; %0.18;

for ai = 1:nA
    srcIdx = plotOrder(ai); % index into box_pro/box_con/etc, in their original collection order

    % Pro
    boxchart((x(ai)-dx)*ones(size(box_pro{srcIdx})), box_pro{srcIdx}, ...
        'BoxFaceColor', color_pro(srcIdx,:), ...
        'LineWidth', 4, 'BoxWidth', 0.6);
    hold on

    errorbar(x(ai)-dx, mean_pro(srcIdx), ...
        errlow_pro(srcIdx), errhigh_pro(srcIdx), ...
        'LineStyle','none', 'LineWidth', 2, ...
        'Color', color_pro(srcIdx,:));

    % Con
    boxchart((x(ai)+dx)*ones(size(box_con{srcIdx})), box_con{srcIdx}, ...
        'BoxFaceColor', color_con(srcIdx,:), ...
        'LineWidth', 4, 'BoxWidth', 0.6);
    hold on

    errorbar(x(ai)+dx, mean_con(srcIdx), ...
        errlow_con(srcIdx), errhigh_con(srcIdx), ...
        'LineStyle','none', 'LineWidth', 2, ...
        'Color', color_con(srcIdx,:));
end

yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2)

ylim([-0.4 0.4])
xlim([0.5 nA+0.5])

set(gca,'XTick',[])
box off
set(gca,'linewidth',2, 'YColor',[0 0 0], 'XColor',[0 0 0]);
set(gca,'FontName','Arial','FontSize',20);

xticks(x)
xticklabels(asymLabel(plotOrder))
xtickangle(25)

ylabel('zscored BOLD psc', 'FontSize', 20);

% Legend (dummy handles, same as you do)
h1 = plot(nan, nan, 'Color', color_pro(1,:), 'LineWidth', 3);
h2 = plot(nan, nan, 'Color', color_con(1,:), 'LineWidth', 3);
%legend([h1 h2], {'pro','con'}, 'Location', 'best');

f = gcf;
f.Position = [298 843 651 494];

%print(fullfile(figureDir, sprintf('MASTER_ROI1_%s_%s', comparisonName, projectName)), ...
%    '-dtiff', '-r300');
%ylim([-0.4 0.4]) % for non zscore
print(fullfile(figureDir, sprintf('MASTER_ROI1_%s_%s', comparisonName, projectName)), '-dpdf', '-bestfit');


%% Per-ROI figures: same layout/colors as the master figure above (4
% asymmetries, pro vs con), but one figure per ROI instead of ROI 1 only.

for roi = 1:length(rois)

    figure; hold on

    for ai = 1:nA
        srcIdx = plotOrder(ai);

        % Pro
        boxchart(x(ai)-dx, mean_pro_allroi(srcIdx,roi), ...
            'BoxFaceColor', color_pro(srcIdx,:), ...
            'LineWidth', 4, 'BoxWidth', 0.6);
        hold on

        errorbar(x(ai)-dx, mean_pro_allroi(srcIdx,roi), ...
            errlow_pro_allroi(srcIdx,roi), errhigh_pro_allroi(srcIdx,roi), ...
            'LineStyle','none', 'LineWidth', 2, ...
            'Color', color_pro(srcIdx,:));

        % Con
        boxchart(x(ai)+dx, mean_con_allroi(srcIdx,roi), ...
            'BoxFaceColor', color_con(srcIdx,:), ...
            'LineWidth', 4, 'BoxWidth', 0.6);
        hold on

        errorbar(x(ai)+dx, mean_con_allroi(srcIdx,roi), ...
            errlow_con_allroi(srcIdx,roi), errhigh_con_allroi(srcIdx,roi), ...
            'LineStyle','none', 'LineWidth', 2, ...
            'Color', color_con(srcIdx,:));
    end

    yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2)

    ylim([-0.4 0.4])
    xlim([0.5 nA+0.5])

    set(gca,'XTick',[])
    box off
    set(gca,'linewidth',2, 'YColor',[0 0 0], 'XColor',[0 0 0]);
    set(gca,'FontName','Arial','FontSize',20);

    xticks(x)
    xticklabels(asymLabel(plotOrder))
    xtickangle(25)

    ylabel('zscored BOLD psc', 'FontSize', 20);
    title(projectSettings.roi_plotnames{roi}, 'FontSize', 20, 'Interpreter', 'none');

    % Legend (dummy handles, same as the master figure)
    h1 = plot(nan, nan, 'Color', color_pro(1,:), 'LineWidth', 3);
    h2 = plot(nan, nan, 'Color', color_con(1,:), 'LineWidth', 3);
    %legend([h1 h2], {'pro','con'}, 'Location', 'best');

    f = gcf;
    f.Position = [298 843 651 494];

    drawnow;

    print(fullfile(figureDir, sprintf('MASTER_%s_%s_%s', rois{roi}, comparisonName, projectName)), '-dpdf', '-bestfit');
end
