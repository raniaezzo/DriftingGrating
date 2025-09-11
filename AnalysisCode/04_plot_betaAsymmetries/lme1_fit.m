clc; clear all; close all

% set up
addpath(genpath(pwd));
projectName = 'dg';
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
               % (condition, polar angle, roi, subject)
               reshaped_mat(index, 1) = filtered_meanBOLDpa(md, pa, roi, subject);
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

    % initialize value for each asymmetry
    mainCardinal = zeros(length(finalMat),1);
    derivedCardinal = zeros(length(finalMat),1);
    polRadial = zeros(length(finalMat),1);

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
        polRadial_idx = (abs(finalMat(:,2)-finalMat(:,3)) == 0 | abs(finalMat(:,2)-finalMat(:,3)) == 180);
        polTangential_idx = abs(finalMat(:,2)-finalMat(:,3)) == 90 | abs(finalMat(:,2)-finalMat(:,3)) == 270;
    elseif strcmp(projectName, 'da')
        polRadial_idx = ismember(finalMat(:,2), [90, 270]); % for da, this is in/out
        polTangential_idx = ismember(finalMat(:,2), [0, 180]); % for da, this is c/cc
    end

    polRadial(polRadial_idx) = 1;
    polRadial(polTangential_idx) = -1;

    finalMat = [finalMat, mainCardinal, derivedCardinal, polRadial];

    %%

    variable_names = {'bold', 'motiondir', 'polarangle', 'sub', 'mainCardinal', 'derivedCardinal', 'polRadial'};
    modeldata = array2table(finalMat, 'VariableNames', variable_names);

    modeldata.subject = categorical(modeldata.sub);

    lme = fitlme(modeldata, 'bold ~ mainCardinal + derivedCardinal + polRadial + (1|sub)');

    anova(lme, 'dfmethod', 'satterthwaite')
    
    global_idx = find(contains(cellstr(lme.Coefficients(:,1)), 'Intercept'));
    maincard_idx = find(contains(cellstr(lme.Coefficients(:,1)), 'mainCardinal'));
    derivedcard_idx = find(contains(cellstr(lme.Coefficients(:,1)), 'derivedCardinal'));
    radial_idx = find(contains(cellstr(lme.Coefficients(:,1)), 'polRadial'));

    global_est = double(lme.Coefficients(global_idx,2));
    maincard_est = double(lme.Coefficients(maincard_idx,2));
    derivedcard_est = double(lme.Coefficients(derivedcard_idx,2));
    radial_est = double(lme.Coefficients(radial_idx,2));

    estimates = [global_est, maincard_est, derivedcard_est, radial_est];

    save(fullfile(saveDir,strcat('LME_',metric)), 'estimates');

    disp(lme)
    save(strcat(saveDir, '/modeldata'), 'modeldata');
end

%% bootstrap data

subjects = {'ALL'};
 
rng('default'); rng(1);

subjID = 1:nSubjs; %10;
bootN = 1000;

for roi=1:length(rois)

    roi

    saveDir = fullfile(glmResultsfolder,'LME_results', comparisonName, rois{roi});
    
    saveboot = {};
    coeffs = nan(4,bootN);
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

        lme = fitlme(randomsample,'bold ~ mainCardinal + derivedCardinal + polRadial + (1|sub)');

        saveboot{bi} = lme;

        global_idx = find(contains(cellstr(lme.Coefficients(:,1)), 'Intercept'));
        maincard_idx = find(contains(cellstr(lme.Coefficients(:,1)), 'mainCardinal'));
        derivedcard_idx = find(contains(cellstr(lme.Coefficients(:,1)), 'derivedCardinal'));
        radial_idx = find(contains(cellstr(lme.Coefficients(:,1)), 'polRadial'));

        global_est = double(lme.Coefficients(global_idx,2));
        maincard_est = double(lme.Coefficients(maincard_idx,2));
        derivedcard_est = double(lme.Coefficients(derivedcard_idx,2));
        radial_est = double(lme.Coefficients(radial_idx,2));

        estimates = [global_est, maincard_est, derivedcard_est, radial_est];
        coeffs(:,bi) = estimates;
        clear lme
        disp(bi)

    end

    save(strcat(saveDir, '/boot'), 'saveboot','coeffs');
end



%%

asymmetryNames = {'mainCardinalVsMainOblique', 'derivedCardinalVsDerivedOblique', 'radialVsTangential'};

%all_labels = {{'Main Cardinal', 'Main Oblique'}, {'Derived Cardinal', 'Derived Oblique'}, {'Radial', 'Tangential'}};

ci_level = 68; %68;
meanRelative=1;
%colors = [[127 191 123]/255; [166 97 26]/255; [146 197 222]/255];
%colors2 = [[175 141 195]/255; [64 176 166]/255; [202 0 32]/255];

for ai=1:numel(asymmetryNames)

    x = ai; % this is condition 1, 2, or 3

    asymmetryName = asymmetryNames{ai};

    colors = colors_data.conditions.(projectName).(asymmetryName).color_pro';
    colors2 = colors_data.conditions.(projectName).(asymmetryName).color_con';

    %labelnames = all_labels(x);
    labelnames = lower(strsplit(asymmetryNames{ai}, 'Vs'));

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
        radial_est = Gintercept + estimates(4);
    
        CIFcn = @(x,p)prctile(x, [100-p, p]); p = ci_level;
    
        % coefficient order is: intercept, abs_cardinality, rel_cardinality, radiality
        CI_maincardinality = CIFcn(coeffs(2,:)+Gintercept,p);
        CI_derivedcardinality = CIFcn(coeffs(3,:)+Gintercept,p);
        CI_radiality = CIFcn(coeffs(4,:)+Gintercept,p);
    
    
        y = [estimates(2) estimates(3) estimates(4)]; %[0.11446, 0.033442, 0.015738];
        %errlow = [0.0053434, 0.0053434, 0.0075567]; % output from model
        %errhigh = [0.0053434, 0.0053434, 0.0075567];
    
        errlow = [CI_maincardinality(1) CI_derivedcardinality(1) CI_radiality(1)];
        errhigh = [CI_maincardinality(2) CI_derivedcardinality(2) CI_radiality(2)];
    
        y1 = Gintercept + y;
        y2 = Gintercept - y;
    
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
    xlim([0 8])
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
    
    ylim([-0.5 0.5])
    ax1 = gca;
    ax1.YTick = [-0.5, -0.25 0, .25, 0.5];
    text(ylabelPosition(1), mean(y_pos), '\Delta standardized BOLD response', 'FontSize', 20, 'Rotation', 90, 'HorizontalAlignment', 'center');
    
    xticks(1:length(rois))
    
    roinamesEdit = rois;
    roinamesEdit{5} = "hMT+"; % edit the MT complex ; it was too long as 'MTcomplex'
    
    % Create a custom legend with the dummy plot objects
    legend([h1, h2], labelnames, 'Location', 'best');
    
    xticklabels(roinamesEdit); % Set x-axis tick labels
    
    f1 = gcf;
    f1.Position = [298 843 651 494];
    
    % Save the figure as a TIFF file with specific options
    print(fullfile(figureDir, sprintf('LME_%s_%s_%s', comparisonName, projectName, asymmetryName)), '-dtiff', '-r300'); % '-r300' specifies a resolution of 300 DPI

end





