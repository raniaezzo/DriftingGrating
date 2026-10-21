clc; clear all; close all

% set up
addpath(genpath(pwd));
projectName = 'da';
subj = 'ALL';
metric = 'meanBOLD';
bidsDir =  '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
githubDir = '~/Documents/GitHub';

% can be 'motion_minus_orientation' ; 'motion_minus_baseline' ; 'orientation_minus_baseline'
comparisonName = 'orientation_minus_baseline';

projectSettings = loadConfig(githubDir);

rois = projectSettings.rois;
axes_limits = projectSettings.axes_limits;
pairaxes_limits = projectSettings.pairaxes_limits;
pairaxes_PAew_limits = projectSettings.pairaxes_PAew_limits;
colors_data = projectSettings.colors_data;
contrasts_dict = projectSettings.contrasts_dict;

figureDir = [strrep(bidsDir, 'data_bids', 'figures'), projectName];

% temporary:
%path2data = fullfile('/Volumes','EXTERNAL_US','Project_dg','data_bids','derivatives', ...
path2data = fullfile(bidsDir,'derivatives', ...
    strcat(projectName,'GLM'),'hRF_glmsingle','LME_results', ...
    comparisonName, 'V1');

load(fullfile(path2data, 'modeldata.mat'))
load(fullfile(path2data, 'LME_bold.mat'))

% %% SIMULATE DATA
% 
% %Define location and direction values
% angles = 0:45:315;
% 
% % Generate all combinations of location and direction
% [locGrid, dirGrid] = ndgrid(angles, angles);
% 
% % Flatten the grids
% location = locGrid(:);
% direction = dirGrid(:);
% 
% % Simulate meanBOLD and medianBOLD values between 0 and 1
% numEntries = numel(location);
% 
% mu = 0.5;
% sigma = 0.1;  % Adjust for how tightly you want values around 0.5
% vals = mu + sigma * randn(numEntries,1);
% 
% % Clip to stay within [0, 1]
% vals(vals < 0) = 0;
% vals(vals > 1) = 1;
% 
% meanBOLD = vals;
% medianBOLD = vals;
% 
% % Create the table
% realdata = table(meanBOLD, medianBOLD, location, direction);
% 
% % Display the table
% disp(realdata);

%% combine table with all subjects into a metasubject

% Group by motiondir and polarangle, average/median across subjects
groupVars = {'motiondir','polarangle','mainCardinal','derivedCardinal','mainSubset', 'derivedSubset'}; %'polRadial'};

% Remove duplicate rows based on selected columns (important for
% orientation)
vars = {'bold','motiondir','polarangle','sub'}; % look for rows with these values identical
[~, ia] = unique(modeldata(:, vars), 'rows', 'stable'); % find unique rows w.r.t. those vars
proc_data = modeldata(ia, :); % keep only the first occurrence of each duplicate set
 
meanTbl = varfun(@mean, proc_data, ...
    'InputVariables','bold', ...
    'GroupingVariables',groupVars);

medianTbl = varfun(@median, proc_data, ...
    'InputVariables','bold', ...
    'GroupingVariables',groupVars);

% Merge the two results into one table
newTbl = join(meanTbl, medianTbl);

% Clean up variable names
newTbl.Properties.VariableNames{'mean_bold'} = 'meanBOLD';
newTbl.Properties.VariableNames{'median_bold'} = 'medianBOLD';

% change name do it applied to orientation and motion dir
newTbl.Properties.VariableNames{'motiondir'} = 'direction';
newTbl.Properties.VariableNames{'polarangle'} = 'location';

% if orientation, need to duplicate axial data
if strcmp(comparisonName, 'orientation_minus_baseline')
    dupData = newTbl; % Duplicate the whole table
    dupData.direction = dupData.direction + 180; % Add 180 to the dir values in the duplicate
    realdata = [newTbl; dupData]; % Append duplicate back to original
else
    realdata = newTbl;
end

%%

% load in a table AVEmodeldata that contains columns for average BOLD,
% location, and direction
%avemodeldata = readtable(strcat('/Users/rania/Desktop/RadialBias_pilot1/Analysis/Modeldata/AVE_modeldata.csv'));

figure
set(gcf, 'Position', [1118 87 1429 1250])
plot(0,0,'+k','MarkerSize',12, 'linewidth',3)
axis square
hold on
xlim([-4 4])
ylim([-4 4])
xticks([]) 
yticks([]) 
xticklabels({}) 
yticklabels({})
p = get(gca, 'Position');

box on
set(gca,'linewidth',1, 'YColor', [0 0 0]);
set(gca,'linewidth',1, 'XColor', [0 0 0]); 


w = p(3)-p(3)*.835; h = p(4)-p(4)*.8;
left = p(1)+p(3)/2-(w/2); bottom = p(2)+p(4)/2-(h/2);
origin = [left bottom w h];

straight = 1.6; diag = 1.1;
shifts = {[w*straight 0 0 0], [w*diag h*diag 0 0], [0 h*straight 0 0], ...
    [-w*diag h*diag 0 0], [-w*straight 0 0 0], [-w*diag -h*diag 0 0], ...
    [0 -h*straight 0 0], [w*diag -h*diag 0 0]};

N = [0 45 90 135 180 225 270 315];

%%

% initialize matrix to save out 8x8 model fits rho
modelPlotVals = nan(8,8); % each row is a location (0, 45, etc.) each col is a dir (0, 45, etc.)

dotsize = 8;

if strcmp(comparisonName, 'orientation_minus_baseline') || ...
    strcmp(comparisonName, 'motion_minus_orientation')
    globalMin = -0.5; %min(realdata.(metric));
    globalMax = 1; %max(realdata.(metric));
elseif strcmp(comparisonName, 'motion_minus_baseline')
    globalMin = 0; %min(realdata.(metric));
    globalMax = 2; %max(realdata.(metric));
else
    error('No condition matches this comparisonName')
end

for pos = 1:length(N) % location

    shift = shifts{pos};

    pax = polaraxes(gcf);
    set(pax, 'Position', origin+shift) %0.1300 0.1100
    
    % modelled data (retrieve design matrix
    % filter positions only using this for the dummy variables
    rows = realdata.location == N(pos);
    subTable = realdata(rows, :);

    % get the Cartesian-converted directions from da
    if strcmp(projectName, 'da')
        % input = 
        % direction 
        % polarangles (N = [0, 45, 90, 135, 180, 225, 270, 315])
        for i = 1:height(subTable)
            remappingDirIdx = deriveLocalMotionfromUVM(subTable.direction(i), N);
            subTable.direction(i) = remappingDirIdx(pos);
        end
    end

    convSaved = subTable.direction(:); 
    subTable = sortrows(subTable, 'direction'); % sort by direction in case not already

    M = [ones(height(subTable),1), ...
         subTable.mainCardinal, ...
         subTable.derivedCardinal, ...
         subTable.mainSubset, ...
         subTable.derivedSubset];
         %subTable.polRadial];

    estimates = estimates(:); % ensure its a col

    temp = M * estimates;

    %temp = [.5 .5 .5 .5 .5 .5 .5 .5]'; %recon.sum(idx(pos):idx(pos)+7); 
    rho = temp';
    modelPlotVals(pos,:) = rho;
    rho = [rho rho(1)];
    
    polarplot(deg2rad([0 45 90 135 180 225 270 315 0]), rho, 'r', 'linewidth', 2)
    hold on
    q=polarplot(deg2rad([0 45 90 135 180 225 270 315 0]), rho, 'o', 'LineWidth', 2, 'MarkerFaceColor', 'r', 'MarkerSize', dotsize);
    q.MarkerEdgeColor = 'w';
    hold on
    
    % Get the BOLD values for those 8 directions at those locations
    targetLocation = N(pos);  % Replace with the location you want (e.g. 90)
    
    % Find rows matching the target location
    rows = realdata.location == targetLocation;
    
    % Extract direction values and sort them
    dirVals = realdata.direction(rows);
    boldVals = realdata{rows, sprintf(metric)};

    % Sort by direction
    [sortedDir, sortIdx] = sort(dirVals);
    bold_sortedbyDir = boldVals(sortIdx);

    if strcmp(projectName, 'da')
        [~, idx] = ismember(sortedDir, convSaved);
        bold_sortedbyDir = bold_sortedbyDir(idx); % re-organize to Cartesian. Just the BOLD values, leave direction as sorted.
    end
    
    rho = bold_sortedbyDir';
    %rho = [rho(d_idx(1)), rho(d_idx(2)), rho(d_idx(3)), rho(d_idx(4)), ...
    %    rho(d_idx(5)), rho(d_idx(6)), rho(d_idx(7)), rho(d_idx(8))];
    
    p = polarplot(deg2rad([0 45 90 135 180 225 270 315]), rho, 'ok', 'LineWidth', 1, 'MarkerFaceColor', 'black', 'MarkerSize', dotsize/1.5); %8); %4);
    p.MarkerEdgeColor = 'w';
    pax.FontSize = 6;
    pax.RTickLabel = {''};
    pax.ThetaTickLabel = {''};
    thetaticks(0:45:315);
    rlim([globalMin globalMax])
    hold on

    
    
end
 
sgtitle(['Allsubs AbscardRelcardRadial', sprintf(' est %s %s %s',projectName, strrep(comparisonName, '_', '-'), metric)])

% % This is to compute R^2 of model predictions for mean data
% pred = recon.sum;
% actual = avemodeldata.mean_sensitivity;
% MSE = (1/(64))*(sum((actual-pred).^2));
% R2 = 1 - (MSE/var(actual)); % only for mean data.
% 
% % marginal vs. conditional R2 (per datapoint)
% marg = fitted(lme, 'Conditional', false);
% cond = fitted(lme);
% actual_all = modeldata.(metric);
% 
% % actual r^2 values
% r_squared_from_cond = 1-(sum((actual_all-cond).^2)/sum((actual_all-mean(actual_all)).^2));
% r_squared_from_marg = 1-(sum((actual_all-marg).^2)/sum((actual_all-mean(actual_all)).^2));
% 
% cond_var_ratio = var(cond)/var(actual_all);
% marg_var_ratio = var(marg)/var(actual_all);


% saveas(gcf,strcat(fullfile(figureDir, projectName), ['Allsubs_AbscardRelcardRadial_4params', sprintf('_est_%s_%s_%s',projectName, ...
%     comparisonName, metric)]), 'png') %'epsc')

print(gcf,strcat(fullfile(figureDir, projectName), ['Allsubs_AbscardRelcardRadial_4params', sprintf('_est_%s_%s_%s',projectName, ...
    comparisonName, metric)]), '-dpdf', '-bestfit');

%% Save the plotted data, and save the condition tags (1,-1, 0) in abs dir/loc
% note; the tag matrix will be the same for da and dg, since this is saving
% the rho

save(fullfile(path2data, 'modelPlotVals.mat'), 'modelPlotVals')

tags = nan(8,8, 4);

% ROW
% per dir/orientation, for abs 0 to 45, to ... 315

% COL
% per location, from 0 to 45, to ... 315

% ASYM (1-4)
% hor (1), ver (-1), other (0)
tmp = [1,0,-1,0,1,0,-1,0];
M = repmat(tmp, 8, 1);
tags(:,:,1) = int8(M);

% card (1), obli (-1), other (0)
M = int8(M ~= 0);
M(M == 0) = -1;
tags(:,:,2) = M;

% radial (1), tang (-1), other (0)
tmp = [1,0,-1,0,1,0,-1,0];
M = zeros(8, numel(tmp));

for ii = 1:8
    M(ii,:) = circshift(tmp, ii-1);
end
tags(:,:,3) = int8(M);

% polar card (1), polar oblique (-1), other (0)
M = int8(M ~= 0);
M(M == 0) = -1;
tags(:,:,4) = M;

genpath = fullfile(bidsDir,'derivatives', ...
    strcat(projectName,'GLM'),'hRF_glmsingle','LME_results', ...
    comparisonName);
save(fullfile(genpath, 'modelPlotValsTags.mat'), 'tags')

