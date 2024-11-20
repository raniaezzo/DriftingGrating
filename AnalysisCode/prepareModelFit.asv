clc; clear all; close all

%%


% motion directions
    % 8 - m0_v_s90
    % 9 - m90_v_s0
    % 10 - m180_v_s90
    % 11 - m270_v_s0
    % 12 - m45_v_s135
    % 13 - m135_v_s45
    % 14 - m225_v_s135
    % 15 - m315_v_s45

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

% polar angle
    anglevals = [90; 45; 0; 315; 270; 225; 180; 135];
% rois
% subjects

cardinalmDir = [0; 90; 180; 270];
primaryMeridians = [90; 0; 270; 180];

% motion - orientation
%filtered_meanBOLDpa = meanBOLDpa(8:15, :, :, :);
%mdirvals = [0; 90; 180; 270; 45; 135; 225; 315];

% motion - baseline
% filtered_meanBOLDpa = meanBOLDpa(18:25, :, :, :);
% mdirvals = [0; 180; 90; 270; 45; 225; 135; 315];

% orientation - baseline
filtered_meanBOLDpa = meanBOLDpa(26:29, :, :, :);
mdirvals = [0; 90; 45; 135];

[nMotDirs, nPAs, nROIs, nSubjs] = size(filtered_meanBOLDpa);

metric = 'bold';
roinames = {'V1','V2','V3','hV4','MTcomplex', 'MT', 'MST'};

for roi=1:length(roinames)
    
    saveDir = fullfile('/Users','rje257','Desktop','LME_results', roinames{roi});
    
    if ~isfolder(saveDir)
        mkdir(saveDir)
    end
    
    % Loop through each element in the new matrix
    index = 1;

    % Preallocate the new matrix
    reshaped_mat = zeros(nMotDirs*nPAs*nSubjs, 1); % just 1 ROI (makes interprettability easier)

    dirCol = repmat(mdirvals, nPAs*nSubjs, 1);
    paRep = repelem(anglevals, nPAs);
    paCol = repmat(paRep, nSubjs, 1);
    subCol = repelem(1:nSubjs, 64)';

    for subject = 1:nSubjs
        for pa = 1:nPAs
            for md = 1:nMotDirs
               reshaped_mat(index, 1) = filtered_meanBOLDpa(md, pa, roi, subject);
               index = index + 1;
            end
        end
    end

    %%
    % now determine which are cartesian cardinal, polar cardinal, and radial:

    % just for orientation
    reshaped_mat = repelem(reshaped_mat, 2);
    dirCol = repelem(dirCol, 2);

    finalMat = [reshaped_mat, dirCol, paCol, subCol];

    cartCardinal = zeros(length(finalMat),1);
    polCardinal = zeros(length(finalMat),1);
    polRadial = zeros(length(finalMat),1);

    cartCardinal_idx = ismember(finalMat(:,2), cardinalmDir);
    cartCardinal(cartCardinal_idx) = 1;
    cartCardinal(~cartCardinal_idx) = -1;

    polCardinal_idx = (ismember(finalMat(:,2), cardinalmDir) & ismember(finalMat(:,3), primaryMeridians)) ...
        | ((~ismember(finalMat(:,2), cardinalmDir) & (~ismember(finalMat(:,3), primaryMeridians))));
    polCardinal(polCardinal_idx) = 1;
    polCardinal(~polCardinal_idx) = -1;

    polRadial_idx = (abs(finalMat(:,2)-finalMat(:,3)) == 0 | abs(finalMat(:,2)-finalMat(:,3)) == 180);
    polRadial(polRadial_idx) = 1;
    polTangential = abs(finalMat(:,2)-finalMat(:,3)) == 90 | abs(finalMat(:,2)-finalMat(:,3)) == 270;
    polRadial(polTangential) = -1;

    finalMat = [finalMat, cartCardinal, polCardinal, polRadial];

    %%

    variable_names = {'bold', 'motiondir', 'polarangle', 'sub', 'cartCardinal', 'polCardinal', 'polRadial'};
    modeldata = array2table(finalMat, 'VariableNames', variable_names);

    modeldata.subject = categorical(modeldata.sub);

    lme = fitlme(modeldata, 'bold ~ cartCardinal + polCardinal + polRadial + (1|sub)');

    anova(lme, 'dfmethod', 'satterthwaite')
    
    global_idx = find(contains(cellstr(lme.Coefficients(:,1)), 'Intercept'));
    abscard_idx = find(contains(cellstr(lme.Coefficients(:,1)), 'cartCardinal'));
    relcard_idx = find(contains(cellstr(lme.Coefficients(:,1)), 'polCardinal'));
    radial_idx = find(contains(cellstr(lme.Coefficients(:,1)), 'polRadial'));

    global_est = double(lme.Coefficients(global_idx,2));
    abscard_est = double(lme.Coefficients(abscard_idx,2));
    relcard_est = double(lme.Coefficients(relcard_idx,2));
    radial_est = double(lme.Coefficients(radial_idx,2));

    estimates = [global_est, abscard_est, relcard_est, radial_est];

    save(fullfile(saveDir,strcat('LME_',metric)), 'estimates');

    disp(lme)
    save(strcat(saveDir, '/modeldata'), 'modeldata');
end

%% bootstrap data

subjects = {'ALL'};
 
rng('default'); rng(1);

subjID = 1:10;
bootN = 1000;


roinames = {'V1','V2','V3','hV4','MTcomplex', 'MT', 'MST'};

for roi=1:length(roinames)
    
    saveboot = {};
    coeffs = nan(4,bootN);
    saveDir = fullfile('/Users','rje257','Desktop','LME_results', roinames{roi});
    load(strcat(saveDir, '/modeldata'), 'modeldata');

    for bi=1:bootN

        y = datasample(subjID,10); % take a random sample of subjects

        %modeldata = readtable(strcat('/Users/rania/Desktop/RadialBias_pilot1/Data/ALLSUBSwFULLDATA/Analysis_1_PF/datatable.csv'));
        
        randomsample = [];
        for i=1:length(y)
            selectsub = y(i);
            temp = modeldata(modeldata.sub == selectsub,:);
            randomsample = [randomsample ; temp];
        end

        lme = fitlme(randomsample,'bold ~ cartCardinal + polCardinal + polRadial + (1|sub)');

        saveboot{bi} = lme;

        global_idx = find(contains(cellstr(lme.Coefficients(:,1)), 'Intercept'));
        abscard_idx = find(contains(cellstr(lme.Coefficients(:,1)), 'cartCardinal'));
        relcard_idx = find(contains(cellstr(lme.Coefficients(:,1)), 'polCardinal'));
        radial_idx = find(contains(cellstr(lme.Coefficients(:,1)), 'polRadial'));

        global_est = double(lme.Coefficients(global_idx,2));
        abscard_est = double(lme.Coefficients(abscard_idx,2));
        relcard_est = double(lme.Coefficients(relcard_idx,2));
        radial_est = double(lme.Coefficients(radial_idx,2));

        estimates = [global_est, abscard_est, relcard_est, radial_est];
        coeffs(:,bi) = estimates;
        clear lme
        disp(bi)

    end

    save(strcat(saveDir, '/boot'), 'saveboot','coeffs');
end

%%

ci_level = 68; %68;
colors = [[127 191 123]/255; [166 97 26]/255; [146 197 222]/255];
colors2 = [[175 141 195]/255; [64 176 166]/255; [202 0 32]/255];
meanRelative=0; %1;

for roi=1:length(roinames)
    % load('bootLME_sensitivity.mat')
    % load('LME_sensitivity.mat')
    saveDir = fullfile('/Users','rje257','Desktop','LME_results', roinames{roi});
    load(strcat(saveDir, '/modeldata'), 'modeldata');
    load(fullfile(saveDir,strcat('LME_',metric)), 'estimates');
    load(strcat(saveDir, '/boot'), 'saveboot','coeffs');

    Gintercept = estimates(1);
    abs_cardinal_est = Gintercept + estimates(2);
    rel_cardinal_est = Gintercept + estimates(3);
    radial_est = Gintercept + estimates(4);

    CIFcn = @(x,p)prctile(x, [100-p, p]); p = ci_level;

    % coefficient order is: intercept, abs_cardinality, rel_cardinality, radiality
    CI_abscardinality = CIFcn(coeffs(2,:)+Gintercept,p);
    CI_relcardinality = CIFcn(coeffs(3,:)+Gintercept,p);
    CI_radiality = CIFcn(coeffs(4,:)+Gintercept,p);

    
    %
    
%     figure
% 
%     subplot(3,1,1)
%     xline(abs_cardinal_est, 'k--', 'linewidth', 4)
%     hold on
%     histogram(coeffs(2,:)+Gintercept, 'FaceColor', colors(1,:))
%     hold on
%     xline(CI_abscardinality(1), 'r--', 'linewidth', 4)
%     hold on
%     xline(CI_abscardinality(2), 'r--', 'linewidth', 4)
%     hold on
%     xline(Gintercept, 'y--', 'linewidth', 4)
%     title('abs cardinal')
%     %xlim([xmin xmax])
% 
%     subplot(3,1,2)
%     xline(rel_cardinal_est, 'k--', 'linewidth', 4)
%     hold on
%     histogram(coeffs(3,:)+Gintercept, 'FaceColor', colors(2,:))
%     hold on
%     xline(CI_relcardinality(1), 'r--', 'linewidth', 4)
%     hold on
%     xline(CI_relcardinality(2), 'r--', 'linewidth', 4)
%     hold on
%     xline(Gintercept, 'y--', 'linewidth', 4)
%     title('rel cardinal')
%     %xlim([xmin xmax])
% 
%     subplot(3,1,3)
%     xline(radial_est, 'k--', 'linewidth', 4)
%     hold on
%     histogram(coeffs(4,:)+Gintercept, 'FaceColor', colors(3,:))
%     hold on
%     xline(CI_radiality(1), 'r--', 'linewidth', 4)
%     hold on
%     xline(CI_radiality(2), 'r--', 'linewidth', 4)
%     hold on
%     xline(Gintercept, 'y--', 'linewidth', 4)
%     title('radial')
%     %xlim([xmin xmax])
% 
%     sgtitle(sprintf('mean betas + %s CIs', num2str(ci_level)))

    % sensitivity (manuscript figure)

    y = [estimates(2) estimates(3) estimates(4)]; %[0.11446, 0.033442, 0.015738];
    %errlow = [0.0053434, 0.0053434, 0.0075567]; % output from model
    %errhigh = [0.0053434, 0.0053434, 0.0075567];

    errlow = [CI_abscardinality(1) CI_relcardinality(1) CI_radiality(1)];
    errhigh = [CI_abscardinality(2) CI_relcardinality(2) CI_radiality(2)];

    y1 = Gintercept + y;
    y2 = Gintercept - y;

    errlow1 = y1-errlow;
    errhigh1 = errhigh -y1;

    errlow2 = errhigh1; %errlow1;
    errhigh2 = errlow1; %errhigh1;

    x = 1:3;
    figure();
    ax = axes();
    hold(ax);
    if meanRelative
        baselineSub = Gintercept;
    else
        baselineSub = 0;
    end
    for i=1:length(x)
        boxchart(x(i)*ones(size(y1(:,i))), y1(:,i)-baselineSub, 'BoxFaceColor', colors(i,:), 'LineWidth', 3, 'BoxWidth', 1)
        plot(x(i)*ones(size(y1(:,i))), y1(:,i)-baselineSub, 'Color', colors(i,:), 'Marker', '.', 'MarkerSize', 10, 'LineStyle','none')
        hold on
        errorbar(x(i),y1(i)-baselineSub,errlow1(i)-baselineSub, errhigh1(i)-baselineSub, 'LineStyle','none', 'LineWidth', 2, 'Color', colors(i,:));
        hold on
        boxchart(x(i)*ones(size(y2(:,i))), y2(:,i)-baselineSub, 'BoxFaceColor', colors2(i,:), 'LineWidth', 3, 'BoxWidth', 1)
        hold on
        errorbar(x(i),y2(i)-baselineSub,errlow2(i)-baselineSub, errhigh2(i)-baselineSub, 'LineStyle','none', 'LineWidth', 2, 'Color', colors2(i,:));
        hold on
    end
    if meanRelative
        yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2)
        ylim([-0.05 0.05])
    else
        yline(Gintercept, '--', 'Color', [0 0 0], 'LineWidth', 2)
        ylim([-0.05+Gintercept 0.05+Gintercept])
    end
    xlim([0.5 3.5])
    %ylim([ymin ymax])
    set(gca,'XTick',[])
    box off
    set(gca,'linewidth',2, 'YColor', [0 0 0]);
    set(gca,'linewidth',2, 'XColor', [0 0 0]);
    title(roinames(roi))
    set(gca, 'FontName', 'Arial', 'FontSize', 12);
    
    
end
hold off

%%

ci_level = 68; %68;
colors = [[127 191 123]/255; [166 97 26]/255; [146 197 222]/255];
colors2 = [[175 141 195]/255; [64 176 166]/255; [202 0 32]/255];
meanRelative=1;
condNames = {'cart cardinal', 'polar cardinal', 'radial'};

x = 3; %:3; % condition

all_labels = {{'Cart Cardinal', 'Cart Oblique'}, {'Pol Cardinal', 'Pol Oblique'}, {'Radial', 'Tangential'}};

labelnames = all_labels(x);

figure
xlim([0 8]);
ylim([0 8]);
for roi=1:length(roinames)
    
    % load('bootLME_sensitivity.mat')
    % load('LME_sensitivity.mat')
    saveDir = fullfile('/Users','rje257','Desktop','LME_results', roinames{roi});
    load(strcat(saveDir, '/modeldata'), 'modeldata');
    load(fullfile(saveDir,strcat('LME_',metric)), 'estimates');
    load(strcat(saveDir, '/boot'), 'saveboot','coeffs');

    Gintercept = estimates(1);
    abs_cardinal_est = Gintercept + estimates(2);
    rel_cardinal_est = Gintercept + estimates(3);
    radial_est = Gintercept + estimates(4);

    CIFcn = @(x,p)prctile(x, [100-p, p]); p = ci_level;

    % coefficient order is: intercept, abs_cardinality, rel_cardinality, radiality
    CI_abscardinality = CIFcn(coeffs(2,:)+Gintercept,p);
    CI_relcardinality = CIFcn(coeffs(3,:)+Gintercept,p);
    CI_radiality = CIFcn(coeffs(4,:)+Gintercept,p);


    y = [estimates(2) estimates(3) estimates(4)]; %[0.11446, 0.033442, 0.015738];
    %errlow = [0.0053434, 0.0053434, 0.0075567]; % output from model
    %errhigh = [0.0053434, 0.0053434, 0.0075567];

    errlow = [CI_abscardinality(1) CI_relcardinality(1) CI_radiality(1)];
    errhigh = [CI_abscardinality(2) CI_relcardinality(2) CI_radiality(2)];

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
        boxchart(roi*1, y1(:,x)-baselineSub, 'BoxFaceColor', colors(x,:), 'LineWidth', 4, 'BoxWidth', .75)
        hold on
        plot(roi*1, y1(:,x)-baselineSub, 'Color', colors(x,:), 'Marker', '.', 'MarkerSize', 10, 'LineStyle','none')
        hold on
        errorbar(roi*1,y1(x)-baselineSub,errlow1(x), errhigh1(x), 'LineStyle','none', 'LineWidth', 2, 'Color', colors(x,:));
        hold on
        boxchart(roi*1, y2(:,x)-baselineSub, 'BoxFaceColor', colors2(x,:), 'LineWidth', 4, 'BoxWidth', .75)
        hold on
        errorbar(roi*1,y2(x)-baselineSub,errlow2(x), errhigh2(x), 'LineStyle','none', 'LineWidth', 2, 'Color', colors2(x,:));
        hold on
    %end
    
    hold on

end

hold on
% Create dummy plot objects for legend
h1 = plot(nan, nan, 'Color', colors(x,:), 'LineWidth', 3); % Dummy red line
hold on;
h2 = plot(nan, nan, 'Color', colors2(x,:), 'LineWidth', 3); % Dummy blue line

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

xticks(1:length(roinames))

roinamesEdit = roinames;
roinamesEdit{5} = "hMT+"; % edit the MT complex ; it was too long as 'MTcomplex'

% Create a custom legend with the dummy plot objects
legend([h1, h2], labelnames{1}, 'Location', 'best');

xticklabels(roinamesEdit); % Set x-axis tick labels

f1 = gcf;
f1.Position = [298 843 651 494];

[parentDirectory, ~, ~] = fileparts(saveDir);
figDir = fullfile(parentDirectory, 'figures');

% Save the figure as a TIFF file with specific options
print(fullfile(figDir, condNames{x}), '-dtiff', '-r300'); % '-r300' specifies a resolution of 300 DPI