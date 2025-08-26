%comparisonName = 'motion_minus_orientation';
comparisonName = 'orientation_minus_baseline'; %'motion_minus_orientation';
projectName = 'dg';
glmResultsfolder = fullfile('/Volumes/server/Projects/Project_dg/data_bids/derivatives', ...
    strcat(projectName, 'GLM'), 'hRF_glmsingle');

asymmetryNames = {'mainCardinalVsMainOblique', 'derivedCardinalVsDerivedOblique', 'radialVsTangential'};

%all_labels = {{'Main Cardinal', 'Main Oblique'}, {'Derived Cardinal', 'Derived Oblique'}, {'Radial', 'Tangential'}};

ci_level = 68; %68;
meanRelative=0;
%colors = [[127 191 123]/255; [166 97 26]/255; [146 197 222]/255];
%colors2 = [[175 141 195]/255; [64 176 166]/255; [202 0 32]/255];

for roi=1:length(rois)
    figure
    xlim([0 8]);
    ylim([0 8]);



    for ai=1:numel(asymmetryNames)

        x = ai; % this is condition 1, 2, or 3

        asymmetryName = asymmetryNames{ai};

        colors = colors_data.conditions.(projectName).(asymmetryName).color_pro';
        colors2 = colors_data.conditions.(projectName).(asymmetryName).color_con';

        %labelnames = all_labels(x);
        labelnames = lower(strsplit(asymmetryNames{ai}, 'Vs'));

% here is where the for loop used to be
        
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
            boxchart(ai*1, y1(:,x)-baselineSub, 'BoxFaceColor', colors, 'LineWidth', 4, 'BoxWidth', 1)
            hold on
            %plot(roi*1, y1(:,x)-baselineSub, 'Color', colors, 'Marker', '.', 'MarkerSize', 10, 'LineStyle','none')
            %hold on
            errorbar(ai*1,y1(x)-baselineSub,errlow1(x), errhigh1(x), 'LineStyle','none', 'LineWidth', 2, 'Color', colors);
            hold on
            boxchart(ai*1, y2(:,x)-baselineSub, 'BoxFaceColor', colors2, 'LineWidth', 4, 'BoxWidth', 1)
            hold on
            errorbar(ai*1,y2(x)-baselineSub,errlow2(x), errhigh2(x), 'LineStyle','none', 'LineWidth', 2, 'Color', colors2);
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
        ylim([-0.3+Gintercept 0.3+Gintercept])
    end
    xlim([0.5 3.5])
    %ylim([]) % motion - baseline (DA)
    %ylim([0.2 0.7]) % orientation (DG)
    set(gca,'XTick',[])
    box on
    set(gca,'linewidth',1, 'YColor', [0 0 0]);
    set(gca,'linewidth',1, 'XColor', [0 0 0]);
    %title(condNames{x})
    set(gca, 'FontName', 'Arial', 'FontSize', 20);
    % ax1 = gca;
    % %ax1.YTick = [-0.03, -0.015, 0, 0.015, 0.03];
    % 
    % % yticks = get(gca, 'ytick'); % Get current y-axis tick positions
    % % text(0.5, mean(yticks), '\Delta', 'VerticalAlignment', 'baseline', 'FontSize', 20); % Add triangle symbol at (0.5, y) position
    % 
    % % Get current y-axis label position
    % ylabelHandle = ylabel('temp', 'FontSize', 20);
    % ylabelPosition = get(ylabelHandle, 'Position');
    % % Remove y-axis label
    % delete(ylabelHandle);
    % 
    % % Add triangle symbol in place of y-axis label
    % y_pos = ylim;
    % %text(ylabelPosition(1), mean(y_pos), '\Delta BOLD signal (%)', 'FontSize', 20, 'Rotation', 90, 'HorizontalAlignment', 'center');
    % 
    % ylim([-0.5 0.5])
    % ax1 = gca;
    % ax1.YTick = [-0.5, -0.25 0, .25, 0.5];
    % text(ylabelPosition(1), mean(y_pos), '\Delta standardized BOLD response', 'FontSize', 20, 'Rotation', 90, 'HorizontalAlignment', 'center');
    % 
    % xticks(1:length(rois))
    % 
    % roinamesEdit = rois;
    % roinamesEdit{5} = "hMT+"; % edit the MT complex ; it was too long as 'MTcomplex'
    % 
    % % Create a custom legend with the dummy plot objects
    % legend([h1, h2], labelnames, 'Location', 'best');
    % 
    % xticklabels(roinamesEdit); % Set x-axis tick labels
    % 
    % f1 = gcf;
    % f1.Position = [298 843 651 494];
    % 
    % % Save the figure as a TIFF file with specific options
    % %print(fullfile(figureDir, sprintf('LME_%s_%s_%s', comparisonName, projectName, asymmetryName)), '-dtiff', '-r300'); % '-r300' specifies a resolution of 300 DPI
    % 
    % disp('')
end

