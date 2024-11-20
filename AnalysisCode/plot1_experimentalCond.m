function plot1_experimentalCond(medianBOLDpa, asymmetryName, projectSettings, varargin)

    % this function will (polar plot) the averages of experimental conditions that
    % do not need to be derived from polar angle.
    % For dg experiment, this will include cardinal v oblique
    % For da experiment, this will include polar cardinal v polar oblique
    %                                   and radial v tangential

    % Check project name & request
    if strcmp(projectSettings.projectName, 'da')
        if strcmp(asymmetryName, 'mainCardinalVsMainOblique')
            derivedVals = 0; % this is needed because radialVsTang occurs for both conditions
            % Ensure a third argument is provided
            if nargin < 3 || isempty(varargin{1})
                error(sprintf('A third input is required when projectName is "%s" and asymmetryName is "%s".', projectSettings.projectName, asymmetryName));
            else
                radialvstang = varargin{1};
                if radialvstang == 1
                    asymmetryName = 'radialVsTangential';
                end
            end
        elseif strcmp(asymmetryName, 'derivedCardinalVsDerivedOblique')
             % Third argument is ignored
             derivedVals = 1;
             radialvstang = 0;
        end
    elseif strcmp(projectSettings.projectName, 'dg')
        if strcmp(asymmetryName, 'mainCardinalVsMainOblique')
            % Third argument is ignored
            derivedVals = 0; 
            radialvstang = 0;
        elseif strcmp(asymmetryName, 'derivedCardinalVsDerivedOblique')
            derivedVals = 1;
            % Ensure a third argument is provided
            if nargin < 3 || isempty(varargin{1})
                error(sprintf('A third input is required when projectName is "%s" and asymmetryName is "%s".', projectSettings.projectName, asymmetryName));
            else
                radialvstang = varargin{1};
                if radialvstang == 1
                    asymmetryName = 'radialVsTangential';
                end
            end
        end
    else
        error('Unknown projectName. Expected "da" or "dg".');
    end

    projectName = projectSettings.projectName;
    colors_data = projectSettings.colors_data;
    rois = projectSettings.rois;
    comparisonName = projectSettings.comparisonName;
    axes_limits = projectSettings.axes_limits;
    figureDir = projectSettings.figureDir;
    
    colors = {colors_data.conditions.(projectName).(asymmetryName).color_pro', ...
        colors_data.conditions.(projectName).(asymmetryName).color_con'};

    markerC = colors_data.conditions.(projectName).plotSettings.markerC;

    % retrieve the indices for specific asymmetries (e.g., motion -
    % orientation for main cardinal v main oblique)
    if ~derivedVals
        [proConditions, conConditions, ~] = retrieveProConIdx(projectName, comparisonName, radialvstang);
    else
        proConditions = 1; conConditions = 2;
    end
    
    figure
    gap = [.04 .01]; % spacing between the subplots vertical gap - horizontal gap
    marg_h = [.015 .13]; % margins of bottom - top of the figure
    marg_w = [.01 .01]; % margina - left, right
    [ha, pos] = tight_subplot(2, 4, gap, marg_h, marg_w);
    % Hide the 8th subplot
    emptyPlots = 8 - numel(rois);
    for ep=1:emptyPlots
        empty_idx = 8+1-ep;
        set(ha(empty_idx), 'Visible', 'off');
    end
    for ri=1:length(rois)
        % Specify the region index
        regionIndex = ri;
        
        % Extract the relevant conditions for the specified region
        conditions1 = medianBOLDpa(proConditions, :, regionIndex, :);
        conditions2 = medianBOLDpa(conConditions, :, regionIndex, :);
        conditions1 = squeeze(conditions1);
        conditions2 = squeeze(conditions2);
        
        % Average across the conditions within subjects
        if strcmp(projectName, 'da') && strcmp(comparisonName, 'orientation_minus_baseline') && radialvstang == 1 || ...
                strcmp(asymmetryName, 'derivedCardinalVsDerivedOblique')
            avgConditions1 = conditions1;  % for orientation DA - there is only 1 value
            avgConditions2 = conditions2; 
        else
            avgConditions1 = nanmean(conditions1, 1);
            avgConditions2 = nanmean(conditions2, 1);
            avgConditions1 = squeeze(avgConditions1);
            avgConditions2 = squeeze(avgConditions2);
        end
        
        % Extract polar angles
        %anglevals = [90, 135, 180, 225, 270, 315, 0, 45];
        anglevals = [90, 45, 0, 315, 270, 225, 180, 135]; % <-- these were manually converted based on the order of polarAngles above (Noah's convention)
        
        vals_1 = nanmean(avgConditions1,2)';
        sem1 = nanstd(avgConditions1,0,2)' ./ sqrt(sum(~isnan(avgConditions1),2)');
        vals_2 = nanmean(avgConditions2,2)';
        sem2 = nanstd(avgConditions2,0,2)' ./ sqrt(sum(~isnan(avgConditions2),2)');
 
        % using this method instead of subplot for TIGHT AXES
        axes(ha(ri)); 
        
        
        % plot average (connecting the last line)
        polarplot([deg2rad(anglevals(end)), deg2rad(anglevals(1))],[vals_1(end), vals_1(1)], 'o-', 'Color', colors{1}, 'MarkerSize', 12, 'LineWidth',1.75,  'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', markerC)
        hold on
        polarplot([deg2rad(anglevals(end)), deg2rad(anglevals(1))],[vals_2(end), vals_2(1)], 'o-', 'Color', colors{2}, 'MarkerSize', 12, 'LineWidth',1.75,  'MarkerFaceColor', colors{2}, 'MarkerEdgeColor', markerC)
        hold on
        
        % plot average (all the data)
        polarplot(deg2rad(anglevals),vals_1, 'o-', 'Color', colors{1},  'MarkerSize', 12, 'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', markerC,'LineWidth',1.75)
        hold on
        polarplot(deg2rad(anglevals),vals_2, 'o-', 'Color', colors{2},  'MarkerSize', 12, 'MarkerFaceColor', colors{2}, 'MarkerEdgeColor', markerC,'LineWidth',1.75)
        hold on

        % Plot SEM per point
        p1 = polarplot([deg2rad(anglevals); deg2rad(anglevals)], [vals_1 - sem1; vals_1 + sem1], '-', 'Color', colors{1}, 'LineWidth',1.75);
        hold on
        p2 = polarplot([deg2rad(anglevals); deg2rad(anglevals)], [vals_2 - sem2; vals_2 + sem2], '-', 'Color', colors{2}, 'LineWidth',1.75);
        hold on

       for subjectIndex = 1:size(medianBOLDpa, 4)
            %p3  = polarplot(deg2rad(anglevals), avgConditions1(:, subjectIndex)', 'o', 'Color', [127/255, 191/255, 123/255]);
            %hold on;
            %p4 = polarplot(deg2rad(anglevals), avgConditions2(:, subjectIndex)', 'o', 'Color', [175/255, 141/255, 195/255]);
            %hold on
        %     hold on
            asymm = avgConditions1(:, subjectIndex) - avgConditions2(:, subjectIndex);
            if (sum(asymm<0)) ~= 0
                sprintf('WaRNING: %s points are below 0', num2str(sum(asymm<0)))
            end
            
        end
        
        thetaticks(0:45:315);
    
    
        if ri==1
            if radialvstang==1 && strcmp(projectName, 'da')
                hLegend = legend('Radial', 'Tangential', 'Location', 'northwest', 'Box', 'off', 'FontSize', 18);
            else
                hLegend = legend('Cardinal', 'Oblique', 'Location', 'northwest', 'Box', 'off', 'FontSize', 18);
            end
            % Adjust the position of the legend
            newPosition = get(hLegend, 'Position'); % Get current position
            newPosition(1) = newPosition(1) + 0.1; % Shift the legend to the right by 0.1 normalized units
            newPosition = [0.7858 0.8972 0.1929 0.1131];
            nlegend = set(hLegend, 'Position', newPosition);
        end
        
        ax = gca;
            
        if ri == 5 || ri == 6 || ri == 7 %|| ri==1 || ri==2
            ax.RLim = [0 2]; %[0 2];%[0 2]; %[1 3]; %
            ax.RLim = [axes_limits.(projectName).(comparisonName).ROIs_motion.min ...
                axes_limits.(projectName).(comparisonName).ROIs_motion.max];
    %     elseif ri==3  % remove later (only for orientation
    %         ax.RLim = [-0.25 2];
        else
            ax = gca;
            ax.RLim = [-.25 1.75]; %[-1 1]; %[-.25 1.75]; %[0 2]; %
            ax.RLim = [axes_limits.(projectName).(comparisonName).ROIs_early.min ...
                axes_limits.(projectName).(comparisonName).ROIs_early.max];
        end
    
    
        ax = gca;
        ax.LineWidth = 3;  % Set the line width (adjust as needed)
        ax.GridColor = [0.25 0.25 0.25];
        ax.ThetaTickLabel = {};
        ax.Box = 0;
        %ax.RTickLabel = [];
        title(rois{ri}, 'FontSize', 18)
    end
    
    fig1 = gcf;
    sgtitle(sprintf('%s: %s', projectName, strrep(comparisonName, "_", " ")), 'FontSize', 40)
    fig1.Position = [34 228 1210 924]; %[152 569 2143 619];
    fig1.Color = 'w';
    hold off;
    
    filename = fullfile(figureDir,sprintf('polarangle_%s_%s_%s', comparisonName, projectName, asymmetryName));
 
    gcf_edit = fitFig2Page(gcf);
    
    % Save as PDF
    print(gcf_edit, filename, '-dpdf');
    close all;

end