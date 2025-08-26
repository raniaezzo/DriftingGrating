function plotSepDirs(medianBOLDpa, asymmetryName, grouping, projectSettings, varargin)


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

    % re-assign asymmetry name within this plotting function
    if strcmp(grouping, 'Abs')
        if strcmp(comparisonName, 'orientation_minus_baseline')
            asymmetryName = 'eachAbsorientation';
            legendNames = {'horizontal', 'upperright', 'vertical', 'upperleft'};
        elseif strcmp(comparisonName, 'motion_minus_baseline')
            asymmetryName = 'eachAbsmotiondir';
            legendNames = {'rightwards', 'upperrightwards', 'upwards', 'upperleftwards', ...
                'leftwards', 'lowerleftwards', 'downwards', 'lowerrightwards'};
        elseif strcmp(comparisonName, 'motion_minus_orientation')
            asymmetryName = 'eachAbsmotiondir';
            legendNames = {'rightwards', 'upperrightwards', 'upwards', 'upperleftwards', ...
                'leftwards', 'lowerleftwards', 'downwards', 'lowerrightwards'};
        end
    elseif strcmp(grouping, 'Rel')
       if strcmp(comparisonName, 'orientation_minus_baseline')
            asymmetryName = 'eachRelorientation';
            legendNames = {'tangential', 'radial', 'cclock', 'clock'};
        elseif strcmp(comparisonName, 'motion_minus_baseline')
            asymmetryName = 'eachRelmotiondir';
            legendNames = {'inwards', 'outwards', 'tangclock', 'tangcclock', ...
                'outclock', 'outcclock', 'inclock', 'incclock'};
        elseif strcmp(comparisonName, 'motion_minus_orientation')
            asymmetryName = 'eachRelmotiondir';
            legendNames = {'inwards', 'outwards', 'tangclock', 'tangcclock', ...
                'outclock', 'outcclock', 'inclock', 'incclock'};
        end
    end

    colors = cellfun(@(name) colors_data.conditions.(projectName).(asymmetryName).(['color_' name]).', legendNames, 'UniformOutput', false);
            
    markerC = colors_data.conditions.(projectName).plotSettings.markerC;
    
    % retrieve the indices for specific asymmetries (e.g., motion -
    % orientation for main cardinal v main oblique)
    if ~derivedVals
        if strcmp(grouping, 'Abs')
            allConditions = retrieveAbsDirIdx(projectName, comparisonName);
        elseif strcmp(grouping, 'Rel')
            allConditions = retrieveRelDirIdx(projectName, comparisonName);
        end
    else
        nConds = size(medianBOLDpa,1);
        allConditions = 1:nConds;
    end
    
    figure
    gap = [.04 .01]; % spacing between the subplots vertical gap - horizontal gap
    marg_h = [.015 .13]; % margins of bottom - top of the figure
    marg_w = [.01 .01]; % margina - left, right
    [ha, pos] = tight_subplot(2, 1, gap, marg_h, marg_w);
        

    % if motion, average the to colinear directions
    if strcmp(comparisonName, 'motion_minus_orientation') || strcmp(comparisonName, 'motion_minus_baseline')
        
        if strcmp(grouping, 'Abs')
            pairedConds = {["rightwards", "leftwards"], ["upperrightwards", "lowerleftwards"], ["upwards", "downwards"], ["upperleftwards", "lowerrightwards"]};
            % colors stays the same its colors{1:4}, ignore 5:8
        elseif strcmp(grouping, 'Rel')
            pairedConds = {["inwards", "outwards"], ["tangclock", "tangcclock"], ["outclock", "incclock"], ["outcclock", "inclock"]};
            % colors are ordered interleaved, so extract
            colors = colors(1:2:end);
        end
        
        pairedConditionIndices = cell(size(pairedConds));
        for i = 1:length(pairedConds)
            pair = pairedConds{i};
            idx1 = find(strcmp(legendNames, pair(1)));
            idx2 = find(strcmp(legendNames, pair(2)));
            pairedConditionIndices{i} = [allConditions(idx1), allConditions(idx2)];
        end
        allConditions = pairedConditionIndices;
        allConditions = cell2mat(allConditions')';

        if strcmp(grouping, 'Abs')
            legendNames = {"right-leftwards", "upperright-lowerleftwards", "up-downwards", "upperleft-lowerrightwards"};
        elseif strcmp(grouping, 'Rel')
            legendNames = {"in-out", "clock-cclock", "outclock-incclock", "outcclock-inclock"};
        end
    end

    for ri=1:length(rois)
       % Specify the region index
       regionIndex = projectSettings.roi_idx{ri};

       % using this method instead of subplot for TIGHT AXES
       axes(ha(ri)); 

        % loop across the directions/orientations
        for cc = 1:size(allConditions,2)

            currCond = allConditions(:,cc)';

            % Extract the relevant conditions for the specified region
            conditions1 = medianBOLDpa(currCond, :, regionIndex, :);
            conditions1 = squeeze(conditions1);
            
            if strcmp(comparisonName, 'motion_minus_orientation') || strcmp(comparisonName, 'motion_minus_baseline')
                avgConditions1 = squeeze(mean(conditions1,1)); % average colinear directions
            else
                avgConditions1 = conditions1;  % for orientation DA - there is only 1 value
            end

            % Extract polar angles
            %anglevals = [90, 135, 180, 225, 270, 315, 0, 45];
            anglevals = [90, 45, 0, 315, 270, 225, 180, 135]; % <-- these were manually converted based on the order of polarAngles above (Noah's convention)
            
            vals_1 = nanmean(avgConditions1,2)';
            sem1 = nanstd(avgConditions1,0,2)' ./ sqrt(sum(~isnan(avgConditions1),2)');
        
            
            % plot average (connecting the last line)
            polarplot([deg2rad(anglevals(end)), deg2rad(anglevals(1))],[vals_1(end), vals_1(1)], 'o-', 'Color', colors{cc}, 'MarkerSize', 12, 'LineWidth',1.75,  'MarkerFaceColor', colors{cc}, 'MarkerEdgeColor', markerC)
            hold on
            
            % plot average (all the data)
            polarplot(deg2rad(anglevals),vals_1, 'o-', 'Color', colors{cc},  'MarkerSize', 12, 'MarkerFaceColor', colors{cc}, 'MarkerEdgeColor', markerC,'LineWidth',1.75)
            hold on
        
            % Plot SEM per point
            p1 = polarplot([deg2rad(anglevals); deg2rad(anglevals)], [vals_1 - sem1; vals_1 + sem1], '-', 'Color', colors{cc}, 'LineWidth',1.75);
            hold on
            
            thetaticks(0:45:315);
        end
    
    
        % Create invisible scatter plots just for the legend
        hold on;
        legendHandles = gobjects(1, length(colors)); % Preallocate handles
        for curCol = 1:length(colors)
            legendHandles(curCol) = polarplot([100, 100], [100, 100], 'o', 'Color', colors{curCol}, 'MarkerFaceColor', colors{curCol}, 'MarkerSize', 8);
        end

        if ri==1
            hLegend = legend(legendHandles, legendNames, 'Location', 'northwest', 'Box', 'off', 'FontSize', 15);
            % Adjust the position of the legend
            newPosition = get(hLegend, 'Position'); % Get current position
            newPosition(1) = newPosition(1) + 0.1; % Shift the legend to the right by 0.1 normalized units
            newPosition = [0.7858 0.7972 0.1929 0.1131];
            nlegend = set(hLegend, 'Position', newPosition);
        end
        
        ax = gca;
            
        if strcmp(rois{ri}, 'hMTcomplex') || strcmp(rois{ri}, 'pMT') || strcmp(rois{ri}, 'pMST')
            ax.RLim = [0 2]; %[0 2];%[0 2]; %[1 3]; %
            ax.RLim = [axes_limits.(projectName).(comparisonName).ROIs_motion.min ...
                axes_limits.(projectName).(comparisonName).ROIs_motion.max];
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