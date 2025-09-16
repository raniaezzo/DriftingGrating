function plot0_experimentalCond(condIdx1, condIdx2, medianBOLD, projectSettings)

    projectName = projectSettings.projectName;
    comparisonName = projectSettings.comparisonName;
    subjects = projectSettings.subjects;
    rois = projectSettings.rois;
    figureDir = projectSettings.figureDir;
    contrasts_dict = projectSettings.contrasts_dict;
    pairaxes_limits = projectSettings.pairaxes_limits;

    if strcmp(projectName, 'dots')
        % if dots, just use da, it's the same strings
        contrastnames = {contrasts_dict.contrasts.('dg_contrast_name')};
    else
        contrastnames = {contrasts_dict.contrasts.(strcat(projectName, '_contrast_name'))};
    end

    % compute number of subjects
    nSubjects = size(medianBOLD,3);

    % Color values for each subject
    hue_values = linspace(0, 1, nSubjects+1);
    hue_values = hue_values(1:end-1);  % Exclude the last point to avoid duplication
    saturation = 1;     % Set saturation and value for all colors
    value = .8; %1; % lower brightness
    % Convert HSV to RGB
    colors_hsv = [hue_values' repmat(saturation, nSubjects, 1) repmat(value, nSubjects, 1)];
    colors_rgb = hsv2rgb(colors_hsv);
    colors = colors_rgb; %rand(length(subjects), 3);

    scatterSize = 55;

    figure
    
    % Select relevant rows and calculate the ratio
    condition1 = squeeze(nanmean(medianBOLD(condIdx1, :, :),1));  % First condition
    condition2 = squeeze(nanmean(medianBOLD(condIdx2, :, :),1));  % Second condition
    
    % initalize new mats
    meanMAT = nan(numel(rois), 2);
    semMAT = nan(numel(rois), 2);
    diffMAT = nan(numel(rois), 1);
    diffSEM = nan(numel(rois), 1);
    
    % fill the values for the initialized mats
    for region=1:numel(rois)
        meanMAT(region,1) = nanmean(condition1(region,:));
        semMAT(region,1) = nanstd(condition1(region, :)) / sqrt(sum(~isnan(condition1(region, :))));
        meanMAT(region,2) = nanmean(condition2(region,:));
        semMAT(region,2) = nanstd(condition2(region, :)) / sqrt(sum(~isnan(condition2(region, :))));
        diffMAT(region) = nanmean(condition1(region,:) - condition2(region,:));
        diffSEM(region) = nanstd(condition1(region,:) - condition2(region,:)) / sqrt(sum(~isnan(condition1(region,:) - condition2(region,:))));
    end

    for region = 1:size(meanMAT, 1)
        subplot(2, 5, region);  % Adjust subplot layout as needed
    
        % Plot points with different colors for conditions
        %scatter(1:2, meanMAT(region, :), 'o', 'DisplayName', 'Mean', 'LineWidth', 1.5);
        %hold on;
        scatter(ones(1, length(condition1(region,:))), condition1(region,:), scatterSize, colors, 'filled','MarkerEdgeColor', 'white','LineWidth', 1.5);
        hold on
        scatter(2*ones(1, length(condition2(region,:))), condition2(region,:), scatterSize, colors, 'filled','MarkerEdgeColor', 'white','LineWidth', 1.5);
        hold on
    
        for si=1:nSubjects
            plot([1,2], [condition1(region,si), condition2(region,si)], 'color', [colors(si, :), 0.25], 'LineWidth',1)
            hold on
        end
    
        % Error bars with different colors for conditions
        errorbar(1:2, meanMAT(region, :), semMAT(region, :), '.', 'Color', 'k', 'DisplayName', 'SEM', 'LineWidth',2);
        hold on
    
        errorbar(1.5, diffMAT(region), semMAT(region), '.', 'Color', [.5 .5 .5], 'DisplayName', 'SEM', 'LineWidth',2);
    
        xlim([0.5,2.5])

        if ismember(rois(region), {'pMT', 'pMST', 'hMTcomplex'})
            ROI_category = 'ROIs_motion';
        else
            ROI_category = 'ROIs_early';
        end

        ylim([pairaxes_limits.(projectName).(comparisonName).(ROI_category).min ...
                     pairaxes_limits.(projectName).(comparisonName).(ROI_category).max])


        title([rois(region)]);
        xlabel('Contrast');
        ylabel('zscored PSC');
        f1 = gcf;
        f1.Position = [1 86 1508 971]; %[-3 137 1508 660];
        ax = gca;
        ax.XTick = [1,2];
        condName1 = strjoin(arrayfun(@(idx) strrep(contrastnames{idx}, 'sep', ''), condIdx1, 'UniformOutput', false), '_');
        %strrep(contrastnames{condIdx1}, 'sep', '');
        condName2 = strjoin(arrayfun(@(idx) strrep(contrastnames{idx}, 'sep', ''), condIdx2, 'UniformOutput', false), '_');
        %strrep(contrastnames{condIdx2}, 'sep', '');
        ax.XTickLabel = {condName1,condName2}; %{contrastnames{1},'',contrastnames{2},''};

        % Plot invisible scatter points to use in the legend
        if region == numel(rois)
            scatterHandles = gobjects(1, nSubjects);
            for i = 1:nSubjects
                scatterHandles(i) = scatter(ones(1,1), ones(1,1), 1, colors(i, :), 'filled');
            end

            leg = legend(scatterHandles, subjects, 'Location', 'eastoutside');
            leg.Position = [0.7988, 0.2443, 0.0666, 0.1273];
        end

        hold off;
    end
    
    sgtitle(sprintf('%s: zscored PSC: %s and %s', strrep(comparisonName, '_', ' '), condName1,condName2));

    filename = fullfile(figureDir,sprintf('pairwise_%s_%s_%s_v_%s', comparisonName, projectName, condName1,condName2));


    gcf_edit = fitFig2Page(gcf);
    
    % Save as PDF
    print(gcf_edit, filename, '-dpdf');
    close all;
    
end