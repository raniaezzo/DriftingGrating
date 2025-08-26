function plot2_experimentalCond(medianBOLDpa, asymmetryName, projectSettings, varargin)

    rng(0)
    % Plot mean across polar angles
    % keep in mind that this equally weighs each PA, whereas there could be
    % differential # of voxels representing the PAs

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
    comparisonName = projectSettings.comparisonName;
    colors_data = projectSettings.colors_data;
    rois = projectSettings.rois;
    %pairaxes_limits = projectSettings.pairaxes_limits;
    pairaxes_PAew_limits = projectSettings.pairaxes_PAew_limits;
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
    
    for ii=1:length(rois)
        rois{ii}
    
        % Specify the region index
        regionIndex = projectSettings.roi_idx{ii};
    
        % Extract the relevant conditions for the specified region
        conditions1 = medianBOLDpa(proConditions, :, regionIndex, :);
        conditions2 = medianBOLDpa(conConditions, :, regionIndex, :);
        conditions1 = squeeze(conditions1);
        conditions2 = squeeze(conditions2);
        
        % Average across the conditions within subjects
        avgConditions1 = nanmean(conditions1, 1);
        avgConditions2 = nanmean(conditions2, 1);
        avgConditions1 = squeeze(avgConditions1);
        avgConditions2 = squeeze(avgConditions2);
        
        % Already averaged across the conditions within subjects < -- already did this in
        % the loop above
        
        % Extract polar angles
        %anglevals = [90, 135, 180, 225, 270, 315, 0, 45];
        anglevals = [90, 45, 0, 315, 270, 225, 180, 135]; % <-- these were manually converted based on the order of polarAngles above (Noah's convention)
        
        vals_1 = nanmean(avgConditions1,1)';
        vals_2 = nanmean(avgConditions2,1)';
    
        vals_1_overall = nanmean(avgConditions1,'all')';
        vals_2_overall = nanmean(avgConditions2,'all')';
        %sem1 = nanstd(avgConditions1,0,2)' ./ sqrt(sum(~isnan(avgConditions1),2)');
        %sem2 = nanstd(avgConditions2,0,2)' ./ sqrt(sum(~isnan(avgConditions2),2)');
        
        % Plot the data on a polar plot
        subplot(2,5,ii)
        
        for subjectIndex = 1:size(medianBOLDpa, 4)
            %scatter(1, vals_1(subjectIndex),  30, 'MarkerFaceColor', colors_rgb(subjectIndex,:), 'MarkerEdgeColor', 'none'); % [127/255, 191/255, 123/255]
            scatter(1, vals_1(subjectIndex),  30, 'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'none'); % 
            hold on
            %scatter(2, vals_2(subjectIndex), 30, 'MarkerFaceColor', colors_rgb(subjectIndex,:), 'MarkerEdgeColor', 'none'); %  [175/255, 141/255, 195/255]
            scatter(2, vals_2(subjectIndex), 30, 'MarkerFaceColor', colors{2}, 'MarkerEdgeColor', 'none'); 
            plot([1 2], [vals_1(subjectIndex) vals_2(subjectIndex)], 'Color', 'k'); %colors_rgb(subjectIndex,:))
            xlim([0 3])
    %         ylim([-0.15 0.25])
        end
        
        %errDiff = std(vals_1-vals_2)/(sqrt((size(medianBOLDpa, 4))));
        %errorbar(1.5, mean(vals_1-vals_2), errDiff, 'k', 'LineWidth', 3);
        %[h, p] = ttest(vals_1-vals_2, 0)

        fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        % Assume: vals_1 and vals_2 are column vectors [nSubjects x 1]
        differences = vals_1 - vals_2;
        nBoot = 10000;
        nSubs = length(differences);
        
        %% Bootstrap Mean
        bootMeans = zeros(nBoot, 1);
        for i = 1:nBoot
            resample = datasample(differences, nSubs);
            bootMeans(i) = mean(resample);
        end
        meanDiff = mean(differences);
        ci_mean = prctile(bootMeans, [2.5 97.5]);
        errLower_mean = meanDiff - ci_mean(1);
        errUpper_mean = ci_mean(2) - meanDiff;
        
        %% Bootstrap Median
        bootMedians = zeros(nBoot, 1);
        for i = 1:nBoot
            resample = datasample(differences, nSubs);
            bootMedians(i) = median(resample);
        end
        medianDiff = median(differences);
        ci_median = prctile(bootMedians, [2.5 97.5]);
        errLower_median = medianDiff - ci_median(1);
        errUpper_median = ci_median(2) - medianDiff;
        
        %% Print to console
        fprintf('Mean difference: %.4f, 95%% CI: [%.4f, %.4f]\n', meanDiff, ci_mean(1), ci_mean(2));
        fprintf('Median difference: %.4f, 95%% CI: [%.4f, %.4f]\n', medianDiff, ci_median(1), ci_median(2));
        differences
        fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        

        plot([1 2], [vals_1_overall vals_2_overall], 'k', 'LineWidth', 3)
        hold on
        scatter(1, vals_1_overall, 70, 'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', markerC, 'LineWidth',1.75); %, 'MarkerFaceAlpha', 0.5);
        hold on
        scatter(2, vals_2_overall,  70, 'MarkerFaceColor', colors{2}, 'MarkerEdgeColor', markerC, 'LineWidth',1.75); %, 'MarkerFaceAlpha', 0.5);
        
        title(rois{ii});
        ylabel('zscored PSC')
        set(gca, 'XTick', []);

        if ismember(rois{ii}, {'pMT', 'pMST', 'hMTcomplex'})
            ROI_category = 'ROIs_motion';
        else
            ROI_category = 'ROIs_early';
        end

        ylim([pairaxes_PAew_limits.(projectName).(comparisonName).(ROI_category).min ...
                     pairaxes_PAew_limits.(projectName).(comparisonName).(ROI_category).max])
    
        if ii==1
            if radialvstang
                lg1 = legend('Radial', 'Tangential', 'Location', 'northeast');
            else
                lg1 = legend('Card', 'Obl', 'Location', 'northeast');
            end
            lg1.Position = [0.7946 0.2614 0.1089 0.0631];
        end
    
        hold off;
    end
    
    sgtitle(sprintf('%s: %s', projectName, strrep(comparisonName, '_', ' ')), 'FontSize', 40)

%     fa = gcf;
%     fa.Position = [1000 555 1514 782];

    filename = fullfile(figureDir,sprintf('pairwise_PAequalweight_%s_%s_%s', comparisonName, projectName, asymmetryName));


    gcf_edit = fitFig2Page(gcf);
    
    % Save as PDF
    print(gcf_edit, filename, '-dpdf');
    close all;

end
