%%

% roi=1;
% for si=1:numel(subjects)
%     subjects(si)
%     rois(roi)
%     cardinaldir_cardloc = medianBOLDpa(1,[1,3,5,7],roi,si);
%     obliquedir_cardloc = medianBOLDpa(2,[1,3,5,7],roi,si);
%     cardinaldir_obliqloc = medianBOLDpa(1,[2,4,6,8],roi,si);
%     obliquedir_obliqloc = medianBOLDpa(2,[2,4,6,8],roi,si);
%     
%     disp('Cardinal advantage at cardinal locs:')
%     %cardinaldir_cardloc-obliquedir_cardloc
%     nanmean(cardinaldir_cardloc-obliquedir_cardloc)
%     disp('Cardinal advantage at oblique locs:')
%     %cardinaldir_obliqloc-obliquedir_obliqloc
%     nanmean(cardinaldir_obliqloc-obliquedir_obliqloc)
% 
% %     disp('Cardinal effect at cardinal locs vs oblique locs:')
% %     cardinaldir_cardloc-cardinaldir_obliqloc
% 
% %     disp('Oblique effect at Oblique locs vs Cardinal locs:')
% %     obliquedir_obliqloc-obliquedir_cardloc
% end




%% Polar plot (POLAR CARDINAL vs POLAR OBLIQUE)

% if strcmp(projectName, 'da')
%     markerC = 'k';
%     colors = {colors_all{1}, colors_all{2}};
% elseif strcmp(projectName, 'dg')
%     markerC = 'w';
%     colors = {colors_all{3}, colors_all{4}};
% end
% 
% figure
% for ri=1:7 %length(rois)
% 
%     % Specify the region index
%     regionIndex = ri;
%     
%     % Average the values in the first two rows along the first dimension
%     averagedMatrix = mean(newMatrix(1:2, :, :, :), 1); % radial and tangential
%     
%     % Combine the averaged values with the third row
%     polarcardMatrix = cat(1, averagedMatrix, newMatrix(3, :, :, :));
%     
%     % Extract the relevant conditions for the specified region
%     conditions1 = polarcardMatrix(1, :, regionIndex, :);
%     conditions2 = polarcardMatrix(2, :, regionIndex, :);
%     avgConditions1 = squeeze(conditions1);
%     avgConditions2 = squeeze(conditions2);
%     
% %     min_axis1 = min(min(squeeze(min(avgConditions1))));
% %     min_axis2 = min(min(squeeze(min(avgConditions2))));
% %     axMIN = customRound(min([min_axis1, min_axis2]));
% %     max_axis1 = max(max(squeeze(max(avgConditions1))));
% %     max_axis2 = max(max(squeeze(max(avgConditions2))));
% %     axMAX = customRound(max([max_axis1, max_axis2]));
% %     axMIN = axMIN - abs(axMAX-axMIN).*0.05;
% %     axMAX = axMAX + abs(axMAX-axMIN).*0.05;
%     
%     % Already averaged across the conditions within subjects < -- already did this in
%     % the loop above
%     
%     % Extract polar angles
%     %anglevals = [90, 135, 180, 225, 270, 315, 0, 45];
%     anglevals = [90, 45, 0, 315, 270, 225, 180, 135]; % <-- these were manually converted based on the order of polarAngles above (Noah's convention)
%     
%     vals_1 = nanmean(avgConditions1,2)';
%     sem1 = nanstd(avgConditions1,0,2)' ./ sqrt(sum(~isnan(avgConditions1),2)');
%     vals_2 = nanmean(avgConditions2,2)';
%     sem2 = nanstd(avgConditions2,0,2)' ./ sqrt(sum(~isnan(avgConditions2),2)');
%     
%     % Plot the data on a polar plot - 
%     subplot(2,4,ri)
%     
%     for subjectIndex = 1:size(medianBOLDpa, 4)
%         polarplot([deg2rad(anglevals(end)), deg2rad(anglevals(1))],[vals_1(end), vals_1(1)], 'o-', 'Color', colors{1}, 'LineWidth',1.75,  'MarkerFaceColor', colors{1})
%         hold on
%         polarplot([deg2rad(anglevals(end)), deg2rad(anglevals(1))],[vals_2(end), vals_2(1)], 'o-', 'Color', colors{2}, 'LineWidth',1.75,  'MarkerFaceColor', colors{2})
%         hold on
%         p1 = polarplot([deg2rad(anglevals); deg2rad(anglevals)], [vals_1 - sem1; vals_1 + sem1], '-', 'Color', colors{1}, 'LineWidth',1.75);
%         hold on
%         p2 = polarplot([deg2rad(anglevals); deg2rad(anglevals)], [vals_2 - sem2; vals_2 + sem2], '-', 'Color', colors{2}, 'LineWidth',1.75);
%         hold on
% %         p3  = polarplot(deg2rad(anglevals), avgConditions1(:, subjectIndex)', 'o', 'Color', [166 97 26]/255);
% %         hold on;
% %         p4 = polarplot(deg2rad(anglevals), avgConditions2(:, subjectIndex)', 'o', 'Color', [64 176 166]/255);
% %         hold on
%         polarplot(deg2rad(anglevals),vals_1, 'o-', 'Color', colors{1}, 'MarkerSize', 12,  'LineWidth',1.75,  'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', markerC)
%         hold on
%         polarplot(deg2rad(anglevals),vals_2, 'o-', 'Color', colors{2}, 'MarkerSize', 12,  'LineWidth',1.75,  'MarkerFaceColor', colors{2}, 'MarkerEdgeColor', markerC)
%         hold on
%     end
%     
%     rlim([axMIN, axMAX])
%     thetaticks(0:45:315);
%     
%     if ri==1
%         hLegend = legend('Cardinal', 'Oblique', 'Location', 'northeast', 'Box', 'off', 'FontSize', 18);
%         % Adjust the position of the legend
%         newPosition = get(hLegend, 'Position'); % Get current position
%         newPosition(1) = newPosition(1) + 0.1; % Shift the legend to the right by 0.1 normalized units
%         set(hLegend, 'Position', newPosition);
%     end
%     
%     ax = gca;
%         
%     if ri == 5 || ri == 6 || ri == 7 || ri == 1 || ri == 2
%         %ax.RLim = [0 2]; %[0 2]; %[0 2]; %[1 3]; %
%         ax.RLim = [axes_limits.(projectName).(comparisonName).ROIs_motion.min ...
%             axes_limits.(projectName).(comparisonName).ROIs_motion.max];
% %     elseif ri == 3
% %         ax.RLim = [-.25 2];
%     else
%         ax = gca;
%         %ax.RLim = [-.25 1.75]; %[-1 1];   %[-.25 1.75]; %[0 2]; %
%         ax.RLim = [axes_limits.(projectName).(comparisonName).ROIs_early.min ...
%             axes_limits.(projectName).(comparisonName).ROIs_early.max];
%     end
% 
% %     if ri == 5 || ri == 6 || ri == 7
% %         ax.RLim = [0 2];
% %         temp = linspace(0, 2, 5);
% %         ax.RTick =  temp(1:end-1);
% %         ax = gca;
% %         ax.RLim = [0 2];
% %     else
% %         ax.RLim = [-0.75 1.25];  % Set the radial limits (adjust as needed)
% %     %     ax.RTick =  [-0.0500, 0, 0.0500]; %lineArray(1, 3, 5);
% %         temp = linspace(-0.75, 1.25, 5);
% %         ax.RTick =  temp(1:end-1);
% %         ax = gca;
% %         ax.RLim = [-0.75 1.25];
% %     end
%     
%     ax = gca;
%     ax.LineWidth = 3;  % Set the line width (adjust as needed)
%     ax.GridColor = [0.25 0.25 0.25];
%     ax.ThetaTickLabel = {};
%     %lineArray = -0.05:0.025:.075; 
%     ax.Box = 0;
%     %ax.RTickLabel = [];
%     title(rois{ri}, 'FontSize', 18)
% end
% 
% fig1 = gcf;
% fig1.Position = [152 810 477 378]; %[152 247 1702 941];
% 
% fig1 = gcf;
% fig1.Position = [152 569 2143 619];
% hold off;
% sgtitle(projectName, 'FontSize', 40)
% 
% % fig1 = gcf;
% % fig1.Position = [152 247 1702 941];
% % hold off;
% 
% filename = fullfile(figureDir,sprintf('PApolarcardvpolaroblmotionVSstatic%s', rois{regionIndex}));
% saveas(gcf, filename, 'pdf');
% 
% 
% %% POLAR CARDINAL vs POLAR OBLIQUE - Combine Polar Angle
% % Plot mean across polar angles 
% % keep in mind that this equally weighs each PA, whereas there could be
% % differential # of voxels representing the PAs
% 
% figure
% 
% for ii=1:length(rois)
% 
%     % Specify the region index
%     regionIndex = ii;
% 
%     % Average the values in the first two rows along the first dimension
%     averagedMatrix = mean(newMatrix(1:2, :, :, :), 1); % radial and tangential
%     
%     % Combine the averaged values with the third row
%     polarcardMatrix = cat(1, averagedMatrix, newMatrix(3, :, :, :));
%     
%     % Extract the relevant conditions for the specified region
%     conditions1 = polarcardMatrix(1, :, regionIndex, :);
%     conditions2 = polarcardMatrix(2, :, regionIndex, :);
%     avgConditions1 = squeeze(conditions1);
%     avgConditions2 = squeeze(conditions2);
%     
%     min_axis1 = min(min(squeeze(min(avgConditions1))));
%     min_axis2 = min(min(squeeze(min(avgConditions2))));
%     axMIN = customRound(min([min_axis1, min_axis2]));
%     max_axis1 = max(max(squeeze(max(avgConditions1))));
%     max_axis2 = max(max(squeeze(max(avgConditions2))));
%     axMAX = customRound(max([max_axis1, max_axis2]));
%     axMIN = axMIN - abs(axMAX-axMIN).*0.05;
%     axMAX = axMAX + abs(axMAX-axMIN).*0.05;
%     
%     % Already averaged across the conditions within subjects < -- already did this in
%     % the loop above
%     
%     % Extract polar angles
%     %anglevals = [90, 135, 180, 225, 270, 315, 0, 45]; % just commented
%     %these out 11/19/2024
%     anglevals = [90, 45, 0, 315, 270, 225, 180, 135]; % <-- these were manually converted based on the order of polarAngles above (Noah's convention)
%     
%     vals_1 = nanmean(avgConditions1,1)';
%     vals_2 = nanmean(avgConditions2,1)';
% 
%     vals_1_overall = nanmean(avgConditions1,'all')';
%     vals_2_overall = nanmean(avgConditions2,'all')';
%     %sem1 = nanstd(avgConditions1,0,2)' ./ sqrt(sum(~isnan(avgConditions1),2)');
%     %sem2 = nanstd(avgConditions2,0,2)' ./ sqrt(sum(~isnan(avgConditions2),2)');
%     
%     % Plot the data on a polar plot
%     subplot(2,5,ii)
%     
%     for subjectIndex = 1:size(medianBOLDpa, 4)
%         scatter(1, vals_1(subjectIndex),  30, 'MarkerFaceColor', colors_rgb(subjectIndex,:), 'MarkerEdgeColor', 'none'); %, [139/255, 69/255, 19/255]
%         hold on
%         scatter(2, vals_2(subjectIndex), 30, 'MarkerFaceColor', colors_rgb(subjectIndex,:), 'MarkerEdgeColor', 'none'); %, [0/255, 139/255, 139/255]
%         plot([1 2], [vals_1(subjectIndex) vals_2(subjectIndex)], 'Color', colors_rgb(subjectIndex,:));  %; %'k')
%         xlim([0 3])
%         %ylim([-0.15 0.25])
%     end
%     
%     scatter(1, vals_1_overall, 70, 'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'black');
%     hold on
%     scatter(2, vals_2_overall,  70, 'MarkerFaceColor', colors{2}, 'MarkerEdgeColor', 'black');
%     plot([1 2], [vals_1_overall vals_2_overall], 'k', 'LineWidth', 3)
%     
%     title(sprintf('BOLD per polar angle in %s',rois{regionIndex}));
%     ylabel('PSC')
%     set(gca, 'XTick', []);
% 
% %     if ii==1
% %         legend('PolarCard', 'PolarObl', 'Location', 'northeast');
% %     end
% 
%     hold off;
% end
% 
% sgtitle(projectName, 'FontSize', 40)
% fa = gcf;
% fa.Position = [1000 555 1514 782];
% 
% %% Polar plot (RADIAL / TANGENTIAL)
% 
% if strcmp(projectName, 'dg')
%     % Plot the data on a polar plot
%     figure;
%     
%     for ri=1:length(rois)
%     
%         % Specify the region index
%         regionIndex = ri;
%         
%         % Extract the relevant conditions for the specified region
%         conditions1 = newMatrix(1, :, regionIndex, :);
%         conditions2 = newMatrix(2, :, regionIndex, :);
%         avgConditions1 = squeeze(conditions1);
%         avgConditions2 = squeeze(conditions2);
%         
%         min_axis1 = min(min(squeeze(min(avgConditions1))));
%         min_axis2 = min(min(squeeze(min(avgConditions2))));
%         axMIN = customRound(min([min_axis1, min_axis2]));
%         max_axis1 = max(max(squeeze(max(avgConditions1))));
%         max_axis2 = max(max(squeeze(max(avgConditions2))));
%         axMAX = customRound(max([max_axis1, max_axis2]));
%         axMIN = axMIN - abs(axMAX-axMIN).*0.05;
%         axMAX = axMAX + abs(axMAX-axMIN).*0.05;
%         
%         % Already averaged across the conditions within subjects < -- already did this in
%         % the loop above
%         
%         % Extract polar angles
%         %anglevals = [90, 135, 180, 225, 270, 315, 0, 45];
%         anglevals = [90, 45, 0, 315, 270, 225, 180, 135]; % <-- these were manually converted based on the order of polarAngles above (Noah's convention)
%         
%         vals_1 = nanmean(avgConditions1,2)';
%         %vals_1 = avgConditions1(:,1);  for only subj1
%         sem1 = nanstd(avgConditions1,0,2)' ./ sqrt(sum(~isnan(avgConditions1),2)');
%         vals_2 = nanmean(avgConditions2,2)';
%         %vals_2 = avgConditions2(:,1);  for only subj1
%         sem2 = nanstd(avgConditions2,0,2)' ./ sqrt(sum(~isnan(avgConditions2),2)');
%         
%     
%         subplot(2,4,ri)
%     %     for subjectIndex = 1:size(medianBOLDpa, 4)
%     %         p1  = polarplot(deg2rad(anglevals), avgConditions1(:, subjectIndex)', 'o', 'Color', [176/255, 224/255, 230/255]);
%     %         hold on;
%     %         p2 = polarplot(deg2rad(anglevals), avgConditions2(:,
%     %         subjectIndex)', 'o', 'Color', 'r');
%     %         hold on
%     %     
%     %     end
%         
%     %     hold on
%         polarplot([deg2rad(anglevals(end)), deg2rad(anglevals(1))],[vals_1(end), vals_1(1)], 'o-', 'Color', [146 197 222]/255, 'LineWidth',1.75,  'MarkerFaceColor', [146 197 222]/255)
%         hold on
%         polarplot([deg2rad(anglevals(end)), deg2rad(anglevals(1))],[vals_2(end), vals_2(1)], 'o-', 'Color', [202 0 32]/255, 'LineWidth',1.75,  'MarkerFaceColor', [202 0 32]/255)
%         hold on
%         p3 = polarplot([deg2rad(anglevals); deg2rad(anglevals)], [vals_1 - sem1; vals_1 + sem1], '-', 'Color', [146 197 222]/255, 'LineWidth',1.75);
%         hold on
%         p4 = polarplot([deg2rad(anglevals); deg2rad(anglevals)], [vals_2 - sem2; vals_2 + sem2], '-', 'Color', [202 0 32]/255, 'LineWidth',1.75);
%         hold on
%         polarplot(deg2rad(anglevals),vals_1, 'o-','MarkerSize', 12,  'Color', [146 197 222]/255, 'LineWidth',1.75,  'MarkerFaceColor', [146 197 222]/255, 'MarkerEdgeColor', 'w')
%         hold on
%         polarplot(deg2rad(anglevals),vals_2, 'o-','MarkerSize', 12,  'Color', [202 0 32]/255, 'LineWidth',1.75,  'MarkerFaceColor', [202 0 32]/255,  'MarkerEdgeColor', 'w')
%         hold on 
%         
%     %     rlim([axMIN, axMAX])
%         thetaticks(0:45:315);
%         
%     %     title(sprintf('BOLD per polar angle in %s (PSC)',rois{regionIndex}));
%     % 
%         if ri==1
%             hLegend = legend('Radial', 'Tangential', 'Location', 'northeast', 'Box', 'off', 'FontSize', 18);
%             % Adjust the position of the legend
%             newPosition = get(hLegend, 'Position'); % Get current position
%             newPosition(1) = newPosition(1) + 0.1; % Shift the legend to the right by 0.1 normalized units
%             set(hLegend, 'Position', newPosition);
%         end
%         
%         ax = gca;
%     
%         if ri == 5 || ri == 6 || ri == 7
%             ax.RLim = [-0.25 2]; %[0 2];%[-0.25 2]; %[1 3]; %
%         %elseif ri == 1
%         %    ax.RLim = [-.75 1.25]; % get rid of this later
%         elseif ri == 3
%             ax.RLim = [0 2.25]; % orientation only
%         else
%             ax = gca;
%             ax.RLim = [-.25 2]; %[-1 1]; %  [-.25 2]; %[0 2]; %
%         end
%     
%     %     if ri == 5 || ri == 6 || ri == 7
%     %         ax.RLim = [0 2];
%     %         temp = linspace(0, 2, 5);
%     %         ax.RTick =  temp(1:end-1);
%     %         ax = gca;
%     %         ax.RLim = [0 2];
%     %     else
%     %         ax.RLim = [-0.75 1.25];  % Set the radial limits (adjust as needed)
%     %     %     ax.RTick =  [-0.0500, 0, 0.0500]; %lineArray(1, 3, 5);
%     %         temp = linspace(-0.75, 1.25, 5);
%     %         ax.RTick =  temp(1:end-1);
%     %         ax = gca;
%     %         ax.RLim = [-0.75 1.25];
%     %     end
%         
%         ax = gca;
%         ax.LineWidth = 3;  % Set the line width (adjust as needed)
%         ax.GridColor = [0.25 0.25 0.25];
%         ax.ThetaTickLabel = {};
%         %lineArray = -0.05:0.025:.075; 
%         ax.Box = 0;
%         %ax.RTickLabel = [];
%         title(rois{ri}, 'FontSize', 18)
%     
%     
%     end
%     
%     fig1 = gcf;
%     fig1.Position = [152 810 477 378]; %[152 247 1702 941];
%     
%     fig1 = gcf;
%     fig1.Position = [152 569 2143 619];
%     hold off;
%     
%     
%     % fig1 = gcf;
%     % fig1.Position = [152 247 1702 941];
%     % hold off;
%     sgtitle(projectName, 'FontSize', 40)
%     
%     filename = fullfile(figureDir,sprintf('PAradialvtangmotionVSstatic%s', rois{regionIndex}));
%     saveas(gcf, filename, 'pdf');
% end
% 
% 
% %% RADIAL / TANGENTIAL - Combine Polar Angle
% % Plot mean across polar angles 
% % keep in mind that this equally weighs each PA, whereas there could be
% % differential # of voxels representing the PAs
% 
% if strcmp(projectName, 'dg')
%     figure
%     
%     for ii=1:length(rois)
%     
%         % Specify the region index
%         regionIndex = ii;
%     
%         % Extract the relevant conditions for the specified region
%         conditions1 = newMatrix(1, :, regionIndex, :);
%         conditions2 = newMatrix(2, :, regionIndex, :);
%         avgConditions1 = squeeze(conditions1);
%         avgConditions2 = squeeze(conditions2);
%         
%         min_axis1 = min(min(squeeze(min(avgConditions1))));
%         min_axis2 = min(min(squeeze(min(avgConditions2))));
%         axMIN = customRound(min([min_axis1, min_axis2]));
%         max_axis1 = max(max(squeeze(max(avgConditions1))));
%         max_axis2 = max(max(squeeze(max(avgConditions2))));
%         axMAX = customRound(max([max_axis1, max_axis2]));
%         axMIN = axMIN - abs(axMAX-axMIN).*0.05;
%         axMAX = axMAX + abs(axMAX-axMIN).*0.05;
%         
%         % Already averaged across the conditions within subjects < -- already did this in
%         % the loop above
%         
%         % Extract polar angles
%         %anglevals = [90, 135, 180, 225, 270, 315, 0, 45];
%         anglevals = [90, 45, 0, 315, 270, 225, 180, 135]; % <-- these were manually converted based on the order of polarAngles above (Noah's convention)
%         
%         vals_1 = nanmean(avgConditions1,1)';
%         vals_2 = nanmean(avgConditions2,1)';
%     
%         vals_1_overall = nanmean(avgConditions1,'all')';
%         vals_2_overall = nanmean(avgConditions2,'all')';
%         %sem1 = nanstd(avgConditions1,0,2)' ./ sqrt(sum(~isnan(avgConditions1),2)');
%         %sem2 = nanstd(avgConditions2,0,2)' ./ sqrt(sum(~isnan(avgConditions2),2)');
%         
%         % Plot the data on a polar plot
%         subplot(2,5,ii)
%         
%         for subjectIndex = 1:size(medianBOLDpa, 4)
%             scatter(1, vals_1(subjectIndex),  30, 'MarkerFaceColor', colors(subjectIndex,:), 'MarkerEdgeColor', 'none'); %, [176/255, 224/255, 230/255]
%             hold on
%             scatter(2, vals_2(subjectIndex), 30, 'MarkerFaceColor', colors(subjectIndex,:), 'MarkerEdgeColor', 'none'); %, 'r'
%             plot([1 2], [vals_1(subjectIndex) vals_2(subjectIndex)], 'Color', colors(subjectIndex,:));   %'k')
%             xlim([0 3])
%             %ylim([-0.15 0.26])
%         end
%         
%         scatter(1, vals_1_overall, 70, 'MarkerFaceColor', [176/255, 224/255, 230/255], 'MarkerEdgeColor', 'black');
%         hold on
%         scatter(2, vals_2_overall,  70, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'black');
%         plot([1 2], [vals_1_overall vals_2_overall], 'k', 'LineWidth', 3)
%         
%         title(sprintf('BOLD per polar angle in %s',rois{regionIndex}));
%         ylabel('PSC')
%         set(gca, 'XTick', []);
%     
%         if ii==1
%             legend('Radial', 'Tangential', 'Location', 'northeast');
%         end
%     
%         hold off;
%     end
%     
%     fa = gcf;
%     fa.Position = [1000 555 1514 782];
%     
%     % filename = fullfile(figureDir,sprintf('PAradialvtangmotionVSstatic%s', rois{regionIndex}));
%     % saveas(gcf, filename, 'pdf');
% end
