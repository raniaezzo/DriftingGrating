%% Individual ROI plots (I think I can delete?)

% ci_level = 68; %68;
% colors = [[127 191 123]/255; [166 97 26]/255; [146 197 222]/255];
% colors2 = [[175 141 195]/255; [64 176 166]/255; [202 0 32]/255];
% meanRelative=0; %1;
% 
% for roi=1:length(rois)
%     % load('bootLME_sensitivity.mat')
%     % load('LME_sensitivity.mat')
% 
%     %saveDir = fullfile(glmResultsfolder,'LME_results', rois{roi});
% 
%     load(strcat(saveDir, '/modeldata'), 'modeldata');
%     load(fullfile(saveDir,strcat('LME_',metric)), 'estimates');
%     load(strcat(saveDir, '/boot'), 'saveboot','coeffs');
% 
%     Gintercept = estimates(1);
%     abs_cardinal_est = Gintercept + estimates(2);
%     rel_cardinal_est = Gintercept + estimates(3);
%     radial_est = Gintercept + estimates(4);
% 
%     CIFcn = @(x,p)prctile(x, [100-p, p]); p = ci_level;
% 
%     % coefficient order is: intercept, abs_cardinality, rel_cardinality, radiality
%     CI_abscardinality = CIFcn(coeffs(2,:)+Gintercept,p);
%     CI_relcardinality = CIFcn(coeffs(3,:)+Gintercept,p);
%     CI_radiality = CIFcn(coeffs(4,:)+Gintercept,p);
% 
%     
%     %
%     
% %     figure
% % 
% %     subplot(3,1,1)
% %     xline(abs_cardinal_est, 'k--', 'linewidth', 4)
% %     hold on
% %     histogram(coeffs(2,:)+Gintercept, 'FaceColor', colors(1,:))
% %     hold on
% %     xline(CI_abscardinality(1), 'r--', 'linewidth', 4)
% %     hold on
% %     xline(CI_abscardinality(2), 'r--', 'linewidth', 4)
% %     hold on
% %     xline(Gintercept, 'y--', 'linewidth', 4)
% %     title('abs cardinal')
% %     %xlim([xmin xmax])
% % 
% %     subplot(3,1,2)
% %     xline(rel_cardinal_est, 'k--', 'linewidth', 4)
% %     hold on
% %     histogram(coeffs(3,:)+Gintercept, 'FaceColor', colors(2,:))
% %     hold on
% %     xline(CI_relcardinality(1), 'r--', 'linewidth', 4)
% %     hold on
% %     xline(CI_relcardinality(2), 'r--', 'linewidth', 4)
% %     hold on
% %     xline(Gintercept, 'y--', 'linewidth', 4)
% %     title('rel cardinal')
% %     %xlim([xmin xmax])
% % 
% %     subplot(3,1,3)
% %     xline(radial_est, 'k--', 'linewidth', 4)
% %     hold on
% %     histogram(coeffs(4,:)+Gintercept, 'FaceColor', colors(3,:))
% %     hold on
% %     xline(CI_radiality(1), 'r--', 'linewidth', 4)
% %     hold on
% %     xline(CI_radiality(2), 'r--', 'linewidth', 4)
% %     hold on
% %     xline(Gintercept, 'y--', 'linewidth', 4)
% %     title('radial')
% %     %xlim([xmin xmax])
% % 
% %     sgtitle(sprintf('mean betas + %s CIs', num2str(ci_level)))
% 
%     % sensitivity (manuscript figure)
% 
%     y = [estimates(2) estimates(3) estimates(4)]; %[0.11446, 0.033442, 0.015738];
%     %errlow = [0.0053434, 0.0053434, 0.0075567]; % output from model
%     %errhigh = [0.0053434, 0.0053434, 0.0075567];
% 
%     errlow = [CI_abscardinality(1) CI_relcardinality(1) CI_radiality(1)];
%     errhigh = [CI_abscardinality(2) CI_relcardinality(2) CI_radiality(2)];
% 
%     y1 = Gintercept + y;
%     y2 = Gintercept - y;
% 
%     errlow1 = y1-errlow;
%     errhigh1 = errhigh -y1;
% 
%     errlow2 = errhigh1; %errlow1;
%     errhigh2 = errlow1; %errhigh1;
% 
%     x = 1:3;
%     figure();
%     ax = axes();
%     hold(ax);
%     if meanRelative
%         baselineSub = Gintercept;
%     else
%         baselineSub = 0;
%     end
%     for i=1:length(x)
%         boxchart(x(i)*ones(size(y1(:,i))), y1(:,i)-baselineSub, 'BoxFaceColor', colors(i,:), 'LineWidth', 3, 'BoxWidth', 1)
%         plot(x(i)*ones(size(y1(:,i))), y1(:,i)-baselineSub, 'Color', colors(i,:), 'Marker', '.', 'MarkerSize', 10, 'LineStyle','none')
%         hold on
%         errorbar(x(i),y1(i)-baselineSub,errlow1(i)-baselineSub, errhigh1(i)-baselineSub, 'LineStyle','none', 'LineWidth', 2, 'Color', colors(i,:));
%         hold on
%         boxchart(x(i)*ones(size(y2(:,i))), y2(:,i)-baselineSub, 'BoxFaceColor', colors2(i,:), 'LineWidth', 3, 'BoxWidth', 1)
%         hold on
%         errorbar(x(i),y2(i)-baselineSub,errlow2(i)-baselineSub, errhigh2(i)-baselineSub, 'LineStyle','none', 'LineWidth', 2, 'Color', colors2(i,:));
%         hold on
%     end
%     if meanRelative
%         yline(0, '--', 'Color', [0 0 0], 'LineWidth', 2)
%         ylim([-0.05 0.05])
%     else
%         yline(Gintercept, '--', 'Color', [0 0 0], 'LineWidth', 2)
%         ylim([-0.05+Gintercept 0.05+Gintercept])
%     end
%     xlim([0.5 3.5])
%     %ylim([ymin ymax])
%     set(gca,'XTick',[])
%     box off
%     set(gca,'linewidth',2, 'YColor', [0 0 0]);
%     set(gca,'linewidth',2, 'XColor', [0 0 0]);
%     title(rois(roi))
%     set(gca, 'FontName', 'Arial', 'FontSize', 12);
%     
%     
% end
% hold off