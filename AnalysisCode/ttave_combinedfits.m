%ttave_combinefits.m

clear all; close all; clc;

% must be in DriftingGrating directory to run

% setup path
addpath(genpath(pwd));
bidsDir =  '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
figDir = strrep(bidsDir, 'data_bids', 'figures');
%bidsDir = '/Volumes/server/Projects/Project_dg/data_bids/';
githubDir = '~/Documents/GitHub';
fsDir = '/Applications/freesurfer/7.2.0';
addpath(genpath(fullfile(githubDir, 'wpToolbox')));
addpath(genpath(fullfile(githubDir, 'atlasmgz')));
setup_user('rania', bidsDir)

projectSettings = loadConfig(githubDir);
rois = projectSettings.rois;
clear projectSettings;
comparisonName = 'motion_minus_baseline';

projects = {'dg', 'da'};
hRF_setting = 'glmsingle'; %
asymmetries = {'cartesianCardinalvsCartesianOblique', 'polarCardinalvsPolarOblique', ...
    'radialVsTangential'};

% if strcmp(projectName, 'dg')
%     % to leave out one subject
%     subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-0397', ...
%         'sub-0442', 'sub-wlsubj121', ...
%         'sub-wlsubj123', 'sub-wlsubj124', 'sub-wlsubj127', 'sub-0395', 'sub-0426', ...
%         'sub-0250'};
% elseif strcmp(projectName, 'da')
    % load subjects
subjects = {'sub-0037', 'sub-0201', 'sub-0255', 'sub-wlsubj123', 'sub-wlsubj124', ...
    'sub-0395', 'sub-0426'};
% end

% make a three 3D arrays (1 for each asymmetry)
% rows (4x) 
% cols (x25) timepoints

nCurves = 4; % conditions (motion pro, motion con, orientation pro, orientation con)
nSubjects = numel(subjects);
nAsymmetries = numel(asymmetries);
nTimepoints = 25;

%% Aggregate

for ri=1:numel(rois)  % specific to the ROI

    roiName = rois{ri};

    for si=1:numel(subjects)
    
        subj = subjects{si};
    
        ttave_combinedPath = fullfile(bidsDir, 'derivatives', 'ttaveCombinedfits', subj);
    
        % this overrides 'projectSettings'
        % save the variables: 'ttave_arrays', 'projectSettings', 'subj'
        load(fullfile(ttave_combinedPath,sprintf('ttaveCombined_%s',comparisonName)));

        for pi=1:numel(projects)    % will take the mean over projects
        
            projectName = projects{pi};
    
            for as=1:numel(asymmetries)

                asymmetryName = asymmetries{as};
    
                % initialize group array if does not exist
                try
                 retrived_fits = group_ttave_arrays.(roiName).(asymmetryName).(strcat(projectName, '_fits'));
                 retrived_rawvals = group_ttave_arrays.(roiName).(asymmetryName).(projectName);
                catch
                 group_ttave_arrays.(roiName).(asymmetryName).(strcat(projectName, '_fits')) = nan(nCurves, nTimepoints, nSubjects);
                 
                 group_ttave_arrays.(roiName).(asymmetryName).(projectName) = nan(nCurves, nTimepoints, nSubjects);
                end

                fits = ttave_arrays.(roiName).(asymmetryName).(strcat(projectName, '_fits'));
                group_ttave_arrays.(roiName).(asymmetryName).(strcat(projectName, '_fits'))(:,:,si) = fits;

                rawvals = ttave_arrays.(roiName).(asymmetryName).(projectName);
                group_ttave_arrays.(roiName).(asymmetryName).(projectName)(:,:,si) = rawvals;

                group_ttave_arrays.(roiName).(asymmetryName).colors = ttave_arrays.(roiName).(asymmetryName).colors;
                group_ttave_arrays.(roiName).(asymmetryName).legendLabels = ttave_arrays.(roiName).(asymmetryName).legendLabels;

            end
        
        end
    end
end

%% Now average and plot

% asymmetries

ttave_combinedGroupPath = fullfile(figDir, 'combined');

if ~isfolder(ttave_combinedGroupPath)
    mkdir(ttave_combinedGroupPath)
end

for ri=1:numel(rois)  % specific to the ROI

    figure
    roiName = rois{ri};

    iter=0;

    for ri=1:numel(asymmetries)

        asymmetryName = asymmetries{ri};

        colors = group_ttave_arrays.(roiName).(asymmetryName).colors;
        legendLabels = group_ttave_arrays.(roiName).(asymmetryName).legendLabels;

        for pi=1:numel(projects)    % will take the mean over projects

            projectName = projects{pi};
        
            iter=iter+1;
            subplot(3,2,iter)

            % copied over for now; condense with ttave_fit LATER
            rawVals = mean(group_ttave_arrays.(roiName).(asymmetryName).(projectName),3);
            fittedVals = mean(group_ttave_arrays.(roiName).(asymmetryName).(strcat(projectName, '_fits')),3);

            advMotionVals = rawVals(1,:);
            disadvMotionVals = rawVals(2,:);
            advStaticVals = rawVals(3,:);
            disadvStaticVals = rawVals(4,:);
    
            fitted_advMotionVals = fittedVals(1,:);
            fitted_disadvMotionVals = fittedVals(2,:);
            fitted_advStaticVals = fittedVals(3,:);
            fitted_disadvStaticVals = fittedVals(4,:);
    
            hold on
            p1=plot(fitted_advMotionVals, 'b-', 'Linewidth',1, 'Color', colors{1});
            p2=plot(fitted_disadvMotionVals, 'b-', 'Linewidth',1, 'Color', colors{2});
            p3=plot(fitted_advStaticVals, 'b--', 'Linewidth',1, 'Color', colors{1});
            p4=plot(fitted_disadvStaticVals, 'b--', 'Linewidth',1, 'Color', colors{2});

            scatter(1:length(advMotionVals), advMotionVals, 25, 'filled', 'Marker', 'd', 'MarkerFaceColor', colors{1});
            scatter(1:length(disadvMotionVals), disadvMotionVals, 25, 'filled', 'Marker', 'd', 'MarkerFaceColor', colors{2});
            scatter(1:length(advStaticVals), advStaticVals, 25, 'Marker', 'd', 'MarkerEdgeColor', colors{1}, 'LineWidth',1);
            scatter(1:length(disadvStaticVals), disadvStaticVals, 25, 'Marker', 'd', 'MarkerEdgeColor', colors{2}, 'LineWidth',1);
            hold off

            if ismember(iter, [1,3,5])
                legend([p1,p2,p3,p4],legendLabels, 'Location', 'northeast', 'Fontsize', 6)
                ylabel('BOLD (psc)')
            end

            if strcmp(projectName, 'dg')
                projectDescription = 'Cartesian grating stimulus';
            elseif strcmp(projectName, 'da')
                projectDescription = 'Polar grating stimulus';
            end

            if ismember(iter, [1,2])
                title(sprintf('%s', projectDescription))
            end

            if ismember(iter, [5,6])
                xlabel('TRs')
                current_ticks = xticks();
                new_labels = current_ticks - projectSettings.eventTRs_prior;
                xticklabels(string(new_labels));
            else
                xticklabels([]);
            end

            %plot(grandMean, '-', 'Linewidth',3, 'Color', 'k')
            ylim([-.5 1.75])
            set(gca, 'FontName', 'Arial', 'FontSize', 12)



        end

    end

    sgtitle(roiName)
    fig = gcf;
    fig.Position = [1 306 1211 1031];
    gcf_edit = fitFig2Page(gcf);

    % Save as PDF
    print(gcf_edit, fullfile(ttave_combinedGroupPath,sprintf('ttaveCombined_%s_%s',roiName,projectSettings.comparisonName)), '-dpdf');

    close all;

end
