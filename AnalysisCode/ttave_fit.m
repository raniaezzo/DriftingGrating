% for a given subject, for a given ROI, compute the mean signal across
% the asymmetry components. Include both dg and da?

% call to plot trial triggered average

clear all; close all; clc;

% must be in DriftingGrating directory to run

% setup path
addpath(genpath(pwd));
bidsDir =  '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
%bidsDir = '/Volumes/server/Projects/Project_dg/data_bids/';
githubDir = '~/Documents/GitHub';
fsDir = '/Applications/freesurfer/7.2.0';
addpath(genpath(fullfile(githubDir, 'wpToolbox')));
addpath(genpath(fullfile(githubDir, 'atlasmgz')));
setup_user('rania', bidsDir)

projectSettings = loadConfig(githubDir);

hRF_setting = 'glmsingle'; %
subj = 'sub-0395'; %'sub-0426';
%ses = 'ses-03'; %'ses-01'; %'ses-nyu3t02'; %'ses-01';

projects = {'dg','da'};
comparisonName = 'motion_minus_baseline'; % technically will read in the same raw data as 'motion_minus_orientation'
asymmetries = {'mainCardinalVsMainOblique', 'derivedCardinalVsDerivedOblique', 'radialVsTangential'};

plonOn=0;

ttave_combinedPath = fullfile(bidsDir, 'derivatives', 'ttaveCombinedfits', subj);

if ~isfolder(ttave_combinedPath)
    mkdir(ttave_combinedPath)
end

for ri=1:numel(projectSettings.rois)  % and specific to the ROI

    roiName = projectSettings.rois{ri};

    for ai=1:numel(asymmetries)   % will take the mean for that specific asymmetry
    
        asymmetryName = asymmetries{ai};

         for pi=1:numel(projects)    % will take the mean over projects

            projectName = projects{pi};

            % get ses name
            subjectDir = fullfile(bidsDir,'derivatives',strcat(projectName, 'GLM'), strcat('hRF_',hRF_setting), subj);
            contents = dir(subjectDir);
            sesNames = {};
            for i = 1:length(contents)
                % Check if the item is a directory and starts with 'ses-'
                if contents(i).isdir && startsWith(contents(i).name, 'ses-')
                    % Add the subfolder name to the cell array
                    sesNames{end+1} = contents(i).name;
                end
            end
            if length(sesNames)>1
                error('Multiple session folder.')
            end
            %

            % note: can just load in .mat file from motion_minus_baseline b/c all raw
            % conditions are saved as ttaveOutput
            filename = sprintf('ttaveSignal_%s_%s_%s_%s.mat', subj, roiName, asymmetryName, comparisonName);
            load(fullfile(subjectDir, sesNames{1}, 'ttaveData', comparisonName, filename))

            % recompute the values of interest
            base = ttaveOutput{1};
            blank = ttaveOutput{2};
            advMotion = ttaveOutput{3};
            disadvMotion = ttaveOutput{4};
            advStatic = ttaveOutput{5};
            disadvStatic = ttaveOutput{6};
            
            shift = nanmean(blank); %mean(base,'all');
            
            advMotionVals = nanmean(advMotion)-shift;
            disadvMotionVals = nanmean(disadvMotion)-shift;
            advStaticVals = nanmean(advStatic)-shift;
            disadvStaticVals = nanmean(disadvStatic)-shift;
            
            % motion - orientation (will probably not use b/c hRFs are
            % flatter) and the whole point is to summarize wach with beta
            % weight anyway
            if ~projectSettings.radialvstang
                advMotionVals_minus_orientation = advMotionVals - (nanmean(advStatic)-shift);
                disadvMotionVals_minus_orientation = disadvMotionVals - (nanmean(disadvStatic)-shift);
            elseif projectSettings.radialvstang  % different b/c orthogonal
                advMotionVals_minus_orientation = advMotionVals - (nanmean(disadvStatic)-shift);
                disadvMotionVals_minus_orientation = disadvMotionVals - (nanmean(advStatic)-shift);
            end
            
            % plot mean motion (either rawValues or minusOrientation) - across
            % directions and/or subjects
            colors = {projectSettings.colors_data.conditions.(projectName).(asymmetryName).color_pro', ...
                    projectSettings.colors_data.conditions.(projectName).(asymmetryName).color_con'};

            if plonOn == 1
                figure
                plot(advMotionVals, 'b-', 'Linewidth',3, 'Color', colors{1})
                hold on
                plot(disadvMotionVals, '-', 'Linewidth',3, 'Color', colors{2})
                hold on
                plot(advStaticVals, 'b--', 'Linewidth',3, 'Color', colors{1})
                hold on
                plot(disadvStaticVals, '--', 'Linewidth',3, 'Color', colors{2})
                title(projectSettings.projectName)
            end

            if (strcmp(asymmetryName, 'mainCardinalVsMainOblique') && strcmp(projectName, 'dg')) || ...
                (strcmp(asymmetryName, 'derivedCardinalVsDerivedOblique') && strcmp(projectName, 'da'))
                renamedAsymmetry = 'cartesianCardinalvsCartesianOblique';
            elseif (strcmp(asymmetryName, 'derivedCardinalVsDerivedOblique') && strcmp(projectName, 'dg')) || ...
                (strcmp(asymmetryName, 'mainCardinalVsMainOblique') && strcmp(projectName, 'da'))
                renamedAsymmetry = 'polarCardinalvsPolarOblique';
            elseif (strcmp(asymmetryName, 'radialVsTangential'))
                renamedAsymmetry = asymmetryName;
            end

            ttave_arrays.(roiName).(renamedAsymmetry).(projectName) = [advMotionVals; disadvMotionVals; advStaticVals; disadvStaticVals];

         end
    end

    renamedAsymmetries = {'cartesianCardinalvsCartesianOblique', 'polarCardinalvsPolarOblique', 'radialVsTangential'};
    
    iter=0;

    for ri=1:numel(renamedAsymmetries)

        renamedAsymmetry = renamedAsymmetries{ri};

        % frame color scheme in this new combined (absolute) condition
        if strcmp(renamedAsymmetry, 'cartesianCardinalvsCartesianOblique')
            colors = {projectSettings.colors_data.conditions.('dg').('mainCardinalVsMainOblique').color_pro', ...
                    projectSettings.colors_data.conditions.('dg').('mainCardinalVsMainOblique').color_con'};
            legendLabels = {'cartCard Motion', 'cartObl Motion', 'cartCard Orientation', 'cartObl Orientation'};
        elseif strcmp(renamedAsymmetry, 'polarCardinalvsPolarOblique')
            colors = {projectSettings.colors_data.conditions.('dg').('derivedCardinalVsDerivedOblique').color_pro', ...
                    projectSettings.colors_data.conditions.('dg').('derivedCardinalVsDerivedOblique').color_con'};
             legendLabels = {'polCard Motion', 'polObl Motion', 'polCard Orientation', 'polObl Orientation'};
        elseif strcmp(renamedAsymmetry, 'radialVsTangential')
            colors = {projectSettings.colors_data.conditions.('dg').('radialVsTangential').color_pro', ...
                    projectSettings.colors_data.conditions.('dg').('radialVsTangential').color_con'};
            legendLabels = {'radial Motion', 'tangential Motion', 'radial Orientation', 'tangential Orientation'};
        end

        ttave_arrays.(roiName).(renamedAsymmetry).grandMean = nanmean([ttave_arrays.(roiName).(renamedAsymmetry).dg; ttave_arrays.(roiName).(renamedAsymmetry).da]);

        % additionally save colors, legendLabels
        ttave_arrays.(roiName).(renamedAsymmetry).colors = colors;
        ttave_arrays.(roiName).(renamedAsymmetry).legendLabels = legendLabels;

        for pi=1:numel(projects)    % will take the mean over projects

            projectName = projects{pi};
        
            iter=iter+1;
            subplot(3,2,iter)

            grandMean = ttave_arrays.(roiName).(renamedAsymmetry).grandMean;
            
            advMotionVals = ttave_arrays.(roiName).(renamedAsymmetry).(projectName)(1,:);
            disadvMotionVals = ttave_arrays.(roiName).(renamedAsymmetry).(projectName)(2,:);
            advStaticVals = ttave_arrays.(roiName).(renamedAsymmetry).(projectName)(3,:);
            disadvStaticVals = ttave_arrays.(roiName).(renamedAsymmetry).(projectName)(4,:);
    
            fitted_advMotionVals = findScalar(grandMean, advMotionVals);
            fitted_disadvMotionVals = findScalar(grandMean, disadvMotionVals);
            fitted_advStaticVals = findScalar(grandMean, advStaticVals);
            fitted_disadvStaticVals = findScalar(grandMean, disadvStaticVals);
    
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

            % prepare to save
            fits = [fitted_advMotionVals; fitted_disadvMotionVals; fitted_advStaticVals; fitted_disadvStaticVals];
            ttave_arrays.(roiName).(renamedAsymmetry).(strcat(projectName, '_fits')) = fits;

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
            ylim([-.5 2.5])
            set(gca, 'FontName', 'Arial', 'FontSize', 12)

        end
    end

    sgtitle(roiName)
    fig = gcf;
    fig.Position = [1 306 1211 1031];
    gcf_edit = fitFig2Page(gcf);
    
    % Save as PDF
    print(gcf_edit, fullfile(ttave_combinedPath,sprintf('ttaveCombined_%s_%s',roiName,projectSettings.comparisonName)), '-dpdf');

    close all;
end

% save the variables
save(fullfile(ttave_combinedPath,sprintf('ttaveCombined_%s',projectSettings.comparisonName)), ...
    'ttave_arrays', 'projectSettings', 'subj'); %, 'legendLabels', 'roiName','renamedAsymmetry', 'colors');





