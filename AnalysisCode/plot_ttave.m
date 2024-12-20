function plot_ttave(ttaveOutput, projectSettings, asymmetryName)

    subj = projectSettings.subject;
    roiName = projectSettings.roiName;
    comparisonName = projectSettings.comparisonName;
    eventTRs_prior = projectSettings.eventTRs_prior;
    eventTRs_after = projectSettings.eventTRs_after;
    projectName = projectSettings.projectName;
    colors_data = projectSettings.colors_data;
    radialvstang = projectSettings.radialvstang;
    legendOn=0;

    colors = {colors_data.conditions.(projectName).(asymmetryName).color_pro', ...
        colors_data.conditions.(projectName).(asymmetryName).color_con'};

    if strcmp(comparisonName, 'motion_minus_baseline') || strcmp(comparisonName, 'orientation_minus_baseline')
        if ~radialvstang
            legendLabels = {'Cardinal Motion', 'Oblique Motion', 'Cardinal Static', 'Oblique Static', 'Mean baseline (padding)'};
        else
            legendLabels = {'Radial Motion', 'Tangential Motion', 'Radial Static', 'Tangential Static', 'Mean baseline (padding)'};
        end
    elseif strcmp(comparisonName, 'motion_minus_orientation') 
        if ~radialvstang
            legendLabels = {'Cardinal Motion - Card Ori', 'Oblique Motion - Card Ori', 'Mean baseline (padding)'};
        else
            legendLabels = {'Radial Motion - Tang Ori', 'Tangential Motion - Rad Ori', 'Mean baseline (padding)'};
        end
    end

    base = ttaveOutput{1};
    blank = ttaveOutput{2};
    advMotion = ttaveOutput{3};
    disadvMotion = ttaveOutput{4};
    advStatic = ttaveOutput{5};
    disadvStatic = ttaveOutput{6};

    figure
    shift = nanmean(blank); %mean(base,'all');

    advMotionVals = nanmean(advMotion)-shift;
    disadvMotionVals = nanmean(disadvMotion)-shift;
    advStaticVals = nanmean(advStatic)-shift;
    disadvStaticVals = nanmean(disadvStatic)-shift;

    % recompute by subtracting orientation (opposite for radialtang)
    if strcmp(comparisonName, 'motion_minus_orientation') && ~radialvstang
        advMotionVals = advMotionVals - (nanmean(advStatic)-shift);
        disadvMotionVals = disadvMotionVals - (nanmean(disadvStatic)-shift);
    elseif strcmp(comparisonName, 'motion_minus_orientation') && radialvstang  % different b/c orthogonal
        advMotionVals = advMotionVals - (nanmean(disadvStatic)-shift);
        disadvMotionVals = disadvMotionVals - (nanmean(advStatic)-shift);
    end

    % plot mean motion (either rawValues or minusOrientation) - across
    % directions and/or subjects
    plot(advMotionVals, 'b-', 'Linewidth',3, 'Color', colors{1})
    hold on
    plot(disadvMotionVals, '-', 'Linewidth',3, 'Color', colors{2})
    hold on

    % repeat for orientation if needed
    if ~strcmp(comparisonName, 'motion_minus_orientation') % plot the raw values of orientation unless plotting motion-ori
        plot(advStaticVals, ':', 'Linewidth',3, 'Color', colors{1})
        hold on
        plot(disadvStaticVals, ':', 'Linewidth',3, 'Color', colors{2})
        hold on
    end

    %plot(mean(blank)-shift, 'k', 'Linewidth',3, 'Color', [175, 175, 175]/255)
    %hold on
    %yline(mean(base,'all')-shift, '--')
    yline(nanmean(base,'all'), '--', 'Linewidth',2, 'Color', [200, 200, 200]/255)
    hold on

    % compute and plot SEM
    if strcmp(subj,'allsubjects')
        [nSubj, nTimepoints] = size(advMotion);

        % need to reccompute pairwise difference (per subject)
        if strcmp(comparisonName, 'motion_minus_orientation')
            if ~radialvstang
                advMotionPairwise = advMotion - advStatic;
                disadvMotionPairwise = disadvMotion - disadvStatic;
            elseif radialvstang
                advMotionPairwise = advMotion - disadvStatic;
                disadvMotionPairwise = disadvMotion - advStatic;
            end
            advStaticPairwise = []; disadvStaticPairwise = [];
        elseif strcmp(comparisonName, 'motion_minus_baseline') || strcmp(comparisonName, 'orientation_minus_baseline')
            advMotionPairwise = advMotion;
            disadvMotionPairwise = disadvMotion;
            advStaticPairwise = advStatic;
            disadvStaticPairwise = disadvStatic;
        end

        advMotion_sem = std(advMotionPairwise, 0, 1) / sqrt(nSubj);
        disadvMotion_sem = std(disadvMotionPairwise, 0, 1) / sqrt(nSubj);
        advStatic_sem = std(advStaticPairwise, 0, 1) / sqrt(nSubj);
        disadvStatic_sem = std(disadvStaticPairwise, 0, 1) / sqrt(nSubj);

        fill([1:nTimepoints, fliplr(1:nTimepoints)], ...
            [advMotionVals + advMotion_sem, fliplr(advMotionVals - advMotion_sem)], ...
            colors{1}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            hold on
        fill([1:nTimepoints, fliplr(1:nTimepoints)], ...
            [disadvMotionVals + disadvMotion_sem, fliplr(disadvMotionVals - disadvMotion_sem)], ...
            colors{2}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                hold on

        if strcmp(comparisonName, 'motion_minus_baseline') || strcmp(comparisonName, 'orientation_minus_baseline')
            fill([1:nTimepoints, fliplr(1:nTimepoints)], ...
                [advStaticVals + advStatic_sem, fliplr(advStaticVals - advStatic_sem)], ...
                colors{1}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                            hold on
             fill([1:nTimepoints, fliplr(1:nTimepoints)], ...
                [disadvStaticVals + disadvStatic_sem, fliplr(disadvStaticVals - disadvStatic_sem)], ...
                colors{2}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end

    end


    title(sprintf('%s %s', subj, roiName))
    xlim([1, eventTRs_prior+eventTRs_after])
    xticks = get(gca, 'XTick');  
    new_xticklabels = xticks - eventTRs_prior;
    set(gca, 'XTickLabel', new_xticklabels);
    xlabel('time (s)')
    ylabel('% signal change')

    %[subjectyMin, subjectyMax] = getSubjectTTaxis(subj, comparisonName);

    ylim([projectSettings.subjectwise_tta_limits.(strrep(subj, '-', '_')).main.min, ...
        projectSettings.subjectwise_tta_limits.(strrep(subj, '-', '_')).main.max])
    ax = gca;
    ax.FontSize = 20;
    if legendOn || strcmp(subj, 'sub-0395')
        legend(legendLabels, 'Location', 'northeast')
    end
    f1 = gcf;
    f1.Position = [72 712 897 441];
    
    % also save figure as png & pdf
    filename = sprintf('ttaveSignal_%s_%s_%s_%s', subj, roiName, asymmetryName, comparisonName);
    ttaveSave = fullfile(projectSettings.ttaveSavePath,filename);
    saveas(gcf, sprintf('%s.png', ttaveSave));

    f1.Position = [72*0.7 712*.7 897*.7 441*.7];
    ax.FontSize = 20;
    if legendOn || strcmp(subj, 'sub-0395')
        legend(legendLabels, 'Location', 'northeast', 'Fontsize', 12)
    end
    saveas(gcf, sprintf('%s.pdf', ttaveSave));

    % save the variables
    save(ttaveSave, 'ttaveOutput', 'projectSettings', 'legendLabels', 'roiName', 'subj');


end