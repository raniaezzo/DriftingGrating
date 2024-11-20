function [legendLabels] = plot_ttave(ttaveOutput, ttaveTime, subj, roiName, ttaveSave, condName, plotValue)

    legendOn=0;

    if strcmp(plotValue, 'rawValues')
        if strcmp(condName, 'cardinaloblique')
            legendLabels = {'Cardinal Motion', 'Oblique Motion', 'Cardinal Static', 'Oblique Static', 'Mean baseline (padding)'};
            colors = {[127, 191, 123]/255, [175, 141, 195]/255};
        elseif strcmp(condName, 'radialtangential')
            legendLabels = {'Radial Motion', 'Tangential Motion', 'Radial Static', 'Tangential Static', 'Mean baseline (padding)'};
            colors = {[146 197 222]/255, [202 0 32]/255};
        end
    elseif strcmp(plotValue, 'subtractOri')
        if strcmp(condName, 'cardinaloblique')
            legendLabels = {'Cardinal Motion - Card Ori', 'Oblique Motion - Card Ori', 'Mean baseline (padding)'};
            colors = {[127, 191, 123]/255, [175, 141, 195]/255};
        elseif strcmp(condName, 'radialtangential')
            legendLabels = {'Radial Motion - Tang Ori', 'Tangential Motion - Rad Ori', 'Mean baseline (padding)'};
            colors = {[146 197 222]/255, [202 0 32]/255};
        end
    end


    eventTRs_prior = ttaveTime{1};
    eventTRs_after = ttaveTime{2};

    base = ttaveOutput{1};
    blank = ttaveOutput{2};
    advMotion = ttaveOutput{3};
    disadvMotion = ttaveOutput{4};
    advStatic = ttaveOutput{5};
    disadvStatic = ttaveOutput{6};

    % clean ROI name
    if strcmp(roiName, 'V1_REmanual')
        roiName = 'V1';
    elseif strcmp(roiName, 'pMT_REmanual')
        roiName = 'MT';
    elseif strcmp(roiName, 'hMTcomplex_REmanual')
        roiName = 'MTcomplex';
    elseif strcmp(roiName, 'pMST_REmanual')
        roiName = 'MST';
    elseif strcmp(roiName, 'V2_REmanual')
        roiName = 'V2';
    elseif strcmp(roiName, 'V3_REmanual')
        roiName = 'V3';
    end

    figure
    shift = nanmean(blank); %mean(base,'all');

    advMotionVals = nanmean(advMotion)-shift;
    disadvMotionVals = nanmean(disadvMotion)-shift;
    advStaticVals = nanmean(advStatic)-shift;
    disadvStaticVals = nanmean(disadvStatic)-shift;

    % recompute by subtracting orientation
    if strcmp(plotValue, 'subtractOri') && strcmp(condName, 'cardinaloblique')
        advMotionVals = advMotionVals - (nanmean(advStatic)-shift);
        disadvMotionVals = disadvMotionVals - (nanmean(disadvStatic)-shift);
    elseif strcmp(plotValue, 'subtractOri') && strcmp(condName, 'radialtangential')  % different b/c orthogonal
        advMotionVals = advMotionVals - (nanmean(disadvStatic)-shift);
        disadvMotionVals = disadvMotionVals - (nanmean(advStatic)-shift);
    end

    % plot mean motion (either rawValues or minusOrientation) - across
    % directions and/or subjects
    plot(advMotionVals, 'b-', 'Linewidth',3, 'Color', colors{1})
    hold on
    plot(disadvMotionVals, '-', 'Linewidth',3, 'Color', colors{2})
    hold on

    % repear for orientation if needed
    if strcmp(plotValue, 'rawValues') % plot the raw values of orientation unless plotting motion-ori
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
        if strcmp(plotValue, 'subtractOri')
            if strcmp(condName, 'cardinaloblique')
                advMotionPairwise = advMotion - advStatic;
                disadvMotionPairwise = disadvMotion - disadvStatic;
            elseif strcmp(condName, 'radialtangential')
                advMotionPairwise = advMotion - disadvStatic;
                disadvMotionPairwise = disadvMotion - advStatic;
            end
            advStaticPairwise = []; disadvStaticPairwise = [];
        elseif strcmp(plotValue, 'rawValues')
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

        if strcmp(plotValue, 'rawValues')
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
    [subjectyMin, subjectyMax] = getSubjectTTaxis(subj, plotValue);

    ylim([subjectyMin, subjectyMax])
    ax = gca;
    ax.FontSize = 20;
    if legendOn || strcmp(subj, 'sub-0395')
        legend(legendLabels, 'Location', 'northeast')
    end
    f1 = gcf;
    f1.Position = [72 712 897 441];
    
    % also save figure as png & pdf
    saveas(gcf, sprintf('%s.png', ttaveSave));

    f1.Position = [72*0.7 712*.7 897*.7 441*.7];
    ax.FontSize = 20;
    if legendOn || strcmp(subj, 'sub-0395')
        legend(legendLabels, 'Location', 'northeast', 'Fontsize', 12)
    end
    saveas(gcf, sprintf('%s.pdf', ttaveSave));


end