function ttave(matrices,datafiles, surfaceROI, tr_s, stimlength)

% cardinal motion
% cardinal static
% oblique motion
% oblique static
% blank
% then plot baseline

padding = 10; % TRs
trialcut = 12; % TRs

% initialize matrices for each event type: c_motion, o_motion, c_static,
% o_static, blank
% 10 sec 
cMotion = []; oMotion = []; cStatic = []; oStatic = []; blank = []; base = [];

for mi=1:length(matrices)
    m = matrices{mi};
    df = datafiles{mi};

    % convert to psc?
    flatrate = mean(df, 2);
    timeseries_psc = (((df-flatrate)./flatrate).*100);

    df_roi = nanmean(timeseries_psc.*surfaceROI);
    %baseline = mean([df_roi(1:padding), df_roi(end-padding+1:end)]);
    
    pad = [df_roi(1:padding) ; df_roi(end-padding+1:end)];
    base = [base ; pad];

    % get unique conditions
    [~, cond_n] = size(m);

    for ci=1:cond_n

        idx1 = find(m(:,ci)==1); 
        idxStart1 = idx1(1:stimlength:end,:); % every third element (to get stim start)
    
        for ti=1:length(idxStart1)
            trial = df_roi(idxStart1(ti):idxStart1(ti)+trialcut-1);

            if (ci <= 4)
                cMotion = [cMotion ; trial];
            elseif (ci <=8) && (ci > 4)
                oMotion = [oMotion ; trial];
            elseif (ci <=10) && (ci > 8)
                cStatic = [cStatic ; trial];
            elseif (ci <=12) && (ci > 10)
                oStatic = [oStatic ; trial];
            elseif ci <=13
                blank = [blank ; trial];
            end

        end
    end

end

figure
shift = 0; %mean(base,'all');
plot(mean(cMotion)-shift, 'b-', 'Linewidth',2, 'Color', [127, 191, 123]/255)
hold on
plot(mean(oMotion)-shift, '-', 'Linewidth',2, 'Color', [175, 141, 195]/255)
hold on
plot(mean(cStatic)-shift, ':', 'Linewidth',2, 'Color', [127, 191, 123]/255)
hold on
plot(mean(oStatic)-shift, ':', 'Linewidth',2, 'Color', [175, 141, 195]/255)
hold on
plot(mean(blank)-shift, 'k', 'Linewidth',2, 'Color', [175, 175, 175]/255)
hold on
yline(mean(base,'all')-shift, '--')
title('Trial triggered average')
xlim([1, 12])
xlabel('Time (TR)')
ylabel('% Signal Change')
ax = gca;
ax.FontSize = 14;
legend({'Cardinal Motion', 'Oblique Motion', 'Cardinal Static', 'Oblique Static', 'Blank', 'Baseline'}, 'Location', 'southwest')

disp('end')

end