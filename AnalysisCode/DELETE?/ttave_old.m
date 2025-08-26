function ttave(matrices,datafiles, surfaceROI, tr_s, stimlength, subj)

% cardinal motion
% cardinal static
% oblique motion
% oblique static
% blank
% then plot baseline

padding = 10; % TRs
eventTRs_after = 20; % TRs
eventTRs_prior = 5;

% initialize matrices for each event type: c_motion, o_motion, c_static,
% o_static, blank
% 10 sec 
cMotion = []; oMotion = []; cStatic = []; oStatic = []; blank = []; base = [];

for mi=1:length(matrices)
    m = matrices{mi};
    df = datafiles{mi};

    % convert to psc?
    %flatrate = mean(df, 2);
    %timeseries_psc = (((df-flatrate)./flatrate).*100);

    % abs (to prevent flipping sign) - done in GLMestimatemodel
    mean_overT = mean(abs(df),2); 
    % doing PSC here - subtract 1 to remove baseline signal
    timeseries_psc = ((df./mean_overT)-1)*100;

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

            try
                % start -5 TRs before event
                trial = df_roi(idxStart1(ti)-eventTRs_prior:idxStart1(ti)+eventTRs_after-1);
            catch
                disp('last trial') % adding extra nans for the trial at the very end b/c eventTRs_after is too long
                trial = df_roi(idxStart1(ti)-eventTRs_prior:end);
                trial = [trial, nan(1,eventTRs_prior+eventTRs_after-length(trial))];
            end

            if (ci <= 4)
                cMotion = [cMotion ; trial];
            elseif (ci <=8) && (ci > 4)
                oMotion = [oMotion ; trial];
            elseif (ci <=10) && (ci > 8)
                cStatic = [cStatic ; trial];
            elseif (ci <=12) && (ci > 10)
                oStatic = [oStatic ; trial];
            elseif ci ==13
                blank = [blank ; trial];
            end

        end
    end

end



end