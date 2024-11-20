function [ttaveOutput, ttaveTime] = ttave_compute(matrices,datafiles, surfaceROI, tr_s, stimlength, subj, eventTRs_prior, eventTRs_after, varargin)

if isempty(varargin)
    prfData = [];  % empty (when not using polar angle)
    polarAngleConditions = {[0:45:315]};
    comparison = 'cardinalOblique'
else
    prfData = varargin{1};  % matrix of polar angle BINS (when using polar angle)
    polarAngleConditions = {[0,180], [45, 225], [90, 270], [135, 315]};
    comparison = 'radialTang'
end

% this is the order of stimulus cols (1-13) -- first 8 are motion, next 4
% are static
stimulusOrder = {0, 90, 180, 270, 45, 135, 225, 315, 0, 90, 45, 135, nan};

padding = 10; % TRs

% initialize matrices for each event type: c_motion, o_motion, c_static,
% o_static, blank
% 10 sec 
advMotion = []; disadvMotion = []; advStatic = []; disadvStatic = []; blank = []; base = [];

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

    % update surface ROI to inlcude only 1 polar angle at a time
    for pa=1:length(polarAngleConditions) 

        % retrieve polar angle value pairs (or all)
        angleVals = polarAngleConditions{pa};

        if ~isempty(prfData)
            surfacePA = nan(length(prfData),1);
            paIdx = find(ismember(prfData, angleVals));
            surfacePA(paIdx) = 1;
            surfaceSelection = surfaceROI.*surfacePA;
        else
            surfaceSelection = surfaceROI;
        end

        df_roi = nanmean(timeseries_psc.*surfaceSelection);
        %baseline = mean([df_roi(1:padding), df_roi(end-padding+1:end)]);
        
        pad = [df_roi(1:padding) ; df_roi(end-padding+1:end)];
        base = [base ; pad];
    
        % get unique conditions
        [~, cond_n] = size(m);
    
        for ci=1:cond_n
    
            idx1 = find(m(:,ci)==1); 
            %idxStart1 = idx1(1:stimlength:end,:); % every third element (to get stim start)
            idxStart1 = idx1;
        
            for ti=1:length(idxStart1)
    
                try
                    % start -5 TRs before event
                    trial = df_roi(idxStart1(ti)-eventTRs_prior:idxStart1(ti)+eventTRs_after-1);
                catch
                    disp('last trial') % adding extra nans for the trial at the very end b/c eventTRs_after is too long
                    trial = df_roi(idxStart1(ti)-eventTRs_prior:end);
                    trial = [trial, nan(1,eventTRs_prior+eventTRs_after-length(trial))];
                end
    
                if ci ==13
                    blank = [blank ; trial]; % this is the same across all conditions
                end

                if strcmp(comparison, 'cardinalOblique')
                    if (ci <= 4)
                        advMotion = [advMotion ; trial];
                    elseif (ci <=8) && (ci > 4)
                        disadvMotion = [disadvMotion ; trial];
                    elseif (ci <=10) && (ci > 8)
                        advStatic = [advStatic ; trial];
                    elseif (ci <=12) && (ci > 10)
                        disadvStatic = [disadvStatic ; trial];
                    end
                elseif strcmp(comparison, 'radialTang')
                    if ci <= 8 % motion
                        if all(ismember(abs(angleVals-stimulusOrder{ci}), [0 180])) % only looking at a given polar anlge (RADIAL)
                            advMotion = [advMotion ; trial];
                        elseif all(ismember(abs(angleVals-stimulusOrder{ci}), [90 270])) % TANGENTIAL 
                            disadvMotion = [disadvMotion ; trial];
                        end
                    elseif ci > 8 && ci <= 13 % static
                        if all(ismember(abs(angleVals-stimulusOrder{ci}), [0 180])) % only looking at a given polar anlge (RADIAL)
                            advStatic = [advStatic ; trial];
                        elseif all(ismember(abs(angleVals-stimulusOrder{ci}), [90 270])) % TANGENTIAL 
                            disadvStatic = [disadvStatic ; trial];
                        end
                    end
                end
    
            end
        end
    end
end

ttaveOutput = {base, blank, advMotion, disadvMotion, advStatic, disadvStatic};
ttaveTime = {eventTRs_prior, eventTRs_after};

end