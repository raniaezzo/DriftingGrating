function ttave_compute(matrices, datafiles, asymmetryName, surfaceROI, projectSettings, varargin)

eventTRs_prior = projectSettings.eventTRs_prior;
eventTRs_after = projectSettings.eventTRs_after;
padding = projectSettings.padding; % TRs
comparisonName = projectSettings.comparisonName;
projectName = projectSettings.projectName;

% Check project name & request
if strcmp(projectSettings.projectName, 'da')
    if strcmp(asymmetryName, 'mainCardinalVsMainOblique')
        derivedVals = 0; % this is needed because radialVsTang occurs for both conditions
        % Ensure a sixth argument is provided
        if nargin < 6 || isempty(varargin{1})
            error(sprintf('A sixth input is required when projectName is "%s" and asymmetryName is "%s".', projectSettings.projectName, asymmetryName));
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
        % Sixth argument is ignored
        derivedVals = 0; 
        radialvstang = 0;
    elseif strcmp(asymmetryName, 'derivedCardinalVsDerivedOblique')
        derivedVals = 1;
        % Ensure a sixth argument is provided
        if nargin < 6 || isempty(varargin{1})
            error(sprintf('A sixth input is required when projectName is "%s" and asymmetryName is "%s".', projectSettings.projectName, asymmetryName));
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

% just to get code to run, remove or generalize?
% if strcmp(projectSettings.projectName, 'dg') && strcmp(asymmetryName, 'mainCardinalVsMainOblique')
%     prfData = [];  % empty (when not using polar angle)
%     polarAngleConditions = {[0:45:315]};
% else %if strcmp(projectSettings.projectName, 'dg') && strcmp(asymmetryName, 'radialVsTangential')
    prfData = projectSettings.filteredPrfBins;  % matrix of polar angle BINS (when using polar angle)
    polarAngleConditions = {[0,180], [45, 225], [90, 270], [135, 315]};
% end


% this is the order of stimulus cols (1-13) -- first 8 are motion, next 4
% are static
stimulusOrder = {0, 90, 180, 270, 45, 135, 225, 315, 0, 90, 45, 135, nan};

if derivedVals
    if radialvstang
        derivedProOffsets = [0, 180]; derivedConOffsets = [90, 270]; % for radial vs tangential
    else
        derivedProOffsets = [0, 180, 90, 270]; derivedConOffsets = [45, 135, 225, 315]; % for other
    end
end


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
            error('ensure polar angle values are provided')
        end

        df_roi = nanmean(timeseries_psc.*surfaceSelection);
        %baseline = mean([df_roi(1:padding), df_roi(end-padding+1:end)]);
        
        pad = [df_roi(1:padding) ; df_roi(end-padding+1:end)];
        base = [base ; pad];
    
        % get unique conditions
        [~, cond_n] = size(m);
    
        for ci=1:cond_n
    
            idx1 = find(m(:,ci)==1); 
            %idxStart1 = idx1(1:projectSettings.stimlength:end,:); % every third element (to get stim start)
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

                if ~derivedVals % for DG: mainCard; DA: mainCard; radialTang
                    % vv if NOT derived

                    % this is for DG, or DA: main Cardinal vs main Oblique
                    if strcmp(asymmetryName, 'mainCardinalVsMainOblique')
                        if (ci <= 4)
                            advMotion = [advMotion ; trial];
                        elseif (ci <=8) && (ci > 4)
                            disadvMotion = [disadvMotion ; trial];
                        elseif (ci <=10) && (ci > 8)
                            advStatic = [advStatic ; trial];
                        elseif (ci <=12) && (ci > 10)
                            disadvStatic = [disadvStatic ; trial];
                        end
                    % this is for DA (radial vs tangential)
                    elseif strcmp(asymmetryName, 'radialVsTangential') && strcmp(projectName, 'da')
                        if (ci == 2) || (ci == 4)
                            advMotion = [advMotion ; trial];
                        elseif (ci == 1) || (ci == 3)
                            disadvMotion = [disadvMotion ; trial];
                        elseif (ci == 10)
                            advStatic = [advStatic ; trial];
                        elseif (ci == 9)
                            disadvStatic = [disadvStatic ; trial];
                        end
                    end
                    % ^^ if NOT derived

                elseif derivedVals
                        % vv if derived

                        % DG: polar Cardinal & oblique; radial tangential
                        if strcmp(projectName, 'dg')
                            if ci <= 8 % motion
                                if all(ismember(abs(angleVals-stimulusOrder{ci}), derivedProOffsets)) % only looking at a given polar anlge (e.g., RADIAL or polar card)
                                    advMotion = [advMotion ; trial];
                                elseif all(ismember(abs(angleVals-stimulusOrder{ci}), derivedConOffsets)) % e.g., TANGENTIAL or Polar Oblique
                                    disadvMotion = [disadvMotion ; trial];
                                else
                                    if ci ~= 13
                                        warning('Is this supposed to not meet a condition?')
                                    end
                                end
                            elseif ci > 8 && ci <= 13 % static
                                if all(ismember(abs(angleVals-stimulusOrder{ci}), derivedProOffsets)) % only looking at a given polar angle (e.g., RADIAL or polar card)
                                    advStatic = [advStatic ; trial];
                                elseif all(ismember(abs(angleVals-stimulusOrder{ci}), derivedConOffsets)) % e.g. TANGENTIAL or Polar Oblique
                                    disadvStatic = [disadvStatic ; trial];
                                else
                                    if ci ~= 13
                                        warning('Is this supposed to not meet a condition?')
                                    end
                                end
                            end
                        % this is da: cartesian cardinal vs oblique
                        elseif strcmp(projectName, 'da')
                            if ci <= 8 % motion
                                if (all(ismember(angleVals, [0,90,180,270])) && (ci <= 4)) || ... % if on the meridians and in/out/c/cc
                                        (all(ismember(angleVals, [45,135,225,315])) && ((ci <=8) && (ci > 4))) % if on the diagonals and spiral
                                    advMotion = [advMotion ; trial];
                                elseif (all(ismember(angleVals, [0,90,180,270])) && ((ci <=8) && (ci > 4))) || ... % if on the meridians and spiral
                                        (all(ismember(angleVals, [45,135,225,315])) && (ci <= 4)) % if on the diagonals and in/out/c/cc
                                    disadvMotion = [disadvMotion ; trial];
                                else
                                    if ci ~= 13
                                        warning('Is this supposed to not meet a condition?')
                                    end
                                end
                            elseif ci > 8 && ci <= 13 % static
                                if (all(ismember(angleVals, [0,90,180,270])) && ((ci ==9) || (ci == 10))) || ... % if on the meridians and annulus/pinwheel
                                        (all(ismember(angleVals, [45,135,225,315])) && ((ci ==11) || (ci ==12))) % if on the diagonals and spiral
                                    advStatic = [advStatic ; trial];
                                elseif (all(ismember(angleVals, [0,90,180,270])) && ((ci ==11) || (ci == 12))) || ... % if on the meridians and spiral
                                        (all(ismember(angleVals, [45,135,225,315])) && ((ci ==9) || (ci ==10))) % if on the diagonals and annulus/pinwheel
                                    disadvStatic = [disadvStatic ; trial];
                                else
                                    if ci ~= 13
                                        warning('Is this supposed to not meet a condition?')
                                    end
                                end
                            end
                        end
                        % ^^ if derived

                end % if derived vs. not
    
            end % of of # of trials
        end % end of cond 1:13
    end % end of polar angle condition PAIRS
end % end of run

ttaveOutput = {base, blank, advMotion, disadvMotion, advStatic, disadvStatic};
%ttaveTime = {eventTRs_prior, eventTRs_after};

projectSettings.radialvstang = radialvstang;

plot_ttave(ttaveOutput, projectSettings, asymmetryName);

end