%% analyze drifting grating eye data

% set up

clc; clear all; close all;

% set up
addpath(genpath(pwd));
projectName = 'dg';

bidsDir =  '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
%bidsDir = '/Volumes/server/Projects/Project_dg/data_bids';
githubDir = '~/Documents/GitHub';
%fsDir = '/Applications/freesurfer'; %/7.2.0';
designDir = fullfile(strrep(bidsDir, 'data_bids', 'experimentalOutput'), projectName);
%addpath(genpath(fullfile(githubDir, 'retinotopy-nyuad')));
%addpath(genpath(fullfile(githubDir, 'wpToolbox')));
addpath(genpath(fullfile(githubDir, 'atlasmgz')));
%setup_user(projectName,bidsDir,githubDir,fsDir);
setup_user('rania',bidsDir);

tr_s = 1;
trialsPerRun = 52;
timeCutOff = 3000;

% for metadata
headings = {'SamplingRate', 'CalibrationMethod', 'Eye', 'CalibrationQuality', ...
            'ValidationQuality', 'AverageError', 'MaxError', 'OffsetDeg', 'Offpix', ...
            'NoTrackingTrials', 'Filename', 'SubjectName','RunID'};
types = {'double', 'string', 'string', 'string', 'string', 'double', 'double', ...
         'double', 'string', 'double', 'string', 'string','double'};

% for timestamps table
varNames = {'runNum', 'trialNum','sTrial', 'eTrial', 'sStim', 'eStim', 'sSacc', 'eSacc', 'sBlink', 'eBlink'}; % Define the variable names and types
    
%%

subjectslist = {'sub-wlsubj123'}; %; 'sub-wlsubj124'; 'sub-0426'; 'sub-0427'; 'sub-0395';
    %'sub-0037'; 'sub-0397'; 'sub-0255'; 'sub-0201'; 'sub-0442'; 'sub-wlsubj127'; 'sub-0250'};

% check sub-0397 & sub-0201 eye data time dimension (is 20,072?) why? -->
% 2,000 sampling rate.

for si=1:1 %12:12 %length(subjectslist)
    
    subj = subjectslist{si}; 
    
    if strcmp(subj, 'sub-wlsubj121')
        continue % NOTE: no eye data for subj121
    end

    
    % retrieve design matrices for a subject (returns 1x8 cell) that has 280
    % rows
    [matrices_init, dg_id, dg_numID, runs] = format_desmats(bidsDir, designDir, subj, tr_s);

    % handle exceptions (e.g., sub-wlsubj127 run9 EDF is corrupt
    if strcmp(subj, 'sub-wlsubj127')
        runs = runs(1:end-1);
    end

    totalRuns = numel(runs);

    eyeanalysisFolder = fullfile(bidsDir, 'derivatives', 'eyetracking', subj, sprintf('ses-%s', projectName));
    if ~isfolder(eyeanalysisFolder)
        mkdir(eyeanalysisFolder)
    end
    
    % PREALLOCATE STORAGE FOR METADATA

    metaDataTable = table('Size', [totalRuns, length(headings)], 'VariableTypes', types, 'VariableNames', headings);

    % PREALLOCATE STORAGE FOR TIMESTAMPS
    numTrials = trialsPerRun*numel(runs);
    % Initialize each column with NaN
    runNum = NaN(numTrials, 1); trialNum = NaN(numTrials, 1);
    sTrial = NaN(numTrials, 1); eTrial = NaN(numTrials, 1); sStim = NaN(numTrials, 1); eStim = NaN(numTrials, 1);
    sSacc = cell(numTrials, 1); eSacc = cell(numTrials, 1); sBlink = cell(numTrials, 1); eBlink = cell(numTrials, 1);
    timestampTable = table(runNum, trialNum, sTrial, eTrial, sStim, eStim, sSacc, eSacc, sBlink, eBlink, 'VariableNames', varNames); % Initialize the table

    % initialize designMat, define the variable names for the table
    designvarNames = {'BaselineStaticMotion', 'staticDir', 'motionDir', 'cardinalGlobal'};
    designAllruns = nan(trialsPerRun*numel(runs), 4);

    % initialize alldat matrix - which will contain all the data data per
    % trial (trial x time=5000 (with a little buffer) x dim = x, y, pupil), where 
    allDat = nan(numTrials, 5500, 3); 
    
    % iterate across each run
    for ri=runs % run is the array 1:X

        path_to_edf = fullfile(designDir, dg_numID, sprintf('Run%i', ri), ...
            'eyedata');

         % Get a list of all .edf files in the directory
        edfFiles = dir(fullfile(path_to_edf, '*.edf'));

        % Extract the names of the .edf files
        edfFileNames = {edfFiles.name};

        if length(edfFileNames) > 1
            warning('More than 1 edf file for this run. Check..')
        else
            edfName = edfFileNames{1};
        end

        edfPath = fullfile(path_to_edf, edfName);

        msgfilePath = strrep(edfPath, '.edf', '.msg');
        datfilePath = strrep(edfPath, '.edf', '.dat');

        if ~isfile(msgfilePath) || ~isfile(datfilePath)
            convert_dgEDF(edfPath);
        end

        % extract meta data into a table, and append the metadata to the main table
        metaDataStruct = extractmsgMetadata(msgfilePath, subj);
        metaDataStruct.RunID = ri;
        metaData = struct2table(metaDataStruct);
        metaDataTable(ri, :) = metaData;

        % preprocess the run (correct for low freq drift and interpolate blinks)
        % where to index the dataTable 
        trialnumInit = (trialsPerRun*ri) - trialsPerRun; % the new "0" relative to this run.
        timestampTable = extractTimestamps(msgfilePath, timestampTable, trialnumInit, trialsPerRun);

        % get trial info from design file
        participant_file = struct2table(tdfread(fullfile(bidsDir,'participants.tsv')));
        rowInd = contains(cellstr(participant_file.participant_id), subj);
        dg_id = participant_file.alternate_id(rowInd,:);
        path_to_dmat = fullfile(designDir, dg_numID, sprintf('Run%i', ri), ...
                sprintf('%s_design_Run%i.mat', dg_id, ri));
        tmp = load(path_to_dmat, 'expDes');
        designAllruns(trialnumInit+1:trialnumInit+trialsPerRun,1:4) = tmp.expDes.trialMat(:,2:end);
        
        % just get params from 1st run; all runs were in the same session so params will be the same for each
        % e.g. params (pixels per degree; scaling of stimulus; offset of
        % center from display etc.
        if ri==1 
            path_to_const = fullfile(designDir, dg_numID, sprintf('Run%i', ri), ...
                sprintf('%s_const_file_Run%i*.mat', dg_id, ri));
            constFiles = dir(path_to_const); path_to_const = fullfile(constFiles.folder, constFiles.name);
            path_to_scr = fullfile(designDir, dg_numID, sprintf('Run%i', ri), ...
                sprintf('%s_scr_file_Run%i.mat', dg_id, ri));
            tmp = load(path_to_const); const = tmp.const;
            tmp = load(path_to_scr); scr = tmp.scr;
        end

        % get the .dat and put it all in a matrix (52 x Time)
        allDat = extractEyeposition(allDat, datfilePath, timestampTable, trialnumInit, trialsPerRun);
    end

    designTable = array2table(designAllruns, 'VariableNames', designvarNames); % Convert the array to a table
    
    % initialize columns to log microsaccades (all empty)
    msLabels = {'sMSacc', 'eMSacc'};
    emptyMSColumn = repmat({[]}, size(timestampTable,1), 2);
    msTable = array2table(emptyMSColumn, 'VariableNames', msLabels);

    manualCorrectionLabels = {'sManualCorr', 'eManualCorr'};
    emptyMCColumn = repmat({[]}, size(timestampTable,1), 2);
    mcTable = array2table(emptyMCColumn, 'VariableNames', manualCorrectionLabels);
        
    finalTable = [timestampTable, msTable, mcTable, designTable];

    % save the metadata table as a CSV
    filename = 'metaEyedata';
    save(fullfile(eyeanalysisFolder, filename), 'metaDataTable');

    % load in scr and const to get exact params
    % convert units to deg: screen is 1080 = y ; 1920 = X
    xCenter = scr.windCenter_px(1); %1920/2;
    yCenter = scr.windCenter_px(2); %1080/2;
    pxPerDeg = const.fixationRadius_px/const.fixationRadius_deg;

%     % NY subjects  were centered on display (no
%     % eyedata for 'sub-wlsubj121') in NY the scaling is perfect
%     if 1 || strcmp(subj, 'sub-wlsubj127') 
%         disp('no scaling')
%         pxPerDeg = 1080/(12.2*2);
%         xCenter = scr.windCenter_px(1); %1920/2;
%         yCenter = scr.windCenter_px(2); %1080/2;
%     else % using newer dataset for 'sub-0201'
%         maxDiam = 0.888; % percent scaled (this makes the stimulus smaller in pixel units)
%         yOffset_percent = 0.9136;
%         pxPerDeg = (1080*0.88)/(12.2*2); % in AD, I have to scale the stimulus b/c the projector takes up more than the bore screen
%         xCenter = 1920/2;
%         %if strcmp(subj, 'sub-wlsubj123')
%         yCenter = (1080/2)+((1-yOffset_percent)*1080/2); %(1080/2)*(yOffset_percent);
%         %else
%         %    yCenter=yOffset_percent*(1080/2);
%         %end
%     end
    
    
    % % 0-center the data
    allDat_dva = nan(size(allDat));
    allDat_dva(:,:,1) = allDat(:,:,1);
    allDat_dva(:,:,2) = (allDat(:,:,2)-xCenter)./pxPerDeg; % 2 = x
    allDat_dva(:,:,3) = (allDat(:,:,3)-yCenter)./pxPerDeg; % 3 = y
    
    % fill in data for msTable (uses same filter as in function below)
    finalTable = findMSs(allDat_dva, finalTable, msLabels, metaDataTable);

    % blinks into NaNs
    cushion = 300; % intrpolate 300 ms before and after blink
    prePeriod = cushion; postPeriod = cushion;
    allDat_dva = filter2nans(allDat_dva, finalTable, 'blinks', prePeriod, postPeriod);

    % eyelink saccades into NaNs
    cushion = 100; % intrpolate 300 ms before and after blink
    prePeriod = cushion; postPeriod = cushion;
    allDat_dva = filter2nans(allDat_dva, finalTable, 'saccades', prePeriod, postPeriod);

    % algorithm saccades into NaNs (includes microsaccades and saccades)
    cushion = 100; % intrpolate 300 ms before and after blink
    prePeriod = cushion; postPeriod = cushion;
    allDat_dva = filter2nans(allDat_dva, finalTable, 'microsaccades', prePeriod, postPeriod);

    % bandpass eye data, and remove small segments
    allDat_dva = smoothEyeData(allDat_dva, metaDataTable);

    allDat_dva = removeShortsegs(allDat_dva, 100); % min duration = 100 ms

    % organize each trial into row (total of 52xNruns rows, 5000 ms + nans for padding)
    eyedataname = 'eyedat_dva_removedblinksaccmsacc';
    save(fullfile(eyeanalysisFolder, eyedataname), 'allDat_dva');

    % save the timestamps table as a CSV
    timetable = 'timestampEyeata';
    save(fullfile(eyeanalysisFolder, timetable), 'finalTable');


    %%
    % filteredMatrix = allDat_dva;
    
    selectRuns = runs; 
    
    % find the latestfile
    datpath = fullfile(eyeanalysisFolder, 'eyedat_dva_removedblinksaccmsacc_manualcorr_*.mat');
    tabpath = fullfile(eyeanalysisFolder, 'timestampEyeData_wManualCorr_*.mat');
    datfiles = dir(datpath); tabfiles = dir(tabpath);
    if ~isempty(datfiles) && ~isempty(tabfiles)
        [~, idx] = max([datfiles.datenum]); % Get index of the most recent file
        datfilepath = fullfile(datfiles(idx).folder, datfiles(idx).name); % Get full path

        [~, idx] = max([tabfiles.datenum]); % Get index of the most recent file
        tabfilepath = fullfile(tabfiles(idx).folder, tabfiles(idx).name); % Get full path

        disp('READING IN MANUALLY CORRECTED INPUT...')
        tmp = load(datfilepath); allDat_dva = tmp.allDat_dva;
        tmp = load(tabfilepath); finalTable = tmp.finalTable;
    end

    % filter data based on horiz or vertical motion

    % rows of the MetaMatrix that have upward or downward motion
    rowIndicesV = intersect(find((finalTable.motionDir == 90 | finalTable.motionDir == 270)), find(ismember(finalTable.runNum, selectRuns))); %finalTable.runNum == selectRun);
    trialIDX_Vmotion = allDat_dva(rowIndicesV,1:timeCutOff,1);
    xPos_Vmotion = allDat_dva(rowIndicesV,1:timeCutOff,2);
    yPos_Vmotion = allDat_dva(rowIndicesV,1:timeCutOff,3);
    rowIndicesH = intersect(find((finalTable.motionDir == 0 | finalTable.motionDir == 180)), find(ismember(finalTable.runNum, selectRuns)));
    trialIDX_Hmotion = allDat_dva(rowIndicesH,1:timeCutOff,1);
    xPos_Hmotion = allDat_dva(rowIndicesH,1:timeCutOff,2);
    yPos_Hmotion = allDat_dva(rowIndicesH,1:timeCutOff,3);
        
    if ~(~isempty(datfiles) && ~isempty(tabfiles))
        % Quality Control (MANUAL + SAVE)
        
        % HORIZ MOTION there are 64 trials each (2 dirs x 4 repeats x 8 runs)
        
        output_hMotion = queryUserInput(xPos_Hmotion, yPos_Hmotion, trialIDX_Hmotion);
        finalTable.sManualCorr(rowIndicesH) = output_hMotion(:,1);
        finalTable.eManualCorr(rowIndicesH) = output_hMotion(:,2);
        
        % VERTICAL MOTION
        output_vMotion = queryUserInput(xPos_Vmotion, yPos_Vmotion, trialIDX_Vmotion);
        finalTable.sManualCorr(rowIndicesV) = output_vMotion(:,1);
        finalTable.eManualCorr(rowIndicesV) = output_vMotion(:,2);
        
        % save the timestamps table as a CSV
        nowstamp = datestr(now, 'yyyymmdd_HHMMSS');
        timetable = sprintf('timestampEyeData_wManualCorr_%s', nowstamp);
        save(fullfile(eyeanalysisFolder, timetable), 'finalTable');
        
        allDat_dva = filter2nans(allDat_dva, finalTable, 'manualcorrection', 0, 0);
        
        eyedataname = sprintf('eyedat_dva_removedblinksaccmsacc_manualcorr_%s', nowstamp);
        save(fullfile(eyeanalysisFolder, eyedataname), 'allDat_dva');
    end

    % rows of the MetaMatrix that have upward or downward motion
    %rowIndices = intersect(find((finalTable.motionDir == 90 | finalTable.motionDir == 270)), find(ismember(finalTable.runNum, selectRuns))); %finalTable.runNum == selectRun);
    %xPos_Vmotion = allDat_dva(rowIndices,:,2);
    %yPos_Vmotion = allDat_dva(rowIndices,:,3);
    
    figure
    subplot(4,1,1)
    plot(1:length(xPos_Vmotion), xPos_Vmotion);
    hold on
    xline(3000, 'k--', 'LineWidth', 2);
    hold on
    text(50, 2.5, 'stimulus ON', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 20, 'Color', 'k');
    text(3050, 2.5, 'stimulus OFF', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 20, 'Color', 'k');
    xlim([0, 5000])
    ylim([-4, 4])
    ax = gca;
    ax.XTick = 0:1000:5000;
    ax.YTick = -12:2:12;
    ax.YTickLabel = {'-12','', '-8', '', '-4', '', '0', '', '4', '', '8', '', '12'};
    ax.XTickLabel = {'0', '1', '2', '3', '4', '5'};
    ax.FontSize = 20; 
    xlabel('time (s)')
    ylabel('X distance (deg)')
    title('Vertical Motion - eye X component')
    subplot(4,1,2)
    plot(1:length(xPos_Vmotion), yPos_Vmotion);
    hold on
    xline(3000, 'k--', 'LineWidth', 2);
    xlim([0, 5000])
    ylim([-4, 4])
    ax = gca;
    ax.XTick = 0:1000:5000;
    ax.YTick = -12:2:12;
    ax.YTickLabel = {'-12','', '-8', '', '-4', '', '0', '', '4', '', '8', '', '12'};
    ax.XTickLabel = {'0', '1', '2', '3', '4', '5'};
    ax.FontSize = 20; 
    xlabel('time (s)')
    ylabel('X distance (deg)')
    title('Vertical Motion - eye Y component')
    
    %rowIndices = intersect(find((finalTable.motionDir == 0 | finalTable.motionDir == 180)), find(ismember(finalTable.runNum, selectRuns)));
    %xPos_Hmotion = allDat_dva(rowIndices,:,2);
    %yPos_Hmotion = allDat_dva(rowIndices,:,3);
    
    subplot(4,1,3)
    plot(1:length(xPos_Hmotion), xPos_Hmotion);
    hold on
    xline(3000, 'k--', 'LineWidth', 2);
    title('Horizontal Motion - eye X component')
    xlim([0, 5000])
    ylim([-4, 4])
    ax = gca;
    ax.XTick = 0:1000:5000;
    ax.YTick = -12:2:12;
    ax.YTickLabel = {'-12','', '-8', '', '-4', '', '0', '', '4', '', '8', '', '12'};
    ax.XTickLabel = {'0', '1', '2', '3', '4', '5'};
    ax.FontSize = 20; 
    ylabel('X distance (deg)')
    xlabel('time (s)')
    subplot(4,1,4)
    plot(1:length(xPos_Hmotion), yPos_Hmotion);
    hold on
    xline(3000, 'k--', 'LineWidth', 2);
    title('Horizontal Motion - eye Y component')
    xlim([0, 5000])
    ylim([-4, 4])
    ax = gca;
    ax.XTick = 0:1000:5000;
    ax.YTick = -12:2:12;
    ax.YTickLabel = {'-12','', '-8', '', '-4', '', '0', '', '4', '', '8', '', '12'};
    xlabel('time (s)')
    ax.XTickLabel = {'0', '1', '2', '3', '4', '5'};
    ax.FontSize = 20; 
    ylabel('Y distance (deg)')
    
    sgtitle(sprintf('%s X and Y gaze components (all runs)',subj), 'FontSize', 20)
    f1 = gcf;
    f1.Position = [155 123 1288 1191];
   
    
    
    %% individual polar plots
    
    % for circle
    radius = 1;
    center_x = 0;
    center_y = 0;
    theta = linspace(0, 2*pi, 100); % Define the range of theta from 0 to 2*pi
    % Calculate the x and y coordinates of the circle
    x = radius * cos(theta) + center_x;
    y = radius * sin(theta) + center_y;
    x2 = 12.2 * cos(theta) + center_x;
    y2 = 12.2 * sin(theta) + center_y;
    
    figure
    subplot(1,2,1)
    plot(x, y, 'b-', 'LineWidth', 1); % 'b-' specifies a blue line
    hold on
    plot(x2, y2, 'b-', 'LineWidth', 1); % 'b-' specifies a blue line
    hold on
    % Overlay semi-transparent scatter points
    scatter(xPos_Hmotion(:), yPos_Hmotion(:), 5, 'k', 'filled', 'MarkerFaceAlpha', 0.01, 'MarkerEdgeAlpha', 0.01);
    title('Horizontal motion', 'FontSize', 20)
    xlim([-12.2 12.2])
    ylim([-12.2, 12.2])
    hold on
    scatter(0,0,50, 'or', 'filled', 'LineWidth', 10)
    ylabel('deg', 'FontSize', 20)
    xlabel('deg', 'FontSize', 20)
    ax = gca; % Get the current axis
    ax.FontSize = 20; 
    
    subplot(1,2,2)
    plot(x, y, 'b-', 'LineWidth', 1); % 'b-' specifies a blue line
    hold on
    plot(x2, y2, 'b-', 'LineWidth', 1); % 'b-' specifies a blue line
    hold on
    % Overlay semi-transparent scatter points
    scatter(xPos_Vmotion(:), yPos_Vmotion(:), 5, 'k', 'filled', 'MarkerFaceAlpha', 0.01, 'MarkerEdgeAlpha', 0.01);
    title('Vertical motion', 'FontSize', 20)
    xlim([-12.2 12.2])
    ylim([-12.2, 12.2])
    hold on
    scatter(0,0,50, 'or', 'filled', 'LineWidth', 10)
    xlabel('deg', 'FontSize', 20)
    ax = gca; % Get the current axis
    ax.FontSize = 20; 
    
    f1 = gcf;
    f1.Position = [1328 861 1012 417];
    sgtitle(sprintf('gaze positions %s during stimulus ON (all runs)', subj), 'FontSize', 20)

end

%%

function allDat = filter2nans(allDat, finalTable, condition, prePeriod, postPeriod)

    if strcmp(condition, 'blinks')
        startlabel = 'sBlink';
        endlabel = 'eBlink';
    elseif strcmp(condition, 'saccades')
        startlabel = 'sSacc';
        endlabel = 'eSacc';
    elseif strcmp(condition, 'microsaccades')
        startlabel = 'sMSacc';
        endlabel = 'eMSacc';
    elseif strcmp(condition, 'manualcorrection')
        startlabel = 'sManualCorr';
        endlabel = 'eManualCorr';
    end

    for tr=1:size(allDat,1)
        currTrialT = allDat(tr,:,1);
        currTrialX = allDat(tr,:,2);
        currTrialY = allDat(tr,:,3);
        currTrialP = allDat(tr,:,4);
        
        % check 
        start_tmp = finalTable.(startlabel)(tr);
        ends_tmp = finalTable.(endlabel)(tr);
        
        condStarts = start_tmp{1};
        condEnds = ends_tmp{1};
        
        if any((condEnds - condStarts)<0)
            warning('special handing: %s occured between trials.', condition)
        else
            % go through each event
            for bi=1:max([numel(condStarts),numel(condEnds)]) 
        
              indexS = find(currTrialT == condStarts(bi));
              indexS = indexS - prePeriod;
              
              indexE = find(currTrialT == condEnds(bi));
              indexE = indexE + postPeriod;
            
                % special handling needed regardless due to cushion
                if indexS <= 0
                    indexS = 1;
                end
                if indexE > length(currTrialT)
                    indexE = length(currTrialT);
                end
    
                %allDat(tr,indexS:indexE,1) = NaN; % leave the idx
                allDat(tr,indexS:indexE,2) = NaN; 
                allDat(tr, indexS:indexE,3) = NaN; 
            end
        end
        
    end
end

%%

function allDat_dva = smoothEyeData(allDat_dva, metaDataTable)
    % Parameters
    numTrials = size(allDat_dva, 1); % Total trials
    numTimepoints = 5000; % Cutoff at 5000 ms
    trialsPerRun = 52; % Number of trials per run

    % get sampling rate (not needed if specificying normalized freq)
    numGroups = ceil(numTrials / trialsPerRun);
    fs = repelem(metaDataTable.SamplingRate(1:numGroups), trialsPerRun, 1);

    %p_cutoff = 0.01; % Low-pass cutoff frequency (normalized)

    % FIR filter uses 36 coefficients (order +1)
    %filter_order = 35; % FIR filter order

    % Design low-pass filter to remove high-frequency drift
    %hp_filter = fir1(filter_order, p_cutoff); 
    
    % buterrworth filter
    [b, a] = butter(2, 30 / (fs(1)/2), 'low');

    trialRange = 1:numTrials;

    % Loop over each run
    %for runIdx = 1:floor(numTrials / trialsPerRun)
        % Get trial indices for this run
    %    trialRange = (runIdx-1) * trialsPerRun + (1:trialsPerRun);

    % Extract X, Y data (reshape into a single time series across trials)
    x_data = reshape(allDat_dva(trialRange, 1:numTimepoints, 2)', [], 1);
    y_data = reshape(allDat_dva(trialRange, 1:numTimepoints, 3)', [], 1);

    % Preserve NaN indices
    nanMaskX = isnan(x_data);
    nanMaskY = isnan(y_data);

    % Fill missing data with linear interpolation
    x_data_filled = fillmissing(x_data, 'linear');
    y_data_filled = fillmissing(y_data, 'linear');

    % Apply mirror padding to avoid edge artifacts
    %x_padded = [flip(x_data_filled(1:filter_order*2)); x_data_filled; flip(x_data_filled(end-filter_order*2+1:end))];
    %y_padded = [flip(y_data_filled(1:filter_order*2)); y_data_filled; flip(y_data_filled(end-filter_order*2+1:end))];

    % butter
    x_padded = [flip(x_data_filled(1:100)); x_data_filled; flip(x_data_filled(end-99:end))];
    y_padded = [flip(y_data_filled(1:100)); y_data_filled; flip(y_data_filled(end-99:end))];


    % Apply high-pass filtering across the entire run
    %x_filtered = filtfilt(hp_filter, 1, x_padded);
    %y_filtered = filtfilt(hp_filter, 1, y_padded);

    %butter
    x_filtered = filtfilt(b, a, x_padded);
    y_filtered = filtfilt(b, a, y_padded);

    % Remove padding
    %x_filtered = x_filtered(filter_order*2+1:end-filter_order*2);
    %y_filtered = y_filtered(filter_order*2+1:end-filter_order*2);
    
    % butter
    x_filtered = x_filtered(101:end-100);
    y_filtered = y_filtered(101:end-100);

    % Restore NaNs
    x_filtered(nanMaskX) = NaN;
    y_filtered(nanMaskY) = NaN;

    % Reshape back into (trials × time) format
    x_filtered = reshape(x_filtered, numTimepoints, numTrials)';
    y_filtered = reshape(y_filtered, numTimepoints, numTrials)';

    % Insert back into allDat_dva
    allDat_dva(trialRange, 1:numTimepoints, 2) = x_filtered;
    allDat_dva(trialRange, 1:numTimepoints, 3) = y_filtered;
    %end
end

% % example trials
% trialID = 300;
% tmpX = allDat_dva(trialID,:,2);
% tmpY = allDat_dva(trialID,:,3);
% figure
% plot(tmpX, 'LineWidth',2)
% hold on
% plot(tmpY, 'LineWidth',2)
% legend({'xPos'; 'yPos'});
% %ylim([-6 6])
% xlim([0 5000])
% ax = gca;
% ax.XTick = 0:1000:5000;
% ax.XTickLabel = {'0', '1', '2', '3', '4', '5'};
% xlabel('time (s)')
% 
% % filter the data to remove eye tracker noise
% Fs = 1000; % Sampling frequency (Hz)
% cutoffFreq = 30; % 60 before % Cutoff frequency (Hz)
% [b, a] = butter(4, cutoffFreq/(Fs/2), 'low'); % 4th-order Butterworth filter
% 
% % Initialize the filtered matrix
% [nTrials, nTimestamps, nDim] = size(allDat_dva);
% 
% % %Apply the filter to each row
% % filteredMatrix = zeros(nTrials, nTimestamps, nDim);
% % for di = 1:nDim % dimensions x, y
% %     for i = 1:nTrials
% %         
% %         timSeg = allDat_dva(i, :, di);
% %         
% %         nonNaNIdx = find(~isnan(timSeg));
% %         if length(nonNaNIdx)>12
% %             % Extract non-NaN part
% %             nonNaNData = timSeg(nonNaNIdx);
% % 
% %             % Apply the filter to the non-NaN part
% %             filteredNonNaNData = filtfilt(b, a, nonNaNData);
% % 
% %             % Place the filtered non-NaN part back into the row
% %             timSeg(nonNaNIdx) = filteredNonNaNData;
% %             
% %              % Store the filtered row
% %             filteredMatrix(i, :, di) = timSeg;
% %         else
% %             disp('filtering out data that is too short.')
% %             filteredMatrix(i, :, di) = nan(size(timSeg));
% %         end
% % 
% %     end
% % end
