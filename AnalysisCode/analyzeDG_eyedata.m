%% analyze drifting grating eye data

% set up

clc; clear all; close all;

% set up
addpath(genpath(pwd));
projectName = 'dg';
%bidsDir =  '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/';
bidsDir = '/Volumes/server/Projects/Project_dg/data_bids';
githubDir = '~/Documents/GitHub';
fsDir = '/Applications/freesurfer'; %/7.2.0';
designDir = '/Volumes/server/Projects/Project_dg/experimentalOutput/dg/';
%addpath(genpath(fullfile(githubDir, 'retinotopy-nyuad')));
addpath(genpath(fullfile(githubDir, 'wpToolbox')));
addpath(genpath(fullfile(githubDir, 'atlasmgz')));
%setup_user(projectName,bidsDir,githubDir,fsDir);
setup_user('sawtooth',bidsDir);

tr_s = 1;
trialsPerRun = 52;

% for metadata
headings = {'SamplingRate', 'CalibrationMethod', 'Eye', 'CalibrationQuality', ...
            'ValidationQuality', 'AverageError', 'MaxError', 'OffsetDeg', 'Offpix', ...
            'NoTrackingTrials', 'Filename', 'SubjectName','RunID'};
types = {'double', 'string', 'string', 'string', 'string', 'double', 'double', ...
         'double', 'string', 'double', 'string', 'string','double'};

% for timestamps table
varNames = {'runNum', 'trialNum','sTrial', 'eTrial', 'sStim', 'eStim', 'sSacc', 'eSacc', 'sBlink', 'eBlink'}; % Define the variable names and types
    
%%

subjectslist = {'sub-wlsubj123'; 'sub-wlsubj124'; 'sub-0426'; 'sub-0427'; 'sub-0395';
    'sub-0037'; 'sub-0397'; 'sub-0255'; 'sub-0201'; 'sub-0442'; 'sub-wlsubj127'; 'sub-0250'};

% check sub-0397 & sub-0201 eye data time dimension (is 20,072?) why? -->
% 2,000 sampling rate.

for si=12:12 %length(subjectslist)
    
    subj = subjectslist{si}; 
    
    if strcmp(subj, 'sub-wlsubj121')
        continue % NOTE: no eye data for subj121
    end
    
    % retrieve design matrices for a subject (returns 1x8 cell) that has 280
    % rows
    [matrices_init, dg_id, dg_numID, runs] = format_desmats(bidsDir, designDir, subj, tr_s);

    totalRuns = numel(runs);
    % conversion the edfs for a subject (returns 1x8 cell) for dat.

%     % retrieve alternate ID
%     participant_file = struct2table(tdfread(fullfile(bidsDir,'participants.tsv')));
%     rowInd = contains(cellstr(participant_file.participant_id), subj);
%     dg_id = participant_file.alternate_id(rowInd,:);
%     dg_numID = regexp(dg_id, '\d+', 'match'); dg_numID = dg_numID{1};
% 
%     % compute how many runs exist for this subject
%     logFile = fullfile(designDir, dg_numID, 'runlog.txt');
%     fid = fopen(logFile, 'r');
%     runs = fscanf(fid, '%s');
%     runarray = split(runs, 'Run');
%     totalRuns = str2num(runarray{end});
%     runs = 1:totalRuns;

    eyeanalysisFolder = fullfile(bidsDir, 'derivatives', 'eyetracking', subj, 'ses-dg');
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
        
        % get the .dat and put it all in a matrix (52 x Time)
        allDat = extractEyeposition(allDat, datfilePath, timestampTable, trialnumInit, trialsPerRun);
    end

    designTable = array2table(designAllruns, 'VariableNames', designvarNames); % Convert the array to a table
    finalTable = [timestampTable, designTable];

    % save the metadata table as a CSV
    filename = 'metaEyedata.csv';
    writetable(metaDataTable, fullfile(eyeanalysisFolder, filename));

    % save the timestamps table as a CSV
    timetable = 'timestampEyedata.csv';
    writetable(finalTable, fullfile(eyeanalysisFolder, timetable));
    
    % organize each trial into row (total of 280xNruns rows)
    eyedataname = 'eyedatOrganized.csv';
    writematrix(allDat, fullfile(eyeanalysisFolder, eyedataname));

end

%%



% convert units to deg
% screen is 1080 = y ; 1920 = X
if strcmp(subj, 'sub-wlsubj127') % no eyedata for 'sub-wlsubj121'
    pxPerDeg = 1080/(12.2*2); % in NY the scaling is perfect
    xCenter = 1920/2;
    yCenter = 1080/2;
else
    maxDiam = 0.888; % percent scaled (this makes the stimulus smaller in pixel units)
    yOffset_percent = 0.9136;
    pxPerDeg = (1080*0.88)/(12.2*2); % in AD, I have to scale the stimulus b/c the projector takes up more than the bore screen
    xCenter = 1920/2;
    if strcmp(subj, 'sub-wlsubj123')
        yCenter = (1080/2)+((1-yOffset_percent)*1080/2); %(1080/2)*(yOffset_percent);
    else
        yCenter=yOffset_percent*(1080/2);
    end
end



% blinks into NaNs (later interpolate blinks)

cushion = 300; % intrpolate 150 ms before and after blink

for tr=1:size(allDat,1)
    currTrialT = allDat(tr,:,1);
    currTrialX = allDat(tr,:,2);
    currTrialY = allDat(tr,:,3);
    currTrialP = allDat(tr,:,4);
    
    % check 
    blinkStart_tmp = finalTable.sBlink(tr);
    blinkEnds_tmp = finalTable.eBlink(tr);
    
    blinkStarts = blinkStart_tmp{1};
    blinkEnds = blinkEnds_tmp{1};
    
    if any((blinkEnds - blinkStarts)<0)
        warning('special handing: blink occured between trials.')
    else
        blinkStarts
        blinkEnds
        % go through each blink
        for bi=1:max([numel(blinkStarts),numel(blinkEnds)]) 
    
          indexS = find(currTrialT == blinkStarts(bi));
          indexS = indexS - cushion;
          
          indexE = find(currTrialT == blinkStarts(bi));
          indexE = indexE + cushion;
        
            % special handling needed regardless due to cushion
            if indexS <= 0
                indexS = 1;
            end
            if indexE > length(currTrialT)
                indexE = length(currTrialT);
            end

            allDat(tr,indexS:indexE,1) = NaN; 
            allDat(tr,indexS:indexE,2) = NaN; 
            allDat(tr, indexS:indexE,3) = NaN; 
        end
    end
    
end
%%
%
% % 0-center the data
allDat_dva = nan(size(allDat));
allDat_dva(:,:,1) = (allDat(:,:,2)-xCenter)./pxPerDeg; % 2 = x
allDat_dva(:,:,2) = (allDat(:,:,3)-yCenter)./pxPerDeg; % 3 = y

% example trials
trialID = 50;
tmpX = allDat_dva(trialID,:,1);
tmpY = allDat_dva(trialID,:,2);
figure
plot(tmpX, 'LineWidth',2)
hold on
plot(tmpY, 'LineWidth',2)
legend({'xPos'; 'yPos'});
%ylim([-6 6])
xlim([0 5000])
ax = gca;
ax.XTick = 0:1000:5000;
ax.XTickLabel = {'0', '1', '2', '3', '4', '5'};
xlabel('time (s)')


%

% filter the data to remove eye tracker noise
Fs = 1000; % Sampling frequency (Hz)
cutoffFreq = 30; % 60 before % Cutoff frequency (Hz)
[b, a] = butter(4, cutoffFreq/(Fs/2), 'low'); % 4th-order Butterworth filter

% Initialize the filtered matrix
[nTrials, nTimestamps, nDim] = size(allDat_dva);

% %Apply the filter to each row
% filteredMatrix = zeros(nTrials, nTimestamps, nDim);
% for di = 1:nDim % dimensions x, y
%     for i = 1:nTrials
%         
%         timSeg = allDat_dva(i, :, di);
%         
%         nonNaNIdx = find(~isnan(timSeg));
%         if length(nonNaNIdx)>12
%             % Extract non-NaN part
%             nonNaNData = timSeg(nonNaNIdx);
% 
%             % Apply the filter to the non-NaN part
%             filteredNonNaNData = filtfilt(b, a, nonNaNData);
% 
%             % Place the filtered non-NaN part back into the row
%             timSeg(nonNaNIdx) = filteredNonNaNData;
%             
%              % Store the filtered row
%             filteredMatrix(i, :, di) = timSeg;
%         else
%             disp('filtering out data that is too short.')
%             filteredMatrix(i, :, di) = nan(size(timSeg));
%         end
% 
%     end
% end

filteredMatrix = allDat_dva;

%%

selectRuns = 1:totalRuns; 

% rows of the MetaMatrix that have upward or downward motion
rowIndices = intersect(find((finalTable.motionDir == 90 | finalTable.motionDir == 270)), find(ismember(finalTable.runNum, selectRuns))); %finalTable.runNum == selectRun);
xPos_Vmotion = filteredMatrix(rowIndices,:,1);
yPos_Vmotion = filteredMatrix(rowIndices,:,2);

figure
subplot(4,1,1)
plot(1:5500, xPos_Vmotion);
hold on
xline(3000, 'k--', 'LineWidth', 2);
hold on
text(50, 6, 'stimulus ON', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 20, 'Color', 'k');
text(3050, 6, 'stimulus OFF', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 20, 'Color', 'k');
xlim([0, 5000])
ylim([-12.2, 12.2])
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
plot(1:5500, yPos_Vmotion);
hold on
xline(3000, 'k--', 'LineWidth', 2);
xlim([0, 5000])
ylim([-12.2, 12.2])
ax = gca;
ax.XTick = 0:1000:5000;
ax.YTick = -12:2:12;
ax.YTickLabel = {'-12','', '-8', '', '-4', '', '0', '', '4', '', '8', '', '12'};
ax.XTickLabel = {'0', '1', '2', '3', '4', '5'};
ax.FontSize = 20; 
xlabel('time (s)')
ylabel('X distance (deg)')
title('Vertical Motion - eye Y component')

rowIndices = intersect(find((finalTable.motionDir == 0 | finalTable.motionDir == 180)), find(ismember(finalTable.runNum, selectRuns)));
xPos_Hmotion = filteredMatrix(rowIndices,:,1);
yPos_Hmotion = filteredMatrix(rowIndices,:,2);

subplot(4,1,3)
plot(1:5500, xPos_Hmotion);
hold on
xline(3000, 'k--', 'LineWidth', 2);
title('Horizontal Motion - eye X component')
xlim([0, 5000])
ylim([-12.2, 12.2])
ax = gca;
ax.XTick = 0:1000:5000;
ax.YTick = -12:2:12;
ax.YTickLabel = {'-12','', '-8', '', '-4', '', '0', '', '4', '', '8', '', '12'};
ax.XTickLabel = {'0', '1', '2', '3', '4', '5'};
ax.FontSize = 20; 
ylabel('X distance (deg)')
xlabel('time (s)')
subplot(4,1,4)
plot(1:5500, yPos_Hmotion);
hold on
xline(3000, 'k--', 'LineWidth', 2);
title('Horizontal Motion - eye Y component')
xlim([0, 5000])
ylim([-12.2, 12.2])
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

selectRuns = 1:8;
timeCutOff = 3000;

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

% rows of the MetaMatrix that have upward or downward motion
rowIndices = intersect(find((finalTable.motionDir == 90 | finalTable.motionDir == 270)), find(ismember(finalTable.runNum, selectRuns))); %finalTable.runNum == selectRun);
xPos_Vmotion = filteredMatrix(rowIndices,1:timeCutOff,1);
yPos_Vmotion = filteredMatrix(rowIndices,1:timeCutOff,2);
rowIndices = intersect(find((finalTable.motionDir == 0 | finalTable.motionDir == 180)), find(ismember(finalTable.runNum, selectRuns)));
xPos_Hmotion = filteredMatrix(rowIndices,1:timeCutOff,1);
yPos_Hmotion = filteredMatrix(rowIndices,1:timeCutOff,2);

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

%%

% A = [1,2,3; 10, 11, 12]; % Example matrix, replace with your actual matrix
% 
% % Transpose the matrix to switch rows and columns
% A_transposed = A';
% 
% % Reshape the transposed matrix into a single row
% B = reshape(A_transposed, 1, []);
