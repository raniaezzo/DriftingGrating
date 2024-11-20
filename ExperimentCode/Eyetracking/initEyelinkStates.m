function [output, exitFlag] = initEyelinkStates(command, window, input)

% adapted from rd_eyeLink

% Possible commands and their ins & outs:
% 1. 'eyestart'
%   in = eyeFile
%   out = el
%
% 2. 'calibrate'
%   in = el
%   out = cal
%
% 3. 'trialstart'
%   in = {el, trialNum, cx, cy, rad}
%   out = []
%
% 4. 'fixholdcheck'
%   in = {cx, cy, rad}
%   out = fixation
%
% 5. 'fixcheck'
%   in = {cx, cy, rad}
%   out = fixation
%
% 6. 'driftcorrect'
%   in = {el, cx, cy}
%   out = driftCorrection
%
% 7. 'trialstop'
%   in = []
%   out = []
%
% 8. 'eyestop'
%   in = {eyeFile, eyeDataDir}
%   out = []

%% Initializations
% assume no output unless some is given
output = [];

% assume everything goes ok (exitFlag=0) until proven otherwise        
exitFlag = 0;
        
%% take action based on command
switch command
    case 'eyestart'
        %% start eyetracker
        eyeFile = input{1};
        screen = input{2};
        
        % First check if we can get a connection
        if EyelinkInit() ~= 1
            fprintf('\nCouldn''t initialize connection with eyetracker! Exiting ...\n');
            return
        end
        
        % Set up the eyetracker
        EL = initEyelinkDefaultVars(window, screen); % do not unify the keyboard!!
        
        Eyelink('Command', 'file_sample_data = LEFT,RIGHT,GAZE,AREA');
        Eyelink('Command', 'calibration_type = HV5');
        %Eyelink('Command', 'enable_automatic_calibration = YES');
        Eyelink('Command', 'calibration_area_proportion 0.5 0.5'); % 0.83 0.88
        Eyelink('Command', 'validation_area_proportion 0.5 0.5'); % 0.83 0.88
        
        Eyelink('Command', 'sample_rate = 1000');
        
        [~, vs] = Eyelink('GetTrackerVersion');
        fprintf('\nRunning experiment on a %s tracker.\n', vs );
        
        % Start the eye file
        edfFile = sprintf('%s.edf', eyeFile);
        edfFileStatus = Eyelink('OpenFile', edfFile);
        if edfFileStatus == 0
            fprintf('\nEye file opened in:\n\n')
            fprintf('%s\n\n', pwd)
        else
            fprintf('\nCannot open eye file (check if the eyelink is really shut down!).\n');
            Screen('CloseAll')
            Eyelink( 'Shutdown');
            exitFlag = 1;
            return
        end
        
        output = EL; % return the el structure as output
        
    case 'calibrate'
        %% calibrate eyetracker
        EL = input;
        
        cali_string = sprintf('Eye tracker calibration:\n\nPlease fixate on the markers.');
        %\n\nPress ''space'' to start or ''q'' to quit');
        DrawFormattedText(window, cali_string, 'center', 'center', 1, []);
        Screen('Flip', window, 0, 1); 
        
        pause(2)

%         % eyetracker instructions until space/q is pressed
%         contKey = '';
%         while isempty(find(strcmp(contKey,'space'), 1))
%             keyIsDown = 0;
%             while ~keyIsDown
%                 [keyIsDown, ~, keyCode] = KbCheck(-1); %% listen to all keyboards
%                 disp(find(keyCode))
%             end
%             contKey = KbName(find(keyCode));
%         end
%         
%         % clear keyBoard events
%         while KbCheck; end
%         FlushEvents('KeyDown');
%         clear KbCheck;
%         
%         if strcmp(contKey,'q')
%             ListenChar(0);
%             ShowCursor;
%             Screen('CloseAll')
%             fclose('all');
%             fprintf('User ended program');
%             exitFlag = 1;
%             return
%         end
%         Screen('Flip', window, 0, 1);
        
        cali = EyelinkDoTrackerSetup(EL);    % this fails --
        if cali == EL.TERMINATE_KEY, exitFlag = 1;return, end
        
        output = cali;
        
    case 'startrecording'
        EL = input;
        
        record = 0;
        while ~record
            Eyelink('StartRecording');	% start recording
            % start recording 100 msec before just to be safe
            WaitSecs(0.1);
            key=1;
            while key~=0, key = EyelinkGetKey(EL); end % dump any pending local keys
            
            
            err = Eyelink('CheckRecording'); 	% check recording status
            if err == 0, record = 1; Eyelink('Message', 'RECORD_START');
            else, record = 0; Eyelink('Message', 'RECORD_FAILURE');% results in repetition of fixation check
            end
        end
        
    case 'trialstart'
        %% trial start
        % start only when we are recording and the subject is fixating
        % initEyelinkStates('trialstart', window, {EL, run.itrial, screen.centerX, screen.centerX, screen.rad})
        EL = input{1};
        itrial = input{2};
        cx = input{3};
        cy = input{4};
        rad = input{5};
        
        driftCorrected = 0;
        
        % Displays a title at the bottom of the eye tracker display
        % Start the trial only when it is 'ready': (1) eyetracker is recording, (2) subject is fixating
%         WaitSecs(.1);
        ready = 0; 
        while ~ready
            % Check that we are recording
            err = Eyelink('CheckRecording'); % report 0 is recording in progress
            if err ~= 0, initEyelinkStates('startrecording', window, EL); end
            
            % Verify that the subject is holding fixation for some set
            % time before allowing the trial to start. A
            % timeout period is built into this function.
            fixation = initEyelinkStates('fixholdcheck', window, {cx, cy, rad});
            
            % Drift correct if fixation timed out

            if ~fixation
                initEyelinkStates('driftcorrect', window, {EL, cx, cy});
                driftCorrected = 1;
                ready = 0;
            else
                ready = 1;
            end
        end
        
        output = driftCorrected;
        
        Eyelink('Message', 'TRIAL_START %d', itrial);
        Eyelink('Message', 'SYNCTIME');		% zero-plot time for EDFVIEW
        
    case 'fixholdcheck'
        %% check that fixation is held for some amount of time
        cx = input{1}; % x coordinate of screen center
        cy = input{2}; % y coordinate of screen center
        rad = input{3}; % acceptable fixation radius %%% in px?
        
        timeout = 3.00; % 3.00 % maximum fixation check time
        tFixMin = 0.30; % 0.10 % minimum correct fixation time
        
        Eyelink('Message', 'FIX_HOLD_CHECK');
        
        tstart = GetSecs;
        fixation = 0; % is the subject fixating now?
        fixStart = 0; % has a fixation already started?
        tFix = 0; % how long has the current fixation lasted so far?
        
        t = tstart;
        while (((t-tstart) < timeout) && (tFix <= tFixMin))
            % get eye position
            evt = Eyelink('newestfloatsample');
            domEye = find(evt.gx ~= -32768);
            if numel(domEye)>1 , domEye = domEye(1); end % if tracking binocularly
            x = evt.gx(domEye);
            y = evt.gy(domEye);

            % check for blink
            if isempty(x) || isempty(y)
                fixation = 0;
            else % check eye position
                if sqrt((x-cx)^2+(y-cy)^2)<rad, fixation = 1;else, fixation = 0;end
            end
            
            % update duration of current fixation
            if fixation==1 && fixStart==0
                tFix = 0;
                tFixStart = GetSecs;
                fixStart = 1;
            elseif fixation==1 && fixStart==1
                tFix = GetSecs-tFixStart;
            else
                tFix = 0;
                fixStart = 0;
            end
            
            t = GetSecs;
        end
        
        output = fixation;
        
    case 'fixcheck'
        %% check fixation at one point in time
        cx = input{1}; % x coordinate of screen center
        cy = input{2}; % y coordinate of screen center
        rad = input{3}; % acceptable fixation radius %%% in px?
        
        % determine recorded eye
        evt = Eyelink('newestfloatsample');
        domEye = find(evt.gx ~= -32768);
        
        % if tracking binocularly, just select one eye to be dominant
        if numel(domEye)>1, domEye = domEye(1); end
        
        Eyelink('Message', 'FIX_CHECK');
        
        % get eye position
        x = evt.gx(domEye);
        y = evt.gy(domEye);
        
        % check for blink
        if isempty(x) || isempty(y), fixation = 0;
        else % check eye position
            if sqrt((x-cx)^2+(y-cy)^2)<rad, fixation = 1; else, fixation = 0;end
        end
        
        if fixation==0, Eyelink('Message', sprintf('BROKE_FIXATION')); end
        
        output = fixation;
        
    case 'driftcorrect'
        %% do a drift correction
        EL = input{1};
        cx = input{2};
        cy = input{3};
        
        Eyelink('Message', 'DRIFT_CORRECTION');
        driftCorrection = EyelinkDoDriftCorrect(EL, cx, cy, 1, 1);
        
        output = driftCorrection;
        
    case 'stoprecording'
        %% stop recording
        Eyelink('StopRecording');
        Eyelink('Message','RECORD_STOP');
        
    case 'eyestop'
        %% get the eye file and close down the eye tracker
        eyeFile = input{1};
        eyeDataDir = input{2};
        
        % if still recording, stop recording
        err = Eyelink('CheckRecording');
        if err==0, initEyelinkStates('stoprecording'); end
        
        fprintf('\n\nSaving file %s/%s ...\n', eyeDataDir, eyeFile)
        
        Eyelink('ReceiveFile', eyeFile, eyeDataDir, 1); 
        Eyelink('CloseFile'); 
        Eyelink('Shutdown');
        
    otherwise
        error('[initEyelinkStates]: ''command'' argument not recognized. See help for available commands.')
end
