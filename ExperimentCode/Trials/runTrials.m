function [task, trial_onsets] = runTrials(scr,const,expDes,my_key,textExp)
% ----------------------------------------------------------------------
% runTrials(scr,const,expDes,my_key,textExp,button)
% ----------------------------------------------------------------------
% Goal of the function :
% Main trial function, display the trial function and save the experi-
% -mental data in different files.
% ----------------------------------------------------------------------
% Input(s) :
% scr : window pointer
% const : struct containing all the constant configurations.
% expDes : struct containing all the variable design and configurations.
% my_key : keyboard keys names
% textExp : struct contanining all instruction text.
% ----------------------------------------------------------------------
% Output(s):
% none
% ----------------------------------------------------------------------

%% General instructions:

disp('Starting runTrials')

% eyetracking (to do: need to add vpixx)
if const.EL_mode
    if (const.run==1) || rem(const.run,4) == 0
        % removed lines below since redundant
        %keyCode = instructions(scr,const,my_key,textExp.eyeinstruction);
        %if keyCode(my_key.escape), Eyelink('Shutdown'); vpixxShutdown(const); return; end
        %FlushEvents('KeyDown');
        % expects T input to start (change this later)
        [~, exitFlag] = initEyelinkStates('calibrate', const.window, const.EL);
        if exitFlag, vpixxShutdown(const); return, end
    end
    err = Eyelink('CheckRecording');
    if err ~= 0
        initEyelinkStates('startrecording', const.window, const.EL);
        disp('Eyelink now recording .. ')
    end
end

HideCursor(scr.scr_num);
keyCode = instructions(scr,const,my_key,textExp.instruction);

task = my_task(expDes, scr);

if keyCode(my_key.escape), vpixxShutdown(const); return, end

FlushEvents('KeyDown');

%% Main Loop
frameCounter=1;
const.expStop = 0;
paddingX = 0;

tic
vbl = Screen('Flip',const.window);
t0=vbl;
trial_onsets = nan(1,(expDes.nb_trials));

while ~const.expStop
    
    % PADDING PERIOD
    [paddingX, task, frameCounter, vbl] = my_padding(my_key, scr,const, expDes, task, frameCounter, paddingX, vbl);
    
    % TRIAL PERIOD
    if paddingX < 2
        for ni=1:expDes.nb_trials
            if const.EL_mode, Eyelink('message', 'TRIAL %d START', ni); end
            trial_onsets(ni) = vbl-t0; % log the onset of each
            [task, frameCounter, vbl] = my_stim(my_key, scr,const,expDes,task, frameCounter,ni,vbl);
            if const.EL_mode, Eyelink('message', 'TRIAL %d END', ni); end
        end
    end
  
    if paddingX >= 2
        const.expStop = 1;
    end
    
end
toc

vpixxShutdown(const);

% save eyetracking file
if const.EL_mode
    disp("Please wait, saving EYELINK file..")
    if ~exist(const.eyeDataDir, 'dir'), mkdir(const.eyeDataDir); end
    initEyelinkStates('eyestop', scr.scr_num, {const.eyeFileName, const.eyeDataDir})  
end


end