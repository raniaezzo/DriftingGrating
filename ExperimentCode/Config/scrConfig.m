function [scr, const]=scrConfig(const)
% ----------------------------------------------------------------------
% [scr]=scrConfig(const)
% ----------------------------------------------------------------------
% Goal of the function :
% Give all information about the screen and the monitor.
% ----------------------------------------------------------------------
% Input(s) :
% const : struct containing initial constant params.
% ----------------------------------------------------------------------
% Output(s):
% scr : struct containing all screen configuration.
% ----------------------------------------------------------------------

%% General

% Instructions
const.text_size = 20;
const.text_font = 'Helvetica';

% Color Configuration :
%const.red =     [200,   0,   0];
%const.green =   [  0, 200,   0];
%const.blue =    [  0,   0, 200];
%const.gray =    [127, 127, 127];
%const.white =   [255, 255, 255]; 
%const.darkgray =[ 64,  64,  64];
const.white =    [ 1,  1,  1];
const.gray =    [ .5,  .5,  .5];
const.black =   [  0,   0,   0];
const.lightgray= [ .8, .8,  .8];

% Time
const.my_clock_ini = clock;

%% Screen
computerDetails = Screen('Computer');  % check computer specs

scr.scr_num = max(Screen('Screens')); % use max screen
disp('Screen Number')
scr.scr_num

scr.maxDiam_percent = 1; % maximum % of radius without eyetracker obstruction

const.vpixx = 0; % assume no vpixx unless set below (overridden by debug==0)

% Size of the display (mm) - will be overwritten if provided below:
[scrX_mm, scrY_mm] = Screen('DisplaySize',scr.scr_num);

yOffset_percent = 1;

% find screen details
if ~computerDetails.windows
    switch computerDetails.localHostName
        case 'Ranias-MacBook-Pro-2'
            scr.experimenter = 'Rania';
            scr.scrViewingDist_cm = 50;
            const.keyboard = 'Apple Internal Keyboard / Trackpad';
        case 'ADUAE08550LP-MX-4'
            scr.experimenter = 'NYUADScanner';
            scrX_mm = 720; scrY_mm = 405;
            scr.scrViewingDist_cm = 83.5; %88;
            disp('ACCOUNTING FOR SCALING OFFSET BETWEEN PROJECTED IMAGE AND DISPLAY')
            scr.maxDiam_percent = 0.89; 
            yOffset_percent = 0.914;
            const.keyboard = 'Magic Keyboard with Numeric Keypad';
            if const.DEBUG == 1
                const.vpixx = 0;
            else
                const.vpixx = 1;
            end
        case 'Stimulus-Mac-2'
            scr.experimenter = 'NYUNYScanner';
            scr.scrViewingDist_cm = 83.5;
            %(check projected image measurements in NYU NY
            %scrX_mm = 600; scrY_mm = 362; % based on documentation 
            % readout is 677 x 380 when automated (keeping this b/c I think there is 
            % a small part of the projected image off screen).
            const.vpixx = 0;
            if const.DEBUG == 1
                const.keyboard = 'Magic Keyboard'; % local experimental mac
            else
                const.keyboard = '932'; % buttonbox / trigger (different device)
            end
        otherwise
            disp('SET screen configuration & viewing distance in scrConfig.m')
            scr.experimenter = 'Unknown';
            scr.scrViewingDist_cm = 57;
            const.DEBUG = 1; 
    end
else    % PC (field names are different)
    switch computerDetails.system
        case 'NT-10.0.9200 - '
            scr.experimenter = 'NYUADMockScanner';
            scr.scrViewingDist_cm = 88;
            scr.maxDiam_percent = 1; %15.82/18.56;
            const.vpixx = 0;
            scr.scr_num = 1;
        otherwise
            disp('SET screen configuration & viewing distance in scrConfig.m')
            scr.experimenter = 'Unknown';
            scr.scrViewingDist_cm = 57;
            const.DEBUG = 1; 
    end
end


scr.scrX_cm = scrX_mm/10; scr.scrY_cm = scrY_mm/10;

% Resolution of the display (pixels):
resolution = Screen('Resolution',scr.scr_num);
scr.scrX_px = resolution.width;
scr.scrY_px = resolution.height;
scr.scrPixelDepth_bpp = resolution.pixelSize; % bits per pixel

[scr.windX_px, scr.windY_px]=Screen('WindowSize', scr.scr_num);

if const.DEBUG == 1
    %PsychDebugWindowConfiguration(0, 0.5)
    % Window resolution (pixel): [small window]
    scr.windX_px = scr.windX_px/2;
    scr.windY_px = scr.windY_px/2;
    scr_dim = [0, 0, scr.windX_px, scr.windY_px];
    const.vpixx = 0;
    [~] = Screen('Preference', 'SkipSyncTests', 1); % change to 0 for real exp
else
    % Window resolution is Screen resolution: [fullscreen]
    scr.windX_px = scr.windX_px;
    scr.windY_px = scr.windY_px;
    scr_dim = []; % PTB says better precision when empty
end

% vpixx
if const.vpixx == 1
    Datapixx('Open')
    Datapixx('StopAllSchedules');
    Datapixx('RegWrRd');    % Synchronize DATAPixx registers to local register cache
end

PsychDefaultSetup(2); % assert OpenGL, setup unifiedkeys and unit color range
PsychImaging('PrepareConfiguration'); % First step in starting pipeline
PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
ListenChar(1);                        % Listen for keyboard input

% open a grey window (size defined above)
[const.window, const.windowRect] = PsychImaging('OpenWindow', scr.scr_num, ...
    [.5 .5 .5], scr_dim, [], [], [], [], []);
% [window, windowRect] = PsychImaging('OpenWindow', scr.scr_num, ...
%     [.5 .5 .5], [0, 0, scr.windX_px, scr.windY_px], 32, 2, [], [], kPsychNeed32BPCFloat);


% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(const.windowRect);

if strcmp(scr.experimenter, 'NYUADScanner') ||  strcmp(scr.experimenter, 'Rania') % account for projector offset at NYUAD
    disp('ACCOUNTING FOR VERTICAL OFFSET BETWEEN PROJECTED IMAGE AND DISPLAY')
    yCenter = yCenter*yOffset_percent;
end

scr.windCenter_px = [xCenter, yCenter];

% Flip to clear
scr.vbl = Screen('Flip', const.window);

% Query the frame duration
scr.ifi = Screen('GetFlipInterval', const.window);

% Enable alpha-blending, set it to a blend equation useable for linear
% additive superposition. This allows to linearly
% superimpose gabor patches in the mathematically correct manner, should
% they overlap. Alpha-weighted source means: The 'globalAlpha' parameter in
% the 'DrawTextures' can be used to modulate the intensity of each pixel of
% the drawn patch before it is superimposed to the framebuffer image, ie.,
% it allows to specify a global per-patch contrast value:
%Screen('BlendFunction', const.window, GL_ONE, GL_ONE);
Screen('BlendFunction', const.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% Set drawing to maximum priority level
topPriorityLevel = MaxPriority(const.window);
Priority(topPriorityLevel);

% .mat file
save(const.scr_fileMat,'scr');

end
