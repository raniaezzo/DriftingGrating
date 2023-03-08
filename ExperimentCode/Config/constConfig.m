function [const]=constConfig(scr,const,expDes)
% ----------------------------------------------------------------------
% [const]=constConfig(const)
% ----------------------------------------------------------------------
% Goal of the function :
% Compute all constant data of this experiment.
% ----------------------------------------------------------------------
% Input(s) :
% scr : window pointer
% const : struct containg previous constant configurations.
% ----------------------------------------------------------------------
% Output(s):
% const : struct containing all constant data.
% ----------------------------------------------------------------------

%% Stimulus Properties             

% stimulus spatial frequency
const.stimSF_cpd = 1;                                   % cycle per degree
const.stimSF_cpp = const.stimSF_cpd/vaDeg2pix(1, scr);  % cycle per pixel
const.stimSF_radians = const.stimSF_cpp*(2*pi);         % in radians
const.stimSF_ppc = ceil(1/const.stimSF_cpp);            % pixel per cycle

% stimulus speed
const.stimSpeed_cpd = 8;                                    % cycles per degree
const.stimSpeed_cps = const.stimSpeed_cpd*const.stimSF_cpd; % cycles per sec
const.stimSpeed_ppc = 1/const.stimSF_cpp;                   % pixel per cycle (without ceil, for precise speed)

% stimulus size: full field
const.stimRadius_xpix = round(scr.windY_px/2 * scr.maxDiam_percent); % constrain X by Y
const.stimRadius_ypix = round(scr.windY_px/2 * scr.maxDiam_percent);
const.stimRadius_deg =  pix2vaDeg(const.stimRadius_xpix,scr);

% stimulus location
const.stimCenterEcc_deg = 0; % 0 degrees eccentricity (center)

% stimulus contrast
const.contrast = .5;

% Initialize displaying of grating (to save time for initial build):
const.grating_halfw= scr.windY_px / 2;
const.visiblesize=2*const.grating_halfw+1;

%gratingtex=ones(2*(const.grating_halfw) +const.stimSF_ppc+1, ...
%    2*(const.grating_halfw) +const.stimSF_ppc+1, 2) * 0.5;
x = meshgrid(-const.grating_halfw:const.grating_halfw + const.stimSF_ppc, 1);
signal=0.5 + 0.5.*cos(const.stimSF_radians*x);
signal = signal.*0.5+0.25; % new
gratingtex = repmat(signal, [length(signal),1]);
%gratingtex(:,:,1) = grating;
%gratingtex(:,:,2) = grating.*const.contrast;
const.gratingtex=Screen('MakeTexture', const.window, gratingtex);

% const.gratingtex = CreateProceduralGabor(const.window, scr.windX_px, scr.windY_px, [],...
%     [0.5 0.5 0.5 0.0], 1, 0.5);

filterparam = const.stimRadius_ypix; % for circular aperature diam
%filterparam = grating_halfw/2; % figure out sigma thresh? (not using anyway)

% Compute outer boundary:
% maskOuter=ones(2*(const.grating_halfw) +const.stimSF_ppc+1, ...
%     2*(const.grating_halfw) +const.stimSF_ppc+1, 2) * 0.5;
% [x,y] = meshgrid(-const.stimSF_ppc/2-const.grating_halfw:const.grating_halfw +const.stimSF_ppc/2, ...
%     -const.stimSF_ppc/2-const.grating_halfw:const.grating_halfw + const.stimSF_ppc/2);
Y_halfdiam = round((scr.windY_px)/2);
maskOuter=ones(Y_halfdiam*2+1, Y_halfdiam*2+1, 2) * 0.5;
[x,y] = meshgrid(-Y_halfdiam:Y_halfdiam, -Y_halfdiam:Y_halfdiam);
%maskOuter=ones(scr.windX_px+1, scr.windY_px+1, 2) * 0.5;
%[x,y] = meshgrid(-scr.windY_px/2:scr.windY_px/2, -scr.windX_px/2:scr.windX_px/2);
maskOuter(:, :, 2)= round(1 * (1 - exp(-((x/filterparam).^2)-((y/filterparam).^2)))); % round for circular aperature
% add cols to fill horizontal space
addcol1 = repmat([0.5], [Y_halfdiam*2+1, round((scr.windX_px-(Y_halfdiam*2+1))/2),1]);
addcol2 = repmat([1], [Y_halfdiam*2+1, round((scr.windX_px-(Y_halfdiam*2+1))/2),1]);
addcols = cat(3, addcol1, addcol2);
maskOuter = [addcols, maskOuter, addcols];
const.maskOutertex=Screen('MakeTexture', const.window, maskOuter);

filterparam = 0.08*const.stimRadius_ypix;

% Compute inner boundary:
maskInner=ones(2*(const.grating_halfw) +const.stimSF_ppc+1, ...
    2*(const.grating_halfw) +const.stimSF_ppc+1, 2) * 0.5;
[x,y] = meshgrid(-const.stimSF_ppc/2-const.grating_halfw:const.grating_halfw +const.stimSF_ppc/2, ...
    -const.stimSF_ppc/2-const.grating_halfw:const.grating_halfw + const.stimSF_ppc/2);
maskInner(:, :, 2)= round((exp(-((x/filterparam).^2)-((y/filterparam).^2))));
const.maskInnertex=Screen('MakeTexture', const.window, maskInner);

% prepare input for stimulus
const.phaseLine = rand(1, expDes.nb_trials) .* 360;
% const.aspectRatio = 1;
% const.gaussianSigma = 0;

%% Fixation Properties

%const.fixationRadius_deg = 0.1;
%const.fixationRadius_px = vaDeg2pix(const.fixationRadius_deg,scr);
const.fixationRadius_px = 0.03*const.stimRadius_ypix;
const.fixationRadius_deg = pix2vaDeg(const.fixationRadius_px, scr);

const.expStart = 0;

%% PTB orientation/direction conversion

orientationids = 0:45:315; ptborientation = {90, 45, 0, 135, 90, 45, 0, 135};
const.maporientation = containers.Map(orientationids,ptborientation);

directionids = 0:45:315; ptbdirection = {180, 135, 90, 45, 0, 315, 270, 225};
const.mapdirection = containers.Map(directionids,ptbdirection);

%% Saving procedure :

% .mat file
save(const.const_fileMat,'const');

end