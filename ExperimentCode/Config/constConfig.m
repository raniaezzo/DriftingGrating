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

% 
const.stimCosEdge_deg = 1; %1.5;
const.stimCosEdge_pix = vaDeg2pix(const.stimCosEdge_deg, scr);

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

% Compute pink noise background with grey circle in the center
Y_halfdiam = round((scr.windY_px)/2); X_halfdiam = round((scr.windX_px)/2);
const.pinknoise_ampSpec = 1.1;
pinknoise = oneoverf(const.pinknoise_ampSpec, const.windowRect(4)+1, const.windowRect(3)+1); %Y_halfdiam*2+1, X_halfdiam*2+1);
% set gray circle with width (filterparam)
masknan = createmask(pinknoise, const.stimRadius_ypix);
background = pinknoise; %.*masknan; % BACK
%background(isnan(background)) = 0.5; % BACK FOR GRAY CIRCLE

const.pinknoise = pinknoise;
const.pinknoiseTex = Screen('MakeTexture', const.window, pinknoise);

maskOuter = cat(3, background, background);
[finalY, finalX, ~] = size(maskOuter);

% make subwindow for cosine ramp for stimulus
imsize = filterparam*2+1;

[x, y] = meshgrid(-imsize/2+0.5:imsize/2-0.5, -imsize/2+0.5:imsize/2-0.5);
[~, r] = cart2pol(x,y);

alpha = zeros(imsize,imsize);

inner_radius = filterparam - const.stimCosEdge_pix; % 384 - 177 (or 44)
outer_radius = filterparam;

% ust added
cosX = linspace(-pi, 0, 1001);
cosY = (cos(cosX)+1)/2;
aa = r - inner_radius;
test = aa(aa < outer_radius-inner_radius);
test = test(test > inner_radius-inner_radius);
maxR = max(test); minR = min(test);

for ii = 1:imsize
    for jj = 1:imsize
        if r(ii,jj) < inner_radius
            alpha(ii,jj) = 0;
        elseif r(ii,jj) < outer_radius
            dist = r(ii,jj)-inner_radius;
            pick = (dist - minR) / (maxR - minR) *1000;
            choose = round(pick);
            alpha(ii,jj) = cosY(choose+1);
            %alpha(ii,jj) = (1-cosd(r(ii,jj)-inner_radius))/2;
        else
            alpha(ii,jj) = 1;
        end
    end
end

[aRows, aCols] = size(alpha);
addcols = repmat([1], [aRows, round((finalX-(aCols))/2),1]);
alpha = [addcols, alpha, addcols];
[aRows, aCols] = size(alpha);
addrows = repmat([1], [round((finalY-(aRows))/2), aCols,1]);
alpha = [addrows; alpha; addrows];
maskOuter(:,:,2) = alpha;

const.maskOutertex=Screen('MakeTexture', const.window, maskOuter);
const.maskOuter = maskOuter; % save in case needed later

filterparam = 0.08*const.stimRadius_ypix;

inner_radius = filterparam + const.stimCosEdge_pix;
outer_radius = filterparam;

% again (for inner)
imsize = round(filterparam*2)+1;

[x, y] = meshgrid(-imsize/2+0.5:imsize/2-0.5, -imsize/2+0.5:imsize/2-0.5);
[~, r] = cart2pol(x,y);

alpha = zeros(imsize,imsize);

cosX = linspace(0, pi, 1001);
cosY = (cos(cosX)+1)/2;
aa = r - inner_radius;
test = aa(aa < outer_radius-inner_radius);
test = test(test > inner_radius-inner_radius);
maxR = max(test); minR = min(test);

for ii = 1:imsize
    for jj = 1:imsize
        if r(ii,jj) < inner_radius
            alpha(ii,jj) = 0;
        elseif r(ii,jj) < outer_radius
            dist = r(ii,jj)-inner_radius;
            pick = (dist - minR) / (maxR - minR) *1000;
            choose = round(pick);
            alpha(ii,jj) = cosY(choose+1);
            %(1-cosd(r(ii,jj)-inner_radius))/2
            %cosY(round(r(ii,jj)-inner_radius));
            %(1-cosd(r(ii,jj)-inner_radius))/2;
        else
            alpha(ii,jj) = 1;
        end
    end
end
%

% Compute inner boundary:
tempSize = (const.grating_halfw) +const.stimSF_ppc; %2*(const.grating_halfw) +const.stimSF_ppc;
tempSize = round(tempSize);
if rem(tempSize,2)==0
    tempSize = tempSize + 1; % make odd to perfectly center
end
maskInner=ones(tempSize, tempSize, 2) * 1; %0.5;
[pn_size_Y, pn_size_X] = size(pinknoise);
pinknoiseCrop = pinknoise((pn_size_Y/2)-(floor(tempSize/2)):(pn_size_Y/2)+(floor(tempSize/2)), ...
    (pn_size_X/2)-(floor(tempSize/2)):(pn_size_X/2)+(floor(tempSize/2)));
maskInner(:,:,1) = pinknoiseCrop;

[x1,y1] = meshgrid(ceil(-const.stimSF_ppc/2-const.grating_halfw):floor(const.grating_halfw +const.stimSF_ppc/2), ...
    ceil(-const.stimSF_ppc/2-const.grating_halfw):floor(const.grating_halfw + const.stimSF_ppc/2));

% [x, y] = meshgrid((pn_size_Y/2)-(floor(tempSize/2)):(pn_size_Y/2)+(floor(tempSize/2)), ...
%     (pn_size_X/2)-(floor(tempSize/2)):(pn_size_X/2)+(floor(tempSize/2)));

x = imresize(x1,[tempSize tempSize]); y = imresize(y1,[tempSize tempSize]); 

maskInner(:, :, 2)= (exp(-((x/filterparam).^2)-((y/filterparam).^2))); % removed round to make gaussian
const.maskInnertex=Screen('MakeTexture', const.window, maskInner);
const.maskInner = maskInner; % save in case needed later

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