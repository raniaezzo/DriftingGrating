function [vaDeg] =cm2vaDeg (cm,scr)
% ----------------------------------------------------------------------
% [vaDeg] = cm2vaDeg(cm,scr)
% ----------------------------------------------------------------------
% Goal of the function :
% Convert cm to visual angle (degree)
% ----------------------------------------------------------------------
% Input(s) :
% cm = size in cm                                   ex : 2 cm target
% scr = scr.scrViewingDist_cm                       ex : scr.scrViewingDist_cm = 60
% ----------------------------------------------------------------------
% Output(s):
% vaDeg = size in visual angle (degree)             ex : vaDeg = 1.7799
% ----------------------------------------------------------------------

if length(cm)>1
    cm=cm(1);   % only use X dimension specs (since isotropic)
end

pix = cm2pix(cm, scr);
vaDeg = pix2vaDeg(pix, scr);

end
