function [cm] = vaDeg2cm (vaDeg,scr)
% ----------------------------------------------------------------------
% [cm] = vaDeg2cm (vaDeg,dist)
% ----------------------------------------------------------------------
% Goal of the function :
% Convert visual angle (degree) in cm
% ----------------------------------------------------------------------
% Input(s) :
% vaDeg = size in visual angle (degree)             ex : deg = 2.0
% scr  = scr.scrViewingDist_cm                      ex : = 60
% ----------------------------------------------------------------------
% Output(s):
% cm    = size in cm                                ex : cm = 12
% ----------------------------------------------------------------------

pix = vaDeg2pix(vaDeg,scr);
cm = pix2cm(pix, scr);

end