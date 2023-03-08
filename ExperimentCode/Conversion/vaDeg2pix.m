function [pix]= vaDeg2pix(vaDeg,scr)
% ----------------------------------------------------------------------
% [pixX, pixY] = vaDeg2pix(vaDeg,scr)
% ----------------------------------------------------------------------
% Goal of the function :
% Convert visual angle (degree) to pixel ( x and y )
% ----------------------------------------------------------------------
% Input(s) :
% vaDeg = size in visual angle (deg)                     ex : = 1
% scr   = screen configurations : scr.scrX_px (pix)      ex : = 2880
%                                 scr.scrY_px (pix)      ex : = 1800
%                                 scr.scrX_cm (cm)       ex : = 28.5
%                                 scr.scrY_cm (cm)       ex : = 17.9
%                                 scr.scrViewingDist_cm  ex : = 60
% ----------------------------------------------------------------------
% Output(s):
% pix  = size in pixel                                   ex : = 101
% ----------------------------------------------------------------------

deg_per_px = rad2deg(atan2(.5 * scr.scrY_cm, scr.scrViewingDist_cm)) / ...
    (.5 * scr.scrY_px);

pix = round(vaDeg / deg_per_px);

end
