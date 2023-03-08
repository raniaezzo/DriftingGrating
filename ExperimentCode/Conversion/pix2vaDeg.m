function [vaDeg]= pix2vaDeg (pix,scr)
% ----------------------------------------------------------------------
% [vaDeg]= pix2vaDeg (pix,scr)
% ----------------------------------------------------------------------
% Goal of the function :
% Convert pixel in visual angle
% ----------------------------------------------------------------------
% Input(s) :
% pix   = size in pixel                                     ex : = 101
% scr   = screen configurations : scr.scrX_px (pix)       ex : = 2880
%                                 scr.scrY_px (pix)       ex : = 1800
%                                 scr.scrX_cm (cm)          ex : = 28.5
%                                 scr.scrY_cm (cm)          ex : = 17.9
%                                 scr.scrViewingDist_cm     ex : = 60
% ----------------------------------------------------------------------
% Output(s):
% vaDeg = size in visual angle (degree)                     ex : =  1.7598
% ----------------------------------------------------------------------

deg_per_px = rad2deg(atan2(.5 * scr.scrY_cm, scr.scrViewingDist_cm)) / ...
    (.5 * scr.scrY_px);

vaDeg = pix * deg_per_px;

end
