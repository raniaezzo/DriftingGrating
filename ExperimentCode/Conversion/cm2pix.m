function [pixX, pixY] = cm2pix(cm,scr)
% ----------------------------------------------------------------------
% [pixX, pixY] = cm2pix(cm,scr)
% ----------------------------------------------------------------------
% Goal of the function :
% Convert cm to pixel ( in x and y screen dimensions )
% ----------------------------------------------------------------------
% Input(s) :
% cm    = size in cm                                 ex : cm = 3
% scr   = screen configurations : scr.scrX_px (pix)  ex : = 2880
%                                 scr.scrY_px (pix)  ex : = 1800
%                                 scr.scrX_cm (cm)   ex : = 28.5
%                                 scr.scrY_cm (cm)   ex : = 17.9
% ----------------------------------------------------------------------
% Output(s):
% pixX  = size in pixel(X)                             ex : = 303
% pixY  = size in pixel(Y)                             ex : = 303  
% ----------------------------------------------------------------------

% round to whole pixels
pix_p_cmX = round(scr.scrX_px/scr.scrX_cm);
pix_p_cmY = round(scr.scrY_px/scr.scrY_cm);

pixX = cm*pix_p_cmX;
pixY = cm*pix_p_cmY;

end