function [cmX,cmY] = pix2cm(pix,scr)
% ----------------------------------------------------------------------
% [cmX,cmY] = pix2cm(pix,scr)
% ----------------------------------------------------------------------
% Goal of the function :
% Convert pix to cm ( x and y )
% ----------------------------------------------------------------------
% Input(s) :
% pix   = size in pixel                                ex : pix = 100
% scr   = screen configurations : scr.scrX_px (pix)    ex : = 2880
%                                 scr.scrY_px (pix)    ex : = 1800
%                                 scr.scrX_cm (cm)     ex : = 28.5
%                                 scr.scrY_cm (cm)     ex : = 17.9
% ----------------------------------------------------------------------
% Output(s):
% cmX  = size in cm(X)                             ex : = 0.99
% cmY  = size in cm(Y)                             ex : = 0.99  
% ----------------------------------------------------------------------

cmX = (pix*(scr.scrX_cm))/scr.scrX_px;
cmY = (pix*(scr.scrY_cm))/scr.scrY_px;


cmX = (pix*(28.48))/scr.scrX_px;
cmY = (pix*(17.80))/scr.scrY_px;


cmX = round(cmX, 2); % round to 2 decimals for minute differences
cmY = round(cmY, 2);

end