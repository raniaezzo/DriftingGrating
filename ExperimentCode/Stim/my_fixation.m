function my_fixation(scr,const,color)
% ----------------------------------------------------------------------
% my_fixationCross(scr,const,color)
% ----------------------------------------------------------------------
% Goal of the function :
% Draw a fixation cross in the center of the screen
% ----------------------------------------------------------------------
% Input(s) :
% scr = scr variable with contains Window Pointer   ex : scr
% color = color of the circle in RBG or RGBA        ex : color = [0 0 0]
% const = structure containing constant configurations
% ----------------------------------------------------------------------
% Output(s):
% ----------------------------------------------------------------------

% Draw the fixation point
Screen('DrawDots', const.window, scr.windCenter_px, ...
    const.fixationRadius_px, color, [], 2);

end