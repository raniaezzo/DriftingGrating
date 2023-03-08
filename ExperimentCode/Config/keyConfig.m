function [my_key]=keyConfig
% ----------------------------------------------------------------------
% [my_key]=keyConfig
% ----------------------------------------------------------------------
% Goal of the function :
% Unify key names and return a structure containing each key names.
% ----------------------------------------------------------------------
% Input(s) :
% none
% ----------------------------------------------------------------------
% Output(s):
% my_key : structure containing all keyboard names.
% ----------------------------------------------------------------------

my_key.escape       = KbName('ESCAPE');
my_key.space        = KbName('Space');
my_key.rightArrow   = KbName('RightArrow');
my_key.leftArrow    = KbName('LeftArrow');
my_key.Trigger      = KbName('T');

end