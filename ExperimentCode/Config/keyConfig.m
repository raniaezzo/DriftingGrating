function [my_key]=keyConfig(const)
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

% keyboard
[keyboardIndices, productNames, ~] = GetKeyboardIndices;
for i=1:length(productNames)                                               % for each possible devicesca
    if strcmp(productNames{i},const.keyboard)                                % compare the name to the name you want
        my_key.keyboardID=keyboardIndices(i);                                   % grab the correct id, and exit loop
        break;
    end
end

my_key.escape       = KbName('ESCAPE');
my_key.space        = KbName('Space');
my_key.rightArrow   = KbName('RightArrow');
my_key.leftArrow    = KbName('LeftArrow');
my_key.Trigger      = KbName('T');

end