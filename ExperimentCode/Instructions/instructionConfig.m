function [textExp] = instructionConfig
% ----------------------------------------------------------------------
% [textExp] = instructionConfig
% ----------------------------------------------------------------------
% Goal of the function :
% Write text of calibration and general instruction for the experiment.
% ----------------------------------------------------------------------
% Input(s) :
% (none)
% ----------------------------------------------------------------------
% Output(s):
% textExp : struct containing all text of general instructions.
% ----------------------------------------------------------------------

%% Main instruction :

eyeinstruction = '-------------  Prepare for eyetracking (T)  -------------';

instruction = '-----------------  Waiting for trigger (T)  -----------------';

textExp.eyeinstruction= {eyeinstruction};
textExp.instruction= {instruction};

end