function main(const)
% ----------------------------------------------------------------------
% main(const)
% ----------------------------------------------------------------------
% Goal of the function :
% Main code of experiment
% ----------------------------------------------------------------------
% Input(s) :
% const : struct containing subject information and saving files.
% ----------------------------------------------------------------------
% Output(s):
% none
% ----------------------------------------------------------------------

% File directory :
[const] = dirSaveFile(const);
% Screen configuration :
[scr, const] = scrConfig(const);
disp('Finished scrConfig')

% Keyboard configuration :
[my_key] = keyConfig(const);
disp('Finished keyConfig')

% Experimental design configuration :
[expDes] = designConfig(scr,const);
disp('Finished designConfig')

% Experimental constant :
[const] = constConfig(scr,const, expDes);
disp('Finished constConfig')

% Instruction file :
[textExp] = instructionConfig;
disp('Finished instrcutionConfig')

% Initialize eyetracking
if const.EL_mode
    if strcmp(scr.experimenter, 'NYUNYScanner')
        Eyelink('SetAddress','192.168.1.5');
    end
    const.EL = initEyetracking(const, const.window);
else
    const.EL = [];
end

% Main part :
if const.expStart;ListenChar(2);end
[responses, trial_onsets] = runTrials(scr,const,expDes,my_key,textExp);

% End
overDone(const, responses, trial_onsets)

end