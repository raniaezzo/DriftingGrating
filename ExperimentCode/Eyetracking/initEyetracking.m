function EL = initEyetracking(const, w_ptr)

%% Initialize Eyetracker

[EL, exitFlag] = initEyelinkStates('eyestart', w_ptr, {const.eyeFileName, const});
if exitFlag, return, end

EL.eyeDataDir = const.eyeDataDir;
EL.eyeFile = const.eyeFileName; 