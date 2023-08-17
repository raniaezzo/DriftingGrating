function [expDes]=designConfig(scr,const)
% ----------------------------------------------------------------------
% [expDes]=designConfig(const)
% ----------------------------------------------------------------------
% Goal of the function :
% Load trial sequence matrix containing condition labels
% used in the experiment.
% ----------------------------------------------------------------------
% Input(s) :
% const : struct containing all constant configurations.
% ----------------------------------------------------------------------
% Output(s):
% expDes : struct containg all trial sequence data.
% ----------------------------------------------------------------------

% save random number generator / seed for each run
expDes.rng = rng(const.run);

%% Experimental sequence

expDes.nb_repeat = 4; % number of repeats

% load in sequence
trialtypes = readtable('Config/trialtypes.csv');
designColNames = trialtypes.Properties.VariableNames;
[expDes.uniqueTrials, ~] = size(trialtypes);

% randomize the unique trials expDes.nb_repeat Xs, then concatenate for
% full run
trialtypesMAT = table2array(trialtypes); trialsequenceMAT = [];
for tt=1:expDes.nb_repeat
    trialsequenceMAT = [trialsequenceMAT; trialtypesMAT(randperm(size(trialtypesMAT,1)), :)];
end

trialsequence = array2table(trialsequenceMAT,'VariableNames',designColNames);

% save attributes
expDes.mainStimTypes = unique(trialsequence.mainCondition); % blank, static, motion
[expDes.nb_trials, ~] = size(trialsequence);
%expDes.nb_repeat = sum(trialsequence.mainCondition == 0); % repeats

% Experimental matrix
trialIDs = 1:expDes.nb_trials;
expDes.trialMat = [trialIDs', table2array(trialsequence)];

%% Experiental timing settings

expDes.runPadding_s = 10;    % 10 sec padding at beginning & end
expDes.stimDur_s  = 3;       % 3 sec stimulus duration
expDes.itiDur_s  = 2;      % 1.5 inter-trial interval
expDes.total_s = (2*expDes.runPadding_s) + (expDes.nb_trials*(expDes.stimDur_s+expDes.itiDur_s));

expDes.runPadding_nFrames  =  round(expDes.runPadding_s/scr.ifi); % # frames
expDes.stimDur_nFrames  =     round(expDes.stimDur_s/scr.ifi); % # frames
expDes.itiDur_nFrames  =      round(expDes.itiDur_s/scr.ifi); % # frames

expDes.totalframes = expDes.runPadding_nFrames+ ...
    expDes.stimDur_nFrames*expDes.nb_trials+ ...
    expDes.itiDur_nFrames*expDes.nb_trials+ ...
    expDes.runPadding_nFrames;
 
expDes.runPadding1Start_frame = 1;
expDes.runPadding1End_frame = expDes.runPadding_nFrames;

expDes.TrialsStart_frame = expDes.runPadding1End_frame+1;
expDes.TrialsEnd_frame = -1 + expDes.TrialsStart_frame+ ...
    expDes.stimDur_nFrames*expDes.nb_trials+ ...
    expDes.itiDur_nFrames*expDes.nb_trials;

expDes.runPadding2Start_frame = expDes.TrialsEnd_frame+1;
expDes.runPadding2End_frame = -1 + expDes.runPadding2Start_frame+expDes.runPadding_nFrames;

%% Saving procedure

save(const.design_fileMat,'expDes');

end
