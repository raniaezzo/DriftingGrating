%% General experimenter launcher %%
%  =============================  %

% Clean up:

sca; Screen('CloseAll'); 
clear functions; clear mex;
close all; clear all; clc;
rng('default'); % run-specific seed set later (this just re-initializes)
                % differs between runs: timing for fixation task & grating phase order
                % constant across all runs: stimulus condition order

[~] = Screen('Preference', 'SkipSyncTests', 1); % change to 0 for real exp

% Initialization
warning('off');        % do not print warnings
const.DEBUG = 0;       % if debug==1, vPixx and eyetracking off
const.EL_mode = 0;     % will later be forced to 0 if const.DEBUG = 1

% Try to find correct path:
dir = (which('expLauncher'));  % find main script
cd(fileparts(fileparts(dir))); % go to general experiment directory

% Verify that path is correct
[MainDirectory, ~] = fileparts(pwd);
[~, MainFolder] = fileparts(MainDirectory);
if ~strcmp(MainFolder, 'DriftingGrating')
    disp('Not in correct directory. Please run code from ExperimentCode directory.')
else
    % Ensure all folders are on path
    addpath(genpath(pwd)); % add folders with experimental code
    
    % Main experimental code
    main(const);
    clear expDes
end
