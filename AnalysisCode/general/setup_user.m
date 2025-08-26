function projectDir = setup_user(username, projectDir)

% check current path
currentDir = pwd; % Get current working directory
[parentDir, currentFolder] = fileparts(currentDir); % Get parent directory

if strcmp(currentFolder, 'AnalysisCode')
    disp('Confirmed in correct working directory');
else
    % Check if 'AnalysisCode' folder is in the path
    if contains(parentDir, 'AnalysisCode')
        disp('Setting current directory to AnalysisCode...');
        cd(parentDir); % Change to AnalysisCode folder
    else
        error('AnalysisCode folder containing setup.json must be the parent folder.');
    end
end

% add all analysis scripts to path
addpath(genpath(pwd));

% add user specific softward paths
username
switch(username)
    case {'rania'}
        fslDir = '/usr/local/fsl/';
        freesurferDir = '/Applications/freesurfer/7.2.0';
        githubDir = '/Users/rje257/Documents/GitHub/';
    case{'sawtooth'}
        fslDir = '/usr/local/fsl/';
        freesurferDir = '/Applications/freesurfer/';
        githubDir = '/Users/rje257/Documents/GitHub/';
    case{'JWLAB08M'}
        fslDir = '/Users/rje257/fsl/';
        freesurferDir = '/Applications/freesurfer/7.4.1/';
        githubDir = '/Users/rje257/Documents/GitHub/';
end

% setup toolboxes
% TODO: Consider using toolboxtoolbox instead - https://github.com/ToolboxHub/ToolboxToolbox
addpath(genpath(fullfile(githubDir, 'MRI_tools'))); % https://github.com/WinawerLab/MRI_tools
addpath(genpath(fullfile(githubDir, 'GLMdenoise'))); % https://github.com/cvnlab/GLMdenoise % required instead of analyzePRF and cvncode when using GLMdenoise
addpath(genpath(fullfile(githubDir, 'GLMsingle'))); % https://github.com/cvnlab/GLMsingle
addpath(genpath(fullfile(githubDir, 'analyzePRF'))); % https://github.com/cvnlab/analyzePRF
addpath(genpath(fullfile(githubDir, 'cvncode'))); % https://github.com/cvnlab/cvncode
addpath(genpath(fullfile(githubDir, 'gifti'))); % https://github.com/gllmflndn/gifti
addpath(genpath(fullfile(githubDir, 'reToolbox'))); 
addpath(genpath(fullfile(githubDir, 'vistasoft'))); % for model
addpath(genpath(fullfile(githubDir, 'prfVista'))); % for model
addpath(genpath(fullfile(githubDir, 'bidsVistapRF'))); % just to save from .mat to surfs etc.
addpath(genpath(fullfile(githubDir, 'atlasmgz')));

% FSL settings
PATH = getenv('PATH'); setenv('PATH', [PATH ':' fslDir '/bin']); % add freesurfer/bin to path
setenv('FSLDIR', fslDir);

% freesurfer settings
PATH = getenv('PATH'); setenv('PATH', [PATH ':' freesurferDir '/bin']); % add freesurfer/bin to path
setenv('FREESURFER_HOME', freesurferDir);
addpath(genpath(fullfile(freesurferDir, 'matlab')));
addpath(genpath(fullfile(freesurferDir, 'fsfast')));
setenv('SUBJECTS_DIR', [projectDir '/derivatives/freesurfer']);
