function projectDir = setup_user(username, projectDir)

% restoredefaultpath
username

% user specific locations
switch(username)
    case 'm1'
        projectDir=  '';
    case 'm2'
        projectDir=  '';
    case 'Server'
        % projectDir=  '/Volumes/Vision/MRI/Decoding';
    case {'Bas', 'bas'}
        freesurferDir = '/Applications/freesurfer/7.2.0';
        fslDir = '/usr/local/fsl';
        githubDir = '/Users/rokers/Documents/GitHub';
    case 'Puti'
        % projectDir=  '/Users/pw1246/Desktop/MRI/Retinotopy';
        freesurferDir = '/Applications/freesurfer/7.2.0';
        githubDir = '/Users/pw1246/Documents/GitHub';
    case 'omnia'
        % projectDir='/Volumes/Vision/MRI/Retinotopic_Mapping';
        fslDir = '/usr/local/fsl';
        % projectDir=  '/Volumes/Vision/MRI/Retinotopy_Class';
        githubDir = '/Users/omniahassanin/Documents/GitHub';
        freesurferDir = '/Applications/freesurfer/7.2.0';
    case {'class'}
        % projectDir=  '/Users/nyuad/Dropbox/fMRI/Retinotopy';
        freesurferDir = '/Applications/freesurfer/7.2.0';
        githubDir = '/Users/nyuad/Dropbox/GitHub';
    case {'caterina'}
        % projectDir=  '/Volumes/Vision/MRI/Retinotopy_Class';
        %projectDir = '/Users/cp3488/Data/Projects/Retinotopy_NYUAD';
        freesurferDir = '/Applications/freesurfer';
        githubDir = '/Users/cp3488/Data/GitHub';
    case {'classbas'}
        % projectDir=  '/Users/br87/Dropbox/Courses_Dropbox/Lab_Visual_Neuroscience/fMRI/Retinotopy';
        % freesurferDir = '/Applications/freesurfer/7.2.0';
        freesurferDir = '/Applications/freesurfer';
        githubDir = '/Users/br87/Dropbox/Courses_Dropbox/Lab_Visual_Neuroscience/GitHub';
    case {'rania'}
        % projectDir=  '/Users/br87/Dropbox/Courses_Dropbox/Lab_Visual_Neuroscience/fMRI/Retinotopy';
        % freesurferDir = '/Applications/freesurfer/7.2.0';
        fslDir = '/usr/local/fsl/';
        freesurferDir = '/Applications/freesurfer/7.2.0';
        githubDir = '/Users/rje257/Documents/GitHub/';
    case{'sawtooth'}
        fslDir = '/usr/local/fsl/';
        freesurferDir = '/Applications/freesurfer/';
        githubDir = '/Users/rje257/Documents/GitHub/';
end

% setup toolboxes
% TODO: Consider using toolboxtoolbox instead - https://github.com/ToolboxHub/ToolboxToolbox
%addpath(genpath(fullfile(githubDir, 'jsonlab'))); % https://github.com/fangq/jsonlab % put first, to resolve mergestruct.m collision
addpath(genpath(fullfile(githubDir, 'MRI_tools'))); % https://github.com/WinawerLab/MRI_tools
addpath(genpath(fullfile(githubDir, 'GLMdenoise'))); % https://github.com/cvnlab/GLMdenoise % required instead of analyzePRF and cvncode when using GLMdenoise
addpath(genpath(fullfile(githubDir, 'GLMsingle')));
addpath(genpath(fullfile(githubDir, 'analyzePRF'))); % https://github.com/cvnlab/analyzePRF
addpath(genpath(fullfile(githubDir, 'cvncode'))); % https://github.com/cvnlab/cvncode
% addpath(genpath(fullfile(githubDir, 'knkutils'))); % https://github.com/cvnlab/knkutils
addpath(genpath(fullfile(githubDir, 'gifti'))); % https://github.com/gllmflndn/gifti
addpath(genpath(fullfile(githubDir, 'reToolbox'))); % https://github.com/gllmflndn/gifti

% add /helper functions - eventually need to be removed, as those functions
% are modified from the toolboxes above
%addpath(genpath(fullfile(githubDir, 'retinotopy-nyuad')));

addpath(genpath(fullfile(githubDir, 'vistasoft'))); % for model
addpath(genpath(fullfile(githubDir, 'prfVista'))); % for model
addpath(genpath(fullfile(githubDir, 'bidsVistapRF'))); % just to save from .mat to surfs etc.

% FSL settings
PATH = getenv('PATH'); setenv('PATH', [PATH ':' fslDir '/bin']); % add freesurfer/bin to path
setenv('FSLDIR', fslDir);
% setenv('FSLOUTPUTTYPE','NIFTI_GZ'); %added to tell where to save the fsl outputs

% freesurfer settings
PATH = getenv('PATH'); setenv('PATH', [PATH ':' freesurferDir '/bin']); % add freesurfer/bin to path
setenv('FREESURFER_HOME', freesurferDir);
addpath(genpath(fullfile(freesurferDir, 'matlab')));
addpath(genpath(fullfile(freesurferDir, 'fsfast')));
setenv('SUBJECTS_DIR', [projectDir '/derivatives/freesurfer']);
