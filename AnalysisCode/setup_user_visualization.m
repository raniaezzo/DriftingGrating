function [projectDir, githubDir, freesurferDir] = setup_user_visualization(username)

% user specific locations
switch(username)
    case {'rania'}
       projectDir=  '/Users/rje257/Desktop/MRI_data/dg/'; %'/Volumes/Vision/UsersShare/Rania/LVN_2023/'; 
       freesurferDir = '/Applications/freesurfer/7.2.0';
       githubDir = '/Users/rje257/Documents/GitHub/'; %'/Volumes/Vision/UsersShare/Rania/LVN_2023/Github/'; 
    case {'classbas'}
        projectDir=  '/Users/br87/Dropbox/Courses_Dropbox/Lab_Visual_Neuroscience/fMRI/Retinotopy';
        freesurferDir = '/Applications/freesurfer';
        githubDir = '/Users/br87/Dropbox/Courses_Dropbox/Lab_Visual_Neuroscience/GitHub';
end

% setup toolboxes
% TODO: Consider using toolboxtoolbox instead - https://github.com/ToolboxHub/ToolboxToolbox
addpath(genpath(fullfile(githubDir, 'cvncode'))); % https://github.com/cvnlab/cvncode
addpath(genpath(fullfile(githubDir, 'knkutils'))); % https://github.com/cvnlab/knkutils
addpath(genpath(fullfile(githubDir, 'gifti')));  % rania added
addpath(genpath(fullfile(githubDir, 'GLMdenoise'))); 

% for localizer only
addpath(genpath(fullfile(githubDir, 'wpToolbox')));
addpath(genpath(fullfile(githubDir, 'atlasmgz')));

%addpath(genpath(fullfile(githubDir, 'GLMsingle')));

% add /helper functions - eventually need to be removed, as those functions
% are modified from the toolboxes above
addpath(genpath(fullfile(githubDir, 'retinotopy-nyuad'))); 

% freesurfer settings
PATH = getenv('PATH'); setenv('PATH', [PATH ':' freesurferDir '/bin']); % add freesurfer/bin to path
setenv('FREESURFER_HOME', freesurferDir);
addpath(genpath(fullfile(freesurferDir, 'matlab')));
setenv('SUBJECTS_DIR', [projectDir '/derivatives/freesurfer']); 


fslDir = '/usr/local/fsl';
PATH = getenv('PATH'); setenv('PATH', [PATH ':' fslDir '/bin']); % add fsl/bin to path
setenv('FSLDIR', fslDir);
setenv('PATH', sprintf('/usr/local/bin:%s', getenv('PATH'))); % add /usr/local/bin to PATH
%setenv('PATH', [fsDir '/bin:' getenv('PATH')]);
setenv('PATH', [getenv('PATH') ':/usr/local/fsl/bin']);
