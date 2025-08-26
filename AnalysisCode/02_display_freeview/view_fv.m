function freeview_cmd = view_fv(subject, bidsDir, varargin)

% view overlay ontop of inflated surface in freeview

% subject - e.g. 'sub-0248' - subject ID
% bidsDir - e.g. '~/Documents/MRI/bids' - BIDS directory
% hemi - e.g. - 'l' - optional, if not defined then we will plot both hemisphere
% varargin - e.g. 'maps/motion:rainbow' - 'subfolderName/mgzFileName:colorMapName' - look for matching files only in the defined subfolder within derivatives
%          - e.g. 'motion:rainbow' - 'mgzFileName:colorMapName' - look for matching files across all subfolders in the derivatives directory
%          - e.g. 'motion' - no colormap

% example usage:
% subject = 'sub-0248'
% bidsDir = '~/Documents/MRI/bids';
% view_fv(subject, bidsDir, 'maps/motion', 'MyelinMap',);
% % or
% view_fv(subject, bidsDir, 'lh', 'mt+2', 'cw:rainbow','prfvista_mov/eccen');

bidsDir = [bidsDir '/derivatives'];
%% check if BIDS directory exsits
if ~exist(bidsDir,'dir')
    error(['BIDS derivatives directory does not exist <' bidsDir '>'])
end
if contains(bidsDir,'~') 
    bidsDir = strrep(bidsDir, '~', getenv('HOME')); % somehow colormaps cannot be find when using ~
end
%% check if subject directory exsits
if ~contains(subject, 'sub-')
else
    tmp = strsplit(subject, '-');
    subject = tmp{2};
end
subfolder = dir(sprintf('%s/freesurfer/*%s*',bidsDir,subject)); % in freesurfer folder check for any subject folder matches our subject ID
subfolderName = subfolder([subfolder.isdir]).name; % get the folder name 
subjectDir = sprintf('%s/freesurfer/%s',bidsDir,subfolderName); % build the path for subject directory
if ~exist(subjectDir,'dir') % check if this subject directory exists
    error(['subject freesurfer directory for <' subject '> does not exist: <' subjectDir '>'])
end

%% check if this subject has more than one names
folderParts=strsplit(subfolderName, '_'); % break the folder name in parts
tmpID=strsplit(folderParts{1}, '-'); % check first part
subject = []; 
subject{1} = ['sub-' tmpID{2}]; % get AD ID
if contains(subfolderName, 'NY')
    tmpID=strsplit(folderParts{2}, '-'); % check second part
    subject{2} = ['sub-' tmpID{2}]; % get NY ID
else
    subject{2} = 'NA'; % no NY ID
end
if contains(subfolderName, 'XXXX') % no AD ID
    subject{1} = 'NA';
end
%% check if freeview exsits
tmpDir = dir(fullfile('/Applications/freesurfer/*'));
for ii = 1:length(tmpDir)
    if ~ismember(tmpDir(ii).name,{'.','..'})
        fvDir = [tmpDir(ii).folder '/' tmpDir(ii).name];
        setenv('FREESURFER_HOME', fvDir); 
        fvCmd = [fvDir '/bin/freeview'];
    end
end
setenv('SUBJECTS_DIR', [bidsDir '/freesurfer']);
%% check which hemisphere to plot
if ismember(lower(varargin{1}), {'l','lh','left'})
    hemi = {'l'};
    varargin = varargin(2:end);
elseif ismember(lower(varargin{1}), {'r','rh','right'})
    hemi = {'r'};
    varargin = varargin(2:end);
else
    hemi = {'l','r'};
end
%% loop through each hemisphere
cmd = [];
for whichHemi = 1:numel(hemi)

    % find the inflated surface
    inflated = sprintf('%s/surf/%sh.inflated',subjectDir,hemi{whichHemi});
    % check if the inflated surf file exists
    if ~exist(inflated,'file')
        error(['<' hemi{whichHemi} 'h.inflate> not found in <' subjectDir '/surf>']);
    end

    % loop through each overlay name provided in varargin
    overlayCmd = ''; % prepare the command for overlay
    for whichOverlay = 1:length(varargin)
        % parse the input to get the folder, overlay, and colormap name
        [folder, fileParts] = parse_input(varargin{whichOverlay},bidsDir,subject,hemi{whichHemi});
        overlayName = fileParts{2};
        colormapName = fileParts{3}; % rename it so we don't forget what this fileParts variable is
        for whichFolder = 1:numel(folder)
            whereOverlay = sprintf('%s/%sh.%s.mgz', folder{whichFolder}, hemi{whichHemi}, overlayName);
            overlayCmd = [overlayCmd, sprintf(':overlay=%s', whereOverlay)];
            % add the customized color map if we have one
            overlayCmd = custom_colormap(bidsDir,overlayName,colormapName,overlayCmd,hemi{whichHemi});
        end
    end
    if isempty(overlayCmd)      
        warning(['none of the overlay files are found for ' hemi{whichHemi} 'h hemi'])
    end
    whichLabel = sprintf('%s/label/Glasser2016/%sh.MT.label', subjectDir,hemi{whichHemi});
    if isfile(whichLabel)
        cmd = sprintf('%s -f %s%s:label=%s:label_outline=1:label_color=black:label_opacity=0.5',cmd,inflated,overlayCmd,whichLabel);
    else
        cmd = sprintf('%s -f %s%s',cmd,inflated,overlayCmd);
    end
end

freeview_cmd = sprintf('%s%s &',fvCmd, cmd);
system(freeview_cmd);

end

%% sub functions
%% which folders are my overlay files in
function [folder, fileParts] = parse_input(whichFile,bidsDir,subject,iHemi)

    % create temp variable fileParts that is {whichFolder} {whichOverlay} {whichColormap} {whichHemi}
    fileParts = regexp(whichFile, '[:/]', 'split'); % if folder and colormap are specified
    if ~contains(whichFile,'/') % if folder is not specified
        fileParts = ['**' fileParts];
    end
    if ~contains(whichFile,':') % if colormap is not specified
        fileParts = [fileParts 'NA'];
    end

    fileParts = [fileParts iHemi]; % which hemi
    
    folder = {};
    for iSub = 1:numel(subject)
        tmpFolder = find_my_files(subject{iSub},bidsDir,fileParts); % look for the overlay files
        if isfolder(tmpFolder)
            folder = [folder tmpFolder];
        end
    end

    if isempty(folder)
        warning(['<' iHemi 'h.' fileParts{2} '.mgz> not found in any subfolders for this subject']);
        warning(['skipping ' fileParts{2} ' for ' iHemi 'h']);
    end

end

%% return all sub folders that contains wanted overlay mgz file
function folder = find_my_files(whichSub,bidsDir,fileParts)
    %tmpFiles = dir(sprintf('%s/%s/%s/%s%sh.%s.mgz', bidsDir, fileParts{1}, whichSub, '**/', fileParts{4}, fileParts{2}));

    subDirs = dir(fullfile(bidsDir, '*'));
    
    if strcmp(fileParts{1},'**')
        subDirs = subDirs([subDirs.isdir] & ~ismember({subDirs.name}, {'.', '..', 'derivatives','freesurfer'}));
    else
        subDirs = subDirs([subDirs.isdir] & ismember({subDirs.name}, {fileParts{1}}));
    end

    tmpFiles = [];
    for subDir = subDirs'
        tmpFiles = [tmpFiles; dir(fullfile(bidsDir, subDir.name, whichSub, '/**/', [fileParts{4} 'h.' fileParts{2} '.mgz']))];
    end

    folder = cell(numel(tmpFiles),1);
    for iFile = 1:length(tmpFiles)
        if ~tmpFiles(iFile).isdir
            folder{iFile}  = tmpFiles(iFile).folder;
        end
    end
end

%% load customized color maps if we have
function overlayCmd = custom_colormap(bidsDir,overlayName,colormapName, overlayCmd, iHemi)
    
    switch overlayName
        case 'myelin'
            whichMap = 'mye';
        case 'eccen'
            whichMap = 'eccentricity_color_scale';
        case'angle_adj'
           whichMap  = ['angle_corr_' iHemi 'h_color_scale'];
        otherwise
            whichMap = [];
    end
    
    if ismember(overlayName,{'mt+1','mt+2','mt+12','bio1','cw','ccw','cwccw','in','out','inout','upper','lower'})
        whichMap = 'rainbow'; % default rainbow color map for bold beta weights 0:0.1:0.7
    end
    
    if ~strcmp(colormapName,'NA') % mannual input overwrites default
        whichMap = colormapName;
    end
    
    whereMap = sprintf('%s/freesurfer/%s',bidsDir,whichMap); % find the colormap   
    if isfile(whereMap)
        overlayCmd = [overlayCmd,sprintf(':overlay_custom=%s',whereMap)];
    end
    
    if ~isempty(whichMap) && ~exist(whereMap) % if it is specified but doesn't exist
        disp(['customized color map <' whichMap '> not found in <' bidsDir '/freesurfer/>'])
        % check if we have a default map for it
        if ismember(overlayName,{'mt+1','mt+2','mt+12','bio1','cw','ccw','cwccw','in','out','inout','upper','lower'})
            disp(['we will use default color map <rainbow> for ' iHemi 'h.' overlayName '.mgz'])
            whereMap = sprintf('%s/freesurfer/%s',bidsDir,'rainbow'); % find the colormap
            if isfile(whereMap)
                overlayCmd = [overlayCmd,sprintf(':overlay_custom=%s',whereMap)];
            else
                disp(['could not locate default color map <rainbow> in <' bidsDir '/freesurfer/>'])
            end
        end
    end

end
