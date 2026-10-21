function viewGLMsingleR2(subject, project, hemi, filtered, varargin)
% viewGLMsingleR2 - display a GLMsingle full-model R^2 .mgz on the inflated
% cortical surface in FreeSurfer's freeview, with the V1 label boundary
% overlaid.
%
%   viewGLMsingleR2('sub-0037', 'dg', 'lh', false)
%   viewGLMsingleR2('sub-0037', 'dg', 'rh', true)
%
% <subject>  e.g. 'sub-0037'
% <project>  'dg' or 'da'
% <hemi>     'lh' or 'rh' -- which hemisphere to display
% <filtered> true  -> fullmodel_rsquared_filteredvertices.mgz (only the
%                      V1+eccen+vexpl vertices used in the summary table
%                      have nonzero values; see computeGLMsingleR2.m)
%            false -> fullmodel_rsquared.mgz (unfiltered, every vertex)
%
% Optional name/value pairs:
%   'bidsDir'   (default: the standard Project_dg data_bids path)
%   'overlayMax' (default: 30, R2 is a percentage; the low threshold is
%                fixed just above 0 -- see note below -- so this is the
%                only bound worth adjusting)
%
% NOTE on the overlay threshold: the low threshold is fixed at 0.01, not 0.
% Real R2 values in this dataset are never exactly 0 (observed minimum
% ~0.8), so this has no effect on the unfiltered map, but it correctly
% hides the filtered map's explicit "0 = outside the mask" vertices
% (freeview doesn't render values below the overlay's low threshold),
% rather than showing them as if they were real near-zero data.
%
% Requires FreeSurfer's freeview to be on the system PATH (added by
% setup_user.m) and $FREESURFER_HOME to be set.
%
% See also COMPUTEGLMSINGLER2

p = inputParser;
p.addParameter('bidsDir', '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids', @ischar);
p.addParameter('overlayMax', 30, @(x) isnumeric(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;

if ~ismember(hemi, {'lh', 'rh'})
    error('hemi must be ''lh'' or ''rh''');
end

if isempty(getenv('FREESURFER_HOME'))
    error('FREESURFER_HOME is not set. Run setup_user first.');
end

fsRoot = fullfile(opt.bidsDir, 'derivatives', 'freesurfer');
surfFile = fullfile(fsRoot, subject, 'surf', sprintf('%s.inflated', hemi));
labelFile = fullfile(fsRoot, subject, 'label', 'retinotopy_RE', sprintf('%s.V1_REmanual.label', hemi));

derivativesFolder = fullfile(opt.bidsDir, 'derivatives', sprintf('%sGLM', project), 'hRF_glmsingle');
sesDirs = dir(fullfile(derivativesFolder, subject, 'ses-*'));
sesDirs = sesDirs([sesDirs.isdir]);
if numel(sesDirs) ~= 1
    error('expected exactly 1 session folder under %s, found %d', ...
        fullfile(derivativesFolder, subject), numel(sesDirs));
end
sesDir = fullfile(sesDirs(1).folder, sesDirs(1).name);

if filtered
    fileTag = 'fullmodel_rsquared_filteredvertices';
else
    fileTag = 'fullmodel_rsquared';
end
overlayFile = fullfile(sesDir, sprintf('%s.%s.mgz', hemi, fileTag));

requiredFiles = {surfFile, labelFile, overlayFile};
for k = 1:numel(requiredFiles)
    if ~exist(requiredFiles{k}, 'file')
        error('file not found: %s', requiredFiles{k});
    end
end

cmd = sprintf(['freeview -f "%s":overlay="%s":overlay_threshold=0.01,%g:' ...
    'label="%s":label_outline=1:label_color=red &'], ...
    surfFile, overlayFile, opt.overlayMax, labelFile);

fprintf('%s\n', cmd);
[status, out] = system(cmd);
if status ~= 0
    error('freeview failed to launch (status %d): %s', status, out);
end

return
