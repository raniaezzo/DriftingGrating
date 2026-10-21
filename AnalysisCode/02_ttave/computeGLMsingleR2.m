function [perHemiTable, combinedTable] = computeGLMsingleR2(varargin)
% computeGLMsingleR2 - full-model GLMsingle R^2 on the cortical surface, dg vs da
%
%   [perHemiTable, combinedTable] = computeGLMsingleR2()
%   [perHemiTable, combinedTable] = computeGLMsingleR2('subjects', {'sub-0037'})
%
% For each subject shared between the dg and da projects (both projects
% live under the same bidsDir/derivatives, in dgGLM and daGLM), this:
%
%   1. Loads modelOutput.mat's modelOut{4} (the full GLMsingle model --
%      TYPED, FITHRF_GLMdenoise_RR) and takes its R2 field (percent
%      variance explained, one value per surface vertex, left hemisphere
%      then right -- same vertex order as get_surfsize/writeMGZfile).
%   2. Saves that R2 vector as lh.fullmodel_rsquared.mgz /
%      rh.fullmodel_rsquared.mgz in the same session folder, via the
%      existing writeMGZfile helper.
%   3. Computes the mean R2 within a V1 ROI, restricted to vertices with
%      eccenMin <= eccentricity <= eccenMax and pRF variance explained >
%      vexplMin -- the SAME mask construction as dg_computeGain.m (V1
%      label + prfvista_mov's vexpl/eccen maps, since that vertex set is
%      what all subsequent analyses were selected on), just with this
%      script's own eccenMin/eccenMax defaults (4-6 deg here, vs. 4-8 deg
%      in dg_computeGain.m). The ROI mask is anatomical, not
%      project-specific, so it's built once per subject and reused for
%      both dg and da.
%
% Returns (and saves to outDir) two tables:
%   perHemiTable  - subject x project x hemi (lh, rh), mean R2 in the ROI
%   combinedTable - subject x project, mean R2 in the ROI across both
%                   hemispheres pooled together
%
% See also DG_COMPUTEGAIN, WRITEMGZFILE, GET_SURFSIZE

p = inputParser;
p.addParameter('bidsDir', '/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids', @ischar);
p.addParameter('projects', {'dg', 'da'}, @iscell);
p.addParameter('movRoot', ['/Volumes/Vision/UsersShare/Rania/Project_dg/' ...
    'data_bids/derivatives/prfvista_mov'], @ischar);
p.addParameter('roiDir',  'retinotopy_RE', @ischar);
p.addParameter('roiName', 'V1_REmanual', @ischar);
p.addParameter('vexplMin', 0.1, @(x) isnumeric(x) && isscalar(x));
p.addParameter('eccenMin', 4,   @(x) isnumeric(x) && isscalar(x));
p.addParameter('eccenMax', 6,   @(x) isnumeric(x) && isscalar(x));
p.addParameter('outFileName', 'fullmodel_rsquared', @ischar);
p.addParameter('outDir', ['/Volumes/Vision/UsersShare/Rania/Project_dg/' ...
    'data_bids/derivatives/summaryTables'], @ischar);
p.addParameter('subjects', {}, @iscell);
p.parse(varargin{:});
opt = p.Results;

fsRoot = fullfile(opt.bidsDir, 'derivatives', 'freesurfer');
setenv('SUBJECTS_DIR', fsRoot);

if ~isfolder(opt.outDir)
    mkdir(opt.outDir);
end

if isempty(opt.subjects)
    opt.subjects = findSharedSubjects(opt);
end
fprintf('%d subject(s) shared across %s\n\n', numel(opt.subjects), strjoin(opt.projects, ', '));

perHemiRows = struct('subject', {}, 'project', {}, 'hemi', {}, ...
    'meanR2', {}, 'nVertices', {});
combinedRows = struct('subject', {}, 'project', {}, 'meanR2', {}, 'nVertices', {});

for ii = 1:numel(opt.subjects)
    sub = opt.subjects{ii};
    fprintf('--- %s\n', sub);

    try
        [maskLH, maskRH] = buildROIMask(sub, opt);
        fprintf('    V1 ROI, %g<=eccen<=%g, vexpl>%.2g (from prfvista_mov): lh %d, rh %d vertices\n', ...
            opt.eccenMin, opt.eccenMax, opt.vexplMin, sum(maskLH), sum(maskRH));
    catch ME
        fprintf('    FAILED to build ROI mask: %s\n', ME.message);
        continue
    end

    for pp = 1:numel(opt.projects)
        proj = opt.projects{pp};
        derivativesFolder = fullfile(opt.bidsDir, 'derivatives', sprintf('%sGLM', proj), 'hRF_glmsingle');
        try
            [R2, ses] = loadFullModelR2(derivativesFolder, sub);
        catch ME
            fprintf('    %s: FAILED: %s\n', proj, ME.message);
            continue
        end

        nL = numel(maskLH);
        if numel(R2) ~= nL + numel(maskRH)
            fprintf('    %s: SKIPPED, R2 has %d vertices but ROI mask has %d (surface mismatch)\n', ...
                proj, numel(R2), nL + numel(maskRH));
            continue
        end

        % save mgz (lh.<outFileName>.mgz / rh.<outFileName>.mgz)
        writeMGZfile(opt.bidsDir, sub, ses, R2, derivativesFolder, opt.outFileName);

        % also save a "filtered" derivative: same values, but zeroed outside
        % the V1+eccen+vexpl mask -- i.e. exactly the vertices that went
        % into the summary table's mean R2, and nothing else.
        R2filtered = zeros(size(R2));
        fullMask = [maskLH; maskRH];
        R2filtered(fullMask) = R2(fullMask);
        writeMGZfile(opt.bidsDir, sub, ses, R2filtered, derivativesFolder, ...
            [opt.outFileName '_filteredvertices']);

        R2lh = R2(1:nL);
        R2rh = R2(nL+1:end);

        rLH = R2lh(maskLH);
        rRH = R2rh(maskRH);

        perHemiRows(end+1) = struct('subject', sub, 'project', proj, 'hemi', 'lh', ...
            'meanR2', mean(rLH), 'nVertices', numel(rLH)); %#ok<AGROW>
        perHemiRows(end+1) = struct('subject', sub, 'project', proj, 'hemi', 'rh', ...
            'meanR2', mean(rRH), 'nVertices', numel(rRH)); %#ok<AGROW>

        rBoth = [rLH(:); rRH(:)];
        combinedRows(end+1) = struct('subject', sub, 'project', proj, ...
            'meanR2', mean(rBoth), 'nVertices', numel(rBoth)); %#ok<AGROW>

        fprintf('    %s (%s): wrote %s.mgz, mean R2 in ROI: lh %.2f, rh %.2f, combined %.2f\n', ...
            proj, ses, opt.outFileName, mean(rLH), mean(rRH), mean(rBoth));
    end
end

perHemiTable = mergeRows(struct2table(perHemiRows), ...
    fullfile(opt.outDir, 'glmsingleR2_V1_perHemi.mat'), 'perHemiTable', {'subject', 'project', 'hemi'});
combinedTable = mergeRows(struct2table(combinedRows), ...
    fullfile(opt.outDir, 'glmsingleR2_V1_combined.mat'), 'combinedTable', {'subject', 'project'});

save(fullfile(opt.outDir, 'glmsingleR2_V1_perHemi.mat'), 'perHemiTable');
writetable(perHemiTable, fullfile(opt.outDir, 'glmsingleR2_V1_perHemi.csv'));
save(fullfile(opt.outDir, 'glmsingleR2_V1_combined.mat'), 'combinedTable');
writetable(combinedTable, fullfile(opt.outDir, 'glmsingleR2_V1_combined.csv'));

fprintf('\nwrote summary tables to:\n  %s\n', opt.outDir);
return


% ------------------------------------------------------------------------
function merged = mergeRows(newRows, matFile, varName, keyCols)
% Resumable: if a table from a previous (possibly partial/interrupted) run
% already exists, merge new rows into it -- new rows replace old ones with
% the same key (subject/project[/hemi]), so re-running is idempotent.
if exist(matFile, 'file')
    S = load(matFile, varName);
    old = S.(varName);
    oldKey = strcat(old.(keyCols{1}));
    newKey = strcat(newRows.(keyCols{1}));
    for kk = 2:numel(keyCols)
        oldKey = strcat(oldKey, '|', old.(keyCols{kk}));
        newKey = strcat(newKey, '|', newRows.(keyCols{kk}));
    end
    keep = ~ismember(oldKey, newKey);
    merged = [old(keep, :); newRows];
else
    merged = newRows;
end
return


% ------------------------------------------------------------------------
function subs = findSharedSubjects(opt)
sets = cell(numel(opt.projects), 1);
for pp = 1:numel(opt.projects)
    d = dir(fullfile(opt.bidsDir, 'derivatives', sprintf('%sGLM', opt.projects{pp}), 'hRF_glmsingle', 'sub-*'));
    sets{pp} = {d([d.isdir]).name};
end
subs = sets{1};
for pp = 2:numel(sets)
    subs = intersect(subs, sets{pp});
end
return


% ------------------------------------------------------------------------
function [maskLH, maskRH] = buildROIMask(sub, opt)
% Same construction as dg_computeGain.m: V1 label (union of ventral/dorsal
% not needed for V1 itself) intersected with prfvista_mov's vexpl/eccen
% criteria. Anatomical/subject-level, so it's shared across projects.
sesDirs = dir(fullfile(opt.movRoot, sub, 'ses-*'));
sesDirs = sesDirs([sesDirs.isdir]);
if isempty(sesDirs)
    error('no prfvista_mov session found for %s', sub);
end
sesDir = fullfile(sesDirs(1).folder, sesDirs(1).name);

vexplLH = readMgzVol(fullfile(sesDir, 'lh.vexpl.mgz'));
vexplRH = readMgzVol(fullfile(sesDir, 'rh.vexpl.mgz'));
eccenLH = readMgzVol(fullfile(sesDir, 'lh.eccen.mgz'));
eccenRH = readMgzVol(fullfile(sesDir, 'rh.eccen.mgz'));

lh_label = read_label(sub, sprintf('%s/lh.%s', opt.roiDir, opt.roiName));
rh_label = read_label(sub, sprintf('%s/rh.%s', opt.roiDir, opt.roiName));

roiLH = false(numel(vexplLH), 1); roiLH(lh_label(:,1) + 1) = true;
roiRH = false(numel(vexplRH), 1); roiRH(rh_label(:,1) + 1) = true;

critLH = vexplLH > opt.vexplMin & eccenLH >= opt.eccenMin & eccenLH <= opt.eccenMax;
critRH = vexplRH > opt.vexplMin & eccenRH >= opt.eccenMin & eccenRH <= opt.eccenMax;

maskLH = roiLH & critLH;
maskRH = roiRH & critRH;
return


% ------------------------------------------------------------------------
function vol = readMgzVol(mgzPath)
if ~exist(mgzPath, 'file')
    error('mgz file not found: %s', mgzPath);
end
mri = MRIread(mgzPath);
vol = double(mri.vol(:));
return


% ------------------------------------------------------------------------
function [R2, ses] = loadFullModelR2(derivativesFolder, sub)
sesDirs = dir(fullfile(derivativesFolder, sub, 'ses-*'));
sesDirs = sesDirs([sesDirs.isdir]);
if isempty(sesDirs)
    error('no session folder found under %s', fullfile(derivativesFolder, sub));
elseif numel(sesDirs) > 1
    error('%d session folders found under %s, expected exactly 1', ...
        numel(sesDirs), fullfile(derivativesFolder, sub));
end
ses = sesDirs(1).name;

modelFile = fullfile(sesDirs(1).folder, ses, 'modelOutput.mat');
if ~exist(modelFile, 'file')
    error('modelOutput.mat not found: %s', modelFile);
end

S = load(modelFile, 'modelOut');
if ~isfield(S, 'modelOut') || ~iscell(S.modelOut) || numel(S.modelOut) ~= 4
    error('modelOutput.mat did not contain a 1x4 cell array "modelOut" as expected');
end

m4 = S.modelOut{4};
if ~isfield(m4, 'R2')
    error('modelOut{4} has no R2 field');
end
R2 = double(m4.R2(:));
return
