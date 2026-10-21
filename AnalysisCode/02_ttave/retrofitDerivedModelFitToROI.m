function status = retrofitDerivedModelFitToROI(bidsDir, projectName, subj, ses, hRF_setting)
% retrofitDerivedModelFitToROI - trim an EXISTING full-surface
% derivedModelFit.mat down to the ROI-union format (see
% ensureDerivedModelFit.m), without recomputing anything from
% modelOutput.mat.
%
% Saves to local disk first, verifies the local copy loads cleanly, THEN
% copies over the network file (rather than writing to the network mount
% in place) -- direct in-place network writes of large .mat files have
% repeatedly corrupted files earlier in this project's session history.
%
% Returns status: 'trimmed', 'already_trimmed', 'missing', or 'corrupt'
% (corrupt = pre-existing unreadable file, needs regeneration via
% ensureDerivedModelFit.m; not caused by this function, since it bails out
% before writing anything in that case).

derivativesFolder = fullfile(bidsDir, 'derivatives', sprintf('%sGLM',projectName), ...
    sprintf('hRF_%s', hRF_setting), subj, ses);
modelFitFile = fullfile(derivativesFolder, 'derivedModelFit.mat');

if ~exist(modelFitFile, 'file')
    fprintf('  %s/%s: no derivedModelFit.mat found, skipping\n', projectName, subj);
    status = 'missing';
    return
end

info = whos('-file', modelFitFile);
varNames = {info.name};
if ismember('modelfitROI', varNames)
    fprintf('  %s/%s: already in trimmed format, skipping\n', projectName, subj);
    status = 'already_trimmed';
    return
end
if ~ismember('modelfit', varNames)
    fprintf('  %s/%s: CORRUPT (unrecognized format) -- needs regeneration, skipping\n', projectName, subj);
    status = 'corrupt';
    return
end

try
    S = load(modelFitFile, 'modelfit');
catch ME
    fprintf('  %s/%s: CORRUPT (load failed: %s) -- needs regeneration, skipping\n', projectName, subj, ME.message);
    status = 'corrupt';
    return
end

oldBytes = getfield(dir(modelFitFile), 'bytes'); %#ok<GFLD>

nRuns = numel(S.modelfit);
nVerticesTotal = size(S.modelfit{1}, 1);

hSize = get_surfsize(subj);
roisJson = jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')), '..', 'general', 'jsons', 'ROIS_ALL.json')));
roiNames = {roisJson.ROIs.filename};
roiVertexIdx = [];
for i = 1:numel(roiNames)
    roiVertexIdx = union(roiVertexIdx, getROIidxs(subj, roiNames{i}, hSize));
end
roiVertexIdx = sort(roiVertexIdx(:));

modelfitROI = cell(1, nRuns);
for r = 1:nRuns
    modelfitROI{r} = S.modelfit{r}(roiVertexIdx, :);
end
clear S % free the large full-surface array before the local save

% --- safe save: local disk first ---
localTemp = fullfile(tempdir, sprintf('derivedModelFit_%s_%s_retrofit.mat', projectName, subj));
save(localTemp, 'modelfitROI', 'roiVertexIdx', 'nVerticesTotal', '-v7.3');

try
    checkS = load(localTemp);
    if ~isequal(size(checkS.modelfitROI), size(modelfitROI)) || numel(checkS.roiVertexIdx) ~= numel(roiVertexIdx)
        error('local copy failed shape check');
    end
catch ME
    fprintf('  %s/%s: LOCAL SAVE VERIFICATION FAILED (%s) -- network file NOT touched\n', projectName, subj, ME.message);
    delete(localTemp);
    status = 'corrupt';
    return
end

% --- verified locally, now copy to network (overwrites the original) ---
copyfile(localTemp, modelFitFile);
delete(localTemp);

try
    load(modelFitFile); %#ok<NASGU> % final read-back verification
    newBytes = getfield(dir(modelFitFile), 'bytes'); %#ok<GFLD>
    fprintf('  %s/%s: %.2f GB -> %.1f MB (%.1fx), verified\n', projectName, subj, ...
        oldBytes/1e9, newBytes/1e6, oldBytes/newBytes);
    status = 'trimmed';
catch ME
    fprintf('  %s/%s: WARNING -- network copy failed verification after copy (%s)\n', projectName, subj, ME.message);
    status = 'corrupt';
end

end
