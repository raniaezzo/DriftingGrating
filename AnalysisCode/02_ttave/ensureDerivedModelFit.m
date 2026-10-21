function modelFitFile = ensureDerivedModelFit(bidsDir, projectName, subj, ses, hRF_setting, stimdur_s, tr_s)
% ensureDerivedModelFit - generate derivedModelFit.mat for a subject/project
% if it doesn't already exist, using the GLMsingle full model (modelOut{4})
% and GLMpredictresponses to reconstruct the predicted time series per run.
% Same logic as plottingGLMsingle.m, factored out so run_runTimeseries.m can
% call it directly instead of requiring a separate manual step.
%
% Only vertices inside the union of the ROIs listed in
% general/jsons/ROIS_ALL.json are stored (~13x smaller on disk than the
% full cortical surface) -- every actual analysis in this project masks to
% one specific ROI within that union, so nothing outside it is ever used.
% Load with loadDerivedModelFit.m, which reconstructs the full-surface
% array (NaN outside the stored ROI union) that downstream code expects.

derivativesFolder = fullfile(bidsDir, 'derivatives', sprintf('%sGLM',projectName), ...
    sprintf('hRF_%s', hRF_setting), subj, ses);
modelFitFile = fullfile(derivativesFolder, 'derivedModelFit.mat');

if exist(modelFitFile, 'file')
    return
end

fprintf('  derivedModelFit.mat not found for %s/%s -- generating (this can take several minutes)...\n', ...
    subj, projectName);

load(fullfile(derivativesFolder,'modelOutput.mat'), 'modelOut');
designSINGLE = load(fullfile(derivativesFolder, 'DESIGNINFO.mat'));

unique_hrfs = unique(modelOut{4}.HRFindex);
numtrials = size(modelOut{4}.modelmd, 4);
hrflibrary = getcanonicalhrflibrary(stimdur_s, tr_s);
n_runs = numel(designSINGLE.designSINGLE);
n_vertices_total = numel(modelOut{4}.HRFindex);
numtimepoints = zeros(1, n_runs);
for n = 1:n_runs
    numtimepoints(n) = size(designSINGLE.designSINGLE{n}, 1);
end

modelfit = cell(1, n_runs);
for r = 1:n_runs
    modelfit{r} = nan(n_vertices_total, numtimepoints(r));
end

for hh = unique_hrfs'
    vox_mask = (modelOut{4}.HRFindex == hh);
    hrf = hrflibrary(hh,:)';
    betas = reshape(modelOut{4}.modelmd .* vox_mask, [], numtrials);
    mf = GLMpredictresponses({hrf, betas}, designSINGLE.designSINGLE, tr_s, numtimepoints, 1);
    for r = 1:numel(mf)
        valid = ~all(mf{r} == 0, 2);
        modelfit{r}(valid,:) = mf{r}(valid,:);
    end
end

% trim to the union of vertices in the project's ROI set before saving --
% see header comment. nVerticesTotal is kept so loadDerivedModelFit.m can
% reconstruct a full-surface-length array on load.
hSize = get_surfsize(subj);
roisJson = jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')), '..', 'general', 'jsons', 'ROIS_ALL.json')));
roiNames = {roisJson.ROIs.filename};
roiVertexIdx = [];
for i = 1:numel(roiNames)
    roiVertexIdx = union(roiVertexIdx, getROIidxs(subj, roiNames{i}, hSize));
end
roiVertexIdx = sort(roiVertexIdx(:));

modelfitROI = cell(1, n_runs);
for r = 1:n_runs
    modelfitROI{r} = modelfit{r}(roiVertexIdx, :);
end
nVerticesTotal = n_vertices_total; %#ok<NASGU>

% safe save: local disk first, verify, then copy to network -- direct
% in-place network writes of large .mat files have repeatedly corrupted
% files earlier in this project's session history.
localTemp = fullfile(tempdir, sprintf('derivedModelFit_%s_%s_generate.mat', projectName, subj));
save(localTemp, 'modelfitROI', 'roiVertexIdx', 'nVerticesTotal', '-v7.3');
checkS = load(localTemp);
if ~isequal(size(checkS.modelfitROI), size(modelfitROI)) || numel(checkS.roiVertexIdx) ~= numel(roiVertexIdx)
    delete(localTemp);
    error('ensureDerivedModelFit:localSaveVerificationFailed', ...
        'local save verification failed for %s/%s -- network file not written', projectName, subj);
end
copyfile(localTemp, modelFitFile);
delete(localTemp);
load(modelFitFile); % final read-back verification; throws if the network copy is bad

fprintf('  saved %s (%d of %d vertices)\n', modelFitFile, numel(roiVertexIdx), n_vertices_total);

end
