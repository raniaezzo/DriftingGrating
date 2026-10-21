function modelfit = loadDerivedModelFit(modelFitFile)
% loadDerivedModelFit - load derivedModelFit.mat and return the same
% full-surface-length {1 x nRuns} cell array of [nVertices x nTimepoints]
% predicted-response matrices that every downstream script (ttave_compute.m,
% runTimeseries_compute.m, ...) expects -- regardless of whether the file
% on disk is the original full-surface format or the ROI-trimmed format.
%
% ROI-trimmed files only store rows for vertices inside the union of the
% ROIs in general/jsons/ROIS_ALL.json (see ensureDerivedModelFit.m); every
% vertex outside that union is reconstructed here as NaN. Since all actual
% analyses mask to one specific ROI within that union (and NaN * anything
% = NaN), this reconstruction is mathematically identical, for any such
% analysis, to loading the original full-surface file -- verified directly
% against the original derivedModelFit.mat for sub-0037/dg (bit-for-bit
% identical downstream V1 trace, max abs diff = 0) before this trimmed
% format was adopted.

S = load(modelFitFile);

if isfield(S, 'modelfit')
    % original full-surface format -- nothing to reconstruct
    modelfit = S.modelfit;
    return
end

% ROI-trimmed format
nRuns = numel(S.modelfitROI);
modelfit = cell(1, nRuns);
for r = 1:nRuns
    full = nan(S.nVerticesTotal, size(S.modelfitROI{r}, 2));
    full(S.roiVertexIdx, :) = S.modelfitROI{r};
    modelfit{r} = full;
end

end
