function subjTraces = computeSubjectRunTraces_GLMsingle(bidsDir, projectName, subj, ses, projectSettingsBase, ...
    stimdur_s, tr_s, matrices_onset, nRunsCheck, roiName)
% computeSubjectRunTraces_GLMsingle - per-run, per-subject, continuous ROI
% time series, BOTH the raw observed trace and the GLMsingle-based
% predicted trace, for one subject, separately for EACH of the 8
% canonical 45-deg polar-angle wedges (0:45:315). Vertex selection (ROI +
% eccentricity 4-8 deg + pRF R^2 >= 0.1) is identical to
% computeSubjectTraces.m -- this is that function's counterpart, extended
% to ALSO return the GLMsingle model prediction (not just the observed
% data), and to return all 8 wedges from ONE pass over the raw data and
% GLMsingle reconstruction (rather than reloading/re-reconstructing once
% per wedge, which would be 8x the I/O and compute for no reason -- the
% wedge assignment only affects which vertices get reduced together at
% the very last step).
%
% Both sides use the SAME conventions createTTaveTable.m's per-condition
% TTA extraction already established (computeConditionTTA_rawdata /
% computeConditionTTA_modelfit), just kept as full continuous per-run
% traces instead of reduced to peristimulus epochs:
%   - observed: raw BOLD converted to %-signal-change FIRST using
%     GLMsingle's OWN meanvol reference, THEN polynomial-denoised (same
%     degree formula GLMestimatesingletrial.m itself uses). This order
%     matters: denoising (which removes the degree-0/DC term) BEFORE
%     dividing by a fixed meanvol reference produces a spurious ~-100%
%     constant offset -- found and fixed 2026-09-01, see git history for
%     the numeric trace before/after.
%   - predicted: GLMsingle's actual fitted single-trial betas
%     (results.modelmd) convolved with each vertex's OWN GLMsingle-
%     selected HRF (results.HRFindex, into getcanonicalhrflibrary) via
%     GLMpredictresponses, batched by HRF-library index across the WHOLE
%     ROI's vertices at once (one call per unique HRF, covering every
%     run and every wedge). This is GLMsingle's real fit -- no fresh/
%     simplified refit anywhere.
%
% Both are reduced to one scalar per TR per run PER WEDGE via nanmedian
% across that wedge's own selected vertices (an exact median over that
% wedge's vertex subset of the full-ROI per-vertex reconstruction -- not
% an approximate combination of other groups' medians).
%
% <roiName> defaults to 'V1' if omitted/empty.
%
% Returns subjTraces with fields:
%   .wedgeAngles     1 x 8, the canonical wedge centers (0:45:315)
%   .nIncludedVerts  1 x 8, vertex count per wedge
%   .observed        1 x 8 cell, each a 1 x nRunsUse cell of row vectors
%   .predicted        1 x 8 cell, each a 1 x nRunsUse cell of row vectors
%
% CACHING (added 2026-09-01): this is the single most expensive step every
% caller (runRawTimeseriesPermutationTest.m, plotRawTimeseriesWithinAcross_
% GLMsingle.m, plotTTA_GLMsingle.m) repeats identically for the same
% (projectName, subj, roiName) -- raw BOLD .mgh loading + the full
% GLMsingle reconstruction, ~4-5s/subject. Result is cached to
% <bidsDir>/derivatives/summaryTables/GLMsingleRunTracesCache/<projectName>/
% <roiName>/<subj>_nRuns<nRunsCheck>.mat and reused on subsequent calls --
% the actual computation below is untouched, this only wraps it with a
% load-if-present / save-after-compute check.

if nargin < 10 || isempty(roiName)
    roiName = 'V1';
end

cacheDir = fullfile(bidsDir, 'derivatives', 'summaryTables', 'GLMsingleRunTracesCache', projectName, roiName);
cacheFile = fullfile(cacheDir, sprintf('%s_nRuns%d.mat', subj, nRunsCheck));
if isfile(cacheFile)
    loaded = load(cacheFile, 'subjTraces');
    subjTraces = loaded.subjTraces;
    return
end

wedgeAngles = 0:45:315;
nWedges = numel(wedgeAngles);

hSize = get_surfsize(subj);
nRuns_subj = numel(matrices_onset);
nRunsUse = min(nRuns_subj, nRunsCheck);

projectSettings = projectSettingsBase;
projectSettings.projectName = projectName;
projectSettings.subject = subj;
projectSettings.bidsDir = bidsDir;
projectSettings.retFolder = 'prfvista_mov';
projectSettings.ses = ses;
projectSettings.minECC = 4;
projectSettings.maxECC = 8;
projectSettings.minVAREXP = .1;
projectSettings.stimdur_s = stimdur_s;
projectSettings.tr_s = tr_s;
projectSettings.roiName = roiName;
projectSettings.polarAngleBinWidth = 45;

primaryROIvertices = getROIidxs(subj, projectSettings.roiName, hSize);
surfaceROI = nan(sum(hSize), 1);
surfaceROI(primaryROIvertices) = 1;

filteredPrfBins = retriveRetData(projectSettings); % 1 x nFullSurface, wedge center or NaN per vertex
surfaceSelection = surfaceROI .* (~isnan(filteredPrfBins))'; % full-ROI, any wedge, ecc/R^2-valid

selVertIdx = find(~isnan(surfaceSelection)); % the FULL ROI's included vertices, all wedges pooled
wedgeOfSelVert = filteredPrfBins(selVertIdx)'; % which of the 8 wedges each selected vertex belongs to

nIncludedVerts = arrayfun(@(w) sum(wedgeOfSelVert == w), wedgeAngles);

%% load GLMsingle output once, reused by both the observed and predicted sides

derivativesFolder = fullfile(bidsDir, 'derivatives', sprintf('%sGLM', projectName), 'hRF_glmsingle', subj, ses);
modelData = load(fullfile(derivativesFolder, 'modelOutput.mat'));
results = modelData.modelOut{1,4};
designInfo = modelData.designSINGLE;

%% ----- observed: raw BOLD -> %sc (GLMsingle meanvol) -> polynomial-denoise -----

meanvol = results.meanvol(selVertIdx, 1);

run_ = 1:nRunsUse;
datafiles = load_data(bidsDir, projectName, 'fsnative', '.mgh', subj, ses, run_);

observedPerVertex = cell(1, nRunsUse); % each: nSel x nTRsRun
for r = 1:nRunsUse
    df = datafiles{r}(selVertIdx, :); % nSel x nTRsRun

    psc = ((df ./ meanvol) - 1) * 100; % %sc FIRST (see header note)

    nTRsRun = size(df, 2);
    maxpolydeg = round(((nTRsRun * tr_s) / 60) / 2);
    pmatrix = constructpolynomialmatrix(nTRsRun, 0:maxpolydeg);
    polymatrix = projectionmatrix(pmatrix);
    psc = (polymatrix * psc')'; % denoise the %sc trace, then back to voxels x time

    observedPerVertex{r} = psc;
end

%% ----- predicted: GLMsingle betas x per-vertex HRF via GLMpredictresponses, whole ROI at once -----

modelmd = results.modelmd;            % nVertFull x 1 x 1 x nTrials, in %
HRFindex = results.HRFindex;          % nVertFull x 1
nTrials = size(modelmd, 4);

hrflibrary = getcanonicalhrflibrary(stimdur_s, tr_s); % nHRFs x time

designRuns = designInfo.designSINGLE; % 1 x nRuns cell, TR x nTrials
nRunsDesign = numel(designRuns);
numtimepoints = cellfun(@(x) size(x,1), designRuns);
nRunsPredict = min(nRunsUse, nRunsDesign);

hrfii_sel = HRFindex(selVertIdx);
uniqueHRFs = unique(hrfii_sel)';

vertRowsByRun = cell(1, nRunsPredict);
for r = 1:nRunsPredict
    vertRowsByRun{r} = {};
end
vertOrderAccum = [];

for hh = uniqueHRFs
    selPos = find(hrfii_sel == hh);      % rows (1..nSel) in selVertIdx
    voxFull = selVertIdx(selPos);        % absolute vertex ids into modelmd

    betas = reshape(modelmd(voxFull,1,1,:), numel(voxFull), nTrials); % nVoxSubset x nTrials
    hrf = hrflibrary(hh,:)';             % time x 1

    mf = GLMpredictresponses({hrf, betas}, designRuns, tr_s, numtimepoints, 1); % 1 x nRuns cell, one call covering every run

    for r = 1:nRunsPredict
        signalMatrix = mf{r};
        if size(signalMatrix, 1) ~= numel(voxFull)
            signalMatrix = signalMatrix(voxFull, :);
        end
        vertRowsByRun{r}{end+1} = signalMatrix; %#ok<AGROW>
    end
    vertOrderAccum = [vertOrderAccum; selPos(:)]; %#ok<AGROW> % selVertIdx-relative positions, in the order rows were appended

    clear mf
end

predictedPerVertex = cell(1, nRunsUse); % each: nSel x nTRsRun, rows in selVertIdx order (via vertOrderAccum)
for r = 1:nRunsPredict
    fullMat = cell2mat(vertRowsByRun{r}(:)); % rows in vertOrderAccum order
    reordered = nan(numel(selVertIdx), size(fullMat,2));
    reordered(vertOrderAccum, :) = fullMat;
    predictedPerVertex{r} = reordered;
end
for r = (nRunsPredict+1):nRunsUse
    predictedPerVertex{r} = nan(numel(selVertIdx), size(observedPerVertex{r}, 2));
end

%% ----- reduce to one trace per TR per run, PER WEDGE -----

observed = cell(1, nWedges);
predicted = cell(1, nWedges);
for w = 1:nWedges
    wedgeRows = find(wedgeOfSelVert == wedgeAngles(w));
    observed{w} = cell(1, nRunsUse);
    predicted{w} = cell(1, nRunsUse);
    for r = 1:nRunsUse
        if isempty(wedgeRows)
            observed{w}{r} = nan(1, size(observedPerVertex{r}, 2));
            predicted{w}{r} = nan(1, size(predictedPerVertex{r}, 2));
        else
            observed{w}{r} = nanmedian(observedPerVertex{r}(wedgeRows, :), 1);
            predicted{w}{r} = nanmedian(predictedPerVertex{r}(wedgeRows, :), 1);
        end
    end
end

subjTraces.wedgeAngles = wedgeAngles;
subjTraces.nIncludedVerts = nIncludedVerts;
subjTraces.observed = observed;
subjTraces.predicted = predicted;

if ~isfolder(cacheDir), mkdir(cacheDir); end
save(cacheFile, 'subjTraces');

end
