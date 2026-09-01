# AGENTS.md

Orientation notes for anyone (human or AI agent) picking up this repository. This
file summarizes what the project is, how the pieces fit together, and the
non-obvious gotchas that have already cost real debugging time.

## What this project is

This repo supports a manuscript in progress:

> **"Local orientation asymmetries in V1 depend on global stimulus properties"**
> Ezzo, Carrasco, Rokers, Winawer — draft at `Support/draft.pdf`

The study measures fMRI (V1) BOLD responses to two families of grating stimuli:
- **"dg" = drifting grating** → the **Cartesian** experiment (horizontal, vertical,
  and two diagonal gratings)
- **"da" = drifting annulus** → the **Polar** experiment (pinwheel, annulus,
  clockwise/counterclockwise spiral gratings)

These two protocol names (`dg`, `da`) are used throughout the MATLAB code and
folder/derivative naming — whenever you see `dg`/`da` as a `projectName`
variable, that's Cartesian vs. Polar.

The central finding: four orientation asymmetries (vertical>horizontal,
oblique>cardinal, radial>tangential, polar-cardinal>polar-oblique) are strongest
when the asymmetry's reference frame matches the global stimulus's reference
frame (Cartesian asymmetries largest for Cartesian gratings, polar asymmetries
largest for polar gratings), even though local orientation/spatial-frequency
content is matched across the two stimulus types at 6° eccentricity.

Read `Support/draft.pdf` for full scientific context before making analysis
decisions — it explains *why* the pipeline is structured the way it is (e.g.
why there's a subject-wise gain correction, why V1 is binned into 8 polar-angle
sectors, why four specific asymmetries were chosen).

### Methods, summarized

- **Observers:** 13 recruited (NYU / NYUAD); 8 completed both the Cartesian
  and polar experiments and are the only ones used in the main-text analyses
  (so contextual comparisons aren't confounded by sample differences) — this
  is the source of the hardcoded `[1:8]` subject subset noted below.
- **Stimuli:** each experiment (`dg`=Cartesian, `da`=polar) had an
  event-related design — 4 stationary-grating conditions, 8 drifting
  conditions, 1 pink-noise condition, 3 s stimulus + 2 s ITI, 52 trials/run,
  50% Michelson contrast, 12.2° aperture radius. **Only the stationary
  conditions are analyzed in this paper**; drifting-stimulus results are left
  for a future study (but drifting trials are still included when fitting the
  GLM, to improve model fit). Cartesian and polar spatial frequency are
  matched only at 6° eccentricity (mid-aperture) — matching breaks down at
  other eccentricities because polar spatial frequency scales with radius.
- **fMRI acquisition/preprocessing:** 3T Siemens Prisma, multiband EPI (TR=1s,
  2mm iso); fMRIPrep + FreeSurfer surface-based pipeline (all analysis is
  vertex-wise on the cortical surface, not voxel-wise).
- **pRF mapping & ROI definition:** a separate retinotopy protocol (bar +
  wedge/ring apertures, stationary- and moving-carrier variants) gives each
  vertex a pRF (x, y, σ) and R² via `vistasoft`. V1 is drawn manually per
  hemisphere, then subdivided into **8 polar-angle sectors** (centered
  0/45/90/.../315°, ±22.5° wide, restricted to 4-8° eccentricity) — this
  "annulus sector" binning is what lets the same local patch be described in
  either Cartesian (e.g. vertical) or polar (e.g. radial) terms, and is the
  origin of the `meanBOLD(pa)` / 8-ROI-column structure everywhere downstream.
  Only vertices with pRF R² ≥ 0.1 are included. Other ROIs (V2, V3, V3A, V3B,
  hV4, MT, MST) are defined the same way but used only in supplementary
  analysis.
- **GLM:** GLMsingle's `GLMestimatesingletrial` (HRF library + noise-pool
  denoising + fractional ridge regression) gives one beta per condition per
  vertex; the pink-noise condition's beta is subtracted from each orientation
  condition's beta as a baseline throughout the paper.
- **Observer gain correction:** each observer's pRF-model amplitude, averaged
  across all 12 retinotopy runs and the same V1 vertices used in the main
  analysis, is used as a single per-observer gain scalar (see *Observer gain
  correction* below for the exact `subjectScale` formula used in code).
- **Two ways the four asymmetries are quantified**, both implemented in
  `04_plot_betaAsymmetries/`:
  1. **Independently** — median % signal change per observer/ROI/condition,
     then bootstrapped (1,000 resamples) confidence intervals on each
     pairwise difference (`plot1_experimentalCond.m`/`plot2_experimentalCond.m`).
  2. **Jointly** — a linear mixed-effects model (`lme1_fit.m`) fits all four
     asymmetries simultaneously (256 datapoints = 8 observers × 4 orientations
     × 8 polar angles per experiment), with each asymmetry coded as a ±1 (or
     0, for the two orientations not in that pairing) predictor, an
     observer-level random intercept, fit by ML, with bootstrapped CIs.
- **Normalization models** (`Models/`): steerable-pyramid energy responses
  (4 orientations × 6 spatial-frequency scales, via `plenoptic`) to the same
  Cartesian/polar stimuli, with cardinal-channel energy artificially doubled
  to simulate V1's cardinal overrepresentation, then passed through 4 divisive
  normalization variants — see `Models/` section below for which succeed.

### Results, summarized

Magnitudes below are the jointly-fit LME estimates (Δ = % signal change,
95% CI), which are what the main text emphasizes:

| Asymmetry | Cartesian (`dg`) experiment | Polar (`da`) experiment |
|---|---|---|
| horizontal − vertical | **−0.55%** [−0.65, −0.44] | −0.22% [−0.34, −0.11] (<half the size) |
| cardinal − oblique | **−0.22%** [−0.32, −0.12] | −0.03% [−0.06, 0.00] (essentially gone) |
| radial − tangential | 0.10% [0.05, 0.17] | **0.15%** [0.03, 0.25] (~1.5x larger) |
| polar cardinal − polar oblique | 0.04% [0.01, 0.07] | 0.04% [−0.03, 0.11] (same magnitude, not reliable here) |

The pattern (bolded = the larger/more reliable of the pair): **Cartesian
asymmetries are largest for Cartesian stimuli; the radial asymmetry is larger
for polar stimuli; the polar-cardinal asymmetry is unaffected by context but
also weak/unreliable in both.** In short, an asymmetry is strongest when its
own reference frame matches the global stimulus's reference frame.

The more surprising/counterintuitive result is the **sign** of the Cartesian
effects: horizontal and cardinal orientations are known to have *more*
V1 neural representation, natural-scene frequency, and (for cardinal)
behavioral sensitivity than vertical/oblique — yet BOLD is **lower**, not
higher, for horizontal and cardinal. Polar asymmetries go the other way:
radial and polar-cardinal orientations, which are also neurally
overrepresented/behaviorally advantaged, show **higher** BOLD. The paper's
interpretation is that Cartesian-preferred orientations, being more
numerous/more mutually-similar in the surround of a spatially coherent
Cartesian grating, are suppressed via context-dependent normalization more
than they're boosted by their raw overrepresentation — an effect that a
plain "more neurons → more BOLD" model would get backwards.

This is what the `Models/` normalization simulation is testing: whether
divisive normalization can turn overrepresentation into net suppression.
Untuned and orientation-tuned normalization (models 1-2) fail — cardinal
responses stay elevated in the model regardless of context. Superlinear
(squared-suppression) and anisotropy-based normalization (models 3-4)
succeed: cardinal responses are suppressed relative to oblique specifically
for Cartesian, not polar, stimuli — mirroring the fMRI asymmetry reversal.
As of this draft, this is presented only as a candidate mechanism (Fig. 9 in
the draft is still a placeholder), not a settled explanation.

## Repository layout

```
AnalysisCode/     MATLAB pipeline: GLM fitting -> ROI/polar-angle asymmetry
                  stats -> linear mixed-effects models -> plotting.
Models/           Python/Jupyter: image-computable normalization models
                  (steerable-pyramid energy + divisive normalization) that
                  attempt to explain the *sign inversion* of Cartesian
                  asymmetries. Implements the paper's Methods subsection
                  "Normalization as a description of contextual suppression".
ExperimentCode/   Psychtoolbox stimulus presentation code (not analysis).
Simulations/      Standalone motion-vector demo scripts.
Reproduction/     Scripts for pulling data from a remote server.
Support/          draft.pdf (the manuscript) and a summary CSV table.
FIGURES.md        Script-by-script data lineage for every manuscript figure
                  (polar, pairwise, subjectwise, context-comparison,
                  directions-per-location, across-cortical-area, TTA, and
                  model-vs-observed time-series plots) -- read this before
                  touching any figure-generating script.
```

### AnalysisCode/ subfolder pipeline order

Numeric prefixes indicate pipeline stage, but they are not strictly run
top-to-bottom for every analysis — the two chains that matter most this year:

1. **`01_process_singlesubjectGLM/`** — per-subject GLMsingle fit
   (`main_singlesub.m`), produces per-subject `results.mat` contrast files.
2. **`01_calculate_observer_gain/`** — computes each observer's mean pRF gain
   from the retinotopy protocols (`dg_computeGain.m`, `dg_compareGainROI.m`),
   used later to correct for cross-subject differences in raw BOLD scale (see
   *Observer gain correction* below). Output: `gainSummary.mat`/`.csv` under
   `derivatives/summaryTables/`.
3. **`03_process_groupBetas/meanWithinLabel.m`** — reads every subject's
   per-ROI FreeSurfer labels + GLM betas, bins by polar angle, and writes the
   group-level matrices `meanBOLD(pa)`/`medianBOLD(pa)` that everything
   downstream consumes. **This is the one script that actually defines the ROI
   dimension's order** (see ROI ordering gotcha below).
4. **`04_plot_betaAsymmetries/`** — the main analysis + figures. See
   `FIGURES.md` for the full script-by-script breakdown of every figure type
   and exactly where its data comes from; the short version is that
   `plot_NeuralAsymmetries.m` orchestrates most of it, `lme1_fit.m` is the
   original (legacy) joint model, and `fitAsymmetryRegression.m`'s cache has
   since become the shared source for most current figures.

Data (`meanBOLD*.mat`, GLM results, FreeSurfer labels, etc.) lives **outside
this git repo**, on a network volume at
`/Volumes/Vision/UsersShare/Rania/Project_dg/data_bids/`. The repo only holds
code and small support tables.

## Known gotchas (read before editing the MATLAB pipeline)

### ROI ordering is positional, not declarative
`AnalysisCode/general/jsons/ROIS_ALL.json` defines the ROI list (`filename`,
`plotname`, `index`). **The `index` field is only meaningful if
`meanWithinLabel.m` was run with that exact ROI list, in that exact order** —
it does not derive from an existing canonical mapping. `meanWithinLabel.m`
writes `meanBOLD(ci, ri, si)` etc. using the ROI list's *position* (`ri`), so
the Nth entry in the JSON becomes the Nth column of every saved
`meanBOLD*`/`medianBOLD*` matrix. Downstream, `plot1_experimentalCond.m` /
`plot2_experimentalCond.m` correctly re-derive the right column via
`projectSettings.roi_idx{ri}`; `lme1_fit.m` was fixed to do the same (previously
it used the raw loop counter, which only ever worked by coincidence). **If you
add/remove/reorder ROIs in `ROIS_ALL.json`, you must re-run
`meanWithinLabel.m`** for every affected project — there is no way to
reconcile a changed ROI list against already-saved matrices. The current
`ROIS_ALL.json` (8 ROIs: V1, V2, V3, V3A, V3B, hV4, pMT, pMST) was restored
from a truncated version (see `ROIS_ALL_2026-08-05.json` backup in the same
folder) — if you need the original 9-ROI canonical order (which included
`hMTcomplex`), it's recoverable from git history (`ROIS_ALL copy.json`,
commit `e37e1be` and earlier).

### Observer gain correction
Each subject's BOLD values get multiplied by `subjectScale(i) = groupGain /
gain_i`, where `gain_i` is that subject's mean pRF gain (from
`retrieveObserverGainWeights.m` / `gainSummary.mat`) and `groupGain =
exp(mean(log(gainWeights)))` (**geometric**, not arithmetic, mean of the
whole subject set — implemented via `exp(mean(log(.)))` rather than
`geomean()` to avoid a Statistics/ML Toolbox dependency). This:
- **down-weights high-gain subjects, up-weights low-gain subjects** — the
  substantive correction.
- multiplying back by `groupGain` afterward is a **pure unit-restoration
  step** — mathematically inert for every inferential statistic (t/F-stats,
  p-values, CIs); any positive constant here rescales every number by the
  same factor and changes no significance conclusion.
- **why geometric mean specifically**: `subjectScale` is itself a
  multiplicative scale factor (a ratio), not an additive quantity. The
  *geometric* mean of the resulting `subjectScale` values across observers
  is exactly 1 only when `groupGain` is the geometric mean of `gainWeights`
  — i.e. this is the choice under which the restoration step is
  magnitude-neutral in the sense that matches what `subjectScale` actually
  is. (Arithmetic mean of `gainWeights` is the corresponding exact choice if
  you instead judge neutrality via the *arithmetic* mean of `subjectScale`;
  neither is more "correct" in the abstract, but geometric mean is the
  natural fit given `subjectScale` is a ratio.) This only affects the
  reported magnitude of the numbers, never any inferential conclusion — see
  previous bullet.
- For `dg` (13 subjects total), only the **first 8** are used in the main
  analysis scripts (`lme1_fit.m`, `plot_NeuralAsymmetries.m`) — hardcoded
  index `[1:8]` — to match the 8 subjects who completed the `da` experiment.
  `lme2_ploteachDirLoc.m` has no subject-selection logic of its own; it
  inherits whatever subset `lme1_fit.m` used when it saved `modeldata.mat`.

### MATLAB working-directory requirement
`AnalysisCode/general/setup_user.m` requires MATLAB's current folder to be
`AnalysisCode/` itself (or a path containing it) — it explicitly checks
`pwd` and errors otherwise (`"AnalysisCode folder containing setup.json must
be the parent folder."`). Scripts that call `setup_user()`
(`meanWithinLabel.m`, GLM-fitting scripts) must be run with that working
directory. Scripts that don't call `setup_user()` (e.g. `lme1_fit.m`) still
rely on `addpath(genpath(pwd))`, so make sure `AnalysisCode/` (and all
subfolders) are on the MATLAB path before running anything, regardless of
which tool/method is used to launch the script.

### da/sub-wlsubj124: runs 1 and 2 share a duplicated trial sequence

Confirmed 2026-08-17: for **da / sub-wlsubj124 only**, `S02_design_Run1.mat`
and `S02_design_Run2.mat` in `experimentalOutput/da/02/` have byte-identical
`expDes.trialMat` content (52 trials, multiple randomized conditions —
confirmed not a coincidence), and this is already baked into
`matrices_onset{1}`/`matrices_onset{2}` in that subject's `rawInfo.mat`. No
other subject/project combination has this issue (checked all 16
subject×project combos pairwise).

**This is already correctly resolved at the GLM level — do not re-run
GLMsingle or "fix" `modelOutput.mat`/`derivedModelFit.mat` for this.** A
cross-run data-vs-model correlation check (real BOLD data for both runs
correlates best with their own run's model, and both runs correlate with
each other's model/data well above the cross-run baseline) shows the same
trial sequence really was presented twice in the scanner and captured
correctly — this is a repeated-stimulus quirk in the experiment itself, not
a mislabeling/misalignment bug. `format_desmats.m` reads each run strictly
by index (`Run<i>` design paired with BIDS `run-<i>` BOLD), and that
index-for-index correspondence is intact; the two runs are just not
independent trials of each other.

**It only matters for pipelines that average across runs/subjects as if
every run were an independent repeat** (e.g. the cross-subject group-average
run time series) — averaging both runs in as independent would double-count
that one trial sequence relative to every other (genuinely unique) run.
`findDuplicateDesignRuns.m` (`AnalysisCode/02_ttave/`) detects this
generically (pairwise `isequaln` on `matrices_onset` per subject) and is
wired into `run_groupAverageRunTimeseries.m`, which excludes the later run
of any duplicate pair from the group average. Per-run/per-subject-only
outputs (e.g. `run_runTimeseries.m` / `plot_runTimeseries.m`, which plot
each run separately and never average across runs) are unaffected and need
no changes.

## Models/ (Python) — normalization simulations

`model_DriftingGratings.ipynb` implements the paper's Methods subsection
*"Normalization as a description of contextual suppression"* and its
(currently incomplete) Results counterpart *"Normalization model describes
context-dependent suppression"* (draft.pdf, Methods p.13, Results p.25 — the
results section still has a placeholder figure, "Fig 9. Image computable
model? ... schematic?").

Pipeline: build Cartesian/polar stimulus images → run through a Steerable
Pyramid (4 orientations × 6 spatial frequencies, via the `plenoptic` library)
→ artificially double the cardinal-orientation channels' energy (simulating
V1's overrepresentation of cardinal-tuned neurons) → test whether various
divisive-normalization formulations can turn that overrepresentation into a
*suppressed* net response (matching the "inverted" BOLD asymmetry seen for
cardinal-vs-oblique in the Cartesian experiment, but not the polar one).

Models implemented in `helper_functions/utils_image.py`:
1. `div_normalization(..., tuned=False)` — untuned divisive normalization.
2. `div_normalization(..., tuned=True)` — feature/orientation-tuned normalization.
3. `div_normalization(..., tuned=False, q_exp=2)` — untuned + exponentiated
   (superlinear) suppression.
4. `normalization_byAnisotropy(...)` — suppression scales with the standard
   deviation of energy across orientation channels (adapted from Fang et al.).
5. `normalization_byStimHomogeneity(...)` — **not yet in the manuscript text**;
   a newer, exploratory model based on cosine similarity between center/surround
   feature vectors. Treat as unpublished/in-progress.

Per the draft: models 1–2 fail to reproduce the suppression; models 3–4
succeed, and only for Cartesian (not polar) stimuli — mirroring the fMRI data.
Beyond the manuscript's 4 named models, the notebook's own internal numbering
(see below) also includes an anisotropy variant that skips spatial collapse
and the not-yet-published `normalization_byStimHomogeneity` model.

**The notebook's title cell (cell 0) is a running DONE / open-questions / TO-DO
log kept by the author — read it first, it is more current than any static
summary of the code.** Highlights as of this writing: radial-orientation
overrepresentation has also been added (in addition to cardinal) with
behaviorally-derived weights; there's an open question about how to choose the
4D Gaussian convolution bandwidth for the untuned/tuned models; and the TO-DO
list includes adding a polar-cardinal asymmetry version, a literature review of
suppressive-exponent normalization, and testing predictions for contrast/
adaptation/natural-scene manipulations.

The notebook otherwise has no markdown documentation of the code itself and
ends in empty cells — it is a working draft, not a finished analysis. See the
installation instructions at the top of the notebook and in `README.md`
before running it.

### Performance: the normalization models are GPU-accelerated, the pyramid step is not

Profiling (2026-08-07) found the notebook's ~12-minute runtime was ~99% spent
in `div_normalization` (Models 1–3) and `normalization_byStimHomogeneity`
(Model 6) — specifically in their `scipy.ndimage.gaussian_filter` calls with
a ~90px-radius kernel — not in the steerable-pyramid step, which is <1% of
total time. `utils_image.py` now has `gaussian_filter_gpu()` (backed by
`get_fast_device()`: CUDA > MPS > CPU), a drop-in torch-based replacement for
`scipy.ndimage.gaussian_filter` used inside those two functions, verified to
match scipy's output to ~1e-7 relative error. This cut total runtime from
~715s to ~44s on this Intel Mac (MPS backend) — a ~16x speedup — with no
change to model behavior. `div_normalization` and
`normalization_byStimHomogeneity` both take an optional `device=` kwarg
(defaults to auto-detect).

Implementation note: `gaussian_filter_gpu` cannot just move the SciPy call's
axis-length-based reasoning to torch — SciPy's per-axis cost scales with
*kernel radius*, not axis length, so even the tiny 4–6-element SF/orientation
axes are expensive to filter at this sigma. It instead builds a dense `(L,
L)` linear operator per axis (matching scipy's default half-sample-symmetric
`mode='reflect'` boundary behavior, including the periodic folding needed
when the kernel radius exceeds the axis length) and applies it via matmul —
padding-based conv1d was tried first and blew MPS's buffer limit (~9GB) when
padding a small axis while other axes stayed huge.

The steerable pyramid step is intentionally left on CPU regardless of
`DEVICE`: `plenoptic`'s `SteerablePyramidFreq` calls `torch.fft.fft2` on
complex tensors internally, and PyTorch's MPS backend does not support
complex dtypes for FFT ops as of `torch==2.2.2` (confirmed by direct test:
`torch.fft.fft2` on a real MPS tensor raises "Trying to convert ComplexFloat
to the MPS backend but it does not have support for that dtype"). This may
be fixed in newer torch releases (worth re-checking on Apple Silicon setups
pinned to `torch==2.9.0`, e.g. via `Models/requirements.txt`), but isn't
worth chasing further given the step is <1% of runtime either way.

### Dependency pins are platform-specific — do not `pip install` unpinned

`Models/requirements.txt` (Apple Silicon / arm64) and
`Models/requirements-intel-mac.txt` (Intel / x86_64) exist because a plain
`pip install plenoptic torch ...` silently grabs versions that break the
notebook, in two independent ways that were both hit and diagnosed in this
repo's history:

1. **`plenoptic>=2.0` removed the `plenoptic.simulate` module** that this
   notebook imports (`from plenoptic.simulate import SteerablePyramidFreq`) —
   this fails with `ModuleNotFoundError: No module named 'plenoptic.simulate'`,
   not an obvious version-mismatch error. Both requirements files pin
   `plenoptic==1.3.1`, matching a real working environment (confirmed via
   `pip show`/`pkgutil.iter_modules` that `2.1.0` genuinely lacks this
   submodule; `1.3.1` has it).
2. **PyTorch dropped Intel-macOS (x86_64) wheels after the `2.2.x` series.**
   On an Intel Mac, no `torch` version at all is installable for Python
   ≥3.12, and the newest installable version for Python 3.11 is `2.2.2`.
   `torch==2.2.2` was compiled against NumPy 1.x's ABI, so installing NumPy
   2.x alongside it produces `UserWarning: Failed to initialize NumPy:
   _ARRAY_API not found` and real crash risk (NumPy's own message: "may
   crash"). This further constrains Intel Macs to `numpy<2` (verified pin:
   `1.26.4`) and, transitively, to an older `opencv-python` release
   (`4.9.0.80`) since `opencv-python>=4.10ish` requires `numpy>=2`, which
   would directly conflict with the `numpy<2` pin. None of this applies on
   Apple Silicon, which has current wheels for everything.

If you see either failure mode again (missing `plenoptic.simulate`, or the
NumPy ABI warning cascading into an unrelated `ModuleNotFoundError`), the fix
is almost always "an unpinned `pip install` grabbed a too-new package" — check
`pip show <package>` for the installed version against the pinned
`requirements*.txt`, don't just reinstall the same unpinned command again.

When updating these pins (e.g. to allow a newer `plenoptic` once its API
stabilizes), verify on **both** architectures before merging — `python3 -c
"import platform; print(platform.machine())"` distinguishes them
(`arm64` vs `x86_64`), and a clean import check looks like:
```python
import torch, numpy
from plenoptic.simulate import SteerablePyramidFreq  # must not raise
```
with no NumPy ABI warning printed.

## When picking up new work here

- Check `Support/draft.pdf` for the current state of the manuscript text
  before assuming what's "done" — sections can be explicitly marked
  incomplete (placeholder figures, missing quantitative comparisons).
- Anything under `AnalysisCode/04_plot_betaAsymmetries/` that reads
  `meanBOLD*`/`medianBOLD*` assumes those matrices already exist for the
  current ROI list; regenerate via `meanWithinLabel.m` first if in doubt.
- Prefer extending the existing per-project (`dg`/`da`) branching pattern
  (`if strcmp(projectName, 'dg') ... elseif strcmp(projectName, 'da')`) rather
  than introducing a new convention — it's used consistently across every
  script in `04_plot_betaAsymmetries/`.
