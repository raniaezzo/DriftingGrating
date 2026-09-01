# FIGURES.md

Script-by-script data lineage for every figure-generating script in this
repo's MATLAB pipeline. See `AGENTS.md` first for project background,
methods, and repo-wide gotchas — this file assumes that context and only
covers *how each figure gets made and what data it's made from*.

Every asymmetry-related figure ultimately traces back to one of two
group-level arrays built by `meanWithinLabel.m`
(`AnalysisCode/03_process_groupBetas/`) from raw per-vertex GLMsingle betas:
**`meanBOLDpa`** (contrast × 8 polar-angle locations × roi × subject,
restricted to 4-8° eccentricity, pRF R²≥0.1) or the newer **`meanBOLDnative`**
(see *Across-cortical-area plots* below). The per-project cached fits
(`fitAsymmetryRegression.m` / `fitAsymmetryRegressionNative.m`) sit between
the raw arrays and most of the plotting scripts, so most figures never touch
`meanBOLD*` directly.

All paths below are relative to `AnalysisCode/` unless noted.

## Polar plots
**`04_plot_betaAsymmetries/plot1_experimentalCond.m`** (called from
`plot_NeuralAsymmetries.m`) — polar-coordinate plots of pro/con condition
values around the visual field. Data: `meanBOLDpa`, either the 4 fixed
"native" contrast indices (`retrieveProConIdx.m`, contrasts 26-29 =
`s0/s90/s45/s135_v_b`) or a location-recombined "derived" matrix built by
`compute_derivativeDirections.m` (da's derived terms need
`deriveLocalMotionfromUVM.m`'s rotation; dg's don't — raw dg codes are
already Cartesian). Per-subject values are gain-corrected
(`retrieveObserverGainWeights2.m`, ROI-specific). The smooth harmonic overlay
(cos(2θ)/cos(4θ) curves) reads `fitAsymmetryRegression.m`'s cached
coefficients directly, not `meanBOLDpa` — it's the same underlying fit as the
across-cortical-area plots, just rendered as a continuous curve instead of a
bar/point.

## Pairwise plots
**`04_plot_betaAsymmetries/plot2_experimentalCond.m`**, default `model`
argument (`'pairwise'`) — one pro dot + one con dot per asymmetry,
group-level, at a fixed physical axis scale. Data: same `meanBOLDpa`-derived
native/derived pro/con values as the polar plots, but averaged EQUALLY across
all 8 polar-angle locations (so every visual-field location counts the same
regardless of vertex count), gain- and precision-weighted
(`retrieveObserverPrecisionWeights.m` — currently a uniform-1s placeholder).
CI is a per-condition bootstrap over subjects.

## Subjectwise plots
**`04_plot_betaAsymmetries/plot2_experimentalCond.m`**, called with
`'model','subjectwiseDiff'` — same underlying per-subject
`meanBOLDpa`-derived pro/con values as the pairwise mode above (same
gain-weighting, same equal-location-averaging), but displayed differently:
each subject's own pro-minus-con difference as a jittered dot at x=1, plus
the group mean ± 95% CI at x=2, instead of two separate pro/con dots. Fixed
y-range `[-0.8, 0.4]` across every asymmetry and project.

## Context comparison plots (dg vs. da)
**`04_plot_betaAsymmetries/fitAsymmetryRegression_dgVsDa.m`** is the shared
data source for both scripts below: reads `meanBOLDpa` directly (independent
of `fitAsymmetryRegression.m`'s per-project cache), restricted to the 7
subjects who completed *both* dg and da, and fits all 4 asymmetries for dg
and da in one **paired** bootstrap (the same resampled subject-index draw
applied to both projects in the same iteration, so the dg-minus-da CI
reflects a genuine paired difference, not two independent CIs subtracted
afterward). Reindexes dg's and da's 4 raw predictor slots into one shared
"conceptual" order first, since which physical asymmetry each slot means
swaps between projects.
- **`plotDgVsDaComparison.m`** — original, `boxchart`-based: 3 rows per
  cortical area (dg per-subject bars, da per-subject bars, main dg-minus-da
  point with CI/asterisks); V1 only.
- **`plotDgVsDaSubjectwise.m`** — same `fitAsymmetryRegression_dgVsDa.m`
  cache, dots-and-lines convention instead of bars: row 1 = per-subject
  paired dg/da dots + grey line (observer) plus paired group means + bold
  line (group); row 2 = the dg-minus-da delta per subject, with the group
  CI/asterisks. V1 only, `orientation_minus_baseline` only.

## Directions-per-location plots
**`04_plot_betaAsymmetries/plotEachDirLocRegression.m`** (successor to the
legacy `lme2_ploteachDirLoc.m`, which instead reads `lme1_fit.m`'s own
`modeldata.mat`/`LME_bold.mat`) — 8 polar subplots, one per visual-field
location, each showing all 4 raw stimulus directions' data + model curve.
Data: the model curve/line is built directly from `fitAsymmetryRegression.m`'s
cached coefficients (`M * [grandInterceptFE; estimates/2]`); the empirical
dots are `meanBOLDpa`, gain- and precision-weighted across subjects, for the
same 4-direction × 8-location grid. For da, each point's *display* angle
(not its value) is rotated to local Cartesian angle (`plotShift = pa-90`)
since raw da stimulus codes are only Cartesian at the UVM (upper vertical
meridian) — dg needs no such rotation.

## Across-cortical-area plots
Two complementary layouts, both sourced from a per-ROI cached fit rather than
`meanBOLDpa` directly — successors to `lme1_fit.m`'s "master figure"/
"across-ROI" figure sections:
- **`04_plot_betaAsymmetries/fitAsymmetryRegression.m`** — the shared fit: a
  joint WLS regression of all 4 asymmetries simultaneously from `meanBOLDpa`,
  gain- and precision-weighted (both looked up fresh per ROI inside the fit
  loop, since gain/precision genuinely vary by cortical area), with a
  1000-draw subject-resample bootstrap for CIs. Cached per ROI under
  `derivatives/summaryTables/regressionResults/<fitLabel>/<roi>.mat`
  (`fitLabel` = `dg`/`dgMatched7`/`da`). This is also what the polar-plot
  harmonic overlay and directions-per-location model curve read.
- **`computeMeanBOLDNative.m`** (`03_process_groupBetas/`) +
  **`fitAsymmetryRegressionNative.m`** (`04_plot_betaAsymmetries/`) — an
  alternate, genuinely different analysis for the 2 asymmetries that are
  "native" to each project (dg: Horizontal-vs-Vertical, Cardinal-vs-Oblique;
  da: Radial-vs-Tangential, Polar-Cardinal-vs-Polar-Oblique — each is a
  single fixed stimulus/pair that means the same physical thing at every
  location, so no polar-angle binning is needed). `computeMeanBOLDNative.m`
  averages vertices within eccentricity **0.5-12°** (explicitly **no** varexp
  filter, no location restriction — a deliberate departure from
  `meanBOLDpa`'s ecc[4,8]+varexp≥0.1+location-binning), producing
  `meanBOLDnative.mat`. `fitAsymmetryRegressionNative.m` then does *not* fit
  a regression (native terms don't depend on location, so there's no design
  matrix) — just a gain-corrected, precision-weighted mean of each subject's
  own pro-minus-con scalar, with the same subject-bootstrap CI convention.
  Cached under
  `derivatives/summaryTables/regressionResultsNative/<fitLabel>/<roi>.mat`.
- **`04_plot_betaAsymmetries/plotROISummary.m`** — one figure per cortical
  area, asymmetries (4, or 2 under `'nativeOnly',true`) along the x-axis.
- **`04_plot_betaAsymmetries/plotAsymmetryAcrossROIs.m`** — the inverse
  layout: one figure per asymmetry (4, or 2 under `'nativeOnly',true`),
  cortical areas along the x-axis, with the axis box doubled in height
  relative to `plotROISummary.m`.
- Both read `regressionResults/` by default, or `regressionResultsNative/`
  when `nativeOnly=true`.

## Trial-triggered-average (TTA) plots
**`02_ttave/run_ttave.m`** (per-subject) → **`ttave_compute.m`** →
**`plot_ttave.m`**, and **`run_ttaveGroup.m`** (group average, reads back the
per-subject `ttaveSignal_*.mat` files already saved rather than reloading raw
data) — peristimulus-time-locked, condition-averaged (cardinal-motion,
oblique-motion, cardinal-static, oblique-static) mean % signal change vs.
time (TR, centered on stimulus onset), with cross-subject SEM shading,
plotted separately for "data" and "model". Data: raw preprocessed BOLD
(`.mgh`, converted to % signal change) for the data curve; `derivedModelFit.mat`
(GLMsingle single-trial betas convolved with the fitted HRF via
`GLMpredictresponses`) for the model curve.

## Model vs. observed time-series plots
Two parallel pathways exist for this figure type — they differ in whether the
group-level "model" trace is GLMsingle's own fit or a simplified stand-in, and
that distinction matters for what the figure can claim.

**Legacy pathway (not GLMsingle at the group level):**
**`02_ttave/run_runTimeseries.m`** (per-subject, per-run) →
`runTimeseries_compute.m` → **`plot_runTimeseries.m`**, and
**`run_groupAverageRunTimeseries.m`** (group average) →
**`plot_groupAverageRun.m`** — continuous (not trial-averaged) % signal
change vs. time (seconds) for a full run, observed overlaid with predicted,
with shaded motion/static stimulus periods and R² in the title. Data:
observed = raw BOLD → % signal change (same conversion as TTA); predicted =
`derivedModelFit.mat` for the per-subject version, or (since per-vertex
GLMsingle HRF fits can't be averaged across subjects) a freshly refit
canonical-HRF GLM on the group-mean observed trace (`refitGroupMeanTrace.m`)
for the group version — a materially different, simpler model (generic
canonical HRF, plain OLS) than GLMsingle's own fit, found and flagged
2026-09-01. Also, `run_groupAverageRunTimeseries.m`'s subject list is
`intersect(dgSubs, daSubs)` for *both* projects, so dg's group average there
only ever uses the ~7-8 subjects shared with da, never the full 13.

**Current pathway (genuinely GLMsingle-based, this is what produces the
manuscript figure):**
**`02_ttave/computeSubjectRunTraces_GLMsingle.m`** (per subject, all 8
polar-angle wedges in one pass) → **`runRawTimeseriesPermutationTest.m`**
(the aggregate stat) and **`plotRawTimeseriesWithinAcross_GLMsingle.m`** (the
figure). Both "observed" and "predicted" are genuinely GLMsingle-based:
observed = raw BOLD → % signal change using GLMsingle's own per-vertex
`meanvol` reference **first**, *then* polynomial-denoised (getting this order
backwards silently reproduces the legacy pathway's ~-100% offset bug — see
below); predicted = GLMsingle's actual fitted single-trial betas
(`results.modelmd`) convolved with each vertex's own GLMsingle-selected HRF
(`results.HRFindex`) via `GLMpredictresponses`, batched by HRF-library index
across the *whole* ROI at once (not refit, not simplified). Both reduced to
one trace per (location, run) via an exact median across all selected
vertices, then a mean across subjects (13 for dg / 7 for da, sub-0395
excluded either way; `findDuplicateDesignRuns.m` still excludes duplicate-
design runs like da/sub-wlsubj124's runs 1-2, same convention as the legacy
pathway). Session-folder discovery filters to folders that actually contain
`rawInfo.mat` (a stray empty `ses-*` folder broke naive `numel(sesNames)>1`
checks for dg/sub-wlsubj121, a subject this pipeline had never processed
before 2026-09-01).

- **`runRawTimeseriesPermutationTest.m`** — the quantitative claim: this is a
  **temporal-specificity test, not a held-out generalization test** (GLMsingle
  fits one set of betas per subject pooling all runs, not separately per run,
  so "run j's fit" is the same betas reconstructed through run j's own
  trial-onset design, not an independently-trained model). True statistic:
  concatenate every valid run's own ("i=i") observed/predicted pair across
  *all 8 polar-angle locations and all valid runs* into one long vector each,
  compute one centered R² (`1-SSres/SStot`, matching the legacy pathway's
  formula). Null: every *derangement* of the valid run indices (a permutation
  with zero fixed points) — with ≤8 runs this is exhaustive (14,833
  derangements enumerated exactly, not Monte Carlo sampled), same derangement
  applied to all 8 locations within one null draw (a run's trial-timing
  mismatch is a property of the run itself, shared across locations). Report
  as: observed R² + exact permutation p (`p < 1/(nDerangements+1)` when the
  count is exactly 0, as it was for both dg and da) — NOT as a bootstrap CI;
  nothing here is resampled, the null is the complete, exact population of
  mismatches, and the true R² itself has no separate uncertainty band unless
  a subject-resampling bootstrap is layered on top (not currently done).
- **`plotRawTimeseriesWithinAcross_GLMsingle.m`** — the illustrative figure:
  one run's data against its own ("within") fit and a mismatched ("across")
  fit overlaid, for a chosen pair of locations and one specific run pair
  (e.g. runs 2 and 5, both directions → 4 plots for 2 locations). Line/shading
  convention: observed = solid black, within fit = solid dark red, across fit
  = dashed light grey; condition-colored shaded blocks at 15% opacity behind
  the lines — grey for all 8 moving conditions, and the 4 stationary
  conditions colored by raw stimulus code (0°→cyan, 90°→yellow, 45°→dark
  green, 135°→magenta; project-invariant — dg's codes are
  horizontal/vertical/right-/left-leaning, da's are the same codes as
  annulus(tangential)/pinwheel(radial)/CW-spiral/CCW-spiral, verified against
  `dg_stimNames`/`da_stimNames` in `createTTaveTable.m`). Axes sized to a
  fixed 15.37×3.37cm (excluding the title), shared y-axis limits across every
  plot one call produces.

"Model" throughout both pathways means a GLM reconstruction (single-trial
betas × fitted/canonical HRF), not a separate mechanistic model — that's what
`Models/` (see `AGENTS.md`) is for.

### R² surface visualization
**`02_ttave/computeGLMsingleR2.m`** (precompute) + **`02_ttave/viewGLMsingleR2.m`**
(`viewGLMsingleR2(subject, project, hemi, filtered, ...)`, display) — renders
GLMsingle's own per-vertex R² (`modelOut{4}.R2`, the full
`TYPED_FITHRF_GLMDENOISE_RR` model) on that subject's inflated FreeSurfer
surface, with the V1 label outlined in red. Not a MATLAB `patch`/`trisurf`
plot — `computeGLMsingleR2.m` writes the R² as an `.mgz` overlay file first
(`fullmodel_rsquared.mgz`, every vertex, or `fullmodel_rsquared_filteredvertices.mgz`,
zeroed outside the V1+eccentricity+pRF-R²≥0.1 mask used elsewhere in this
pipeline), and `viewGLMsingleR2.m` shells out to FreeSurfer's `freeview` to
display it (`system(...)`, backgrounded) — needs `freeview` on `PATH` and
`$FREESURFER_HOME` set. Strictly per-subject, per-hemisphere, per-project;
there is no group-average surface rendering. No new file is saved by the
viewer itself — it only opens an interactive `freeview` window.

## Implementation gotcha: figure-export loops need `drawnow` before `print()`
`lme1_fit.m` previously had a real bug where two figures generated
back-to-back in the same loop (`mainSubset`/`derivedSubset` asymmetries) were
saved with *identical* content — colors, legend, and data all from the wrong
figure — because `print()` grabbed a stale render before MATLAB finished
drawing the legend/boxcharts. Fixed by calling `drawnow;` immediately before
every `print(...)` inside a loop that creates multiple figures. Any new
figure-generation loop added to this pipeline should do the same.
