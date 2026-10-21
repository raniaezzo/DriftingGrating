# Regression / weighting reorganization — handoff summary

Context for a fresh agent picking this up. Manuscript: "Local orientation
asymmetries in V1 depend on global stimulus properties." Working directory:
`AnalysisCode/04_plot_betaAsymmetries/`.

## Terminology note (read this first)

The paper's own vocabulary uses "ROI" for a *polar-angle location* (one of 8
visual-field wedges). The codebase — inherited, predates this work, not
renamed — uses `rois`/`roiname`/`roi_idx` for **cortical area** (V1, V2, V3,
V3a, V3b, hV4, pMT, pMST) instead. Every file below, and this document,
follows the *code's* convention: "ROI" = cortical area. Where the
paper's-sense location matters (the 8-wedge polar angle), it's called
"location," never "ROI," to keep the two apart.

## What this reorg is

A parallel, independent successor to `lme1_fit.m`'s LME-based fitting and
plotting. `lme1_fit.m` and everything under `LME_results/` are untouched and
remain independently runnable — nothing here reads or writes there. The new
pipeline instead: (1) fits one joint weighted-least-squares regression per
(project, cortical area) via `fitAsymmetryRegression.m`, caching the result;
(2) every downstream plot reads that cache rather than refitting.

Subjects: dg = all 13 (or the 7-subject `'matched'` subset shared with da);
da = 7 (sub-0395 excluded, mismatched pilot stimulus); dg-vs-da comparison =
the same 7 matched subjects.

## The regression model

For a given (project, cortical area), one row per (subject, location,
orientation-direction) — 8 locations × 4 directions × nSubj rows, some
missing (NaN, dropped by `fitlm`) if a subject/location has no data:

```
y_im = B0 + B1·X1,im + B2·X2,im + B3·X3,im + B4·X4,im + C_m + e_im,  Σ_m C_m = 0
```

`i` = (location, direction) observation, `m` = subject. `X1..X4` are the 4
±1/0-coded asymmetry contrasts (mainCardinal, derivedCardinal, mainSubset,
derivedSubset — physical meaning of each swaps between dg/da, see
`fitAsymmetryRegression_dgVsDa.m`'s `dgTermForConcept`/`daTermForConcept`).
`C_m` is a **fixed**, sum-to-zero subject intercept (not random/REML) — see
"Subject intercept" below for why fixed-not-random, and why it's
mathematically guaranteed not to change B1–B4 regardless of missingness.
Fit via `fitlm(..., 'Weights', w)`, `w` = gain-corrected-then-precision-
weighted per subject-row (see "Weighting" below). `estimates` in the saved
`.mat` = `2×[B1..B4]` (pro-minus-con scale, matching `lme1_fit.m`'s own
convention). Bootstrap CIs resample *subjects* (not rows), refitting the
whole model each draw (1000 draws) — this is what makes the CI a paired,
subject-level statistic rather than a naive trial-level one (more in Q&A
below).

### Per-subject contribution decomposition

`β̂ = (X'WX)⁻¹X'Wy` always decomposes as `Σ_m contribution_m`, where
`contribution_m = (X'WX)⁻¹ Xₘ'Wₘyₘ` uses subject m's own rows but the
*shared* (whole-pooled-dataset) `(X'WX)⁻¹` — pure algebra (distributing a
fixed matrix over a sum), true unconditionally. For a subject with
**complete** data (all 8 locations), this reduces *exactly* to
`(1/nSubj) × that subject's own precision-weighted marginal pro-minus-con
difference`, because the 4 predictors are mutually orthogonal within one
complete subject's own 32 rows. For an incomplete subject it does not reduce
that cleanly — it instead corrects for confounding between the 4 asymmetries
that missingness introduces, using whatever partial data they have. Saved as
`subjectContributions` (nSubj×4); **must be rescaled by ×nSubj** before
comparing to a raw individual difference (plot2's grey lines had a real,
now-fixed bug from skipping this rescale).

**Consequence**: the arithmetic mean of plot2's grey lines (each subject's
own raw, gain-corrected-only value) equals the model's `estimates` *exactly*
when (a) precision weights are uniform (true today — still placeholder) and
(b) every subject has complete data in that cortical area. Neither condition
is guaranteed in general; both currently hold, which is why the visual match
seems automatic but isn't unconditional.

### Subject intercept — final decision

Fixed, sum-to-zero `C_m`, chosen over random/REML and over omitting it
entirely (a per-subject constant is provably orthogonal to B1–B4 in this
design, for any missingness pattern, so all three choices give numerically
identical B1–B4 — fixed was chosen for simplicity/no convergence risk, not
because it changes the slopes). Implemented as a full-dummy refit
(`fitlm(X,y,'Weights',w,'Intercept',false)`, one indicator column per subject
with data). `grandInterceptFE` = unweighted mean of per-subject intercepts;
`subjectFixedEffects` = each subject's deviation from it (`NaN` for a subject
with zero data in that area, not 0). Bootstrap is **not** re-run for this —
the orthogonality is a structural property of any row subset, so B1–B4 under
the fixed-effects spec is guaranteed identical to the reduced-model
bootstrap already computed; re-deriving it would be pure added compute for a
guaranteed-identical number.

## Weighting architecture

Two independent corrections, applied in sequence, doing different jobs:

**Gain** (rescales each subject's raw *values*, before anything else divides
or averages): `subjectScale_m = groupGain / gain_m`, `groupGain` = geometric
mean of `gain_m` across subjects (geometric, not arithmetic, because the
applied factor is itself multiplicative and geometric mean is what makes
*that factor's* mean exactly 1). **Now cortical-area-specific** (as of this
session — see "Today's work" below); previously one V1-only value per
subject applied identically everywhere.

**Precision** (changes how much each already-gain-corrected subject's value
*counts* when combining across subjects — inside the WLS fit as a per-row
weight, and in plot1/plot2's `weightedNanMean` for the descriptive group
dot/across-subject line). **Still a placeholder**: `precisionWeightsSource =
[]` everywhere → uniform weight 1 for every (subject, cortical area), i.e.
every weighted fit/mean is currently mathematically an ordinary unweighted
one. Real values exist (`observerPrecisionWeights.csv`, computed by
`01_calculate_observer_precisionWeights/computeObserverPrecisionWeights.m`,
bootstrap-over-runs reliability, ecc 4–8°/R²≥0.1 vertex criteria) but are not
yet wired in — swapping them in is the next deferred item (see Agenda).

Both weight sources use the identical lookup contract:
`retrieveObserver{Gain,Precision}Weights2(subjects, roiname, source)` → 1×nSubj
row vector. **Difference in missing-data handling**: precision's `[]` source
is a legitimate no-op placeholder; gain has no such no-op (it's not
optional — it always rescales real data), so a missing (subject, cortical
area) row returns `NaN`, not 1 and not an error. Every caller that computes
`groupGain = exp(mean(log(gainWeights)))` **must** use `'omitnan'`, or one
missing subject silently NaNs out the correction factor for every other
subject in that cortical area (this was a real bug caught and fixed this
session — see "Today's work").

## Script inventory

| Script | What it does |
|---|---|
| `fitAsymmetryRegression.m` | THE fit, per (project, cortical area). Single source every plot below reads. |
| `fitAsymmetryRegression_dgVsDa.m` | Paired dg-vs-da fit, 7 matched subjects, one shared bootstrap draw applied to both simultaneously (genuine paired comparison, not two independently-bootstrapped estimates subtracted). |
| `plotROISummary.m` | One figure per cortical area, 4 asymmetries on x-axis. Pro/con dots = model-derived, zero-centered (`±β`, not `Gintercept±β` — matches what `lme1_fit.m` actually plots once its `meanRelative=1` setting is accounted for). |
| `plotAsymmetryAcrossROIs.m` | Inverse layout: one figure per asymmetry, cortical areas on x-axis. Same conventions as above. |
| `plotEachDirLocRegression.m` | 8-location polar "compass" grid, model prediction vs. gain+precision-weighted empirical data. Offshoot of `lme2_ploteachDirLoc.m`. |
| `plotDgVsDaComparison.m` | dg-vs-da comparison, V1-only by explicit instruction. 3-row layout: dg per-subject bars / da per-subject bars / DG−DA main point, all 4 asymmetries. See "Visual conventions" for its sign-flip and fill rules. |
| `retrieveObserverPrecisionWeights.m` | `(subjects, roiname, source)` → 1×nSubj weights. `source=[]` = uniform placeholder. |
| `retrieveObserverGainWeights2.m` | Same contract, gain. No placeholder — `NaN` for a missing (subject, area) row, not a no-op. Supersedes the old V1-only `retrieveObserverGainWeights.m` (still used by untouched `lme1_fit.m`; do not delete). |
| `computeObserverGainWeightsByROI.m` (in `01_calculate_observer_gain/`) | Computes the per-(subject, cortical area) gain table (`gainSummaryByROI.mat`/`.csv`). Generalizes `dg_computeGain.m` (V1-only) to all 8 areas, same eccentricity/R² vertex criteria as precision's script. |
| `potentiallyUseful/` | Dead/superseded scripts relocated here this session: `REDUNDANTplotAbsDirs.m`, `plot0_experimentalCond.m`, `lme1_fit_perSubjectContextDiagnostic.m`. |

## Visual conventions (settled, do not re-litigate without being asked)

- `plotROISummary.m`/`plotAsymmetryAcrossROIs.m`: `boxchart` markers (not
  scatter — tried and explicitly reverted). No separate edge color
  (`BoxEdgeColor=BoxFaceColor`). Whichever of pro/con is bigger gets the
  full-strength color, the smaller gets a 50%-white blend — reusing the same
  two fixed colors regardless of which one ends up bigger.
- Full-width 68% CI (of the bootstrapped *difference*, `coeffs`), applied
  symmetrically to both pro and con dots — approximates the conventional
  ~95% "non-overlap" significance heuristic. This is deliberately NOT
  `lme1_fit.m`'s own convention (mirroring the raw undoubled coefficient's
  CI onto each side), which is a more liberal ~68%-of-the-difference
  boundary — the two differ by a factor of ~4 in effective significance
  threshold, confirmed both analytically and against real bootstrap numbers.
- Asterisks at fixed y=0.35: `**` = 95% CI of the difference excludes 0,
  `*` = only 68% does.
- `plotDgVsDaComparison.m`: Horizontal-vs-Vertical and Cardinal-vs-Oblique
  are ×-1'd for display so all 4 asymmetries' "expected" direction reads
  positive. Labels stay in ORIGINAL word order with an explicit "(x-1)"
  marker — never relabel to match the flip. Fill state (bars and main-plot
  point) tracks the ORIGINAL pre-flip pro>con sign, not the displayed sign
  — computed as `vals(si)*sf >= 0` (sf=±1 flip factor; `sf²=1` recovers the
  original sign). "Minus" not "vs" throughout. All 8 per-subject bar panels
  share one symmetric-about-0 ylim. Main plot's y:x display ratio is
  increased (Position height ×1.6) for a squarer look — axes/ticks
  untouched, pure physical resize.
- Con condition's "transparency" everywhere in this pipeline is actually a
  flat 50%-white color blend, never true alpha — vector PDF export via
  `painters` doesn't support real alpha.
- Known failure mode: reading fill-state/color-shade off a downsampled PDF
  preview is unreliable (got it wrong twice this session). Verify against
  the underlying numbers (which of pro/con is actually bigger) before
  trusting a visual read.

## dg=7 (matched) vs dg=13 (all) infrastructure

`fitAsymmetryRegression.m`'s `'dgSubjectMode'` (`'all'`|`'matched'`) saves to
separate cache folders (`regressionResults/dg/` vs.
`regressionResults/dgMatched7/`) — never overwrites the 13-subject cache.
Downstream scripts read via `projectSettings.fitLabel` (defaults to
`projectName`) or, for `plotEachDirLocRegression.m`, its own
`'dgSubjectMode'` parameter directly.

## Today's work (this session)

1. Relocated 3 dead/superseded scripts to `potentiallyUseful/` (git mv).
2. Confirmed gain was NOT cortical-area-specific (`dg_computeGain.m` only
   ever computed one V1-based value per subject) — this was the first item
   on the "come back to this later" list from a previous session.
3. Built `computeObserverGainWeightsByROI.m`: same eccentricity(4–8°)/R²(≥0.1)
   vertex criteria as precision's script, same ROI-label reading (V2/V3
   dorsal+ventral union), one `rmModelGain` call per subject per protocol
   (mov/stat) over the union of all 8 areas' vertices, then split/aggregated
   per area. Ran it: all 13 subjects, `gainSummaryByROI.mat`/`.csv` saved.
   V1 values cross-checked exactly against the old V1-only `gainSummary.mat`
   (bit-for-bit match on displayed precision) — confirms the generalization
   didn't change V1's own numbers. Two subjects have one missing area each
   (`sub-wlsubj123`→pMT, `sub-wlsubj127`→pMST — presumably insufficient
   retinotopic coverage there).
4. Built `retrieveObserverGainWeights2.m`, mirroring precision's
   `(subjects, roiname, source)` contract exactly.
5. Found and fixed a real bug during validation: a missing (subject, area)
   gain row must return `NaN`, not error — and every `groupGain =
   exp(mean(log(gainWeights)))` call needed `'omitnan'`, or the one NaN
   subject would propagate to `Inf`/`NaN` and silently corrupt every other
   subject's correction factor for that area. Fixed in all 5 call sites
   (verified: `pMT` gain vector for dg's 13 subjects has exactly 1 NaN,
   `groupGain` stays finite, `subjectScale` isolates the NaN to just that
   one subject).
6. Wired the new gain lookup through all 5 places that previously read the
   old flat, non-ROI-aware value, **reloading it inside the per-ROI loop in
   each case** (previously, both gain and — before an earlier fix —
   precision were computed once outside the loop, silently applying one
   area's weight to every other area):
   - `fitAsymmetryRegression.m`, `fitAsymmetryRegression_dgVsDa.m`: table
     loaded once, `retrieveObserverGainWeights2` called per-area inside the
     fitting loop.
   - `plot_NeuralAsymmetries.m`: new `projectSettings.gainWeightsSource`
     field (parallel to `precisionWeightsSource`) replaces the old
     precomputed `projectSettings.observerGain`.
   - `plot1_experimentalCond.m`, `plot2_experimentalCond.m`: gain now
     looked up per-area inside their existing ROI loops, same place
     precision was already being looked up.
   - `plotEachDirLocRegression.m`: single-ROI-per-call script, updated to
     use the ROI-aware lookup for its one `roiname` instead of the old flat
     value.
   - `plotROISummary.m`/`plotAsymmetryAcrossROIs.m`/`plotDgVsDaComparison.m`
     needed NO changes — they only read the cached fit's `estimates`/
     `coeffs`, never raw data, so they inherit whatever gain (and
     precision) went into producing that fit automatically.
   - `lme1_fit.m` was deliberately NOT touched (hard constraint — see
     terminology note above; it still uses the old flat
     `retrieveObserverGainWeights.m`, which is why that file was kept
     rather than deleted).
7. Validation refit — **completed and confirmed clean**:
   `fitAsymmetryRegression('dg')`, `('dg','dgSubjectMode','matched')`,
   `('da')`, and `fitAsymmetryRegression_dgVsDa()`, all with
   `overwrite=true`, across all 8 cortical areas each. All 4 runs finished
   with no errors, including the two areas with a NaN-gain subject
   (`pMT`/`pMST`). V1 estimates match the previously-validated numbers
   exactly: dg `[-0.1863, 0.0348, -0.4646, 0.0950]`, da
   `[0.0595, -0.0281, 0.1950, -0.1994]` — confirming the ROI-specificity
   change is numerically a no-op for V1 (its own gain value didn't change,
   it just became independently addressable per cortical area) while now
   producing correct, area-specific gain everywhere else. All 4 cached-fit
   directories (`regressionResults/{dg,dgMatched7,da,dgVsDa7}/`) are
   current as of this run.

## Agenda — remaining, in the order the user wants to revisit them

1. ~~Make gain cortical-area-specific~~ — done this session (validation in
   flight, see above).
2. **Swap placeholder precision weights for the real computed values**
   (`observerPrecisionWeights.csv`, already computed, not yet wired in —
   same `retrieveObserverPrecisionWeights.m` contract already accepts a
   real table, just needs `precisionWeightsSource` populated instead of
   `[]` everywhere it's set).
3. **Run the real production pipeline** — regenerate the actual saved
   figures in `/Volumes/Vision/UsersShare/Rania/Project_dg/figures/`.
   Explicitly tabled by the user until everything else is done; nothing
   there has been touched all session (verified via file timestamps more
   than once).
4. **Pooled-WLS vs. subject-wise two-stage fit discussion** — user wants to
   revisit whether the current one-stage pooled approach (never explicitly
   fitting any single subject alone) is the right call vs. fitting each
   subject's own regression first and combining afterward. Partially
   discussed this session (see "Q&A" below) but user has not yet decided;
   treat as still open.
5. Repo-wide directory reorganization.
6. Methods write-up assistance (the subject-intercept equation above is
   one piece already settled; the rest of the methods section is not yet
   drafted).

Also explicitly declined for now: lifting the V1-only (`rois(1)`)
restriction in `plot_NeuralAsymmetries.m` — user judged it redundant with
`plotROISummary.m`/`plotAsymmetryAcrossROIs.m` already covering all 8 areas.

## Q&A this session worth preserving (methodological reasoning, not just code)

**Is the pooled fit "similar to a paired t-test" if no subject's own effect
is ever explicitly estimated first?** Point estimate: for complete-data
subjects, yes exactly (see "contribution decomposition" above) — the pooled
coefficient IS the arithmetic mean of individual differences in that case.
Inference: the bootstrap resamples *subjects*, not rows, which is what makes
it paired in spirit — but it's an approximation to a paired t-test's
closed-form SE (`std(dᵢ)/√n`), not an implementation of it, and the gap
matters most for (a) small n (7–13 — bootstrap-with-replacement has real
small-sample roughness at that n) and (b) incomplete-data subjects, who have
no well-defined "own difference" for a t-test to use in the first place —
which is also precisely where the joint model's advantage lives (it uses
partial data rather than discarding those subjects). A two-stage
fit-each-subject-then-combine alternative would make the paired-t-test
analogy exact and auditable, at the cost of an explicit, possibly lossy
decision about subjects whose own partial data isn't full rank.

**How does the model "know" how much each subject contributes?** It doesn't,
in any decision-making sense — `β̂=(X'WX)⁻¹X'Wy` decomposes into per-subject
pieces purely because matrix multiplication distributes over a sum
(`X'Wy=Σₘ Xₘ'Wₘyₘ`, true for any row partition). The same shared `(X'WX)⁻¹`
is applied to every subject's own partial sum; there is no per-subject
matrix inversion.

**Does the decomposition assume subjects have equal true effects?** No —
each complete subject's own observed difference enters the combination at
full strength (only scaled by their precision weight, which tracks
measurement *reliability*, not effect-size similarity). There is no
shrinkage toward a common value, unlike a genuine random-effects/hierarchical
model (which `lme1_fit.m`'s original LME is, and this pooled-WLS approach
deliberately is not) — worth keeping in mind for item 4 above, which is
really a three-way choice (fixed pooled / hierarchical with shrinkage /
two-stage), not a two-way one.

**Precision weighting, mechanically**: applied in exactly two places, never
to a subject's own raw value in isolation. (1) Inside the joint fit — every
one of a subject's rows gets the same per-(subject, area) weight, WLS
minimizes `Σwᵢ(yᵢ-Xβ)²`. (2) In plot1/plot2's `weightedNanMean` — group dot
at each location = `Σₛwₛxₛ/Σₛwₛ`, NaN subjects dropped from both sums.
Grey lines (per-subject values) are always plotted raw/unweighted —
weighting only ever enters when subjects are combined into one number.
