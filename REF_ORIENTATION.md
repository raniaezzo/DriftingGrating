# Orientation / reference-frame logic across the `04_plot_betaAsymmetries` pipeline

This document exists because the polar-vs-Cartesian relationship between stimulus
condition labels and what they *mean* at a given visual-field location is the
single most error-prone part of this pipeline. Read this before touching any
script that builds `mainCardinal`/`derivedCardinal`/`mainSubset`/`derivedSubset`,
reads `meanBOLDpa`, or plots anything indexed by direction/orientation AND
location together.

## 1. The two experiments

- **dg** ("drifting grating"): stimuli are gratings with a genuinely fixed
  **Cartesian** orientation. A grating nominally labeled "0°" *is* a horizontal
  grating, physically, at every location in the visual field. Raw condition
  codes (0/45/90/135) are already absolute Cartesian orientation — no
  location-dependent reinterpretation is ever needed to know what a dg
  condition means visually.
- **da**: stimuli are annulus / pinwheel / spiral patterns, defined in a
  **polar** reference frame relative to the stimulus's own center. The same
  nominal code (0/45/90/135) is a *different local Cartesian orientation*
  depending on which visual-field location it's shown at. This is the source
  of essentially every bug this pipeline has had involving da.

## 2. The core geometric fact (da only)

Verified exhaustively (all 32 (stimulus, location) pairs) against ground truth
provided by the PI:

```
localCartesianAngle(S, L) = mod(S + L - 90, 180)
```

where `S` is the raw da stimulus code (0=annulus, 45=CW spiral, 90=pinwheel,
135=CCW spiral) and `L` is the location's own polar angle (0/45/.../315).
Cartesian angle convention: 0=horizontal, 45=right-leaning, 90=vertical,
135=left-leaning (all mod 180, since orientation is axial).

Consequences:
- At **L=90 (the UVM, upper vertical meridian)**, `localCartesianAngle = S`
  exactly — the polar and Cartesian frames coincide there by construction.
  This is why UVM-anchored code exists: da's raw codes are defined so that at
  the UVM they read directly as Cartesian codes.
- **Pinwheel (S=90) is always locally radial**, and **annulus (S=0) is always
  locally tangential**, at *every* location — this is a rotation-invariant
  property of the stimulus pattern itself (a radiating pattern is radial
  everywhere; a ring pattern is tangential everywhere). CW/CCW spirals
  (S=45/135) are always the two patterns exactly *between* radial and
  tangential (45° off both).
- Everywhere else (Cartesian horizontal/vertical/oblique for a *specific*
  location), the mapping genuinely rotates with `L`, and code that treats a
  raw da code as if it were already a Cartesian angle will silently be wrong
  at every location except the UVM.

dg has no analog of this rotation — raw code IS Cartesian angle, always,
`plotShift`/UVM-correction never applies to dg.

## 3. `meanBOLDpa`: what it is

4-D array: `(contrast, location, roi, subject)`. `location`'s 8 columns are
ordered `anglevals = [90, 45, 0, 315, 270, 225, 180, 135]` (both dg and da use
this same location ordering/convention — location is a physical visual-field
position, project-independent). `contrast` indexes into
`projectSettings.contrasts_dict.contrasts`, whose `.dg_contrast_name` field
(named `dg_contrast_name` for historical reasons, but shared by both
projects) gives the human-readable condition name. Full 29-entry table:

```
 1: cardMsep        8: m0_v_s90    15: m315_v_s45  22: m45_v_b     29: s135_v_b
 2: oblMsep          9: m90_v_s0   16: cardsVblank  23: m225_v_b
 3: allmValls       10: m180_v_s90 17: oblsVblank   24: m135_v_b
 4: allsVblank      11: m270_v_s0  18: m0_v_b       25: m315_v_b
 5: allmVblank      12: m45_v_m135 19: m180_v_b     26: s0_v_b
 6: cardmVblank     13: m135_v_m45 20: m90_v_b      27: s90_v_b
 7: oblmVblank      14: m225_v_m135 21: m270_v_b    28: s45_v_b
```
*(indices 1-17 are motion-related contrasts, not used by the
`orientation_minus_baseline` pathway this whole pipeline runs; listed for
completeness — verify against `loadConfig(...).contrasts_dict` if precision
matters for a motion-comparison script.)*

**The four that matter for everything in this directory** (`comparisonName =
'orientation_minus_baseline'`):

```
26: s0_v_b    (raw code 0)
27: s90_v_b   (raw code 90)
28: s45_v_b   (raw code 45)
29: s135_v_b  (raw code 135)
```

These are looked up by name, never hardcoded, via e.g.
`s0_idx = find(strcmp(contrastnames,'s0_v_b'))` — always resolves to 26 for
the current `COLORS.json`/config, but don't hardcode 26-29 in new code either;
look them up.

**Critically: these same 4 raw per-direction contrasts are the *only* raw
material available.** There is no separate set of "already-corrected"
contrasts distinct from these — indices 26-29 ARE the 4 raw individual
stimulus-type responses. Whether reading them "natively" (single contrast,
correct everywhere with no further work) or needing them "derived" (combined
across the 4 via a location-aware computation) depends entirely on which
asymmetry concept you're after and which project you're in — see §4.

## 4. Native vs. derived asymmetries (and why the pair swaps between projects)

Four asymmetry concepts exist. For a given project, two are **native**
(directly = one raw contrast, or one fixed pairing of two raw contrasts —
correct at every location with zero further computation) and two are
**derived** (require combining across the 4 raw directions in a
location-aware way to construct).

| | dg native | dg derived | da native | da derived |
|---|---|---|---|---|
| | cardinal vs oblique | radial vs tangential | radial vs tangential | cardinal vs oblique |
| | vertical vs horizontal | polar cardinal vs polar oblique | polar cardinal vs polar oblique | vertical vs horizontal |

**Why native is native**: for dg, "cardinal" (0/90) vs "oblique" (45/135) and
"horizontal" (0) vs "vertical" (90) are directly Cartesian properties of the
raw codes themselves — true at every location, no rotation. For da, "radial"
(always S=90) vs "tangential" (always S=0), and "polar cardinal" (radial-or-
tangential, i.e. S∈{0,90}) vs "polar oblique" (S∈{45,135}), are rotation-
invariant properties of the *stimulus pattern* — also true at every location,
no rotation, for the reason in §2.

**Why derived needs a location mapping**: dg's "radial vs tangential" asks
"does this grating's Cartesian orientation point at/away from THIS location"
— inherently a location-relative question, needs `abs(direction - location)`.
da's "cardinal vs oblique" (true Cartesian horizontal/vertical/oblique) asks
"what is this stimulus's local Cartesian orientation at THIS location" —
needs the §2 rotation, since S alone doesn't determine it.

## 5. Where the correction lives in the pipeline — the native path never needs it, the derived path always does

**Native path — `retrieveProConIdx.m`** (used by `plot1_experimentalCond.m`
and the non-derived branch of `plot2_experimentalCond.m`/production runners):
returns FIXED contrast indices into `meanBOLDpa`'s 26-29 block, no location
loop, no rotation:
```matlab
% orientation_minus_baseline branch:
subset==1 & da:  pro=27(radial=s90),   con=26(tangential=s0)
subset==1 & dg:  pro=26(horizontal=s0), con=27(vertical=s90)
otherwise (both): pro=26:27 (cardinal/polar-cardinal), con=28:29 (oblique/polar-oblique)
```
This is correct *by construction* for whichever pair is native to that
project — verified numerically this session (contrast 27 alone reproduces
the UVM-corrected reconstruction of da's radial value to 4 decimal places at
every one of the 8 locations).

**Derived path — `compute_derivativeDirections.m`**: the ONLY place in this
pipeline that calls `deriveLocalMotionfromUVM`, and only for da
(`if strcmp(projectName,'da')...`). For dg it uses a direct
`abs(polarAngles(polarIndex) - motionValue)` distance comparison — no UVM
call needed, because dg's raw codes are already Cartesian (§1). Produces one
3-way split per project (dg: radial/tangential/neut; da: horizontal/
vertical/neut — internally labeled provals/convals/neutvals, printed as
"is H"/"is V"/"is NEUT" for da). Both of a project's derived concepts come
out of this single split:
- `radialVsTangential` (dg) / `verticalVsHorizontal` (da) = rows 1 & 2
  directly (subset=1 callers, `n_derivedConditions{2} = {1, 2}` in
  `plot_NeuralAsymmetries.m`).
- `polarCardinalVsPolarOblique` (dg) / `cardinalVsOblique` (da) = mean(rows
  1,2) vs row 3 (subset=0 callers, `n_derivedConditions{1} = {1:2, 3}`) —
  "polar cardinal"/"cardinal" is "aligned with EITHER native axis", "polar
  oblique"/"oblique" is "aligned with neither" (the neut bucket).

**Regression fit — `fitAsymmetryRegression.m`**: builds ALL FOUR terms
(`mainCardinal`, `derivedCardinal`, `mainSubset`, `derivedSubset`)
simultaneously from the same 4 raw per-direction contrasts, via a per-
`(direction, location)` closed-form formula (not a call to
`deriveLocalMotionfromUVM`) — but the formula's own branching (the
`isCardMd`/`isOblMd`/`proH`/`conV` logic for da) already encodes the correct
location-dependent relationship; it's an algebraically-derived equivalent of
the UVM rotation, not a missing correction. **Confirmed correct**: the design
matrix is location-based and already accounts for each orientation condition
in each experiment (verified independently against the native-path ground
truth in §5's `retrieveProConIdx` block: matches to numerical precision at
every location for the terms that were checked). Do not "fix" this file to
call `deriveLocalMotionfromUVM` — it doesn't need it, and doing so would
likely double-apply the correction.

**Continuous harmonic overlays (`plot1_experimentalCond.m`'s polar plots,
`plotEachDirLocRegression.m`'s smooth curve)**: both are built directly from
`fitAsymmetryRegression.m`'s coefficients using verified closed-form
trigonometric formulas (`cos(2θ)`, `cos(4θ)`, `cos(2(θ-θᵥ))`, `cos(4(θ-θᵥ))`
and da's product-form `derivedSubset`) — since the underlying fit is already
correct (previous paragraph), these need no additional rotation for their
*values*. `plotEachDirLocRegression.m` needed a separate, purely cosmetic
**plotting-angle** rotation — see §6.

## 6. The one place a *plotting* rotation was needed: `plotEachDirLocRegression.m`

This script is unique in the pipeline: it's the only one that displays a full
sweep across all 4 raw directions *within a single location's own subplot*,
arranged spatially around a circle. Every other script only ever plots a
scalar pro/con value per location — never "where within this location's
own compass should raw direction code X land." That question only has an
answer once you fix what "0° on this subplot's circle" *means* — and the
natural, useful answer is "local Cartesian angle," per §2, not "raw code."

The values (data and model) were never wrong — `weightedBold`,
`mainCardinal`/`derivedCardinal`/`mainSubset`/`derivedSubset`, and the
regression betas are all built exactly as in §5, no rotation applied to any
of them. The fix (see the file's own comments around `plotShift`) is: when
drawing a point/curve sample for raw direction `θ_stim` in the subplot for
location `pa`, plot it at angle `mod(θ_stim + (pa - 90), 360)` for da, and at
angle `θ_stim` unchanged for dg. `pa - 90` is exactly the rotation needed to
turn "raw code" into "local Cartesian angle" per the §2 formula (since
`localCartesianAngle = θ_stim + (pa - 90) mod 180`, extended to the full
mod-360 mirrored display). Verified: after this fix, "radial" (the point at
plotting-angle = `pa`, i.e. pointing away from the figure's center) is > 
"tangential" (plotting-angle = `pa+90`) at 6/8 locations for da, with the
same two near-tied exceptions (pa=0, pa=180, sub-0.01 magnitude) that also
appear in the independent native-contrast ground truth from §5 — convergent
confirmation the rotation is right, not a coincidence.

**Rule of thumb for any new plot**: if you are ever arranging multiple raw
orientation codes around a circle (any polar/compass-style layout, any
project), the display angle should be the *local Cartesian angle*
(§2 formula for da, raw code unchanged for dg) — never the raw stimulus
code directly, except for dg or for da specifically at the UVM location.

## 7. Complete (stimulus, location) → local Cartesian angle table (da)

Cartesian labels: H=horizontal(0°), R=right-leaning(45°), V=vertical(90°),
L=left-leaning(135°). Rows = raw stimulus code, columns = location (mod 180,
since orientation is axial — location 0 and 180 give identical results, etc).

| raw code | L=0 | L=45 | L=90 (UVM) | L=135 |
|---|---|---|---|---|
| 0 (annulus)      | V | L | H | R |
| 45 (CW spiral)   | L | H | R | V |
| 90 (pinwheel)    | H | R | V | L |
| 135 (CCW spiral) | R | V | L | H |

Read a column top-to-bottom and you'll see it's the same H/R/V/L cycle,
just started at a different point — that rotation-by-90°-per-45°-of-location
*is* the `(S+L-90) mod 180` formula. dg has no equivalent table: raw code IS
the Cartesian label at every location, by definition (§1).

## 8. Script-by-script: is the transform relevant, and where

All paths relative to `AnalysisCode/04_plot_betaAsymmetries/` unless noted.

| Script | Relevant? | Where / how |
|---|---|---|
| `retrieveProConIdx.m` | Yes — native path | No transform. Returns fixed indices into `meanBOLDpa`'s 26-29 block (§5). Correct by construction; never touch to add rotation logic. |
| `deriveLocalMotionfromUVM.m` | Yes — **is** the transform | The rotation utility itself (§2). Called only from `compute_derivativeDirections.m` (da branch) in this session's verified pipeline. |
| `compute_derivativeDirections.m` | Yes — derived path | Calls `deriveLocalMotionfromUVM` explicitly, da only (§5). dg branch does a direct `abs(direction-location)` comparison, no UVM call, no rotation needed. |
| `compute_derivativeEachAbsDirection.m` | Likely yes | Also calls `deriveLocalMotionfromUVM` per a repo-wide grep this session, but its internals were **not examined** — read it fully and verify against §2/§5's pattern before modifying or trusting it. |
| `fitAsymmetryRegression.m` | Yes, but implicitly | Builds all 4 predictor terms directly from raw `(direction, location)` pairs via its own closed-form branching (§5) — mathematically equivalent to the UVM rotation, not a separate call to it. **Confirmed correct; do not add a `deriveLocalMotionfromUVM` call here** — would double-apply the correction. |
| `plot1_experimentalCond.m` | Indirect only | Consumes already-correct data — either `retrieveProConIdx`'s native indices, or a `proconMatrix` the *caller* (`plot_NeuralAsymmetries.m`/production runners) already built via `compute_derivativeDirections.m`. Its own harmonic overlay reads `fitAsymmetryRegression.m` coefficients directly; theta there is *location*, not orientation, so §2's rotation doesn't apply to it at all. |
| `plot2_experimentalCond.m` | Indirect only | Same pattern as `plot1_experimentalCond.m` — native vs. derived data arrives pre-built from the caller; this file does no transform itself. |
| `plot_NeuralAsymmetries.m` (and the `run_production_*.m` runner scripts) | Yes — orchestration | Decides whether to call `compute_derivativeDirections.m` before invoking `plot1`/`plot2` for a "derived" asymmetry call (`n_derivedConditions` grouping, §5). This is where native-vs-derived routing actually happens; get this wrong and the wrong raw material reaches the plotting scripts regardless of how correct they are internally. |
| `plotEachDirLocRegression.m` | Yes — plotting-angle only | Reads the 4 raw contrasts directly (no native/derived routing at all — always builds all 4 terms itself, mirroring `fitAsymmetryRegression.m`'s formula). Needed its own **plotting-angle-only** rotation (`plotShift = pa-90` for da, 0 for dg), since it's the only script displaying a full orientation sweep within one location (§6). Data/model *values* were never wrong. |
| `lme2_ploteachDirLoc.m` | Yes — legacy predecessor of the above | Different data source (`modeldata.mat`/`LME_bold.mat`, its own `fitlme`), but applies the same rotation idea explicitly via `deriveLocalMotionfromUVM`, called per-row for da to determine each point's plotting angle. Superseded by `plotEachDirLocRegression.m`; kept for reference/legacy figures only. |
| `plotROISummary.m` / `plotAsymmetryAcrossROIs.m` | No | Plot one fixed asymmetry value per ROI (`pro=+β`, `con=-β` from the cached fit) — no per-orientation or per-location display, so §2's rotation is not applicable. |

**If a new script needs to add here**: ask first whether it needs the *native*
value (read via `retrieveProConIdx.m`, zero transform), the *derived* value
(route through `compute_derivativeDirections.m`, transform already applied),
or a *raw per-orientation* sweep like `plotEachDirLocRegression.m` (transform
needed only for display angle, per §6's rule of thumb) — these are three
different questions with three different correct answers, and conflating them
is how this bug class keeps recurring.
