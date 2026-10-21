# `model_DriftingGratings.ipynb` — summary

## relevant functions / scripts:
`model_DriftingGratings.ipynb`
`helper_functions/pyrmodel.py` +
`helper_functions/utils_image.py`

Bugs found and fixed in this code (plus open issues still being tracked) are
logged separately in `BUG_FIXES.md`, not duplicated here.

## What it does, end to end

1. **Build a steerable pyramid** (`pyrmodel.SteerablePyramidSF`, wrapping
   `plenoptic.simulate.SteerablePyramidFreq`) with 4 orientations × 6 spatial
   frequencies, `tight_frame=True`.
2. **Generate matched Cartesian vs. polar grating stimuli**
   (`pyrmodel.InputStimuli`) at a single SF (1 cyc/deg) and 4 orientations
   (0/45/90/135°), soft-masked to an annulus. Cartesian gratings are plain
   sinusoids in rotated `(x,y)`; polar gratings are logarithmic spirals
   (`cos(w_r·log(r) + w_a·θ)`), so "orientation" for the polar set means
   pinwheel/spiral/annulus rather than a literal line orientation.
3. **Run both stimulus sets through the pyramid** and take **energy**
   (squared magnitude of the complex coefficients, via
   `get_energy_of_filters(..., to_magnitude=False)`) → a
   `(set, stimulus) -> (n_sf, n_ori, 896, 896)` dict, `stim_energy`.
4. **Impose a cardinal/radial overrepresentation** (`imposeImbalance`): a
   flat multiplier on the two cardinal orientation channels, and a
   cosine-ish radial gain (full weight at each channel's radial polar
   angle, tapering to 1× at the orthogonal angle) — simulating V1's known
   over-representation of cardinal- and radial-tuned neurons, with the
   radial weighting now sourced from human behavioral data (per the DONE
   log) rather than Fang et al.'s macaque numbers, which are left as a
   comment for reference.
5. **Feed the imbalanced energy through five normalization models**
   (below), each trying to turn that overrepresentation into a *suppressed*
   net response — the model-side analogue of the fMRI finding that
   cardinal/radial responses are *smaller*, not larger, despite more
   neurons being tuned to them.
6. **Visualize/quantify** each model's output with a shared set of plotting
   helpers (filter maps, per-channel scatter/bar, polar profiles, an 8-location
   polar-angle grid) and a local `plot_quickView` closure (Cartesian-set vs.
   polar-set mean energy, and a cardinal-minus-oblique asymmetry bar) defined
   inline in the notebook, not in `utils_image.py`.

## The five normalization models (`utils_image.py`)

All five share the same divisive-normalization skeleton,
`response = energy^p_exp / (sigma + suppression_drive)`, and differ only in
how `suppression_drive` pools energy across orientation, SF, and space.
Numbering below follows the notebook's own section numbering (Sections
7–11), which does **not** match the order the models were introduced in
`AGENTS.md`'s prose — mapped explicitly here to avoid confusion. (An earlier,
6th model — a spatially-collapsed anisotropy model, originally "Model 4" —
was deleted as redundant with the one below; see `BUG_FIXES.md` for what
else has changed in this area since this table was first written.)

| # (notebook) | Function | Pooling | Manuscript status |
|---|---|---|---|
| 1 | `div_normalization(tuned=False)` | 4D Gaussian blur (SF, ori, x, y all pooled with one isotropic σ) | fails to reproduce suppression |
| 2 | `div_normalization(tuned=False, q_exp=2)` | same, with superlinear (`q_exp=2`) suppressive drive | succeeds, Cartesian only |
| 3 | `div_normalization(tuned=True)` | per-orientation spatial surround + untblurred cross-orientation sum | fails to reproduce suppression |
| 4 | `normalization_byAnisotropy_NOA` | SF-averaged; output stays a full `(n_ori, X, Y)` map, with a single global (non-spatial) std-across-orientations scalar as the suppressive drive | not discussed in `AGENTS.md`'s model list; exploratory |
| 5 | `normalization_byStimHomogeneity` | cosine similarity between small-σ "center" and large-σ "surround" 24-D (SF×ori) energy vectors, mapped to a per-pixel exponent `q ∈ [1,2]` (1 = dissimilar/no extra suppression, 2 = identical distributions → superlinear suppression) | newest, not yet in manuscript text; exploratory |

Caveats worth flagging for anyone extending these:

- **`canonical_normalization`** (above `div_normalization` in the file) looks
  like an earlier/alternate version of model 1/3 and is **not called
  anywhere in the notebook** — it operates on a single representative SF
  channel rather than all of them. Likely dead code kept for reference.
- **`normalization_byAnisotropy_NOA`**'s `p_exp = 1,` and `sigma = 0.01,`
  (trailing commas) are tuples, not scalars — harmless
  here only because `p_exp` is used as a plain exponent via `**` (broadcasts
  fine against a 1-tuple) and `sigma` only ever appears in `+`, but worth
  fixing if this function is extended, since it's easy to introduce a real
  bug by copying this pattern elsewhere.
- **`normalization_byStimHomogeneity`** computes `P_c = C` / `P_s = S`
  (the actual L2-sum normalization, `C / C_sum`, `S / S_sum`, is commented
  out) — so despite the docstring's claim that this step "makes the
  homogeneity measure contrast-independent," the function as currently
  written is *not* contrast-normalized; cosine similarity on the raw
  (non-renormalized) vectors is still scale-invariant per-pixel, so the
  *homogeneity* value itself is unaffected, but anyone relying on the
  stated contrast-independence rationale should check this before citing it.
- **`normalization_byAnisotropy_NOA`** takes `pixpdeg` as an argument
  but never uses it — std is computed directly in energy units, not scaled
  by eccentricity/pixel geometry. Fine for the current fixed-size stimuli,
  but a latent trap if this is reused with a different `img_size`.
- Output shapes differ across models: 1/2/3/5 return full `(n_sf, n_ori, X, Y)`
  maps; 4 also returns `(n_sf, n_ori, X, Y)` but with every pixel sharing one
  global scalar suppression value (no real spatial map). The downstream
  plotting helpers (`retrieve_dim`, `plot_energy_per_channel`, etc.) branch
  on `arr.ndim == 2` for a fully-collapsed case that no longer exists after
  the Model-4-deletion above — now dead code in those helpers, not a bug.

## GPU-accelerated Gaussian filtering

`gaussian_filter_gpu` / `_torch_gaussian_1d` are a drop-in replacement for
`scipy.ndimage.gaussian_filter`, built because models 1/3 and 5 convolve
with a ~90px-radius kernel across axes as short as 4–6 elements (SF, ori),
where scipy's per-axis cost scales with kernel *radius* rather than axis
*length*. They reconstruct scipy's default `mode='reflect'` boundary
behavior (including periodic folding when the kernel radius exceeds the
axis length) as a dense `(L,L)` matmul per axis, auto-selecting
CUDA > MPS > CPU via `get_fast_device()`. Verified against scipy to ~1e-7
relative error; cut this notebook's runtime from ~715s to ~44s on an Intel
Mac. Used inside `div_normalization` (models 1–3) and
`normalization_byStimHomogeneity` (model 5's surround term only — the
center term, `sigma_c=1`, is cheap enough to stay on plain
`scipy.ndimage.gaussian_filter`). Model 4 (`normalization_byAnisotropy_NOA`)
does no spatial convolution at all, so it doesn't use this path.

## Biological / theoretical grounding of the suppression exponents

Two distinct "exponent" ideas sit in these models' denominators, and it's
worth being precise about which parts are grounded in established
physiology/modeling literature and which are free, unconstrained choices.

**The energy step itself** (`|filter response|²`, computed once for every
model in `get_energy_of_filters`) is well-grounded: Heeger's (1992) canonical
normalization model uses the same half-squaring nonlinearity on linear
filter outputs before pooling, standing in for the roughly-expansive
contrast-response of simple/complex cells.

**A fixed, asymmetric `q_exp` on top of that** (Model 2: suppression pool is
`energy^2` while the numerator stays at `energy^1`) is a real, documented
modeling move, not an arbitrary one. The classical Naka-Rushton /
contrast-gain-control formulation (Albrecht & Geisler, 1991, "Motion
selectivity and the contrast-response function of simple cells in the
visual cortex," *Vis. Neurosci.* 7:531–546; Heeger, 1992) normally uses the
*same* exponent on numerator and pooled denominator, which saturates but
never declines at high contrast. Making the suppressive exponent larger
than the excitatory one is specifically how models produce
**supersaturation** — response rising then falling at high drive. **Verified
citation for this:** Peirce (2007, "The potential importance of saturating
and supersaturating contrast response functions in visual cortex,"
*J. Vision* 7(6):13, doi:10.1167/7.6.13) directly estimates supersaturation
prevalence in macaque visual cortex (~17% of V1 neurons, ~26% of V2) using
an extended Naka-Rushton that can fit non-monotonic curves — this is solid
support for "supersaturation is a real, documented V1/V2 phenomenon."
(Correction: an earlier pass of this document cited "Albrecht, Geisler &
Crane, 2003" alongside Peirce for this point — that citation could not be
verified on a direct check and was likely a misremembered/conflated
reference; the closest real paper by overlapping authors is Albrecht,
Geisler, Frazor & Crane, 2002, "Visual cortex neurons of monkeys and cats:
temporal dynamics of the contrast response function," *J. Neurophysiol.*
88:888–913 — confirmed to exist and topically related to contrast-response
nonlinearity, but **not independently confirmed to make the specific
supersaturation claim**, so treat it as background, not as verified
support, until checked directly.) Whether *this* dataset's asymmetries
require supersaturation specifically (versus just stronger-than-linear
suppression that doesn't invert the response) is an open, testable
question, not something the current fixed `q_exp=2` choice settles on its
own.

**A context-dependent (non-fixed) exponent** — Model 5's entire premise — has
a real precedent, but a looser one than it might first appear. Coen-Cagli,
Kohn & Schwartz (2015, *Nat. Neurosci.* 18(11):1648–1655, "Flexible gating of
contextual influences in natural vision") showed macaque V1
surround-suppression strength is modulated by whether the surround is
*statistically redundant* with the center. The comparison below is checked
directly against the paper's Online Methods (not just the abstract), term
by term against what Model 5 actually computes, so the two can be compared
apples to apples rather than just at the level of "both involve homogeneity."

*Coen-Cagli et al.'s standard normalization model* (their equation 1):

$$R = \alpha \cdot \frac{E_{RF}^{\,n}}{\sigma^n + E_{RF}^{\,n} + \beta \, E_k^{\,n} + \gamma \, E_s^{\,n}}$$

Three separately-defined, separately-weighted pool terms sit in the
denominator, each from a geometrically distinct set of filters: **$E_{RF}$**
(drive at the RF's own location/orientation/SF — also the numerator, via
standard self-normalization); **$E_k$**, the "within-RF" pool (4 filter
pairs, same size/position/SF as $E_{RF}$, spanning all 4 orientations —
orientation-untuned but spatially co-located with the RF); and **$E_s$**,
the "surround" pool (8 filter pairs arranged in a ring around the RF, all
sharing *one* preferred orientation — matched per-image to whatever
orientation dominates inside the RF, so spatially displaced *and*
orientation-selective, not an all-orientations average). $n$, $\sigma$, and
$\alpha$ (exponent, semi-saturation, output gain) are shared across every
term — the symmetric Naka-Rushton form. **The flexible model changes exactly
one thing:** $\gamma E_s^n \to \gamma \cdot \text{gate} \cdot E_s^n$ — a
binary 0/1 switch (thresholded posterior probability of homogeneity) on the
surround term alone. $E_{RF}$ and $E_k$ are never touched by the gate.

*Model 5* (`normalization_byStimHomogeneity`), the same way:

```python
C = gaussian_filter(E_full, sigma=(0,0,sigma_c,sigma_c))     # used only to compute homogeneity
S = gaussian_filter_gpu(E_full, sigma=(0,0,sigma_s,sigma_s)) # used only to compute homogeneity
H_map = cosine_similarity(C, S)          # across all 24 (6 SF x 4 ori) channels
q_exp = H_map + 1
Z = gaussian_filter(E_full ** q_exp, sigma=(sigma_c,)*4)     # the ONLY pool term
norm_energy = E_full ** p_exp / (sigma + Z)
```

| | Coen-Cagli et al. (2015) | Model 5 |
|---|---|---|
| **Pool structure** | 3 separate additive terms ($E_{RF}$, $E_k$, $E_s$), independently weighted | 1 single pool term (`Z`), no additive separation |
| **What homogeneity modulates** | Whether $E_s$ (the surround term specifically) is *included at all* | The *exponent* applied to the *entire* pool before it's formed |
| **Is "center" pool ever touched?** | No — $E_{RF}$, $E_k$ are always fully present regardless of homogeneity | Yes, implicitly — there's no structurally separate "center-only" term left untouched; the adaptive exponent acts on all of `E_full` before the single pooling step |
| **Surround geometry** | 8 discrete points in a ring, same orientation as the RF's dominant content | Smooth isotropic Gaussian blur (`sigma_s`), *all* orientations/SF pooled equally |
| **Fate of the "surround" quantity itself** | $E_s$ is a literal term summed into the denominator | `S` (the blurred surround) is used only to *derive* `q_exp`, then discarded — never itself enters `Z` |
| **What the exponent $n$ does** | Fixed (fit per neuron), applied equally everywhere — entirely separate from gating | *Is* the thing that varies with homogeneity — the whole mechanism |
| **Role of $\sigma$** | Single additive constant, same everywhere | Single additive constant, same everywhere — this part *is* structurally analogous |
| **Source of the homogeneity signal** | Bayesian posterior from a fitted generative model (multivariate Gaussian scale mixture), with covariance/prior learned on an *independent* natural-image database | Raw cosine similarity computed directly from the *same* test image's own center/surround |
| **Free-parameter parity checked?** | Yes — paper states explicitly "the flexible and standard models used the same free parameters" | Not checked against the simpler models |
| **Binary vs. continuous** | Fitted model uses a binary gate, but the paper separately confirms (via synthesized intermediate-probability images) that the underlying neural effect is continuous | Continuous by construction — arguably closer to what the paper says the *biology* does than their own fitted (binary) model |

**Pool size, quantified.** Coen-Cagli et al.'s own filter geometry gives explicit
numbers for the "surround geometry" row above: the untuned pool ($E_k$) uses
filters of the **same size and position as the RF itself** (ratio ≈1×, not
spatially displaced at all), while the tuned pool ($E_s$) is **2–3× the RF's
spatial extent**, arranged in a ring outside it. An independent data point
corroborates the tuned/surround side specifically: Cavanaugh, Bair & Movshon
(2002, *J. Neurophysiol.* 88:2530–2546), using annular-grating stimuli on
macaque V1, report the suppressive surround's diameter is on average
**3.2× ± 2.0 (geometric mean ± SD) the diameter of the excitatory center**,
fairly constant across eccentricities. Worth being precise about what each
number does and doesn't cover: Cavanaugh et al.'s 3.2× is measured from the
*overall* surround (their annular-stimulus method doesn't separate tuned
from untuned contributions), so it isn't itself a tuned-vs-untuned split —
it corroborates the *tuned/surround* side of Coen-Cagli's 2–3× figure, not
their 1× untuned figure independently. Their companion paper (Cavanaugh,
Bair & Movshon 2002b, *J. Neurophysiol.* 88:2547–2556) characterizes the
surround's orientation-tuning *properties* (bandwidth, preference match to
the center) but does not report a separately-measured size for an
untuned-only component — so the clean "untuned ≈1×, tuned ≈2–3×" split used
above is Coen-Cagli's own model specification, not independently replicated
by a size measurement in the classic literature.

**The sharpest point:** Coen-Cagli's gate never touches an exponent.
Homogeneity controls *term inclusion* — a switch on one additive component.
Model 5's homogeneity controls *the exponent itself* — a continuously
varying power applied to the whole, undifferentiated pool. These are
different classes of nonlinearity: switching a linear term on/off changes
suppression's overall *magnitude*; changing an exponent changes how
*superlinearly* suppression scales with the underlying energy — closer in
spirit to the supersaturation mechanism above (Model 2's `q_exp`) than to
anything in the 2015 paper's own gating mechanism. Both being triggered by a
"homogeneity"-like signal is a conceptual parallel, not a structural one:
there is no term in Coen-Cagli's equation that plays the role Model 5's
`q_exp` plays, and no exponent in Model 5 that plays the role their fixed
$n$ plays.

**Can the denominator's exponent legitimately differ from the numerator's,
and can a homogeneity-like signal act on suppression more than excitation?**
Three separate claims, worth not conflating:

1. *Can `q` (denominator) differ from `p` (numerator) at all?* Established —
   this is the supersaturation literature already covered above (Albrecht,
   Geisler & Crane 2003; Peirce 2007). Worth noting the irony, though:
   Coen-Cagli et al.'s own equation uses the *same* `n` everywhere
   ($E_{RF}^n$, $E_k^n$, $E_s^n$ all share one exponent) — the paper being
   leaned on for homogeneity-gated suppression does not itself use an
   asymmetric exponent. The asymmetric-exponent precedent and the
   homogeneity-gating precedent come from different sources, not the same
   one, so citing both together for Model 5 is citing two separate ideas,
   not one unified result.
2. *Does homogeneity currently act only on suppression, not excitation, in
   Model 5?* Yes, already, by construction: `p_exp=1` is fixed and untouched
   by homogeneity; only `q_exp` (feeding `Z`) varies with it. So this part
   of the question is already answered by the existing code, not an open
   design choice.
3. *Is there precedent for homogeneity controlling the degree of asymmetry
   itself* (the size of the `q`-vs-`p` gap), rather than just picking a
   fixed `q`? This is the part with no clean citation found so far. Nothing
   pulled from Coen-Cagli, Festa, Weiss, or Farzmahdi ties the *magnitude*
   of a q-vs-p gap to a context/statistics signal the way Model 5 does —
   their gating is binary inclusion of an additive term, not a continuously
   scaled exponent gap. This specific move — letting a homogeneity-like
   signal set how much more superlinear suppression is than excitation —
   looks like a genuine extension beyond what's been published here, not an
   implementation of a documented mechanism. Worth stating plainly if this
   design is kept: it combines two separately-precedented ideas (asymmetric
   exponents; homogeneity-gated suppression) in a way that is itself novel,
   not independently supported by either source on its own.

Schwartz & Simoncelli (2001) is the earlier, more general ancestor of the
idea (normalization pool strength adapting to natural-scene feature
co-occurrence statistics). Three more recent papers from the same group
extend the broader program, though none targets this exact
exponent-vs-gate question: Festa, Aschner, Davila, Kohn & Coen-Cagli (2021,
*Nat. Commun.* 12:3635) apply the same natural-image-statistics framework to
trial-to-trial response *variability* rather than suppression strength;
Weiss, Bounds, Adesnik & Coen-Cagli (2023, *PLoS Comput. Biol.*
19(11):e1011667) build a *stochastic*, pairwise divisive-normalization model
(mouse V1) examining whether normalization signals are shared between
neurons; and Farzmahdi, Kohn & Coen-Cagli (2025, *Nat. Commun.* 16:6757)
relate normalization-pool size/overlap to response covariability in the same
macaque V1 preparation as the 2015 paper — relevant to this project's own
open question about how to principledly choose the pooling radius
(`sigma_c`/`sigma_s`), which is currently hand-set rather than derived or
fit.

**Bottom line:** the *qualitative* claim Model 5 leans on — suppression
strength should depend, continuously, on inferred center-surround
homogeneity — is well-supported by this literature. The *specific
mechanism* is not, on essentially every structural axis in the table above.
The honest framing is that Model 5 is *inspired by* or *in the spirit of*
Coen-Cagli et al. (2015), not that it *implements* or is *directly based on*
that model — citing it as direct justification for "the suppressive
exponent should be a cosine-similarity-driven power in `[1,2]`" would
overstate the connection.

This is also why the semi-saturation-constant sensitivity noted in
`BUG_FIXES.md` (Model 5's radial-vs-tangential sign flipping with `sigma`)
matters more than it might otherwise: with both `sigma` and the exponent's
functional form unconstrained (and, per the above, not even matching the
specific mechanism of their closest literature precedent), the model has
multiple free knobs that could each independently flip the qualitative
result. Closing that gap would mean either deriving the exponent from an
explicit generative/inference model the way Coen-Cagli et al. did (fit
independently of the test images and of `sigma`), or at minimum fitting
`sigma` and the homogeneity-to-exponent mapping jointly against some
external target before reporting the sign of an asymmetry as a finding.

This is the literature review the notebook's own TODO list already flags as
outstanding ("Literature review of suppressive drive as exponent... Check
out Kough, Poggio for suppression exp.") — the papers above are the natural
next citations to fold into that review, alongside whatever Kough & Poggio
turns out to say about the general exponent question.

## Known open questions (from the notebook's own DONE/Qs/TODO log, cell 1)

The notebook's title cell is the author's running log and is more current
than this file — re-read it before assuming anything below is still
accurate. As of the version reviewed:

- No principled way yet to pick the 4D Gaussian's σ (currently one
  symmetric `std_deg=3` reused across SF/ori/x/y axes in `div_normalization`).
- Radial overrepresentation currently increases local radial energy even in
  the *balanced* population check, which isn't expected — flagged as unresolved.
- Open question on whether to treat cardinal-overrepresentation-as-suppression
  and polar-overrepresentation-as-facilitation as two distinct mechanisms
  (more neurons vs. feedback/lateral connections).
- Anisotropy models currently use a binary large/small-Gaussian split
  (Fang et al. 2023 style) rather than a graded orientation Gaussian.
- TODO list includes: plotting difference-from-mean, varying the
  normalization pool size, switching the radial profile from a linear ramp
  to a von Mises (per Fang et al.), testing that extra orientation channels
  propagate correctly, a polar-cardinal asymmetry version of these models,
  a literature review on suppressive exponents, and testing predictions for
  contrast/adaptation/natural-scene manipulations.

## Scope not covered by this review

`Steerable_Pyramid_tutorial.ipynb` and `matlab_code/` (in this same
`Models/` folder) are a separate pyramid tutorial and unrelated MATLAB
filter-design scripts, respectively — neither is called by
`model_DriftingGratings.ipynb` and neither was reviewed for this summary.
