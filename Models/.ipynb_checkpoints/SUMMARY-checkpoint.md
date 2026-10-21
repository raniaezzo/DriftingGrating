# `model_DriftingGratings.ipynb` — summary

## relevant functions / scripts:
`model_DriftingGratings.ipynb`
`helper_functions/pyrmodel.py` +
`helper_functions/utils_image.py`

## What it does, end to end

1. **Build a steerable pyramid** (`pyrmodel.SteerablePyramidSF`, wrapping
   `plenoptic.simulate.SteerablePyramidFreq`) with 4 orientations × 6 spatial
   frequencies, `tight_frame=True`, kept on CPU (MPS can't FFT complex
   tensors as of `torch==2.2.2`).
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
5. **Feed the imbalanced energy through six normalization models**
   (below), each trying to turn that overrepresentation into a *suppressed*
   net response — the model-side analogue of the fMRI finding that
   cardinal/radial responses are *smaller*, not larger, despite more
   neurons being tuned to them.
6. **Visualize/quantify** each model's output with a shared set of plotting
   helpers (filter maps, per-channel scatter/bar, polar profiles, an 8-location
   polar-angle grid) and a local `plot_quickView` closure (Cartesian-set vs.
   polar-set mean energy, and a cardinal-minus-oblique asymmetry bar) defined
   inline in the notebook, not in `utils_image.py`.

## The six normalization models (`utils_image.py`)

All six share the same divisive-normalization skeleton,
`response = energy^p_exp / (sigma + suppression_drive)`, and differ only in
how `suppression_drive` pools energy across orientation, SF, and space.
Numbering below follows the notebook's own DONE-log numbering (cell 1),
which does **not** match the order the models were introduced in
`AGENTS.md`'s prose — mapped explicitly here to avoid confusion:

| # (notebook) | Function | Pooling | Manuscript status |
|---|---|---|---|
| 1 | `div_normalization(tuned=False)` | 4D Gaussian blur (SF, ori, x, y all pooled with one isotropic σ) | fails to reproduce suppression |
| 2 | `div_normalization(tuned=False, q_exp=2)` | same, with superlinear (`q_exp=2`) suppressive drive | succeeds, Cartesian only |
| 3 | `div_normalization(tuned=True)` | per-orientation spatial surround + untblurred cross-orientation sum | fails to reproduce suppression |
| 4 | `normalization_byAnisotropy` | SF-averaged, then **spatially collapsed** (mean per channel) before taking std across orientations | succeeds, Cartesian only |
| 5 | `normalization_byAnisotropy_NOA` | same idea as #4 but **does not collapse space** first — output stays a full `(n_ori, X, Y)` map, with a single global (non-spatial) std scalar as the suppressive drive | not discussed in `AGENTS.md`'s model list; exploratory variant of #4 |
| 6 | `normalization_byStimHomogeneity` | cosine similarity between small-σ "center" and large-σ "surround" 24-D (SF×ori) energy vectors, mapped to a per-pixel exponent `q ∈ [1,2]` (1 = dissimilar/no extra suppression, 2 = identical distributions → superlinear suppression) | newest, not yet in manuscript text; exploratory |

Caveats worth flagging for anyone extending these:

- **`canonical_normalization`** (above `div_normalization` in the file) looks
  like an earlier/alternate version of model 1/3 and is **not called
  anywhere in the notebook** — it operates on a single representative SF
  channel rather than all of them. Likely dead code kept for reference.
- **`normalization_byAnisotropy_NOA`**'s `p_exp = 1,` and `sigma = 0.01,`
  (lines 1181–1182) are tuples (trailing commas), not scalars — harmless
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
- Models 4/5 (`normalization_byAnisotropy*`) take `pixpdeg` as an argument
  but never use it — std is computed directly in energy units, not scaled
  by eccentricity/pixel geometry. Fine for the current fixed-size stimuli,
  but a latent trap if these are reused with a different `img_size`.
- Output shapes differ across models: 1/2/3/6 return full `(n_sf, n_ori, X, Y)`
  maps; 4 collapses to `(n_sf, n_ori)` scalars (no spatial map at all,
  spatial dims never existed in the output); 5 returns `(n_sf, n_ori, X, Y)`
  but with every pixel sharing one global scalar suppression value. The
  downstream plotting helpers (`retrieve_dim`, `plot_energy_per_channel`,
  etc.) branch on `arr.ndim == 2` to handle model 4's collapsed case, so this
  is accounted for, but a 7th model with yet another shape would need the
  same branch added.

## GPU-accelerated Gaussian filtering

`gaussian_filter_gpu` / `_torch_gaussian_1d` are a drop-in replacement for
`scipy.ndimage.gaussian_filter`, built because models 1/3 and 6 convolve
with a ~90px-radius kernel across axes as short as 4–6 elements (SF, ori),
where scipy's per-axis cost scales with kernel *radius* rather than axis
*length*. They reconstruct scipy's default `mode='reflect'` boundary
behavior (including periodic folding when the kernel radius exceeds the
axis length) as a dense `(L,L)` matmul per axis, auto-selecting
CUDA > MPS > CPU via `get_fast_device()`. Verified against scipy to ~1e-7
relative error; cut this notebook's runtime from ~715s to ~44s on an Intel
Mac. Used inside `div_normalization` (models 1–3) and
`normalization_byStimHomogeneity` (model 6's surround term only — the
center term, `sigma_c=1`, is cheap enough to stay on plain
`scipy.ndimage.gaussian_filter`). Models 4/5 (`normalization_byAnisotropy*`)
do no spatial convolution at all, so they don't use this path.

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
