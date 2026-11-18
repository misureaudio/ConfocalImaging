# Confocal Imaging — Richards–Wolf PSF & Richardson–Lucy Deconvolution

Companion notes for the two scripts in this directory (run with the workspace
`.venv`):

| file | what it does | output |
|---|---|---|
| `richards_wolf_axial.py` | on-axis high-NA axial + lateral PSF, paraxial vs full | `richards_wolf_axial.png` |
| `rl_deconvolution.py` | 3D Poisson (Richardson–Lucy) deconvolution of a simulated confocal stack | `rl_deconvolution.png` |

Both use the workspace interpreter:
`D:/Source/hermes-dir/.venv/Scripts/python.exe` (numpy 2.4.6, scipy 1.17.1,
matplotlib 3.11.2).

Parameters throughout: `λ = 500 nm`, `n0 = 1.5` (oil), `NA = 1.4`
(→ `α = arcsin(NA/n0) = 68.96°`, `cos α = 0.3590`, `k_med = 1.885e7 m⁻¹`).

---

## Part (a) — Richards–Wolf high-NA PSF vs the paraxial (Abbe) picture

### Model

On-axis field of a high-NA focused beam (direction-cosine / angular-spectrum form),
`u = cos θ ∈ [cos α, 1]`, `β = k_med z`:

- **scalar** Richards–Wolf:  `U_s(z) = -i k_med ∫ u^{1/2} e^{i β u} du`
- **vectorial** on-axis (`E_x`):  `E_x(z) = -i k_med/2 ∫ (u^{1/2}+u^{5/2}) e^{i β u} du`
  (the extra `u^{5/2}` is the `(1+cos²θ)` polarisation weighting)
- **paraxial (Abbe)** limit of the same integral, in terms of the measurable NA
  (`α → NA/n0`):  `|U_par(z)|² = (4/z²) sin²(k_med (NA/n0)² z/4)`

The integrals are evaluated on a fixed `u`-grid (4001 pts) + trapezoidal rule,
vectorised over the whole `z`-grid. **Self-check at z=0** (closed form
`∫ u^p du = (1-cosα^{p+1})/(p+1)`): quad vs closed differ by `1.28e-4` (`p=1/2`)
and `8.63e-5` (`p=5/2`) — trapezoidal error only.

### Measured FWHMs

| quantity | paraxial | full scalar RW | full vectorial RW |
|---|---|---|---|
| **axial FWHM** | 678.0 nm | 467.9 nm | 478.9 nm |
| **lateral FWHM** | 91.9 nm | 88.4 nm | — |

**Findings**

- The high-NA focus is **tighter than the paraxial (Abbe) prediction**: the
  vectorial axial FWHM is `478.9/678.0 = 0.706 ≈ 1/√2` of the paraxial value
  (the scalar is `0.690`). This is the same `1/√2` factor that appears for the
  confocal lateral resolution — a consistent narrowing, not a coincidence.
- The **lateral** PSF is a close Airy match (`88.4` vs `91.9` nm, ~4% narrower);
  the paraxial Airy pattern is a good lateral model even at NA = 1.4, whereas the
  paraxial *axial* picture over-estimates the depth extent by ~30%.
- The **vectorial (polarisation) correction** is small on-axis (`478.9` vs `467.9`
  nm, ~2%) but is what pushes the axial ratio to almost exactly `1/√2`. Off-axis it
  is large (a significant `E_z`), which is the basis of dipole-orientation
  super-resolution.

---

## Part (b) — 3D Richardson–Lucy deconvolution

### Forward model

3D **confocal** PSF from the scalar solid-immersion (Debye) integral

```
U(r,z) = ∫_0^{kmax} J0(kρ r) e^{i kz z} kρ dkρ ,   kz = sqrt(kmed² - kρ²)
h(r,z) = |U(r,z)|⁴        (illumination × detection, symmetric)
```

unit-normalised (`sum = 1`, flux-conserving), `kmax = 2π NA/λ`. Computed on an
`(r,z)` grid (201×401) and interpolated onto the 3D grid.

Grid: `64×64×32`, `dx = dy = 10 nm`, `dz = 20 nm`. The forward model and the RL
cross-correlation use explicit FFTs with the PSF at its natural centre in a
2×-padded array and a `(N//2)` crop, so the blur response is exactly centred on
each object pixel. **Self-test:** a point source blurred through the forward model
reproduces the PSF to `max|diff| = 2.2e-19` (machine precision). (A naive
`fftconvolve(mode='same')` is off by one pixel for even-sized arrays — avoid it
when the forward model and the correction must be exactly consistent.)

### Data

Object = 4 Gaussian blobs (`σ = 3 px`) at `(16,16,6) (48,16,6) (16,48,26)
(48,48,26)`; blurred with `h`, scaled to `N_photons = 200000`, corrupted by
Poisson noise. Observed `200532` photons (expected `200000`).

### Richardson–Lucy update (Poisson ML)

From the Poisson likelihood `I_k ~ Pois((h*O)_k)`, the EM / complete-data log-
likelihood reduces (using flux conservation `Σ_k h_k = 1`) to
`Σ_j [S_j ln O_j − O_j]`, giving the multiplicative update

```
O ← O · [ h ⋆ ( I / (h*O) ) ]        (h symmetric ⇒ h⋆R = h*R)
```

which is positivity-preserving and monotonically increases the likelihood.

### Stopping rule & measured convergence

Smooth init `O0 = h*I` (a noisy init `O0 = I` makes the TV dip first and the
L-curve degenerate). The **L-curve corner** (max distance from the chord in
normalised log-log residual-vs-TV) and the **min-NMSE** (the U-curve, which needs
the ground truth) **both land at iteration 10** — the model-free L-curve agrees
with the (oracle) minimum-error iteration.

| iteration | NMSE | residual ‖I−h*O‖ | TV | peak |
|---|---|---|---|---|
| 0  | 0.732 | 881.1 | 223234 | 13.5 |
| 10 | **0.398** | 554.6 | 272583 | 82.6 |
| 50 | 0.839 | 524.9 | 373464 | 236.7 |
| 159| 1.659 | 523.9 | 541391 | 425.5 |

- **NMSE: blurred = 0.722 → deconvolved = 0.398 → 1.81× improvement.**
- The NMSE falls to a minimum at it 10 then **rises** (it 50 → 0.84, it 159 →
  1.66): past the knee the algorithm fits the Poisson noise (over-fitting), which
  is exactly what early stopping / the L-curve prevents.

### Pitfalls hit (and fixed)

1. **Even-sized-array centring** — `fftconvolve(mode='same')` is off by one pixel
   here; use an explicit FFT with the PSF at its natural centre and a `(N//2)` crop.
2. **Negative blur values** — FFT round-off puts ~`1e-5`-relative negatives in the
   blur; clip to `≥ 0` before `poisson()` (which rejects `lam < 0`).
3. **Degenerate L-curve** — a noisy init makes the TV non-monotonic; use a smooth
   init `O0 = h*I`.
4. **Chord-distance bug** — the L-curve corner constant must be
   `y0(x1−x0) − x0(y1−y0)`, not `x0·y1 − y0·x1` (the latter is its negation and
   measures distance to the wrong line).
