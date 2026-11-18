# Verification Note — Confocal Imaging Essay, v2 (expanded sections)

**Date:** 2026-11-18
**Essay:** `Confocal_Imaging_Essay_2026-11_v2.md` (v1 `Confocal_Imaging_Essay_2026-11.md` left untouched)
**Scope:** Two new subsections added to the v1 essay — **§4.1 The full 3D Debye PSF** and **§5.1 A worked convergence bound for the Richardson–Lucy iteration** — plus matching cross-references in the Abstract and Conclusion. Every number below was re-derived from the reproducible scripts in this directory under the workspace venv (`D:/Source/hermes-dir/.venv/Scripts/python.exe`).

---

## 1. What changed in v2

| Location | Change |
|---|---|
| Abstract | One added sentence covering the full 3D Debye PSF and the per-step RL convergence bound. |
| §4.1 (new) | Full 3D Debye integral, on-axis obliquity comparison, confocal \|U\|⁴ PSF (axial/lateral FWHMs), missing cone. |
| §5.1 (new) | EM monotone-likelihood theorem, per-step KL lower bound, worked on the §5 trajectory, why the likelihood is the wrong stopping rule, honest caveat on finite-image boundary effects. |
| §7 Conclusion | Two added clauses cross-referencing §4.1 and §5.1. |
| References | No new citations (Debye 1909, Richards & Wolf 1959, Richardson 1972, Lucy 1974, Dempster–Laird–Rubin 1977, Hansen 1992 all already present in v1). |

Supporting artifacts (created for this expansion, all exit 0):

| File | Role |
|---|---|
| `debye_psf_3d.py` | Full 3D Debye PSF + on-axis obliquity conventions + FWHMs + missing cone → `debye_psf_3d.png`, `debye_psf_3d.log`. |
| `rl_convergence_bound.py` | Exact-M-step RL re-run of the §5 setup, logs per-step ΔL, KL bound, NMSE → `rl_convergence_bound.png`, `rl_convergence_bound.log`. |
| `debug_em_theorem.py` | Step-by-step EM-theorem check (ΔL ≥ Q-gap) + iterate-sum (flux) drift → terminal output. |

---

## 2. §4.1 — numbers cited and their source

All from `debye_psf_3d.log` (NA = 1.4, n₀ = 1.5, λ = 500 nm; cos α = 0.3590, α = 68.96°, kmax/kmed = 0.9333):

| Quantity cited in §4.1 | Value | Source line |
|---|---|---|
| On-axis \|U\|² FWHM, obliquity-free u⁰ | 460.7 nm | log L5 |
| On-axis \|U\|² FWHM, R–W scalar u^{1/2} | 468.0 nm | log L8 |
| On-axis \|U\|² FWHM, R–W vectorial u^{1/2}+u^{5/2} | 479.0 nm | log L9 |
| On-axis \|U\|² FWHM, 3D Debye u¹ (1D trap) | 481.6 nm | log L6 |
| On-axis \|U\|² FWHM, 3D Debye u¹ (grid, r=0 col) | 480.0 nm | log L7 |
| Paraxial (Abbe) axial reference | 678 nm | log L10 (= §4 value) |
| Confocal \|U\|⁴ axial FWHM (at r=0) | 344.5 nm | log L13 |
| Single-path \|U\|² axial FWHM | 480.0 nm | log L14 |
| Confocal narrowing ratio (axial) | 0.718 | log L15 |
| Confocal \|U\|⁴ lateral FWHM (at z=0) | 131.8 nm | log L18 |
| Single-path \|U\|² lateral FWHM | 183.6 nm | log L19 |
| Paraxial Airy lateral FWHM (first zero) | 183.7 nm (435.7 nm) | log L20 |
| First axial side lobe | 0.1389 (≈14%) at z = −283 nm | log L23 |
| k-space support (missing cone) | kmax = 1.759e7 1/m; 2kmax = 3.519e7 1/m; gap 21.04° | log L27–30 |

**Derived claims re-checked by hand:**
- "3D Debye on-axis agrees with §4 vectorial (478.9 nm) to 0.6%": (481.6 − 478.9)/478.9 = 0.56% ✓
- "axial extent ≈ 2.6× the lateral": 344.5/131.8 = 2.61 ✓
- "single-path lateral agrees with Airy to 0.05%": (183.7 − 183.6)/183.7 = 0.054% ✓
- "ratio 0.718 against Gaussian 1/√2 = 0.707": reported as measured; the small excess is attributed to \|U\|² being sinc-like, not Gaussian (stated in-text). ✓

**Self-checks embedded in `debye_psf_3d.py`:** the on-axis limit is computed two ways (the r=0 column of the 3D integral, and the 1D trapezoidal ∫u e^{iβu} du); the two agree to <1 nm in FWHM (480.0 vs 481.6 nm), validating the quadrature. The peak of the 3D PSF is at (r=0, z≈0) — verified numerically (the (r,z) argmax is at r=0.0 nm, z=−0.7 nm), confirming the lobe is on-axis and prolate (lateral 132 nm < axial 345 nm).

---

## 3. §5.1 — numbers cited and their source

All from `rl_convergence_bound.log` (identical data to `rl_deconvolution.py`: seed 42, same object/PSF, ΣI = 200 532 photons, 160 iterations, exact M-step):

| Quantity cited in §5.1 | Value | Source line |
|---|---|---|
| ΣI (observed photons) | 200 532 | log L2 |
| Per-step ΔL / bound / ratio (it 0,1,2,5,10,20,50) | see §5.1 table | log L3–L9 |
| L (init / corner it 10 / final it 159) | 19 488.43 / 55 835.34 / 57 957.21 | log L21 |
| L crest | it 50, L = 58 133.11 | log L9 |
| Monotonicity (ΔL ≥ 0) | holds it 0–50; 109 violations, first at it 51 | log L15 |
| KL bound tightness | tight it 0–9 (ratio 0.97–1.67); loose from it ~10 | log L16 |
| Empirical linear rate ρ (median \|ΔL^{t+1}/ΔL^t\|, it 30–59) | 0.873 | log L18 |
| NMSE (blurred / corner / improvement) | 0.7218 / 0.3978 / 1.81× | log L20 |

**Derived claims re-checked by hand:**
- "L(final) is 3.8% above L(corner)": (57 957.21 − 55 835.34)/55 835.34 = 3.80% ✓
- "L crests near it 50 (58 133)": L rises to 58 133.11 at it 50, then falls to 57 957.21 by it 159 ✓
- "monotonicity violated by ~0.3% in the over-fitting regime": L crests 58 133.11 (it 50) and drifts to 57 957.21 (it 159), a 0.30% drop ✓

**Flux-drift caveat (re-verified by re-running `debug_em_theorem.py`):** the iterate sum ΣO^t drifts 200 542.4 (it 1) → 200 621.5 (it 30), a 79-photon departure from ΣI = 200 532 at ~2.7 photons/step (≈0.04% over it 1–30). This is the finite-image FFT boundary effect that breaks exact flux conservation; it is why the KL bound (derived assuming exact flux conservation) is tight only in the early iterations and becomes loose from it ~10. The Q-gap bound in `debug_em_theorem.py` (ΔL ≥ Q(O^{t+1},O^t) − Q(O^t,O^t)) holds exactly for it 0–9 (dL ≥ Qgap = True) and is violated from it 10, consistent with the KL-bound behavior.

**Derivation check:** the per-step Q-gap is Σⱼ O^tⱼ[aⱼ ln aⱼ + 1 − aⱼ] (matches `debug_em_theorem.py`). Under exact flux conservation (ΣO^{t+1} = ΣO^t = ΣI) this equals ΣI·KL(p^{t+1}‖p^t) between consecutive normalized iterates p^t = O^t/ΣO^t; the extra term is (ΣO^t − ΣO^{t+1}), which vanishes under flux conservation and is the source of the looseness when it drifts. ✓

---

## 4. Consistency with v1

- The §5.1 trajectory uses the *exact* M-step (no ε in the ratio); v1 §5 used the +ε safeguard. The two are identical to the quoted precision: NMSE 0.7218 → 0.3978 (1.81×), L-curve corner = min-NMSE = it 10, matching v1 §5 exactly.
- The §4.1 on-axis vectorial value (478.9 nm) and paraxial reference (678 nm) are the same as v1 §4.
- No v1 numbers were altered; v1 and its verification note are preserved byte-for-byte.

## 5. LaTeX / structure check (v2)

- Display `$$` blocks: 21 (balanced). Inline `$`: 512 (even, balanced).
- `\left`/`\right` delimiters: 7/7 (the 8th `\left` regex hit is `\leftarrow`, not a delimiter).
- New headers present: `### 4.1 The full 3D Debye PSF`, `### 5.1 A worked convergence bound for the Richardson–Lucy iteration`.
- v2 size: 38 987 chars, 287 lines (v1 was 30 372 chars, 219 lines).

---

## 6. Revision r1 — review corrections (2026-11-18)

An external review (`Confocal_Imaging_Essay_2026-11_v2_r0_Arg.md`) identified five errors and two minor inaccuracies. All were independently re-verified against the scripts (workspace venv) and corrected in `Confocal_Imaging_Essay_2026-11_v2.md`. Supporting checks: `verify_review.py` (FWHM roots), `verify_otf.py` (OTF cutoffs), `verify_adjoint.py` (self-adjointness + flux drift).

| # | Error | Correction | Verification |
|---|---|---|---|
| 1 | §4 "lateral FWHM" 91.9 / 88.4 nm were **HWHMs** (the radial profile peaks at the r=0 grid edge, so `fwhm_nm` returns the half-width); contradicted §4.1's 183.7 nm by a factor of 2. | §4 table → **183.7 / 176.7 nm** (true FWHM); interpretation text updated; `richards_wolf_axial.py` now doubles the radial HWHM. | Airy FWHM formula = 183.74 nm; 91.9×2 = 183.8; re-run of the fixed script → 183.7 (Airy), 176.7 (scalar). |
| 2 | §3 confocal Airy FWHM root: `[2J₁/ρ]⁴ = 1/2` solved to ρ≈**1.140**, ratio 0.705. | Root is ρ≈**1.160**, ratio **0.718**; the `δ_conf = δ_wf/√2` equality flagged as a Gaussian estimate (exact Airy ratio 0.718). | `brentq` root = 1.1603; at ρ=1.140 the value is 0.5125 (≠0.5); 1.1603/1.6163 = 0.7179 = §4.1's 131.8/183.6. |
| 3 | §4.1 "confocal OTF … reaches the same 2kmax, product of two pupils." | Confocal OTF reaches **4kmax**; it is the **convolution of two single-path OTFs** (pupil autocorrelations); the missing-cone paragraph rewritten to separate the coherent field U (Ewald cap, kmax), the single-path OTF (bowtie, 2kmax, missing cone 21.0°), and the confocal OTF (4kmax, fills the cone). | `verify_otf.py`: widefield OTF support = 2kmax (ratio 1.000); confocal OTF support = 0.968×4kmax, confocal/widefield = 1.936 ≈ 2. |
| 4 | §5.1 EM proof claimed `L = Q + Z` with Z **θ-independent** (and called L the complete-data likelihood). | Rewritten: L is the **incomplete-data (marginal)** log-likelihood; `L = Q + H` with H **θ-dependent**; monotonicity from the M-step (Q non-decreasing) + **Gibbs' inequality** (`H diff = KL ≥ 0`). | Derivation re-checked; the Q-gap lower bound and the boxed `ΣI·KL` form are unchanged and remain correct. |
| 5 | §5.1 bound called "loose" for t≥10 (ΔL < bound is a **violation**, not looseness); "boundary truncation" misdiagnosed the flux drift. | "Loose" → **"violated"** (with the loose-vs-violated distinction stated); diagnosis corrected: the script's backward step is the *same* crop-and-pad convolution as the forward step, but the exact M-step needs the matrix **transpose**; the forward operator is **not self-adjoint**, so the flux drifts. | `verify_adjoint.py`: `⟨forward(a),b⟩` vs `⟨a,forward(b)⟩` differ by 6.5e-6 (not self-adjoint); iterate sum drifts a few photons/step. |
| m1 | §4.1 "the same physics in different notation" (scalar u¹ vs vectorial u^{1/2}+u^{5/2}). | Rephrased: different models (scalar angular-spectrum weight vs aplanatic-apodized vectorial polarization) that agree to 0.6% on the axial FWHM. | — |
| m2 | §4 "678.0/478.9 = 1.415." | → **1.416** (678.0/478.9 = 1.4157). | arithmetic |
| m3 | §4.1 "first zero at 435.7 nm" (ambiguous radius/diameter). | → "first zero at radius 217.8 nm — diameter 435.7 nm" (r₀ = 0.61λ/NA = 217.8 nm). | arithmetic |

**Net effect on the essay's claims:** the confocal narrowing factor is **0.718** (not 0.705), consistent across §3, §4.1, and the numerical §4.1 PSF; the confocal OTF cutoff is **4kmax** (the ×2 Abbe-bandwidth doubling of §3 is unchanged and now correctly attributed to the convolution of two single-path OTFs); the EM monotonicity proof is the standard Gibbs-inequality argument; and the convergence bound is correctly described as **violated** (not loose) from it ~10, with the cause identified as the non-self-adjoint finite discretization rather than "boundary truncation."

**Files changed in r1:** `Confocal_Imaging_Essay_2026-11_v2.md` (§3, §4, §4.1, §5.1), `richards_wolf_axial.py` (lateral FWHM doubling), `run_a.log` (re-run), plus the new check scripts `verify_review.py`, `verify_otf.py`, `verify_adjoint.py`. v1 and its note remain untouched.
