# Verification Note — `Confocal_Imaging_Essay_2026-11.md`

Companion audit trail for the essay. The essay itself carries no process scaffolding; this file records the verification standard, the item-by-item audit (queries run, primary sources confirmed, candidates rejected), the re-derivation of every self-measured number from the final run artifacts, and the honest residual-gaps list.

Date: 2026-11-18. Interpreter: `D:/Source/hermes-dir/.venv/Scripts/python.exe` (numpy 2.4.6, scipy 1.17.1, matplotlib 3.11.2). Scripts re-run fresh for this note; both exited 0.

---

## 1. Verification standard

- **Citations.** Every reference in the essay was checked against a primary source (publisher page, DOI resolution, or the paper's own record) for existence, author list, year, volume, and pages — not against search-result similarity. A search that returns several real-but-different papers was treated as *not* confirming the specific citation. Two would-be chimera citations were caught and corrected (see §3).
- **Self-measured numbers.** Every number in the essay's tables and prose was re-derived from the final run artifacts (`run_a.log`, `run_b.log`) produced by re-running the two scripts immediately before writing the note. No number was carried from memory or from an earlier, superseded run.
- **Math.** Key claims (the confocal PSF product, the OTF convolution / doubled support, the EM derivation of the Richardson–Lucy update, the FWHM ratios) were re-derived from scratch; sign and adjoint conventions were checked (the back-projection is the cross-correlation $h\star(\cdot)$, equal to $h*(\cdot)$ only because the Debye PSF is real and symmetric).

---

## 2. Citation audit (all confirmed against primary sources)

| # | Reference | Primary source checked | Confirmed fields |
|---|---|---|---|
| 1 | Abbe (1873) | DOI 10.1007/bf02956173 (Springer); abbeverein.de reprint | *M. Schultze's Archiv für mikroskopische Anatomie* **9**, 413–468 (1873) |
| 2 | Debye (1909) | DOI 10.1002/andp.19093351406 (Wiley); multiple citing records | "Das Verhalten von Lichtwellen in der Nähe eines Brennpunktes oder einer Brennlinie", *Ann. Phys.* **335**(14), 755–776 (1909) |
| 3 | Richards & Wolf (1959) | DOI 10.1098/rspa.1959.0200 (Royal Society) | "Electromagnetic diffraction in optical systems. II. Structure of the image field in an aplanatic system", *Proc. R. Soc. A* **253**(1274), 358–379 (1959) |
| 4 | Sheppard & Wilson (1978) | DOI 10.1080/713819784 (Taylor & Francis) | "Image formation in scanning microscopes with partially coherent source and detector", *Optica Acta* **25**(4), 315–325 (1978) |
| 5 | Gu & Sheppard (1993) | DOI 10.1117/12.2308503 (SPIE); full text | "Three-dimensional coherent transfer functions in confocal scanning imaging", *Proc. SPIE* **1983** (1993) |
| 6 | Richardson (1972) | DOI 10.1364/JOSA.62.000055 (Optica); ADS; full text | "Bayesian-based iterative method of image restoration", *J. Opt. Soc. Am.* **62**(1), 55–59 (1972) |
| 7 | Lucy (1974) | DOI 10.1086/111605 (AAS); ADS; full text | "An iterative technique for the rectification of observed distributions", *Astron. J.* **79**(6), 745–754 (1974) |
| 8 | Dempster, Laird & Rubin (1977) | JSTOR record (Vol. 39, No. 1); MIT-hosted PDF | "Maximum likelihood from incomplete data via the EM algorithm", *J. R. Statist. Soc. B* **39**(1), 1–38 (1977) |
| 9 | Hansen (1992) | DOI 10.1137/1034115 (SIAM); ACM | "Analysis of discrete ill-posed problems by means of the L-curve", *SIAM Review* **34**(4), 561–580 (1992) |
| 10 | Born & Wolf (1999) | Cambridge UP page; Wikipedia; Google Books | *Principles of Optics*, 7th ed., Cambridge UP (1999), ISBN 9780521642224 |
| 11 | Goodman (2005) | Stanford publist; Google Books; library catalogs | *Introduction to Fourier Optics*, 3rd ed., Roberts & Company (2005), 491 pp. |
| 12 | Pawley (2006) | Springer book page; Google Books | *Handbook of Biological Confocal Microscopy*, 3rd ed., Springer (2006), 988 pp. |

Notes on specific confirmations:
- **Abbe (1873):** the journal is *M. Schultze's Archiv für mikroskopische Anatomie* vol. 9, pp. 413–468 (two parts, DOIs 10.1007/bf02956173 and …74). The resolution criterion comes from part IV, "Das optische Vermögen des Mikroskops."
- **Gu & Sheppard (1993):** the full text explicitly states "the spatial frequency pass band of the confocal scanning microscope is twice as large, in both axial and transverse directions, as that in the conventional microscope" — the direct primary-source support for the essay's §3 doubled-bandwidth claim.
- **Richardson (1972) / Lucy (1974):** full texts retrieved; both confirm the iterative scheme and (for Lucy) the monotone maximum-likelihood convergence that the essay's §5 relies on.

---

## 3. Candidates rejected (would-be chimera citations)

These were considered and **not** used, because verification failed or returned a different work:

1. **"Transfer function in confocal microscopy", J. Opt. Soc. Am. 68, 782 (1978), Sheppard & Wilson.** Repeated searches for this exact (title, journal, volume, page) kept returning *different* Sheppard–Wilson papers. The paper as stated could not be surfaced from any primary source — the classic fabricated-citation signature. **Replaced by** the verified Optica Acta 1978 paper (ref. 4) for the scanning/confocal transfer function, and the verified Gu & Sheppard 1993 SPIE paper (ref. 5) for the explicit doubled-pass-band statement.
2. **Debye (1909), "Beugung, Auflösungsvermögen und Dispersion in starken Objektiven und Kondensern", Poggendorff's Ann. 352, 1–62.** This (title, journal, volume, pages) could not be confirmed; searches returned Debye's unrelated 1909 radiation-pressure paper (*Ann. Phys.* 335, 57–136) instead. **Corrected to** the verified "Das Verhalten von Lichtwellen…" (*Ann. Phys.* 335, 755–776, DOI 10.1002/andp.19093351406), which is the actual origin of the Debye angular-spectrum integral (ref. 2).
3. **Abbe (1873), "Annalen der Physik 248, 1–106".** Wrong journal/volume (a chimera: correct title + year, wrong venue). **Corrected to** *M. Schultze's Archiv für mikroskopische Anatomie* 9, 413–468 (ref. 1).
4. **Goodman (2005), "Wiley".** Wrong publisher for the 3rd edition. **Corrected to** Roberts & Company (ref. 11). (The 1st/2nd editions were McGraw-Hill; the 3rd is Roberts & Company.)

---

## 4. Self-measured numbers, re-derived from the final run artifacts

Both scripts were re-run fresh (`run_a.log`, `run_b.log`); the values below are read from those logs.

**Part (a) — `richards_wolf_axial.py` → `run_a.log`:**
- Parameters: λ = 500 nm, n₀ = 1.5, NA = 1.4 → α = 68.96°, cos α = 0.3590, k_med = 1.8850e7 m⁻¹.
- z = 0 quadrature self-check (vs closed form): scalar diff = 1.28e-04, vectorial diff = 8.63e-05. (Essay quotes 1.3e-4 and 8.6e-5.)
- Axial FWHM: paraxial (Abbe) **678.0 nm**, full scalar **467.9 nm**, full vectorial **478.9 nm**.
- Lateral FWHM: paraxial (Airy) **91.9 nm**, full scalar **88.4 nm**.
- Derived ratios (recomputed): 478.9/678.0 = 0.7063 ≈ 1/√2 = 0.7071; 678.0/478.9 = 1.4154 ≈ √2; (678.0−478.9)/678.0 = 29.4%; (91.9−88.4)/91.9 = 3.8%.

**Part (b) — `rl_deconvolution.py` → `run_b.log`:**
- PSF: sum = 1.000000 (unit flux), peak = 2.822e-04.
- Point-source forward-model test: max|diff| = **2.17e-19** (essay quotes 2.2e-19).
- Data: observed = 200 532 photons, expected = 200 000.
- Per-iteration diagnostics (it : NMSE, residual, TV, peak):
  - 0: 0.73181, 881.1, 223 234, 13.5
  - 10: 0.39784, 554.6, 272 583, 82.6
  - 20: 0.47470, 534.4, 304 014, 134.4
  - 50: 0.83946, 524.9, 373 464, 236.7
  - 100: 1.28059, 523.9, 459 377, 335.7
  - 159: 1.65946, 523.9, 541 391, 425.5
- Stopping: L-curve corner @ iter **10**; min-NMSE @ iter **10** (the two agree).
- Headline: NMSE blurred = **0.7218**, deconvolved = **0.3978**, improvement = 0.7218/0.3978 = **1.81×**.

All table cells and prose values in the essay match these log lines; no value was carried from an earlier run.

---

## 5. Math re-derivation (claims checked from scratch)

- **Confocal PSF = product of two PSFs** (§3): for a point emitter, detected signal = (illumination intensity at the emitter) × (detection efficiency at the emitter). Standard; consistent with Sheppard & Wilson 1978 and Gu & Sheppard 1993.
- **Confocal OTF = convolution of two pupil autocorrelations; support 4k⊥** (§3): F{h_ill·h_det} = F{h_ill} ⊗ F{h_det}; each F{h} is the pupil autocorrelation with support 2k⊥; the Minkowski sum of two radius-2k⊥ disks is radius 4k⊥. Cutoff 4 NA/λ = 2× the widefield 2 NA/λ. Matches Gu & Sheppard 1993's explicit "twice as large" statement.
- **√2 vs ×2 are different measures** (§3): FWHM (main-lobe width) vs OTF support (total bandwidth incl. sidelobes). Stated as distinct, not conflated.
- **EM derivation of the Richardson–Lucy update** (§5): complete-data log-likelihood Q = Σ_j[S_j ln O_j − O_j] (the −O_j term uses unit-sum flux conservation); E-step S_j = Σ_k E[S_kj|I] = O_j·(h⋆R)_j with R = I/(h*O); M-step ∂Q/∂O_j = S_j/O_j − 1 = 0 ⇒ O_j ← S_j = O_j·(h⋆R)_j. Positivity-preserving, monotone likelihood ascent. Sign/adjoint checked: the back-projection is the cross-correlation h⋆(·), equal to h*(·) only for the real-symmetric Debye PSF.
- **FWHM ratios** (§4): recomputed from the log values (0.7063, 1.4154, 29.4%, 3.8%); the "≈1/√2" and "≈√2" characterizations are within 0.1–0.2% and are labeled as near-coincidences, not identities.

---

## 6. Residual gaps (honest)

- **Pixel size (dx, dy) of the 64×64×32 RL grid** was not independently re-confirmed for this note; the essay deliberately does not quote a lateral pixel size, only the (re-confirmed) NMSE/residual/TV/peak diagnostics and the stopping iteration, so no unverified number enters the essay. (The z-spacing maps slice index 6→60 nm and 26→260 nm, i.e. dz = 10 nm, which is consistent with the plotted slices.)
- **Dempster–Laird–Rubin (1977) page range:** the JSTOR record and the MIT-hosted PDF both give **pp. 1–38**; some secondary bibliographies cite 1–37. The essay uses the value confirmed from the JSTOR-hosted PDF (1–38). Flagged here because the off-by-one is a known bibliographic ambiguity.
- **Paraxial axial reference formula:** the essay quotes only the measured paraxial axial FWHM (678.0 nm) and does not print the closed-form paraxial axial profile, to avoid committing to a convention-dependent prefactor. The shape (a sinc²-like axial profile) is standard (Born & Wolf; Goodman).
- **Richards–Wolf integrals written as ∝** (proportional): the overall prefactor is a convention-dependent normalization; only the FWHM (normalization-independent) is quoted. The z = 0 closed-form self-check validates the quadrature, not the absolute scale.
- **The two √2's are a near-coincidence, not a theorem:** the essay states this explicitly (§4, item 1). The measured ratio 0.7063 is 0.1% below 1/√2; it is presented as "≈1/√2" with the physical distinction made, not as an exact identity.

---

## 7. File set (this task, in `D:\Source\hermes-dir\ConfocalImaging`)

- `Confocal_Imaging_Essay_2026-11.md` — the essay (ends at its reference list; no process scaffolding).
- `Confocal_Imaging_Essay_2026-11_Verification_Note_2026-11.md` — this note.
- `richards_wolf_axial.py` / `.png` and `run_a.log` — Part (a) script, figure, log.
- `rl_deconvolution.py` / `.png` and `run_b.log` — Part (b) script, figure, log.
- `Confocal_Imaging_Notes_2026-11.md` — pre-existing companion notes (left untouched).
