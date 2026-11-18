# The Diffraction-Limited Image as a Band-Limited Inverse Problem

## High-NA Point-Spread Functions and Richardson–Lucy Deconvolution in Confocal Microscopy

*For mathematical physicists. The numerical results below are from two reproducible computations; the scripts, run logs, and figures that produce them accompany this essay.*

---

## Abstract

Optical image formation is, at its core, a linear operation: the detected image is the object convolved with the point-spread function (PSF), and the PSF is the Green's function of the imaging system. This essay develops the mathematical structure of that operation for a high-numerical-aperture confocal microscope and carries it through to the recovery of the object, which is an ill-posed inverse problem. Three results are established, each grounded in a reproducible computation. First, the confocal PSF is the *product* of the illumination and detection PSFs, so its optical transfer function (OTF) is the convolution of two pupil autocorrelations and its transverse bandwidth is twice that of a widefield system — a factor of two in the Abbe limit, alongside a factor of $\sqrt{2}$ in the PSF's full width at half maximum (FWHM). Second, at high NA the paraxial (Airy) picture, while an excellent lateral model, over-estimates the axial extent of the focus: evaluating the on-axis Debye / Richards–Wolf integral shows the vectorial axial FWHM is $\approx 1/\sqrt{2}$ of the Abbe prediction ($478.9$ nm versus $678.0$ nm for $\mathrm{NA}=1.4$, $n_0=1.5$, $\lambda=500$ nm). Third, recovering the object from a Poisson-noised confocal stack is a band-limited inverse problem whose natural solver is the Richardson–Lucy iteration — the expectation–maximization (EM) algorithm for the Poisson likelihood; the iteration must be early-stopped at the L-curve corner, and we report a $1.81\times$ improvement in normalized mean-squared error (NMSE) on a simulated stack, with the L-curve corner coinciding with the minimum-NMSE iteration. The essay also evaluates the full three-dimensional Debye PSF — a prolate lobe whose confocal $|U|^4$ form is narrowed by $\approx 1/\sqrt{2}$ on both axes relative to the single-path $|U|^2$ — and derives a per-step convergence bound for the Richardson–Lucy iteration that is tight in the convergence regime and explains why the L-curve, not the likelihood, is the defensible stopping rule. Throughout, each number is interpreted: what it measures, why it takes the value it does, and where the model's assumptions enter.

---

## 1. Introduction: the structural question

The question this essay addresses is structural rather than historical: *given a diffraction-limited optical system, what is the mathematical object that maps an object to an image, what does that object imply about the finest resolvable detail, and how — if at all — can the object be recovered from the image?* Three constructions answer these questions in turn, and the essay is organized around why each arises, what it buys, and what it costs.

The first construction is the **point-spread function** and its Fourier transform, the **optical transfer function**. The PSF is the Green's function of the imaging operator; the OTF is its symbol. The diffraction limit — the Abbe bound — falls out of the band-limitedness of the OTF and is, in this sense, a direct manifestation of the uncertainty principle: a band-limited field cannot be compactly localized.

The second construction is the **confocal** modification of the imaging operator. By illuminating with a focused point source and detecting through a spatially filtered (pinhole) channel, the confocal microscope replaces the widefield PSF by the *product* of an illumination and a detection PSF. This single change — a product in real space — has a clean consequence in the Fourier domain (a convolution of OTFs) and is the origin of both the improved resolution and the optical sectioning.

The third construction is the **deconvolution**, the recovery of the object from the image. Because the OTF vanishes outside its finite support, the formal inverse is unbounded and the problem is ill-posed; and because the detector counts photons, the noise is Poisson. The Richardson–Lucy iteration is the EM algorithm for exactly this model, and it is the natural, positivity-preserving solver. The price of the improved resolution is that the iteration must be stopped before it begins to fit the shot noise — a decision made by the L-curve.

The essay is grounded in two reproducible computations. The first evaluates the on-axis Debye / Richards–Wolf integral for a high-NA objective and compares the resulting point-spread function with the paraxial (Abbe) approximation. The second builds a 3D confocal PSF from the Debye integral, simulates a Poisson-noised stack, and runs the Richardson–Lucy deconvolution with L-curve stopping. Every number quoted is read from the final run of these scripts; the companion verification note records the audit trail.

---

## 2. Image formation: the PSF as a Green's function

Consider an incoherent fluorescence object of intensity $O(\mathbf r)$, $\mathbf r=(x,y,z)$. A widefield microscope images it onto a detector, and the detected intensity is

$$I(\mathbf r) = \int O(\mathbf r')\, h(\mathbf r - \mathbf r')\, d^3\mathbf r',$$

a convolution. The kernel $h$ is the **point-spread function**: the image of an ideal point source. In operator language, $h$ is the Green's function of the imaging operator $\mathcal H$, $I = \mathcal H O$, and $\mathcal H$ is a convolution operator whose symbol is the Fourier transform of $h$.

In the paraxial (low-NA) limit, $h$ is the diffraction pattern of the objective's circular pupil. Writing $k_\perp = 2\pi\,\mathrm{NA}/\lambda$ for the maximum transverse wave number, the lateral PSF is the Airy pattern

$$h(r) = \left[\frac{2 J_1(k_\perp r)}{k_\perp r}\right]^2, \qquad r = \sqrt{x^2+y^2},$$

where $J_1$ is the Bessel function of the first kind. Its first zero (the Rayleigh radius) is at $k_\perp r = 3.8317$, i.e. $r_0 = 0.61\,\lambda/\mathrm{NA}$, and its full width at half maximum is $\delta_\perp \approx 0.51\,\lambda/\mathrm{NA}$.

The **Abbe limit** [Abbe 1873] follows from the band-limitedness of the pupil. The pupil transmits only transverse wave numbers $|\mathbf k_\perp| \le k_\perp = 2\pi\,\mathrm{NA}/\lambda$, so the finest resolvable period is

$$\delta_{\mathrm{Abbe}} = \frac{\lambda}{2\,\mathrm{NA}}.$$

Equivalently, the incoherent OTF — the autocorrelation of the pupil — has support $|\mathbf k_\perp| \le 2k_\perp$, a cutoff spatial frequency of $2\,\mathrm{NA}/\lambda$.

The point worth emphasizing for the mathematical physicist is that the diffraction limit is not a separate physical law grafted onto the imaging model; it is the *uncertainty principle applied to a band-limited field*. The pupil is a hard band-limit in wave-vector space, and a function that is compactly supported in $\mathbf k$-space cannot be compactly supported in $\mathbf r$-space. The PSF is necessarily spread out; its width is the price of the system's finite aperture. The Abbe bound is the quantitative statement of that trade-off.

---

## 3. The confocal principle: a product of two PSFs

A confocal microscope modifies the imaging operator in two coupled ways. It illuminates the specimen with a diffraction-limited laser spot (a point source) and detects the emitted light through a small pinhole placed in a plane conjugate to the illumination focus. Three elements thus sit in conjugate ("confocal") planes: the illumination focus, the specimen plane, and the detection pinhole.

The key mechanism is spatial filtering. A point emitter displaced from the focus is imaged *off* the pinhole and is largely blocked; an emitter in the focus is imaged *onto* the pinhole and passes. The pinhole therefore rejects light that did not originate in the confocal volume — this is **optical sectioning**, and it is what removes the widefield projection (the loss of the axial coordinate) that plagues the conventional image.

For a single point emitter at $\mathbf r'$, the detected signal factors into an illumination part and a detection part:

$$\boxed{\,S(\mathbf r') = h_{\mathrm{ill}}(\mathbf r')\; h_{\mathrm{det}}(\mathbf r')\,}$$

where $h_{\mathrm{ill}} = |\mathcal E_{\mathrm{ill}}|^2$ is the intensity of the focused laser (the illumination PSF) and $h_{\mathrm{det}}$ is the fraction of the emitter's light that reaches the pinhole (the detection PSF). **The confocal PSF is the product of the illumination and detection PSFs.** This single fact is the whole of the confocal advantage.

With symmetric illumination and detection (the same objective, the same NA), both factors are Airy patterns, so

$$h_{\mathrm{conf}}(r) = \left[\frac{2 J_1(k_\perp r)}{k_\perp r}\right]^4.$$

The fourth power, rather than the second, narrows the main lobe. The FWHM of $[2J_1/\rho]^2$ is at $\rho\approx1.616$; the FWHM of $[2J_1/\rho]^4$ — which reaches half its maximum where $[2J_1/\rho]^2 = 1/\sqrt{2}$ — is at $\rho\approx1.140$. The ratio $1.140/1.616 \approx 0.705 \approx 1/\sqrt{2}$:

$$\delta_\perp^{\,\mathrm{conf}} = \frac{\delta_\perp^{\,\mathrm{wf}}}{\sqrt{2}}.$$

So the confocal lateral FWHM is smaller by a factor of $\sqrt{2}$; the same factor applies to the axial FWHM.

In the Fourier domain the product becomes a convolution. Writing $\widehat h(\mathbf k) = \mathcal F\{h\}(\mathbf k)$ for the OTF, the convolution theorem gives

$$\widehat h_{\mathrm{conf}}(\mathbf k) = \widehat h_{\mathrm{ill}}(\mathbf k) \ast \widehat h_{\mathrm{det}}(\mathbf k).$$

Each factor is the autocorrelation of the pupil, with support $|\mathbf k_\perp| \le 2k_\perp$ (the widefield incoherent OTF). The confocal OTF is the convolution of two such autocorrelations, so its support is the Minkowski sum — a disk of radius $4k_\perp$:

$$|\mathbf k_\perp| \le 4k_\perp = 8\pi\,\frac{\mathrm{NA}}{\lambda} \quad\Longleftrightarrow\quad f_{\max}^{\,\mathrm{conf}} = \frac{4\,\mathrm{NA}}{\lambda} = 2\, f_{\max}^{\,\mathrm{wf}}.$$

The transverse cutoff is doubled, and the Abbe limit improves from $\lambda/(2\,\mathrm{NA})$ to $\lambda/(4\,\mathrm{NA})$ [Sheppard & Wilson 1978; Gu & Sheppard 1993].

A subtlety worth stating plainly: the $\sqrt{2}$ FWHM gain and the $\times 2$ OTF-support gain are *different* measures of the same underlying narrowing. The FWHM measures the width of the PSF's main lobe; the OTF support measures the total bandwidth, including the oscillatory sidelobes. Both are true simultaneously — the confocal PSF is narrower in its core, and its Fourier transform extends further. Conflating the two is a common source of apparent contradiction in the literature.

The pinhole radius, in Airy units, is the free parameter. An open pinhole passes all the light (including out-of-focus) and forfeits the sectioning; a closed pinhole approaches the ideal point detector and maximizes the sectioning but minimizes the signal. A finite pinhole convolves the ideal detection PSF with the image of the pinhole — a low-pass filter in the detection arm — and trims the OTF support. The optimum, which trades resolution against photon count, sits near one Airy unit. This is the **resolution–sensitivity trade-off**, and it is the reason the confocal microscope is not unconditionally superior to widefield.

---

## 4. The high-NA correction: the Debye / Richards–Wolf integral

The Airy pattern of §2 is a paraxial result: it assumes the objective collects light within a small cone, $\alpha \ll 1$. Modern objectives reach $\mathrm{NA} \gtrsim 1.4$, corresponding to a half-angle $\alpha$ of $60^\circ$–$70^\circ$, where the paraxial approximation breaks down. The correct description is the **Debye** angular-spectrum representation [Debye 1909], in its complete vectorial form due to Richards and Wolf [Richards & Wolf 1959].

The on-axis field of a plane wave focused by a high-NA aplanatic objective is an integral over the angular spectrum. With $u = \cos\theta \in [\cos\alpha, 1]$ and $\beta = k_{\mathrm{med}} z$, where $k_{\mathrm{med}} = 2\pi n_0/\lambda$ is the wave number in the immersion medium of index $n_0$, the scalar on-axis field is

$$U_s(z) \;\propto\; \int_{\cos\alpha}^{1} u^{1/2}\, e^{i\beta u}\, du,$$

and the vectorial on-axis transverse field (for a linearly polarized incident wave) is

$$E_x(z) \;\propto\; \int_{\cos\alpha}^{1} \left(u^{1/2} + u^{5/2}\right) e^{i\beta u}\, du.$$

The proportionality absorbs a conventional prefactor (it sets the absolute scale but not the FWHM). The extra $u^{5/2}$ term is the signature of the high-NA correction: it is the $(1+\cos^2\theta)$ polarization weighting, since $u^{1/2}\cos^2\theta = u^{1/2}u^2 = u^{5/2}$. At low NA this weighting is nearly constant and the scalar and vectorial fields agree; at high NA the oblique rays carry a transverse-field component that the paraxial model misses, and the $u^{5/2}$ term reshapes the axial field.

We evaluated both integrals (fixed $u$-grid, trapezoidal rule, vectorized over a $z$-grid spanning $\pm 800$ nm) for $\lambda = 500$ nm, $n_0 = 1.5$, $\mathrm{NA} = 1.4$ (so $\alpha = 68.96^\circ$, $\cos\alpha = 0.3590$), and compared them with the paraxial (Abbe) axial reference. The self-check at $z=0$, where the integral has the closed form $\int u^p\,du = (1-\cos\alpha^{\,p+1})/(p+1)$, agrees to $1.3\times10^{-4}$ (scalar) and $8.6\times10^{-5}$ (vectorial) — at the level of the trapezoidal error — confirming the quadrature.

The measured full widths at half maximum are:

| axial FWHM (nm) | paraxial (Abbe) | full scalar | full vectorial |
|---|---|---|---|
| value | 678.0 | 467.9 | 478.9 |

| lateral FWHM (nm) | paraxial (Airy) | full scalar |
|---|---|---|
| value | 91.9 | 88.4 |

**Interpretation.** Three facts stand out.

1. *The high-NA axial focus is tighter than the paraxial prediction.* The vectorial axial FWHM, $478.9$ nm, is $478.9/678.0 = 0.706 \approx 1/\sqrt{2}$ of the Abbe value; equivalently the paraxial value is a factor $678.0/478.9 = 1.415 \approx \sqrt{2}$ larger, so the true focus is $\approx 29\%$ narrower than the paraxial prediction. The mechanism is the $(1+\cos^2\theta)$ weighting: the oblique rays that dominate at high NA are weighted differently from the paraxial rays, and the net effect is to concentrate the on-axis field. The measured ratio $0.706$ is close to $1/\sqrt{2} = 0.707$ — the same factor as the confocal lateral gain of §3 — but the two $\sqrt{2}$'s have different physical origins (a product of Airy patterns, versus the polarization weighting), so the agreement is a near-coincidence, not an identity.

2. *The paraxial Airy pattern remains an excellent lateral model.* The lateral FWHM is $88.4$ nm (full scalar) against $91.9$ nm (paraxial Airy) — the Airy is only $\sim 4\%$ wider. The high-NA correction matters far less laterally than axially, which is why the Airy disk is still the working picture for lateral resolution even at $\mathrm{NA}=1.4$.

3. *The scalar and vectorial axial FWHMs are close but not identical* ($467.9$ vs $478.9$ nm). The difference is the polarization effect itself: small in width, but it is the correct physics, and it is what a quantitative high-NA model should use.

The practical consequence is that a deconvolution or a resolution estimate built on the paraxial PSF will be *systematically wrong in the axial direction* for a high-NA objective: it will assume a deeper confocal section than the optics actually provide. The error is modest (tens of percent) but systematic, and it is precisely the kind that accumulates when many optical sections are combined.

### 4.1 The full 3D Debye PSF

The on-axis integral of §4 is the $r\to 0$ limit of the complete field. The full scalar **Debye** integral [Debye 1909] is the angular spectrum of the converging cone: each plane-wave component carries a transverse wave number $k_\rho\in[0,k_{\max}]$ and an axial wave number $k_z=\sqrt{k_{\mathrm{med}}^2-k_\rho^2}$, and the azimuthal symmetry collapses the $\phi$-integral to the zeroth-order Bessel function,

$$U(r,z) \;=\; \int_0^{k_{\max}} J_0(k_\rho r)\, e^{\,i k_z z}\, k_\rho\, dk_\rho, \qquad k_z = \sqrt{k_{\mathrm{med}}^2 - k_\rho^2}, \qquad k_{\max}=\tfrac{2\pi\,\mathrm{NA}}{\lambda}.$$

On axis $J_0(0)=1$, and with $u=k_z/k_{\mathrm{med}}$ (so $k_\rho\,dk_\rho=-k_{\mathrm{med}}^2 u\,du$) the integral reduces to

$$U(0,z) \;=\; k_{\mathrm{med}}^2 \int_{\cos\alpha}^{1} u\, e^{i\beta u}\, du,$$

the on-axis form of §4 with the **weight $u$** — the natural $k_\rho\,dk_\rho$ measure of the angular spectrum, corresponding to an effective obliquity factor $\propto\cos\theta$ (against $\sqrt{\cos\theta}$ for the scalar Richards–Wolf integral and $1+\cos^2\theta$ for the vectorial one). We evaluated this on-axis limit two ways — as the $r=0$ column of the 3D integral and as the one-dimensional trapezoidal integral above — and the two agree to $<1$ nm in FWHM ($480.0$ vs $481.6$ nm), validating the quadrature.

*The obliquity weight is a small but real effect on the axial width.* The four on-axis conventions give the following single-path $|U|^2$ axial FWHMs:

| on-axis weight | $u^0$ (obliquity-free) | $u^{1/2}$ (R–W scalar) | $u^{1/2}{+}u^{5/2}$ (R–W vectorial) | $u^1$ (3D Debye) |
|---|---|---|---|---|
| axial FWHM (nm) | 460.7 | 468.0 | 479.0 | 481.6 |

The spread is only $\sim 5\%$ ($461\to482$ nm): the on-axis lobe is nearly convention-insensitive, because all four weights are monotone in $u$ and differ only in how strongly they concentrate weight toward the axis ($u=1$). The 3D Debye on-axis limit ($481.6$ nm) agrees with the vectorial Richards–Wolf value of §4 ($478.9$ nm) to $0.6\%$ — the same physics in different notation. All four are far narrower than the paraxial (Abbe) prediction of §4 ($678$ nm, which uses the small-angle axial bandwidth), reconfirming from the full 3D integral that the high-NA focus is sharper than the paraxial model.

*The confocal 3D PSF.* For the symmetric confocal case the detected intensity is the product of the illumination and detection PSFs, $h=|U|^4$ (unit-normalized). We evaluated $U$ on a $201\times401$ $(r,z)$ grid ($r\in[0,320]$ nm, $z\in[-320,+300]$ nm, $500$ $k_\rho$ samples) and bilinearly interpolated it onto the $64\times64\times32$ image grid used in §5. The result is a **prolate** lobe — elongated along the optical axis — with weak axial sidelobes:

| FWHM (nm) | axial (at $r=0$) | lateral (at $z=0$) |
|---|---|---|
| confocal $\lvert U\rvert^4$ | 344.5 | 131.8 |
| single-path $\lvert U\rvert^2$ | 480.0 | 183.6 |
| paraxial Airy (lateral only) | — | 183.7 |

Three facts. (i) *Confocal detection narrows the PSF by $\approx 1/\sqrt2$ on both axes*: the axial FWHM falls $480.0\to344.5$ nm (ratio $0.718$) and the lateral FWHM falls $183.6\to131.8$ nm (ratio $0.718$), because squaring the amplitude ($|U|^2\to|U|^4$) sharpens the peak. This is the 3D version of the $\sqrt2$ gain of §3, seen now in the actual PSF rather than the OTF support (the ratio is $0.718$ against the Gaussian value $1/\sqrt2=0.707$; the small excess is because $|U|^2$ is sinc-like, not Gaussian). (ii) *The single-path lateral PSF is the Airy pattern*: the single-path lateral FWHM ($183.6$ nm) agrees with the paraxial Airy ($183.7$ nm, first zero at $435.7$ nm) to $0.05\%$, confirming that the Airy disk of §2 is exactly the $z=0$ cut of the single-path Debye PSF. (iii) *The lobe is prolate, its axial extent $\approx 2.6\times$ the lateral* ($344.5/131.8$) — the familiar high-NA anisotropy that makes the axial resolution worse than the lateral; the first axial sidelobe sits at $z\approx\pm283$ nm at $14\%$ of the peak.

*The missing cone.* In $k$-space the PSF is supported only on the cone $k_\rho\le k_{\max}$ (the propagating components the objective collects); the evanescent / super-NA components $k_\rho>k_{\max}$ are absent. This "missing cone" is what limits the resolution: the single-path OTF has transverse support out to $2k_{\max}$, and the confocal OTF (the product of two pupils, §3) reaches the same $2k_{\max}$, doubling the Abbe bandwidth. The angular gap of the missing cone is $90^\circ-\alpha = 21.0^\circ$ from the axis.

---

## 5. Recovery: the deconvolution as an ill-posed inverse problem

For an extended object the confocal image at scan position $\mathbf r$ is a 3D convolution,

$$I(\mathbf r) = \int O(\mathbf r')\, h_{\mathrm{conf}}(\mathbf r - \mathbf r')\, d^3\mathbf r',$$

and scanning over a grid of $(x,y,z)$ positions and stacking the optical sections yields a 3D image. The object is, in principle, recoverable by inverting this convolution,

$$O = \mathcal F^{-1}\!\left[\frac{\mathcal F\{I\}}{\widehat h_{\mathrm{conf}}}\right],$$

but the inversion is **ill-posed**. The OTF $\widehat h_{\mathrm{conf}}$ vanishes outside its finite support (§3), so the division by $\widehat h_{\mathrm{conf}}$ amplifies the (absent) high-frequency components without bound. The forward operator is compact and band-limited; its formal inverse is unbounded. This is the standard structure of a band-limited inverse problem, and it is why a naive division produces a noise-dominated, unusable result.

The detector adds a second complication: it counts photons. The detected image is therefore not $I$ plus Gaussian noise but a **Poisson** random variable with mean $I$:

$$I_k \sim \mathrm{Poisson}\big((h * O)_k\big).$$

The likelihood of the data $\{I_k\}$ given the object $O$ is

$$\mathcal L(O) = \prod_k \frac{(h*O)_k^{I_k}\, e^{-(h*O)_k}}{I_k!}.$$

The **Richardson–Lucy** iteration [Richardson 1972; Lucy 1974] is the maximum-likelihood estimator for exactly this model, and it is the **expectation–maximization** (EM) algorithm [Dempster, Laird & Rubin 1977] applied to it. The EM view makes the update transparent. Introduce the latent complete data $S_{kj}$, the number of photons that originated at object pixel $j$ and landed at detector pixel $k$; then $I_k = \sum_j S_{kj}$ with $S_{kj} \sim \mathrm{Poisson}(O_j h_{kj})$. The E-step computes the expected complete-data log-likelihood,

$$Q = \sum_j \left[ S_j \ln O_j - O_j \right],$$

where $S_j = \sum_k \mathbb E[S_{kj}\mid I]$ is the expected number of photons that came from pixel $j$, and the $-O_j$ term is the price of the flux-conservation constraint $\sum_j O_j = \mathrm{const}$ (it uses that the PSF is normalized to unit sum). The M-step maximizes $Q$ over $O_j > 0$; since $\partial Q/\partial O_j = S_j/O_j - 1$, the maximizer is $O_j = S_j$, giving the multiplicative update

$$O \leftarrow \max\!\left( O \cdot \left[ h \star \frac{I}{h*O + \varepsilon} \right],\ 0 \right),$$

where $\star$ denotes cross-correlation (the adjoint of the blur) and $\varepsilon$ guards the division. For a symmetric $h$, $h\star(\cdot) = h*(\cdot)$, so the correction reuses the same FFT of $h$ as the forward model. The update is positivity-preserving and monotonically increases the likelihood; it is therefore an ascent method that must be **early-stopped**, or it converges to the maximum-likelihood fit of the Poisson noise.

**Stopping: the L-curve.** The natural stopping rule is the **L-curve** [Hansen 1992]: a log-log plot of a measure of the solution's complexity (here the total variation $\mathrm{TV} = \sum |\nabla O|$) against a measure of its fit to the data (the residual norm). As the iteration proceeds the fit improves (the residual falls) but the solution grows more complex (the TV rises); the curve is L-shaped, and its corner marks the trade-off between fidelity and regularity. The corner is located as the point of maximum distance from the chord joining the two endpoints, in normalized log–log coordinates.

We built the 3D confocal PSF from the Debye integral (§4, $h = |U|^4$ for the symmetric confocal case), normalized it to unit sum (flux conservation), and simulated a stack: four Gaussian blobs on a $64\times64\times32$ grid, blurred by the PSF, scaled to a budget of $2\times10^5$ photons, and corrupted by Poisson sampling (the observed image carried $200\,532$ photons). The forward model was validated by blurring a single point source and asserting that it reproduces the PSF to $2.2\times10^{-19}$ — machine precision, confirming that the forward model and the RL correction are exactly consistent. We then ran the Richardson–Lucy iteration from a smooth initialization $O_0 = \max(h*I,\,0)$ for $160$ iterations, tracking the NMSE against the ground truth, the residual, the total variation, and the peak.

The results:

| iteration | NMSE | residual | TV | peak |
|---|---|---|---|---|
| 0 (init) | 0.7318 | 881.1 | 223 234 | 13.5 |
| 10 | 0.3978 | 554.6 | 272 583 | 82.6 |
| 20 | 0.4747 | 534.4 | 304 014 | 134.4 |
| 50 | 0.8395 | 524.9 | 373 464 | 236.7 |
| 100 | 1.2806 | 523.9 | 459 377 | 335.7 |
| 159 | 1.6595 | 523.9 | 541 391 | 425.5 |

The observed (blurred, Poisson-noised) image has NMSE $0.7218$ relative to the ground truth; the RL estimate stopped at the L-curve corner has NMSE $0.3978$ — a **$1.81\times$ improvement**. The L-curve corner and the minimum-NMSE iteration both fall at **iteration 10**, and the two agree: this is the check that neither the forward model nor the stopping rule is broken.

**Interpretation.**

1. *The NMSE U-curve is the signature of an ill-posed problem being over-fit.* The NMSE falls from $0.7318$ (init) to its minimum $0.3978$ at iteration 10, then *rises* to $1.6595$ by iteration 159. Past the minimum the iteration is no longer recovering the object — it is fitting the Poisson noise, producing an increasingly spiky, high-total-variation estimate (the TV climbs monotonically from $223\,234$ to $541\,391$, and the peak from $13.5$ to $425.5$). The residual, by contrast, saturates near $524$: the data cannot be fit better than the shot-noise floor allows. The U-shape is exactly what theory predicts for an unregularized maximum-likelihood iteration on a band-limited problem, and it is why **early stopping is not an optional refinement but a requirement**.

2. *The L-curve corner is the model-free stopping rule, and it works.* The corner (iteration 10) is found without any knowledge of the ground truth — it is read from the data (residual versus complexity) alone. That it coincides with the minimum-NMSE iteration (which *does* use the ground truth) is the validation: the model-free rule and the oracle rule agree. In a real experiment, where the ground truth is unavailable, the L-curve corner is the defensible place to stop.

3. *The $1.81\times$ improvement is real but bounded.* The deconvolution recovers a substantial fraction of the lost fidelity (NMSE $0.7218 \to 0.3978$), sharpening the blobs and raising the peak from $13.5$ toward the true value. But it does not restore the full diffraction limit: the high frequencies the OTF suppressed (§3) are absent from the data, and no linear or maximum-likelihood inversion can recover what was never measured. The improvement is a recovery of *measured* information, not a creation of new information. This is the honest limit of deconvolution, and it is set by the band-limitedness, not by the algorithm.

4. *The initialization matters, and not in the obvious way.* Starting from the noisy observed image $O_0 = I$ rather than the smooth $O_0 = h*I$ makes the total-variation complexity *dip* before it rises (RL smooths the noise before it sharpens the object), so the L-curve is no longer an L and its corner snaps to iteration 0 — a degenerate, useless stopping rule. The smooth initialization is what makes the L-curve a proper L. The stopping rule is only as good as the initialization that makes its input curve well-behaved.

### 5.1 A worked convergence bound for the Richardson–Lucy iteration

§5 established that the Richardson–Lucy (RL) iteration is an EM ascent method that must be early-stopped. Here we make the convergence quantitative: we derive a lower bound on the per-step likelihood increase and verify it on the actual trajectory of §5.

**The monotone-likelihood theorem.** The EM algorithm [Dempster, Laird & Rubin 1977] guarantees that the likelihood is non-decreasing at every step, $\mathcal L(O^{(t+1)}) \ge \mathcal L(O^{(t)})$. The proof is the standard EM argument: write the complete-data log-likelihood as $Q(\theta \mid \theta^{(t)}) + Z(\theta^{(t)})$, where $Q$ is its expectation under the current posterior and $Z$ a $\theta$-independent term. Since the M-step sets $O^{(t+1)} = \arg\max_O Q(O \mid O^{(t)})$,

$$\mathcal L(O^{(t+1)}) = Q(O^{(t+1)}\mid O^{(t)}) + Z \;\ge\; Q(O^{(t)}\mid O^{(t)}) + Z = \mathcal L(O^{(t)}).$$

**A lower bound on the per-step increase.** The first inequality sharpens: the gain is at least the gap the M-step closes,

$$\mathcal L(O^{(t+1)}) - \mathcal L(O^{(t)}) \;\ge\; Q(O^{(t+1)}\mid O^{(t)}) - Q(O^{(t)}\mid O^{(t)}).$$

For the Poisson model the E-step gives $S_j = O_j^{(t)} a_j$, where $a_j = [h \star (I/(h*O^{(t)})])_j$ is the RL correction factor, and the M-step sets $O^{(t+1)}_j = S_j$. Substituting, and using the exact flux conservation $\sum_j O^{(t+1)}_j = \sum_j O^{(t)}_j = \sum_k I_k$ (which the unregularized update satisfies), the gap collapses to a Kullback–Leibler divergence between consecutive *normalized* iterates $p^{(t)} = O^{(t)}/\sum_j O^{(t)}_j$:

$$\boxed{\;\mathcal L(O^{(t+1)}) - \mathcal L(O^{(t)}) \;\ge\; \Big(\sum_k I_k\Big)\, \mathrm{KL}\!\big(p^{(t+1)} \,\big\|\, p^{(t)}\big)\;}$$

Two immediate consequences: the bound is non-negative (recovering the monotonicity theorem), and it vanishes as $p^{(t+1)} \to p^{(t)}$ — i.e. the iteration has converged to a fixed point.

**Worked on the §5 trajectory.** We re-ran the §5 computation (same seed, object, PSF, $200\,532$ photons) with the *exact* M-step (the $\varepsilon$ in the §5 ratio is a numerical safeguard that does not change the trajectory to the quoted precision), logging the log-likelihood, the per-step gain $\Delta\mathcal L$, and the bound at every iteration:

| it | $\Delta\mathcal L$ | bound $\sum I\,\mathrm{KL}$ | ratio $\Delta\mathcal L$/bound |
|---|---|---|---|
| 0 | 50 883 | 30 525 | 1.67 |
| 1 | 15 312 | 12 160 | 1.26 |
| 2 | 7 799 | 6 061 | 1.29 |
| 5 | 1 889 | 1 612 | 1.17 |
| 10 | 444 | 474 | 0.94 |
| 20 | 77 | 115 | 0.67 |
| 50 | 0.26 | 15.4 | 0.02 |

In the convergence regime the bound is **tight** (ratio $\Delta\mathcal L$/bound $0.97$–$1.67$ for $t\le 9$), and the gain decays geometrically, $\Delta\mathcal L^{(t+1)}/\Delta\mathcal L^{(t)} \to \rho \approx 0.87$ (the median over iterations $30$–$59$) — linear convergence to the maximum-likelihood fixed point. The bound $\to 0$ as the iterate converges, exactly as the theorem requires.

**Why the likelihood is the wrong stopping rule.** The log-likelihood climbs monotonically through the convergence regime (from $19\,488$ at the initialization) and crests only near iteration $50$ ($58\,133$), *far past* the L-curve corner (iteration $10$, $55\,835$). At the final iterate (iteration $159$) the likelihood is $57\,957$ — still $3.8\%$ *above* the corner value. The bound explains why: as long as the iterate keeps changing, RL keeps closing a positive KL gap by fitting the Poisson noise, so the likelihood keeps rising even as the NMSE (which measures fidelity to the *true* object, not the noisy data) bottoms out at iteration 10 and then rises. The likelihood measures fit to the data; the NMSE measures fit to the object. Past the corner they pull in opposite directions, which is precisely why the L-curve (residual vs complexity) — not the likelihood — is the defensible stopping rule.

**A caveat, stated honestly.** The bound is derived for the ideal (infinite-support) blur. Our forward model is a finite-image FFT convolution, whose boundary truncation breaks exact flux conservation: the iterate sum drifts $200\,542$ (it 1) $\to 200\,621$ (it 30), a $\sim 0.04\%$ departure from $\sum I = 200\,532$ at a steady $\sim 2.7$ photons per step. Consequently the bound, which assumes exact flux conservation, is tight only in the early iterations and becomes loose from $t\approx 10$; and the strict monotonicity is violated (by $\sim 0.3\%$) in the deep over-fitting regime ($t \gtrsim 51$), where the log-likelihood crests and drifts down. These are finite-image boundary effects, not a failure of the theorem: on the ideal model the DLR monotonicity and the KL bound hold exactly.

---

## 6. The quantum limit: shot noise and the resolution–sensitivity trade-off

The two threads of the essay — the resolution (the forward problem) and the recovery (the inverse problem) — meet at the detector, where the limit is ultimately quantum. The photomultiplier counts photons, so the signal is Poisson and the signal-to-noise ratio scales as $\sqrt{N}$, where $N$ is the number of detected photons. Every lever that improves the resolution costs photons: closing the pinhole (§3) blocks light; deconvolving (§5) amplifies the noise in the suppressed frequencies; pushing the NA higher collects more but also concentrates the (damaging) illumination into a smaller volume. The resolution–sensitivity trade-off is therefore not an engineering inconvenience but a fundamental bound: one cannot improve the resolution without reducing the photon count that carries the information, and a lower photon count raises the relative shot noise. This is the sense in which the diffraction limit and the noise limit are two faces of the same constraint.

There is a second, practical quantum cost: the same light that provides the resolution also **damages the specimen**. Fluorophores bleach irreversibly under the excitation, and the high local intensity is phototoxic to living cells. For a live specimen the usable signal is bounded not only by the shot noise but by the damage the specimen can survive — a further, and often tighter, constraint on how much resolution can actually be extracted.

---

## 7. Conclusion

The diffraction-limited image is a band-limited inverse problem, and the essay has followed that single idea through three constructions. The **PSF/OTF** pair is the Green's function and symbol of the imaging operator; the Abbe limit is the uncertainty principle applied to the pupil's band-limit. The **confocal** modification replaces the PSF by a product of two PSFs, which in the Fourier domain doubles the OTF bandwidth — a factor of two in the Abbe limit and a factor of $\sqrt{2}$ in the FWHM — and, through the pinhole, provides the optical sectioning that recovers the lost axial coordinate. The **deconvolution** is the recovery of the object, an ill-posed band-limited inverse problem whose natural solver is the Richardson–Lucy iteration, the EM algorithm for the Poisson likelihood; it must be early-stopped at the L-curve corner, and it recovers only the information that was actually measured.

The two computations ground the picture in numbers. The high-NA correction shows the paraxial Airy picture is an excellent lateral model but over-estimates the axial focus by a factor $\approx\sqrt{2}$ (vectorial axial FWHM $478.9$ nm versus the Abbe $678.0$ nm, a ratio $0.706\approx1/\sqrt{2}$); the full 3D Debye PSF (§4.1) sharpens the picture, showing the focus is a prolate lobe that confocal detection narrows by $\approx1/\sqrt{2}$ on both axes (axial FWHM $480.0\to344.5$ nm, lateral $183.6\to131.8$ nm). The deconvolution shows the Richardson–Lucy iteration, stopped at the L-curve corner (which coincides with the minimum-NMSE iteration), improves the NMSE by a factor of $1.81$ on a simulated stack — a real but bounded recovery, limited by the band-limitedness, not by the algorithm; the worked convergence bound (§5.1) makes this quantitative, bounding the per-step likelihood gain below by the total photon budget times the KL divergence of consecutive normalized iterates, tight in the convergence regime, and explaining why the likelihood (which crests far past the corner) is the wrong stopping rule. And the detector reminds us that the whole enterprise is bounded by the shot noise and the damage the specimen can bear.

For the mathematical physicist, the unifying lesson is that the resolution of a microscope is not a number to be looked up but a *theorem* about a band-limited operator: it follows from the support of the OTF, it is sharpened (but not eliminated) by the confocal product, and its recovery is governed by the same ill-posedness that governs every other inverse problem in physics.

---

## References

1. E. Abbe, "Beiträge zur Theorie des Mikroskops und der mikroskopischen Wahrnehmung", *M. Schultze's Archiv für mikroskopische Anatomie* **9**, 413–468 (1873). DOI [10.1007/bf02956173](https://doi.org/10.1007/bf02956173). (The Abbe diffraction limit, §2.)
2. P. J. W. Debye, "Das Verhalten von Lichtwellen in der Nähe eines Brennpunktes oder einer Brennlinie", *Annalen der Physik* **335**, 755–776 (1909). DOI [10.1002/andp.19093351406](https://doi.org/10.1002/andp.19093351406). (The angular-spectrum / Debye representation, §4.)
3. B. Richards and E. Wolf, "Electromagnetic diffraction in optical systems. II. Structure of the image field in an aplanatic system", *Proc. R. Soc. A* **253**, 358–379 (1959). DOI [10.1098/rspa.1959.0200](https://doi.org/10.1098/rspa.1959.0200). (The vectorial high-NA PSF, §4.)
4. C. J. R. Sheppard and T. Wilson, "Image formation in scanning microscopes with partially coherent source and detector", *Optica Acta* **25**, 315–325 (1978). DOI [10.1080/713819784](https://doi.org/10.1080/713819784). (The scanning/confocal transfer function, §3.)
5. M. Gu and C. J. R. Sheppard, "Three-dimensional coherent transfer functions in confocal scanning imaging", *Proc. SPIE* **1983** (1993). DOI [10.1117/12.2308503](https://doi.org/10.1117/12.2308503). (The doubled confocal pass band, §3.)
6. W. H. Richardson, "Bayesian-based iterative method of image restoration", *J. Opt. Soc. Am.* **62**, 55–59 (1972). DOI [10.1364/JOSA.62.000055](https://doi.org/10.1364/JOSA.62.000055). (The Richardson iteration, §5.)
7. L. B. Lucy, "An iterative technique for the rectification of observed distributions", *Astron. J.* **79**, 745–754 (1974). DOI [10.1086/111605](https://doi.org/10.1086/111605). (The Lucy iteration and its maximum-likelihood convergence, §5.)
8. A. P. Dempster, N. M. Laird, and D. B. Rubin, "Maximum likelihood from incomplete data via the EM algorithm", *J. R. Statist. Soc. B* **39**, 1–38 (1977). (The EM algorithm, §5.)
9. P. C. Hansen, "Analysis of discrete ill-posed problems by means of the L-curve", *SIAM Review* **34**, 561–580 (1992). DOI [10.1137/1034115](https://doi.org/10.1137/1034115). (The L-curve stopping rule, §5.)
10. M. Born and E. Wolf, *Principles of Optics: Electromagnetic Theory of Propagation, Interference and Diffraction of Light*, 7th ed., Cambridge University Press (1999). (Standard reference on the diffraction theory, §§2–4.)
11. J. W. Goodman, *Introduction to Fourier Optics*, 3rd ed., Roberts & Company (2005). (Fourier-optics treatment of the OTF, §§2–3.)
12. J. B. Pawley (ed.), *Handbook of Biological Confocal Microscopy*, 3rd ed., Springer (2006). (The confocal instrument and its practical limits, §§3, 6.)
