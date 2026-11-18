Yes — there are **five clear errors** (mathematical, numerical, and self-contradictory) in the essay, along with a few minor physical/terminological inaccuracies. Here is the breakdown:

---

### 1. §4 vs. §4.1: Factor-of-2 error in the lateral FWHMs (HWHM mislabeled as FWHM)
In **§4**, the lateral FWHM table reports:
* **paraxial (Airy):** $91.9\text{ nm}$
* **full scalar:** $88.4\text{ nm}$

For $\lambda = 500\text{ nm}$ and $\mathrm{NA} = 1.4$, the Airy FWHM is
$$\delta_\perp = \frac{1.6163}{\pi}\,\frac{\lambda}{\mathrm{NA}} = 0.5145 \times \frac{500}{1.4} = 183.74\text{ nm},$$
which **§4.1** correctly reports as $183.7\text{ nm}$ (with single-path $|U|^2$ at $183.6\text{ nm}$). The §4 value $91.9\text{ nm} = 183.74/2$ is the **HWHM** (half-width / radius at half-maximum), not the FWHM; likewise $88.4\text{ nm}$ is the HWHM of the scalar lateral profile ($176.8\text{ nm}$ FWHM). Consequently, §4 and §4.1 contradict each other by a factor of 2 on the same paraxial Airy lateral FWHM.

*(Related minor note in §4.1: "first zero at $435.7\text{ nm}$" gives the diameter $2r_0$, whereas the Rayleigh first-zero radius is $r_0 = 0.61\lambda/\mathrm{NA} = 217.8\text{ nm}$. Calling it the "first-zero diameter" avoids ambiguity.)*

---

### 2. §3: Wrong root for the confocal Airy FWHM ($1.140$ vs. $1.160$)
In **§3**, the text states:
> *"the FWHM of $[2J_1/\rho]^4$ — which reaches half its maximum where $[2J_1/\rho]^2 = 1/\sqrt{2}$ — is at $\rho \approx 1.140$. The ratio $1.140/1.616 \approx 0.705 \approx 1/\sqrt{2}$"*

Solving $[2J_1(\rho)/\rho]^4 = 1/2$ (i.e., $2J_1(\rho)/\rho = 2^{-1/4} \approx 0.8409$) actually gives **$\rho \approx 1.160$**, not $1.140$ (at $\rho = 1.140$, $[2J_1/\rho]^4 \approx 0.513$). The exact Airy-fourth-power ratio is
$$\frac{1.1603}{1.6163} = 0.718,$$
not $0.705$ (and the equality $\delta_\perp^{\,\mathrm{conf}} = \delta_\perp^{\,\mathrm{wf}}/\sqrt{2}$ is only a Gaussian approximation, not exact). This also fixes an internal inconsistency with **§4.1**, which numerically finds $131.8/183.6 = 0.718$ and explicitly notes the excess above $1/\sqrt{2} = 0.707$.

---

### 3. §4.1: Wrong confocal OTF cutoff ($2k_{\max}$ vs. $4k_{\max}$) and confused OTF/missing-cone statement
In the final paragraph of **§4.1**:
> *"In $k$-space the PSF is supported only on the cone $k_\rho \le k_{\max}$ ... the single-path OTF has transverse support out to $2k_{\max}$, and the confocal OTF (the product of two pupils, §3) reaches the same $2k_{\max}$, doubling the Abbe bandwidth."*

There are three problems here:
1. **Cutoff error:** As **§3** correctly derives, the single-path incoherent OTF reaches $2k_{\max}$ (where $k_{\max} = k_\perp = 2\pi\mathrm{NA}/\lambda$), and the **confocal** OTF reaches **$4k_{\max}$**, not "the same $2k_{\max}$."
2. **Operation error:** The confocal OTF is the **convolution of two single-path OTFs** (each of which is a pupil autocorrelation), not "the product of two pupils."
3. **Missing-cone conflation:** It is the 3D Fourier transform of the coherent *field* $U$ (supported on the Ewald cap $k_\rho^2 + k_z^2 = k_{\mathrm{med}}^2,\; k_\rho \le k_{\max}$) that has transverse radius $k_{\max}$; the single-path OTF ($|U|^2$) has the bowtie support missing the cone of half-angle $90^\circ - \alpha = 21.0^\circ$ around the $k_z$-axis, while the 3D confocal OTF (convolution of two single-path OTFs) actually fills in the missing cone along $k_z$.

---

### 4. §5.1: Flawed proof of EM monotonicity
In **§5.1**, the proof states:
> *"write the complete-data log-likelihood as $Q(\theta \mid \theta^{(t)}) + Z(\theta^{(t)})$, where $Q$ is its expectation under the current posterior and $Z$ a $\theta$-independent term. Since the M-step sets $O^{(t+1)} = \arg\max_O Q(O \mid O^{(t)})$,*
> $$\mathcal L(O^{(t+1)}) = Q(O^{(t+1)}\mid O^{(t)}) + Z \;\ge\; Q(O^{(t)}\mid O^{(t)}) + Z = \mathcal L(O^{(t)})."$$

This is mathematically incorrect:
1. $\mathcal L(\theta)$ is the **incomplete-data (marginal) log-likelihood**, not the complete-data log-likelihood.
2. In the EM decomposition $\mathcal L(\theta) = Q(\theta \mid \theta^{(t)}) + H(\theta \mid \theta^{(t)})$, where $H(\theta \mid \theta^{(t)}) = -\mathbb E_{S\mid I,\theta^{(t)}}[\ln P(S \mid I,\theta)]$, the remainder $H(\theta \mid \theta^{(t)})$ is **not $\theta$-independent**. (If it were, maximizing $Q$ would maximize $\mathcal L$ in a single step.)

Monotonicity — and the per-step lower bound — holds because Gibbs' inequality gives
$$H(O^{(t+1)} \mid O^{(t)}) - H(O^{(t)} \mid O^{(t)}) = \mathrm{KL}\!\big(P(S \mid I, O^{(t)}) \,\big\|\, P(S \mid I, O^{(t+1)})\big) \ge 0,$$
so $\mathcal L(O^{(t+1)}) - \mathcal L(O^{(t)}) \ge Q(O^{(t+1)}\mid O^{(t)}) - Q(O^{(t)}\mid O^{(t)})$.

---

### 5. §5.1: Lower bound is *violated* (not "loose") for $t \ge 10$, and imprecise FFT boundary explanation
In the **§5.1** table:
* At $t=10$: $\Delta\mathcal L = 444 < 474$ (bound)
* At $t=20$: $\Delta\mathcal L = 77 < 115$ (bound)
* At $t=50$: $\Delta\mathcal L = 0.26 < 15.4$ (bound)
* For $t \gtrsim 51$: $\Delta\mathcal L < 0$ while the KL bound remains $> 0$.

1. **Terminology error:** When the actual gain $\Delta\mathcal L$ falls *below* a lower bound, the lower bound is **violated** (because an assumption of the bound failed), not "loose" (a lower bound is loose when $\Delta\mathcal L \gg \text{bound}$).
2. **Inaccurate diagnosis:** A periodic FFT convolution has no boundary truncation and conserves flux identically. If boundary-cropped linear convolution (`mode='same'`) is used instead, exact flux conservation $\sum_j O_j^{(t+1)} = \sum_k I_k$ (for $\varepsilon = 0$) *still* holds whenever the backward step uses the exact matrix transpose $h^T$; what fails without dividing by the boundary sensitivity factor $\sum_k h_{kj} < 1$ is the EM M-step condition (breaking monotonicity). Conversely, the observed drift in $\sum_j O_j$ ($200\,532 \to 200\,621$) indicates that the backward operator in the script is not the exact transpose of the forward operator (e.g., due to boundary cropping or parity/centering mismatch of an even-sized PSF).

---

### Minor inaccuracies worth cleaning up
* **§4.1 ("the same physics in different notation"):** Comparing the scalar Debye integral with weight $u^1$ ($481.6\text{ nm}$) to the vectorial Richards–Wolf integral with weight $u^{1/2}+u^{5/2}$ ($478.9\text{ nm}$), calling them "the same physics in different notation" is slightly inaccurate: one is scalar ($k_\rho\,dk_\rho$ angular-spectrum weight) and the other includes aplanatic apodization ($\sqrt{\cos\theta}$) plus vectorial polarization ($1+\cos^2\theta$). They are different models that happen to agree to $0.6\%$ on the axial FWHM.
* **§4 ("$678.0/478.9 = 1.415$"):** From the rounded table values, $678.0/478.9 = 1.4157 \approx 1.416$ (minor rounding discrepancy if computed from unrounded numbers).