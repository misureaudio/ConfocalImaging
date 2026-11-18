#!/usr/bin/env python
"""
Part (a): Richards-Wolf high-NA axial PSF vs the paraxial (Abbe) approximation.

On-axis field of a high-NA focused beam (direction-cosine / angular-spectrum form):

  scalar  Richards-Wolf :  U_s(z) = -i k_med ∫_{cosα}^{1} u^{1/2} e^{i k_med z u} du
  vectorial on-axis     :  Ex(z)  = -i k_med/2 ∫_{cosα}^{1} (u^{1/2}+u^{5/2}) e^{i k_med z u} du
  paraxial (Abbe) limit : |U_par(z)|^2 = (4/z^2) sin^2(k_med z α^2 / 4)

The integrals are evaluated with a fixed u-grid + trapezoidal rule, vectorised
over the whole z-grid (fast).  The z->0 values are checked against the closed
form ∫ u^p du = (1-cosα^{p+1})/(p+1).

Output: richards_wolf_axial.png  (axial + lateral PSF comparison)
"""
import numpy as np
from scipy.special import j0, j1
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------- parameters -----------------------------------
lam   = 500e-9                 # vacuum wavelength [m]
n0    = 1.5                    # immersion (oil) refractive index
NA    = 1.4                    # numerical aperture
k0    = 2.0*np.pi/lam          # vacuum wavenumber
kmed  = n0*k0                  # wavenumber in the medium
alpha = np.arcsin(NA/n0)       # semi-aperture angle [rad]
cosA  = np.cos(alpha)

print(f"lambda = {lam*1e9:.1f} nm | n0 = {n0} | NA = {NA}")
print(f"alpha  = {np.degrees(alpha):.2f} deg | cos(alpha) = {cosA:.4f}")
print(f"k_med  = {kmed:.4e} 1/m | k0*NA = {k0*NA:.4e} 1/m")

# --------------------- on-axis integrals (vectorised) ------------------------
Nu = 4001
u  = np.linspace(cosA, 1.0, Nu)
du = u[1] - u[0]

zmax = 800e-9
Nz   = 801
z    = np.linspace(-zmax, zmax, Nz)
beta = kmed*z                       # (Nz,)

def I_p(p):
    """I_p(z) = ∫_{cosα}^{1} u^p e^{i k_med z u} du  for every z (vectorised)."""
    phase = np.exp(1j*np.outer(beta, u))          # (Nz, Nu)
    return np.sum((u**p)[None, :] * phase, axis=1) * du

I05 = I_p(0.5)
I25 = I_p(2.5)
Us  = -1j*kmed*I05                 # scalar on-axis field
Ev  = -0.5j*kmed*(I05 + I25)       # vectorial on-axis Ex field

# closed-form check at z = 0  (z[0] is -zmax; the on-axis point is z[Nz//2])
z0 = Nz//2
I05_0 = (1.0 - cosA**1.5)/1.5     # ∫ u^{1/2} du = (1-cosα^{3/2})/(3/2)
I25_0 = (1.0 - cosA**3.5)/3.5     # ∫ u^{5/2} du = (1-cosα^{7/2})/(7/2)
print(f"\n[check z=0]  I(1/2): quad={I05[z0]:.6f}  closed={I05_0:.6f}  "
      f"diff={abs(I05[z0]-I05_0):.2e}")
print(f"[check z=0]  I(5/2): quad={I25[z0]:.6f}  closed={I25_0:.6f}  "
      f"diff={abs(I25[z0]-I25_0):.2e}")

# paraxial (Abbe) axial intensity — small-alpha limit of the SAME integral,
# expressed in terms of the measurable NA (alpha -> NA/n0 in the paraxial limit)
bpar = kmed*(NA/n0)**2
Ipar = np.empty_like(z)
nz   = np.abs(z) > 1e-15
Ipar[nz]  = 4.0/z[nz]**2 * np.sin(bpar*z[nz]/4.0)**2
Ipar[~nz] = bpar**2/4.0

# normalised intensities
Is  = np.abs(Us)**2
Iv  = np.abs(Ev)**2
Is_n, Iv_n, Ipar_n = Is/Is.max(), Iv/Iv.max(), Ipar/Ipar.max()

# ----------------------------- lateral PSF ----------------------------------
Nr = 501
r  = np.linspace(0.0, 500e-9, Nr)
Nth = 2001
th = np.linspace(0.0, alpha, Nth)
dth = th[1] - th[0]

# full scalar lateral field: U(r) = -i kmed ∫_0^α sqrt(cosθ) J0(kmed r sinθ) sinθ dθ
arg   = np.outer(r, np.sin(th))*kmed          # (Nr, Nth)
J0m   = j0(arg)
Ulat  = -1j*kmed*np.sum(J0m*(np.sqrt(np.cos(th))*np.sin(th))[None, :], axis=1)*dth
Ilat_full = np.abs(Ulat)**2

# paraxial (Airy) lateral intensity: [2 J1(rho)/rho]^2, rho = k0 NA r
rho = k0*NA*r
Ilat_par = np.ones_like(r)
m = np.abs(rho) > 1e-12
Ilat_par[m] = (2*j1(rho[m])/rho[m])**2

Ilat_full_n = Ilat_full/Ilat_full.max()
Ilat_par_n  = Ilat_par/Ilat_par.max()

# ----------------------------- FWHM helper ----------------------------------
def fwhm_nm(xnm, y):
    y = y/y.max()
    ipk = int(np.argmax(y))
    half = 0.5
    xl = xnm[ipk]
    for i in range(ipk-1, -1, -1):
        if y[i] < half:
            xl = np.interp(half, [y[i], y[i+1]], [xnm[i], xnm[i+1]]); break
    xr = xnm[ipk]
    for i in range(ipk+1, len(xnm)):
        if y[i] < half:
            xr = np.interp(half, [y[i], y[i-1]], [xnm[i], xnm[i-1]]); break
    return xr - xl

znm = z*1e9
rnm = r*1e9
print("\n--- axial FWHM (nm) ---")
print(f"  paraxial (Abbe) : {fwhm_nm(znm, Ipar_n):8.1f}")
print(f"  full scalar RW  : {fwhm_nm(znm, Is_n):8.1f}")
print(f"  full vectorial  : {fwhm_nm(znm, Iv_n):8.1f}")
print("--- lateral FWHM (nm) ---")
# NOTE: the lateral profile is radial (r >= 0) with its peak at r = 0 (the grid
# edge), so fwhm_nm returns the HWHM (radius at half-max); the full width is 2x.
print(f"  paraxial (Airy) : {2*fwhm_nm(rnm, Ilat_par_n) :8.1f}")
print(f"  full scalar RW  : {2*fwhm_nm(rnm, Ilat_full_n):8.1f}")

# ----------------------------- plot -----------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))

ax1.plot(znm, Ipar_n, lw=2.0, color="0.45", ls="--", label="paraxial (Abbe)  $|U_{par}|^2$")
ax1.plot(znm, Is_n,   lw=2.0, color="tab:blue",  label="full scalar Richards-Wolf  $|U_s|^2$")
ax1.plot(znm, Iv_n,   lw=2.0, color="tab:red",   label="full vectorial on-axis  $|E_x|^2$")
ax1.set_xlabel("axial position $z$  [nm]")
ax1.set_ylabel("normalised intensity")
ax1.set_title(f"Axial PSF   (NA={NA}, $n_0$={n0}, $\\lambda$={lam*1e9:.0f} nm)")
ax1.set_xlim(-450, 450)
ax1.axvline(0, color="k", lw=0.5, alpha=0.4)
ax1.grid(alpha=0.3)
ax1.legend(fontsize=9)

ax2.plot(rnm, Ilat_par_n,  lw=2.0, color="0.45", ls="--", label="paraxial (Airy)  $[2J_1/\\rho]^2$")
ax2.plot(rnm, Ilat_full_n, lw=2.0, color="tab:blue", label="full scalar Richards-Wolf")
ax2.set_xlabel("radial position $r$  [nm]")
ax2.set_ylabel("normalised intensity")
ax2.set_title("Lateral PSF (same parameters)")
ax2.set_xlim(0, 350)
ax2.grid(alpha=0.3)
ax2.legend(fontsize=9)

fig.suptitle("Richards-Wolf high-NA PSF vs paraxial (Abbe) approximation", y=1.02, fontsize=12)
fig.tight_layout()
out = "richards_wolf_axial.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"\nsaved {out}")
