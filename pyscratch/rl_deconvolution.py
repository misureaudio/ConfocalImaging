#!/usr/bin/env python
"""
Part (b): 3D Richardson-Lucy deconvolution of a simulated confocal stack.

  * 3D confocal PSF from the scalar solid-immersion (Debye) integral
        U(r,z) = ∫_0^{kmax} J0(kρ r) e^{i kz z} kρ dkρ ,  kz = sqrt(kmed²-kρ²)
    with the confocal (illumination × detection) PSF  h = |U|⁴,  unit-normalised.
  * A 3D object (4 point sources) is blurred with h and corrupted by Poisson noise.
  * The 3D Richardson-Lucy (Poisson maximum-likelihood) update is run:
        O ← O · [ h  ⋆  ( I / (h * O) ) ]
  * We track the L-curve (residual norm vs total-variation) and the true NMSE
    (a U-shape), locate the L-curve corner, and plot before/after z-slices.

Output: rl_deconvolution.png
"""
import numpy as np
from scipy.special import j0
from scipy.interpolate import RegularGridInterpolator
from scipy.signal import fftconvolve
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rng = np.random.default_rng(42)

# ----------------------------- parameters -----------------------------------
lam  = 500e-9
n0   = 1.5
NA   = 1.4
k0   = 2.0*np.pi/lam
kmed = n0*k0
kmax = 2.0*np.pi*NA/lam            # max transverse wavenumber

Nx = Ny = 64
Nz = 32
dx = dy = 10e-9
dz = 20e-9
cx, cy, cz = Nx//2, Ny//2, Nz//2
N_photons = 200000
max_iter  = 160

# --------------------- 3D confocal PSF (Debye integral) ----------------------
Nr, Nzp, Nk = 201, 401, 500
r_grid = np.linspace(0.0, cx*dx, Nr)
z_grid = np.linspace(-cz*dz, (Nz-1-cz)*dz, Nzp)
k_rho  = np.linspace(0.0, kmax, Nk)
dk     = k_rho[1]-k_rho[0]
kz     = np.sqrt(np.maximum(kmed**2 - k_rho**2, 0.0))

U = np.zeros((Nr, Nzp), dtype=complex)
for k in range(Nk):
    U += np.outer(j0(k_rho[k]*r_grid), np.exp(1j*kz[k]*z_grid))*k_rho[k]*dk

h_conf = np.abs(U)**4
h_conf /= h_conf.sum()

interp = RegularGridInterpolator((r_grid, z_grid), h_conf,
                                 method='linear', bounds_error=False, fill_value=0.0)
XX, YY = np.meshgrid(np.arange(Nx), np.arange(Ny), indexing='ij')
r3d = np.sqrt((XX-cx)**2*dx**2 + (YY-cy)**2*dy**2)
R3 = np.tile(r3d[:,:,None], (1,1,Nz))
Z3 = np.tile(((np.arange(Nz)-cz)*dz)[None,None,:], (Nx,Ny,1))
h3d = interp(np.stack([R3.ravel(), Z3.ravel()], axis=1)).reshape(Nx, Ny, Nz)
h3d = np.maximum(h3d, 0.0)
h3d /= h3d.sum()
print(f"PSF: sum={h3d.sum():.6f}  peak={h3d.max():.3e}")

# ---- explicit FFT forward-model / cross-correlation (exact centering) ----
# Keep the PSF at its natural centre in a 2x-padded array; the circular conv
# then equals the linear conv (64+64-1=127 < 128).  Crop with offset (N//2) so
# the blur response is centred on each object pixel.  For the symmetric PSF the
# RL correction (h ⋆ R) equals the convolution (h * R), so both use H_fft.
fft_shape = (2*Nx, 2*Ny, 2*Nz)
h3d_pad = np.zeros(fft_shape); h3d_pad[:Nx,:Ny,:Nz] = h3d
H_fft   = np.fft.fftn(h3d_pad)
crop    = (slice(Nx//2, Nx//2+Nx), slice(Ny//2, Ny//2+Ny), slice(Nz//2, Nz//2+Nz))

def _pad(a):
    p = np.zeros(fft_shape); p[:Nx,:Ny,:Nz] = a; return p

def _conv(a):
    return np.real(np.fft.ifftn(np.fft.fftn(_pad(a)) * H_fft))[crop]

def forward(O):
    """(h * O) cropped to O's size — the forward (blur) model."""
    return _conv(O)

def xcorr(R):
    """(h ⋆ R) cropped to R's size — the RL correction factor (h symmetric)."""
    return _conv(R)

# --------------------- forward-model self-test -------------------------------
O_t = np.zeros((Nx, Ny, Nz)); O_t[cx, cy, cz] = 1.0
blurred_t = forward(O_t)
print(f"[test] point-source blur vs PSF:  max|diff| = {np.max(np.abs(blurred_t-h3d)):.2e}")

# ----------------------------- object ----------------------------------------
# A realistic confocal "object": 4 small Gaussian blobs (sigma ~ 3 px), like
# fluorescent beads / protein clusters — sharp enough to deconvolve, smooth
# enough that the noisy RL can actually approach them.
sigma_pix = 3.0
XXg, YYg, ZZg = np.meshgrid(np.arange(Nx), np.arange(Ny), np.arange(Nz), indexing='ij')
O_true = np.zeros((Nx, Ny, Nz))
for (px, py, pz) in [(16,16,6),(48,16,6),(16,48,26),(48,48,26)]:
    O_true += np.exp(-(((XXg-px)**2 + (YYg-py)**2 + (ZZg-pz)**2))/(2*sigma_pix**2))

B     = np.maximum(forward(O_true), 0.0)
scale = N_photons / B.sum()
O_true_s = O_true*scale          # true object on the data (photon) scale
B_s     = B*scale                # = H O_true_s  (noiseless mean)
I       = rng.poisson(B_s)       # noisy data
print(f"data: observed={I.sum():.0f}  expected={B_s.sum():.0f} photons")

# ----------------------------- RL deconvolution ------------------------------
eps = 1e-12
O2s = float((O_true_s**2).sum())

def tv3d(x):
    return (np.abs(np.diff(x,1)).sum()+np.abs(np.diff(x,0)).sum()+np.abs(np.diff(x,2)).sum())

def rl_step(O):
    blurred = forward(O)
    R       = I/(blurred+eps)
    corr    = xcorr(R)
    return np.maximum(O*corr, 0.0), blurred

# pass 1: track the diagnostics
# Smooth init O0 = h*I (blurred data) so the complexity rises monotonically
# and the L-curve is a proper L (a noisy init O0=I makes the TV dip first).
resid = np.zeros(max_iter); tv = np.zeros(max_iter); nmse = np.zeros(max_iter)
O = np.maximum(forward(I), 0.0)
for t in range(max_iter):
    O, blurred = rl_step(O)
    resid[t] = np.linalg.norm(I - blurred)
    tv[t]    = tv3d(O)
    nmse[t]  = float(((O-O_true_s)**2).sum())/O2s
    if t in (0,1,2,5,10,20,50,100,max_iter-1):
        print(f"  it {t:3d}:  NMSE={nmse[t]:.5f}  resid={resid[t]:9.1f}  "
              f"TV={tv[t]:10.0f}  peak={O.max():9.1f}  sum={O.sum():9.0f}")

# L-curve corner: max distance from the chord in normalised log-log space
x  = np.log(tv);  y = np.log(resid)
xn = (x-x.min())/(x.max()-x.min()); yn = (y-y.min())/(y.max()-y.min())
x0, y0 = xn[0], yn[0];  x1, y1 = xn[-1], yn[-1]
a = y1-y0;  b = x0-x1;  c = y0*(x1-x0) - x0*(y1-y0)   # chord: a x + b y + c = 0
L  = np.hypot(a, b)
dist = np.abs(a*xn + b*yn + c)/L
t_corner = int(np.argmax(dist))
t_best   = int(np.argmin(nmse))
print(f"L-curve corner @ iter {t_corner}   min-NMSE @ iter {t_best}")

# pass 2: recover the deconvolved object at the corner
O = np.maximum(forward(I), 0.0); O_dec = None
for t in range(t_corner+1):
    O, _ = rl_step(O)
    if t == t_corner:
        O_dec = O.copy()

nmse_blur  = float(((I - O_true_s)**2).sum())/O2s
nmse_dec   = float(((O_dec - O_true_s)**2).sum())/O2s
print(f"NMSE  blurred={nmse_blur:.4f}   deconvolved={nmse_dec:.4f}   "
      f"improvement={nmse_blur/nmse_dec:.2f}x")

# ----------------------------- plots -----------------------------------------
def imshow(ax, arr, title, cmap='inferno'):
    ax.imshow(arr, cmap=cmap)
    ax.set_title(title, fontsize=9)
    ax.axis('off')

fig = plt.figure(figsize=(13, 11))
gs  = fig.add_gridspec(4, 2, hspace=0.55, wspace=0.25)

# (0,0) L-curve
axL = fig.add_subplot(gs[0,0])
axL.plot(tv, resid, '-o', ms=2.5, lw=1.5, color='tab:blue')
axL.plot(tv[t_corner], resid[t_corner], 'r*', ms=16, label=f'L-curve corner (it {t_corner})')
axL.set_xscale('log'); axL.set_yscale('log')
axL.set_xlabel('total variation  (solution complexity)')
axL.set_ylabel('residual  $\|I - h*O\|_2$')
axL.set_title('Richardson-Lucy L-curve')
axL.legend(fontsize=8); axL.grid(alpha=0.3, which='both')

# (0,1) NMSE U-curve
axN = fig.add_subplot(gs[0,1])
axN.plot(nmse, '-o', ms=2.5, lw=1.5, color='tab:green')
axN.plot(t_best,   nmse[t_best],   'g^', ms=12, label=f'min NMSE (it {t_best})')
axN.plot(t_corner, nmse[t_corner], 'r*', ms=16, label=f'L-curve corner (it {t_corner})')
axN.set_xlabel('RL iteration')
axN.set_ylabel('NMSE vs true object')
axN.set_title('True error (U-shape) — overfitting after the knee')
axN.legend(fontsize=8); axN.grid(alpha=0.3)

# rows 1-3: true / blurred+noisy / deconvolved, at z=6 and z=26
rows = [
    (1, 'ground truth (4 point sources)'),
    (2, 'blurred + Poisson noise  (the data $I$)'),
    (3, f'Richardson-Lucy deconvolved (it {t_corner})'),
]
zs = [6, 26]
data = [O_true_s, I, O_dec]
for ri,(row, rtitle) in enumerate(rows):
    for ci, z in enumerate(zs):
        ax = fig.add_subplot(gs[row, ci])
        arr = data[ri][:, :, z]
        arr = arr/arr.max() if arr.max() > 0 else arr
        imshow(ax, arr, f'{rtitle}\n  z = {z*dx*1e9:.0f} nm')
        if ri == 0:
            for (px,py,pz) in [(16,16,6),(48,16,6),(16,48,26),(48,48,26)]:
                if pz == z:
                    ax.plot(py, px, 'c+', ms=10, mew=1.5)

fig.suptitle('3D Richardson-Lucy deconvolution of a simulated confocal stack '
             f'(NA={NA}, $\\lambda$={lam*1e9:.0f} nm)', y=0.995, fontsize=13)
out = 'rl_deconvolution.png'
fig.savefig(out, dpi=150, bbox_inches='tight')
print(f'saved {out}')
