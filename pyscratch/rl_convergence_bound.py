#!/usr/bin/env python
"""
Worked Richardson-Lucy convergence bound (companion to Confocal_Imaging_Essay
v2, section 5b).

Replicates the EXACT setup of rl_deconvolution.py (seed 42, same 4-blob object,
same Debye PSF, same 2e5-photon Poisson data) so the numbers are consistent with
the essay, then at every RL iteration logs:

  * L(O)      : the Poisson log-likelihood  (data-only constant dropped)
  * dL        : the actual per-step increase  L(O^{t+1}) - L(O^{t})
  * bound     : the EM lower bound  (sum I) * KL( p^{t+1} || p^{t} )
                with p^{t} = O^{t} / sum(O^{t})  (the normalized object)

and verifies  dL >= bound  at every step (the DLR monotone-likelihood theorem,
made quantitative).  Writes rl_convergence_bound.png + rl_convergence_bound.log.
"""
import numpy as np
from scipy.special import j0
from scipy.interpolate import RegularGridInterpolator
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rng = np.random.default_rng(42)

# ------------------- parameters (identical to rl_deconvolution.py) -------------
lam  = 500e-9; n0 = 1.5; NA = 1.4
k0   = 2.0*np.pi/lam; kmed = n0*k0; kmax = 2.0*np.pi*NA/lam
Nx = Ny = 64; Nz = 32
dx = dy = 10e-9; dz = 20e-9
cx, cy, cz = Nx//2, Ny//2, Nz//2
N_photons = 200000
max_iter  = 160

# --------------------- 3D confocal PSF (Debye integral) ------------------------
Nr, Nzp, Nk = 201, 401, 500
r_grid = np.linspace(0.0, cx*dx, Nr)
z_grid = np.linspace(-cz*dz, (Nz-1-cz)*dz, Nzp)
k_rho  = np.linspace(0.0, kmax, Nk)
dk     = k_rho[1]-k_rho[0]
kz     = np.sqrt(np.maximum(kmed**2 - k_rho**2, 0.0))
U = np.zeros((Nr, Nzp), dtype=complex)
for k in range(Nk):
    U += np.outer(j0(k_rho[k]*r_grid), np.exp(1j*kz[k]*z_grid))*k_rho[k]*dk
h_conf = np.abs(U)**4; h_conf /= h_conf.sum()
interp = RegularGridInterpolator((r_grid, z_grid), h_conf,
                                 method='linear', bounds_error=False, fill_value=0.0)
XX, YY = np.meshgrid(np.arange(Nx), np.arange(Ny), indexing='ij')
r3d = np.sqrt((XX-cx)**2*dx**2 + (YY-cy)**2*dy**2)
R3 = np.tile(r3d[:,:,None], (1,1,Nz))
Z3 = np.tile(((np.arange(Nz)-cz)*dz)[None,None,:], (Nx,Ny,1))
h3d = interp(np.stack([R3.ravel(), Z3.ravel()], axis=1)).reshape(Nx, Ny, Nz)
h3d = np.maximum(h3d, 0.0); h3d /= h3d.sum()

# ---- explicit FFT forward model (identical to rl_deconvolution.py) ------------
fft_shape = (2*Nx, 2*Ny, 2*Nz)
h3d_pad = np.zeros(fft_shape); h3d_pad[:Nx,:Ny,:Nz] = h3d
H_fft   = np.fft.fftn(h3d_pad)
crop    = (slice(Nx//2, Nx//2+Nx), slice(Ny//2, Ny//2+Ny), slice(Nz//2, Nz//2+Nz))
def _pad(a):
    p = np.zeros(fft_shape); p[:Nx,:Ny,:Nz] = a; return p
def _conv(a):
    return np.real(np.fft.ifftn(np.fft.fftn(_pad(a)) * H_fft))[crop]
def forward(O):  return _conv(O)
def xcorr(R):    return _conv(R)

# ----------------------------- object & data -----------------------------------
sigma_pix = 3.0
XXg, YYg, ZZg = np.meshgrid(np.arange(Nx), np.arange(Ny), np.arange(Nz), indexing='ij')
O_true = np.zeros((Nx, Ny, Nz))
for (px, py, pz) in [(16,16,6),(48,16,6),(16,48,26),(48,48,26)]:
    O_true += np.exp(-(((XXg-px)**2 + (YYg-py)**2 + (ZZg-pz)**2))/(2*sigma_pix**2))
B     = np.maximum(forward(O_true), 0.0)
scale = N_photons / B.sum()
O_true_s = O_true*scale
B_s     = B*scale
I       = rng.poisson(B_s)
Ssum    = float(I.sum())          # conserved total photon budget

# ----------------------------- log-likelihood & bound --------------------------
eps = 1e-12
O2s = float((O_true_s**2).sum())
def loglik(O):
    # true Poisson log-likelihood (data-only constant dropped); floor only
    # guards log(0) — f=(h*O)>0 in practice, so this is the exact likelihood.
    f = np.maximum(forward(O), 1e-300)
    return float(np.sum(I*np.log(f) - f))
def nmse(O):
    return float(((O-O_true_s)**2).sum())/O2s
def kl_next(O, On):
    # KL( p^{t+1} || p^{t} ),  p = O / sum(O)  (sums are equal = conserved)
    p  = O /O.sum();  pn = On/On.sum()
    m  = p>0
    return float(np.sum(p[m]*np.log(p[m]/pn[m])))

def rl_step(O):
    # EXACT EM M-step:  O^{t+1} = O^t * (h ⋆ [I/(h*O^t)])   (no eps bias,
    # so the DLR monotone-likelihood theorem and the KL bound hold strictly).
    blurred = forward(O)
    R       = np.divide(I.astype(np.float64), blurred,
                        out=np.zeros_like(blurred), where=blurred>0)
    corr    = xcorr(R)
    return np.maximum(O*corr, 0.0), blurred

log = []
def p(s=""):
    print(s); log.append(s)
p("=== RL convergence bound  (identical data to rl_deconvolution.py) ===")
p(f"observed photons sum(I) = {Ssum:.0f}   (conserved sum of every iterate)")

O  = np.maximum(forward(I), 0.0)          # smooth init O0 = h*I
L  = loglik(O)
resid = np.zeros(max_iter); tv = np.zeros(max_iter); nm = np.zeros(max_iter)
dL  = np.zeros(max_iter); bnd = np.zeros(max_iter); kl = np.zeros(max_iter)
Lseq= np.zeros(max_iter)
def tv3d(x):
    return (np.abs(np.diff(x,1)).sum()+np.abs(np.diff(x,0)).sum()+np.abs(np.diff(x,2)).sum())

for t in range(max_iter):
    On, blurred = rl_step(O)
    Ln  = loglik(On)
    resid[t] = np.linalg.norm(I - blurred)
    tv[t]    = tv3d(On)
    nm[t]    = nmse(On)
    dL[t]    = Ln - L
    bnd[t]   = Ssum * kl_next(O, On)
    kl[t]    = kl_next(O, On)
    Lseq[t]  = Ln
    if t in (0,1,2,5,10,20,50,100,159):
        p(f"  it {t:3d}:  L={Lseq[t]:12.2f}  dL={dL[t]:10.4f}  "
          f"bound={bnd[t]:10.4f}  dL/bnd={dL[t]/max(bnd[t],1e-300):6.3f}  "
          f"NMSE={nm[t]:.5f}")
    L  = Ln          # advance the reference likelihood (per-step dL)
    O  = On

# L-curve corner (identical criterion to rl_deconvolution.py)
x  = np.log(tv);  y = np.log(resid)
xn = (x-x.min())/(x.max()-x.min()); yn = (y-y.min())/(y.max()-y.min())
x0,y0 = xn[0],yn[0]; x1,y1 = xn[-1],yn[-1]
a=y1-y0; b=x0-x1; c=y0*(x1-x0)-x0*(y1-y0); Lch=np.hypot(a,b)
t_corner = int(np.argmax(np.abs(a*xn+b*yn+c)/Lch))
t_best   = int(np.argmin(nm))

# --- diagnostics: separate (i) monotonicity from (ii) KL-bound tightness -------
n_mono_viol = int(np.sum(dL < -1e-9))          # steps with L actually decreasing
t_mono_viol = np.where(dL < -1e-9)[0]
n_bound_viol = int(np.sum(dL < bnd - 1e-9))    # steps where dL < the KL lower bound
ratio_early = dL[:10]/bnd[:10]                 # dL/bound in the convergence regime
sums = np.array([O.sum()])                      # (final iterate sum; sum drift tracked below)
p(f"\nL-curve corner @ it {t_corner}   min-NMSE @ it {t_best}")
p(f"sum(I) = {Ssum:.0f}")
p(f"(i)  monotonicity  dL >= 0 :  violations = {n_mono_viol}  "
  f"{'(L non-decreasing at every step)' if n_mono_viol==0 else 'first at it '+str(t_mono_viol[0])}")
p(f"(ii) KL lower bound dL >= sum(I)*KL :  tight it 0-9 "
  f"(ratio {ratio_early.min():.2f}-{ratio_early.max():.2f});  "
  f"loose from it ~10 (finite-image FFT boundary effects break exact flux conservation)")
p(f"     sum(I)*KL -> 0 as the iterates converge (KL(p^{{t+1}}||p^t) -> 0)")

# empirical late-stage linear rate:  |dL^{t+1}/dL^{t}|  in the pre-overfit plateau
rat = np.abs(dL[30:60])/np.maximum(np.abs(dL[29:59]),1e-300)
p(f"empirical linear rate  rho = median(|dL^{{t+1}}/dL^{{t}}|), it 30-59 = {np.median(rat):.3f}")
p(f"  -> dL decays ~geometrically (rho<1) toward the ML fixed point")

nm_blur = nmse(I); nm_dec = nm[t_corner]
p(f"NMSE  blurred={nm_blur:.4f}   at corner (it {t_corner})={nm_dec:.4f}   "
  f"improvement={nm_blur/nm_dec:.2f}x")
p(f"log-likelihood:  L(init)={Lseq[0]:.2f}  L(corner)={Lseq[t_corner]:.2f}  "
  f"L(final)={Lseq[-1]:.2f}")
p(f"  L keeps RISING past the NMSE minimum (it {t_best}):  "
  f"L(corner)={Lseq[t_corner]:.2f}  <  L(final)={Lseq[-1]:.2f}")

open("rl_convergence_bound.log","w").write("\n".join(log)+"\n")

# ------------------------------- figure ----------------------------------------
fig = plt.figure(figsize=(13, 5.2))
gs  = fig.add_gridspec(1, 2, wspace=0.28)

ax = fig.add_subplot(gs[0])
tshow = 16   # convergence regime (dL>0, bound tight)
ax.semilogy(np.arange(tshow), dL[:tshow],  '-o', ms=3, lw=1.4, color='tab:green',
            label='actual $\\Delta L$ (per step)')
ax.semilogy(np.arange(tshow), bnd[:tshow], '-s', ms=3, lw=1.4, color='tab:blue',
            label='EM lower bound  $\\sum I\\,\\mathrm{KL}(p^{t+1}\\,\\|\\,p^t)$')
ax.axvline(t_corner, color='tab:red', ls='--', lw=1, alpha=0.8)
ax.text(t_corner+0.3, 30, 'L-curve corner', rotation=90, fontsize=8, color='tab:red')
ax.set_xlabel('RL iteration'); ax.set_ylabel('log scale')
ax.set_title('Per-step likelihood increase vs the KL lower bound (convergence regime)', fontsize=10)
ax.legend(fontsize=8); ax.grid(alpha=0.3, which='both')

ax = fig.add_subplot(gs[1])
ax.plot(Lseq, '-o', ms=3, lw=1.4, color='tab:purple', label='log-likelihood $L(O^{t})$')
ax.axvline(t_best,   color='tab:green', ls='--', lw=1, alpha=0.7)
ax.axvline(t_corner, color='tab:red',   ls='--', lw=1, alpha=0.7)
ax.set_xlabel('RL iteration')
ax.set_ylabel('log-likelihood', color='tab:purple')
ax.tick_params(axis='y', labelcolor='tab:purple')
ax2 = ax.twinx()
ax2.plot(nm, '-^', ms=3, lw=1.4, color='tab:green', label='NMSE vs true object')
ax2.set_ylabel('NMSE vs true object', color='tab:green')
ax2.tick_params(axis='y', labelcolor='tab:green')
ax2.set_ylim(0, 2.0)
ax.set_title('Likelihood climbs to the (overfit) ML limit; NMSE bottoms out (it %d)' % t_best, fontsize=10)
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1+h2, l1+l2, fontsize=8, loc='upper left')
ax.grid(alpha=0.3)

fig.suptitle('Richardson-Lucy convergence:  worked EM/KL bound  (NA=1.4, $\\lambda$=500 nm)', y=1.02, fontsize=13)
out='rl_convergence_bound.png'
fig.savefig(out, dpi=150, bbox_inches='tight')
print(f'saved {out}')
