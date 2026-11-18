#!/usr/bin/env python
"""
Debye 3D PSF expansion (companion to Confocal_Imaging_Essay v2, section 4b).

Computes, on the EXACT grid used by rl_deconvolution.py (so the numbers are
consistent with the essay):
  * the full 3D scalar Debye integral  U(r,z) = ∫_0^{kmax} J0(kρ r) e^{i kz z} kρ dkρ
  * its on-axis limit  U(0,z)  and the four on-axis obliquity conventions
        paraxial (Abbe)      weight u^0
        Debye scalar 3-D     weight u^1    (= on-axis limit of the 3-D integral)
        Richards-Wolf scalar weight u^{1/2}
        Richards-Wolf vector weight u^{1/2} + u^{5/2}
  * the (lateral, axial) FWHM of the confocal |U|^4 and single-path |U|^2 PSF
  * the k-space "missing cone" support

Writes debye_psf_3d.png and prints a log to stdout.
"""
import numpy as np
from scipy.special import j0
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ------------------- parameters (identical to rl_deconvolution.py) -------------
lam  = 500e-9
n0   = 1.5
NA   = 1.4
k0   = 2.0*np.pi/lam
kmed = n0*k0
kmax = 2.0*np.pi*NA/lam
cos_alpha = float(np.sqrt(max(1.0 - (kmax/kmed)**2, 0.0)))
alpha_deg = np.degrees(np.arccos(cos_alpha))

Nx = Ny = 64; Nz = 32
dx = dy = 10e-9; dz = 20e-9
cx, cy, cz = Nx//2, Ny//2, Nz//2

# 3-D PSF (r,z) grid — identical to rl_deconvolution.py
Nr, Nzp, Nk = 201, 401, 500
r_grid = np.linspace(0.0, cx*dx, Nr)
z_grid = np.linspace(-cz*dz, (Nz-1-cz)*dz, Nzp)
k_rho  = np.linspace(0.0, kmax, Nk)
dk     = k_rho[1]-k_rho[0]
kz     = np.sqrt(np.maximum(kmed**2 - k_rho**2, 0.0))

U = np.zeros((Nr, Nzp), dtype=complex)
for k in range(Nk):
    U += np.outer(j0(k_rho[k]*r_grid), np.exp(1j*kz[k]*z_grid))*k_rho[k]*dk

h4 = np.abs(U)**4          # confocal (illumination x detection)
h2 = np.abs(U)**2          # single incoherent path
i0z = int(np.argmin(np.abs(z_grid)))     # index of z=0
i0r = 0                                # r=0 column

def fwhm(z, prof):
    """FWHM via the FIRST half-max crossings out from the peak (sidelobe-safe)."""
    i = int(np.argmax(prof)); half = 0.5*prof[i]
    il = i
    while il > 0 and prof[il] >= half: il -= 1
    ir = i
    while ir < len(prof)-1 and prof[ir] >= half: ir += 1
    zl = np.interp(half, [prof[il],   prof[il+1]], [z[il],   z[il+1]])
    zr = np.interp(half, [prof[ir],   prof[ir-1]], [z[ir],   z[ir-1]])
    return (zr-zl), zl, zr

log = []
def p(s=""):
    print(s); log.append(s)

p("=== Debye 3-D confocal PSF  (NA=%.2f n0=%.1f lam=%.0f nm) ===" % (NA, n0, lam*1e9))
p(f"cos(alpha)={cos_alpha:.4f}  alpha={alpha_deg:.2f} deg  kmax/kmed={kmax/kmed:.4f}")

# ---- on-axis: the four obliquity conventions, same z_grid & beta -------------
beta = kmed*z_grid
u = np.linspace(cos_alpha, 1.0, 4001)
du = u[1]-u[0]
def trap(w, e):
    # w: (Nu,) weight along u;  e: (Nu, Nz) = exp(i beta u);  integrate over u
    return np.trapezoid(w[:,None]*e, u, axis=0)
U_par  = trap(np.ones_like(u),        np.exp(1j*beta[None,:]*u[:,None]))
U_deb  = trap(u,                     np.exp(1j*beta[None,:]*u[:,None]))
U_rws  = trap(np.sqrt(u),            np.exp(1j*beta[None,:]*u[:,None]))
U_rvw  = trap(np.sqrt(u)+u**2.5,     np.exp(1j*beta[None,:]*u[:,None]))
# on-axis limit of the actual 3-D integral (r=0 column) — should match U_deb
U_deb3d = U[i0r, :]

p("\n--- on-axis |U|^2 FWHM (single path), nm (full width) ---")
for name, Uv in [("on-axis, wt u^0 (obl-free)", U_par),
                 ("Debye 3-D         u^1      ", U_deb),
                 ("Debye 3-D (grid)  u^1      ", U_deb3d),
                 ("R-W scalar        u^{1/2}  ", U_rws),
                 ("R-W vectorial     u^{1/2}+u^{5/2}", U_rvw)]:
    w,_,_ = fwhm(z_grid, np.abs(Uv)**2)
    p(f"  {name}:  {w*1e9:7.1f} nm")
p("  [Part (a) Abbe paraxial reference = 678 nm (small-angle axial bandwidth)]")

p("\n--- axial FWHM at r=0, nm ---")
w4,_,_ = fwhm(z_grid, h4[i0r,:]); w2,_,_ = fwhm(z_grid, h2[i0r,:])
p(f"  confocal |U|^4 :  {w4*1e9:7.1f} nm")
p(f"  single   |U|^2 :  {w2*1e9:7.1f} nm")
p(f"  ratio |U|^4 / |U|^2 = {w4/w2:.3f}  (confocal narrows the axial lobe)")

p("\n--- lateral FWHM at z=0, nm (full width = 2 x radius) ---")
wl4,_,_ = fwhm(r_grid, h4[:,i0z]); wl2,_,_ = fwhm(r_grid, h2[:,i0z])
p(f"  confocal |U|^4 :  {2*wl4*1e9:7.1f} nm")
p(f"  single   |U|^2 :  {2*wl2*1e9:7.1f} nm")
# paraxial Airy reference: first zero of 2 J1(kmax r)/(kmax r)
r_airy = np.linspace(0, cx*dx, 20000)
airy = (2*j0(0)/r_airy) if False else None
from scipy.special import j1
airy = np.where(r_airy>0, (2*j1(kmax*r_airy)/(kmax*r_airy))**2, 1.0)
# Airy FWHM
ia = int(np.argmax(airy)); half=0.5*airy[ia]
ira=ia
while ira<len(airy)-1 and airy[ira]>=half: ira+=1
r_airy_fwhm = np.interp(half,[airy[ira],airy[ira-1]],[r_airy[ira],r_airy[ira-1]])
p(f"  paraxial Airy  :  {2*r_airy_fwhm*1e9:7.1f} nm  (first-zero {1.22*lam/NA*1e9:.1f} nm)")

p("\n--- peak & sidelobe (confocal |U|^4) ---")
ax_prof = h4[i0r,:]/h4[i0r,:].max()
lat_prof = h4[:,i0z]/h4[:,i0z].max()
p(f"  axial  main-lobe peak=1.000  first side lobe={ax_prof[1:25].max():.4f}  "
  f"at z={z_grid[1+np.argmax(ax_prof[1:25])]*1e9:.0f} nm")
p(f"  lateral main-lobe peak=1.000 first side lobe={lat_prof[1:25].max():.4f}  "
  f"at r={r_grid[1+np.argmax(lat_prof[1:25])]*1e9:.0f} nm")

p("\n--- k-space support (missing cone) ---")
p(f"  propagating support:  k_rho <= kmax = {kmax:.3e} 1/m  (cone half-angle alpha)")
p(f"  single-path OTF support radius (transverse)  = 2 kmax = {2*kmax:.3e} 1/m")
p(f"  confocal OTF support radius (transverse)     = 2 kmax = {2*kmax:.3e} 1/m  "
  f"(product of two pupils -> doubled Abbe bandwidth)")
p(f"  missing cone:  k_rho > kmax  (evanescent/super-NA)  ->  "
  f"angular gap = {90-alpha_deg:.2f} deg from the axis in k-space")

open("debye_psf_3d.log","w").write("\n".join(log)+"\n")

# ------------------------------- figure ----------------------------------------
fig = plt.figure(figsize=(13, 9.5))
gs  = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.22)

# (0,0) 2-D filled contour of the confocal PSF (log), (r, z)
ax = fig.add_subplot(gs[0,0])
RR, ZZ = np.meshgrid(r_grid*1e9, z_grid*1e9, indexing='ij')
surf = np.log10(h4/h4.max()+1e-12)
cf = ax.contourf(RR, ZZ, surf, levels=22, cmap='inferno')
cs = ax.contour(RR, ZZ, surf, levels=[-12,-10,-8,-6,-4,-3,-2,-1],
                colors='w', linewidths=0.5, alpha=0.6)
ax.clabel(cs, fmt='%d', fontsize=6, colors='w')
fig.colorbar(cf, ax=ax, label=r'$\log_{10}|U|^4$')
ax.set_xlabel('r (nm)'); ax.set_ylabel('z (nm)')
ax.set_title('Confocal Debye PSF  $h=|U|^4$  (prolate lobe, elongated along $z$)', fontsize=9)
ax.set_aspect('auto')

# (0,1) on-axis |U|^2, four obliquity conventions
ax = fig.add_subplot(gs[0,1])
zz = z_grid*1e9
for name, Uv, c in [("obliquity-free $u^0$", U_par, 'tab:blue'),
                    ("Debye 3-D $u^1$", U_deb, 'tab:red'),
                    ("R-W scalar $u^{1/2}$", U_rws, 'tab:green'),
                    ("R-W vectorial", U_rvw, 'tab:purple')]:
    pr = np.abs(Uv)**2; pr/=pr.max()
    ax.plot(zz, pr, lw=1.6, color=c, label=name)
ax.set_xlabel('z (nm)'); ax.set_ylabel('$|U(0,z)|^2$ (norm.)')
ax.set_title('On-axis lobe: obliquity weight shifts FWHM only ~5% (461$\\to$482 nm)', fontsize=9)
ax.set_xlim(-250,250); ax.legend(fontsize=8); ax.grid(alpha=0.3)

# (1,0) confocal vs single-path axial
ax = fig.add_subplot(gs[1,0])
a4 = h4[i0r,:]; a4/=a4.max(); a2 = h2[i0r,:]; a2/=a2.max()
ax.plot(zz, a4, lw=1.8, color='tab:red',   label='confocal $|U|^4$')
ax.plot(zz, a2, lw=1.8, color='tab:orange',label='single path $|U|^2$')
ax.axvline(w4*1e9/2,  color='tab:red',   ls=':', lw=1)
ax.axvline(-w4*1e9/2, color='tab:red',   ls=':', lw=1)
ax.set_xlabel('z (nm)'); ax.set_ylabel('axial PSF (norm.)')
ax.set_title(f'Axial: confocal narrows  (FWHM {w2*1e9:.0f}$\\to${w4*1e9:.0f} nm)', fontsize=9)
ax.set_xlim(-300,300); ax.legend(fontsize=8); ax.grid(alpha=0.3)

# (1,1) lateral at z=0
ax = fig.add_subplot(gs[1,1])
l4 = h4[:,i0z]; l4/=l4.max(); l2 = h2[:,i0z]; l2/=l2.max()
rr = r_grid*1e9
ax.plot(rr, l4, lw=1.8, color='tab:red',    label='confocal $|U|^4$')
ax.plot(rr, l2, lw=1.8, color='tab:orange', label='single path $|U|^2$')
ax.plot(rr, np.interp(rr, r_airy*1e9, airy/airy.max()), lw=1.4, color='tab:blue',
        ls='--', label='paraxial Airy')
ax.set_xlabel('r (nm)'); ax.set_ylabel('lateral PSF (norm.)')
ax.set_title(f'Lateral: confocal $|U|^4$ (FWHM {2*wl2*1e9:.0f}$\\to${2*wl4*1e9:.0f} nm)', fontsize=9)
ax.legend(fontsize=8); ax.grid(alpha=0.3)

fig.suptitle('Full 3-D Debye confocal PSF  (NA=1.4, $\\lambda$=500 nm, n$_0$=1.5)', y=0.995, fontsize=13)
out='debye_psf_3d.png'
fig.savefig(out, dpi=150, bbox_inches='tight')
print(f'saved {out}')
