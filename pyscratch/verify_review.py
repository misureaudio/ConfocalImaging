#!/usr/bin/env python
"""Verify the review's numerical claims (errors 1, 2, 3)."""
import numpy as np
from scipy.special import j1, j0
from scipy.optimize import brentq

lam = 500e-9; NA = 1.4; n0 = 1.5
k0 = 2*np.pi/lam; kmed = n0*k0; kmax = 2*np.pi*NA/lam
alpha = np.arcsin(NA/n0)

print("=== CHECK 1: Airy FWHM vs §4 'lateral FWHM' (91.9 / 88.4) ===")
airy_fwhm = (1.6163/np.pi) * lam/NA
print(f"Airy FWHM (formula 1.6163/pi * lam/NA) = {airy_fwhm*1e9:.2f} nm")
print(f"  91.9 * 2 = {91.9*2:.1f} nm   <- §4 value doubled = true FWHM")
print(f"  88.4 * 2 = {88.4*2:.1f} nm   <- §4 scalar value doubled")

print("\n=== CHECK 2: confocal Airy FWHM root  [2J1/rho]^4 = 1/2 ===")
def f(rho): return (2*j1(rho)/rho)**4 - 0.5
rho_conf = brentq(f, 1.0, 1.5)
print(f"root rho_conf = {rho_conf:.4f}")
print(f"  at rho=1.140 (the §3 claim): [2J1/rho]^4 = {(2*j1(1.140)/1.140)**4:.4f}  (NOT 0.5)")
print(f"  ratio rho_conf/1.6163 = {rho_conf/1.6163:.4f}")
print(f"  §3 claim:  1.140/1.616  = {1.140/1.616:.4f}")
print(f"  review:    1.160/1.6163 = {1.160/1.6163:.4f}")

print("\n=== standard Airy FWHM root  [2J1/rho]^2 = 1/2 (sanity) ===")
def g(rho): return (2*j1(rho)/rho)**2 - 0.5
rho_airy = brentq(g, 1.0, 2.0)
print(f"root rho_airy = {rho_airy:.4f}  (should be ~1.6163)")

print("\n=== §4.1 numerical ratio (131.8/183.6) ===")
print(f"  131.8/183.6 = {131.8/183.6:.4f}   <- matches rho_conf/1.6163 = {rho_conf/1.6163:.4f}")

print("\n=== CHECK 3: confocal OTF transverse cutoff (3D FT of |U|^4) ===")
# Debye field U(r,z) on a cylindrical grid
Nr, Nzp, Nk = 241, 241, 400
rmax = 600e-9; zmax = 600e-9
r_grid = np.linspace(0.0, rmax, Nr)
z_grid = np.linspace(-zmax, zmax, Nzp)
k_rho = np.linspace(0.0, kmax, Nk); dk = k_rho[1]-k_rho[0]
kz = np.sqrt(np.maximum(kmed**2 - k_rho**2, 0.0))
U = np.zeros((Nr, Nzp), dtype=complex)
for k in range(Nk):
    U += np.outer(j0(k_rho[k]*r_grid), np.exp(1j*kz[k]*z_grid)) * k_rho[k]*dk
h_wf = np.abs(U)**2
h_conf = np.abs(U)**4

# map to a Cartesian grid and take the 3D FT; inspect the k_z=0 transverse slice
N = 192
dx = 2*rmax/N
X = (np.arange(N)-N//2)*dx
xx, yy = np.meshgrid(X, X, indexing='ij')
rr = np.sqrt(xx**2+yy**2)
from scipy.interpolate import RegularGridInterpolator
interp = RegularGridInterpolator((r_grid, z_grid), h_conf, method='linear',
                                 bounds_error=False, fill_value=0.0)
ZZ = np.tile(z_grid[None, None, :], (N, N, 1))   # (N, N, Nzp)
RR = np.tile(rr[:, :, None], (1, 1, Nzp))        # (N, N, Nzp)
pts = np.stack([RR.ravel(), ZZ.ravel()], axis=1) # (N*N*Nzp, 2)
cart = interp(pts).reshape(N, N, Nzp)
# zero-pad to 2N to avoid wrap-around aliasing of the support
P = np.zeros((2*N, 2*N, 2*Nzp))
P[:N, :N, :Nzp] = cart
H = np.fft.fftn(P)
H = np.fft.fftshift(H, axes=(0,1,2))
dkx = 1.0/(2*N*dx)
kx = (np.arange(2*N)-N)*dkx
# k_z = 0 plane
kz0 = N  # index of k_z=0 after fftshift of 2*Nzp -> index Nzp
slice_kz0 = np.abs(H[:, :, Nzp])
KX, KY = np.meshgrid(kx, kx, indexing='ij')
KR = np.sqrt(KX**2 + KY**2)
# find the transverse radius where the OTF falls to 1e-3 of peak (support edge)
prof = slice_kz0.max(axis=1)  # max over y for each x  -> crude; use radial profile instead
# radial profile of the k_z=0 slice
bins = np.linspace(0, kmax*1.2, 400)
rads = KR.ravel(); vals = slice_kz0.ravel()
import numpy as np
order = np.argsort(rads)
rads_s = rads[order]; vals_s = vals[order]
# bin-average
idx = np.digitize(rads_s, bins)
prof_r = np.array([vals_s[(idx==i)].mean() if np.any(idx==i) else 0 for i in range(1,len(bins))])
rc = 0.5*(bins[1:]+bins[:-1])
pn = prof_r/prof_r.max()
# support edge: last radius where pn > 1e-3
edge = rc[pn > 1e-3]
print(f"  kmax = {kmax:.4e} 1/m")
print(f"  2 kmax = {2*kmax:.4e} 1/m  (widefield OTF cutoff)")
print(f"  4 kmax = {4*kmax:.4e} 1/m  (confocal OTF cutoff, predicted)")
print(f"  confocal OTF transverse support edge (k_z=0 slice) = {edge.max():.4e} 1/m")
print(f"  ratio edge/(4 kmax) = {edge.max()/(4*kmax):.3f}   (expect ~1.0 if cutoff is 4kmax)")
print(f"  ratio edge/(2 kmax) = {edge.max()/(2*kmax):.3f}   (would be ~2 if cutoff were 2kmax)")
