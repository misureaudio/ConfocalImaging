#!/usr/bin/env python
"""Check whether the forward operator in rl_convergence_bound.py is self-adjoint.

The script uses  xcorr(R) = _conv(R) = forward(R)  (the SAME crop-and-pad
convolution for the forward and the backward step).  The exact EM M-step
requires the backward step to be the matrix TRANSPOSE of the forward operator.
If forward is not self-adjoint, using it as its own adjoint breaks the exact
M-step condition, and the flux is not exactly conserved -> the observed drift.
"""
import numpy as np
from scipy.special import j0
from scipy.interpolate import RegularGridInterpolator

rng = np.random.default_rng(0)
lam=500e-9; n0=1.5; NA=1.4
k0=2*np.pi/lam; kmed=n0*k0; kmax=2*np.pi*NA/lam
Nx=Ny=64; Nz=32; dx=dy=10e-9; dz=20e-9
cx,cy,cz=Nx//2,Ny//2,Nz//2
Nr,Nzp,Nk=201,401,500
r_grid=np.linspace(0.0,cx*dx,Nr); z_grid=np.linspace(-cz*dz,(Nz-1-cz)*dz,Nzp)
k_rho=np.linspace(0.0,kmax,Nk); dk=k_rho[1]-k_rho[0]
kz=np.sqrt(np.maximum(kmed**2-k_rho**2,0.0))
U=np.zeros((Nr,Nzp),dtype=complex)
for k in range(Nk):
    U+=np.outer(j0(k_rho[k]*r_grid),np.exp(1j*kz[k]*z_grid))*k_rho[k]*dk
h_conf=np.abs(U)**4; h_conf/=h_conf.sum()
interp=RegularGridInterpolator((r_grid,z_grid),h_conf,method='linear',bounds_error=False,fill_value=0.0)
XX,YY=np.meshgrid(np.arange(Nx),np.arange(Ny),indexing='ij')
r3d=np.sqrt((XX-cx)**2*dx**2+(YY-cy)**2*dy**2)
R3=np.tile(r3d[:,:,None],(1,1,Nz)); Z3=np.tile(((np.arange(Nz)-cz)*dz)[None,None,:],(Nx,Ny,1))
h3d=interp(np.stack([R3.ravel(),Z3.ravel()],axis=1)).reshape(Nx,Ny,Nz)
h3d=np.maximum(h3d,0.0); h3d/=h3d.sum()
fft_shape=(2*Nx,2*Ny,2*Nz); h3d_pad=np.zeros(fft_shape); h3d_pad[:Nx,:Ny,:Nz]=h3d
H_fft=np.fft.fftn(h3d_pad)
crop=(slice(Nx//2,Nx//2+Nx),slice(Ny//2,Ny//2+Ny),slice(Nz//2,Nz//2+Nz))
def _pad(a):
    p=np.zeros(fft_shape); p[:Nx,:Ny,:Nz]=a; return p
def _conv(a): return np.real(np.fft.ifftn(np.fft.fftn(_pad(a))*H_fft))[crop]
def forward(a): return _conv(a)

# adjointness test: <forward(a), b>  vs  <a, forward(b)>
a = rng.random((Nx,Ny,Nz))
b = rng.random((Nx,Ny,Nz))
lhs = np.sum(forward(a)*b)
rhs = np.sum(a*forward(b))
print(f"<forward(a), b> = {lhs:.6f}")
print(f"<a, forward(b)> = {rhs:.6f}")
print(f"relative difference = {abs(lhs-rhs)/abs(lhs):.3e}")
print(f"forward is self-adjoint? {abs(lhs-rhs)/abs(lhs) < 1e-9}")

# flux conservation of the RL update (exact M-step, eps=0)
XXg,YYg,ZZg=np.meshgrid(np.arange(Nx),np.arange(Ny),np.arange(Nz),indexing='ij')
O_true=np.zeros((Nx,Ny,Nz))
for (px,py,pz) in [(16,16,6),(48,16,6),(16,48,26),(48,48,26)]:
    O_true+=np.exp(-(((XXg-px)**2+(YYg-py)**2+(ZZg-pz)**2))/(2*3.0**2))
B=np.maximum(forward(O_true),0.0); scale=200000.0/B.sum()
I=rng.poisson(B*scale)
def rl_step(O):
    blurred=forward(O)
    R=np.divide(I.astype(float),blurred,out=np.zeros_like(blurred),where=blurred>0)
    return np.maximum(O*forward(R),0.0), blurred
O=np.maximum(forward(I),0.0)
print(f"\nsum(I) = {I.sum()}")
for t in range(6):
    s0=O.sum()
    O,blurred=rl_step(O)
    print(f"  it {t}: sumO {s0:.1f} -> {O.sum():.1f}   (drift {O.sum()-s0:+.1f})")
