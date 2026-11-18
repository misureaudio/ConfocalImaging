#!/usr/bin/env python
"""Debug: verify the EM monotone-likelihood theorem step by step for the
Poisson RL deconvolution, and pinpoint any non-monotonicity."""
import numpy as np
from scipy.special import j0
from scipy.interpolate import RegularGridInterpolator

rng = np.random.default_rng(42)
lam=500e-9; n0=1.5; NA=1.4
k0=2*np.pi/lam; kmed=n0*k0; kmax=2*np.pi*NA/lam
Nx=Ny=64; Nz=32; dx=dy=10e-9; dz=20e-9
cx,cy,cz=Nx//2,Ny//2,Nz//2; N_photons=200000
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
def forward(O): return _conv(O)
def xcorr(R): return _conv(R)
sigma_pix=3.0
XXg,YYg,ZZg=np.meshgrid(np.arange(Nx),np.arange(Ny),np.arange(Nz),indexing='ij')
O_true=np.zeros((Nx,Ny,Nz))
for (px,py,pz) in [(16,16,6),(48,16,6),(16,48,26),(48,48,26)]:
    O_true+=np.exp(-(((XXg-px)**2+(YYg-py)**2+(ZZg-pz)**2))/(2*sigma_pix**2))
B=np.maximum(forward(O_true),0.0); scale=N_photons/B.sum()
O_true_s=O_true*scale; B_s=B*scale; I=rng.poisson(B_s)
Ssum=float(I.sum())

def loglik(O, floor):
    f=np.maximum(forward(O),floor)
    return float(np.sum(I*np.log(f)-f))

def rl_step(O):
    blurred=forward(O)
    R=np.divide(I.astype(np.float64),blurred,out=np.zeros_like(blurred),where=blurred>0)
    corr=xcorr(R)
    return np.maximum(O*corr,0.0), blurred, R

print("floor      L(O0)        L(O1)        dL0       monotone?")
for floor in [0.0, 1e-300, 1e-12, 1e-6]:
    O=np.maximum(forward(I),0.0)
    L0=loglik(O,floor)
    On,blurred,R=rl_step(O)
    L1=loglik(On,floor)
    print(f"{floor:<9} {L0:12.4f} {L1:12.4f} {L1-L0:10.4f}   {L1>=L0}")

# Now: run 30 iterations with floor=0 (no floor, but guard log(0) via where)
# and check the EM theorem  dL >= Q-gap  at each step.
print("\n--- EM theorem check (no floor, log guarded) ---")
print("  it   L(O^t)       dL        Q-gap     dL>=Qgap   sumO^t  sumO^t+1")
O=np.maximum(forward(I),0.0)
def loglik_g(O):
    f=forward(O)
    # log only where f>0 and I>0; where I==0 the I*log term is 0
    term=np.zeros_like(f)
    m=(I>0)&(f>0)
    term[m]=I[m]*np.log(f[m])
    return float(np.sum(term-f))
L=loglik_g(O)
prev_sum=O.sum()
for t in range(30):
    On,blurred,R=rl_step(O)
    a=xcorr(R)                       # correction factor (h ⋆ R)
    # Q-gap = sum_j O^t_j [a_j ln a_j + 1 - a_j]
    am=np.maximum(a,1e-300)
    qgap=float(np.sum(O*(a*np.log(am)+1-a)))
    Ln=loglik_g(On)
    dL=Ln-L
    print(f"  {t:3d}  {Ln:11.4f}  {dL:9.5f}  {qgap:9.5f}   {dL>=qgap-1e-9}   {prev_sum:.1f}  {On.sum():.1f}")
    L=Ln; O=On; prev_sum=O.sum()
