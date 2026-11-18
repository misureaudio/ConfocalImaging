#!/usr/bin/env python
"""Clean 2D verification of the OTF cutoffs via the pupil autocorrelation.

Widefield incoherent OTF = autocorrelation of the pupil (a disk of radius kmax)
  -> support |k| <= 2 kmax.
Confocal OTF = convolution of two widefield OTFs
  -> support |k| <= 4 kmax.
"""
import numpy as np

lam = 500e-9; NA = 1.4
kmax = 2*np.pi*NA/lam

# 2D pupil: disk of radius kmax
N = 1024
kspan = 8*kmax                      # grid half-width in k-space (well beyond 4kmax)
dk = 2*kspan/N
K = (np.arange(N)-N//2)*dk
KX, KY = np.meshgrid(K, K, indexing='ij')
KR = np.sqrt(KX**2+KY**2)
P = (KR <= kmax).astype(float)      # pupil: 1 inside the NA disk

def otf_support(otf, thresh=1e-6):
    o = np.abs(otf)
    o /= o.max()
    r = KR[o > thresh]
    return r.max() if r.size else 0.0

# widefield OTF = autocorrelation of the pupil
P_fft = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(P)))
otf_wf = np.abs(np.fft.ifft2(P_fft * np.conj(P_fft)))
otf_wf = np.fft.fftshift(otf_wf)
sup_wf = otf_support(otf_wf)

# confocal OTF = convolution of two widefield OTFs (in k-space)
#   (f * g) = ifft( fft(f) . fft(g) )
fft_wf = np.fft.fft2(np.fft.ifftshift(otf_wf))
otf_conf = np.fft.fftshift(np.abs(np.fft.ifft2(fft_wf * fft_wf)))
sup_conf = otf_support(otf_conf)

print(f"kmax   = {kmax:.4e} 1/m")
print(f"2 kmax = {2*kmax:.4e} 1/m  (widefield OTF cutoff)")
print(f"4 kmax = {4*kmax:.4e} 1/m  (confocal OTF cutoff)")
print()
print(f"widefield OTF support edge  = {sup_wf:.4e} 1/m   ratio/(2 kmax) = {sup_wf/(2*kmax):.3f}")
print(f"confocal  OTF support edge  = {sup_conf:.4e} 1/m   ratio/(4 kmax) = {sup_conf/(4*kmax):.3f}")
print(f"confocal/widefield support  = {sup_conf/sup_wf:.3f}   (expect 2.0)")
