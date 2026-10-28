"""Numerical check of the Procrustes identities in procrustes_history_mathematics.md, section 4.

(a) The complex full fit w*z/(w*w) is no worse than a brute-force search over 4,001 turns at
    the fitted scale, and E/|z|^2 = d_F^2.
(b) For the real 2x2 cross-product A, |w*z| = s1 + sign(det A) s2, |w^T z| = s1 - sign(det A) s2,
    and det A = (|w*z|^2 - |w^T z|^2)/4.
(c) d_P = 2 sin(rho/2) and d_F = sin(rho), with rho the Procrustes (Riemannian) distance.
(d) The shape-sphere identity cos^2 rho = (1 + n.n')/2, and a reflection flips the height.
Here `rho` is Dryden's distance, not the scale field of js/procrustes.js (now `beta`).
2,000 random configurations per identity. Needs numpy. Originally scratchpad/scripts/check.py.
Run from the repo root: python3 docs/research/scripts/identities_check.py
"""
import numpy as np
rng = np.random.default_rng(1)
def c(v): return v - v.mean()
ok = True
for trial in range(2000):
    k = rng.integers(3, 8)
    z = c(rng.normal(size=k) + 1j*rng.normal(size=k)); w = c(rng.normal(size=k) + 1j*rng.normal(size=k))
    # (a) complex full fit vs brute force over theta, beta
    a = np.vdot(w, z) / np.vdot(w, w).real          # w*z / w*w
    th = np.linspace(-np.pi, np.pi, 4001); best = min(np.linalg.norm(z - (abs(a))*np.exp(1j*t)*w)**2 for t in th)
    E = np.linalg.norm(z - a*w)**2
    ok &= E <= best + 1e-6
    dF2 = 1 - abs(np.vdot(w, z))**2 / (np.vdot(w, w).real*np.vdot(z, z).real)
    ok &= abs(E/np.vdot(z,z).real - dF2) < 1e-9
    # (b) real 2x2 SVD vs complex: A = Zt^T W (rows = landmarks)
    Z = np.c_[z.real, z.imag]; W = np.c_[w.real, w.imag]
    A = Z.T @ W; s = np.linalg.svd(A, compute_uv=False)
    p, q = abs(np.vdot(w, z)), abs(np.sum(w*z))   # |w*z|, |w^T z|
    sg = np.sign(np.linalg.det(A))
    ok &= abs(p - (s[0] + sg*s[1])) < 1e-9 and abs(q - (s[0] - sg*s[1])) < 1e-9
    ok &= abs(np.linalg.det(A) - (p**2 - q**2)/4) < 1e-9
    # partial distance and rho relations
    cr = p/np.sqrt(np.vdot(w,w).real*np.vdot(z,z).real); rho = np.arccos(min(1,cr))
    dP = np.linalg.norm(z/np.linalg.norm(z) - np.exp(1j*np.angle(np.vdot(w,z)))*w/np.linalg.norm(w))
    ok &= abs(dP - 2*np.sin(rho/2)) < 1e-9 and abs(np.sqrt(dF2) - np.sin(rho)) < 1e-9
# (c) triangle shape sphere identity (code's Helmert + Hopf map)
def sp(t):
    w1 = (t[1]-t[0])/np.sqrt(2); w2 = (2*t[2]-t[0]-t[1])/np.sqrt(6); n = abs(w1)**2+abs(w2)**2; cc = np.conj(w1)*w2
    return np.array([(abs(w1)**2-abs(w2)**2)/n, 2*cc.real/n, 2*cc.imag/n])
for trial in range(2000):
    A = rng.normal(size=3)+1j*rng.normal(size=3); B = rng.normal(size=3)+1j*rng.normal(size=3)
    z, u = c(A), c(B); cos2 = abs(np.vdot(z, u))**2/(np.vdot(z,z).real*np.vdot(u,u).real)
    ok &= abs(cos2 - (1+sp(A)@sp(B))/2) < 1e-9
    ok &= abs(sp(np.conj(A))[2] + sp(A)[2]) < 1e-12 and abs(sp(np.conj(A))[0]-sp(A)[0])<1e-12  # reflection flips height
eq = np.array([0, 1, np.exp(1j*np.pi/3)]); col = np.array([0, 1, 2.5])
print("all identities hold:", bool(ok)); print("equilateral ->", np.round(sp(eq),6), " collinear ->", np.round(sp(col),6))
