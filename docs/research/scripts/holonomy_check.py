"""Chained Procrustes holonomy on the triangle shape sphere.

Carries a triangle round a closed loop of shapes, fitting each step to the previous one by
the best rotation (mass-weighted). The net rotation after the loop should equal half the
solid angle Omega/2 enclosed on the shape sphere, for any masses and any loop, and the
one-shot fit of last onto first should be zero. Needs numpy.
Run from the repo root: python3 docs/research/scripts/holonomy_check.py
"""
import numpy as np
def tri_from_w(w, m):
    m1,m2,m3=m; M12=m1+m2; M=M12+m3
    mu1=np.sqrt(m1*m2/M12); mu2=np.sqrt(m3*M12/M)
    w=w/np.linalg.norm(w); w1,w2,w3=w
    Z1=np.sqrt((1+w1)/2); Z2=(w2-1j*w3)/(2*Z1)
    c12=-m3*Z2/(mu2*M); q3=Z2/mu2*M12/M
    q1=c12+(m2/M12)*Z1/mu1; q2=c12-(m1/M12)*Z1/mu1
    return np.array([q1,q2,q3])
def chained(loop, m):
    m=np.array(m,float); ref=tri_from_w(loop[0],m); q0=ref.copy()
    for w in loop[1:]:
        q=tri_from_w(w,m); c=np.sum(m*q*np.conj(ref)); ref=np.exp(-1j*np.angle(c))*q
    c=np.sum(m*ref*np.conj(q0)); return np.angle(c)   # net rotation of final vs initial
def one_shot(loop,m):
    m=np.array(m,float); a=tri_from_w(loop[0],m); b=tri_from_w(loop[-1],m)
    c=np.sum(m*b*np.conj(a)); r=np.exp(-1j*np.angle(c))*b; return np.angle(np.sum(m*r*np.conj(a)))
def cap(alpha,n,center='pole'):
    t=np.linspace(0,2*np.pi,n+1)
    if center=='pole':  # circle around w3 axis (equilateral pole for equal masses)
        return np.stack([np.sin(alpha)*np.cos(t),np.sin(alpha)*np.sin(t),np.cos(alpha)*np.ones_like(t)],1)
    # circle around equatorial point (0,1,0), i.e. a collinear shape
    return np.stack([np.sin(alpha)*np.cos(t),np.cos(alpha)*np.ones_like(t),np.sin(alpha)*np.sin(t)],1)
wrap=lambda x:(x+np.pi)%(2*np.pi)-np.pi
for m in [(1,1,1),(1,2,5)]:
  for center in ['pole','equator']:
    for alpha in [0.3,1.0,2.0]:
        pred=np.pi*(1-np.cos(alpha))   # half the solid angle of the cap
        res=[chained(cap(alpha,n,center),m) for n in (12,100,2000)]
        print(f"m={m} loop={center:7s} alpha={alpha:.1f} Omega/2={wrap(pred):+.5f} chained(N=12,100,2000)=" + ", ".join(f"{r:+.5f}" for r in res) + f"  one-shot={one_shot(cap(alpha,2000,center),m):+.1e}")
