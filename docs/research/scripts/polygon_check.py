"""Half the solid angle of a geodesic polygon on the unit sphere, summed from triangles.

Companion to holonomy_check.py: computes Omega/2 for the same circular caps by
Van Oosterom-Strackee triangle areas, as an independent check of the closed form
pi(1 - cos alpha). Needs numpy.
Run from the repo root: python3 docs/research/scripts/polygon_check.py
"""
import numpy as np
def tri_area(a,b,c):
    num=abs(np.dot(a,np.cross(b,c))); den=1+np.dot(a,b)+np.dot(b,c)+np.dot(c,a)
    return 2*np.arctan2(num,den)
for alpha in (0.3,1.0,2.0):
  for n in (12,100):
    t=np.linspace(0,2*np.pi,n+1); V=np.stack([np.sin(alpha)*np.cos(t),np.sin(alpha)*np.sin(t),np.cos(alpha)*np.ones_like(t)],1)
    P=np.array([0,0,1.0]); A=sum(tri_area(P,V[k],V[k+1]) for k in range(n))
    print(f"alpha={alpha} N={n}: half solid angle of geodesic polygon = {((A/2)+np.pi)%(2*np.pi)-np.pi:+.5f}")
