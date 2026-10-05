# Search for an actual closed unit square inside [0,s]^2 capturing weight < 1 (exact check at the end).
import sys, numpy as np
from fractions import Fraction as F
cert, N, k = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
cx0, cy0, rad = float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])
toks=[l.split('#')[0].split() for l in open(cert)]; toks=[t for t in toks if t]
sn,sd=map(int,toks[0]); D=int(toks[1][0]); W=int(toks[2][0]); m=int(toks[3][0])
P=np.array([[int(a) for a in t] for t in toks[4:4+m]],dtype=float)
s=F(sn,sd); g=N*N+k*k; c=F(N*N-k*k,g); si=F(2*k*N,g); w=c+si
X=P[:,0]/D; Y=P[:,1]/D; Wt=P[:,2]
cf,sf=float(c),float(si)
lo=float(w/2); hi=float(s-w/2)
xs=np.linspace(max(lo,cx0-rad),min(hi,cx0+rad),401); ys=np.linspace(max(lo,cy0-rad),min(hi,cy0+rad),401)
best=(1e18,None)
for x in xs:
    dx=X-x
    for y in ys:
        dy=Y-y
        u0=cf*dx+sf*dy; u1=-sf*dx+cf*dy
        cap=Wt[(np.abs(u0)<=0.5)&(np.abs(u1)<=0.5)].sum()
        if cap<best[0]: best=(cap,(x,y))
print("float min capture", best[0]/W, "at", best[1])
# exact re-check at a rational centre near the best point, strictly inside admissible box
x,y=best[1]
qx=F(round(x*10**9),10**9); qy=F(round(y*10**9),10**9)
qx=max(qx,w/2); qy=max(qy,w/2); qx=min(qx,s-w/2); qy=min(qy,s-w/2)
tot=0
for (a,b,wt) in [tuple(int(v) for v in t) for t in toks[4:4+m]]:
    px,py=F(a,D)-qx,F(b,D)-qy
    u0=c*px+si*py; u1=-si*px+c*py
    if abs(u0)<=F(1,2) and abs(u1)<=F(1,2): tot+=wt
inside = (w/2<=qx<=s-w/2) and (w/2<=qy<=s-w/2)
print("exact: angle 2*atan(%d/%d), centre (%s, %s) ~ (%.9f, %.9f), inside container: %s, captured = %d/%d = %.7f" % (k,N,qx,qy,float(qx),float(qy),inside,tot,W,tot/W))
