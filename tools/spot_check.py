# Exact random spot check of the covering claim (sanity only; not a proof).
import numpy as np, random, sys
t=open(sys.argv[1]).read().split()
s_num,s_den,D,W,m=map(int,t[:5]); a=np.array(list(map(int,t[5:])),dtype=np.int64).reshape(m,3)
assert s_den==D  # container side = s_num/D
X,Y,w=a[:,0],a[:,1],a[:,2]
M=1000; DC=D*M; S=s_num*M          # centre coords over DC; container side = S/DC
rng=random.Random(int(sys.argv[2]) if len(sys.argv)>2 else 1)
best=None; n=0; B=1<<10
while n<int(sys.argv[3]) if len(sys.argv)>3 else 20000:
    aa=rng.randrange(0,B+1)              # u = aa/B = tan(theta/2), theta in [0,90]
    cn,sn,G=B*B-aa*aa,2*aa*B,B*B+aa*aa
    # admissible centre: 2G c >= (cn+sn) DC  and  2G c <= 2G S - (cn+sn) DC
    lo=-(-(cn+sn)*DC//(2*G)); hi=(2*G*S-(cn+sn)*DC)//(2*G)
    if lo>hi: continue
    if rng.random()<0.3:   # bias toward walls and corners
        cx=rng.choice([lo,hi,rng.randrange(lo,hi+1)]); cy=rng.choice([lo,hi,rng.randrange(lo,hi+1)])
    else:
        cx=rng.randrange(lo,hi+1); cy=rng.randrange(lo,hi+1)
    dx=X*M-cx; dy=Y*M-cy
    r0=cn*dx+sn*dy; r1=-sn*dx+cn*dy
    inside=(2*np.abs(r0)<=G*DC)&(2*np.abs(r1)<=G*DC)
    cap=int(w[inside].sum()); n+=1
    if best is None or cap<best[0]: best=(cap,aa,cx,cy)
print(f"{n} random admissible poses; min captured = {best[0]}/{W} = {best[0]/W:.7f} at u={best[1]}/{B}, c=({best[2]}/{DC},{best[3]}/{DC})")
print("ALL >= 1" if best[0]>=W else "VIOLATION FOUND")
