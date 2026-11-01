# Zero-guessing hybrid at fixed budget 2^166: guess g coords are zero (prob = C(c0,g)/C(n,g)),
# drop them, run BKZ-beta on dimension n-g+1, beta<=568 (core-SVP<=2^166). Norm unchanged.
import math, json
from scipy import stats
src = open("usvp_2016_check.py", encoding="utf-8").read().split("rows = []")[0]
ns = {}; exec(src, ns); ld = ns["log_delta"]; ROOT = ns["ROOT"]
fz = json.loads((ROOT/"target_s/WP3_C_N1024.json").read_text(encoding="utf-8"))
cnt = fz["counts"]; p=251; n=1024; m=512; c0=cnt[0]
norm2 = sum(min(a,p-a)**2*c for a,c in enumerate(cnt))
def lg_guess(g): return (math.lgamma(c0+1)-math.lgamma(c0-g+1)-math.lgamma(n+1)+math.lgamma(n-g+1))/math.log(2)
best=None; rows=[]
for g in range(0, c0+1, 2):
    d=n-g+1
    for b in range(400, 569, 8):
        reps = 166 - 0.292*b          # log2 of repetitions affordable
        thr2=math.exp(2*((2*b-d-1)*ld(b)+m/d*math.log(p)))
        lp=stats.beta.logcdf(min(1.0,thr2/norm2), b/2,(d-b)/2)/math.log(2)
        tot=min(0.0, lg_guess(g)+lp+reps)
        if best is None or tot>best[0]: best=(tot,g,b,lp,lg_guess(g))
    if g % 40 == 0: rows.append((g, best[0]))
print("best log2 Pr at 2^166: %.1f  (g=%d beta=%d lattice=%.1f guess=%.1f)"%best)
print("assumption F requires <= -136")
