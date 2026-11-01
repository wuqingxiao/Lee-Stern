import math, json, pathlib
from scipy import stats
import importlib.util
spec = importlib.util.spec_from_file_location("u", "usvp_2016_check.py")
src = open("usvp_2016_check.py", encoding="utf-8").read().split("rows = []")[0]
ns = {}; exec(src, ns)
ld = ns["log_delta"]; sigma_of = ns["sigma_of"]; ROOT = ns["ROOT"]
def curve(n, m, q, norm, budget=166):
    d = n + 1; best = (-1e9, None)
    for b in range(60, int(budget/0.292)+1):
        thr2 = math.exp(2*((2*b-d-1)*ld(b) + m/d*math.log(q)))
        lp = stats.beta.logcdf(min(1.0, thr2/norm**2), b/2, (d-b)/2)/math.log(2)
        amp = lp + (budget - 0.292*b)   # repeat on re-randomised lattices, independence heuristic
        amp = min(amp, 0.0)
        if amp > best[0]: best = (amp, b, lp)
    return best
c = json.loads((ROOT/"security/inputs/candidates.json").read_text(encoding="utf-8"))
c = c["candidates"] if isinstance(c, dict) else c
c = list(c.values()) if isinstance(c, dict) else c
fz = json.loads((ROOT/"target_s/WP3_C_N1024.json").read_text(encoding="utf-8"))
items = [(x.get("id") or x.get("name"), 640, 320, x.get("counts") or x.get("C")) for x in c] + [(fz["name"], 1024, 512, fz["counts"])]
print("instance                                   log2Pr(single,beta=568)  log2Pr(repeat,best)  beta*")
for name, n, k, cnt in items:
    s = sigma_of(cnt, 251, n); norm = s*math.sqrt(n)
    d=n+1; b=568; thr2=math.exp(2*((2*b-d-1)*ld(b)+(n-k)/d*math.log(251)))
    single = stats.beta.logcdf(min(1.0,thr2/norm**2), b/2,(d-b)/2)/math.log(2)
    amp, bb, lp = curve(n, n-k, 251, norm)
    print(f"{name:42s} {single:10.1f} {amp:20.1f} {bb:6d}")
