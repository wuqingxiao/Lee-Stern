"""Where is min_beta [0.292 beta - log2 P(<=beta)] attained? (all instances, chi2, CN11, q-ary start)"""
import json, math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import prob_model as pm
insts = json.load(open(os.path.join(HERE, "instances.json")))
for inst in insts:
    n, k = inst["n"], inst["k"]; cnt = {int(v): c for v, c in inst["counts_centered"].items()}
    ex2 = sum(c * v * v for v, c in cnt.items()) / n
    best = None
    for m in range(max(10, (n - k) // 2), n - k + 1):
        d = m + k + 1; b = pm.beta_2016(d, m, math.sqrt(ex2))
        if b and (best is None or b < best[0]):
            best = (b, m, d)
    b16, m, d = best; norm2 = (m + k) * ex2 + 1
    for tours in (1, 2):
        dist, _ = pm.first_success(pm.qary(m, k), norm2, tours, sim="cn11", succ_min=40)
        dd = dist["chi2"]; acc = 0.0; pts = []
        for b in sorted(dd):
            acc += dd[b]
            if acc > 0:
                pts.append((0.292 * b - math.log2(acc), b, math.log2(acc)))
        pts.sort()
        early = [(b, round(math.log2(p), 1)) for b, p in sorted(dd.items()) if p > 0 and b < b16 - 60][:6]
        print(inst["name"][:28], "tours", tours, "best", [(round(c, 1), b, round(l, 1)) for c, b, l in pts[:3]], "early mass", early, flush=True)