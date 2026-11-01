"""Calibration of the probabilistic uSVP model and conservative full-size figures.

(1) PIT check: for every small-scale trial, F_model(observed beta) from its own LLL
    profile (CN11, same tours). Calibrated model -> PIT ~ Uniform(0,1).
    mid-PIT = (F(obs-2) + F(obs)) / 2 handles the even block-size grid.
(2) Full-size instances, chi2 and beta models, CN11 q-ary start, tours 1/2/4/8:
    core-SVP (0.292 beta) at the median and at the 1% quantile of the first-success
    block size, and the cost per unit success probability
        min_beta [0.292*beta - log2 P(success by beta)]
    (one SVP call per block size, core-SVP convention).
    Main column counts success only for beta >= beta_2016/2: below that the last
    beta Gram-Schmidt vectors still come from the unit (identity) part of the q-ary
    basis, the projected target is not unique there, and the Gaussian chi2 tail
    (~2^-230) is an artefact. The column "incl_tail" counts from beta = 40 anyway.
    The real event behind that tail (many zero coordinates at the end) is priced by
    the exact zero-guess hybrid in p2_lattice_struct (optimum g = 0).
(3) Cross-check: BSW18 probabilistic simulator for geo2 (binding instance), 2 tours.
"""
import json, math, os, statistics as st, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import prob_model as pm
from scipy.stats import kstest

def cdf(dist, b):
    return sum(p for bb, p in dist.items() if bb <= b)

def cost_per_success(dd):
    acc = 0.0; best = None
    for b in sorted(dd):
        acc += dd[b]
        if acc > 0:
            v = 0.292 * b - math.log2(acc)
            best = v if best is None or v < best else best
    return round(best, 1) if best is not None else None

def summarize(dd):
    q = pm.quantiles(dd)
    return dict(beta_med=q[0.5], core_med=round(0.292 * q[0.5], 1),
                beta_q01=q[0.01], core_q01=round(0.292 * q[0.01], 1),
                cost_per_success=cost_per_success(dd))

out = {"pit": {}, "scale": []}
rows = [json.loads(l) for l in open(os.path.join(HERE, "bkz_scale_v2.jsonl"))]
for model in pm.MODELS:
    pits = []; per_k = {}
    for t in rows:
        if not t.get("beta"):
            continue
        r0 = [2.0 ** (2 * x) for x in t["lll_profile"]]
        dist, _ = pm.first_success(r0, t["norm2"], t.get("tours_per_beta", 2), sim="cn11")
        lo, hi = cdf(dist[model], t["beta"] - 2), cdf(dist[model], t["beta"])
        mid = 0.5 * (lo + hi); pits.append(mid); per_k.setdefault(t["k"], []).append(round(mid, 3))
    ks = kstest(pits, "uniform")
    out["pit"][model] = dict(n=len(pits), mean=round(st.mean(pits), 3), median=round(st.median(pits), 3),
                             ks_stat=round(ks.statistic, 3), ks_p=round(ks.pvalue, 3), per_k=per_k)
    print(model, json.dumps(out["pit"][model]), flush=True)

insts = json.load(open(os.path.join(HERE, "instances.json")))
for idx, inst in enumerate(insts):
    n, k = inst["n"], inst["k"]
    cnt = {int(v): c for v, c in inst["counts_centered"].items()}
    ex2 = sum(c * v * v for v, c in cnt.items()) / n
    best = None
    for m in range(max(10, (n - k) // 2), n - k + 1):
        d = m + k + 1; b = pm.beta_2016(d, m, math.sqrt(ex2))
        if b and (best is None or b < best[0]):
            best = (b, m, d)
    b16, m, d = best; norm2 = (m + k) * ex2 + 1
    res = dict(name=inst["name"], beta_2016=b16, core_2016=round(0.292 * b16, 1), m=m, d=d)
    for tours in (1, 2, 4, 8):
        dist, _ = pm.first_success(pm.qary(m, k), norm2, tours, sim="cn11", succ_min=b16 // 2)
        tail, _ = pm.first_success(pm.qary(m, k), norm2, tours, sim="cn11", succ_min=40)
        for model in pm.MODELS:
            s = summarize(dist[model]); s["cost_per_success_incl_tail"] = cost_per_success(tail[model])
            res[f"t{tours}_{model}"] = s
    if idx == 1:
        t0 = time.time()
        dist, _ = pm.first_success(pm.qary(m, k), norm2, 2, sim="bsw18", succ_min=b16 // 2)
        for model in pm.MODELS:
            res[f"bsw18_t2_{model}"] = summarize(dist[model])
        res["bsw18_seconds"] = round(time.time() - t0, 1)
    print(json.dumps(res), flush=True); out["scale"].append(res)
json.dump(out, open(os.path.join(HERE, "prob_calib.json"), "w"), indent=1)