"""Most pessimistic stack: l1 proxy (sigma_eff = b_factor * sigma, from p2_lattice_struct)
combined with the calibrated chi2 model (CN11, q-ary start, 8 tours), 1% quantile.
The m (number of samples) is re-optimised for the reduced sigma."""
import json, math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import prob_model as pm
ST = json.load(open(os.path.join(HERE, "..", "p2_lattice_struct", "lattice_struct.json")))
insts = json.load(open(os.path.join(HERE, "instances.json")))
rows = []
for inst in insts:
    n, k = inst["n"], inst["k"]; cnt = {int(v): c for v, c in inst["counts_centered"].items()}
    ex2 = sum(c * v * v for v, c in cnt.items()) / n
    f = ST[inst["name"]]["b_factor"]; ex2s = ex2 * f * f
    best = None
    for m in range(max(10, (n - k) // 2), n - k + 1):
        d = m + k + 1; b = pm.beta_2016(d, m, math.sqrt(ex2s))
        if b and (best is None or b < best[0]):
            best = (b, m, d)
    b16, m, d = best; norm2 = (m + k) * ex2s + 1
    res = dict(name=inst["name"], b_factor=round(f, 3), beta_2016_l1=b16, core_2016_l1=round(0.292 * b16, 1))
    for tours in (2, 8):
        dist, _ = pm.first_success(pm.qary(m, k), norm2, tours, sim="cn11", succ_min=b16 // 2)
        q = pm.quantiles(dist["chi2"])
        res[f"t{tours}"] = dict(beta_med=q[0.5], core_med=round(0.292 * q[0.5], 1), beta_q01=q[0.01], core_q01=round(0.292 * q[0.01], 1))
    print(json.dumps(res), flush=True); rows.append(res)
json.dump(rows, open(os.path.join(HERE, "prob_stack.json"), "w"), indent=1)