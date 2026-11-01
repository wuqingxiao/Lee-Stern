"""Model predictions (chi2 / beta, CN11 and BSW18, q-ary start, 2 tours) for k not yet run."""
import json, os, sys, math
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import bkz_scale_v2 as bk
import prob_model as pm
rows = []
for k in [60, 65, 70, 80]:
    comp = bk.scaled(2 * k)
    norm2 = sum(c * v * v for v, c in comp.items()) + 1
    d = 2 * k + 1
    res = dict(k=k, d=d, norm2=norm2, beta_2016=pm.beta_2016(d, k, math.sqrt((norm2 - 1) / (2 * k))))
    for sim in ("cn11", "bsw18"):
        dist, _ = pm.first_success(pm.qary(k, k), norm2, 2, sim=sim)
        for m in pm.MODELS:
            res[f"{sim}_{m}"] = {str(q): v for q, v in pm.quantiles(dist[m]).items()}
    print(json.dumps(res), flush=True); rows.append(res)
json.dump(rows, open(os.path.join(HERE, "prob_pred.json"), "w"), indent=1)