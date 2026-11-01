"""RC.MATZOV bit-operation cost at the most pessimistic stack (l1 proxy + chi2 model, 8 tours,
1% quantile; prob_stack.json). Lattice dimension d re-derived as in prob_stack.py."""
import json, math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, "/home/gql/lattice-estimator")
import prob_model as pm
from estimator import RC
from sage.all import log
ST = json.load(open(os.path.join(HERE, "..", "p2_lattice_struct", "lattice_struct.json")))
stack = {r["name"]: r for r in json.load(open(os.path.join(HERE, "prob_stack.json")))}
rows = []
for inst in json.load(open(os.path.join(HERE, "instances.json"))):
    n, k = inst["n"], inst["k"]; cnt = {int(v): c for v, c in inst["counts_centered"].items()}
    ex2 = sum(c * v * v for v, c in cnt.items()) / n * ST[inst["name"]]["b_factor"] ** 2
    best = None
    for m in range(max(10, (n - k) // 2), n - k + 1):
        d = m + k + 1; b = pm.beta_2016(d, m, math.sqrt(ex2))
        if b and (best is None or b < best[0]):
            best = (b, m, d)
    _, m, d = best; s = stack[inst["name"]]
    out = dict(name=inst["name"], d=d)
    for tag, b in [("t2_med", s["t2"]["beta_med"]), ("t8_q01", s["t8"]["beta_q01"])]:
        out[tag] = dict(beta=b, adps16=round(float(log(RC.ADPS16(b, d), 2)), 1), matzov=round(float(log(RC.MATZOV(b, d), 2)), 1))
    print(json.dumps(out), flush=True); rows.append(out)
json.dump(rows, open(os.path.join(HERE, "matzov_stack.json"), "w"), indent=1)