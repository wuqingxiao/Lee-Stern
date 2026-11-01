"""Where does the chi2 model put probability mass for geo2 (2 tours, CN11, q-ary start)?"""
import json, math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import prob_model as pm
from fpylll import BKZ
from fpylll.tools.bkz_simulator import simulate
inst = json.load(open(os.path.join(HERE, "instances.json")))[1]
n, k = inst["n"], inst["k"]; cnt = {int(v): c for v, c in inst["counts_centered"].items()}
ex2 = sum(c * v * v for v, c in cnt.items()) / n
m = 318; d = m + k + 1; norm2 = (m + k) * ex2 + 1
r = pm.qary(m, k); surv = 1.0; acc = 0.0
for b in range(10, 600, 2):
    for _ in range(2):
        r = simulate(r, BKZ.Param(block_size=b, max_loops=1))[0]
        p = pm.p_short(r[d - b], norm2, d, b, "chi2"); acc += surv * p; surv *= 1 - p
    if b % 40 == 10 or 520 <= b <= 560:
        flat = sum(1 for x in r if abs(x - 1.0) < 1e-9)
        print(f"b={b:4d} r[d-b]={r[d-b]:.4g} log2P(<=b)={math.log2(acc) if acc > 0 else float('-inf'):8.2f} "
              f"min_r={min(r):.4g} n_exact1={flat} max_r={max(r):.4g}", flush=True)