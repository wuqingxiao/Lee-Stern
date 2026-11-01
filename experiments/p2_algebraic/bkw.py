"""coded-BKW for FCSD instances (lattice-estimator). The estimator cannot amplify samples
for a nonzero-mean distribution, so the composition marginal is used with its mean set to 0
and the same second moment (stddev^2 + mean^2), which only helps the attacker slightly."""
import json, sys
sys.path.insert(0, "/home/gql/lattice-estimator")
sys.path.insert(0, "/mnt/f/codex-pj/ac3/research/next_steps_128_20261101/p1_lattice")
from sage.all import log, sqrt, RR
from run_estimator import FixedComposition, LWE
from estimator import ND
insts = json.load(open("/mnt/f/codex-pj/ac3/research/next_steps_128_20261101/p1_lattice/instances.json"))
out = []
for inst in insts:
    n, k = inst["n"], inst["k"]
    X = FixedComposition(inst["counts_centered"])
    sd = float(RR(sqrt(X.stddev**2 + X.mean**2)))
    G = ND.DiscreteGaussian(float(sd))
    P = LWE.Parameters(n=k, q=251, Xs=G, Xe=G, m=n - k, tag=inst["name"])
    try:
        b = LWE.coded_bkw(P)
        r = dict(name=inst["name"], stddev=float(sd), log2_rop=round(float(log(b["rop"], 2)), 1),
                 log2_m=round(float(log(b["m"], 2)), 1), log2_mem=round(float(log(b["mem"], 2)), 1))
    except Exception as ex:
        r = dict(name=inst["name"], error=repr(ex))
    print(r, flush=True); out.append(r)
json.dump(out, open("/mnt/f/codex-pj/ac3/research/next_steps_128_20261101/p2_algebraic/bkw.json", "w"), indent=1)