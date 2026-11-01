"""Guess-and-determine Arora-Ge hybrid and coded-BKW for FCSD.

Hybrid: guess g of the k free coordinates to be 0 (exact hypergeometric probability
C(c0,g)/C(n,g), the cheapest value to guess), then Groebner on nv = k-g variables with
n-g degree-|supp| equations (semi-regular d_reg). Cost = guesses x binom(nv+D, D)^2.
BKW: lattice-estimator coded_bkw with sample amplification (m = n-k given samples).
"""
import json, sys, math
sys.path.insert(0, "/home/gql/lattice-estimator")
sys.path.insert(0, "/mnt/f/codex-pj/ac3/research/next_steps_128_20261101/p1_lattice")
from sage.all import ZZ, PolynomialRing, binomial, log
from run_estimator import FixedComposition, LWE

insts = json.load(open("/mnt/f/codex-pj/ac3/research/next_steps_128_20261101/p1_lattice/instances.json"))
R = PolynomialRing(ZZ, "z"); z = R.gen()

def dreg(nv, deg, cnt):
    num = (1 - z**deg)**cnt
    top = min(num.degree(), 6 * nv * deg) + nv + 1
    inv = R([binomial(nv - 1 + i, i) for i in range(top + 1)])
    co = (num.truncate(top + 1) * inv).truncate(top + 1).list()
    for D, c in enumerate(co):
        if c <= 0:
            return D
    return None

def lb(a, b):
    return (math.lgamma(a + 1) - math.lgamma(b + 1) - math.lgamma(a - b + 1)) / math.log(2)

out = []
for inst in insts:
    n, k = inst["n"], inst["k"]; cc = {int(v): c for v, c in inst["counts_centered"].items()}
    supp, c0 = len(cc), cc[0]
    best = None
    for g in sorted(set([0, c0 // 4, c0 // 2, 3 * c0 // 4, c0])):
        nv = k - g
        D = dreg(nv, supp, n - g)
        guess = -(lb(c0, g) - lb(n, g))
        gb = 2 * min(lb(nv + D, D), nv * math.log2(supp))
        row = dict(name=inst["name"], g=g, nv=nv, dreg=D, log2_guess=round(guess, 1), log2_gb=round(gb, 1), log2_total=round(guess + gb, 1))
        print(row, flush=True); out.append(row)
    X = FixedComposition(inst["counts_centered"])
    P = LWE.Parameters(n=k, q=251, Xs=X, Xe=X, m=n - k, tag=inst["name"])
    try:
        b = LWE.coded_bkw(P)
        print(inst["name"], "coded_bkw", b, flush=True)
        out.append(dict(name=inst["name"], coded_bkw_log2_rop=float(log(b["rop"], 2)), log2_m=float(log(b["m"], 2)), log2_mem=float(log(b["mem"], 2))))
    except Exception as ex:
        print(inst["name"], "coded_bkw FAILED", repr(ex), flush=True)
json.dump(out, open("/mnt/f/codex-pj/ac3/research/next_steps_128_20261101/p2_algebraic/hybrid_bkw.json", "w"), indent=1)