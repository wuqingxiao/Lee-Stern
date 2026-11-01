"""Semi-regular degree of regularity for the Arora-Ge system of FCSD (Sage/FLINT).

H(z) = prod_i (1 - z^{d_i}) / (1 - z)^nv. Equations: n coordinate equations of degree
|supp(C)| (each coordinate is a variable or an affine form after eliminating e1),
optionally the power sums t = 2..p-2 (one equation of degree t each; t = 1 is linear
and removes a variable). d_reg = first index with coefficient <= 0.
Cost (F5 / linear algebra, omega in {2, 2.81}): binom(nv + d_reg, d_reg)^omega.
"""
import json, sys
from sage.all import ZZ, PolynomialRing, binomial, log, RR

P = 251
insts = json.load(open("/mnt/f/codex-pj/ac3/research/next_steps_128_20261101/p1_lattice/instances.json"))
R = PolynomialRing(ZZ, "z"); z = R.gen()

def dreg(nv, eqs):
    num = R(1)
    for deg, cnt in eqs:
        num *= (1 - z**deg)**cnt
    top = num.degree() + nv + 1
    # (1-z)^-nv truncated: coefficients binom(nv-1+i, i)
    inv = R([binomial(nv - 1 + i, i) for i in range(top + 1)])
    ser = (num * inv).truncate(top + 1)
    co = ser.list()
    for D, c in enumerate(co):
        if c <= 0:
            return D
    return None

rows = []
for inst in insts:
    n, k = inst["n"], inst["k"]; supp = len(inst["counts_centered"])
    for use_ps in (False, True):
        nv = k - 1 if use_ps else k
        eqs = [(supp, n)] + ([(t, 1) for t in range(2, P - 1)] if use_ps else [])
        D = dreg(nv, eqs)
        lm = float(log(binomial(nv + D, D), 2)) if D is not None else None
        row = dict(name=inst["name"], n=n, k=k, support=supp, power_sums=use_ps, nv=nv, dreg=D,
                   log2_monomials=None if lm is None else round(lm, 1),
                   log2_cost_omega2=None if lm is None else round(2 * lm, 1),
                   log2_cost_omega281=None if lm is None else round(2.81 * lm, 1))
        print(row, flush=True); rows.append(row)
json.dump(rows, open("/mnt/f/codex-pj/ac3/research/next_steps_128_20261101/p2_algebraic/arora_ge_dreg.json", "w"), indent=1)