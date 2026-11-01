"""Arora-Ge / Groebner estimate for FCSD with exact equation counts.

Variables: after eliminating e1 = s - A e2, the k coordinates of e2 (k-1 after
using the public linear power sum). Equations:
  - n coordinate equations prod_{v in supp(C)} (x - v) = 0, degree |supp(C)|,
    each coordinate an affine form in the variables;
  - power sums sum_j e_j^t = sum_v c_v v^t (mod p), one equation of degree t,
    for t = 2 .. p-2 (t=1 is linear and removes one variable).
Semi-regular Hilbert series H(z) = prod_i (1 - z^{d_i}) / (1 - z)^nv; d_reg is the
first degree whose coefficient is <= 0. Cost = binom(nv + d_reg, d_reg)^omega
(omega = 2 lower bound and 2.81), memory binom(nv+d_reg, d_reg)^2.
Also reports linearisation (XL at degree D) monomial counts.
"""
import json, math

P = 251
insts = json.load(open(r"F:\codex-pj\ac3\research\next_steps_128_20261101\p1_lattice\instances.json"))

def lbinom(a, b):
    return (math.lgamma(a + 1) - math.lgamma(b + 1) - math.lgamma(a - b + 1)) / math.log(2)

def hilbert_dreg(nv, eqs, dmax):
    """eqs: list of (degree, count). Exact integer series up to dmax."""
    s = [0] * (dmax + 1); s[0] = 1
    for deg, cnt in eqs:
        # multiply by (1 - z^deg)^cnt via binomial expansion
        coeff = [0] * (dmax + 1)
        j = 0; b = 1
        while j * deg <= dmax and j <= cnt:
            coeff[j * deg] = b if j % 2 == 0 else -b
            b = b * (cnt - j) // (j + 1); j += 1
        new = [0] * (dmax + 1)
        nz = [(i, c) for i, c in enumerate(coeff) if c]
        for i, si in enumerate(s):
            if si:
                for o, c in nz:
                    if i + o > dmax: break
                    new[i + o] += si * c
        s = new
    # divide by (1 - z)^nv: multiply by sum binom(nv-1+i, i) z^i
    inv = [math.comb(nv - 1 + i, i) for i in range(dmax + 1)]
    for D in range(dmax + 1):
        val = sum(s[i] * inv[D - i] for i in range(D + 1) if s[i])
        if val <= 0:
            return D
    return None

rows = []
for inst in insts:
    n, k = inst["n"], inst["k"]
    supp = len(inst["counts_centered"])
    for use_ps in (False, True):
        nv = k - 1 if use_ps else k
        eqs = [(supp, n)]
        if use_ps:
            eqs += [(t, 1) for t in range(2, P - 1)]
        dreg = hilbert_dreg(nv, eqs, 4 * nv)
        if dreg is None:
            rows.append(dict(name=inst["name"], power_sums=use_ps, nv=nv, deg=supp, dreg=None)); continue
        lm = lbinom(nv + dreg, dreg)
        rows.append(dict(name=inst["name"], power_sums=use_ps, nv=nv, deg=supp, n_eq=n, dreg=dreg,
                         log2_monomials=round(lm, 1), log2_cost_w2=round(2 * lm, 1),
                         log2_cost_w281=round(2.81 * lm, 1), log2_mem=round(2 * lm, 1)))
    # linearisation lower bound: monomials of degree <= supp in nv vars vs number of equations
    rows.append(dict(name=inst["name"], linearisation_log2_monomials_deg_supp=round(lbinom(k + supp, supp), 1),
                     equations=n))
for r in rows: print(r)
json.dump(rows, open(r"F:\codex-pj\ac3\research\next_steps_128_20261101\p2_algebraic\arora_ge.json", "w"), indent=1)
