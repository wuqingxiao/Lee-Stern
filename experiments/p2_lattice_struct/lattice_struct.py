"""P2 structure-aware lattice checks (lattice-estimator + 2016 criterion).

(a) public sum row 1^T e = a_C appended to H: LWE(n=k-1, m=n-k+1).
(b) l1 / Lee proxy: relative shortness of e in l1 vs l2 against the Gaussian heuristic
    of the Construction-A lattice; pessimistic proxy scales sigma by r1/r2.
(c) exact-probability zero-guessing hybrid: guess g positions are 0 with probability
    C(c0,g)/C(n,g); remaining LWE(k-g, m) with sigma_g = ||e|| / sqrt(n-g).
    cost(g) = 0.292*beta(g) + log2(1/Pr).  Per-coordinate trade-off vs MITM entropy.
(d) dual / dual_hybrid with RC.ADPS16 (core-SVP) for confirmation.
"""
import json, math, sys
sys.path.insert(0, "/home/gql/lattice-estimator")
sys.path.insert(0, "/mnt/f/codex-pj/ac3/research/next_steps_128_20261101/p1_lattice")
from sage.all import log, RR
from run_estimator import FixedComposition, LWE
from estimator import ND, RC

BASE = "/mnt/f/codex-pj/ac3/research/next_steps_128_20261101/"
insts = json.load(open(BASE + "p1_lattice/instances.json"))
P = 251

def ld(b):
    return math.log(((b / (2 * math.pi * math.e)) * (math.pi * b) ** (1 / b)) ** (1 / (2 * (b - 1))))

def beta_2016(ns, m, sigma, q=P):
    best = None
    for mm in range(32, m + 1, 2):
        d = ns + mm + 1
        lo, hi = 40, d
        if math.log(math.sqrt(hi - 1) * sigma) > (2 * (hi - 1) - d - 1) * ld(hi - 1) + mm / d * math.log(q):
            continue
        while lo < hi:
            mid = (lo + hi) // 2
            if math.log(math.sqrt(mid) * sigma) <= (2 * mid - d - 1) * ld(mid) + mm / d * math.log(q):
                hi = mid
            else:
                lo = mid + 1
        best = lo if best is None else min(best, lo)
    return best

def lb(a, b):
    return (math.lgamma(a + 1) - math.lgamma(b + 1) - math.lgamma(a - b + 1)) / math.log(2)

out = {}
for inst in insts:
    name, n, k = inst["name"], inst["n"], inst["k"]; m = n - k
    cc = {int(v): c for v, c in inst["counts_centered"].items()}
    norm2 = sum(v * v * c for v, c in cc.items()); l1 = sum(abs(v) * c for v, c in cc.items())
    sigma = math.sqrt(norm2 / n); c0 = cc[0]
    H = -sum(c / n * math.log2(c / n) for c in cc.values())
    r = dict(sigma=sigma, norm2=math.sqrt(norm2), l1=l1, entropy_bits=H)
    X = FixedComposition(inst["counts_centered"])
    # (a) sum row
    Pa = LWE.Parameters(n=k - 1, q=P, Xs=X, Xe=X, m=m + 1, tag=name + "+sumrow")
    ua = LWE.primal_usvp(Pa, red_cost_model=RC.ADPS16)
    r["a_sumrow_coreSVP"] = float(log(ua["rop"], 2)); r["a_sumrow_beta"] = int(ua["beta"])
    u0 = LWE.primal_usvp(LWE.Parameters(n=k, q=P, Xs=X, Xe=X, m=m), red_cost_model=RC.ADPS16)
    r["base_coreSVP"] = float(log(u0["rop"], 2)); r["base_beta"] = int(u0["beta"])
    # (b) l1 proxy on Construction-A lattice (dim n, det p^m)
    logdet_n = m / n * math.log(P)
    gh2 = math.exp(logdet_n + math.lgamma(n / 2 + 1) / n) / math.sqrt(math.pi)
    gh1 = math.exp(logdet_n + math.lgamma(n + 1) / n) / 2
    r1, r2 = l1 / gh1, math.sqrt(norm2) / gh2
    f = r1 / r2
    G = ND.DiscreteGaussian(sigma * min(1.0, f))
    ub = LWE.primal_usvp(LWE.Parameters(n=k, q=P, Xs=G, Xe=G, m=m), red_cost_model=RC.ADPS16)
    r.update(b_gh2=gh2, b_gh1=gh1, b_ratio_l2=r2, b_ratio_l1=r1, b_factor=f,
             b_proxy_coreSVP=float(log(ub["rop"], 2)), b_proxy_beta=int(ub["beta"]))
    # (c) exact zero-guess hybrid
    rows = []
    b0 = beta_2016(k, m, sigma)
    for g in sorted(set(list(range(0, min(c0, 80) + 1, 4)) + [c0 // 2, c0])):
        sg = math.sqrt(norm2 / (n - g))
        b = beta_2016(k - g, m, sg)
        guess = -(lb(c0, g) - lb(n, g))
        rows.append(dict(g=g, beta=b, coreSVP=0.292 * b, guess=guess, total=0.292 * b + guess))
    best = min(rows, key=lambda x: x["total"])
    sav = (rows[0]["coreSVP"] - [x for x in rows if x["g"] == 40][0]["coreSVP"]) / 40
    r.update(c_beta0_2016=b0, c_best=best, c_saving_per_coord=sav,
             c_zero_guess_bits_first=-math.log2(c0 / n), c_mitm_bits_per_coord=H / 2, c_rows=rows)
    # (d) dual with ADPS16
    for nm, fn in [("dual", LWE.dual), ("dual_hybrid", LWE.dual_hybrid)]:
        try:
            dv = fn(LWE.Parameters(n=k, q=P, Xs=X, Xe=X, m=m), red_cost_model=RC.ADPS16)
            r["d_" + nm + "_ADPS16"] = float(log(dv["rop"], 2))
        except Exception as ex:
            r["d_" + nm + "_ADPS16"] = repr(ex)
    out[name] = r
    print(name, {kk: (round(v, 2) if isinstance(v, float) else v) for kk, v in r.items() if kk != "c_rows"}, flush=True)
json.dump(out, open(BASE + "p2_lattice_struct/lattice_struct.json", "w"), indent=1, default=float)