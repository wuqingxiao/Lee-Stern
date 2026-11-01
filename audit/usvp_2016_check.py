"""Primal uSVP check with the 2016 estimate (ADPS16 / Albrecht et al. 2017).

Success condition for BKZ-beta on the Kannan embedding of dimension d:
    sqrt(beta) * sigma <= delta(beta)^(2*beta - d - 1) * q^(m/d)
Cost: core-SVP 2^(0.292*beta). Sanity anchor: Kyber512 primal gives beta=406.
Probability column: chance that the projection of the fixed secret onto the last
beta Gram-Schmidt directions is short enough, at the largest beta with 0.292*beta <= 166.
"""
import json, math, pathlib
from scipy import stats

ROOT = pathlib.Path(r"F:\codex-pj\ac3\research\lee_stern_paper_restructured_20261028")

def log_delta(b):
    return math.log(((b / (2 * math.pi * math.e)) * (math.pi * b) ** (1 / b)) ** (1 / (2 * (b - 1))))

def primal(ns, m, q, sigma):
    """LWE form: ns secret coords, m samples (systematic H). Returns min beta over m."""
    best = None
    for mm in range(64, m + 1, 4):
        d = ns + mm + 1
        for b in range(50, d):
            if math.log(math.sqrt(b) * sigma) <= (2 * b - d - 1) * log_delta(b) + mm / d * math.log(q):
                if best is None or b < best[0]:
                    best = (b, mm)
                break
    return best

def prob_at_budget(n, m, q, norm, budget=166):
    d = n + 1
    b = int(budget / 0.292)
    thr2 = math.exp(2 * ((2 * b - d - 1) * log_delta(b) + m / d * math.log(q)))
    return b, stats.beta.logcdf(min(1.0, thr2 / norm ** 2), b / 2, (d - b) / 2) / math.log(2)

def sigma_of(counts, p, n):
    return math.sqrt(sum(min(a, p - a) ** 2 * c for a, c in enumerate(counts)) / n)

rows = []
k = primal(512, 512, 3329, math.sqrt(1.5))
rows.append(("anchor Kyber512", 512, k[0], 0.292 * k[0], None))
cands = json.loads((ROOT / "security/inputs/candidates.json").read_text(encoding="utf-8"))
cands = cands if isinstance(cands, list) else list(cands.get("candidates", cands).values()) if isinstance(cands.get("candidates", cands), dict) else cands["candidates"]
frozen = json.loads((ROOT / "target_s/WP3_C_N1024.json").read_text(encoding="utf-8"))
items = [(c.get("id") or c.get("name"), c.get("p", 251), c.get("n", 640), c.get("k", 320), c.get("counts") or c.get("C")) for c in cands]
items.append((frozen["name"], 251, 1024, 512, frozen["counts"]))
for name, p, n, kk, counts in items:
    sg = sigma_of(counts, p, n)
    b, mm = primal(kk, n - kk, p, sg)
    _, lp = prob_at_budget(n, n - kk, p, sg * math.sqrt(n))
    rows.append((name, n, b, 0.292 * b, lp))
print(f"{'instance':40s} {'n':>5s} {'beta':>5s} {'log2 coreSVP':>13s} {'log2 Pr@2^166':>14s}")
for r in rows:
    print(f"{r[0]:40s} {r[1]:5d} {r[2]:5d} {r[3]:13.1f} {'' if r[4] is None else f'{r[4]:14.1f}'}")
