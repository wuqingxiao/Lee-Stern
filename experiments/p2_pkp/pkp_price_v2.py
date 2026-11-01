"""PKP pricing, v2: adds structural k-tree (Wagner) with disjoint coordinate blocks,
and an extreme non-structural k-tree proxy. Reuses seq_counts from pkp_price.py.
Wagner(a): 2^a lists on disjoint blocks of j=(k+u)/2^a coordinates; list size L=N(j);
needs log L >= u*log p/(a+1); time ~ 2^a * 2^(u log p/(a+1)); plus the solution-count term N(k+u)/p^u.
Extreme proxy: list size N(k+u)/2^a in log (no block structure), a<=15. Not a real algorithm;
upper-bounds what any list-merging refinement could gain under this accounting."""
import json
from math import log2
src = open(__file__.replace("pkp_price_v2.py", "pkp_price.py"), encoding="utf8").read().split("res = []")[0]
exec(src)
old = {r["name"]: r for r in json.load(open(OUT + r"\pkp_price.json"))}
res = []
for inst in json.load(open(INST)):
    n, k, p = inst["n"], inst["k"], inst["p"]; m = n - k; lp = log2(p)
    L = seq_counts(list(inst["counts_centered"].values()), n)
    wag = None; ext = None
    for a in range(1, 16):
        for u in range(0, m + 1):
            sol = L[k + u] - u * lp
            j = (k + u) >> a
            if a <= 10 and j >= 1 and L[j] >= u * lp / (a + 1):
                t = max(a + u * lp / (a + 1), sol)
                if wag is None or t < wag[0]: wag = (t, a, u)
            lsz = max(L[k + u] / 2 ** a, u * lp / (a + 1))
            t = max(a + lsz, sol)
            if ext is None or t < ext[0]: ext = (t, a, u)
    r = dict(old[inst["name"]]); r.update(wagner=wag[0], wagner_a=wag[1], wagner_u=wag[2],
                                       extreme=ext[0], extreme_a=ext[1], extreme_u=ext[2])
    res.append(r)
    print(f'{r["name"]:36s} KMP {r["kmp_bitops"]:7.1f}  4-list {r["proxy4"]:6.1f}  8-list {r["proxy8"]:6.1f}  Wagner {wag[0]:6.1f} (a={wag[1]},u={wag[2]})  extreme {ext[0]:6.1f} (a={ext[1]},u={ext[2]})')
json.dump(res, open(OUT + r"\pkp_price_v2.json", "w"), indent=1)
