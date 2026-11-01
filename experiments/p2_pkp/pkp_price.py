"""PKP-type combinatorial attacks on the IPKP form of FCSD.
IPKP: given H (m x n, m=n-k) over F_p, s, and composition C (multiset of n values),
find v with composition C and Hv = s.
N(j) = number of distinct length-j sequences drawn from multiset C (exact integer DP).
KMP-type (ePrint 2019/412, generalised in ePrint 2022/1749): pick u rows of systematic H,
they involve k+u coordinates; split k+u = k1+k2; lists N(k1), N(k2); match on F_p^u.
T2(u,k1) = max(N(k1), N(k2), N(k1)N(k2)/p^u)  (list elements; memory of same order)
4-list proxy (stand-in for k-list / representation refinements such as SBC and 2026/706):
T4(u) = max(N(k+u)^(1/4), N(k+u)/p^u). This is an optimistic-for-attacker proxy, not a proven bound.
"""
import json, math
from math import comb, log2
INST = r"F:\codex-pj\ac3\research\next_steps_128_20261101\p1_lattice\instances.json"
OUT = r"F:\codex-pj\ac3\research\next_steps_128_20261101\p2_pkp"

def lg(x): return x.bit_length() - 1 + log2(x / (1 << (x.bit_length() - 1))) if x.bit_length() < 1000 else (x.bit_length() - 53) + log2(x >> (x.bit_length() - 53))

def seq_counts(counts, n):
    N = [1] + [0] * n
    tot = 0
    for c in counts:
        tot += c
        new = [0] * (n + 1)
        for j in range(tot + 1):
            s = 0
            for t in range(0, min(c, j) + 1):
                if N[j - t]:
                    s += comb(j, t) * N[j - t]
            new[j] = s
        N = new
    return [lg(x) if x > 0 else float('-inf') for x in N]

res = []
for inst in json.load(open(INST)):
    n, k, p = inst["n"], inst["k"], inst["p"]; m = n - k; lp = log2(p)
    cnt = list(inst["counts_centered"].values()); assert sum(cnt) == n
    L = seq_counts(cnt, n)
    best2 = None
    for u in range(0, m + 1):
        J = k + u
        for k1 in range(J // 2, J + 1):
            k2 = J - k1
            t = max(L[k1], L[k2], L[k1] + L[k2] - u * lp)
            if best2 is None or t < best2[0]: best2 = (t, u, k1, k2, L[J] - u * lp)
    best4 = min((max(L[k + u] / 4, L[k + u] - u * lp), u) for u in range(0, m + 1))
    best8 = min((max(L[k + u] / 8, L[k + u] - u * lp), u) for u in range(0, m + 1))
    H = -sum(c / n * log2(c / n) for c in cnt)
    ops = log2((k + best2[1]) * best2[1] * lp)  # field ops per list element, rough bit-op factor
    r = dict(name=inst["name"], n=n, k=k, logSC=L[n], entropy_per_coord=H,
             kmp_elems=best2[0], kmp_u=best2[1], kmp_k1=best2[2], kmp_k2=best2[3],
             kmp_bitops=best2[0] + ops, proxy4=best4[0], proxy4_u=best4[1],
             proxy8=best8[0], proxy8_u=best8[1])
    res.append(r)
    print(f'{r["name"]:36s} log2|S_C|={L[n]:7.1f} H={H:.2f}  KMP 2^{best2[0]:7.1f} elems (u={best2[1]}, k1={best2[2]}, k2={best2[3]}) ~2^{r["kmp_bitops"]:.1f} bitops  4-list proxy 2^{best4[0]:.1f} (u={best4[1]})  8-list proxy 2^{best8[0]:.1f}')
json.dump(res, open(OUT + r"\pkp_price.json", "w"), indent=1)
