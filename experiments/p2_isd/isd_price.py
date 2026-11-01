"""P2: LMMT-ISD pricing in bit operations for the four FCSD instances.

Cost model: Li-Wang LMMT prime-field model, paper eq. (5.4) (bit operations:
Gaussian elimination beta^2 m^2 (n+1), list building, sorting, completion).
Two success factors per parameter point (ell, ell2 = ell-1, w1 = w2 = w):
  P_th : paper's Lee-shell split factor D1 D2 S(j,tau) / S(n,w_C)   (Table 5.1 model)
  P_nec: exact Pr[wt_L(e_I) <= 4w] for the fixed-composition key, I a uniform
         d-subset (d = k+ell). Necessary for any four-leaf success (triangle
         inequality), so C_iter / P_nec is a LOWER bound on this family's cost.
Also: generic MITM / representation scale over S_C, and multi-target adjustment.
"""
import json, math, sys, time
import numpy as np

P = 251; R = P // 2; BETA = math.log2(P); MU = (P * P - 1) / (4 * P)
BASE = r"F:\codex-pj\ac3\research\next_steps_128_20261101"
insts = json.load(open(BASE + r"\p1_lattice\instances.json"))

def lee_shell_logs(n, t):
    """S[len][a] = log2 #vectors of length len, Lee weight a; V = cumulative."""
    S = np.full((n + 1, t + 1), -np.inf); V = np.full((n + 1, t + 1), -np.inf)
    row = [1] + [0] * t
    for length in range(n + 1):
        if length:
            prev = row; row = [0] * (t + 1); window = 0
            for a in range(t + 1):
                if a: window += prev[a - 1]
                if a - R - 1 >= 0: window -= prev[a - R - 1]
                row[a] = prev[a] + 2 * window
        acc = 0
        for a, x in enumerate(row):
            acc += x
            if x: S[length, a] = math.log2(x)
            if acc: V[length, a] = math.log2(acc)
    return S, V

def comp_subset_logs(m_lee, n, wmax):
    """N[j][w] = log2 #(j-subsets of positions with total Lee weight w) for the key."""
    LN2 = math.log(2)
    N = np.full((n + 1, wmax + 1), -np.inf); N[0, 0] = 0.0
    for i, mi in m_lee.items():
        new = np.full_like(N, -np.inf)
        for x in range(mi + 1):
            lc = (math.lgamma(mi + 1) - math.lgamma(x + 1) - math.lgamma(mi - x + 1)) / LN2
            sh = i * x
            if sh > wmax: break
            src = N[: n + 1 - x, : wmax + 1 - sh] + lc
            dst = new[x:, sh:]
            np.logaddexp2(dst, src, out=dst)
        N = new
    return N

def logadd(*v):
    m = max(v); return m + math.log2(sum(2 ** (x - m) for x in v))

def lbin(a, b):
    return (math.lgamma(a + 1) - math.lgamma(b + 1) - math.lgamma(a - b + 1)) / math.log(2)

def price(inst):
    n, k = inst["n"], inst["k"]; m = n - k
    cc = {int(v): c for v, c in inst["counts_centered"].items()}
    t = sum(abs(v) * c for v, c in cc.items())
    m_lee = {}
    for v, c in cc.items(): m_lee[abs(v)] = m_lee.get(abs(v), 0) + c
    t0 = time.time()
    S, V = lee_shell_logs(n, t)
    N = comp_subset_logs(m_lee, n, t)
    Ncum = np.logaddexp2.accumulate(N, axis=1)
    lcn = np.array([lbin(n, d) for d in range(n + 1)])
    G = 2 * math.log2(BETA) + 2 * math.log2(m) + math.log2(n + 1)
    best_th = best_c = None; pts = 0
    for ell in range(2, n - k):
        if (k + ell) % 2: continue
        h = (k + ell) // 2; j = n - k - ell; u = ell - 1; d = 2 * h
        low = max(1, (t - j * R) // 4 + 1); high = min((t - 1) // 4, (h * R - 1) // 2)
        for w in range(low, high + 1):
            tau = t - 4 * w
            A = S[h, w]; D = S[h, 2 * w]; U = V[h, w]
            B = math.log2(BETA) + logadd(math.log2(4 * ell) + U, A, 2 * A - u * BETA)
            srt = 1 + A + math.log2(A)
            Rc = math.log2(BETA * BETA + BETA) + math.log2(k + ell) + math.log2(tau + 1) - math.log2(MU) + 4 * A - (ell + u) * BETA
            Cit = logadd(G, B, srt, Rc)
            Pth = 2 * D + S[j, tau] - S[n, t]
            Pnec = min(0.0, Ncum[d, 4 * w] - lcn[d])
            mem = A + math.log2(h * BETA)
            rec_th = (Cit - Pth, ell, w, Cit, Pth, mem)
            rec_c = (Cit - Pnec, ell, w, Cit, Pnec, mem)
            pts += 1
            if best_th is None or rec_th < best_th: best_th = rec_th
            if best_c is None or rec_c < best_c: best_c = rec_c
    lsc = (math.lgamma(n + 1) - sum(math.lgamma(c + 1) for c in cc.values())) / math.log(2)
    out = dict(name=inst["name"], n=n, k=k, w_C=t, points=pts, seconds=round(time.time() - t0, 1),
               log2_S_C=lsc, log2_lee_shell=float(S[n, t]), log2_p_m=m * BETA,
               lmmt_shell=dict(zip(["log2_cost", "ell", "w", "log2_Citer", "log2_P", "log2_mem_bits"], best_th)),
               lmmt_fcsd_lower=dict(zip(["log2_cost", "ell", "w", "log2_Citer", "log2_P", "log2_mem_bits"], best_c)),
               mitm_half=lsc / 2, representation_0283=0.283 * lsc)
    return out, S, N, lcn

if __name__ == "__main__":
    res = []
    for inst in insts:
        o, S, N, lcn = price(inst); res.append(o)
        print(json.dumps(o), flush=True)
        if inst["name"].endswith("a17a2941ee99"):
            d, wi = 408, 664
            pn = np.logaddexp2.reduce(N[d, : wi + 1]) - lcn[d]
            print("check Table 5.2 (wide700, d=408, w_int=664): exact log2 P_nec =", round(float(pn), 3), "(paper Chernoff bound B_C = -164.860)")
    json.dump(res, open(BASE + r"\p2_isd\isd_price.json", "w"), indent=1)
