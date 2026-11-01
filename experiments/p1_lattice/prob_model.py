"""Probabilistic success model for primal uSVP (DDGR20 / PV21 style).

Validation: every small-scale trial in bkz_scale_v2.jsonl starts from its saved LLL
profile and the same progressive schedule is simulated (beta = 10, 12, ...,
tours_per_beta tours each). After every simulated tour, the probability that the
projection of the embedded target t = (e1, e2, 1) onto the last beta Gram-Schmidt
directions is shorter than ||b*_{d-beta}|| is
  chi2 : F_chi2(beta)( r_{d-beta} * d / ||t||^2 )              Gaussian target (DDGR20/PV21)
  beta : F_Beta(beta/2,(d-beta)/2)( r_{d-beta} / ||t||^2 )    fixed norm, uniform direction
(fixed composition fixes ||t|| exactly). Tours are independent trials (PV21).
Output: distribution of the first-success block size (mixture over trials).
Profile simulators: CN11 (fpylll simulate; what lattice-estimator uses) and BSW18
(fpylll simulate_prob). The canonical q-ary start profile is checked against the
LLL start, then used for the four full-size instances.
Usage: python prob_model.py [small|scale|all]
"""
import json, math, os, statistics as st, sys, time
from fpylll import BKZ
from fpylll.tools.bkz_simulator import simulate, simulate_prob
from scipy.stats import chi2, beta as Beta

HERE = os.path.dirname(os.path.abspath(__file__))
Q = 251
MODELS = ("chi2", "beta")

def log_delta(b):
    return math.log(((b / (2 * math.pi * math.e)) * (math.pi * b) ** (1 / b)) ** (1 / (2 * (b - 1))))

def beta_2016(d, m, sigma):
    """2016 criterion for a Kannan embedding of dimension d with m q-vectors."""
    for b in range(40, d):  # delta(b) formula is meaningless for small b
        if math.log(math.sqrt(b) * sigma) <= (2 * b - d - 1) * log_delta(b) + m / d * math.log(Q):
            return b
    return None

def p_short(x, norm2, d, b, model):
    if model == "chi2":
        return float(chi2.cdf(x * d / norm2, b))
    return float(Beta.cdf(x / norm2, b / 2.0, (d - b) / 2.0))

def first_success(r0, norm2, tours, sim="cn11", b_start=10, step=2, seed=1, eps=1e-7, succ_min=0):
    """succ_min: success events are only counted for block size >= succ_min.
    For full-size instances the Gaussian (chi2) tail at tiny b gives a spurious ~2^-27
    "success" (projection onto the last few unit vectors); it is impossible because the
    Kannan coordinate alone has norm 1, and lifting from such a tiny block fails anyway."""
    d = len(r0); r = list(r0)
    surv = {m: 1.0 for m in MODELS}; dist = {m: {} for m in MODELS}
    b = b_start
    while b < d and max(surv.values()) > eps:
        for _ in range(tours):
            par = BKZ.Param(block_size=b, max_loops=1)
            if sim == "cn11":
                r = simulate(r, par)[0]
            else:
                r = simulate_prob(r, par, prng_seed=seed)[0]; seed += 1
            x = r[d - b]
            for m in MODELS:
                p = p_short(x, norm2, d, b, m) if b >= succ_min else 0.0
                dist[m][b] = dist[m].get(b, 0.0) + surv[m] * p
                surv[m] *= 1.0 - p
        b += step
    return dist, surv

def quantiles(dist, qs=(0.01, 0.1, 0.25, 0.5, 0.75, 0.9)):
    out = {}; acc = 0.0; items = sorted(dist.items()); i = 0
    for q in qs:
        while i < len(items) and acc + items[i][1] < q:
            acc += items[i][1]; i += 1
        out[q] = items[i][0] if i < len(items) else None
    return out

def mixture(dists):
    mix = {}
    for dd in dists:
        for b, p in dd.items():
            mix[b] = mix.get(b, 0.0) + p / len(dists)
    return mix

def qary(m, n):
    return [float(Q * Q)] * m + [1.0] * n + [1.0]

def small():
    rows = [json.loads(l) for l in open(os.path.join(HERE, "bkz_scale_v2.jsonl"))]
    by_k = {}
    for t in rows:
        by_k.setdefault(t["k"], []).append(t)
    out = []
    for k in sorted(by_k):
        ts = [t for t in by_k[k] if t.get("beta")]
        if not ts:
            continue
        d = 2 * k + 1; tours = ts[0].get("tours_per_beta", 2)
        obs = sorted(t["beta"] for t in ts)
        res = dict(k=k, d=d, n_trials=len(ts), obs=obs, obs_median=st.median(obs))
        for sim in ("cn11", "bsw18"):
            for start in ("lll", "qary"):
                per = {m: [] for m in MODELS}
                for j, t in enumerate(ts):
                    r0 = [2.0 ** (2 * x) for x in t["lll_profile"]] if start == "lll" else qary(k, k)
                    dist, _ = first_success(r0, t["norm2"], tours, sim=sim, seed=1 + 1000 * j)
                    for m in MODELS:
                        per[m].append(dist[m])
                for m in MODELS:
                    qv = quantiles(mixture(per[m]))
                    res[f"{sim}_{start}_{m}"] = {str(q): v for q, v in qv.items()}
        res["beta_2016"] = beta_2016(d, k, math.sqrt((ts[0]["norm2"] - 1) / (2 * k)))
        print(json.dumps(res), flush=True); out.append(res)
    # v1 check (bkz_small.txt): 8 tours per block size, only row 0 checked, no saved profile
    v1_obs = {40: [34, 22, 18], 50: [50, 48, 46]}
    for k, obs in v1_obs.items():
        norm2 = by_k[k][0]["norm2"]
        dist, _ = first_success(qary(k, k), norm2, 8, sim="cn11")
        res = dict(k=k, d=2 * k + 1, check="v1_8tours", obs=sorted(obs), obs_median=st.median(obs))
        for m in MODELS:
            res[f"cn11_qary_{m}_t8"] = {str(q): v for q, v in quantiles(dist[m]).items()}
        print(json.dumps(res), flush=True); out.append(res)
    return out

def scale(tours_list=(2,), sims=("cn11",)):
    insts = json.load(open(os.path.join(HERE, "instances.json")))
    out = []
    for inst in insts:
        n, k = inst["n"], inst["k"]
        cnt = {int(v): c for v, c in inst["counts_centered"].items()}
        ex2 = sum(c * v * v for v, c in cnt.items()) / n
        sigma = math.sqrt(ex2)
        best = None
        for m in range(max(10, (n - k) // 2), n - k + 1):
            d = m + k + 1; b = beta_2016(d, m, sigma)
            if b and (best is None or b < best[0]):
                best = (b, m, d)
        b16, m, d = best
        norm2 = (m + k) * ex2 + 1
        res = dict(name=inst["name"], n=n, k=k, sigma=sigma, m=m, d=d, beta_2016=b16, core_2016=round(0.292 * b16, 1))
        for sim in sims:
            for tours in tours_list:
                t0 = time.time()
                dist, surv = first_success(qary(m, k), norm2, tours, sim=sim, b_start=10, step=2)
                for mdl in MODELS:
                    qv = quantiles(dist[mdl])
                    res[f"{sim}_t{tours}_{mdl}"] = {str(q): v for q, v in qv.items()}
                    res[f"{sim}_t{tours}_{mdl}_core_median"] = round(0.292 * qv[0.5], 1)
                    res[f"{sim}_t{tours}_{mdl}_core_q01"] = round(0.292 * qv[0.01], 1)
                res[f"{sim}_t{tours}_time"] = round(time.time() - t0, 1)
        print(json.dumps(res), flush=True); out.append(res)
    return out

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    result = {}
    if mode in ("small", "all"):
        result["small"] = small()
    if mode in ("scale", "all"):
        result["scale"] = scale(tours_list=(1, 2, 4, 8))
    json.dump(result, open(os.path.join(HERE, f"prob_model_{mode}.json"), "w"), indent=1)