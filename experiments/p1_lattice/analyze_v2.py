"""Compare observed success block sizes with two predictions.

(1) 2016 criterion on the full embedding (d = 2k+1, m = k samples):
    sqrt(b) * sigma <= delta(b)^(2b - d - 1) * q^(k/d)        -> beta_2016
(2) Chen-Nguyen simulator (fpylll bkz_simulator) on each trial's actual LLL profile,
    progressive block sizes 10, 12, ... with the same number of tours as the experiment;
    success when the projection of the target onto the last b Gram-Schmidt directions
    (expected length sigma_t * sqrt(b), sigma_t^2 = ||target||^2 / d) is shorter than
    the simulated b_{d-b}^*.                                  -> beta_sim
Both use the same success notion as lattice-estimator's usvp (2016 model).
"""
import json, math, os, statistics as st
from fpylll import BKZ
from fpylll.tools.bkz_simulator import simulate

HERE = os.path.dirname(os.path.abspath(__file__))
Q = 251

def log_delta(b):
    return math.log(((b / (2 * math.pi * math.e)) * (math.pi * b) ** (1 / b)) ** (1 / (2 * (b - 1))))

def beta_2016(k, sigma):
    d = 2 * k + 1
    for b in range(40, d):
        if math.log(math.sqrt(b) * sigma) <= (2 * b - d - 1) * log_delta(b) + k / d * math.log(Q):
            return b
    return None

def beta_sim(prof_log2, norm2, tours):
    d = len(prof_log2); r = [2.0 ** (2 * x) for x in prof_log2]
    sig2 = norm2 / d
    for b in range(10, d + 1, 2):
        r, _ = simulate(r, BKZ.Param(block_size=b, max_loops=tours))
        if sig2 * b <= r[d - b]:
            return b
    return None

rows = [json.loads(l) for l in open(os.path.join(HERE, "bkz_scale_v2.jsonl"))]
EST = {r["k"]: r for r in json.load(open(os.path.join(HERE, "est_small.json")))}
by_k = {}
for t in rows:
    t["beta_sim"] = beta_sim(t["lll_profile"], t["norm2"], t.get("tours_per_beta", 2))
    by_k.setdefault(t["k"], []).append(t)
summary = []
for k in sorted(by_k):
    ts = by_k[k]; obs = [t["beta"] for t in ts if t["beta"]]
    b16 = beta_2016(k, ts[0]["sigma"]); sims = [t["beta_sim"] for t in ts if t["beta_sim"]]
    row = dict(k=k, d=2 * k + 1, trials=len(ts), solved=len(obs), beta_obs=sorted(obs),
               obs_median=st.median(obs) if obs else None, beta_2016=b16, beta_est_gsa=EST[k]["beta_gsa"], beta_est_cn11=EST[k]["beta_cn11"],
               gap_est=(EST[k]["beta_gsa"] - st.median(obs)) if obs else None, ratio_est=(st.median(obs) / EST[k]["beta_gsa"]) if obs else None,
               sim_median=st.median(sims) if sims else None,
               gap_2016=(b16 - st.median(obs)) if obs else None,
               gap_sim=(st.median(sims) - st.median(obs)) if obs and sims else None,
               time_median=st.median(t["time"] for t in ts))
    summary.append(row); print(json.dumps(row))
json.dump(dict(trials=[{k: v for k, v in t.items() if k != "lll_profile"} for t in rows], summary=summary),
          open(os.path.join(HERE, "bkz_scale_v2.json"), "w"), indent=1)