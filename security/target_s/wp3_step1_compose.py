# -*- coding: utf-8 -*-
"""WP3 step 1: freeze the composition C at n=1024 and settle the |S_C| vs p^k regime.

Inputs : ../security/inputs/candidates.json  (the three n=640 frozen compositions)
Outputs: review/WP3_STEP1.json, WP3_STEP1_TABLE.csv

Every sign decision (|S_C| >= p^k ?) uses exact integers; floats are for reporting only.
Counts are indexed by the value a in Z_p; the balanced representative is
b(a) = a for a <= 125, a - p otherwise, so |b(a)| = min(a, p-a).
"""
from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CAND = HERE.parent / "security" / "inputs" / "candidates.json"
OUT_JSON = HERE.parent / "review" / "WP3_STEP1.json"
OUT_CSV = HERE / "WP3_STEP1_TABLE.csv"

P = 251
N0, K0 = 640, 320
N, K = 1024, 512
LOGE = math.log(2.0)
PK = P ** K
FACT_N = math.factorial(N)


def lee_of(a: int) -> int:
    return min(a, P - a)


def log2_fact(m: int) -> float:
    return math.lgamma(m + 1) / LOGE


def lee_weight(counts: list[int]) -> int:
    return sum(c * lee_of(a) for a, c in enumerate(counts))


def l2sq(counts: list[int]) -> int:
    return sum(c * lee_of(a) ** 2 for a, c in enumerate(counts))


def log2_size(counts: list[int], n: int) -> float:
    return log2_fact(n) - sum(log2_fact(c) for c in counts if c)


def exact_ge_pk(counts: list[int], n: int = N, k: int = K) -> bool:
    den = P ** k
    for c in counts:
        if c:
            den *= math.factorial(c)
    return math.factorial(n) >= den


def scale_counts(counts: list[int], n_from: int, n_to: int) -> list[int]:
    ratio = n_to / n_from
    base = [int(math.floor(c * ratio)) for c in counts]
    rem = n_to - sum(base)
    order = sorted(range(len(counts)),
                   key=lambda i: -(counts[i] * ratio - base[i]))
    j = 0
    while rem > 0:
        base[order[j % len(order)]] += 1
        rem -= 1
        j += 1
    assert sum(base) == n_to
    return base


def family_geometric(r: float, n: int = N) -> list[int]:
    """counts ~ r^{|b(a)|} over all p values, largest remainder."""
    w = [r ** lee_of(a) for a in range(P)]
    tot = sum(w)
    target = [x / tot * n for x in w]
    base = [int(math.floor(x)) for x in target]
    rem = n - sum(base)
    order = sorted(range(P), key=lambda i: -(target[i] - base[i]))
    for j in range(rem):
        base[order[j % P]] += 1
    return base


def family_uniform_t(t: int, n: int = N) -> list[int]:
    """t equally likely Lee values 0, +-1, +-2, ...; counts differ by at most 1."""
    vals = [0]
    a = 1
    while len(vals) < t:
        vals.append(a)
        if len(vals) < t:
            vals.append(P - a)
        a += 1
    q, r = divmod(n, t)
    counts = [0] * P
    for i, v in enumerate(vals):
        counts[v] = q + (1 if i < r else 0)
    return counts


def describe(name: str, counts: list[int], n: int = N, k: int = K) -> dict:
    lg = log2_size(counts, n)
    lp = k * math.log2(P)
    return {
        "name": name,
        "zeros": counts[0],
        "zero_fraction": round(counts[0] / n, 6),
        "distinct_values": sum(1 for c in counts if c),
        "lee_weight": lee_weight(counts),
        "l2": math.sqrt(l2sq(counts)),
        "log2_size": round(lg, 4),
        "log2_pk": round(lp, 4),
        "log2_multiplicity": round(lg - lp, 4),
        "size_ge_pk": exact_ge_pk(counts, n, k),
        "t_eff": round(2 ** (lg / n), 3),
    }


def crosses(counts: list[int]) -> bool:
    den = PK
    for c in counts:
        if c:
            den *= math.factorial(c)
    return FACT_N >= den


def delta0(beta: float) -> float:
    return (beta / (2 * math.pi * math.e)) ** (1 / (2 * beta))


def bkz_beta(target_norm: float, n: int, det_root: float) -> float:
    """Smallest beta with delta0(beta)^n * det^(1/n) <= target_norm (G1 model)."""
    lo, hi = 2.0, 5000.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if delta0(mid) ** n * det_root >= target_norm:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")

    cands = json.loads(CAND.read_text(encoding="utf-8"))
    gh = P ** 0.5 * math.sqrt(N / (2 * math.pi * math.e))
    gh0 = P ** 0.5 * math.sqrt(N0 / (2 * math.pi * math.e))

    rows: list[dict] = []
    checks: list[dict] = []

    def chk(label: str, ok: bool, detail: str) -> None:
        checks.append({"check": label, "pass": bool(ok), "detail": detail})

    # ---- 1. replay the three n=640 candidates --------------------------------
    replay = []
    for c in cands:
        full = c["C"]
        assert len(full) == P and sum(full) == c["n"]
        item = describe(c["name"] + "@640", full, N0, K0)
        item["weight_declared"] = c["weight"]
        item["l2_over_GH"] = round(item["l2"] / gh0, 6)
        replay.append(item)
        chk(f"Lee weight matches candidates.json [{c['name'][:20]}]",
            item["lee_weight"] == c["weight"],
            f"{item['lee_weight']} == {c['weight']}")

    # ---- 2. scale to n=1024 --------------------------------------------------
    scaled = []
    for c in cands:
        h = scale_counts(c["C"], N0, N)
        d = describe(c["name"] + "@1024", h)
        d["l2_over_GH"] = round(d["l2"] / gh, 6)
        scaled.append(d)

    # ---- 3. G1 §6 entry point: BKZ cost of the scaled compositions ------------
    det_root = P ** 0.5
    for r, d in zip(replay, scaled):
        b0 = bkz_beta(r["l2"], N0, det_root)
        d["beta_ref_n640"] = round(b0, 1)
        d["log2_cost_n640"] = round(0.292 * b0, 2)
        b1 = bkz_beta(d["l2"], N, det_root)
        d["beta_n1024"] = round(b1, 1)
        d["log2_cost_n1024"] = round(0.292 * b1, 2)

    # ---- 4. balanced reference and frontiers ---------------------------------
    uni = family_uniform_t(P, N)
    d_uni = describe("uniform_all_251", uni)
    d_uni["l2_over_GH"] = round(d_uni["l2"] / gh, 6)

    lo, hi = 2, P
    while lo < hi:
        mid = (lo + hi) // 2
        if crosses(family_uniform_t(mid)):
            hi = mid
        else:
            lo = mid + 1
    t_star = lo
    d_tstar = describe(f"uniform_t{t_star}", family_uniform_t(t_star))
    d_tstar["l2_over_GH"] = round(d_tstar["l2"] / gh, 6)
    d_prev = describe(f"uniform_t{t_star - 1}", family_uniform_t(t_star - 1))
    d_prev["l2_over_GH"] = round(d_prev["l2"] / gh, 6)

    # geometric family: |S_C| is unimodal in r (max at r=1), so the sparse-side
    # frontier is the crossing on (0, 1]; above r=1 the family is anti-sparse.
    glo, ghi = 0.05, 1.0
    for _ in range(90):
        mid = (glo + ghi) / 2
        if crosses(family_geometric(mid)):
            ghi = mid
        else:
            glo = mid
    r_star = (glo + ghi) / 2
    d_rstar = describe(f"geometric_r{r_star:.4f}", family_geometric(r_star))
    d_rstar["l2_over_GH"] = round(d_rstar["l2"] / gh, 6)

    # ---- 5. checks -----------------------------------------------------------
    chk("BKZ model replays G1's n=640 betas 813.7 / 665.6 / 653.7",
        all(abs(d["beta_ref_n640"] - v) < 0.5
            for d, v in zip(scaled, (813.7, 665.6, 653.7))),
        str([d["beta_ref_n640"] for d in scaled]))
    chk("BKZ model replays G1's n=640 costs 237.6 / 194.4 / 190.9 bits",
        all(abs(d["log2_cost_n640"] - v) < 0.15
            for d, v in zip(scaled, (237.6, 194.4, 190.9))),
        str([d["log2_cost_n640"] for d in scaled]))
    chk("n=640 replay reproduces log2|S_C| inside [2147.2, 2408.1]",
        all(2147.2 <= r["log2_size"] <= 2408.1 for r in replay),
        str([r["log2_size"] for r in replay]))
    chk("n=640 replay reproduces G1 l2 norms 72.40 / 92.17 / 94.35",
        all(abs(r["l2"] - v) < 0.01 for r, v in zip(replay, (72.4017, 92.1683, 94.3451))),
        str([r["l2"] for r in replay]))
    chk("all three scaled candidates stay below p^k (unique-preimage regime)",
        all(not d["size_ge_pk"] for d in scaled),
        str([d["log2_multiplicity"] for d in scaled]))
    chk("scaling preserves the l2/GH ratio (this is exactly G1 6's assumption)",
        all(abs(d["l2_over_GH"] - r["l2_over_GH"]) < 5e-3
            for d, r in zip(scaled, replay)),
        str([d["l2_over_GH"] for d in scaled]))
    chk("balanced composition is above p^k (many-witness regime)",
        d_uni["size_ge_pk"] and d_uni["log2_multiplicity"] > 0,
        str(d_uni["log2_multiplicity"]))
    chk(f"t*={t_star} crosses, t*={t_star - 1} does not",
        crosses(family_uniform_t(t_star)) and not crosses(family_uniform_t(t_star - 1)),
        f"{d_tstar['log2_size']} vs {d_prev['log2_size']} (p^k={round(K * math.log2(P), 3)})")
    chk("scaled geometric:2 keeps l2/GH ~ 0.7465 (G1 table row)",
        abs(scaled[0]["l2_over_GH"] - 0.7465) < 5e-3, str(scaled[0]["l2_over_GH"]))
    chk("every scaled candidate keeps Lee weight below 3n (sparse family)",
        all(d["lee_weight"] < 3 * N for d in scaled),
        str([d["lee_weight"] for d in scaled]))
    chk("uniform composition has Lee weight above 0.45*n*(p-1)/2",
        d_uni["lee_weight"] > 0.45 * N * (P - 1) / 2, str(d_uni["lee_weight"]))
    chk("the three n=640 candidates were already below p^k (regime is inherited)",
        all(r["log2_multiplicity"] < 0 for r in replay),
        str([r["log2_multiplicity"] for r in replay]))
    chk("at the frontier the secret is no longer short: l2/GH > 1",
        d_tstar["l2_over_GH"] > 1.0 and d_prev["l2_over_GH"] > 1.0,
        f"t=16: {d_prev['l2_over_GH']}, t=17: {d_tstar['l2_over_GH']}")
    chk("geometric r* is on the sparse side (0, 1] and brackets the crossing",
        crosses(family_geometric(r_star)) and not crosses(family_geometric(r_star - 0.002)),
        str(round(r_star, 5)))

    frozen_counts = scale_counts(cands[0]["C"], N0, N)
    freeze = {
        "name": "p251_n1024_k512_geometric2_scaled",
        "status": "frozen_pending_step2_coverage",
        "p": P, "n": N, "k": K,
        "derivation": "counts = geometric:2 n=640 vector scaled by 1024/640 = 1.6, largest remainder",
        "source_candidate": cands[0]["name"],
        "source_sha256": cands[0].get("C_sha256"),
        "counts": frozen_counts,
        "sum": sum(frozen_counts),
        "lee_weight": lee_weight(frozen_counts),
        "l2": round(math.sqrt(l2sq(frozen_counts)), 6),
        "log2_size": round(log2_size(frozen_counts, N), 6),
        "log2_pk": round(K * math.log2(P), 6),
        "log2_multiplicity": round(log2_size(frozen_counts, N) - K * math.log2(P), 6),
    }
    (HERE / "WP3_C_N1024.json").write_text(
        json.dumps(freeze, ensure_ascii=False, indent=2), encoding="utf-8")
    chk("frozen composition sums to n=1024 over p=251 values",
        sum(frozen_counts) == N and len(frozen_counts) == P,
        f"{sum(frozen_counts)}")
    chk("frozen composition is the scaled geometric:2 shape (zeros=242, Lee=2025)",
        frozen_counts[0] == 242 and lee_weight(frozen_counts) == 2025,
        f"zeros={frozen_counts[0]}, Lee={lee_weight(frozen_counts)}")

    status = "PASS" if all(c["pass"] for c in checks) else "FAIL"
    result = {
        "status": status,
        "work_package": "WP3 step 1",
        "date": "2026-10-30",
        "parameters": {"p": P, "n": N, "k": K, "log2_pk": round(K * math.log2(P), 4),
                       "GH_n1024": round(gh, 4), "GH_n640": round(gh0, 4)},
        "replay_n640": replay,
        "scaled_n1024": scaled,
        "balanced_reference": d_uni,
        "frontier": {
            "t_star_uniform": t_star,
            "at_t_star": d_tstar,
            "below_t_star": d_prev,
            "geometric_r_star": round(r_star, 5),
            "at_geometric_r_star": d_rstar,
        },
        "checks": checks,
    }
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")

    with OUT_CSV.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["composition", "zeros", "zero_fraction", "distinct_values",
                    "lee_weight", "l2", "l2_over_GH", "log2_size", "log2_pk",
                    "log2_multiplicity", "size_ge_pk", "t_eff"])
        for name, d in ([(r["name"], r) for r in replay]
                        + [(d["name"], d) for d in scaled]
                        + [(d["name"], d) for d in (d_prev, d_tstar, d_rstar, d_uni)]):
            w.writerow([name, d["zeros"], d["zero_fraction"], d["distinct_values"],
                        d["lee_weight"], round(d["l2"], 4), d.get("l2_over_GH"),
                        d["log2_size"], d["log2_pk"], d["log2_multiplicity"],
                        d["size_ge_pk"], d["t_eff"]])

    print(json.dumps({"status": status,
                      "scaled": [{k: d[k] for k in ("name", "zeros", "lee_weight",
                                                    "l2_over_GH", "log2_size",
                                                    "log2_multiplicity", "size_ge_pk",
                                                    "t_eff")} for d in scaled],
                      "uniform": {k: d_uni[k] for k in ("zeros", "lee_weight",
                                                        "l2_over_GH", "log2_size",
                                                        "log2_multiplicity", "t_eff")},
                      "t_star": t_star, "geometric_r_star": round(r_star, 5)},
                     ensure_ascii=False, indent=2))
    print(f"{sum(c['pass'] for c in checks)}/{len(checks)} checks pass")


if __name__ == "__main__":
    main()
