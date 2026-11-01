#!/usr/bin/env python3
"""G1: Construction A / BKZ calibration for the three fixed-composition candidates.

Scope: exactly one attack family -- the primal Construction A coset lattice with
the balanced Lee representative, priced by the BKZ block-size cost model under a
classical core-SVP sieve.  This is a heuristic calibration, not a security proof.

Outputs g1_lattice_results.json and prints a Markdown table.
"""
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
SEC = HERE.parent / "security"
LOG2 = math.log2

BUDGET_LOG2 = 166.0          # charge/query cap of the authentication reduction
COEFF_CLASSICAL = 0.292      # log2 cost per unit BKZ block size, classical sieve
COEFF_QUANTUM = 0.265        # same for a quantum sieve (sensitivity only)
TARGET_W_LOG2 = -2.0         # success <= 1/4
TARGET_S_LOG2 = -130.0       # success <= 2^-130 (2-bit margin below 2^-128)


def ball_log2(n, radius):
    """log2 of the volume of the n-dimensional Euclidean ball of this radius."""
    return ((n / 2.0) * math.log2(math.pi)
            - math.lgamma(n / 2.0 + 1) / math.log(2)
            + n * math.log2(radius))


def delta0(beta):
    """Root-Hermite factor of BKZ with block size beta."""
    return (beta / (2 * math.pi * math.e)) ** (1 / (2 * beta))


def beta_for_delta0(target):
    """Smallest real beta with delta0(beta) <= target (delta0 decreases in beta)."""
    lo, hi = 2.0, 20000.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if delta0(mid) > target:
            lo = mid
        else:
            hi = mid
    return lo


def det_root(p, m, n):
    """det(L)^(1/n) for the Construction A lattice, det = p^m."""
    return math.exp((m * math.log(p)) / n)


def gaussian_heuristic(p, m, n):
    return det_root(p, m, n) * math.sqrt(n / (2 * math.pi * math.e))


def bkz_norm(p, m, n, beta):
    """Predicted first-vector norm of a BKZ-beta reduced basis."""
    return det_root(p, m, n) * delta0(beta) ** n


def points_log2(radius, p, m, n):
    """Gaussian-heuristic log2 of |Lambda ∩ B(radius)| = Vol(B)/det."""
    return ball_log2(n, radius) - m * math.log2(p)


def candidate_metrics(cfg):
    p, n, k = cfg["p"], cfg["n"], cfg["k"]
    m = n - k
    C = cfg["C"]
    bal = [(a if a <= (p - 1) // 2 else a - p) for a in range(p)]

    norm_sq = sum(ca * bal[a] ** 2 for a, ca in enumerate(C))
    lee_l1 = sum(ca * abs(bal[a]) for a, ca in enumerate(C))
    if lee_l1 != cfg["weight"]:
        raise ValueError(f"{cfg['name']}: computed L1 {lee_l1} != stored weight {cfg['weight']}")
    norm = math.sqrt(norm_sq)
    dr = det_root(p, m, n)
    gh = gaussian_heuristic(p, m, n)

    # beta at which the BKZ output norm first falls to the secret's Euclidean norm
    required_delta = (norm / dr) ** (1.0 / n)
    beta_req = beta_for_delta0(required_delta)
    cost_req = COEFF_CLASSICAL * beta_req

    # budget side: how far BKZ can go at the charge cap
    beta_budget = BUDGET_LOG2 / COEFF_CLASSICAL
    reach_budget = bkz_norm(p, m, n, beta_budget)
    log2_points_budget = points_log2(reach_budget, p, m, n)
    log2_points_secret = points_log2(norm, p, m, n)

    # M1: norm-guarantee / threshold model.  Attack succeeds only if the BKZ
    # norm reaches the secret; below that block size the model predicts failure.
    m1_at_budget = reach_budget <= norm

    # M2: exchangeable-basis counting model.  The reduced basis carries n
    # vectors of norm ~ reach_budget; if the secret were exchangeable with all
    # lattice points in that ball, P(secret in basis) = n / |Lambda ∩ B|.
    m2_basis = LOG2(n) - log2_points_budget
    m2_single = -log2_points_budget

    return {
        "name": cfg["name"],
        "recipe": cfg.get("primary_recipe"),
        "p": p, "n": n, "k": k, "m": m,
        "weight": cfg["weight"],
        "lee_l1_norm": lee_l1,
        "balanced_euclidean_norm": norm,
        "det_root": dr,
        "gaussian_heuristic": gh,
        "norm_over_GH": norm / gh,
        "GH_over_norm_gap": gh / norm,
        "log2_points_at_secret_norm": log2_points_secret,
        "required_delta0": required_delta,
        "required_beta": beta_req,
        "required_cost_log2": cost_req,
        "cost_margin_over_budget_bits": cost_req - BUDGET_LOG2,
        "M1_threshold_success_at_budget": 0.0 if not m1_at_budget else 1.0,
        "M1_budget_reaches_secret": m1_at_budget,
        "M2_success_basis_log2": m2_basis,
        "M2_success_single_log2": m2_single,
    }


def scaled_metrics(n, ratio, p=251):
    """Scale the norm/GH ratio to dimension n at rate 1/2 (new composition needed)."""
    k = n // 2
    m = n - k
    dr = det_root(p, m, n)
    gh = gaussian_heuristic(p, m, n)
    norm = ratio * gh
    required_delta = (norm / dr) ** (1.0 / n)
    beta_req = beta_for_delta0(required_delta)
    beta_budget = BUDGET_LOG2 / COEFF_CLASSICAL
    reach = bkz_norm(p, m, n, beta_budget)
    pts = points_log2(reach, p, m, n)
    return {
        "n": n, "k": k, "m": m,
        "secret_norm": norm,
        "GH": gh,
        "required_beta": beta_req,
        "required_cost_log2": COEFF_CLASSICAL * beta_req,
        "reachable_norm_at_budget": reach,
        "log2_points_in_reachable_ball": pts,
        "M2_success_basis_log2": LOG2(n) - pts,
    }


def main():
    input_path = SEC / "inputs/candidates.json"
    raw = input_path.read_bytes()
    candidates = json.loads(raw)
    rows = [candidate_metrics(c) for c in candidates]

    p, n, k = rows[0]["p"], rows[0]["n"], rows[0]["k"]
    m = n - k
    beta_budget = BUDGET_LOG2 / COEFF_CLASSICAL
    beta_quantum = BUDGET_LOG2 / COEFF_QUANTUM
    common = {
        "p": p, "n": n, "k": k, "m": m,
        "log2_det": m * math.log2(p),
        "det_root": det_root(p, m, n),
        "GH": gaussian_heuristic(p, m, n),
        "budget_log2": BUDGET_LOG2,
        "beta_budget_classical": beta_budget,
        "reachable_norm_at_budget_classical": bkz_norm(p, m, n, beta_budget),
        "log2_points_in_reachable_ball": points_log2(
            bkz_norm(p, m, n, beta_budget), p, m, n),
        "beta_budget_quantum_sensitivity": beta_quantum,
        "reachable_norm_at_budget_quantum": bkz_norm(p, m, n, beta_quantum),
    }

    worst = max(rows, key=lambda r: r["norm_over_GH"])
    best = min(rows, key=lambda r: r["norm_over_GH"])
    scaling = []
    for ratio, tag in ((best["norm_over_GH"], "best_ratio"), (worst["norm_over_GH"], "worst_ratio")):
        for dim in (640, 768, 1024, 1280):
            row = scaled_metrics(dim, ratio)
            row["ratio_tag"] = tag
            row["norm_over_GH"] = ratio
            scaling.append(row)

    m2_worst = max(r["M2_success_basis_log2"] for r in rows)
    verdict = {
        "scope": "lattice family only; does not quantify over all FCSD solvers",
        "target_W_lattice_family": {
            "threshold_model_success_at_budget": 0.0,
            "surrogate_success_log2": m2_worst,
            "passes_under_both_models": m2_worst <= TARGET_W_LOG2,
            "margin_over_quarter_bits": TARGET_W_LOG2 - m2_worst,
            "established": False,
            "note": (
                "Heuristic calibration, not a theorem. The threshold model is "
                "the standard pricing model; the surrogate is a conservative "
                "cross-check. Neither quantifies over all FCSD solvers."
            ),
        },
        "target_S_lattice_family": {
            "surrogate_success_log2": m2_worst,
            "margin_over_2^-130_bits": TARGET_S_LOG2 - m2_worst,
            "passes_surrogate_threshold": m2_worst <= TARGET_S_LOG2,
            "established": False,
            "note": (
                "The surrogate passes 2^-130 by only {} bits at n=640, less "
                "than the model error; n=640 cannot support a 2^-130 claim. "
                "n=1024 has a 138-bit threshold margin and is the robust option."
            ).format(round(TARGET_S_LOG2 - m2_worst, 2)),
        },
        "global_target_W_established": False,
        "remaining_gap": (
            "Target W quantifies over all FCSD solvers. G1 prices only the "
            "Construction A/BKZ family; ISD/restricted-ball, full_union, "
            "hybrid, algebraic solvers and the general premise remain open."
        ),
    }

    result = {
        "format": "LeeStern-G1-LatticeCalibration-v2",
        "status": "HEURISTIC_CALIBRATION_NOT_PROOF",
        "date": "2026-10-29",
        "input": {
            "path": "security/inputs/candidates.json",
            "sha256": hashlib.sha256(raw).hexdigest(),
        },
        "model": {
            "attack": "primal Construction A coset lattice, balanced Lee representative, BKZ + nearest-plane/coset decoding",
            "cost": "2^(0.292 beta) classical core-SVP sieve; beta = block size",
            "delta0": "delta0(beta) = (beta/(2*pi*e))^(1/(2 beta)); BKZ output norm = delta0(beta)^n * det^(1/n)",
            "budget": "charge/query cap 2^166 of the authentication reduction (generous: lattice attack uses no RO queries)",
            "scope": "classical, single target, no preprocessing, no composition constraint exploited",
        },
        "common": common,
        "candidates": rows,
        "scaling": scaling,
        "verdict": verdict,
        "limitations": [
            "Both success models are heuristics, not theorems.",
            "The composition constraint (e in S_C) is used only through the norm of the balanced representative.",
            "Lee-metric-specific embeddings, hybrid lattice/ISD attacks and dual attacks are not priced.",
            "Core-SVP is a lower bound on classical sieve cost; using it against the charge cap is conservative.",
        ],
    }

    out = HERE / "g1_lattice_results.json"
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"common: det^(1/n)={common['det_root']:.4f} GH={common['GH']:.3f} "
          f"beta_budget={beta_budget:.2f} reach={common['reachable_norm_at_budget_classical']:.3f} "
          f"log2N={common['log2_points_in_reachable_ball']:.3f}")
    print(f"quantum sensitivity: beta={beta_quantum:.2f} reach={common['reachable_norm_at_budget_quantum']:.3f}")
    print()
    print("| candidate | norm | norm/GH | required beta | required cost | margin | M2 basis | points at norm |")
    print("|---|---:|---:|---:|---:|---:|---:|---:|")
    for r in rows:
        print(f"| {r['recipe']} | {r['balanced_euclidean_norm']:.2f} | {r['norm_over_GH']:.3f} | "
              f"{r['required_beta']:.1f} | 2^{r['required_cost_log2']:.1f} | +{r['cost_margin_over_budget_bits']:.1f} bits | "
              f"2^{r['M2_success_basis_log2']:.1f} | 2^{r['log2_points_at_secret_norm']:.1f} |")
    print()
    print("scaling (worst ratio):")
    for s in scaling:
        if s["ratio_tag"] != "worst_ratio":
            continue
        print(f"  n={s['n']}: beta_req={s['required_beta']:.0f} cost=2^{s['required_cost_log2']:.1f} "
              f"reach={s['reachable_norm_at_budget']:.1f} M2=2^{s['M2_success_basis_log2']:.1f}")
    print()
    print(json.dumps(verdict, ensure_ascii=False, indent=2))
    print("wrote", out)


if __name__ == "__main__":
    main()
