#!/usr/bin/env python3
"""Scoped parameter-validity metrics for the repaired Lee--Stern protocol.

This script reads public candidate compositions and the saved finite-grid
attack assessment. It does not generate keys, execute attacks, or measure
performance. A high attack-cost upper budget is not a security lower bound.
"""
import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SEC = ROOT / "security"


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def width(bound):
    """Bytes needed for an integer in [0, bound); same convention as protocol.py."""
    return max(1, ((bound - 1).bit_length() + 7) // 8)


def natural_len(value):
    return max(1, (value.bit_length() + 7) // 8)


def is_prime(value):
    if value < 3 or value % 2 == 0:
        return False
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 2
    return True


def log2_multinomial(n, counts):
    """log2(n! / prod(counts!)). lgamma is accurate enough for these sizes."""
    return (math.lgamma(n + 1) - sum(math.lgamma(c + 1) for c in counts)) / math.log(2)


def validate_composition(cfg):
    p, n, k, weight, C = cfg["p"], cfg["n"], cfg["k"], cfg["weight"], cfg["C"]
    checks = {
        "p_is_odd_prime": is_prime(p),
        "dimensions": 0 < k < n,
        "C_length_equals_p": len(C) == p,
        "C_entries_nonnegative_integers": all(type(c) is int and c >= 0 for c in C),
        "C_sum_equals_n": sum(C) == n,
        "C_lee_weight_matches": sum(min(a, p - a) * c for a, c in enumerate(C)) == weight,
        "C_sha256_matches": hashlib.sha256(
            json.dumps(C, separators=(",", ":")).encode()
        ).hexdigest()
        == cfg["C_sha256"],
    }
    return checks, all(checks.values())


def public_field_bytes(cfg):
    p, n, k = cfg["p"], cfg["n"], cfg["k"]
    m = n - k
    bp = width(p)
    bn1 = width(n + 1)
    matrix = m * n * bp
    syndrome = m * bp
    composition = p * bn1
    amplitude = ((p + 1) // 2) * bn1
    return {
        "matrix_H_bytes": matrix,
        "syndrome_s_bytes": syndrome,
        "composition_C_bytes": composition,
        "amplitude_T_bytes": amplitude,
        "C_minus_T_bytes": composition - amplitude,
        "raw_total_bytes": matrix + syndrome + composition,
        "bytes_per_field_element": bp,
        "bytes_per_composition_count": bn1,
        "note": "raw fields only; no seed compression, wire framing, or transport overhead",
    }


def canonical_context_bytes(cfg):
    """Length of protocol.Context's canonical SHA-256 input, without building H."""
    p, n, k = cfg["p"], cfg["n"], cfg["k"]
    parts = [
        len(b"LeeStern-ExactC-v1"),  # VERSION
        natural_len(p),
        natural_len(n),
        natural_len(k),
        p * width(n + 1),  # pack(C, n+1)
        (n - k) * n * width(p),  # pack(H, p)
        (n - k) * width(p),  # pack(s, p)
    ]
    return sum(part + 8 for part in parts)


def communication_metrics(cfg, rounds):
    p, n = cfg["p"], cfg["n"]
    bp, bn = width(p), width(n)
    response_b01 = 64 + n * (bn + bp)
    response_b2 = 64 + 2 * n * bp
    payload_b01 = 96 + 1 + response_b01
    payload_b2 = 96 + 1 + response_b2
    expected = rounds * (2 * payload_b01 + payload_b2) / 3
    return {
        "rounds": rounds,
        "commitment_bytes_per_round": 96,
        "challenge_bytes_per_round": 1,
        "response_bytes_b01": response_b01,
        "response_bytes_b2": response_b2,
        "payload_bytes_b01": payload_b01,
        "payload_bytes_b2": payload_b2,
        "expected_payload_per_authentication_bytes": expected,
        "expected_payload_219_rounds_bytes": 219
        * (2 * payload_b01 + payload_b2)
        / 3,
        "note": "uniform three-challenge expectation; public setup is not repeated",
    }


def non_core_error_log2(auth_file):
    terms = auth_file["authentication"]["terms"]
    total = sum(Fraction(t["exact"]) for t in terms.values())
    return math.log2(float(total))


def evaluate():
    candidates = read_json(SEC / "inputs/candidates.json")
    decision_table = read_json(SEC / "FINAL_PARAMETER_DECISION.json")
    auth_table = read_json(SEC / "inputs/baseline_and_authentication.json")
    config = read_json(SEC / "config.json")
    verification = read_json(SEC / "final_round/VERIFICATION.json")
    saved = {x["candidate"]: x for x in decision_table["candidates"]}

    rows = []
    for cfg in candidates:
        checks, valid = validate_composition(cfg)
        p, n, k = cfg["p"], cfg["n"], cfg["k"]
        log2_orbit = log2_multinomial(n, cfg["C"])
        log2_syndrome = (n - k) * math.log2(p)
        expected_log2 = log2_orbit - log2_syndrome
        saved_row = saved[cfg["name"]]
        auth = auth_table[cfg["name"]]["authentication"]
        coverage = saved_row.get("coverage_statuses", {})
        coverage_gaps = (
            coverage.get("DOMAIN_LIMIT_OR_WORKER_FAILURE", 0)
            + coverage.get("CURRENT_LAYOUT_RESERVATION_LIMIT", 0)
        )
        rows.append(
            {
                "name": cfg["name"],
                "recipe": cfg.get("primary_recipe"),
                "p": p,
                "n": n,
                "k": k,
                "weight": cfg["weight"],
                "C_sha256": cfg["C_sha256"],
                "C_checks": checks,
                "C_valid": valid,
                "composition_shape": {
                    "max_count": max(cfg["C"]),
                    "nonzero_counts": sum(c > 0 for c in cfg["C"]),
                    "zero_counts": sum(c == 0 for c in cfg["C"]),
                },
                "witness_metrics": {
                    "log2_orbit_size": log2_orbit,
                    "log2_syndrome_space": log2_syndrome,
                    "log2_expected_witnesses_per_random_syndrome": expected_log2,
                    "bits_per_coordinate": log2_orbit / n,
                    "uniform_min_entropy_bits": log2_orbit,
                    "note": (
                        "orbit entropy is a search-space size, not a security level; "
                        "the random-syndrome expectation is a diagnostic, not the cost "
                        "of sampling the planted instance"
                    ),
                },
                "simple_generic_attack_costs_log2": {
                    "exhaustive_orbit": log2_orbit,
                    "uniform_sampling_expected_trials_upper": log2_orbit,
                    "random_syndrome_preimage_expected": log2_syndrome,
                    "grover_orbit": log2_orbit / 2,
                    "note": (
                        "generic strategies against the planted instance s=He; the "
                        "random-syndrome preimage figure is a different task and is "
                        "not the protocol attack cost"
                    ),
                },
                "public_fields_raw": public_field_bytes(cfg),
                "canonical_context_bytes": canonical_context_bytes(cfg),
                "communication": communication_metrics(cfg, 219),
                "saved_attack_assessment": {
                    "best_available_constructive_upper_log2": saved_row.get(
                        "best_available_constructive_upper_log2"
                    ),
                    "v7_best_upper_log2": saved_row.get("v7_best_upper_log2"),
                    "v7_best_heuristic_log2": saved_row.get("v7_best_heuristic_log2"),
                    "v7_declared_reserved_GiB": saved_row.get(
                        "v7_declared_reserved_GiB"
                    ),
                    "v7_best_success_lower": saved_row.get("v7_best_success_lower"),
                    "coverage_statuses": coverage,
                    "coverage_gap_count": coverage_gaps,
                    "cheap_attack_certificate_found_in_saved_rows": saved_row.get(
                        "cheap_attack_certificate_found_in_saved_rows"
                    ),
                    "general_FCSD_premise_established": saved_row.get(
                        "general_FCSD_premise_established"
                    ),
                    "supports_complete_classical128": saved_row.get(
                        "supports_complete_classical128"
                    ),
                    "security_bits": saved_row.get("security_bits"),
                    "verdict": saved_row.get("verdict"),
                },
                "auth_ledger": {
                    "rounds": auth.get("N"),
                    "hash_bits": auth.get("h"),
                    "salt_bits": auth.get("rho"),
                    "learning_sessions_cap": auth.get("learning_sessions_cap"),
                    "log2_reduction_charge": auth.get("log2_Q"),
                    "required_FCSD_budget_log2": saved_row.get(
                        "required_FCSD_budget_log2"
                    ),
                    "required_FCSD_success_upper": saved_row.get(
                        "required_FCSD_success_upper"
                    ),
                    "arithmetic_ledger_under_inherited_premise_passes": saved_row.get(
                        "arithmetic_ledger_under_inherited_premise_passes"
                    ),
                },
                "evidence_status": {
                    "necessary_orbit_size_ge_2^128": log2_orbit >= 128,
                    "simple_generic_attacks_above_2^128": log2_orbit > 128,
                    "evaluated_strategy_no_cheap_certificate": not saved_row.get(
                        "cheap_attack_certificate_found_in_saved_rows", True
                    ),
                    "general_FCSD_premise_established": False,
                    "supports_complete_classical128": False,
                },
                "verdict": "NOT_APPROVED_EVIDENCE_INSUFFICIENT",
            }
        )

    summary = {
        "candidates_total": len(rows),
        "candidates_C_valid": sum(r["C_valid"] for r in rows),
        "candidates_orbit_size_ge_2^128": sum(
            r["evidence_status"]["necessary_orbit_size_ge_2^128"] for r in rows
        ),
        "candidates_simple_generic_attacks_above_2^128": sum(
            r["evidence_status"]["simple_generic_attacks_above_2^128"] for r in rows
        ),
        "candidates_with_general_FCSD_premise": sum(
            r["evidence_status"]["general_FCSD_premise_established"] for r in rows
        ),
        "candidates_supporting_complete_classical128": sum(
            r["evidence_status"]["supports_complete_classical128"] for r in rows
        ),
        "approved_candidates": 0,
        "performance_measurement_triggered": False,
        "performance_measurement_reason": (
            "No candidate passed the step-3 security gate; step 4 is not run for a "
            "secure-parameter claim. Existing Chapter 6 timings use a different "
            "scaled composition family and are engineering evidence only."
        ),
    }

    return {
        "format": "LeeStern-ParameterAssessment-v1",
        "status": "SCOPED_ASSESSMENT_NOT_SECURITY_PROOF",
        "date": "2026-10-29",
        "target": {
            "security_bits": 128,
            "classical_only": True,
            "security_notion": (
                "single-user, strictly sequential active learning, one final "
                "impersonation attempt"
            ),
            "FCSD_budget_log2": config["required_FCSD_budget_log2"],
            "FCSD_success_upper": "1/4",
            "reduction_log2_charge": auth_table[
                next(iter(auth_table))
            ]["authentication"]["log2_Q"],
            "non_core_error_log2": non_core_error_log2(
                auth_table[next(iter(auth_table))]
            ),
            "authenticated_rounds": auth_table[
                next(iter(auth_table))
            ]["authentication"]["N"],
            "hash_bits": auth_table[next(iter(auth_table))]["authentication"]["h"],
            "salt_bits": auth_table[next(iter(auth_table))]["authentication"]["rho"],
            "learning_sessions_cap_log2": 20,
            "performance_configuration": {
                "rounds": 219,
                "hash": "SHA-256",
                "salt_bits": 256,
                "security_claim": False,
                "note": "performance settings are not a 128-bit parameter set",
            },
        },
        "sources": {
            "candidates": str(SEC / "inputs/candidates.json"),
            "attack_decision": str(SEC / "FINAL_PARAMETER_DECISION.json"),
            "authentication_ledger": str(SEC / "inputs/baseline_and_authentication.json"),
            "verification": str(SEC / "final_round/VERIFICATION.json"),
            "source_sha256": {
                "security/inputs/candidates.json": sha256(SEC / "inputs/candidates.json"),
                "security/FINAL_PARAMETER_DECISION.json": sha256(
                    SEC / "FINAL_PARAMETER_DECISION.json"
                ),
                "security/inputs/baseline_and_authentication.json": sha256(
                    SEC / "inputs/baseline_and_authentication.json"
                ),
                "security/config.json": sha256(SEC / "config.json"),
                "security/final_round/VERIFICATION.json": sha256(
                    SEC / "final_round/VERIFICATION.json"
                ),
            },
            "verification_status": verification.get("status"),
            "target_attack_executed": verification.get("target_attack_executed"),
            "verification_scope": (
                "model consistency and saved-certificate replay; not a global "
                "security proof"
            ),
        },
        "candidates": rows,
        "summary": summary,
        "claim_boundary": {
            "supported": [
                "all three fixed C pass the composition and hash checks",
                "all three have witness-orbit entropy far above 2^128",
                "simple exhaustive and random-sampling attacks are far above 2^128",
                "the saved finite-grid search found no cheap attack certificate",
                "the conditional authentication arithmetic passes only if the "
                "general FCSD premise is assumed",
            ],
            "not_supported": [
                "any of the three candidates is a concrete classical 128-bit parameter",
                "the best attack cost is at least 2^128",
                "general FCSD hardness at the required budget",
                "full impersonation security",
                "lattice, advanced Hamming/multiset, algebraic, or quantum attacks "
                "have been excluded",
                "the repair is more efficient than the original protocol",
            ],
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "results.json")
    args = parser.parse_args()
    result = evaluate()
    args.output.write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(result["summary"], ensure_ascii=False, indent=2))
    print("wrote", args.output)


if __name__ == "__main__":
    main()
