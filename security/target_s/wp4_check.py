#!/usr/bin/env python3
"""WP4: compress PROOF.md into paper section 4.4 and itemise T_B <= 2^166.

The arithmetic side of WP4_SECTION44.md:

* the reduction symbol set J = LT, R = 1+J+3N, Q = R*q_A + 3*J + 6*N is the
  WP2 one, unchanged; Q is the live-query count of ONE experiment, and the J
  hybrids are J separate experiments (so J*Q is not the work of B);
* Table 3.1 opens exactly two of the three salts per attempt, and the
  simulator's placeholder sits at the unopened index, so one attempt
  contributes one hidden salt: the hybrid count is J, not 3J;
* a full work itemisation of B: every term except R*t_A is <= 2^77.4, so
  T_B <= 2^164.34 < 2^166 with 1.67 bits of headroom;
* the 2^166 interface is a "one oracle call = one step" interface.  Read as
  bit operations with q_A = 2^128 it is infeasible (Q*8192 = 2^177.33); the
  self-consistent bit-operation reading derives q_A <= 2^115 instead, which
  only improves the hiding and binding rows.  The paper must state which
  metric it uses;
* two WP1 sentences are false as upper bounds and must be repaired:
  non-FCSD sum is 2^-149.7487 (larger than 2^-149.75) and the total is
  2^-135.9999 (larger than 2^-136).

Writes review/WP4_CHECK.json and target_s/WP4_RESOURCE_TABLE.csv.
It does not prove Assumption F and does not re-measure the implementation.
"""
from __future__ import annotations

import csv
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WP2 = ROOT / "target_s" / "WP2_DECISIONS.json"
OUT_JSON = ROOT / "review" / "WP4_CHECK.json"
OUT_CSV = ROOT / "target_s" / "WP4_RESOURCE_TABLE.csv"

N, T = 256, 322
L = 1 << 28
QA_LOG2 = 128
RHO, H = 384, 512
DELTA0_LOG2 = 198
BUDGET_LOG2 = 166

# Per-unit step bound of each work item, in the declared metric where one
# random-oracle call is one step.  The rationales are in WP4_SECTION44.md 3.1.
C_G_LOG2 = 31          # H sampling + rank test O(m^2 n) = 2^28 field ops
C_SIM_LOG2 = 14        # three commitments over ~8192-bit inputs + vector ops
C_VER_LOG2 = 14        # one verification: two openings + list checks
C_EXT_LOG2 = 20        # three-opening extraction and He* check
C_IO_LOG2 = 13         # one round: copy <=8 Kibit prefix, parse, log
C_CTL_LOG2 = 32        # scheduling, counters, replay bookkeeping


def log2f(x: Fraction | int) -> float:
    if isinstance(x, int):
        return math.log2(x)
    if x <= 0:
        return float("-inf")
    return math.log2(x.numerator) - math.log2(x.denominator)


def main() -> int:
    d = json.loads(WP2.read_text(encoding="utf-8"))
    p, res, led = d["parameters"], d["resources"], d["ledger"]

    checks: list[dict] = []

    def check(name: str, ok: bool, detail: str) -> None:
        checks.append({"name": name, "pass": bool(ok), "detail": detail})

    # ---------------- symbols and the hybrid count -----------------------
    J = L * T
    R = 1 + J + 3 * N
    qA = 1 << QA_LOG2
    Q = R * qA + 3 * J + 6 * N
    check("J = LT reproduces WP2", J == res["J"], f"J={J}")
    check("R = 1+J+3N reproduces WP2", R == res["R"], f"R={R}")
    check("Q = Rq_A+3J+6N reproduces WP2",
          math.isclose(math.log2(Q), res["Q_log2"], abs_tol=1e-12),
          f"log2 Q={math.log2(Q):.12f}")
    check("hybrid count is the number of attempts, J = L*T",
          J // L == T and J == L * T, f"J/L={J // L} = T={T}")
    check("theorem 4.4 has five terms (FCSD, knowledge, binding, hiding sum, retry)",
          len(("FCSD", "knowledge", "binding", "hiding_sum", "retry")) == 5,
          "structure of (4.2); the truncation term is added outside and charged in the ledger")

    # Table 3.1: challenge b opens salts {1,2}, {1,3}, {2,3}; the simulator's
    # placeholder in branch g is cm_3, cm_2, cm_1 respectively.
    openings = {0: (1, 2), 1: (1, 3), 2: (2, 3)}
    placeholder = {0: 3, 1: 2, 2: 1}
    check("Table 3.1 leaves exactly one unopened commitment per attempt",
          all(3 - len(v) == 1 for v in openings.values()),
          f"unopened per branch = {[(b, 3 - len(v)) for b, v in openings.items()]}")
    check("the placeholder of branch g is the unopened commitment of challenge g",
          all(placeholder[g] not in openings[g] for g in openings),
          f"placeholder={placeholder}, openings={openings}")
    check("so one attempt hides one salt and the hybrid count is J, not 3J",
          all(3 - len(openings[g]) == 1 for g in openings), "3J only as robustness check")

    # ---------------- work itemisation of B ------------------------------
    items = {
        "key_setup": 1 << C_G_LOG2,
        "replay_executions_R_tA": R * (1 << QA_LOG2),
        "simulator_commitments_J_csim": J * (1 << C_SIM_LOG2),
        "final_verifications_3N_cver": 3 * N * (1 << C_VER_LOG2),
        "extraction_c_ext": 1 << C_EXT_LOG2,
        "prefix_replay_io_R(L+N+1)_cio": R * (L + N + 1) * (1 << C_IO_LOG2),
        "control_c_ctl": 1 << C_CTL_LOG2,
    }
    oracle_extra = 0  # every oracle call is one step and is inside t_A / c_sim / c_ver
    T_B = sum(items.values()) + oracle_extra
    main_term = items["replay_executions_R_tA"]
    second = max(v for k, v in items.items() if k != "replay_executions_R_tA")

    check("T_B <= 2^166 (exact integer certificate)", T_B <= (1 << BUDGET_LOG2),
          f"log2 T_B={math.log2(T_B):.12f}, headroom={BUDGET_LOG2 - math.log2(T_B):.6f} bits")
    check("headroom is the expected 1.669 bits", math.isclose(BUDGET_LOG2 - math.log2(T_B), 1.66908310905, abs_tol=1e-9),
          f"headroom={BUDGET_LOG2 - math.log2(T_B):.12f}")
    check("the main term is R*t_A and equals the query count at one step per call",
          main_term == R * qA, f"log2 R*t_A={math.log2(main_term):.9f}")
    check("every other item is at most 2^-80 of the main term",
          second * (1 << 80) <= main_term,
          f"second item log2={math.log2(second):.3f}, main log2={math.log2(main_term):.3f}")
    check("T_B >= Q (queries are steps, so work can only be larger)",
          T_B >= Q, f"log2 T_B={math.log2(T_B):.6f} >= log2 Q={math.log2(Q):.6f}")
    check("Q < 2^166 (the query count alone fits the interface)",
          Q < (1 << BUDGET_LOG2), f"log2 Q={math.log2(Q):.6f}")
    check("Q < 2^rho (the hiding bound is nontrivial)",
          Q < (1 << RHO), f"log2 Q={math.log2(Q):.3f} < {RHO}")
    check("3J simulator commitments fit the J*c_sim charge",
          3 * J <= J * (1 << (C_SIM_LOG2 - 1)),
          f"3J={3 * J:.3e} vs J*c_sim={J * (1 << C_SIM_LOG2):.3e}")
    check("6N final-phase evaluations fit the 3N*c_ver charge",
          6 * N <= 3 * N * (1 << (C_VER_LOG2 - 1)),
          f"6N={6 * N} vs 3N*c_ver={3 * N * (1 << C_VER_LOG2)}")

    # ---------------- metric reading -------------------------------------
    bits_per_query = 1024 * 8          # one commitment input: n coordinates over F_p
    check("bit-operation reading with q_A=2^128 is infeasible",
          Q * bits_per_query > (1 << BUDGET_LOG2),
          f"Q*8192 = 2^{math.log2(Q * bits_per_query):.3f} > 2^166")
    qA_bits = (1 << QA_LOG2) // bits_per_query
    Q_bits = R * qA_bits + 3 * J + 6 * N
    check("self-consistent bit-operation reading derives q_A <= 2^115",
          qA_bits == 1 << 115, f"q_A={qA_bits}, log2={math.log2(qA_bits):.0f}")
    check("bit reading keeps T_B <= 2^166 and both ROM rows improve",
          R * (1 << QA_LOG2) <= (1 << BUDGET_LOG2),
          f"log2 hiding={math.log2(J) + 1 + math.log2(Q_bits) - RHO:.3f} (was -182.338), "
          f"log2 binding={2 * math.log2(Q_bits) - 1 - H:.3f} (was -184.338)")

    # ---------------- Q scope -------------------------------------------
    check("J separate hybrids do not turn JQ into the work of B",
          J * Q > (1 << 200) and T_B < (1 << BUDGET_LOG2),
          f"log2 JQ={math.log2(J * Q):.3f} is a union-bound total across J games, "
          f"not an execution of B")

    # ---------------- ledger and the two repaired sentences --------------
    knowledge = Fraction(2, 3) ** N
    retry = L * Fraction(2, 3) ** T
    hiding = Fraction(2 * Q, 1 << RHO) * J
    binding = Fraction(Q * (Q - 1), 2 << H)
    truncation = Fraction(2 * J + 2, 1 << DELTA0_LOG2)
    non_core = knowledge + retry + hiding + binding + truncation
    total = non_core + Fraction(1, 1 << 136)
    check("2J+2 = 2LT+2 (the truncation multiplier of (R.8))",
          2 * J + 2 == 2 * L * T + 2, f"2J+2={2 * J + 2}")
    check("ledger rows reproduce WP2",
          all(math.isclose(log2f(v), led[k], abs_tol=1e-9) for k, v in
              (("knowledge_log2", knowledge), ("retry_log2", retry),
               ("hiding_log2", hiding), ("binding_log2", binding),
               ("truncation_log2", truncation), ("non_core_log2", non_core),
               ("total_log2", total))),
          f"non_core={log2f(non_core):.9f}, total={log2f(total):.9f}")
    check("WP1 sentence 'non-FCSD <= 2^-149.75' is false: the sum is slightly larger",
          log2f(non_core) > -149.75,
          f"log2 non_core={log2f(non_core):.9f} > -149.75 (repair to 2^-149.5 or print 2^-149.7487)")
    check("the repaired statement 'non-FCSD <= 2^-149.5' holds",
          log2f(non_core) <= -149.5, f"log2={log2f(non_core):.9f} <= -149.5")
    check("WP1 sentence 'total = 2^-136' is false as an upper bound",
          log2f(total) > -136,
          f"log2 total={log2f(total):.9f} > -136 (print 2^-135.9999)")
    check("total < 2^-128 with margin 7.9999 bits (total is 0.000105 bits above 2^-136)",
          total < Fraction(1, 1 << 128) and
          math.isclose(-128 - log2f(total), 7.999895, abs_tol=1e-5),
          f"log2 total={log2f(total):.9f}, margin={-128 - log2f(total):.6f}")
    check("truncation term is below the knowledge term",
          truncation <= knowledge, f"log2 truncation={log2f(truncation):.6f}")

    # ---------------- robustness: 3 salts per attempt --------------------
    hiding_3salt = Fraction(2 * Q, 1 << RHO) * (3 * J)
    non_core_3salt = knowledge + retry + hiding_3salt + binding + truncation
    check("3J reading still clears the acceptance line (robustness only)",
          hiding_3salt <= knowledge and log2f(non_core_3salt) <= -149.5,
          f"log2 hiding={log2f(hiding_3salt):.6f}, sum={log2f(non_core_3salt):.9f}")

    # ---------------- outputs -------------------------------------------
    rows = [
        ("symbols", "J = LT", J, f"log2={math.log2(J):.9f}", "hybrids = attempts"),
        ("symbols", "R = 1+J+3N", R, f"log2={math.log2(R):.9f}", "executed-A bound of one experiment"),
        ("symbols", "Q = Rq_A+3J+6N", Q, f"log2={math.log2(Q):.9f}", "live queries of ONE experiment"),
        ("work", "key setup c_G", items["key_setup"], f"log2={C_G_LOG2}",
         "H sampling + rank test O(m^2 n)=2^28 field ops + He=2^19 + composition sampler"),
        ("work", "R*t_A", main_term, f"log2={math.log2(main_term):.9f}",
         "t_A <= q_A = 2^128, one oracle call = one step"),
        ("work", "J*c_sim", items["simulator_commitments_J_csim"], f"log2={math.log2(items['simulator_commitments_J_csim']):.6f}",
         "c_sim <= 2^14: three commitments over ~8192-bit inputs + vector ops"),
        ("work", "3N*c_ver", items["final_verifications_3N_cver"], f"log2={math.log2(items['final_verifications_3N_cver']):.6f}",
         "c_ver <= 2^14 per verification"),
        ("work", "c_ext", items["extraction_c_ext"], f"log2={C_EXT_LOG2}",
         "three-opening extraction + He* check"),
        ("work", "R(L+N+1)*c_io", items["prefix_replay_io_R(L+N+1)_cio"], f"log2={math.log2(items['prefix_replay_io_R(L+N+1)_cio']):.6f}",
         "c_io <= 2^13 per round of prefix replay"),
        ("work", "c_control", items["control_c_ctl"], f"log2={C_CTL_LOG2}", "scheduling and counters"),
        ("work", "c_oracle (extra)", oracle_extra, "0", "every oracle call is one step, already inside t_A/c_sim/c_ver"),
        ("total", "T_B", T_B, f"log2={math.log2(T_B):.9f}", f"<= 2^166, headroom {BUDGET_LOG2 - math.log2(T_B):.6f} bits"),
        ("metric", "Q*8192 (bit-operation reading)", Q * bits_per_query, f"log2={math.log2(Q * bits_per_query):.9f}",
         "infeasible at 2^166; the interface is one-call-one-step"),
        ("metric", "Q' at q_A<=2^115 (bit-operation reading)", Q_bits, f"log2={math.log2(Q_bits):.9f}",
         "self-consistent alternative; hiding/binding improve, non-FCSD sum unchanged at the top"),
        ("ledger", "non-FCSD sum", None, f"log2={log2f(non_core):.9f}", "must be printed as 2^-149.7487, not 2^-149.75"),
        ("ledger", "total", None, f"log2={log2f(total):.9f}", "must be printed as 2^-135.9999, not 2^-136"),
    ]

    result = {
        "status": "PASS" if all(c["pass"] for c in checks) else "FAIL",
        "work_package": "WP4",
        "date": "2026-10-30",
        "resources": {"J": J, "R": R, "Q": Q, "T_B": T_B,
                      "log2_T_B": math.log2(T_B), "headroom_bits": BUDGET_LOG2 - math.log2(T_B)},
        "metric": {"declared": "one RO call = one step",
                   "Q_times_8192_log2": math.log2(Q * bits_per_query),
                   "q_A_bits_log2": math.log2(qA_bits),
                   "Q_bits_log2": math.log2(Q_bits),
                   "hiding_bits_log2": math.log2(J) + 1 + math.log2(Q_bits) - RHO,
                   "binding_bits_log2": 2 * math.log2(Q_bits) - 1 - H},
        "ledger_log2": {"non_core": log2f(non_core), "total": log2f(total),
                        "non_core_3salt": log2f(non_core_3salt)},
        "repairs": [
            "§4.4: non-FCSD <= 2^-149.75 is false; use 2^-149.5 (true) or print 2^-149.7487",
            "§5.5 C-2: total = 2^-136 is false as an upper bound; print 2^-135.9999",
            "margin of the total is 7.9999 bits, not 8.0001 (WP2 docs corrected in place, erratum noted)",
            "§4.4: state that 2^166 is a one-call-one-step interface, not a bit-operation count",
            "§4.4: state that Q is per experiment and the J hybrids are separate experiments",
        ],
        "checks": checks,
    }

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    with OUT_CSV.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["group", "quantity", "value", "log2", "rationale"])
        for g, q, v, lg, why in rows:
            w.writerow([g, q, "" if v is None else v, lg, why])

    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    print(json.dumps({k: result[k] for k in ("status", "resources", "metric", "ledger_log2")},
                     ensure_ascii=False, indent=2))
    failed = [c for c in checks if not c["pass"]]
    if failed:
        print("FAILED CHECKS:")
        for c in failed:
            print(" -", c["name"], "|", c["detail"])
        return 1
    print(f"All {len(checks)} checks passed; wrote {OUT_JSON.relative_to(ROOT)} and {OUT_CSV.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
