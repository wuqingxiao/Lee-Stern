#!/usr/bin/env python3
"""WP2: exact recomputation of the Target-S ROM ledger.

See WP2_LEDGER.md for the derivations.  This script is the arithmetic side:

* Lemma T.2: one experiment runs A at most R = 1 + J + 3N times, so the number
  of live RO queries against a single fixed random function is
  Q = R*q_A + 3*J + 6*N -- the (S.2) quantity.  Rewinds multiply live queries;
  the per-branch cap q_A does not bound them.
* Lemma T.1: eps_H <= 2*Q / 2^rho per hybrid, J hybrids.  The WP1 candidate
  tightening (Q -> q_A) is tested and rejected, at both rho=320 and rho=384.
* Lemma T.3: beta_bind <= binom(Q,2) / 2^h.  The WP1 candidate tightening
  (q_A^2 / 2^h) is tested and rejected: a final-phase commitment may copy a
  digest already defined earlier in the same experiment.
* Corollary T.4: the ledger, the acceptance certificate
  (non-FCSD sum <= 2^-149.5, checked exactly as sum^2 * 2^299 <= 1 in
  Fractions), the smallest working (rho, h), and the byte sizes.

Writes review/WP2_CHECK.json and target_s/WP2_EXACT_TABLE.csv.
It does not prove Assumption F (WP3) and does not price the reduction time (WP4).
"""
from __future__ import annotations

import csv
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DECISIONS = ROOT / "target_s" / "WP1_DECISIONS.json"
FROZEN = ROOT / "target_s" / "WP2_DECISIONS.json"
OUT_JSON = ROOT / "review" / "WP2_CHECK.json"
OUT_CSV = ROOT / "target_s" / "WP2_EXACT_TABLE.csv"

# WP2 freezes these.  WP1 had rho=320 and treated h=448 as an acceptable floor.
RHO = 384
H = 512
DELTA0_LOG2 = -198.0                 # truncation budget per full reduction execution
ACC = Fraction(-299, 2)              # acceptance line: sum <= 2^-149.5


def log2f(x: Fraction) -> float:
    if x <= 0:
        return float("-inf")
    return math.log2(x.numerator) - math.log2(x.denominator)


def le_pow2(x: Fraction, exp2: Fraction) -> bool:
    """exact test x <= 2**exp2, with half-integer exp2 cleared by squaring"""
    e = 2 * exp2
    if e.denominator != 1:
        raise ValueError("only integer or half-integer exponents are supported")
    e = int(e)
    if e <= 0:
        return x * x * (1 << -e) <= 1
    return x * x <= (1 << e)


def main() -> int:
    d = json.loads(DECISIONS.read_text(encoding="utf-8"))
    p = d["parameters"]
    N, T, L, n = p["N"], p["T"], p["L"], p["n"]
    qA = 1 << p["q_A_log2"]
    rho1, h1 = p["rho"], p["h"]          # WP1 provisional values, for contrast

    checks: list[dict] = []

    def check(name: str, ok: bool, detail: str) -> None:
        checks.append({"name": name, "pass": bool(ok), "detail": detail})

    # ---------------- resources (S.2), reinterpreted as live queries ----
    J = L * T
    R = 1 + J + 3 * N
    Q = R * qA + 3 * J + 6 * N

    check("J = LT reproduces WP1", J == d["resources"]["J"], f"J={J}")
    check("R = 1+LT+3N reproduces WP1", R == d["resources"]["R"], f"R={R}")
    check("log2 Q reproduces WP1",
          math.isclose(math.log2(Q), d["resources"]["Q_log2"], abs_tol=1e-9),
          f"log2 Q={math.log2(Q):.9f} recorded={d['resources']['Q_log2']}")
    # Lemma T.2: R is also the executed-A bound of one experiment, so Q already
    # is the experiment-level live count (1+J+2N full prefixes, N probes).
    # These three are algebraic identities of the (S.2) definitions: they guard
    # against a typo in R/Q; they do not verify Lemma T.2 itself (that is WP4).
    E_exec = 1 + J + 2 * N + N
    check("identity 1+J+2N+N = R (typo guard, not a proof of T.2)", E_exec == R,
          f"1+J+2N+N={E_exec}, R={R}")
    check("identity R >= J+3N+1 (typo guard)",
          R >= J + 3 * N + 1, f"R={R}, J+3N+1={J + 3 * N + 1}")
    check("identity 2J+6N probe commitments fit the 3J+6N charge",
          2 * J + 4 * N + 2 * N <= 3 * J + 6 * N,
          f"probe commitments={2 * J + 6 * N}, charged={3 * J + 6 * N}")

    # WP2 second-pass audit: the hiding union bound prices ONE hidden salt per
    # attempt.  That is a protocol fact of Table 3.1, not an assumption:
    # challenge b opens the pair (1,2), (1,3) or (2,3) of the three salts drawn
    # per attempt, so exactly one commitment per attempt is never opened.  The
    # worst case of three hidden salts per attempt is priced below to show the
    # ledger does not depend on which of the two readings is used.
    openings = {0: (1, 2), 1: (1, 3), 2: (2, 3)}
    check("Table 3.1 leaves exactly one unopened commitment per attempt",
          all(3 - len(v) == 1 for v in openings.values()),
          f"unopened per branch = {[(b, 3 - len(v)) for b, v in openings.items()]}")

    # ---------------- Lemma T.1: hiding ---------------------------------
    two_Q = 2 * Q
    check("guessing budget below the salt space at frozen rho", two_Q < (1 << RHO),
          f"log2 2Q={math.log2(two_Q):.3f} < rho={RHO}")
    epsH = Fraction(two_Q, 1 << RHO)
    hiding = epsH * J
    # WP1 recorded hiding at rho=320 with Q -> q_A; both readings are shown.
    epsH_t = Fraction(two_Q, 1 << rho1)
    hiding_rho320 = epsH_t * J
    hiding_qA_only_rho320 = Fraction(2 * qA, (1 << rho1) - 2 * qA) * J
    hiding_qA_only_frozen = Fraction(2 * qA, (1 << RHO) - 2 * qA) * J
    check("WP1 hiding count (q_A per hybrid) is 36.33 bits below the Lemma T.1 sum",
          hiding_qA_only_rho320 < hiding_rho320,
          f"log2 q_A-only={log2f(hiding_qA_only_rho320):.4f} vs live-query={log2f(hiding_rho320):.4f}"
          f" (gap {log2f(hiding_rho320)-log2f(hiding_qA_only_rho320):.4f}; rejection rests on C.1)")
    check("rho=320 under the Lemma T.1 count misses 2^-128",
          not le_pow2(hiding_rho320, Fraction(-128)),
          f"log2 hiding(rho=320)={log2f(hiding_rho320):.6f}")

    # ---------------- Lemma T.3: binding --------------------------------
    binding = Fraction(Q * (Q - 1), 2 << H)
    binding_qA_only = Fraction(qA * (qA - 1), 2 << H)
    binding_h448 = Fraction(Q * (Q - 1), 2 << 448)
    check("WP1 binding count (q_A^2) is 71.66 bits below the Lemma T.3 bound",
          binding_qA_only < binding,
          f"log2 q_A^2/2^h={log2f(binding_qA_only):.4f} vs binom(Q,2)/2^h={log2f(binding):.4f}"
          f" (gap {log2f(binding)-log2f(binding_qA_only):.4f}; rejection rests on C.2)")
    check("h=448 floor misses the acceptance line",
          binding_h448 > Fraction(1, 1 << 152),
          f"log2 binding(h=448)={log2f(binding_h448):.6f}")

    # ---------------- Corollary T.4: ledger -----------------------------
    knowledge = Fraction(2, 3) ** N
    retry = L * Fraction(2, 3) ** T
    truncation = Fraction(2 * J + 2, 1 << 198)
    terms = {"knowledge": knowledge, "retry": retry, "hiding": hiding,
             "binding": binding, "truncation": truncation}
    non_core = sum(terms.values(), Fraction(0))
    eps_f = Fraction(1, 1 << 136)
    total = non_core + eps_f

    for name in ("hiding", "binding", "retry", "truncation"):
        check(f"{name} <= knowledge term", terms[name] <= knowledge,
              f"log2={log2f(terms[name]):.6f} vs knowledge {log2f(knowledge):.6f}")
    check("WP2 acceptance: non-FCSD sum <= 2^-149.5 (exact certificate)",
          le_pow2(non_core, ACC),
          f"log2 non_core={log2f(non_core):.9f}")
    check("non-FCSD sum <= 2^-128", le_pow2(non_core, Fraction(-128)),
          f"log2 non_core={log2f(non_core):.9f}")
    check("total <= 2^-128 (exact certificate)", le_pow2(total, Fraction(-128)),
          f"log2 total={log2f(total):.9f}")
    check("total keeps >= 7 bits of margin", total * (1 << 135) <= 1,
          f"log2 total={log2f(total):.9f}")
    check("non-FCSD sum < FCSD term", non_core < eps_f,
          f"log2 non_core={log2f(non_core):.6f} < -136")

    # Second-pass audit items.
    # (a) robustness to the hidden-salt count: three salts per attempt instead
    #     of one costs only log2(3)=1.585 bits and stays under the knowledge term.
    hiding_3salt = epsH * (3 * J)
    non_core_3salt = knowledge + retry + hiding_3salt + binding + truncation
    check("ledger unchanged if three salts per attempt were hidden",
          le_pow2(non_core_3salt, ACC) and hiding_3salt <= knowledge,
          f"log2 hiding(3J)={log2f(hiding_3salt):.6f}, sum={log2f(non_core_3salt):.9f}")
    # (b) the 2^-128 margin belongs to the FCSD term, not to the ledger:
    #     total = 2^-136 + (much smaller) so margin = 7.9999 bits (= 8 - 0.000105,
    #     because the total sits 0.000105 bits above 2^-136), while the
    #     non-FCSD sum alone sits 21.7487 bits under 2^-128.
    margin_total = -128 - log2f(total)
    margin_core = -128 - log2f(non_core)
    check("total margin over 2^-128 is 7.9999 bits (assumption-dominated)",
          7.99 < margin_total < 8.01, f"margin_total={margin_total:.6f}")
    check("non-FCSD margin over 2^-128 is 21.75 bits",
          21.7 < margin_core < 21.8, f"margin_non_core={margin_core:.6f}")
    check("FCSD term is the binding constraint on the total",
          eps_f <= total <= 2 * eps_f,
          f"log2 total={log2f(total):.9f} vs eps_F=-136 (2*eps_F = -135.0)")

    # ---------------- sensitivity: smallest working (rho, h) ------------
    def ledger(rho: int, h: int) -> tuple[Fraction, Fraction, Fraction]:
        hid = (Fraction(two_Q, 1 << rho) if two_Q < (1 << rho) else Fraction(1)) * J
        bind = Fraction(Q * (Q - 1), 2 << h)
        return hid, bind, knowledge + retry + hid + bind + truncation

    rho_min = min(rho for rho in range(8, 1025, 8) if le_pow2(ledger(rho, H)[2], ACC))
    h_min = min(h for h in range(8, 1025, 8) if le_pow2(ledger(RHO, h)[2], ACC))
    check("frozen rho clears the acceptance line", RHO >= rho_min,
          f"rho_min={rho_min}, frozen={RHO}, margin={RHO - rho_min} bits")
    check("rho just below the minimum fails",
          not le_pow2(ledger(rho_min - 8, H)[2], ACC),
          f"log2 sum(rho={rho_min - 8})={log2f(ledger(rho_min - 8, H)[2]):.6f}")
    check("frozen h clears the acceptance line", H >= h_min,
          f"h_min={h_min}, frozen={H}, margin={H - h_min} bits")
    check("h just below the minimum fails",
          not le_pow2(ledger(RHO, h_min - 8)[2], ACC),
          f"log2 sum(h={h_min - 8})={log2f(ledger(RHO, h_min - 8)[2]):.6f}")

    # ---------------- cross-check against WP2_DECISIONS.json ------------
    fz = json.loads(FROZEN.read_text(encoding="utf-8"))
    fp, fl = fz["parameters"], fz["ledger"]
    check("frozen rho in WP2_DECISIONS matches the script", fp["rho"] == RHO,
          f"json={fp['rho']}, script={RHO}")
    check("frozen h in WP2_DECISIONS matches the script", fp["h"] == H,
          f"json={fp['h']}, script={H}")
    check("frozen J/R match", (fp.get("J") or fl.get("J")) in (None, J) and
          fz["resources"]["J"] == J and fz["resources"]["R"] == R,
          f"J={J}, R={R}")
    for key, val in [("knowledge_log2", knowledge), ("retry_log2", retry),
                     ("hiding_log2", hiding), ("binding_log2", binding),
                     ("truncation_log2", truncation),
                     ("non_core_log2", non_core), ("total_log2", total)]:
        check(f"WP2_DECISIONS ledger {key}",
              math.isclose(fl[key], log2f(val), abs_tol=1e-9),
              f"json={fl[key]}, computed={log2f(val):.12f}")
    rj = fz["rejected_wp1_tightenings"]
    for key, val in [("hiding_qA_only_at_rho320_log2", hiding_qA_only_rho320),
                     ("hiding_correct_at_rho320_log2", hiding_rho320),
                     ("hiding_correct_at_rho384_log2", hiding),
                     ("binding_qA_only_at_h512_log2", binding_qA_only),
                     ("binding_correct_at_h512_log2", binding),
                     ("binding_correct_at_h448_log2", binding_h448)]:
        check(f"WP2_DECISIONS rejected {key}",
              math.isclose(rj[key], log2f(val), abs_tol=1e-9),
              f"json={rj[key]}, computed={log2f(val):.12f}")
    sen = fz["sensitivity"]
    check("WP2_DECISIONS sensitivity rho_min/h_min",
          sen["rho_min_for_acceptance"] == rho_min and
          sen["h_min_for_acceptance"] == h_min,
          f"json rho={sen['rho_min_for_acceptance']} h={sen['h_min_for_acceptance']}, computed {rho_min}/{h_min}")

    # ---------------- T and delta0 sanity -------------------------------
    T_min = min(t for t in range(1, 4096)
                if L * Fraction(2, 3) ** t <= Fraction(1, 1 << 160))
    check("T=322 is the smallest cap with L(2/3)^T <= 2^-160", T_min == T,
          f"T_min={T_min}, frozen={T}")
    check("WP2_DECISIONS sensitivity T_min",
          json.loads(FROZEN.read_text(encoding="utf-8"))["sensitivity"]["T_min_for_retry_below_2^-160"] == T_min,
          f"computed T_min={T_min}")
    check("delta0=2^-198 keeps the truncation term below knowledge",
          truncation <= knowledge, f"log2 truncation={log2f(truncation):.6f}")

    # ---------------- sizes ---------------------------------------------
    salt_b, dig_b = RHO // 8, H // 8
    salt_b1, dig_b1 = rho1 // 8, h1 // 8
    # Only the unambiguous deltas are priced here; the per-round payload of the
    # response leg depends on how sigma is encoded and is a WP5 measurement.
    comm_a = 3 * dig_b                       # commitment leg: three digests
    salt_delta_per_round = 2 * (salt_b - salt_b1)   # two openings carried per round
    digest_delta_commitment_leg = 3 * (dig_b - 384 // 8)
    # Fixed per-round overhead, counted exactly as Ch.6 of the manuscript counts
    # 161 B: three commitments, one challenge byte and the two salts that Table
    # 3.1 opens.  The vector/permutation payload (n*B_n, n*B_p terms) is WP5.
    ovh_security = 3 * dig_b + 1 + 2 * salt_b
    ovh_measured = 3 * 32 + 1 + 2 * 32      # Ch.6 measurement: SHA-256, 32-byte salts
    check("salts and digests are byte integral", RHO % 8 == 0 and H % 8 == 0,
          f"salt={salt_b}B, digest={dig_b}B")
    check("frozen salts and digests are at least as large as WP1's",
          salt_b > salt_b1 and dig_b >= dig_b1,
          f"salt {salt_b1}->{salt_b}B, digest {dig_b1}->{dig_b}B")
    check("fixed per-round overhead is 289 B against the measured 161 B",
          ovh_security == 289 and ovh_measured == 161 and ovh_security - ovh_measured == 128,
          f"{ovh_security} B vs {ovh_measured} B, +{ovh_security - ovh_measured} B/round,"
          f" +{(ovh_security - ovh_measured) * N // 1024} KiB per {N}-round session")

    # ---------------- diagnostic: size of S_C ---------------------------
    # Second-pass audit correction.  Spreading n=1024 coordinates as evenly as
    # possible over all 251 field values MAXIMISES the multinomial, so 7571.7
    # bits is an UPPER BOUND on log2|S_C| over every composition, not the value
    # of a usable composition.  A positive multiplicity at the bound says nothing
    # about a frozen C.  The protocol compositions actually recorded for n=640
    # (2147.29-2408.04 bits over 640 coordinates) sit at 3.36-3.76 bits per
    # coordinate; the same density at n=1024 with k=512 gives 3435.7-3852.9
    # bits, still below log2 p^k = 4081.4 bits.  So the sign of the multiplicity
    # at n=1024 is undetermined until WP3 freezes C.
    def log2_fact(m: int) -> float:
        return math.lgamma(m + 1) / math.log(2)

    counts = [n // p["p"]] * p["p"]
    for i in range(n - sum(counts)):
        counts[i] += 1
    log_Sc_max = log2_fact(n) - sum(log2_fact(c) for c in counts)
    log_pk = p["k"] * math.log2(p["p"])
    n640 = (2147.29, 2408.04)
    log_Sc_scaled = tuple(v * n / 640 for v in n640)
    check("even spread over 251 values maximises the multinomial (upper bound)",
          log_Sc_max > max(log_Sc_scaled),
          f"bound={log_Sc_max:.1f} vs equal-density n=1024 estimate "
          f"{log_Sc_scaled[0]:.1f}-{log_Sc_scaled[1]:.1f}")
    check("multiplicity sign at n=1024 is undetermined before WP3",
          log_Sc_scaled[0] < log_pk < log_Sc_max,
          f"log2 p^k={log_pk:.1f} lies inside [{log_Sc_scaled[0]:.1f}, {log_Sc_max:.1f}]")
    check("n=640 recorded compositions do sit below p^k",
          max(n640) < 320 * math.log2(p["p"]),
          f"n640 max={max(n640):.1f} < log2 p^320={320 * math.log2(p['p']):.1f}")
    log_Sc, multiplicity = log_Sc_max, log_Sc_max - log_pk

    rows = [
        ("parameters", "p", p["p"], ""),
        ("parameters", "n", n, ""),
        ("parameters", "k", p["k"], ""),
        ("parameters", "N", N, "final verification rounds"),
        ("parameters", "T", T, "retries per learning round"),
        ("parameters", "L", L, "learning rounds = 2^20 sessions x N"),
        ("parameters", "J", J, "hybrids = LT"),
        ("parameters", "R", R, "1+J+3N; also the executed-A bound of Lemma T.2"),
        ("parameters", "q_A", qA, "per-branch query cap"),
        ("parameters", "Q", Q, "live RO queries in one experiment"),
        ("parameters", "h", H, "WP2 frozen; 448 rejected"),
        ("parameters", "rho", RHO, "WP2 frozen; 320 rejected"),
        ("hiding", "eps_H per hybrid = 2Q/2^rho", None, f"log2={log2f(epsH):.9f}"),
        ("hiding", "sum J*eps_H (frozen rho=384)", None, f"log2={log2f(hiding):.9f}"),
        ("hiding", "sum with correct count at rho=320", None, f"log2={log2f(hiding_rho320):.9f}"),
        ("hiding", "sum with WP1 q_A-only at rho=320", None, f"log2={log2f(hiding_qA_only_rho320):.9f}"),
        ("hiding", "sum with WP1 q_A-only at rho=384", None, f"log2={log2f(hiding_qA_only_frozen):.9f}"),
        ("binding", "binom(Q,2)/2^512 (frozen)", None, f"log2={log2f(binding):.9f}"),
        ("binding", "q_A^2/2^512 (rejected)", None, f"log2={log2f(binding_qA_only):.9f}"),
        ("binding", "binom(Q,2)/2^448 (rejected)", None, f"log2={log2f(binding_h448):.9f}"),
        ("ledger", "knowledge (2/3)^256", None, f"log2={log2f(knowledge):.9f}"),
        ("ledger", "retry 2^28(2/3)^322", None, f"log2={log2f(retry):.9f}"),
        ("ledger", "truncation (2J+2)2^-198", None, f"log2={log2f(truncation):.9f}"),
        ("ledger", "non-FCSD sum", None, f"log2={log2f(non_core):.9f}"),
        ("ledger", "FCSD assumption F", None, "-136 (assumption; WP3 evidence)"),
        ("ledger", "total", None, f"log2={log2f(total):.9f}"),
        ("sensitivity", "rho_min for 2^-149.5 (h=512)", rho_min, "multiple of 8"),
        ("sensitivity", "h_min for 2^-149.5 (rho=384)", h_min, "multiple of 8"),
        ("sensitivity", "T_min for L(2/3)^T <= 2^-160", T_min, ""),
        ("communication", "salt bytes", salt_b, f"rho={RHO} (WP1 provisional had {salt_b1})"),
        ("communication", "digest bytes", dig_b, f"h={H} (WP1 provisional had {dig_b1})"),
        ("communication", "commitment leg bytes", comm_a, "three digests only"),
        ("communication", "fixed per-round overhead, security", ovh_security,
         "3 digests + 1 challenge byte + 2 salts, same accounting as the measured 161 B"),
        ("communication", "fixed per-round overhead, measured", ovh_measured,
         "Ch.6 configuration: SHA-256 + 32-byte salts"),
        ("communication", "overhead growth per round", ovh_security - ovh_measured,
         f"+{(ovh_security - ovh_measured) * N // 1024} KiB over the {N}-round session"),
        ("diagnostic", "log2 |S_C| upper bound (even spread over 251 values)",
         round(log_Sc_max, 1), "maximum over all compositions, not a usable C"),
        ("diagnostic", "log2 p^k", round(log_pk, 1), ""),
        ("diagnostic", "equal-density n=1024 estimate from the recorded n=640 range",
         f"{log_Sc_scaled[0]:.1f}-{log_Sc_scaled[1]:.1f}",
         "both below log2 p^k: multiplicity sign undetermined until WP3 freezes C"),
        ("audit", "hiding sum if 3 salts per attempt stay hidden", None,
         f"log2={log2f(hiding_3salt):.9f}; sum={log2f(non_core_3salt):.9f}; still under 2^-149.5"),
        ("audit", "margin of total over 2^-128", round(margin_total, 4),
         f"assumption-dominated; non-FCSD margin={margin_core:.4f}"),
    ]

    result = {
        "status": "PASS" if all(c["pass"] for c in checks) else "FAIL",
        "work_package": "WP2",
        "date": "2026-10-29",
        "frozen": {"rho": RHO, "h": H, "salt_bytes": salt_b, "digest_bytes": dig_b,
                   "supersedes_wp1": {"rho": rho1, "h_floor": 448}},
        "resources": {"J": J, "R": R, "Q": Q, "log2_Q": math.log2(Q)},
        "ledger_log2": {
            "eps_H": log2f(epsH), "hiding": log2f(hiding), "binding": log2f(binding),
            "knowledge": log2f(knowledge), "retry": log2f(retry),
            "truncation": log2f(truncation), "non_core": log2f(non_core),
            "fcsd": -136.0, "total": log2f(total),
        },
        "rejected_tightenings": {
            "hiding_correct_rho320_log2": log2f(hiding_rho320),
            "hiding_qA_only_rho320_log2": log2f(hiding_qA_only_rho320),
            "hiding_qA_only_frozen_log2": log2f(hiding_qA_only_frozen),
            "binding_qA_only_log2": log2f(binding_qA_only),
            "binding_h448_log2": log2f(binding_h448),
        },
        "sensitivity": {"rho_min": rho_min, "h_min": h_min, "T_min": T_min},
        "audit": {"margin_total_over_2^-128_bits": margin_total,
                  "margin_non_core_over_2^-128_bits": margin_core,
                  "hiding_if_3_salts_per_attempt_log2": log2f(hiding_3salt),
                  "non_core_if_3_salts_per_attempt_log2": log2f(non_core_3salt),
                  "hidden_salts_per_attempt": 1,
                  "hidden_salts_source": "Table 3.1: each challenge opens 2 of 3 salts"},
        "sizes": {"salt_bytes": salt_b, "digest_bytes": dig_b,
                  "commitment_leg_bytes": comm_a,
                  "overhead_security_bytes": ovh_security,
                  "overhead_measured_bytes": ovh_measured,
                  "overhead_growth_bytes": ovh_security - ovh_measured,
                  "overhead_growth_per_session_KiB": (ovh_security - ovh_measured) * N // 1024,
                  "salt_growth_per_round_bytes": salt_delta_per_round,
                  "digest_growth_commitment_leg_bytes": digest_delta_commitment_leg,
                  "note": "vector/permutation payload needs the sigma encoding; deferred to WP5"},
        "diagnostics": {"log2_S_C_upper_bound": round(log_Sc_max, 3),
                        "log2_p_k": round(log_pk, 3),
                        "equal_density_n1024_from_n640": [round(v, 1) for v in log_Sc_scaled],
                        "multiplicity_sign_at_n1024": "undetermined until WP3 freezes C",
                        "note": "the even-spread multinomial is a maximum over compositions, "
                                "not the value of a usable low-weight composition"},
        "checks": checks,
        "notes": [
            "Hiding is priced by the live queries of one experiment (Q), because the ROM is one fixed function and rewinds reuse it.",
            "Binding keeps binom(Q,2)/2^h: a final-phase commitment can copy a digest already defined earlier in the same experiment.",
            "A private-placeholder variant could restore rho=320 but needs a min-entropy lemma for the hidden message; recorded, not used.",
            "Assumption F(2^166, 2^-136) is unchanged and still unproved here; L_* does not cover a general solver B.",
            "Second-pass audit: the J prefactor is legitimate because each hybrid experiment re-draws its own salt and re-runs A R times, "
            "so hybrid j is not billed from a shared query budget; the q_A cap bounds one execution only.",
            "Second-pass audit: 1+J+2N+N=R and the two R inequalities are definitional identities, kept only as typo guards.",
            "Second-pass audit: rho=384 is a conservative freeze; the binding constraint at that value is the knowledge term, "
            "and the 7.9999-bit margin of the total comes from Assumption F, not from the ledger.",
        ],
    }

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    with OUT_CSV.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["group", "quantity", "value", "detail"])
        for g, q, v, det in rows:
            w.writerow([g, q, "" if v is None else v, det])

    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    print(json.dumps({k: result[k] for k in
                      ("status", "frozen", "resources", "ledger_log2",
                       "rejected_tightenings", "sensitivity", "audit", "sizes", "diagnostics")},
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
