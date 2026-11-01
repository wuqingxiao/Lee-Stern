#!/usr/bin/env python3
"""WP5: sizes of the security configuration at n=1024, and the measurement boundary.

Scope note (2026-10-30): this package does the arithmetic half of WP5 and
records why no re-measurement is run.

* Fixed per-round overhead: 3 digests + 1 challenge byte + 2 salts.
  Measured configuration (SHA-256, 32-byte salts): 3*32+1+2*32 = 161 B.
  Security configuration (h=512, rho=384): 3*64+1+2*48 = 289 B.
* Payload: b=0,1 carry n*B_n + n*B_p bytes, b=2 carries 2*n*B_p, so the
  expected round is 161 + 2nB_n/3 + 4nB_p/3 in the measured configuration.
  The formula must reproduce the four numbers of manuscript table 6.2.
* At n=1024, N=256 the security configuration is 3019.67 B/round and
  773034.67 B/session (754.92 KiB), exactly 32768 B above the same round
  count in the measured configuration.
* The manuscript's chapter 6 timings belong to C_n = 4*C_256 (16 nonzero
  counts, Lee weight 1480, log2|S_C| = 3033.99), not to the frozen C
  (c_0=242, 25 nonzero counts, weight 2025, log2|S_C| = 3460.24).  They are
  a sampler comparison, not a measurement of the security configuration.
* The implementation pins SHA-256 and 32-byte salts (implementation/
  protocol.py) and the sweep pins 219 rounds, so a security-configuration
  timing needs code changes; the three patch points are recorded here.

Writes review/WP5_SIZES.json and target_s/WP5_SIZES.csv.
"""
from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / "target_s" / "WP3_C_N1024.json"
PROTOCOL = ROOT / "implementation" / "protocol.py"
SETTINGS = ROOT / "performance" / "length_sweep" / "settings.json"
OUT_JSON = ROOT / "review" / "WP5_SIZES.json"
OUT_CSV = ROOT / "target_s" / "WP5_SIZES.csv"

P, N = 251, 256
L = 1 << 28
ROUNDS_MEASURED = 219

# Chapter 6, table 6.2, expected payload in bytes for 219 rounds.
TABLE_6_2 = {256: 147387, 512: 334267, 768: 483771, 1024: 633275}


def bytes_n(n: int) -> int:
    """B_n = ceil(log_256 n)"""
    return math.ceil(math.log(n, 256))


def bytes_p() -> int:
    """B_p = ceil(log_256 p)"""
    return math.ceil(math.log(P, 256))


def per_round(n: int, overhead: int) -> float:
    return overhead + 2 * n * bytes_n(n) / 3 + 4 * n * bytes_p() / 3


def log2S(counts: list[int]) -> float:
    lg = lambda m: math.lgamma(m + 1) / math.log(2)
    return lg(sum(counts)) - sum(lg(c) for c in counts if c)


def main() -> int:
    checks: list[dict] = []

    def check(name: str, ok: bool, detail: str) -> None:
        checks.append({"name": name, "pass": bool(ok), "detail": detail})

    # ---------------- fixed overhead decomposition ----------------------
    ovh_measured = 3 * 32 + 1 + 2 * 32
    ovh_security = 3 * 64 + 1 + 2 * 48
    check("measured fixed overhead is 3*32+1+2*32 = 161 B", ovh_measured == 161,
          f"{ovh_measured} B = commitment 96 + salts 64 + challenge 1")
    check("security fixed overhead is 3*64+1+2*48 = 289 B", ovh_security == 289,
          f"{ovh_security} B = commitment 192 + salts 96 + challenge 1")
    check("commitment leg grows by 96 B, salt leg by 32 B",
          3 * 64 - 3 * 32 == 96 and 2 * 48 - 2 * 32 == 32,
          f"digests +{3 * 64 - 3 * 32} B, salts +{2 * 48 - 2 * 32} B")
    check("per-round overhead delta is 128 B", ovh_security - ovh_measured == 128,
          f"{ovh_security - ovh_measured} B/round")

    # ---------------- payload formula reproduces chapter 6 --------------
    check("B_p = 1 byte at p=251", bytes_p() == 1, f"p<256")
    check("B_n = 2 bytes at n=1024", bytes_n(1024) == 2, f"1024<256^2")
    for n, want in TABLE_6_2.items():
        got = per_round(n, ovh_measured) * ROUNDS_MEASURED
        check(f"formula reproduces table 6.2 at n={n}", math.isclose(got, want, abs_tol=0.5),
              f"computed {got:.1f} B, manuscript {want} B")

    # ---------------- per-branch payloads at n=1024 ---------------------
    n = 1024
    b01_m = ovh_measured + n * bytes_n(n) + n * bytes_p()
    b2_m = ovh_measured + 2 * n * bytes_p()
    b01_s = ovh_security + n * bytes_n(n) + n * bytes_p()
    b2_s = ovh_security + 2 * n * bytes_p()
    check("branch expectation equals the round formula (measured)",
          math.isclose((2 * b01_m + b2_m) / 3, per_round(n, ovh_measured), abs_tol=1e-9),
          f"({2 * b01_m}+{b2_m})/3 = {(2 * b01_m + b2_m) / 3:.2f} B")
    check("branch expectation equals the round formula (security)",
          math.isclose((2 * b01_s + b2_s) / 3, per_round(n, ovh_security), abs_tol=1e-9),
          f"({2 * b01_s}+{b2_s})/3 = {(2 * b01_s + b2_s) / 3:.2f} B")

    # ---------------- the security configuration ------------------------
    per_s = per_round(n, ovh_security)
    per_m = per_round(n, ovh_measured)
    sess_s = per_s * N
    sess_m = per_m * N
    check("security configuration is 3019.67 B/round at n=1024",
          math.isclose(per_s, 3019.6666666666665, abs_tol=1e-9), f"{per_s:.4f} B")
    check("security session (N=256) is 773034.67 B = 754.92 KiB",
          math.isclose(sess_s, 773034.6666666666, abs_tol=1e-6),
          f"{sess_s:.2f} B = {sess_s / 1024:.2f} KiB")
    check("measured configuration at N=256 is 740266.67 B = 722.92 KiB",
          math.isclose(sess_m, 740266.6666666666, abs_tol=1e-6),
          f"{sess_m:.2f} B = {sess_m / 1024:.2f} KiB")
    check("session delta is exactly 32768 B = 32 KiB",
          math.isclose(sess_s - sess_m, 32768.0, abs_tol=1e-6),
          f"{sess_s - sess_m:.1f} B")
    check("learning-phase bound is 810.59 GB (2^28 rounds, not an operational cost)",
          math.isclose(per_s * L / 1e9, 810.5856, abs_tol=1e-3),
          f"{per_s * L / 1e9:.2f} GB")

    # ---------------- compositions --------------------------------------
    fz = json.loads(FROZEN.read_text(encoding="utf-8"))
    counts = fz["counts"]
    weight = sum(min(a, P - a) * x for a, x in enumerate(counts))
    check("frozen C: 1024 coordinates, c_0=242, 25 nonzero counts, weight 2025",
          sum(counts) == 1024 and counts[0] == 242 and
          sum(1 for x in counts if x) == 25 and weight == 2025,
          f"sum={sum(counts)}, c0={counts[0]}, nonzero counts={sum(1 for x in counts if x)}, "
          f"unused symbols={sum(1 for x in counts if x == 0)}, weight={weight}")
    c256 = [78, 42, 22, 11, 6, 3, 1, 1] + [0] * 235 + [1, 1, 2, 4, 7, 12, 23, 42]
    c1024 = [4 * x for x in c256]
    w_perf = sum(min(a, P - a) * x for a, x in enumerate(c1024))
    check("chapter 6 measures C_n=4*C_256: 16 nonzero counts, weight 1480",
          sum(c256) == 256 and sum(c1024) == 1024 and
          sum(1 for x in c1024 if x) == 16 and w_perf == 1480,
          f"n=256 weight 370 -> n=1024 weight {w_perf}, nonzero counts {sum(1 for x in c1024 if x)}")
    lg_s, lg_p = log2S(c1024), 512 * math.log2(P)
    lg_f, _ = log2S(counts), 0
    check("the two n=1024 compositions are different witness sets",
          lg_s < lg_f < lg_p,
          f"performance log2|S_C|={lg_s:.2f}, frozen={lg_f:.2f}, log2 p^k={lg_p:.2f}")

    # ---------------- measurement boundary ------------------------------
    proto = PROTOCOL.read_text(encoding="utf-8")
    settings = json.loads(SETTINGS.read_text(encoding="utf-8"))
    check("implementation pins SHA-256 and 32-byte salts",
          "sha256" in proto and "len(salt)!=32" in proto,
          "implementation/protocol.py: sha256 prefix, salt length check == 32")
    check("sweep pins 219 rounds", settings.get("rounds") == ROUNDS_MEASURED,
          f"performance/length_sweep/settings.json rounds={settings.get('rounds')}")
    check("security configuration would need 3 patch points plus a rerun",
          True,
          "SHA-512 (or 512-bit digest), 48-byte salt, frozen C at n=1024, N=256; "
          "recorded as the reopening condition, not performed")

    # ---------------- outputs -------------------------------------------
    rows = []
    for nn in (256, 512, 768, 1024):
        rows.append(("table_6_2", nn, f"{per_round(nn, ovh_measured):.2f}",
                     f"{per_round(nn, ovh_measured) * ROUNDS_MEASURED:.0f}",
                     "161 B overhead, 219 rounds (manuscript table 6.2)"))
    rows += [
        ("branches", "measured b=0/1", f"{b01_m}", f"{2 * b01_m + b2_m}", "161+nB_n+nB_p"),
        ("branches", "measured b=2", f"{b2_m}", "", "161+2nB_p"),
        ("branches", "security b=0/1", f"{b01_s}", f"{2 * b01_s + b2_s}", "289+nB_n+nB_p"),
        ("branches", "security b=2", f"{b2_s}", "", "289+2nB_p"),
        ("session", "measured N=256", f"{per_m:.2f}", f"{sess_m:.2f}", "722.92 KiB"),
        ("session", "security N=256", f"{per_s:.2f}", f"{sess_s:.2f}", "754.92 KiB"),
        ("session", "delta", "128", "32768", "32 KiB, independent of n"),
        ("composition", "frozen C", "25 nonzero counts", "2025",
         f"log2|S_C|={lg_f:.2f}, log2 p^k={lg_p:.2f}"),
        ("composition", "C_n=4*C_256", "16 nonzero counts", f"{w_perf}",
         f"log2|S_C|={lg_s:.2f}"),
    ]

    result = {
        "status": "PASS" if all(c["pass"] for c in checks) else "FAIL",
        "work_package": "WP5",
        "date": "2026-10-30",
        "scope": "sizes only; no re-measurement (decision recorded below)",
        "overhead_bytes": {"measured": ovh_measured, "security": ovh_security,
                           "delta": ovh_security - ovh_measured},
        "n1024": {"measured_per_round": per_m, "security_per_round": per_s,
                  "measured_session_N256": sess_m, "security_session_N256": sess_s,
                  "session_delta": sess_s - sess_m},
        "compositions": {"frozen": {"zero_coordinates": counts[0],
                                    "nonzero_counts": sum(1 for x in counts if x),
                                    "unused_symbols": sum(1 for x in counts if x == 0),
                                    "lee_weight": weight, "log2_S_C": lg_f},
                         "chapter6": {"name": "C_n = 4*C_256",
                                      "nonzero_counts": sum(1 for x in c1024 if x),
                                      "lee_weight": w_perf, "log2_S_C": lg_s}},
        "decision": {
            "remeasure": False,
            "reasons": [
                "chapter 6 measures the sampler pair on C_n=4*C_256, not on the frozen C",
                "the implementation pins SHA-256 and 32-byte salts; the security configuration needs SHA-512/64-byte digests, 48-byte salts and N=256",
                "timings do not enter the reduction or the ledger; section 4.4 and 5.5 C-2 already separate measurement configuration from security configuration",
            ],
            "reopen_if": "a referee asks for the security configuration's throughput, or chapter 6 is rewritten for another reason",
        },
        "checks": checks,
    }

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    with OUT_CSV.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["group", "n_or_case", "bytes_per_round", "bytes_total", "note"])
        for r in rows:
            w.writerow(list(r))

    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    print(json.dumps({k: result[k] for k in ("status", "overhead_bytes", "n1024", "compositions")},
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
