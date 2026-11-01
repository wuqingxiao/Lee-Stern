#!/usr/bin/env python3
"""WP1 arithmetic and consistency check.

Recomputes the resource interface and the Target-S ledger from the frozen
parameters in WP1_DECISIONS.json, compares them with the recorded values, and
writes review/WP1_CHECK.json. This script does not prove the tightened lemmas
(WP2) or the FCSD assumption (WP3); it only checks that the frozen numbers are
internally consistent.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DECISIONS = ROOT / "target_s" / "WP1_DECISIONS.json"
OUT = ROOT / "review" / "WP1_CHECK.json"

TOL = 1e-6


def log2(x: float) -> float:
    return math.log2(x)


def logsumexp2(values: list[float]) -> float:
    m = max(values)
    return m + log2(sum(2.0 ** (v - m) for v in values))


def close(a: float, b: float, tol: float = TOL) -> bool:
    return abs(a - b) <= tol


def main() -> int:
    d = json.loads(DECISIONS.read_text(encoding="utf-8"))
    p = d["parameters"]
    N, T, h, rho = p["N"], p["T"], p["h"], p["rho"]
    L = p["L"]
    qA = 2 ** p["q_A_log2"]

    checks: list[dict] = []

    def check(name: str, ok: bool, detail: str) -> None:
        checks.append({"name": name, "pass": bool(ok), "detail": detail})

    # Parameter consistency
    check("L = 2^20 sessions * N rounds", L == (2 ** 20) * N,
          f"L={L}, 2^20*N={2 ** 20 * N}")
    check("N=256 (security)", N == 256, f"N={N}")
    check("T=322 (conservative retry cap)", T == 322, f"T={T}")
    check("h=512", h == 512, f"h={h}")
    check("rho=320", rho == 320, f"rho={rho}")

    # Resource interface
    R = 1 + L * T + 3 * N
    Q = R * qA + 3 * L * T + 6 * N
    J = L * T
    rec = d["resources"]
    check("R = 1+LT+3N", R == rec["R"], f"R={R}, recorded={rec['R']}")
    check("log2 R", close(log2(R), rec["R_log2"]),
          f"log2 R={log2(R):.12f}, recorded={rec['R_log2']}")
    check("log2 Q", close(log2(Q), rec["Q_log2"]),
          f"log2 Q={log2(Q):.12f}, recorded={rec['Q_log2']}")
    check("J = LT", J == rec["J"], f"J={J}, recorded={rec['J']}")
    check("log2 J", close(log2(J), rec["J_log2"]),
          f"log2 J={log2(J):.12f}, recorded={rec['J_log2']}")
    check("Assumption F work covers Q", log2(Q) < d["assumption_F"]["W_log2"],
          f"log2 Q={log2(Q):.6f} < W={d['assumption_F']['W_log2']}")

    # Ledger terms
    knowledge = N * log2(2.0 / 3.0)
    retry = p["L_log2"] + T * log2(2.0 / 3.0)
    hiding = log2(J) + 1 + p["q_A_log2"] - rho
    binding = 2 * p["q_A_log2"] - h
    delta0 = -198.0
    truncation = log2(2 * J + 2) + delta0
    non_core = logsumexp2([knowledge, retry, hiding, binding, truncation])
    eps_f = d["assumption_F"]["epsilon_log2"]
    total = logsumexp2([eps_f, non_core])

    led = d["ledger"]
    for name, value, key in [
        ("knowledge", knowledge, "knowledge_log2"),
        ("retry", retry, "retry_log2"),
        ("hiding", hiding, "hiding_log2"),
        ("binding", binding, "binding_log2"),
        ("non_core", non_core, "non_core_log2"),
        ("total", total, "total_log2"),
    ]:
        check(f"ledger {name}", close(value, led[key]),
              f"computed={value:.12f}, recorded={led[key]}")

    check("binding at h=448 is -192", 2 * p["q_A_log2"] - 448 == -192,
          "2*128-448=-192")
    check("non-core below 2^-128", non_core < -128,
          f"non-core log2={non_core:.6f}")
    check("total below 2^-128", total < -128,
          f"total log2={total:.9f}")
    check("total at least 7 bits of margin", total <= -135,
          f"total log2={total:.9f}")

    # Performance separation
    k219 = 219 * log2(2.0 / 3.0)
    ratio = 256 / 219
    check("219-round knowledge is only a reliability figure",
          close(k219, -128.10678765793324, 1e-9),
          f"log2 (2/3)^219={k219:.9f}")
    check("256/219 conversion", close(ratio, 1.1689497716894977, 1e-12),
          f"ratio={ratio:.12f}")

    result = {
        "status": "PASS" if all(c["pass"] for c in checks) else "FAIL",
        "work_package": "WP1",
        "date": d["date"],
        "parameters": {"p": p["p"], "n": p["n"], "k": p["k"],
                       "N": N, "T": T, "h": h, "rho": rho,
                       "L": L, "sessions": 2 ** 20},
        "resources": {"R": R, "log2_R": log2(R), "log2_Q": log2(Q),
                      "J": J, "log2_J": log2(J)},
        "ledger_log2": {
            "knowledge": knowledge, "retry": retry, "hiding": hiding,
            "binding": binding, "truncation": truncation,
            "non_core": non_core, "fcsd": eps_f, "total": total,
        },
        "checks": checks,
        "notes": [
            "The tightened hiding and binding lemmas are WP2 targets, not proved here.",
            "The FCSD bound is an assumption (WP3 evidence), not a theorem.",
            "219 rounds and SHA-256 are performance settings; they are not used in the security ledger.",
        ],
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")

    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    print(json.dumps({k: result[k] for k in
                      ("status", "resources", "ledger_log2")},
                     ensure_ascii=False, indent=2))
    failed = [c for c in checks if not c["pass"]]
    if failed:
        print("FAILED CHECKS:")
        for c in failed:
            print(" -", c["name"], c["detail"])
        return 1
    print(f"All {len(checks)} checks passed; wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
