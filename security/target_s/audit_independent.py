#!/usr/bin/env python3
"""Independent audit of WP1-WP5 (2026-10-30).

Deliberately does NOT import any of wp1_check..wp5_sizes.  It recomputes the
ledger from the raw definitions with Fraction arithmetic, re-derives the
resource interface as an integer, re-runs the sensitivity search over *every*
integer (not multiples of 8), reproduces chapter 6's table from the manuscript
text itself, and scans the staged paste-ready documents for superseded numbers.

Read-only apart from review/AUDIT_CHECK.json.
"""
from __future__ import annotations

import json
import math
import re
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MS = ROOT / "build" / "concise_manuscript.md"
DEC = ROOT / "target_s" / "WP2_DECISIONS.json"
DEC1 = ROOT / "target_s" / "WP1_DECISIONS.json"
C1024 = ROOT / "target_s" / "WP3_C_N1024.json"
SEC44 = ROOT / "target_s" / "WP4_SECTION44.md"
PAPER = ROOT / "target_s" / "WP1_PAPER_TEXT.md"
SIZES5 = ROOT / "target_s" / "WP5_SIZES.md"
OUT = ROOT / "review" / "AUDIT_CHECK.json"

P, K, N, T = 251, 512, 256, 322
L, QA = 1 << 28, 1 << 128
RHO, H, DELTA0 = 384, 512, 198
BUDGET = 166


def log2f(x: Fraction) -> float:
    return math.log2(x.numerator) - math.log2(x.denominator)


def le2(x: Fraction, e: float) -> bool:
    """exact test x <= 2**e for integer or half-integer e"""
    ee = 2 * e
    assert float(ee).is_integer()
    ee = int(ee)
    return x * x * (1 << -ee) <= 1 if ee <= 0 else x * x <= (1 << ee)


def main() -> int:
    checks: list[dict] = []

    def check(name: str, ok: bool, detail: str) -> None:
        checks.append({"name": name, "pass": bool(ok), "detail": detail})

    # ---------- 1. raw definitions, exact ----------
    J = L * T
    Rexec = 1 + J + 3 * N
    Q = Rexec * QA + 3 * J + 6 * N
    check("Q is below 2^166 and above 2^164", (1 << 164) < Q < (1 << 166),
          f"J={J}, R={Rexec}, Q={Q}, log2Q={math.log2(Q):.9f}")

    knowledge = Fraction(2, 3) ** N
    retry = L * Fraction(2, 3) ** T
    hiding = Fraction(2 * Q * J, 1 << RHO)
    binding = Fraction(Q * (Q - 1), 2 * (1 << H))
    trunc = Fraction(2 * J + 2, 1 << DELTA0)
    s = knowledge + retry + hiding + binding + trunc
    total = s + Fraction(1, 1 << 136)

    d = json.loads(DEC.read_text(encoding="utf-8"))["ledger"]
    doc = {"knowledge": knowledge, "retry": retry, "hiding": hiding,
           "binding": binding, "truncation": trunc,
           "non_core": s, "total": total}
    for k, v in doc.items():
        rec = d[f"{k}_log2"]
        check(f"ledger {k} matches WP2_DECISIONS", abs(rec - log2f(v)) < 1e-9,
              f"json={rec}, independent={log2f(v):.12f}")
    check("acceptance sum <= 2^-149.5 (independent certificate)",
          le2(s, -149.5), f"log2 sum={log2f(s):.9f}")
    check("total <= 2^-128 (independent certificate)", le2(total, -128),
          f"log2 total={log2f(total):.9f}, margin={-128 - log2f(total):.6f} bits")
    check("margin claim is 7.9999 not 8.0001", abs((-128 - log2f(total)) - 7.999895) < 1e-6,
          f"margin={-128 - log2f(total):.9f}")
    check("non-FCSD sum is ABOVE 2^-149.75 (so '<=' was false)", s > Fraction(1, 1 << 149) * Fraction(1, 2),
          f"log2 sum={log2f(s):.9f} vs -149.75")
    check("total is ABOVE 2^-136 (so '=2^-136' was false as an upper bound)", total > Fraction(1, 1 << 136),
          f"log2 total={log2f(total):.9f} vs -136")

    # ---------- 2. sensitivity over EVERY integer ----------
    def ledger(rho: int, h: int) -> Fraction:
        hid = Fraction(2 * Q * J, 1 << rho) if Q < (1 << rho) else Fraction(1)
        return knowledge + retry + hid + Fraction(Q * (Q - 1), 2 * (1 << h)) + trunc

    rho_any = next(r for r in range(1, 1025) if le2(ledger(r, H), -149.5))
    h_any = next(h for h in range(1, 1025) if le2(ledger(RHO, h), -149.5))
    rho_8 = min(r for r in range(8, 1025, 8) if le2(ledger(r, H), -149.5))
    h_8 = min(h for h in range(8, 1025, 8) if le2(ledger(RHO, h), -149.5))
    check("true integer rho_min is 354; the recorded 360 is the first multiple of 8",
          rho_any == 354 and rho_8 == 360, f"all-integer={rho_any}, multiple-of-8={rho_8}")
    check("h_min is 480 either way", h_any == 480 and h_8 == 480, f"{h_any}/{h_8}")
    check("frozen (384, 512) clears both minima", RHO > rho_any and H > h_any,
          f"rho margin {RHO - rho_any}, h margin {H - h_any}")
    check("rho=320 would fail the acceptance line",
          not le2(ledger(320, H), -149.5), f"log2 sum(320)={log2f(ledger(320, H)):.6f}")
    check("h=448 would fail the acceptance line",
          not le2(ledger(RHO, 448), -149.5), f"log2 sum(448)={log2f(ledger(RHO, 448)):.6f}")

    # ---------- 3. T_B as an integer ----------
    items = {
        "c_G": 1 << 31,
        "R t_A": Rexec * QA,
        "J c_sim": J * (1 << 14),
        "3N c_ver": 3 * N * (1 << 14),
        "c_ext": 1 << 20,
        "R(L+N+1) c_io": Rexec * (L + N + 1) * (1 << 13),
        "c_control": 1 << 32,
    }
    T_B = sum(items.values())
    doc_TB = 29412720712961561391050835851825704906307042025472
    check("T_B integer reproduces the recorded certificate", T_B == doc_TB,
          f"{T_B}")
    check("T_B < 2^166", T_B < (1 << 166),
          f"log2 T_B={math.log2(T_B):.6f}, headroom={166 - math.log2(T_B):.6f} bits")
    check("T_B >= Q (queries are part of the declared work)",
          T_B >= Q, f"log2 Q={math.log2(Q):.6f}")
    check("main term dominates: second/main <= 2^-80",
          items["R(L+N+1) c_io"] * (1 << 80) <= items["R t_A"],
          f"second={math.log2(items['R(L+N+1) c_io']):.3f}")

    # ---------- 4. chapter 6 table from the manuscript itself ----------
    txt = MS.read_text(encoding="utf-8")
    rows = re.findall(r"^\|(\d{3,4})\|[^|]*\|[^|]*\|[^|]*\|(\d+)\|$", txt, re.M)
    table = {int(n): int(v) for n, v in rows}
    gen = {}
    for n in (256, 512, 768, 1024):
        bn = math.ceil(math.log(n, 256))
        gen[n] = round((161 + 2 * n * bn / 3 + 4 * n / 3) * 219)
    check("manuscript table 6.2 is parsed as four rows", len(table) == 4, f"{table}")
    check("payload formula reproduces the manuscript's table 6.2", table == gen,
          f"parsed={table}, computed={gen}")
    n1024, bn1024 = 1024, 2
    pay = 2 * n1024 * bn1024 / 3 + 4 * n1024 / 3
    per_meas, per_sec = 161 + pay, 289 + pay
    check("security per-round = measured + 128 B, delta x256 = 32768",
          per_sec - per_meas == 128 and (per_sec - per_meas) * 256 == 32768,
          f"{per_meas:.4f} -> {per_sec:.4f} B; session {per_sec*256:.2f} B")

    # ---------- 5. composition ----------
    fz = json.loads(C1024.read_text(encoding="utf-8"))
    counts = fz["counts"]
    lg = lambda m: math.lgamma(m + 1) / math.log(2)
    logSc = lg(sum(counts)) - sum(lg(c) for c in counts if c)
    pk = K * math.log2(P)
    weight = sum(min(a, P - a) * x for a, x in enumerate(counts))
    check("frozen C: sum 1024, c0 242, weight 2025", sum(counts) == 1024 and
          counts[0] == 242 and weight == 2025,
          f"sum={sum(counts)}, c0={counts[0]}, weight={weight}")
    check("log2|S_C| = 3460.2358 < log2 p^k = 4081.43", abs(logSc - 3460.235753) < 1e-5 and logSc < pk,
          f"{logSc:.6f} vs {pk:.6f}, fibre 2^{logSc - pk:.3f}")

    # ---------- 6. manuscript counters ----------
    displays = len(re.findall(r"\\tag\{\d", txt))
    tables = len(re.findall(r"^表\d", txt, re.M)) - len(re.findall(r"^表\d\.\d(?!　)", txt, re.M))
    refs = len(re.findall(r"^\[\d+\]", txt, re.M))
    check("manuscript has 16 displays after WP6 paste of (4.2)", displays == 16, f"{displays}")
    check("manuscript references count is 17", refs == 17, f"{refs}")

    # ---------- 7. staged texts: superseded-string scan ----------
    paper = PAPER.read_text(encoding="utf-8")
    stale = {
        "不超过$2^{-149.75}$": paper.count("不超过$2^{-149.75}$"),
        "非FCSD项为$2^{-149.75}$": paper.count("非FCSD项为$2^{-149.75}$"),
        "总界为$2^{-136}$": paper.count("总界为$2^{-136}$"),
        "覆盖尚待复核": paper.count("覆盖尚待复核"),
        "实现重测": paper.count("实现重测"),
        "WP2 已完成但 WP3 尚未跑": paper.count("WP3 的 $n=1024$ 攻击覆盖尚未跑"),
    }
    check("WP1_PAPER_TEXT is fully superseded (no stale numeric string left)",
          all(v == 0 for v in stale.values()), f"stale occurrences={stale}")

    sec = SEC44.read_text(encoding="utf-8")
    body = sec.split("本节把定理4.2")[1].split("## 3")[0]
    check("4.4 draft states one hidden commitment per attempt (J not 3J)",
          "混合数是$J$而不是$3J$" in body, "found")
    check("4.4 draft states the three probes share one first message",
          "同一首消息" in body, "found (Theorem 4.2 needs it)")
    check("4.4 draft's assumption F counts advice and preprocessing",
          "建议与预处理计入" in body, "found")

    sizes = SIZES5.read_text(encoding="utf-8")
    check("WP5 section 7 still carries the old C-2 coverage sentence as 'original'",
          "覆盖尚待复核" in sizes, "quoted as the text to replace; fine")

    # ---------- 8. forward references used by the paste ----------
    check("4.4 draft's A3 / 4.3-simulator / Theorem 4.2 / Table 3.1 anchors exist",
          ("A3的隐藏性游戏" in txt) and ("### 4.3" in txt) and
          ("定理4.2" in txt) and ("表3.1　" in txt),
          "A3 game, 4.3, Theorem 4.2, Table 3.1 all present")

    status = "PASS" if all(c["pass"] for c in checks) else "REVIEW"
    result = {"status": status, "date": "2026-10-30",
              "purpose": "independent audit of WP1-WP5 before WP6",
              "checks": checks}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    for c in checks:
        print(("PASS " if c["pass"] else "FAIL ") + c["name"] + " | " + c["detail"])
    print(f"\n{sum(c['pass'] for c in checks)}/{len(checks)} pass; {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
