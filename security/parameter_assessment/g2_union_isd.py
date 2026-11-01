"""G2: finish the full_union coverage grid and calibrate the restricted-ball family at 2^166.

Two independent computations:

1. exact_union_size: exact cardinality of the full permutation-closure base set
       U = {v in F_p^h : exists E subset [h], |E|=q,
            v_j low (|v_j|_L <= r) for j not in E, sum_{j not in E} |v_j|_L = a}.
   The v7 vendor.Domain enumerates one histogram per accepted group and stops at
   an arbitrary candidate cap.  Here the same partition is evaluated by a DP over
   the low Lee categories (0 = zero, w = +-w for 1<=w<=r) and the high category
   (|v_j|_L > r).  The state carries the number of low coordinates, the low
   weight, and the exact set of achievable removals (size t, weight delta) as a
   bitset.  The base set is counted once per vector, not once per mask.

2. optimistic_charge: for a domain whose exact size is out of reach, the fixed-E
   subdomain B_E (one prescribed free set) is a subset of U.  Its size is an
   exact, certified lower bound |U| >= |B_E|.  The v7 survival is decreasing in
   |U|, hence the charged budget computed with |B_E| is a lower bound on the
   charged budget of the true union domain.  Replacing the union metadata by a
   single group and dropping the union histogram cache reservation only makes
   the surrogate more optimistic (allows larger lists).  Every surviving point
   is therefore a rigorous lower bound on the declared v7 charge.

The script never executes a target-scale attack and never claims a security
level.  All totals are declared-charge budgets, not measured runtimes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from fractions import Fraction as F
from math import comb
from pathlib import Path
from types import SimpleNamespace

HERE = Path(__file__).resolve().parent
SEC = HERE.parent / 'security'
sys.path.insert(0, str(SEC))
sys.path.insert(0, str(SEC / 'vendor/v61/vendor'))

import final_model as model          # noqa: E402
import run_round as rr               # noqa: E402
from fixed_domain import FixedDomain  # noqa: E402

BUDGET_LOG2 = 166.0
CANDIDATES = json.loads((SEC / 'inputs/candidates.json').read_text(encoding='utf-8'))
CONFIG = json.loads((SEC / 'config.json').read_text(encoding='utf-8'))
POINTS = SEC / 'final_round/points.csv'
COVERAGE = json.loads((SEC / 'final_round/coverage.json').read_text(encoding='utf-8'))


class UnionStateLimit(RuntimeError):
    pass


def _row_mask(q, r, limit):
    d = q * r + 1
    m = 0
    for t in range(q + 1):
        for delta in range(0, min(limit, q * r) + 1):
            m |= 1 << (t * d + delta)
    return m


def exact_union_size(p, h, a, q, r, *, state_cap=6_000_000, time_cap=600.0, verbose=False):
    """Exact |U| by category DP with an exact removal-achievement bitset.

    Validated against vendor.Domain on every completed v7 full_union domain and
    against brute-force enumeration on small parameter sets (see validate()).
    """
    if not (p >= 3 and 1 <= q <= h and 0 <= a and 1 <= r <= p // 2):
        raise ValueError('union domain parameters')
    stride = q * r + 1
    idx = lambda t, d: t * stride + d
    mask = (1 << ((q + 1) * stride)) - 1
    limits = {w: _row_mask(q, r, q * r - w) for w in range(r + 1)}
    # (low_count, low_weight, removal_bitset) -> number of labelled words on the
    # h positions whose non-low positions are still unassigned (all high later).
    dp = {(0, 0, 1): 1}
    started = time.perf_counter()
    for w in range(r + 1):
        alphabet = 1 if w == 0 else 2
        ndp = {}
        seen = 0
        for (used, weight, bits), count in dp.items():
            seen += 1
            if (seen & 0x3FFFF) == 0 and (len(ndp) > state_cap or time.perf_counter() - started > time_cap):
                raise UnionStateLimit(f'state/time cap while expanding {len(ndp)} states (category {w})')
            remaining = h - used
            # After all low categories the remaining positions are high and must
            # not exceed q; the smallest possible final low weight is `weight`.
            for c in range(remaining + 1):
                new_weight = weight + w * c
                if new_weight > a + q * r:
                    break
                new_bits = bits
                for d in range(1, min(c, q) + 1):
                    new_bits |= (bits & limits[w]) << idx(d, w * d)
                new_bits &= mask
                ways = comb(remaining, c) * alphabet ** c
                key = (used + c, new_weight, new_bits)
                ndp[key] = ndp.get(key, 0) + count * ways
        dp = ndp
        if verbose:
            print(f'  category {w}: states={len(dp)} t={time.perf_counter()-started:.1f}s', flush=True)
        if len(dp) > state_cap or time.perf_counter() - started > time_cap:
            raise UnionStateLimit(f'state cap: {len(dp)} states after category {w}')
    total = 0
    for (used, weight, bits), count in dp.items():
        high = h - used
        if high > q:
            continue
        t = q - high
        delta = weight - a
        if delta < 0 or delta > q * r:
            continue
        if (bits >> idx(t, delta)) & 1:
            total += count * (p - 1 - 2 * r) ** high
    return total, time.perf_counter() - started


def brute_union_size(p, h, a, q, r):
    """Independent enumeration for tiny cases only."""
    import itertools
    total = 0
    for v in itertools.product(range(p), repeat=h):
        low = [min(x, p - x) for x in v]
        high = sum(1 for w in low if w > r)
        if high > q:
            continue
        weight = sum(w for w in low if w <= r)
        delta = weight - a
        if delta < 0:
            continue
        achievable = {(0, 0)}
        for w in (x for x in low if x <= r):
            achievable |= {(c + 1, d + w) for (c, d) in list(achievable)}
        if (q - high, delta) in achievable:
            total += 1
    return total


def union_lower_bound(p, h, a, q, r):
    """|B_E| for one fixed free set; exact subset of U."""
    return FixedDomain(p, h, a, q, r, max_cells=10 ** 7).size


def _patched_layout():
    """Keep the v7 list/rank/PGE areas; drop only union histogram metadata.

    The union histogram table and its worst-case cache formula are setup
    reservations of the vendor.Domain enumerator.  A domain of the same size can
    be stored with one explicit group list; dropping these areas permits larger
    L, i.e. it is deliberately favourable to the attacker.
    """
    original = model.old.memory_layout

    def layout(p, n, k, ell, u, B, L, b, domain):
        m = original(p, n, k, ell, u, B, L, b, domain)
        areas = m['areas']
        for key in ('group_prefix_metadata', 'setup_removal_cache_reservation',
                    'setup_factorial_power_tables', 'setup_DP_buffers'):
            areas[key] = 8
        m['reserved_bytes'] = sum(areas.values())
        m['domain_metadata_model'] = 'G2 optimistic: single-group metadata, union cache not charged'
        return m

    return original, layout


def make_stub(p, h, a, q, r, size):
    return SimpleNamespace(Q=p, h=h, a=a, q=q, radius=r, size=size,
                           groups=(None,), stats={}, low=h - q,
                           setup_operand_max_bits=size.bit_length() + h.bit_length())


def evaluate(cfg, ell, radius, q, size, *, optimistic_metadata=True):
    ts = rr.variants(cfg, ell, radius, q, 'full_union', CONFIG['full_union']['asymmetry_shifts'])
    if not ts:
        return None
    domain = make_stub(cfg['p'], (cfg['k'] + ell) // 2, ts[0]['a'], q, radius, size)
    rows, selected, rejected = rr.evaluate_domain(cfg, ts, domain, CONFIG)
    if not rows:
        return dict(status='NO_ROW', rejected=dict(rejected), templates=len(ts), a=ts[0]['a'])
    best = min(rows, key=lambda x: x['declared_bit_upper_log2'])
    best_h = min(rows, key=lambda x: x['heuristic_log2'])
    return dict(status='EVALUATED', templates=len(ts), a=ts[0]['a'], rows=len(rows),
                best=best, best_heuristic=best_h, rejected=dict(rejected))


def union_cases():
    """Every full_union (candidate, ell, r, q) with its v7 coverage status."""
    out = []
    for entry in COVERAGE:
        if entry['family'] != 'full_union':
            continue
        out.append(dict(candidate=entry['candidate'],
                        ell=entry['case'][0], radius=entry['case'][1], q=entry['case'][2],
                        coverage=entry['status']))
    return out


def load_points():
    import csv
    with POINTS.open(encoding='utf-8-sig') as fh:
        return list(csv.DictReader(fh))


def validate(verbose=True):
    """Exact DP vs vendor.Domain on all completed union domains, and brute force."""
    report = {'small': [], 'domains': []}
    small = [(7, 6, 4, 1, 2), (7, 6, 4, 2, 2), (7, 5, 3, 2, 2),
             (11, 4, 3, 1, 2), (11, 5, 5, 1, 2)]
    for (p, h, a, q, r) in small:
        brute = brute_union_size(p, h, a, q, r)
        dp, _ = exact_union_size(p, h, a, q, r)
        report['small'].append(dict(case=[p, h, a, q, r], brute=brute, dp=dp, match=brute == dp))
    from domain import Domain as VendorDomain
    for entry in union_cases():
        if entry['coverage'] != 'DOMAIN_GRID_COMPLETE':
            continue
        cfg = next(c for c in CANDIDATES if c['name'] == entry['candidate'])
        ts = rr.variants(cfg, entry['ell'], entry['radius'], entry['q'], 'full_union',
                         CONFIG['full_union']['asymmetry_shifts'])
        h = (cfg['k'] + entry['ell']) // 2
        a = ts[0]['a']
        domain = VendorDomain(cfg['p'], h, a, entry['q'], entry['radius'], max_candidates=50000)
        dp, elapsed = exact_union_size(cfg['p'], h, a, entry['q'], entry['radius'])
        report['domains'].append(dict(candidate=entry['candidate'][-12:], ell=entry['ell'],
                                      radius=entry['radius'], q=entry['q'], a=a,
                                      vendor_log2=math.log2(domain.size), dp_log2=math.log2(dp),
                                      exact_match=dp == domain.size, seconds=round(elapsed, 1)))
        if verbose:
            print(report['domains'][-1], flush=True)
    report['all_small_match'] = all(x['match'] for x in report['small'])
    report['all_domains_match'] = all(x['exact_match'] for x in report['domains'])
    return report


EXACT_DIR = HERE / 'g2_exact'


def case_parameters(cfg, ell, radius, q):
    ts = rr.variants(cfg, ell, radius, q, 'full_union', CONFIG['full_union']['asymmetry_shifts'])
    h = (cfg['k'] + ell) // 2
    return h, ts[0]['a']


def exact_path(candidate, ell, radius, q):
    return EXACT_DIR / f"{candidate[-12:]}_e{ell}_r{radius}_q{q}.json"


def load_exact():
    out = {}
    if not EXACT_DIR.exists():
        return out
    for path in EXACT_DIR.glob('*.json'):
        value = json.loads(path.read_text(encoding='utf-8'))
        out[(value['candidate'], value['ell'], value['radius'], value['q'])] = value
    return out


def run_exact_one(candidate_name, ell, radius, q):
    cfg = next(c for c in CANDIDATES if c['name'].endswith(candidate_name[-12:]))
    h, a = case_parameters(cfg, ell, radius, q)
    value, elapsed = exact_union_size(cfg['p'], h, a, q, radius,
                                      state_cap=5_000_000, time_cap=300.0, verbose=True)
    out = dict(candidate=cfg['name'], ell=ell, radius=radius, q=q, a=a,
               status='EXACT', size=str(value), log2=math.log2(value),
               seconds=round(elapsed, 2))
    EXACT_DIR.mkdir(exist_ok=True)
    exact_path(cfg['name'], ell, radius, q).write_text(
        json.dumps(out, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(out, ensure_ascii=False), flush=True)


def feasible_exact_cases():
    """Cases where the bitset DP is expected to fit the state/time budget."""
    out = []
    for entry in union_cases():
        if entry['radius'] <= 4:
            out.append((entry['candidate'][-12:], entry['ell'], entry['radius'], entry['q']))
    return out


def isd_metrics():
    """Heuristic success at 2^166 from the saved v7 grid (Poisson proxy only)."""
    rows = load_points()
    per_candidate = {}
    for cfg in CANDIDATES:
        mine = [r for r in rows if r['candidate'] == cfg['name']]
        hbest = min(mine, key=lambda r: float(r['heuristic_log2']))
        cbest = min(mine, key=lambda r: float(r['declared_bit_upper_log2']))
        H = float(hbest['heuristic_log2'])
        per_candidate[cfg['name'][-12:]] = dict(
            heuristic_best_log2=H,
            heuristic_best_point={k: hbest[k] for k in
                                  ('family', 'ell', 'radius', 'q', 'a', 'u', 'L', 'schedule',
                                   'heuristic_log2', 'declared_bit_upper_log2')},
            declared_best_log2=float(cbest['declared_bit_upper_log2']),
            declared_best_point={k: cbest[k] for k in
                                 ('family', 'ell', 'radius', 'q', 'a', 'u', 'L', 'schedule',
                                  'heuristic_log2', 'declared_bit_upper_log2')},
            heuristic_success_at_166_log2=math.log2(math.log(2)) + BUDGET_LOG2 - H,
            declared_gap_bits=float(cbest['declared_bit_upper_log2']) - BUDGET_LOG2)
    H_all = min(v['heuristic_best_log2'] for v in per_candidate.values())
    C_all = min(v['declared_best_log2'] for v in per_candidate.values())
    return dict(per_candidate=per_candidate,
                family_min_heuristic_log2=H_all,
                family_min_declared_log2=C_all,
                family_heuristic_success_at_166_log2=math.log2(math.log(2)) + BUDGET_LOG2 - H_all,
                family_declared_gap_bits=C_all - BUDGET_LOG2,
                note='Poisson/log-normal proxy of the same heuristic columns; '
                     'not an upper bound on FCSD success and not a theorem')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, default=HERE / 'g2_union_isd_results.json')
    ap.add_argument('--validate', action='store_true')
    ap.add_argument('--exact-list', action='store_true')
    ap.add_argument('--exact-one', nargs=4, metavar=('CANDIDATE', 'ELL', 'R', 'Q'))
    ap.add_argument('--opt', action='store_true', help='evaluate all union cases with exact/lower-bound size')
    args = ap.parse_args()

    if args.exact_list:
        print(json.dumps(feasible_exact_cases(), ensure_ascii=False, indent=1))
        return

    if args.exact_one:
        candidate, ell, radius, q = args.exact_one
        run_exact_one(candidate, int(ell), int(radius), int(q))
        return

    started = time.perf_counter()
    result = {'scope': 'restricted-ball / PGE family on the fixed v7 grid; '
                       'declared-charge budgets, not measured runtimes',
              'budget_log2': BUDGET_LOG2,
              'input_sha256': hashlib.sha256((SEC / 'inputs/candidates.json').read_bytes()).hexdigest(),
              'config_sha256': hashlib.sha256((SEC / 'config.json').read_bytes()).hexdigest()}

    if args.validate:
        result['validation'] = validate()
        args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
        print(json.dumps(result['validation'], ensure_ascii=False, indent=1))
        return

    if not args.opt:
        ap.error('choose --validate, --exact-list, --exact-one or --opt')

    original_mem, opt_mem = _patched_layout()
    model.old.memory_layout = opt_mem
    exact = load_exact()

    checks = []
    for entry in union_cases():
        cfg = next(c for c in CANDIDATES if c['name'] == entry['candidate'])
        key = (entry['candidate'], entry['ell'], entry['radius'], entry['q'])
        h, a = case_parameters(cfg, entry['ell'], entry['radius'], entry['q'])
        bound = union_lower_bound(cfg['p'], h, a, entry['q'], entry['radius'])
        row = dict(candidate=entry['candidate'][-12:], ell=entry['ell'], radius=entry['radius'],
                   q=entry['q'], a=a, coverage=entry['coverage'],
                   lower_bound_log2=math.log2(bound), lower_bound=str(bound))
        if key in exact and exact[key]['status'] == 'EXACT':
            row['size_used'] = 'EXACT'
            row['size_log2'] = exact[key]['log2']
            row['exact'] = exact[key]
            size = int(exact[key]['size'])
        else:
            row['size_used'] = 'LOWER_BOUND'
            row['size_log2'] = math.log2(bound)
            size = bound
        if row['size_log2'] + 1e-9 < row['lower_bound_log2']:
            raise ArithmeticError('exact size below certified subset bound')
        out = evaluate(cfg, entry['ell'], entry['radius'], entry['q'], size)
        row.update(out or {})
        checks.append(row)
        print('eval', entry['candidate'][-12:], entry['ell'], entry['radius'], entry['q'],
              row['size_used'], round(row['size_log2'], 2),
              row.get('best', {}).get('declared_bit_upper_log2'), flush=True)

    result['exact'] = {f'{k[0][-12:]}|{k[1]}|{k[2]}|{k[3]}': v for k, v in exact.items()}
    result['union_checks'] = checks
    result['isd_metrics'] = isd_metrics()
    result['elapsed_seconds'] = time.perf_counter() - started
    result['summary'] = dict(
        evaluated=sum(1 for r in checks if r.get('status') == 'EVALUATED'),
        total=len(checks),
        exact=sum(1 for r in checks if r['size_used'] == 'EXACT'),
        lower_bound=sum(1 for r in checks if r['size_used'] == 'LOWER_BOUND'),
        min_optimistic_declared_log2=min(
            (r['best']['declared_bit_upper_log2'] for r in checks if r.get('status') == 'EVALUATED'),
            default=None),
        min_optimistic_heuristic_log2=min(
            (r['best_heuristic']['heuristic_log2'] for r in checks if r.get('status') == 'EVALUATED'),
            default=None))
    args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(result['summary'], ensure_ascii=False, indent=1))
    print(json.dumps(result['isd_metrics'], ensure_ascii=False, indent=1))


if __name__ == '__main__':
    main()
