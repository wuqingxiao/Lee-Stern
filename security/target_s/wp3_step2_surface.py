"""WP3 Step 2: attack-surface pricing at the frozen n=1024 composition.

Three independent blocks:

1. STRUCTURAL.  |S_C|, p^k, the average fibre, and the two blind families
   (enumerate S_C / walk the affine coset), all in exact integer arithmetic.

2. GRID.  Re-run the frozen v7 declared-charge model at n=1024 for the single
   frozen composition, over the same (ell, radius, q) grid as the n=640 final
   round: 100 fixed-E domains and 12 full-union domains.  Where the exact union
   DP fits, the exact cardinality is used; otherwise the fixed-E column is a
   certified subset, hence a lower bound on the union, hence an optimistic
   (attacker-favourable) charge.

3. LATTICE.  Primal Construction A / BKZ in the G1 model, evaluated at n=1024,
   with an explicit check of the model's own domain of validity beta <= n.

Nothing here is executed at target scale, and no total is a measured runtime.
All charges are declared-bit budgets under the frozen v7 conventions.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
SEC = HERE.parent / 'security'
PARAM = HERE.parent / 'parameter_assessment'
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(PARAM))
sys.path.insert(0, str(SEC))
sys.path.insert(0, str(SEC / 'vendor/v61'))

import run_round as rr                       # noqa: E402
import final_model as model                  # noqa: E402
import g2_union_isd as g2                    # noqa: E402
from fixed_domain import FixedDomain         # noqa: E402

P, N, K = 251, 1024, 512
BUDGET_LOG2 = 166.0
LOGE = math.log(2.0)
EXACT_DIR = HERE / 'wp3_exact'
FROZEN = HERE / 'WP3_C_N1024.json'

FIXED_GRID = dict(ells=[2, 8, 16, 32], radii=[2, 3, 4, 6, 8],
                  q=[1, 4, 8, 16, 32], shifts=[-4, 0, 4], dp_cell_limit=500_000)
UNION_GRID = dict(ells=[8, 32], pairs=[[3, 1], [3, 4], [4, 1], [4, 4], [6, 8], [8, 16]],
                  shifts=[-4, 0, 4])
# The frozen grid's cheapest n=1024 point sits at radius 2, the smallest radius
# in the frozen list.  This probe varies ell, radius and q around that point so
# that a boundary artefact cannot be mistaken for the family's minimum.
BOUNDARY = dict(ells=[4, 6, 8, 10, 12, 64], radii=[1, 2, 3], q=[1, 2], shifts=[-4, 0, 4])


def log2_int(value):
    if value <= 0:
        return None
    return math.log2(value)


def lee(a):
    return min(a, P - a)


def load_frozen():
    value = json.loads(FROZEN.read_text(encoding='utf-8'))
    counts = value['counts']
    if len(counts) != P or any(type(x) is not int or x < 0 for x in counts):
        raise ValueError('frozen composition malformed')
    if sum(counts) != N:
        raise ValueError('frozen composition length')
    return value, counts


def cfg_frozen(counts, name):
    return dict(name=name, p=P, n=N, k=K, C=list(counts))


# ---------------------------------------------------------------- structural

def structural(counts):
    pk = P ** K
    size = math.factorial(N)
    for c in counts:
        if c:
            size //= math.factorial(c)
    product = math.factorial(N)
    for c in counts:
        if c:
            product *= math.factorial(c)
    if product % size:                       # |S_C| divides n! exactly
        raise ArithmeticError('n! not divisible by prod c_a!')
    log2_size = (math.lgamma(N + 1) - sum(math.lgamma(c + 1) for c in counts if c)) / LOGE
    if abs(log2_size - log2_int(size)) > 1e-9:
        raise ArithmeticError('lgamma cardinality disagrees with the integer count')
    lee_weight = sum(c * lee(a) for a, c in enumerate(counts))
    l2 = math.sqrt(sum(c * lee(a) ** 2 for a, c in enumerate(counts)))
    nonzero = sum(1 for c in counts if c)
    # Two exact linear relations leaked by the public composition: sum e_j and
    # sum e_j^2 are both determined by the histogram (the latter modulo p).
    walk_dim = (N - K) * math.log2(P)
    walk_dim_leaked = (N - K - 2) * math.log2(P)
    return dict(
        size_exact=str(size), log2_size=log2_size, log2_pk=(N - K) * math.log2(P),
        pk_exact=str(pk), size_below_pk=size < pk,
        log2_fibre=log2_size - (N - K) * math.log2(P),
        lee_weight=lee_weight, l2_norm=l2, distinct_symbols=nonzero,
        zeros=counts[0],
        brute_force_log2=log2_size,
        coset_walk_log2=walk_dim,
        coset_walk_leaked_log2=walk_dim_leaked,
        leaked_linear_relations=2,
        leaked_gain_bits=2 * math.log2(P))


# ------------------------------------------------------------------- lattice

def delta0(beta):
    return (beta / (2 * math.pi * math.e)) ** (1 / (2 * beta))


class Lattice:
    def __init__(self, n, k, p, target):
        self.n, self.k, self.p, self.target = n, k, p, target
        self.det = p ** (n - k)
        self.det_root = math.exp(math.log(self.det) / n)
        self.gh = self.det_root * math.sqrt(n / (2 * math.pi * math.e))
        if abs(math.log(self.det) / n - 0.5 * math.log(p)) > 1e-12:
            raise ArithmeticError('determinant root is not p^(1/2) for k=n/2')

    def norm(self, beta):
        return delta0(beta) ** self.n * self.det_root

    def beta_required(self, lo=2.0, hi=1.0e5):
        if self.norm(hi) > self.target:
            return None
        for _ in range(200):
            mid = math.sqrt(lo * hi)
            if self.norm(mid) > self.target:
                lo = mid
            else:
                hi = mid
        return math.sqrt(lo * hi)

    def report(self, budget_log2=BUDGET_LOG2, c=0.292):
        beta_max = budget_log2 / c
        beta_req = self.beta_required()
        return dict(
            n=self.n, k=self.k, det_log2=self.n * math.log2(self.det_root),
            det_root=self.det_root, gh=self.gh, target_norm=self.target,
            target_over_gh=self.target / self.gh,
            beta_at_budget=beta_max, norm_at_budget=self.norm(beta_max),
            budget_reaches_target=self.norm(beta_max) <= self.target,
            beta_required=beta_req,
            beta_required_exceeds_n=(beta_req is not None and beta_req > self.n),
            required_cost_log2=(None if beta_req is None else c * beta_req),
            norm_at_beta_n=self.norm(self.n),
            gh_identity_error=abs(self.norm(self.n) - self.gh) / self.gh,
            floor_cost_log2=c * self.n,
            floor_margin_bits=c * self.n - budget_log2,
            note='beta>n is outside the root-Hermite model; the reported floor '
                 'is the cheapest charge the model can express with beta<=n')


# ---------------------------------------------------------------------- grid

def exact_path(ell, radius, q):
    return EXACT_DIR / f'e{ell}_r{radius}_q{q}.json'


def load_exact():
    out = {}
    if EXACT_DIR.exists():
        for path in EXACT_DIR.glob('*.json'):
            value = json.loads(path.read_text(encoding='utf-8'))
            if value.get('status') == 'EXACT':
                out[(value['ell'], value['radius'], value['q'])] = value
    return out


def run_exact_one(ell, radius, q, state_cap=4_000_000, time_cap=240.0):
    ts = rr.variants(cfg_frozen(load_frozen()[1], 'probe'), ell, radius, q,
                     'full_union', UNION_GRID['shifts'])
    if not ts:
        raise ValueError('no template')
    h = (K + ell) // 2
    a = ts[0]['a']
    value, elapsed = g2.exact_union_size(P, h, a, q, radius,
                                         state_cap=state_cap, time_cap=time_cap, verbose=True)
    out = dict(ell=ell, radius=radius, q=q, h=h, a=a, status='EXACT',
               size=str(value), log2=math.log2(value), seconds=round(elapsed, 2))
    EXACT_DIR.mkdir(exist_ok=True)
    exact_path(ell, radius, q).write_text(json.dumps(out, ensure_ascii=False, indent=2),
                                          encoding='utf-8')
    print(json.dumps(out, ensure_ascii=False), flush=True)


def fixed_grid(cfg, spec=None):
    spec = spec or FIXED_GRID
    rows_out, coverage = [], []
    for ell in spec['ells']:
        for radius in spec['radii']:
            for q in spec['q']:
                case = [ell, radius, q]
                ts = rr.variants(cfg, ell, radius, q, 'fixed_E', spec['shifts'])
                if not ts:
                    coverage.append(dict(family='fixed_E', case=case,
                                         status='NO_TEMPLATE_IN_THIS_SEED_RULE'))
                    continue
                h = (K + ell) // 2
                try:
                    domain = FixedDomain(P, h, ts[0]['a'], q, radius,
                                         max_cells=spec.get('dp_cell_limit', 500_000))
                    rows, selected, rejected = rr.evaluate_domain(cfg, ts, domain, g2.CONFIG)
                except ValueError as exc:
                    coverage.append(dict(family='fixed_E', case=case,
                                         status='ANALYSIS_LIMIT', reason=str(exc)[:160]))
                    continue
                if not rows:
                    coverage.append(dict(family='fixed_E', case=case, status='NO_ROW',
                                         reason=str(dict(rejected))[:200]))
                    continue
                best = min(rows, key=lambda r: r['declared_bit_upper_log2'])
                bh = min(rows, key=lambda r: r['heuristic_log2'])
                rows_out.append(dict(family='fixed_E', ell=ell, radius=radius, q=q, h=h,
                                     a=ts[0]['a'], domain_log2=math.log2(domain.size),
                                     size_used='EXACT_SUBSET', points=len(rows),
                                     templates=len(ts),
                                     best_declared_log2=best['declared_bit_upper_log2'],
                                     best_heuristic_log2=bh['heuristic_log2'],
                                     best_parameters={k: best[k] for k in
                                                      ('u', 'L', 'block_records', 'tile_records',
                                                       'schedule', 'asymmetric')}))
                coverage.append(dict(family='fixed_E', case=case, status='EVALUATED',
                                     points=len(rows)))
    return rows_out, coverage


def union_grid(cfg, exact):
    rows_out, coverage = [], []
    for ell in UNION_GRID['ells']:
        h = (K + ell) // 2
        for radius, q in UNION_GRID['pairs']:
            ts = rr.variants(cfg, ell, radius, q, 'full_union', UNION_GRID['shifts'])
            if not ts:
                coverage.append(dict(family='full_union', case=[ell, radius, q],
                                     status='NO_TEMPLATE_IN_THIS_SEED_RULE'))
                continue
            bound = g2.union_lower_bound(P, h, ts[0]['a'], q, radius)
            key = (ell, radius, q)
            if key in exact:
                size, used, exact_log2 = int(exact[key]['size']), 'EXACT', exact[key]['log2']
                if size < bound:
                    raise ArithmeticError('exact union below the certified subset bound')
            else:
                size, used, exact_log2 = bound, 'LOWER_BOUND', None
            out = g2.evaluate(cfg, ell, radius, q, size)
            if not out or out.get('status') != 'EVALUATED':
                coverage.append(dict(family='full_union', case=[ell, radius, q],
                                     status='NO_ROW' if out else 'NOT_EVALUATED'))
                continue
            row = dict(family='full_union', ell=ell, radius=radius, q=q, h=h, a=ts[0]['a'],
                       domain_log2=exact_log2 if used == 'EXACT' else math.log2(bound),
                       subset_lower_bound_log2=math.log2(bound),
                       size_used=used, points=out['rows'], templates=out['templates'],
                       best_declared_log2=out['best']['declared_bit_upper_log2'],
                       best_heuristic_log2=out['best_heuristic']['heuristic_log2'],
                       best_parameters={k: out['best'][k] for k in
                                        ('u', 'L', 'block_records', 'tile_records',
                                         'schedule', 'asymmetric')})
            rows_out.append(row)
            coverage.append(dict(family='full_union', case=[ell, radius, q],
                                 status='EVALUATED', points=out['rows'], size_used=used))
    return rows_out, coverage


def evaluate_grid(cfg):
    original, patched = g2._patched_layout()
    model.old.memory_layout = patched
    try:
        exact = load_exact()
        fixed_rows, fixed_cov = fixed_grid(cfg)
        union_rows, union_cov = union_grid(cfg, exact)
        boundary_rows, boundary_cov = fixed_grid(cfg, BOUNDARY)
    finally:
        model.old.memory_layout = original
    return dict(fixed=fixed_rows, union=union_rows, boundary=boundary_rows,
                coverage=fixed_cov + union_cov,
                boundary_coverage=boundary_cov,
                exact_available=sorted([list(k) for k in exact]))


# ---------------------------------------------------------------------- main

def closed_form_subset(p, h, a, q, r):
    """Independent count of one fixed-E domain: p^q * [x^a](1+2 sum_{w<=r} x^w)^(h-q)."""
    m = h - q
    if a < 0 or a > m * r:
        return 0
    dp = [0] * (a + 1)
    dp[0] = 1
    for _ in range(m):
        nd = [0] * (a + 1)
        for w, v in enumerate(dp):
            if not v:
                continue
            nd[w] += v
            for step in range(1, r + 1):
                if w + step <= a:
                    nd[w + step] += 2 * v
        dp = nd
    return p ** q * dp[a]


def validate_subset_sizes(rows):
    """Recompute every fixed-E domain size from the closed form.

    For a full-union row the closed form reproduces the certified subset bound,
    not the union size; the two coincide only when the union rows fall back to
    the lower bound.
    """
    out = []
    for row in rows:
        value = closed_form_subset(P, row['h'], row['a'], row['q'], row['radius'])
        target = row.get('subset_lower_bound_log2', row['domain_log2'])
        out.append(dict(case=[row['ell'], row['radius'], row['q']], h=row['h'], a=row['a'],
                        closed_log2=math.log2(value) if value else None,
                        compared_to=target,
                        match=(abs(math.log2(value) - target) < 1e-9
                               if value else False)))
    return out


def write_csv(result, path):
    import csv
    rows = []
    for family in ('fixed', 'union', 'boundary'):
        for r in result['grid'][family]:
            rows.append(dict(
                family=r['family'] if family != 'boundary' else 'fixed_E_boundary',
                ell=r['ell'], radius=r['radius'], q=r['q'], h=r['h'], a=r['a'],
                domain_log2=round(r['domain_log2'], 3),
                size_used=r.get('size_used', 'EXACT_SUBSET'),
                best_declared_log2=round(r['best_declared_log2'], 3),
                best_heuristic_log2=round(r['best_heuristic_log2'], 3),
                margin_over_2_166=round(r['best_declared_log2'] - BUDGET_LOG2, 3)))
    if not rows:
        return
    with path.open('w', encoding='utf-8-sig', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--exact-one', nargs=3, metavar=('ELL', 'R', 'Q'))
    ap.add_argument('--exact-attempt', action='store_true')
    ap.add_argument('--out', type=Path, default=HERE / 'review/WP3_STEP2.json')
    args = ap.parse_args()

    if args.exact_one:
        ell, radius, q = (int(x) for x in args.exact_one)
        run_exact_one(ell, radius, q)
        return
    if args.exact_attempt:
        for ell in UNION_GRID['ells']:
            for radius, q in UNION_GRID['pairs']:
                if (ell, radius, q) in load_exact():
                    continue
                try:
                    run_exact_one(ell, radius, q, time_cap=200.0)
                except g2.UnionStateLimit as exc:
                    print(json.dumps(dict(ell=ell, radius=radius, q=q,
                                          status='STATE_OR_TIME_CAP', reason=str(exc)[:90])),
                          flush=True)
                except Exception as exc:          # noqa: BLE001 - reported, not hidden
                    print(json.dumps(dict(ell=ell, radius=radius, q=q,
                                          status='FAILED', reason=repr(exc)[:120])), flush=True)
        return

    started = time.perf_counter()
    frozen, counts = load_frozen()
    cfg = cfg_frozen(counts, frozen['name'])
    result = dict(
        scope='n=1024 frozen composition; declared-charge budgets and exact '
              'combinatorics, not measured runtimes and not a proof of assumption F',
        budget_log2=BUDGET_LOG2,
        frozen_composition=frozen['name'],
        frozen_sha256=hashlib.sha256(FROZEN.read_bytes()).hexdigest(),
        input_sha256=hashlib.sha256((SEC / 'inputs/candidates.json').read_bytes()).hexdigest(),
        config_sha256=hashlib.sha256((SEC / 'config.json').read_bytes()).hexdigest())

    st = structural(counts)
    lat = Lattice(N, K, P, st['l2_norm']).report()
    grid = evaluate_grid(cfg)
    result.update(structural=st, lattice=lat, grid=grid)
    subset_check = validate_subset_sizes(grid['fixed'] + grid['union'] + grid['boundary'])
    result['subset_size_validation'] = subset_check

    declared = [r['best_declared_log2'] for r in grid['fixed'] + grid['union']]
    heuristic = [r['best_heuristic_log2'] for r in grid['fixed'] + grid['union']]
    b_declared = [r['best_declared_log2'] for r in grid['boundary']]
    b_heuristic = [r['best_heuristic_log2'] for r in grid['boundary']]

    def point(row):
        return dict((k, row[k]) for k in ('ell', 'radius', 'q', 'h', 'a', 'domain_log2',
                                          'size_used', 'best_declared_log2',
                                          'best_heuristic_log2'))

    best_row = min(grid['fixed'] + grid['union'],
                   key=lambda r: r['best_declared_log2']) if declared else None
    all_rows = grid['fixed'] + grid['union'] + grid['boundary']
    overall = min(all_rows, key=lambda r: r['best_declared_log2']) if all_rows else None
    result['summary'] = dict(
        fixed_domains=len(grid['fixed']), union_domains=len(grid['union']),
        points_total=sum(r['points'] for r in grid['fixed'] + grid['union']),
        min_declared_log2=min(declared) if declared else None,
        min_heuristic_log2=min(heuristic) if heuristic else None,
        min_declared_family=best_row['family'] if best_row else None,
        min_declared_point=point(best_row) if best_row else None,
        points_at_or_below_budget=sum(1 for x in declared if x <= BUDGET_LOG2),
        margin_min_declared_bits=(min(declared) - BUDGET_LOG2) if declared else None,
        proxy_success_at_budget_log2=(math.log2(math.log(2)) + BUDGET_LOG2
                                      - min(heuristic)) if declared else None,
        exact_union_domains=sum(1 for r in grid['union'] if r['size_used'] == 'EXACT'),
        lower_bound_union_domains=sum(1 for r in grid['union']
                                      if r['size_used'] == 'LOWER_BOUND'),
        boundary_domains=len(grid['boundary']),
        boundary_min_declared_log2=min(b_declared) if b_declared else None,
        boundary_min_heuristic_log2=min(b_heuristic) if b_heuristic else None,
        boundary_points_at_or_below_budget=sum(1 for x in b_declared if x <= BUDGET_LOG2),
        overall_min_declared_log2=overall['best_declared_log2'] if overall else None,
        overall_min_point=point(overall) if overall else None,
        overall_margin_bits=(overall['best_declared_log2'] - BUDGET_LOG2) if overall else None,
        elapsed_seconds=round(time.perf_counter() - started, 1))

    checks = []

    def check(name, ok, detail=''):
        checks.append(dict(name=name, ok=bool(ok), detail=str(detail)))

    check('composition sums to n', sum(counts) == N, sum(counts))
    check('S_C below p^k (unique preimage)', st['size_below_pk'], st['log2_fibre'])
    check('fibre equals size minus pk', abs(st['log2_fibre'] -
                                            (st['log2_size'] - st['log2_pk'])) < 1e-9)
    check('lgamma log2 matches integer count',
          abs(st['log2_size'] - math.log2(int(st['size_exact']))) < 1e-9)
    check('brute force equals log2|S_C|', abs(st['brute_force_log2'] - st['log2_size']) < 1e-12)
    check('coset walk equals k log2 p', abs(st['coset_walk_log2'] - K * math.log2(P)) < 1e-12)
    check('brute force cheaper than coset walk',
          st['brute_force_log2'] < st['coset_walk_log2'])
    check('leaked relations do not close the gap',
          st['coset_walk_leaked_log2'] > st['brute_force_log2'])
    check('lattice GH identity at beta=n', lat['gh_identity_error'] < 1e-12,
          lat['gh_identity_error'])
    check('target below Gaussian heuristic', lat['target_over_gh'] < 1,
          lat['target_over_gh'])
    check('required beta exceeds dimension', lat['beta_required_exceeds_n'],
          lat['beta_required'])
    check('budget cannot reach target in model', not lat['budget_reaches_target'],
          lat['norm_at_budget'])
    check('model floor above budget', lat['floor_cost_log2'] > BUDGET_LOG2 + 100,
          lat['floor_cost_log2'])
    check('every union domain evaluated', len(grid['union']) == 12, len(grid['union']))
    check('union size accounting complete',
          result['summary']['exact_union_domains']
          + result['summary']['lower_bound_union_domains'] == len(grid['union']))
    check('no grid point within budget', result['summary']['points_at_or_below_budget'] == 0,
          result['summary']['min_declared_log2'])
    check('grid minimum above 2^166', (result['summary']['min_declared_log2'] or 0) > BUDGET_LOG2,
          result['summary']['margin_min_declared_bits'])
    check('fixed-E grid complete', len(grid['fixed']) == 100, len(grid['fixed']))
    check('boundary probe accounted for', len(grid['boundary_coverage'])
          == len(BOUNDARY['ells']) * len(BOUNDARY['radii']) * len(BOUNDARY['q']),
          len(grid['boundary_coverage']))
    check('boundary probe has no analysis limits',
          all(c['status'] in ('EVALUATED', 'NO_TEMPLATE_IN_THIS_SEED_RULE')
              for c in grid['boundary_coverage'])
          and len(grid['boundary']) >= 30, len(grid['boundary']))
    check('boundary probe stays above budget',
          result['summary']['boundary_points_at_or_below_budget'] == 0,
          result['summary']['boundary_min_declared_log2'])
    check('cheapest priced family beats lattice floor?', st['brute_force_log2'] > lat['floor_cost_log2'],
          [st['brute_force_log2'], lat['floor_cost_log2']])
    check('every fixed-E domain size matches the closed form',
          all(x['match'] for x in subset_check), len(subset_check))
    result['checks'] = checks
    result['status'] = 'PASS' if all(c['ok'] for c in checks) else 'FAIL'
    result['passed'] = sum(c['ok'] for c in checks)
    result['total_checks'] = len(checks)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
    write_csv(result, HERE / 'WP3_STEP2_TABLE.csv')
    print(json.dumps(dict(status=result['status'], passed=result['passed'],
                          total=result['total_checks'],
                          summary=result['summary']), ensure_ascii=False, indent=1))
    bad = [c for c in checks if not c['ok']]
    for c in bad:
        print('FAILED:', c['name'], c['detail'])


if __name__ == '__main__':
    main()
