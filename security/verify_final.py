"""Replay saved certificates and verify the finite-grid verdict, not hardness."""
import argparse
from collections import Counter
import csv
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import final_model as model
from fixed_domain import FixedDomain

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify(out):
    summary = read(out/'SUMMARY.json')
    assert summary['status'] == 'DECLARED_FINAL_ROUND_FINISHED'
    for name, value in summary['code_binding'].items():
        assert sha(HERE/name) == value, name
    for name, value in read(HERE/'SOURCE_BINDINGS.json').items():
        assert sha(HERE/name) == value['sha256'], name
    assert sha(HERE/'inputs/candidates.json') == summary['input_sha256']
    assert sha(HERE/'config.json') == summary['config_sha256']
    assert read(out/'config.json') == read(HERE/'config.json')
    cfgs = {x['name']: x for x in read(HERE/'inputs/candidates.json')}
    assert read(out/'candidates.json') == list(cfgs.values())
    for cfg in cfgs.values():
        C, p = cfg['C'], cfg['p']
        assert sum(C) == cfg['n'] and len(C) == p and min(C) >= 0
        assert sum(min(i,p-i)*c for i,c in enumerate(C)) == cfg['weight']
        assert hashlib.sha256(json.dumps(C,separators=(',',':')).encode()).hexdigest() == cfg['C_sha256']

    with (out/'points.csv').open(encoding='utf-8-sig',newline='') as f:
        points = list(csv.DictReader(f))
    assert len(points) == summary['points'] == 35802
    families = Counter(x['family'] for x in points)
    assert families == {'fixed_E': 34650, 'full_union': 1152}
    assert all(float(x['memory_GiB']) <= 64 for x in points)
    assert all(F(x['success_lower']) >= F(1,2) for x in points)
    best = read(out/'best.json')
    assert len(best) == 6
    domains = {}
    for row in best:
        cfg = cfgs[row['candidate']]; p = row['parameters']
        key = (row['family'],cfg['p'],p['ell'],p['a'],p['q'],p['radius'])
        if key not in domains:
            ctor = FixedDomain if row['family'] == 'fixed_E' else model.old.Domain
            domains[key] = ctor(cfg['p'],(cfg['k']+p['ell'])//2,p['a'],p['q'],p['radius'])
        t = dict(family=row['family'],ell=p['ell'],a=p['a'],q=p['q'],radius=p['radius'],
                 left=row['left'],right=row['right'])
        checked = model.assess(cfg,t,domains[key],int(p['L']),p['u'],p['block_records'],
                               p['schedule'],p['tile_records'],2**36)
        for field in ('half_counts','representation_scope','domain_size','alpha','sigma_lower',
                      'sigma_method','bad_map_mass','success_lower','proposal_cap','resource_caps',
                      'prices','declared_bit_upper','declared_bit_upper_log2','heuristic_log2',
                      'memory','counters_payload','outside'):
            assert checked[field] == row[field], (row['candidate'],field)
        objective = row['final_selection_objective']
        assert min(float(x[objective]) for x in points if x['candidate']==row['candidate']) == row[objective]
        assert sum(row['memory']['areas'].values()) == row['memory']['reserved_bytes'] <= 2**36
        regions = sorted((x['offset'],x['offset']+x['bytes']) for x in row['memory']['regions'].values() if x['bytes'])
        assert all(a[1]<=b[0] for a,b in zip(regions,regions[1:]))
        assert row['counters_payload'] <= row['memory']['areas']['control_and_unrank_workspace']
        assert row['security_bits'] is None and not row['physical_RSS_certified']

    coverage = read(out/'coverage.json')
    status = Counter(x['status'] for x in coverage)
    assert status == summary['coverage_statuses']
    failures = Counter()
    for x in coverage:
        if x['status']=='DOMAIN_LIMIT_OR_WORKER_FAILURE':
            tail = x['reason'].splitlines()[-1]
            assert tail in ('domain.CapacityLimit: candidate histogram limit; no partial domain returned',
                            'domain.CapacityLimit: histogram DFS limit; no partial domain returned'), tail
            failures[tail] += 1
    assert sum(failures.values()) == 22
    assert sum(x['status']=='CURRENT_LAYOUT_RESERVATION_LIMIT' for x in coverage) == 6
    assert sum(x.get('points',0) for x in coverage) == len(points)
    decision = HERE/'FINAL_PARAMETER_DECISION.json'
    if decision.exists():
        table = read(decision)
        assert table['round_status']=='FINISHED'
        assert len(table['candidates']) == 3
        for x in table['candidates']:
            assert x['verdict']=='证据不足' and not x['supports_complete_classical128']
            b = next(r for r in best if r['candidate']==x['candidate'] and r['final_selection_objective']=='declared_bit_upper_log2')
            assert x['v7_best_upper_log2']==b['declared_bit_upper_log2']
    report = dict(status='PASS',saved_certificates_recomputed=6,points_verified=len(points),
                  source_bindings_verified=True,exact_C_verified=3,coverage_statuses=dict(status),
                  all_worker_failures_are_explicit_capacity_limits=True,worker_failure_tails=dict(failures),
                  declared_memory_bound_bytes=2**36,physical_RSS_certified=False,
                  target_attack_executed=False,check_is_model_consistency_not_global_security_proof=True)
    (out/'VERIFICATION.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--results',type=Path,default=HERE/'final_round')
    verify(parser.parse_args().results)
