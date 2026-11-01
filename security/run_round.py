"""Exactly one declared round; unsupported/timeout is NOT attack exclusion."""
import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time

import final_model as model
from fixed_domain import FixedDomain

HERE=Path(__file__).resolve().parent


def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(value,ensure_ascii=False,indent=2),encoding='utf-8');tmp.replace(path)
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def composition_seed(cfg,ell,radius,q,family):
    p=cfg['p'];h=(cfg['k']+ell)//2;length=h-q if family=='fixed_E' else h
    caps=[(x-(2*q if j==0 and family=='fixed_E' else 0))//2 if min(j,p-j)<=radius else 0 for j,x in enumerate(cfg['C'])]
    if min(caps)<0 or sum(caps)<length:return None
    total=sum(caps);D=[length*x//total for x in caps]
    order=sorted(range(p),key=lambda j:(-(length*caps[j]%total),j))
    for j in order:
        if sum(D)==length:break
        if D[j]<caps[j]:D[j]+=1
    weight=sum(min(j,p-j)*x for j,x in enumerate(D))
    if weight%2:
        choices=[]
        for i,x in enumerate(D):
            if not x:continue
            for j,cap in enumerate(caps):
                if D[j]>=cap or abs(min(i,p-i)-min(j,p-j))!=1:continue
                score=math.log(x)-math.log(D[j]+1);choices.append((score,i,j))
        if not choices:return None
        _,i,j=max(choices);D[i]-=1;D[j]+=1
    return D


def variants(cfg,ell,radius,q,family,shifts):
    D=composition_seed(cfg,ell,radius,q,family)
    if D is None:return []
    a=sum(min(j,cfg['p']-j)*x for j,x in enumerate(D))//2;out=[]
    for shift in shifts:
        A=list(D);B=list(D)
        for j,delta in ((0,shift),(1,-2*shift),(2,shift)):A[j]+=delta;B[j]-=delta
        if min(A+B)<0:continue
        if family=='full_union' and min(A[0],B[0])<q:continue
        if not any(min(j,cfg['p']-j)%2 and x for j,x in enumerate(A)):continue
        if not any(min(j,cfg['p']-j)%2 and x for j,x in enumerate(B)):continue
        out.append(dict(family=family,ell=ell,radius=radius,q=q,a=a,left=A,right=B,shift=shift))
    return out


def brief(row):
    return dict(candidate=row['candidate'],family=row['family'],asymmetric=row['asymmetric'],
        **row['parameters'],domain_log2=model.old.lg(row['domain_size']),
        alpha_log2=model.old.lg(row['alpha']),sigma_log2=model.old.lg(row['sigma_lower']),
        declared_bit_upper_log2=row['declared_bit_upper_log2'],heuristic_log2=row['heuristic_log2'],
        memory_GiB=row['memory']['reserved_bytes']/2**30,success_lower=row['success_lower'],
        representation_left=row['representation_scope'][0],representation_right=row['representation_scope'][1])


def evaluate_domain(cfg,templates,domain,config):
    rows=[];winners={};rejected=Counter()
    if not templates:return rows,[],rejected
    ell=templates[0]['ell']
    for u in sorted({1,max(1,ell//3),(ell+2)//3,ell-1}):
        for b in config['blocks']:
            for mode in config['schedules']:
                maximum=domain.size if mode=='tiled_full' else model.maximum_subset(cfg,ell,u,domain,b,mode,config['tile_records'],config['memory_bytes'])
                choices=[domain.size] if mode=='tiled_full' else sorted({maximum}|{2**e for e in config['list_bits'] if 2**e<=maximum}|
                    ({cfg['p']**u} if config.get('include_p_power_lists') and cfg['p']**u<=maximum else set()))
                for L in choices:
                    if L<2:continue
                    for t in templates:
                        try:row=model.assess(cfg,t,domain,L,u,b,mode,config['tile_records'],config['memory_bytes'])
                        except ValueError as exc:rejected[str(exc)]+=1;continue
                        rows.append(brief(row))
                        for obj in ('declared_bit_upper_log2','heuristic_log2'):
                            key=(mode,t['shift']!=0,obj)
                            if key not in winners or row[obj]<winners[key][obj]:winners[key]=row
    selected=[dict(row,selection_objective=key[2],selection_asymmetric=key[1]) for key,row in winners.items()]
    return rows,selected,rejected


def union_worker(job,out):
    cfg=job['candidate'];case=job['case'];config=job['config'];ts=variants(cfg,*case,'full_union',config['full_union']['asymmetry_shifts'])
    if not ts:write(out,dict(status='NO_TEMPLATE_IN_THIS_DECLARED_SEED_RULE',rows=[],selected=[]));return
    ell,radius,q=case
    d=model.old.Domain(cfg['p'],(cfg['k']+ell)//2,ts[0]['a'],q,radius,max_candidates=config['full_union']['candidate_histogram_limit'])
    rows,selected,rejected=evaluate_domain(cfg,ts,d,config)
    write(out,dict(status='DOMAIN_GRID_COMPLETE',rows=rows,selected=selected,rejections=dict(rejected),domain_size=str(d.size),domain_groups=len(d.groups)))


def run(out):
    if out.exists():raise FileExistsError('use a fresh output directory')
    out.mkdir(parents=True);config=read(HERE/'config.json');cfgs=read(HERE/'inputs/candidates.json')
    write(out/'config.json',config);write(out/'candidates.json',cfgs)
    code={p.relative_to(HERE).as_posix():digest(p) for p in list(HERE.glob('*.py'))+list((HERE/'vendor').rglob('*.py'))}
    rows=[];selected=[];coverage=[];start=time.perf_counter();counts=Counter()
    for cfg in cfgs:
        f=config['fixed_E'];index=0
        for ell in f['ells']:
            for radius in f['radii']:
                for q in f['q']:
                    ts=variants(cfg,ell,radius,q,'fixed_E',f['asymmetry_shifts']);index+=1
                    if not ts:
                        coverage.append(dict(candidate=cfg['name'],family='fixed_E',case=[ell,radius,q],status='NO_TEMPLATE_IN_THIS_DECLARED_SEED_RULE'));continue
                    try:
                        d=FixedDomain(cfg['p'],(cfg['k']+ell)//2,ts[0]['a'],q,radius,max_cells=f['DP_cell_limit'])
                        br,sel,reject=evaluate_domain(cfg,ts,d,config)
                        rows.extend(br);selected.extend(sel);counts.update(reject);counts['fixed_E_points']+=len(br)
                        coverage.append(dict(candidate=cfg['name'],family='fixed_E',case=[ell,radius,q],status='DOMAIN_GRID_COMPLETE',points=len(br),templates=len(ts)))
                    except ValueError as exc:coverage.append(dict(candidate=cfg['name'],family='fixed_E',case=[ell,radius,q],status='ANALYSIS_LIMIT_OR_UNSUPPORTED',reason=str(exc)))
                    if index%10==0:print(json.dumps(dict(candidate=cfg['name'],stage='fixed_E',domains_done=index,domains_total=100,points=len(rows),elapsed=round(time.perf_counter()-start,1))),flush=True)
        write(out/(cfg['name']+'_fixed_selected.json'),[r for r in selected if r['candidate']==cfg['name']])
    for cfg in cfgs:
        u=config['full_union']
        for ell in u['ells']:
            for radius,q in u['radius_q_pairs']:
                case=[ell,radius,q];entry=dict(candidate=cfg['name'],family='full_union',case=case)
                cache=(q+1)**(radius+2)*model.old.align((radius+2)*8+(q*radius+8)//8+32)
                if cache>config['memory_bytes']:
                    coverage.append(dict(entry,status='CURRENT_LAYOUT_RESERVATION_LIMIT',reason='declared cache reservation alone exceeds64GiB; not general attack exclusion'));continue
                name=cfg['name']+f'_e{ell}_r{radius}_q{q}';job=out/'union_jobs'/(name+'.json');result=job.with_suffix('.result.json')
                write(job,dict(candidate=cfg,case=case,config=config))
                try:
                    r=subprocess.run([sys.executable,str(HERE/'run_round.py'),'--union-job',str(job),'--union-result',str(result)],capture_output=True,text=True,timeout=u['worker_wall_seconds'])
                    if r.returncode or not result.exists():
                        coverage.append(dict(entry,status='DOMAIN_LIMIT_OR_WORKER_FAILURE',reason=r.stderr[-1800:]));continue
                    value=read(result);rows.extend(value['rows']);selected.extend(value['selected']);counts['full_union_points']+=len(value['rows']);counts.update(value.get('rejections',{}));coverage.append(dict(entry,status=value['status'],points=len(value['rows'])))
                except subprocess.TimeoutExpired:coverage.append(dict(entry,status='ANALYSIS_TIMEOUT',reason='wall limit; no security conclusion'))
                print(json.dumps(dict(candidate=cfg['name'],stage='full_union',case=case,points=len(rows),elapsed=round(time.perf_counter()-start,1))),flush=True)
    write(out/'coverage.json',coverage);write(out/'selected.json',selected)
    with (out/'points.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    best=[]
    for cfg in cfgs:
        candidates=[r for r in selected if r['candidate']==cfg['name']]
        for objective in ('declared_bit_upper_log2','heuristic_log2'):
            best.append(dict(min(candidates,key=lambda r:r[objective]),final_selection_objective=objective))
    write(out/'best.json',best)
    summary=dict(status='DECLARED_FINAL_ROUND_FINISHED',candidates=3,points=len(rows),counts=dict(counts),
        coverage_statuses=dict(Counter(r['status'] for r in coverage)),elapsed_seconds=time.perf_counter()-start,
        code_binding=code,input_sha256=digest(HERE/'inputs/candidates.json'),config_sha256=digest(HERE/'config.json'),
        results_are_upper_budgets_and_labeled_heuristics=True,target_attack_executed=False,security_bits=None)
    if any(digest(HERE/name)!=value for name,value in code.items()):raise RuntimeError('code changed during round')
    write(out/'SUMMARY.json',summary);print(json.dumps(summary,ensure_ascii=True,indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path);p.add_argument('--union-job',type=Path);p.add_argument('--union-result',type=Path)
    a=p.parse_args()
    if a.union_job:union_worker(read(a.union_job),a.union_result)
    elif a.out:run(a.out)
    else:p.error('choose --out or --union-job')
