"""One bounded extension: larger r/q and asymmetric halves; same v6.1 prices."""
from pathlib import Path
from types import SimpleNamespace
import math
import sys
from fractions import Fraction as F

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'vendor/v61'))
import competitive as old
from fixed_domain import FixedDomain


def layout(cfg,ell,u,B,L,b,domain,mode='resident',tile=2**20,budget=64*2**30):
    if not isinstance(domain,FixedDomain):return old.layout(cfg,ell,u,B,L,b,domain,mode,tile,budget)
    # Reuse ONLY repaired record/rank/matrix reservations. Histogram-specific
    # metadata is replaced by exact restricted-sphere DP metadata below.
    proxy=SimpleNamespace(Q=domain.Q,h=domain.h,size=domain.size,q=0,radius=0,groups=(None,))
    # Defer budget gate until the actual metadata has been inserted.
    mem=old.layout(cfg,ell,u,B,L,b,proxy,mode,tile,max(budget,2**1000))
    s=max(1,(domain.setup_operand_max_bits+7)//8);a=mem['areas']
    a['group_prefix_metadata']=old.align(domain.table_bytes+domain.Q*16+4*s)
    a['setup_removal_cache_reservation']=0;a['setup_factorial_power_tables']=0
    a['setup_DP_buffers']=old.align(2*(domain.a+1)*s)
    a['control_and_unrank_workspace']=old.align(65536+64*(domain.h+domain.radius+domain.q+4))
    mem['reserved_bytes']=sum(a.values());mem['budget_bytes']=budget
    if mem['reserved_bytes']>budget:raise ValueError('memory budget')
    cursor=0;mem['regions']={}
    for name,value in a.items():mem['regions'][name]=dict(offset=cursor,bytes=value);cursor+=value
    mem['domain_metadata_model']='fixed-E exact sphere DP; no histogram closure table'
    return mem


def maximum_subset(cfg,ell,u,domain,b,mode,tile,budget):
    lo=2;hi=min(domain.size,2**64-1);best=0
    while lo<=hi:
        mid=(lo+hi)//2
        try:layout(cfg,ell,u,domain.size,mid,b,domain,mode,tile,budget)
        except ValueError:hi=mid-1
        else:best=mid;lo=mid+1
    return best


def allowances(cfg,ell,u,L,b,domain,mem):
    ops,setup,diag=old.allowances(cfg,ell,u,L,b,domain,mem)
    if isinstance(domain,FixedDomain):
        generated=diag['generated_leaves'];partial=mem['schedule']!='tiled_full';height=2*(L+1).bit_length()
        # Direct DP-symbol unranking, NOT histogram/multiset unranking.
        ops['rank_bigint_multiplications']=F(0)
        ops['rank_bigint_divisions']=generated*(domain.q+1)
        ops['rank_bigint_additions_subtractions']=generated*domain.low*(2*domain.radius+1)
        ops['rank_bigint_comparisons']=generated*(2+domain.low*(2*domain.radius+1))+(8*L*height if partial else 0)
    return ops,setup,diag


def half_count(D,domain):
    if isinstance(domain,FixedDomain):return domain.half_count(D),'EXACT_FIXED_E_ALL_REPRESENTATIONS'
    p=domain.Q;old.validate_composition(D,p,length=domain.h,radius=domain.radius)
    if sum(min(i,p-i)*v for i,v in enumerate(D))!=2*domain.a:raise ValueError('half weight')
    if domain.radius<=2:return old.representations(tuple(D),domain.a,domain.q,p,domain.radius)[:2]
    value=old.half_representations(D,domain.a,domain.q,p)
    return value,'GEODESIC_SUBFAMILY_IN_FULL_UNION_NOT_ALL_SIGMA'


def assess(cfg,t,domain,L,u,b,mode='resident',tile=2**20,budget=64*2**30):
    p,n,k=(cfg[x] for x in ('p','n','k'));ell=t['ell'];A=tuple(t['left']);B=tuple(t['right'])
    old.validate_geometry(cfg,ell,u);old.validate_composition(cfg['C'],p,length=n)
    if (domain.Q,domain.h,domain.a,domain.q,domain.radius)!=(p,(k+ell)//2,t['a'],t['q'],t['radius']):raise ValueError('domain geometry')
    expected=domain.low if isinstance(domain,FixedDomain) else domain.h
    old.validate_composition(A,p,length=expected,radius=domain.radius);old.validate_composition(B,p,length=expected,radius=domain.radius)
    outside=[v-x-y for v,x,y in zip(cfg['C'],A,B)]
    if isinstance(domain,FixedDomain):outside[0]-=2*domain.q
    if min(outside)<0 or sum(outside)!=n-k-ell:raise ValueError('secret composition capacities')
    rL,scopeL=half_count(A,domain);rR,scopeR=half_count(B,domain)
    if not min(rL,rR):raise ValueError('no certified representation subevent')
    mem=layout(cfg,ell,u,domain.size,L,b,domain,mode,tile,budget)
    prob=old.survival(p,k+ell,u,domain.size,L,rL,rR)
    alpha=F(old.orbit(A)*old.orbit(B)*old.orbit(outside),old.orbit(cfg['C']))
    ops,setup,diag=allowances(cfg,ell,u,L,b,domain,mem);T=old.ceil(4/(alpha*prob['sigma']));K=len(ops)
    caps={name:old.ceil(10*K*T*v)+old.ceil(setup.get(name,0)) for name,v in ops.items()}
    for name,v in setup.items():
        if name not in caps:caps[name]=old.ceil(v)
    address=max(36,max(v.bit_length() for v in caps.values()));rank=max(diag['rank_operand_bits'],domain.setup_operand_max_bits)
    prices=old.prices(cfg,rank,address)
    if set(caps)-set(prices):raise ValueError('unpriced unit')
    total=sum(v*prices[name] for name,v in caps.items())
    counter_bytes=sum(2*old.ceil(F(max(1,v.bit_length()),8))+64 for v in caps.values())+8*old.ceil(F(T.bit_length(),8))+8192
    if counter_bytes>mem['areas']['control_and_unrank_workspace']:raise ValueError('counter reservation')
    exponent=old.lg(prob['first_moment'])
    pz=exponent if exponent<-40 else 0.0 if exponent>10 else math.log2(-math.expm1(-2**exponent))
    cycle=sum(F(v)*prices[name] for name,v in ops.items());initial=sum(F(v)*prices[name] for name,v in setup.items())
    work=old.lg(cycle)+math.log2(math.log(2))-old.lg(alpha)-math.log2(1-1/(p-1))-pz
    top=max(work,old.lg(initial));heuristic=top+math.log2(2**(work-top)+2**(old.lg(initial)-top))
    return dict(candidate=cfg['name'],family=t['family'],left=list(A),right=list(B),outside=outside,
        asymmetric=A!=B,parameters=dict(ell=ell,u=u,a=t['a'],q=t['q'],radius=t['radius'],
            L=str(L),block_records=mem['root_block_records'],tile_records=mem['tile_records'],schedule=mode),
        half_counts=[str(rL),str(rR)],representation_scope=[scopeL,scopeR],domain_size=str(domain.size),
        alpha=str(alpha),sigma_lower=str(prob['sigma']),sigma_method=prob['method'],bad_map_mass=str(prob['bad_map_mass']),
        success_lower=str(prob['success']),proposal_cap=str(T),resource_caps={k:str(v) for k,v in caps.items()},
        prices=prices,declared_bit_upper=str(total),declared_bit_upper_log2=old.lg(total),heuristic_log2=heuristic,
        heuristic_assumptions=['Poisson representation survival','homogeneous keys','full-cycle cost proxy for first hit'],
        memory=mem,counters_payload=counter_bytes,scope='ordinary uniform full-rank H; exact C; specified finite scheduling',
        base_is_full_permutation_closure=not isinstance(domain,FixedDomain),
        L_equals_p_power_u=L==p**u,three_u_ge_ell=3*u>=ell,
        all_target_sigma_exact=False,physical_RSS_certified=False,security_bits=None)
