"""Finite exact-C competitors; free L, affine targets, and explicit regeneration.

This is an explicitly named finite adaptation, not Li--Wang Algorithm 1 or a
claim to reproduce the original Restricted-Balls authors' finite sampler.
"""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import math
import sys

VENDOR=Path(__file__).resolve().parent/'vendor'
sys.path.insert(0,str(VENDOR))
from domain import Domain
from moments import inclusion,kernel,orbit,half_representations
from model import memory_layout,base_survival,offset_moments
from prior_prices import resource_prices as prior_prices


def ceil(x):
    x=F(x);return -(-x.numerator//x.denominator)


def lg(x):
    x=F(x);return None if x<=0 else math.log2(x.numerator)-math.log2(x.denominator)


def align(x):return (x+7)//8*8


@lru_cache(maxsize=128,typed=True)
def require_prime(p):
    if type(p)!=int or p<3 or p%2==0 or any(p%d==0 for d in range(3,math.isqrt(p)+1,2)):
        raise ValueError('odd prime field required')


def require_ints(**values):
    if any(type(x)!=int for x in values.values()):raise ValueError('integer arguments required')


def validate_composition(D,p,*,length=None,radius=None):
    if len(D)!=p or any(type(x)!=int or x<0 for x in D):raise ValueError('canonical nonnegative composition required')
    if length is not None and sum(D)!=length:raise ValueError('composition length mismatch')
    if radius is not None and any(x and min(i,p-i)>radius for i,x in enumerate(D)):
        raise ValueError('target symbol outside proved radius')


def validate_geometry(cfg,ell,u):
    p,n,k=(cfg[x] for x in ('p','n','k'));require_prime(p);require_ints(n=n,k=k,ell=ell,u=u)
    if not 0<k<n or not 1<=u<ell<=n-k or (k+ell)%2:raise ValueError('PGE/filter dimensions')


def multinomial(n,parts):
    if min(parts,default=0)<0 or sum(parts)!=n:return 0
    out=1
    for x in parts:out*=math.comb(n,x);n-=x
    return out


@lru_cache(maxsize=2048,typed=True)
def full_radius1_representations(p,h,a,q):
    """Exact |{x in B: e-x in B}| for e in {0,+1,-1}^h of weight 2a.

    The membership simplification requires h>=a+q and prime p>=5. The
    category polynomial is (1+2XY+(p-3)UV)^(h-2a) times
    (X+Y+XV+UY+(p-4)UV)^(2a). High counts are <=q; low weights
    lie in [a,a+q-high]. All allowed actual vectors, without mask tags.
    """
    require_prime(p);require_ints(h=h,a=a,q=q)
    if p<5 or not 0<a or 2*a>h or h<a+q or not 1<=q<=h:
        raise ValueError('full radius-1 counting domain')
    zeros=h-2*a;nonzeros=2*a;total=0
    for t in range(min(q,zeros)+1):
        for s in range(min(q-t,zeros-t)+1):
            zero_choices=multinomial(zeros,(t,s,zeros-t-s))*2**s*(p-3)**t
            for u in range(q-t+1):
                for v in range(q-t+1):
                    for z in range(q-t-max(u,v)+1):
                        hx=t+u+z;hy=t+v+z
                        remaining=nonzeros-u-v-z
                        if remaining<0:continue
                        for nx in range(max(0,a-s-v),min(remaining,a+q-hx-s-v)+1):
                            ny=remaining-nx;wx=s+nx+v;wy=s+ny+u
                            if not(a<=wx<=a+q-hx and a<=wy<=a+q-hy):continue
                            total+=zero_choices*multinomial(nonzeros,(nx,ny,v,u,z))*(p-4)**z
    return total


@lru_cache(maxsize=2048,typed=True)
def full_radius2_q1_representations(p,h,a,ones,twos):
    """Exact distinct-domain pair count, using a sparse category polynomial.

    For q=1 and h>=a+3 every low-only admissible word has a zero.
    Membership is: high=1 => weight=a; high=0 => weight=a, or
    weight=a+1 with a Lee-1 entry, or a+2 with a Lee-2 entry.
    Four presence bits retain these conditions for BOTH summands.
    """
    require_prime(p);require_ints(h=h,a=a,ones=ones,twos=twos)
    if p<11 or a<=0 or min(ones,twos)<0 or ones+twos>h or ones+2*twos!=2*a or h<a+3:
        raise ValueError('full radius-2/q1 counting domain')
    dp={(0,0,0,0,0):1};prefix=0;processed=0
    for t,count in ((0,h-ones-twos),(1,ones),(2,twos)):
        transitions={}
        for x in range(p):
            y=(t-x)%p;wx=min(x,p-x);wy=min(y,p-y)
            hx=int(wx>2);hy=int(wy>2);wx=0 if hx else wx;wy=0 if hy else wy
            bits=(1 if wx==1 else 2 if wx==2 else 0)|(4 if wy==1 else 8 if wy==2 else 0)
            key=(hx,hy,wx,wx+wy-t,bits);transitions[key]=transitions.get(key,0)+1
        for _ in range(count):
            prefix+=t;processed+=1;out={}
            for (hx,hy,wx,excess,bits),value in dp.items():
                for (ix,iy,iw,ie,ib),multiplicity in transitions.items():
                    nx=hx+ix;ny=hy+iy;nw=wx+iw;ne=excess+ie
                    if nx>1 or ny>1 or nw>a+2 or ne>4+2*(2-nx-ny):continue
                    other=prefix+ne-nw
                    if other>a+2 or min(nw,other)+2*(h-processed)<a:continue
                    key=(nx,ny,nw,ne,bits|ib)
                    out[key]=out.get(key,0)+value*multiplicity
            dp=out
    total=0
    def member(high,weight,has1,has2):
        return weight==a or (not high and ((weight==a+1 and has1) or (weight==a+2 and has2)))
    for (hx,hy,wx,excess,bits),value in dp.items():
        wy=2*a+excess-wx
        if member(hx,wx,bits&1,bits&2) and member(hy,wy,bits&4,bits&8):total+=value
    return total


def representations(D,a,q,p,radius):
    require_prime(p);require_ints(a=a,q=q,radius=radius)
    if radius not in (1,2) or radius>p//2:raise ValueError('unsupported proved radius')
    validate_composition(D,p,radius=radius);h=sum(D)
    if a<=0 or not 1<=q<=h or sum(min(i,p-i)*x for i,x in enumerate(D))!=2*a:
        raise ValueError('representation weight/free-coordinate constraints')
    if radius==1 and (p<5 or h<a+q):raise ValueError('outside proved full radius-1 domain')
    return _representations_cached(tuple(D),a,q,p,radius)


@lru_cache(maxsize=2048)
def _representations_cached(D,a,q,p,radius):
    geo=half_representations(D,a,q,p)
    if radius==1 and p>=5 and sum(D)>=a+q:
        full=full_radius1_representations(p,sum(D),a,q)
        if full<geo:raise ArithmeticError('full representation count smaller than subfamily')
        return full,'EXACT_FULL_RADIUS1',geo
    if radius==2 and q==1 and p>=11 and sum(D)>=a+3 and D[1]+D[-1]>0:
        full=full_radius2_q1_representations(p,sum(D),a,D[1]+D[-1],D[2]+D[-2])
        if full<geo:raise ArithmeticError('full representation count smaller than subfamily')
        return full,'EXACT_FULL_RADIUS2_Q1',geo
    return geo,'GEODESIC_SUBFAMILY_RADIUS2_NOT_FULL_SIGMA',geo


def survival(p,N,u,B,L,rL,rR,delta=F(1,32)):
    """Two valid alternatives. The bad-map Markov certificate is explicit.

    Conditional first moment is CONSTANT for every fixed map, thanks to tau.
    Markov on its conditional second moment gives bad-map mass <=delta.
    This avoids raising an annealed one-trial probability to a repeat count.
    """
    require_prime(p);require_ints(N=N,u=u,B=B,L=L,rL=rL,rR=rR)
    if not 1<=u<N or not 2<=L<=B or min(rL,rR)<2 or max(rL,rR)>B or rL%2 or rR%2:
        raise ValueError('fixed-point-free representation/subset constraints')
    delta=F(delta)
    if not 0<delta<1:raise ValueError('bad-map allowance')
    pointwise=base_survival(B,L,rL,rR)/p**u
    mean,second,annealed=offset_moments(p,N,u,B,L,rL,rR)
    gamma=1-F(1,p-1)
    alternatives=[]
    for name,sigma,bad in [('POINTWISE_ALL_MAPS',pointwise,F(0)),
                           ('CONDITIONAL_SECOND_MOMENT_MARKOV',delta*annealed,delta)]:
        success=(2*gamma*(1-bad)-1)*F(6,7)-F(1,10)
        if sigma and success>=F(1,2):alternatives.append((sigma,name,bad,success))
    if not alternatives:raise ValueError('no half-success cutoff certificate')
    sigma,name,bad,success=max(alternatives)
    return dict(sigma=sigma,method=name,bad_map_mass=bad,success=success,
        first_moment=mean,second_moment=second,annealed_pz_lower=annealed,
        pointwise_sigma=pointwise)


def layout(cfg,ell,u,B,L,b,domain,mode='resident',tile=2**20,budget=64*2**30):
    validate_geometry(cfg,ell,u);require_ints(B=B,L=L,b=b,tile=tile,budget=budget)
    if min(b,tile,budget)<1 or B!=domain.size or (domain.Q,domain.h)!=(cfg['p'],(cfg['k']+ell)//2):
        raise ValueError('layout domain/positive sizes')
    p,n,k=(cfg[x] for x in ('p','n','k'));b=min(b,L*L)
    if mode=='resident':
        if not 2<=L<=min(B,2**64-1):raise ValueError('resident subset bounds')
        mem=memory_layout(p,n,k,ell,u,B,L,b,domain);tile=L
    elif mode in ('tiled_partial','tiled_full'):
        if not 2<=L<=B or not 1<=tile<2**64:raise ValueError('tile bounds')
        if mode=='tiled_full' and L!=B:raise ValueError('full domain required')
        if mode=='tiled_partial' and L>=2**64:raise ValueError('rank-array index bound')
        tile=min(tile,L)
        mem=memory_layout(p,n,k,ell,u,B,2,min(b,4),domain)
        areas=mem['areas'];rankbytes=max(1,(B.bit_length()+7)//8)
        for key in ('base_left','base_right','base_sort_scratch','rank_sampling_tree','root_block','root_sort_scratch'):
            areas.pop(key)
        # Distinct arenas: L1 can remain suspended while L2 sorts/regenerates.
        # No implicit borrowing from root scratch or an inactive sampling tree.
        for branch in ('L1','L2'):
            areas[branch+'_tile']=tile*mem['base_record_bytes']
            areas[branch+'_sort_scratch']=tile*mem['base_record_bytes']
            areas[branch+'_right_record']=mem['base_record_bytes']
        areas['root_block']=b*mem['mid_record_bytes'];areas['root_sort_scratch']=b*mem['mid_record_bytes']
        # Partial subsets must be replayable: both rank arrays + one reusable AVL.
        # Full domain uses an implicit rank range, with no fake PRP/random oracle.
        areas['replayable_rank_arrays']=2*L*rankbytes if mode=='tiled_partial' else 0
        areas['rank_sampling_tree']=L*mem['rank_node_bytes'] if mode=='tiled_partial' else 0
        mem['reserved_bytes']=sum(areas.values())
    else:raise ValueError('unknown schedule')
    mem['areas']={key:align(value) for key,value in mem['areas'].items()}
    mem['reserved_bytes']=sum(mem['areas'].values());cursor=0;regions={}
    for key,value in mem['areas'].items():regions[key]=dict(offset=cursor,bytes=value);cursor+=value
    mem['regions']=regions
    if b>=2**64 or mem['reserved_bytes']>budget:raise ValueError('memory budget')
    mem.update(schedule=mode,tile_records=tile,root_block_records=b,budget_bytes=budget,
        fits_budget=True,allocator_implemented_at_target=False,process_RSS_certified=False,
        lifetime_contract='v6.1_disjoint_L1_L2_tiles_scratch_and_right_records',
        region_alignment_bytes=8,implicit_scratch_aliasing=False)
    return mem


def maximum_subset(cfg,ell,u,domain,b,mode,tile,budget):
    """Exact layout test, including b=min(b,L^2), never a fixed-b shortcut."""
    lo=2;hi=min(domain.size,2**64-1);best=0
    while lo<=hi:
        mid=(lo+hi)//2
        try:layout(cfg,ell,u,domain.size,mid,b,domain,mode,tile,budget)
        except ValueError:hi=mid-1
        else:best=mid;lo=mid+1
    return best


def allowances(cfg,ell,u,L,b,domain,mem):
    p,n,k=(cfg[x] for x in ('p','n','k'));N=k+ell;h=N//2;m=n-k;j=n-N
    W=L*L;M=W*W;P=p**u;b=min(b,W);B=domain.size
    X=F(W,P);Y=X;XY=(W+(M-W)*kernel(p,N,u))/P
    roots=(W+(M-W)*kernel(p,N,ell))/P
    blocks=min(F(ceil(F(W,b))),X/b+min(F(1),X))
    regenerated=min(ceil(F(W,b))*Y,XY/b+Y)
    passes=1+blocks;tile=mem['tile_records'];tiles=ceil(F(L,tile))
    resident=mem['schedule']=='resident';partial=mem['schedule']!='tiled_full'
    generated=F(2*L) if resident else L*(1+tiles)*passes
    small_sorts=F(L) if resident else L*passes
    queries=L*passes if resident else L*tiles*passes
    base_moves=small_sorts*(tile.bit_length()+2);root_moves=X*(b.bit_length()+2)
    mid=X+regenerated;rankbits=B.bit_length()+h.bit_length();height=2*(L+1).bit_length()
    f=max(1,((p-1).bit_length()+7)//8);width=n+1+m
    # Charge rectangular augmented PGE as in the existing runtime contract.
    ops={
        'field_multiplications':F(j*m*width)+generated*h*ell+roots*(j*N+m*n),
        'field_additions_subtractions':F(j*(m-1)*width+u)+generated*h*ell
            +2*mid*(ell-u)+queries*u+roots*(N+j*(N+1)+m*n),
        'field_inversions':F(j),'PGE_pivot_tests':F(j*m),
        'field_coordinate_comparisons':(base_moves+2*queries*tile.bit_length())*u
            +(X*b.bit_length()+2*regenerated*b.bit_length())*(ell-u)+roots*3*n,
        'histogram_counter_comparisons':roots*2*p,'histogram_counter_updates':roots*(n+p),
        'Lee_integer_additions':roots*n,
        'vector_coordinate_writes':generated*(h+ell)+mid*(N+2*(ell-u))+queries*u+roots*3*(n+1),
        'base_sort_record_moves':base_moves,'root_sort_record_moves':root_moves,
        'rank_bigint_multiplications':generated*h*(domain.radius+2),
        'rank_bigint_divisions':generated*(h*(domain.radius+2)+h+1),
        'rank_bigint_additions_subtractions':generated*(1+h*(domain.radius+2)),
        'rank_bigint_comparisons':generated*(len(domain.groups).bit_length()+h*(domain.radius+2)+2)
            +(8*L*height if partial else 0),
        'rank_tree_node_visits_upper':F(64*L*height if partial else 0),
        'expected_random_bits':F(2*(n-1)*n.bit_length()+2*u*p.bit_length()+(4*L*B.bit_length() if partial else 0)),
        'partition_uniform_integer_calls':F(n-1),'sample_uniform_integer_calls':F(2*L if partial else 0),
        'offset_uniform_field_calls':F(u),'first_join_queries':queries,'first_join_records':mid,
        'root_matches':roots,'streamed_leaf_vectors':generated,
        'join_pass_starts':passes,'base_tile_starts':tiles*passes,
        'allocation_zero_bytes':F(sum(v for key,v in mem['areas'].items() if key not in
            ('group_prefix_metadata','setup_removal_cache_reservation','setup_factorial_power_tables','setup_DP_buffers'))),
        'rank_array_byte_writes':F(2*L*max(1,(B.bit_length()+7)//8) if mem['schedule']=='tiled_partial' else 0),
        'rank_array_byte_reads':generated*max(1,(B.bit_length()+7)//8) if mem['schedule']=='tiled_partial' else F(0),
    }
    ops['sort_payload_bytes_written_upper']=base_moves*mem['base_record_bytes']+root_moves*mem['mid_record_bytes']
    ops['sort_payload_bytes_read_upper']=ops['sort_payload_bytes_written_upper']+2*f*(base_moves*u+X*b.bit_length()*(ell-u))
    setup={'field_multiplications':F(m*m*width),'field_additions_subtractions':F(m*(m-1)*width),
        'field_inversions':F(m),'PGE_pivot_tests':F(m*n),
        'public_setup_steps':F(sum(domain.stats.values())+len(domain.groups)*(domain.radius+10)),
        'allocation_zero_bytes':F(mem['reserved_bytes'])}
    return ops,setup,dict(EX=X,EY=Y,EXY_upper=XY,expected_blocks=blocks,
        regenerated_L2=regenerated,root_matches_upper=roots,generated_leaves=generated,
        join_queries=queries,join_passes=passes,rank_operand_bits=rankbits)


def prices(cfg,rank,address):
    """Same schoolbook declared-charge convention, with new units named."""
    b=(cfg['p']-1).bit_length();a=address;R=max(rank,a,b,cfg['n'].bit_length());loop=64*(a+R)
    base=prior_prices(cfg,rank,address)
    base['public_setup_steps']=16*R*R+4*a+loop
    for name in ('join_pass_starts','base_tile_starts'):base[name]=4*a+loop
    for name in ('allocation_zero_bytes','rank_array_byte_writes','rank_array_byte_reads',
                 ):base[name]=8+2*a+loop
    return base


def assess(cfg,t,domain,L,u,b,mode='resident',tile=2**20,budget=64*2**30):
    p,n,k=(cfg[x] for x in ('p','n','k'));ell=t['ell'];D=tuple(t['D']);a=t['a'];q=t['q'];radius=t['radius']
    validate_geometry(cfg,ell,u);validate_composition(cfg['C'],p,length=n)
    require_ints(a=a,q=q,radius=radius)
    validate_composition(D,p,length=(k+ell)//2,radius=radius)
    if (domain.Q,domain.h,domain.a,domain.q,domain.radius)!=(p,(k+ell)//2,a,q,radius):raise ValueError('domain mismatch')
    if any(2*x>y for x,y in zip(D,cfg['C'])) or sum(D)!=(k+ell)//2:raise ValueError('composition capacity')
    if sum(min(i,p-i)*x for i,x in enumerate(D))!=2*a:raise ValueError('target half weight')
    r,rep_scope,geo=representations(D,a,q,p,radius)
    if not r:raise ValueError('no certified representations')
    mem=layout(cfg,ell,u,domain.size,L,b,domain,mode,tile,budget)
    prob=survival(p,k+ell,u,domain.size,L,r,r)
    outside=[c-2*d for c,d in zip(cfg['C'],D)]
    alpha=F(orbit(D)**2*orbit(outside),orbit(cfg['C']))
    ops,setup,diag=allowances(cfg,ell,u,L,b,domain,mem)
    T=ceil(4/(alpha*prob['sigma']));K=len(ops)
    caps={key:ceil(10*K*T*x)+ceil(setup.get(key,0)) for key,x in ops.items()}
    for key,x in setup.items():
        if key not in caps:caps[key]=ceil(x)
    # Reserve two arbitrary-precision integers (used,limit) per class plus
    # proposal counters, indices and conservative bookkeeping within control.
    counter_bytes=sum(2*ceil(F(max(1,value.bit_length()),8))+64 for value in caps.values())
    counter_bytes+=8*ceil(F(max(1,T.bit_length()),8))+8192
    if counter_bytes>mem['areas']['control_and_unrank_workspace']:
        raise ValueError('cutoff counter workspace exceeds reservation')
    address=max(36,max(v.bit_length() for v in caps.values()));rank=max(diag['rank_operand_bits'],domain.setup_operand_max_bits)
    unit=prices(cfg,rank,address)
    if set(caps)-set(unit):raise ValueError('unknown resource price')
    total=sum(v*unit[key] for key,v in caps.items())
    # Explicitly heuristic: annealed Poisson root-survival and homogeneous keys.
    lambda_log=lg(prob['first_moment'])
    poisson_log=(lambda_log if lambda_log<-40 else
        0.0 if lambda_log>10 else math.log2(-math.expm1(-2**lambda_log)))
    proposal_charge=sum(F(v)*unit[key] for key,v in ops.items())
    setup_charge=sum(F(v)*unit[key] for key,v in setup.items())
    heuristic_log=lg(proposal_charge)+math.log2(math.log(2))-lg(alpha)-math.log2(1-1/(p-1))-poisson_log
    setup_log=lg(setup_charge);maximum=max(heuristic_log,setup_log)
    heuristic_log=maximum+math.log2(2**(heuristic_log-maximum)+2**(setup_log-maximum))
    return dict(candidate=cfg['name'],parameters=dict(ell=ell,u=u,a=a,q=q,radius=radius,L=str(L),
            block_records=min(b,L*L),tile_records=mem['tile_records'],schedule=mode),
        D=list(D),outside=outside,domain_size=str(domain.size),domain_groups=len(domain.groups),
        representation_scope=rep_scope,half_representation_count=str(r),geodesic_half_count=str(geo),
        alpha=str(alpha),sigma_lower=str(prob['sigma']),sigma_method=prob['method'],bad_map_mass_upper=str(prob['bad_map_mass']),
        ensemble_survival_moments={key:str(prob[key]) for key in ('first_moment','second_moment','annealed_pz_lower','pointwise_sigma')},
        proposal_cap=str(T),success_lower=str(prob['success']),memory=mem,
        resource_caps={key:str(v) for key,v in caps.items()},prices=unit,declared_bit_upper_log2=lg(total),
        declared_bit_upper=str(total),resource_projection_log2={key:lg(caps[key]) for key in
            ('field_multiplications','field_coordinate_comparisons','rank_bigint_divisions','streamed_leaf_vectors','first_join_queries','root_matches')},
        heuristic_half_success_log2=heuristic_log,
        heuristic_assumptions=['Poisson representation survival','homogeneous keys / annealed-to-repeated-key work',
            'complete-cycle mean cost as first-hit cost proxy; no exact stopping distribution'],
        diagnostics={key:str(v) for key,v in diag.items()},
        cost_version='v6.1_category_comparisons_disjoint_coroutine_buffers',
        comparison_per_unrank_upper=len(domain.groups).bit_length()+((k+ell)//2)*(radius+2)+2,
        cutoff_contract=dict(resources=K,failure_allowance='1/10',
            counters_payload_upper_bytes=counter_bytes,control_reserved_bytes=mem['areas']['control_and_unrank_workspace'],
            counted_kernel_precharge_monitor_implemented=True,complete_target_multi_proposal_driver_implemented=False),
        public_sum_ignored_by_this_valid_competitor=True,whole_bit_compiler_proved=False,
        target_lists_executed=False,general_FCSD_lower_bound=False,security_bits=None)
