"""Finite memory / free-L analysis. Never materializes target-size lists.

Random offsets are an explicit finite adaptation, not attributed to the
Restricted-Balls authors. Resource vectors are NOT scalar security bits.
"""
from fractions import Fraction as F
from functools import lru_cache
import math
from moments import orbit, half_representations, inclusion, kernel


def ceil(x):
    x=F(x)
    return -(-x.numerator//x.denominator)


def value(x):
    x=F(x)
    return {'exact':str(x),'log2':None if not x else math.log2(x.numerator)-math.log2(x.denominator)}


def align(n):return (n+7)//8*8


@lru_cache(maxsize=256)
def pge_probability(p,n,k,ell):
    N=k+ell
    return math.prod(1-F(1,p**i) for i in range(ell+1,n-k+1))/math.prod(
        1-F(1,p**i) for i in range(N+1,n+1))


def memory_layout(p,n,k,ell,u,B,L,b,domain):
    """Specified packed arrays; not Python heap or measured RSS.

    Sum all reserved areas, even when lifetimes do not overlap. A budget must
    cover setup as well. AVL rank keys and 4 uint64 fields are charged.
    """
    N=k+ell;h=N//2;m=n-k;s=max(1,(B.bit_length()+7)//8)
    field=max(1,((p-1).bit_length()+7)//8)
    if not (1<=u<=ell and 2<=L<=B and 1<=b<=L*L and max(L,b)<2**64):
        raise ValueError('memory/index constraints')
    base=align((h+ell)*field);mid=align((N+ell-u)*field)
    rank_node=align(s+32)
    q,r=domain.q,domain.radius
    # Explicit fixed-width group data plus prefix endpoints. Lookup for rank()
    # is unnecessary in the attack; unrank() uses sorted group-prefix search.
    groups=len(domain.groups)*align((r+2)*8+5*s+32)
    cache_entries=(q+1)**(r+1)*(q+1)
    cache=cache_entries*align((r+2)*8+(q*r+8)//8+32)
    setup_bits=max(math.factorial(h).bit_length(),(p**h).bit_length())
    tables=(h+1)*(r+3)*align((setup_bits+7)//8)
    dp=2*(q+1)*align((q*r+8)//8)
    areas={'base_left':L*base,'base_right':L*base,'base_sort_scratch':L*base,
        'rank_sampling_tree':L*rank_node,'root_block':b*mid,'root_sort_scratch':b*mid,
        'group_prefix_metadata':groups,'setup_removal_cache_reservation':cache,
        'setup_factorial_power_tables':tables,'setup_DP_buffers':dp,
        'matrices_and_completion_workspace':align((2*m*(n+1)+m*N+3*n)*field),
        'histogram':p*8,'permutation_and_inverse':2*n*8,
        'control_and_unrank_workspace':65536+64*(h+r+q+4)}
    return {'areas':areas,'reserved_bytes':sum(areas.values()),
        'base_record_bytes':base,'mid_record_bytes':mid,'rank_node_bytes':rank_node,
        'memory_model':'packed arrays, 64-bit indices, conservative sum of reserved areas; excludes Python/interpreter heap',
        'implemented_target_allocator':False}


def max_L(p,n,k,ell,u,B,b,domain,budget):
    b=max(1,b)
    # For the searched L>=2 and fixed b, a root block bigger than L^2 is
    # unnecessary. Return a bound and recheck each actual min(b,L^2) layout.
    floor=memory_layout(p,n,k,ell,u,B,2,min(b,4),domain)
    unit=3*floor['base_record_bytes']+floor['rank_node_bytes']
    fixed=floor['reserved_bytes']-2*unit-2*min(b,4)*floor['mid_record_bytes']
    return max(0,min(B,2**64-1,(budget-fixed-2*b*floor['mid_record_bytes'])//unit))


@lru_cache(maxsize=20000)
def base_survival(B,L,rL,rR):
    if min(rL,rR)<2:return F(0)
    R=rL*rR;f2=inclusion(B,L,2);f4=inclusion(B,L,4)
    mean=R*f2*f2
    second=4*R*f2*f2+2*R*(rL+rR-4)*f2*f4+R*(rL-2)*(rR-2)*f4*f4
    return mean*mean/second if second else F(0)


def offset_moments(p,N,u,B,L,rL,rR):
    """Ensemble moments for diagnostics, not used as pointwise repetition law."""
    R=rL*rR;f2=inclusion(B,L,2);f4=inclusion(B,L,4)
    mean=F(R,p**u)*f2*f2
    pairs=3*R*f2*f2+2*R*(rL+rR-4)*f2*f4+R*(rL-2)*(rR-2)*f4*f4
    second=mean+kernel(p,N,u)*pairs/p**u
    return mean,second,mean*mean/second if second else F(0)


def offset_list_counts(p,n,k,ell,u,L,b):
    N=k+ell;h=N//2;m=n-k;j=n-N;W=L*L;M=W*W;b=min(b,W)
    P=p**u;q0=kernel(p,N,u);ql=kernel(p,N,ell)
    X=F(W,P);Y=X;XY=(W+(M-W)*q0)/P
    roots=(W+(M-W)*ql)/P
    chunks=ceil(F(W,b));K=min(F(chunks),X/b+min(F(1),X))
    KY=min(chunks*Y,XY/b+Y)
    queries=L*(1+K);lh=L.bit_length();bh=b.bit_length()
    ops={'field_multiplications':F(2*L*h*ell+j*m*(n+1))+roots*(j*N+m*n),
        'field_additions_subtractions':F(2*L*h*ell+j*(m-1)*(n+1)+u)
            +(X+KY)*(ell-u)+queries*u+KY*(ell-u)+roots*(N+j*(N+1)+m*n),
        'field_inversions':F(j),
        'field_coordinate_comparisons':(L*lh+2*queries*lh)*u+(X*bh+2*KY*bh)*(ell-u)
            +roots*3*N,
        'histogram_counter_comparisons':roots*2*p,
        'Lee_integer_additions':roots*N,
        'histogram_counter_updates':roots*n,
        'vector_coordinate_writes':F(2*L*h)+(X+KY)*N+roots*3*n,
        'base_sort_record_moves':F(L*(lh+2)),
        'root_sort_record_moves':X*(bh+2),
        'partition_uniform_integer_calls':F(n-1),
        'offset_uniform_field_calls':F(u),
        'PGE_pivot_tests':F(m*j),
        'first_join_queries':queries,'root_join_queries':KY,
        'first_join_records':X+KY,'root_matches':roots}
    return ops,{'EX':X,'EY':Y,'EXY_upper':XY,'expected_blocks_upper':K,
        'expected_L2_regenerated_upper':KY,'root_matches_upper':roots}


def canonical_counts(p,n,k,ell,u,L,b,domain):
    ops,diag=offset_list_counts(p,n,k,ell,u,L,b)
    B=domain.size;g=domain.generation_allowances(2*L,True);h=(k+ell)//2
    height=2*(L+1).bit_length();rank_bits=B.bit_length()
    ops.update({'rank_bigint_multiplications':F(g['rank_bigint_multiplications']),
        'rank_bigint_divisions':F(g['rank_bigint_divisions']+g['rank_bigint_divmod']),
        'rank_bigint_additions_subtractions':F(g['rank_bigint_subtractions']),
        'rank_bigint_comparisons':F(g['rank_prefix_comparisons']+8*L*height),
        'rank_tree_node_visits_upper':F(64*L*height),
        'expected_random_bits':F(2*(n-1)*n.bit_length()+4*L*rank_bits+2*u*p.bit_length()),
        'sample_uniform_integer_calls':F(2*L)})
    # No scalar price is assigned to these different resource units.
    s=max(1,((p-1).bit_length()+7)//8)
    base=align((h+ell)*s);mid=align((k+ell+ell-u)*s)
    ops['sort_payload_bytes_written_upper']=ops['base_sort_record_moves']*base+ops['root_sort_record_moves']*mid
    ops['sort_payload_bytes_read_upper']=ops['sort_payload_bytes_written_upper']
    return ops,diag


@lru_cache(maxsize=256)
def cached_representations(D,a,q,p):return half_representations(D,a,q,p)


@lru_cache(maxsize=256)
def cached_orbit(D):return orbit(D)


def assess(cfg,template,domain,L,u,b,budget,*,detailed=True):
    p,n,k,C=cfg['p'],cfg['n'],cfg['k'],cfg['C'];ell=template['ell'];a=template['a'];q=template['q']
    N=k+ell;h=N//2
    if p<3 or any(p%d==0 for d in range(2,math.isqrt(p)+1)):raise ValueError('odd prime required')
    if not 0<k<n or N%2 or not 1<=u<ell<=n-k:raise ValueError('dimensions/filter')
    if len(C)!=p or sum(C)!=n or any(type(c)!=int or c<0 for c in C):raise ValueError('exact C')
    if (domain.Q,domain.h,domain.a,domain.q,domain.radius)!=(p,h,a,q,template['radius']):
        raise ValueError('wrong domain')
    left=tuple(template['left'].get(str(i),0) for i in range(p))
    right=tuple(template['right'].get(str(i),0) for i in range(p))
    if min(left+right)<0 or any(c-x-y<0 for c,x,y in zip(C,left,right)):raise ValueError('composition capacity')
    if any(sum(D)!=h or sum(c*min(i,p-i) for i,c in enumerate(D))!=2*a for D in (left,right)):
        raise ValueError('composition length/weight')
    if any((x+y) and min(i,p-i)>domain.radius for i,(x,y) in enumerate(zip(left,right))):
        raise ValueError('target radius')
    if not 1<=q<=h or p%2!=1 or a<=0:raise ValueError('subfamily constraints')
    rL=cached_representations(left,a,q,p);rR=cached_representations(right,a,q,p)
    beta=base_survival(domain.size,L,rL,rR);sigma=beta/p**u
    if not sigma:raise ValueError('empty subfamily')
    outside=tuple(c-x-y for c,x,y in zip(C,left,right))
    alpha=F(cached_orbit(left)*cached_orbit(right)*cached_orbit(outside),cached_orbit(tuple(C)))
    # Keep the proven conservative rational bound in certificates rather than
    # serializing tens of thousands of digits of the exact finite product.
    gamma=1-F(1,p-1)
    success=(2*gamma-1)*F(6,7)-F(1,10)
    layout=memory_layout(p,n,k,ell,u,domain.size,L,min(b,L*L),domain)
    if layout['reserved_bytes']>budget:raise ValueError('memory budget')
    ops,diag=canonical_counts(p,n,k,ell,u,L,b,domain)
    T=ceil(4/(alpha*sigma));K=len(ops)
    caps={name:ceil(10*K*T*x) for name,x in ops.items()}
    if not detailed:
        return {'field_multiplications':caps['field_multiplications'],
            'field_coordinate_comparisons':caps['field_coordinate_comparisons'],
            'rank_bigint_arithmetic':caps['rank_bigint_multiplications']+caps['rank_bigint_divisions'],
            'reserved_bytes':layout['reserved_bytes']}
    return {'family':'restricted_distinct_random_offset_adaptation','template_id':template['id'],
        'p':p,'n':n,'k':k,'ell':ell,'u':u,'L':str(L),'block_records':min(b,L*L),
        'a':a,'q':q,'radius':domain.radius,'domain_size':str(domain.size),
        'half_representations':[str(rL),str(rR)],'alpha':value(alpha),
        'leaf_subfamily_existence_lower':value(beta),'pointwise_survival_given_split_PGE_lower':value(sigma),
        'PGE_ensemble_probability_lower':value(gamma),'ensemble_per_proposal_success_lower':value(gamma*alpha*sigma),
        'proposal_cap':str(T),'ensemble_success_with_counted_cutoffs_lower':value(success),
        'good_key_fraction_lower':value(max(F(0),2*gamma-1)),
        'memory_budget_bytes':budget,'memory':layout,
        'operation_allowances':{name:value(x) for name,x in ops.items()},
        'simultaneous_counted_caps':{name:str(x) for name,x in caps.items()},
        'list_diagnostics':{name:value(x) for name,x in diag.items()},
        'rank_operand_max_bits':domain.size.bit_length()+h.bit_length(),
        'setup_operand_max_bits':domain.setup_operand_max_bits,
        'setup_one_time_counts':dict(domain.stats),
        'cutoff_monitor_implemented':False,'target_lists_generated':False,
        'sampling':'two independent uniform subsets, reused for both first joins; fresh uniform offset each proposal',
        'scope':'ordinary uniform full-row-rank H, exact C, ensemble success; packed-array realization specification',
        'original_author_algorithm_reproduction':False,'whole_bit_cost':None,'security_bits':None}
