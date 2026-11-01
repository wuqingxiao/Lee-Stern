"""Small prime-field reference for the source prose vector domain and joins.

Uniform independent half-subsets without replacement are OUR explicit rule.
The same two subsets are reused for the two first joins. No mask tags survive.
"""
from collections import Counter,defaultdict
from dataclasses import dataclass
import itertools
import math
import secrets
from reference_domain import unique_domain,lee,ReferenceLimit


def mv(matrix,word,p):return tuple(sum(x*y for x,y in zip(row,word))%p for row in matrix)


def rank(matrix,p):
    a=[list(row) for row in matrix];r=0
    for col in range(len(a[0])):
        pivot=next((i for i in range(r,len(a)) if a[i][col]),None)
        if pivot is None:continue
        a[r],a[pivot]=a[pivot],a[r];inv=pow(a[r][col],-1,p)
        a[r]=[x*inv%p for x in a[r]]
        for i in range(r+1,len(a)):
            factor=a[i][col];a[i]=[(x-factor*y)%p for x,y in zip(a[i],a[r])]
        r+=1
        if r==len(a):break
    return r


@dataclass(frozen=True)
class Instance:
    p:int
    H:tuple
    s:tuple
    C:tuple|None=None
    weight:int|None=None

    def __post_init__(self):
        p=self.p
        if type(p)!=int or p<2 or any(p%d==0 for d in range(2,math.isqrt(p)+1)):
            raise ValueError('reference PGE supports prime fields only')
        try:
            object.__setattr__(self,'H',tuple(tuple(row) for row in self.H))
            object.__setattr__(self,'s',tuple(self.s))
            if self.C is not None:object.__setattr__(self,'C',tuple(self.C))
        except TypeError as exc:raise ValueError('invalid instance arrays') from exc
        if not self.H or not self.H[0]:raise ValueError('empty matrix')
        n=len(self.H[0]);r=len(self.H)
        if not 0<r<n or len(self.s)!=r or any(len(row)!=n for row in self.H):raise ValueError('dimensions')
        if any(type(x)!=int or not 0<=x<p for row in self.H for x in row) or any(type(x)!=int or not 0<=x<p for x in self.s):
            raise ValueError('noncanonical field element')
        if rank(self.H,p)!=r:raise ValueError('H must have full row rank')
        if self.C is not None:
            if len(self.C)!=p or any(type(c)!=int or c<0 for c in self.C) or sum(self.C)!=n:
                raise ValueError('invalid C')
            weight=sum(c*lee(i,p) for i,c in enumerate(self.C))
            if self.weight is not None and self.weight!=weight:raise ValueError('weight inconsistent with C')
            object.__setattr__(self,'weight',weight)
        if type(self.weight)!=int or not 0<=self.weight<=n*(p//2):raise ValueError('weight required')

    @property
    def n(self):return len(self.H[0])

    def legal(self,word):
        return (len(word)==self.n and all(type(x)==int and 0<=x<self.p for x in word)
            and sum(lee(x,self.p) for x in word)==self.weight
            and (self.C is None or tuple(word.count(i) for i in range(self.p))==self.C)
            and mv(self.H,word,self.p)==self.s)


@dataclass(frozen=True)
class Prepared:
    instance:Instance
    I:tuple
    J:tuple
    rows:tuple
    B:tuple
    target:tuple
    base:tuple
    v:int
    u:int

    @property
    def h(self):return len(self.I)//2

    def complete(self,inside):
        p=self.instance.p;outside=tuple((row[-1]-sum(row[i]*x for i,x in zip(self.I,inside)))%p
            for row in self.rows[len(self.B):])
        word=[0]*self.instance.n
        for i,x in zip(self.I,inside):word[i]=x
        for i,x in zip(self.J,outside):word[i]=x
        return tuple(word)


def prepare(instance,I,a,q,radius,u,*,max_words=100000):
    I=tuple(I);n=instance.n;r=len(instance.H);N=len(I);ell=r-(n-N)
    if (N%2 or not 0<N<=n or len(set(I))!=N or any(type(i)!=int or not 0<=i<n for i in I)
        or type(u)!=int or not 0<=u<=ell):raise ValueError('invalid partition/filter')
    base=unique_domain(instance.p,N//2,a,q,radius,max_words)
    if not base:raise ValueError('empty base domain')
    if 4*a>instance.weight:raise ValueError('internal weight exceeds total weight')
    J=tuple(i for i in range(n) if i not in I)
    rows=[list(row)+[b] for row,b in zip(instance.H,instance.s)];p=instance.p
    for j,col in enumerate(J):
        pivot=next((i for i in range(j,r) if rows[i][col]),None)
        if pivot is None:return None
        rows[j],rows[pivot]=rows[pivot],rows[j]
        inv=pow(rows[j][col],-1,p);rows[j]=[x*inv%p for x in rows[j]]
        for i in range(r):
            if i!=j:
                factor=rows[i][col];rows[i]=[(x-factor*y)%p for x,y in zip(rows[i],rows[j])]
    rows=tuple(tuple(row) for row in rows[len(J):]+rows[:len(J)])
    # Source =_u means the LAST u syndrome coordinates. Bring exactly those
    # low rows to the front for internal prefix-key indexing, permuting the
    # entire augmented rows, including syndrome, together.
    low_order=tuple(range(ell-u,ell))+tuple(range(ell-u))
    rows=tuple(rows[i] for i in low_order)+rows[ell:]
    B=tuple(tuple(row[i] for i in I) for row in rows[:ell]);target=tuple(row[-1] for row in rows[:ell])
    return Prepared(instance,I,J,rows,B,target,base,4*a,u)


def sample_tables(prepared,size=None,rng=None):
    size=prepared.instance.p**prepared.u if size is None else size
    if type(size)!=int or not 1<=size<=len(prepared.base):raise ValueError('invalid distinct subset size')
    rng=rng if rng is not None else secrets.SystemRandom()
    # Explicit additional convention, not claimed to be specified by the paper.
    return tuple(tuple(prepared.base[i] for i in sorted(rng.sample(range(len(prepared.base)),size))) for _ in range(2))


def select(prepared,left,right):
    universe=frozenset(prepared.base);out=[]
    for given in (left,right):
        words=prepared.base if given is None else tuple(tuple(x) for x in given)
        if not words or len(set(words))!=len(words) or not set(words)<=universe:raise ValueError('invalid vector subset')
        out.append(words)
    return tuple(out)


def run_prepared(prepared,left=None,right=None,*,max_join_pairs=2000000,max_root_pairs=2000000):
    left,right=select(prepared,left,right);p=prepared.instance.p;h=prepared.h;u=prepared.u
    if 2*len(left)*len(right)>max_join_pairs:raise ReferenceLimit('first-join capacity limit')
    keyleft=[mv(tuple(row[:h] for row in prepared.B),x,p) for x in left]
    keyright=[mv(tuple(row[h:] for row in prepared.B),x,p) for x in right]
    buckets=defaultdict(list)
    for x,key in zip(left,keyleft):buckets[key[:u]].append((x,key))
    def first_join(target):
        result=[]
        for y,ky in zip(right,keyright):
            wanted=tuple((a-b)%p for a,b in zip(target,ky[:u]))
            for x,kx in buckets.get(wanted,()):
                result.append((x+y,tuple((a+b)%p for a,b in zip(kx[u:],ky[u:]))))
        return result
    L1=first_join((0,)*u);L2=first_join(prepared.target[:u])
    if len(L1)*len(L2)>max_root_pairs:raise ReferenceLimit('root-pair capacity limit')
    rootbuckets=defaultdict(list)
    for word,key in L1:rootbuckets[key].append(word)
    counters=Counter(L1=len(L1),L2=len(L2),selected_left=len(left),selected_right=len(right))
    outputs=Counter()
    for y,key in L2:
        wanted=tuple((a-b)%p for a,b in zip(prepared.target[u:],key))
        for x in rootbuckets.get(wanted,()):
            counters['raw_root_matches']+=1;inside=tuple((a+b)%p for a,b in zip(x,y))
            if sum(lee(z,p) for z in inside)!=prepared.v:continue
            counters['internal_weight_survivors']+=1;candidate=prepared.complete(inside)
            if prepared.instance.legal(candidate):outputs[candidate]+=1
    return {'status':'FOUND' if outputs else 'NO_SOLUTION','solutions':tuple(sorted(outputs)),
        'valid_root_multiplicities':dict(outputs),'counters':dict(counters),
        'scope':'literal prose unique-vector domain; fixed-C acceptance is a research extension'}


def direct_quadruple_outputs(prepared,left=None,right=None,max_quadruples=2000000):
    """Independent four-loop oracle, with no bucket-join construction."""
    left,right=select(prepared,left,right)
    if (len(left)*len(right))**2>max_quadruples:raise ReferenceLimit('quadruple oracle capacity')
    p=prepared.instance.p;outputs=Counter();raw=0
    for l,r,lc,rc in itertools.product(left,right,left,right):
        z=l+r;y=lc+rc
        if any(mv(prepared.B,z,p)[:prepared.u]):continue
        inside=tuple((x+t)%p for x,t in zip(z,y))
        if mv(prepared.B,inside,p)!=prepared.target:continue
        raw+=1
        if sum(lee(x,p) for x in inside)!=prepared.v:continue
        candidate=prepared.complete(inside)
        if prepared.instance.legal(candidate):outputs[candidate]+=1
    return {'solutions':tuple(sorted(outputs)),'valid_root_multiplicities':dict(outputs),'raw_root_matches':raw}


def legal_witnesses(instance,max_candidates=100000):
    n=instance.n;p=instance.p
    if instance.C is None:
        if p**n>max_candidates:raise ReferenceLimit('witness full-space limit')
        candidates=itertools.product(range(p),repeat=n)
    else:
        count=math.factorial(n)//math.prod(math.factorial(c) for c in instance.C)
        if count>max_candidates:raise ReferenceLimit('witness orbit limit')
        def generate(prefix,counts):
            if len(prefix)==n:yield prefix;return
            for x,c in enumerate(counts):
                if c:
                    rest=list(counts);rest[x]-=1
                    yield from generate(prefix+(x,),tuple(rest))
        candidates=generate((),instance.C)
    return tuple(word for word in candidates if instance.legal(word))
