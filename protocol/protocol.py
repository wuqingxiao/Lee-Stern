"""Exact-composition Lee-Stern research protocol; SHA-256 in classical ROM.
Not a signature implementation, not a production security certification.
"""
import hashlib,hmac,math,secrets
from dataclasses import dataclass,field
from fast_sampling import FastRandom
VERSION=b'LeeStern-ExactC-v1'

def width(bound):return max(1,((bound-1).bit_length()+7)//8)
def vector(v,length,bound):
    return isinstance(v,(list,tuple)) and len(v)==length and all(type(x)==int and 0<=x<bound for x in v)
def pack(v,bound):
    if not all(type(x)==int and 0<=x<bound for x in v):raise ValueError('noncanonical integer')
    size=width(bound)
    return bytes(v) if size==1 else b''.join(x.to_bytes(size,'big') for x in v)
def frame(*parts):
    return b''.join(len(x).to_bytes(8,'big')+x for x in parts)
def natural(x):
    if type(x)!=int or x<0:raise ValueError('invalid natural')
    return x.to_bytes(max(1,(x.bit_length()+7)//8),'big')
def counts(v,p):
    out=[0]*p
    for x in v:out[x]+=1
    return tuple(out)
def permute(v,sigma):return tuple(v[i] for i in sigma)
def add(a,b,p):return tuple((x+y)%p for x,y in zip(a,b))
def matrix_vector(H,v,p):return tuple(sum(a*b for a,b in zip(row,v))%p for row in H)
def rank(H,p):
    A=[list(row) for row in H];r=0
    for col in range(len(A[0])):
        pivot=next((i for i in range(r,len(A)) if A[i][col]),None)
        if pivot is None:continue
        A[r],A[pivot]=A[pivot],A[r];inv=pow(A[r][col],-1,p)
        A[r]=[(x*inv)%p for x in A[r]]
        for i in range(r+1,len(A)):
            f=A[i][col]
            if f:A[i]=[(x-f*y)%p for x,y in zip(A[i],A[r])]
        r+=1
        if r==len(A):break
    return r

@dataclass(frozen=True)
class Context:
    p:int
    n:int
    k:int
    C:tuple
    H:tuple
    s:tuple
    backend:str='python'
    _prefix:object=field(init=False,repr=False,compare=False)
    _matrix:object=field(init=False,repr=False,compare=False,default=None)
    def __post_init__(self):
        p,n,k=self.p,self.n,self.k
        if any(type(x)!=int for x in (p,n,k)) or p<3 or p%2==0 or any(p%d==0 for d in range(3,math.isqrt(p)+1,2)) or not 0<k<n:
            raise ValueError('invalid dimensions or odd prime')
        if not vector(self.C,p,n+1) or sum(self.C)!=n:raise ValueError('invalid exact composition')
        if len(self.H)!=n-k or any(not vector(row,n,p) for row in self.H) or not vector(self.s,n-k,p):raise ValueError('invalid H or s')
        object.__setattr__(self,'C',tuple(self.C));object.__setattr__(self,'H',tuple(tuple(row) for row in self.H));object.__setattr__(self,'s',tuple(self.s))
        if rank(self.H,p)!=n-k:raise ValueError('H is not full rank')
        ctx=frame(VERSION,natural(p),natural(n),natural(k),pack(self.C,n+1),pack(tuple(x for row in self.H for x in row),p),pack(self.s,p))
        # Cached internal hash state after the FULL context; no context digest substitution.
        object.__setattr__(self,'_prefix',hashlib.sha256(frame(b'COM-SHA256-SALT32-v1',ctx)))
        if self.backend=='numpy':
            import numpy as np
            if n*(p-1)**2>2**63-1:raise ValueError('int64 accumulation overflow')
            A=np.array(self.H,dtype=np.int64);A.flags.writeable=False
            object.__setattr__(self,'_matrix',A)
        elif self.backend!='python':raise ValueError('unknown backend')
    def mv(self,v):
        if self._matrix is None:return matrix_vector(self.H,v,self.p)
        import numpy as np
        return tuple(int(x) for x in ((self._matrix@np.asarray(v,dtype=np.int64))%self.p))
    def commitment(self,tag,salt,*parts):
        if type(tag)!=int or tag not in (1,2,3) or type(salt)!=bytes or len(salt)!=32:raise ValueError('invalid commitment input')
        h=self._prefix.copy();h.update(frame(bytes([tag]),*parts,salt));return h.digest()

class Prover:
    """Validate and copy the witness once; it is immutable during this session."""
    def __init__(self,ctx,e):
        if not vector(e,ctx.n,ctx.p) or counts(e,ctx.p)!=ctx.C or ctx.mv(e)!=ctx.s:raise ValueError('invalid witness')
        self.ctx=ctx;self._e=tuple(e)
    def commit(self):
        ctx=self.ctx;p,n=ctx.p,ctx.n;rng=FastRandom(1024)
        try:
            sigma=list(range(n));rng.shuffle(sigma);sigma=tuple(sigma)
            u=tuple(rng.vector(p,n))
        finally:rng.discard()
        salts=tuple(secrets.token_bytes(32) for _ in range(3));z=add(self._e,u,p)
        A=(ctx.commitment(1,salts[0],pack(sigma,n),pack(ctx.mv(u),p)),
           ctx.commitment(2,salts[1],pack(permute(u,sigma),p)),
           ctx.commitment(3,salts[2],pack(permute(z,sigma),p)))
        return A,State(sigma,u,z,self._e,salts)

@dataclass
class State:
    sigma:tuple
    u:tuple
    z:tuple
    e:tuple
    salts:tuple
    used:bool=False

def respond(state,b):
    if state.used:raise ValueError('commitment already answered')
    # Every response attempt consumes the state, including an invalid challenge.
    state.used=True
    if type(b)!=int or b not in (0,1,2):raise ValueError('invalid challenge')
    if b==0:return {'sigma':state.sigma,'x':state.u,'salts':state.salts[:2]}
    if b==1:return {'sigma':state.sigma,'x':state.z,'salts':(state.salts[0],state.salts[2])}
    return {'v':permute(state.e,state.sigma),'x':permute(state.u,state.sigma),'salts':state.salts[1:]}

def verify(ctx,A,b,response):
    p,n=ctx.p,ctx.n
    try:
        if type(b)!=int or b not in (0,1,2) or not isinstance(A,(tuple,list)) or len(A)!=3 or any(type(a)!=bytes or len(a)!=32 for a in A):return False
        if not isinstance(response,dict) or set(response)!=({'sigma','x','salts'} if b<2 else {'v','x','salts'}):return False
        x=response['x'];salts=response['salts']
        if not vector(x,n,p) or not isinstance(salts,(tuple,list)) or len(salts)!=2 or any(type(r)!=bytes or len(r)!=32 for r in salts):return False
        if b<2:
            sigma=response['sigma']
            if not vector(sigma,n,n) or len(set(sigma))!=n:return False
            h=ctx.mv(x)
            if b==1:h=tuple((a-s)%p for a,s in zip(h,ctx.s))
            checks=((ctx.commitment(1,salts[0],pack(sigma,n),pack(h,p)),A[0]),
                    (ctx.commitment(2 if b==0 else 3,salts[1],pack(permute(x,sigma),p)),A[1 if b==0 else 2]))
        else:
            v=response['v']
            if not vector(v,n,p) or counts(v,p)!=ctx.C:return False
            checks=((ctx.commitment(2,salts[0],pack(x,p)),A[1]),
                    (ctx.commitment(3,salts[1],pack(add(v,x,p),p)),A[2]))
        return all(hmac.compare_digest(a,b) for a,b in checks)
    except (KeyError,TypeError,ValueError,OverflowError):return False

def response_bytes(ctx,b,r):
    return b''.join(r['salts'])+pack(r['sigma'] if b<2 else r['v'],ctx.n if b<2 else ctx.p)+pack(r['x'],ctx.p)

def make_instance(cfg,backend='python'):
    p,n,k=cfg['p'],cfg['n'],cfg['k'];C=cfg['C']
    if any(type(x)!=int for x in (p,n,k)) or p<3 or p%2==0 or any(p%d==0 for d in range(3,math.isqrt(p)+1,2)) or not 0<k<n:
        raise ValueError('invalid parameters')
    if not vector(C,p,n+1) or sum(C)!=n:raise ValueError('invalid C')
    if 'weight' in cfg and sum(c*min(a,p-a) for a,c in enumerate(C))!=cfg['weight']:raise ValueError('wrong Lee weight')
    rng=secrets.SystemRandom()
    while True:
        H=tuple(tuple(rng.randrange(p) for _ in range(n)) for _ in range(n-k))
        if rank(H,p)==n-k:break
    e=[a for a,c in enumerate(C) for _ in range(c)];rng.shuffle(e);e=tuple(e)
    s=matrix_vector(H,e,p)
    return Context(p,n,k,tuple(C),H,s,backend),e

def simulate(ctx,b):
    """Offline HVZK candidate using real dummy commitments, not RO programming."""
    if type(b)!=int or b not in (0,1,2):raise ValueError('invalid challenge')
    p,n=ctx.p,ctx.n;rng=FastRandom()
    try:
        sigma=list(range(n));rng.shuffle(sigma);sigma=tuple(sigma)
        x=tuple(rng.vector(p,n))
        if b==2:
            v=[a for a,c in enumerate(ctx.C) for _ in range(c)];rng.shuffle(v);v=tuple(v)
    finally:rng.discard()
    salts=tuple(secrets.token_bytes(32) for _ in range(3));zero=pack((0,)*n,p)
    if b<2:
        h=ctx.mv(x)
        if b==1:h=tuple((a-s)%p for a,s in zip(h,ctx.s))
        A=(ctx.commitment(1,salts[0],pack(sigma,n),pack(h,p)),
           ctx.commitment(2,salts[1],pack(permute(x,sigma),p) if b==0 else zero),
           ctx.commitment(3,salts[2],pack(permute(x,sigma),p) if b==1 else zero))
        r={'sigma':sigma,'x':x,'salts':salts[:2] if b==0 else (salts[0],salts[2])}
    else:
        A=(ctx.commitment(1,salts[0],pack(tuple(range(n)),n),pack((0,)*(n-ctx.k),p)),
           ctx.commitment(2,salts[1],pack(x,p)),ctx.commitment(3,salts[2],pack(add(v,x,p),p)))
        r={'v':v,'x':x,'salts':salts[1:]}
    return A,r
