"""Re-derived GEODESIC SUBFAMILY lower bound inside the unique-vector domain.

Not the full representation set: only this subfamily has the 2/4 footprint
classification. Domain cardinality is the new exact |B|, never a tagged count.
"""
from fractions import Fraction as F
import math


def orbit(C):
    left=sum(C);out=1
    for c in C:out*=math.comb(left,c);left-=c
    return out


def half_representations(D,a,q,p):
    if D[0]<q:return 0
    if sum(c*min(i,p-i) for i,c in enumerate(D))!=2*a:raise ValueError('half target weight')
    dp=[0]*(a+1);dp[0]=1
    for symbol,count in enumerate(D):
        w=min(symbol,p-symbol)
        if not w:continue
        for _ in range(count):
            out=[0]*(a+1)
            for old,c in enumerate(dp):
                for added in range(min(w,a-old)+1):out[old+added]+=c
            dp=out
    return math.comb(D[0],q)*(p-1)**q*dp[a]


def kernel(p,N,rows,d=1):
    if rows>N-d:return F(0)
    return math.prod(F(p**(N-rows)-p**i,p**N-p**i) for i in range(d))


def inclusion(B,L,d):
    if min(B,L)<d:return F(0)
    return math.prod(F(L-i,B-i) for i in range(d))


def subfamily(p,N,u,rL,rR,B,L):
    if p<3 or not 0<L<=B:raise ValueError('odd prime and valid base subset required')
    R=rL*rR
    if not R:return {'R':0,'mean':F(0),'second':F(0),'lower':F(0),'upper':F(0)}
    if min(rL,rR)<2:raise ValueError('fixed point free half family')
    f2=inclusion(B,L,2);f4=inclusion(B,L,4)
    Q1=kernel(p,N,u);Q2=kernel(p,N,u,2)
    counts={(2,2):3*R,(2,4):2*rL*rR*(rR-2),
        (4,2):2*rL*(rL-2)*rR,(4,4):rL*(rL-2)*rR*(rR-2)}
    assert sum(counts.values())==R*(R-1)
    mean=R*f2*f2*Q1
    second=mean+Q2*(counts[2,2]*f2*f2+(counts[2,4]+counts[4,2])*f2*f4+counts[4,4]*f4*f4)
    return {'R':R,'mean':mean,'second':second,'lower':mean*mean/second if second else F(0),
        'upper':min(F(1),mean),'ordered_pair_counts':counts}
