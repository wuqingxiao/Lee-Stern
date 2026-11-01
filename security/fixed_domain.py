"""Literal fixed-free-position list: F_p^q times a restricted Lee sphere.

No mask multiplicities. This is NOT the complete permutation closure.
"""
from collections import Counter
from functools import lru_cache


@lru_cache(maxsize=8192)
def geodesic_coefficient(magnitudes,a):
    dp=[1]+[0]*a
    for weight in magnitudes:
        if not weight:continue
        new=[0]*(a+1);window=0
        for i in range(a+1):
            window+=dp[i]
            if i>weight:window-=dp[i-weight-1]
            new[i]=window
        dp=new
    return dp[a]


class FixedDomain:
    family='fixed_E'
    def __init__(self,p,h,a,q,radius,*,max_cells=500000):
        if any(type(x)!=int for x in (p,h,a,q,radius)):raise ValueError('integer domain arguments')
        if not(p>=3 and 1<=q<h and a>0 and 1<=radius<=p//2):raise ValueError('fixed-E domain')
        self.Q=p;self.h=h;self.a=a;self.q=q;self.radius=radius;self.low=h-q
        if (self.low+1)*(a+1)>max_cells:raise ValueError('fixed-E DP cell analysis limit')
        self.alphabet=tuple(range(radius+1))+tuple(range(p-radius,p))
        self.weights=tuple(min(x,p-x) for x in self.alphabet);self.stats=Counter()
        rows=[[1]+[0]*a];peak=1
        for _ in range(self.low):
            previous=rows[-1];row=[0]*(a+1)
            for w in range(a+1):
                for weight in self.weights:
                    if weight<=w:row[w]+=previous[w-weight];self.stats['setup_bigint_additions']+=1
            rows.append(row);peak=max(peak,max(row))
        self.dp=rows;self.sphere=rows[-1][a]
        if not self.sphere:raise ValueError('empty fixed-E sphere')
        self.size=p**q*self.sphere;self.groups=(None,)
        self.setup_operand_max_bits=max(peak.bit_length(),self.size.bit_length())+h.bit_length()
        self.table_bytes=(self.low+1)*(a+1)*max(1,(peak.bit_length()+7)//8)
        self.stats['DP_cells']=(self.low+1)*(a+1);self.stats['setup_bigint_multiplications']=q

    def unrank(self,rank,*,counters=None):
        def inc(key,x=1):
            if counters is not None:counters[key]+=x
        inc('rank_bigint_comparisons',2)
        if type(rank)!=int or not 0<=rank<self.size:raise ValueError('fixed-E rank')
        inc('rank_bigint_divmod');free,rank=divmod(rank,self.sphere)
        prefix=[0]*self.q
        for i in range(self.q-1,-1,-1):inc('rank_bigint_divmod');free,prefix[i]=divmod(free,self.Q)
        out=[];weight=self.a
        for remaining in range(self.low-1,-1,-1):
            for symbol,w in zip(self.alphabet,self.weights):
                block=self.dp[remaining][weight-w] if w<=weight else 0
                inc('rank_bigint_comparisons')
                if rank<block:out.append(symbol);weight-=w;break
                inc('rank_bigint_subtractions');rank-=block
            else:raise ArithmeticError('invalid fixed-E unranking')
        if rank or weight or free:raise ArithmeticError('fixed-E residue')
        return tuple(prefix+out)

    def rank(self,word):
        if len(word)!=self.h or any(type(x)!=int or not 0<=x<self.Q for x in word):raise ValueError('word')
        if any(x not in self.alphabet for x in word[self.q:]) or sum(min(x,self.Q-x) for x in word[self.q:])!=self.a:
            raise ValueError('word outside fixed-E list')
        free=0
        for x in word[:self.q]:free=free*self.Q+x
        rank=0;weight=self.a
        for remaining,value in zip(range(self.low-1,-1,-1),word[self.q:]):
            for symbol,w in zip(self.alphabet,self.weights):
                if symbol==value:weight-=w;break
                if w<=weight:rank+=self.dp[remaining][weight-w]
        return free*self.sphere+rank

    def half_count(self,D):
        if len(D)!=self.Q or sum(D)!=self.low or any(type(x)!=int or x<0 for x in D):raise ValueError('low composition')
        if self.Q<=4*self.radius:raise ValueError('fixed-E no-wrap proof not applicable')
        if any(x and min(i,self.Q-i)>self.radius for i,x in enumerate(D)):raise ValueError('low target alphabet')
        magnitudes=tuple(min(i,self.Q-i) for i,x in enumerate(D) for _ in range(x))
        if sum(magnitudes)!=2*self.a or not any(w%2 for w in magnitudes):
            raise ValueError('weight or fixed-point-free proof not applicable')
        return self.Q**self.q*geodesic_coefficient(magnitudes,self.a)
