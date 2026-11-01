"""Literal prose permutation closure, with vectors as unique elements.

Supports residue alphabets; the attack kernel separately restricts PGE to fields.
Every enumerator has an explicit preflight capacity check, never truncation.
"""
from functools import lru_cache
import itertools
import math


class ReferenceLimit(ValueError):pass


def check_parameters(Q,h,a,q,radius):
    if any(type(x)!=int for x in (Q,h,a,q,radius)):
        raise ValueError('integer domain parameters required')
    if Q<2 or h<1 or a<0 or not 0<=q<=h or not 0<=radius<=Q//2:
        raise ValueError('invalid domain parameters')


def lee(x,Q):return min(x,Q-x)


def histogram_mask_count(hist,ordinary_count,a,radius):
    """[t^ordinary_count z^a] product_(w<=r) (1+t z^w)^hist[w]."""
    dp={(0,0):1}
    for w,count in enumerate(hist):
        if w>radius:break
        out={}
        for (size,weight),ways in dp.items():
            for take in range(min(count,ordinary_count-size)+1):
                new_weight=weight+take*w
                if new_weight>a:break
                key=(size+take,new_weight)
                out[key]=out.get(key,0)+ways*math.comb(count,take)
        dp=out
    return dp.get((ordinary_count,a),0)


def mask_multiplicity(word,Q,a,q,radius):
    check_parameters(Q,len(word),a,q,radius)
    if any(type(x)!=int or not 0<=x<Q for x in word):raise ValueError('invalid word')
    hist=[0]*(Q//2+1)
    for x in word:hist[lee(x,Q)]+=1
    return histogram_mask_count(hist,len(word)-q,a,radius)


@lru_cache(maxsize=64)
def unique_domain(Q,h,a,q,radius,max_words=100000):
    check_parameters(Q,h,a,q,radius)
    if type(max_words)!=int or max_words<1:raise ValueError('invalid word limit')
    if Q**h>max_words:raise ReferenceLimit('domain word scan exceeds explicit reference limit')
    return tuple(word for word in itertools.product(range(Q),repeat=h)
                 if mask_multiplicity(word,Q,a,q,radius)>0)


def histograms(length,bins):
    if bins==1:yield (length,);return
    for first in range(length+1):
        for tail in histograms(length-first,bins-1):yield (first,)+tail


def histogram_domain_size(Q,h,a,q,radius,max_histograms=100000):
    """Exact distinct-vector count; not a generation, sampling or attack cost."""
    check_parameters(Q,h,a,q,radius)
    bins=Q//2+1;states=math.comb(h+bins-1,bins-1)
    if type(max_histograms)!=int or max_histograms<1:raise ValueError('invalid histogram limit')
    if states>max_histograms:raise ReferenceLimit('histogram count exceeds explicit reference limit')
    total=0;accepted=0
    multiplicity=[sum(lee(x,Q)==w for x in range(Q)) for w in range(bins)]
    for hist in histograms(h,bins):
        if not histogram_mask_count(hist,h-q,a,radius):continue
        accepted+=1;arrangements=math.factorial(h)//math.prod(math.factorial(m) for m in hist)
        total+=arrangements*math.prod(g**m for g,m in zip(multiplicity,hist))
    return {'distinct_vectors':total,'histograms_examined':states,'histograms_accepted':accepted,
        'scope':'exact cardinality only; no target-size generation or attack-cost claim'}
