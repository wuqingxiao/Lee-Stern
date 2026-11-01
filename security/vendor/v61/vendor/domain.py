"""Exact union-domain partition by low Lee histogram plus high-coordinate count.

No mask labels, inclusion-exclusion over masks, or vector-space scan. Prefix
blocks give a bijection between [0,|B|) and distinct vectors in the source-prose
permutation closure. Preprocessing limits are explicit, never truncation.
"""
from collections import Counter
from functools import lru_cache
import math
import secrets


class CapacityLimit(ValueError):pass


class Domain:
    def __init__(self,Q,h,a,q,radius,*,max_candidates=2000000):
        if any(type(x)!=int for x in (Q,h,a,q,radius,max_candidates)):
            raise ValueError('integer parameters required')
        if Q<2 or h<1 or a<0 or not 0<=q<=h or not 0<=radius<=Q//2 or max_candidates<1:
            raise ValueError('invalid domain parameters')
        self.Q=Q;self.h=h;self.a=a;self.q=q;self.radius=radius
        self.alphabets=tuple(tuple(x for x in range(Q) if min(x,Q-x)==w) for w in range(radius+1))
        self.alphabets+=(tuple(x for x in range(Q) if min(x,Q-x)>radius),)
        self.symbol_index={x:(category,i) for category,alphabet in enumerate(self.alphabets) for i,x in enumerate(alphabet)}
        self.stats=Counter();self.groups=[];self.ends=[];self.lookup={};self.size=0
        factorial=[1]
        for i in range(1,h+1):
            factorial.append(factorial[-1]*i);self.stats['setup_bigint_multiplications']+=1
        powers=[]
        for alphabet in self.alphabets:
            values=[1]
            for _ in range(h):
                values.append(values[-1]*len(alphabet));self.stats['setup_bigint_multiplications']+=1
            powers.append(values)

        @lru_cache(None)
        def removal_weights(capped,t):
            self.stats['removal_signature_cache_misses']+=1
            dp=[0]*(t+1)
            for count in range(min(capped[0],t)+1):dp[count]=1
            for w,c in enumerate(capped[1:],1):
                nxt=[0]*(t+1)
                for used,bits in enumerate(dp):
                    if not bits:continue
                    for take in range(min(c,t-used)+1):
                        nxt[used+take]|=bits<<(w*take)
                        self.stats['removal_DP_shift_or_updates']+=1
                dp=nxt
            return dp[t]

        def accept(nonzero,weight,high):
            self.stats['candidate_histograms']+=1
            if self.stats['candidate_histograms']>max_candidates:
                raise CapacityLimit('candidate histogram limit; no partial domain returned')
            zero=h-high-sum(nonzero)
            if zero<0:return
            t=q-high;delta=weight-a
            if delta<0 or delta>t*radius:return
            counts=(zero,)+tuple(nonzero)+(high,)
            self.stats['removal_signature_queries']+=1
            bits=removal_weights(tuple(min(c,t) for c in counts[:-1]),t)
            if not bits>>delta&1:return
            denominator=1;symbols=1
            for c,power in zip(counts,powers):
                denominator*=factorial[c];symbols*=power[c]
                self.stats['setup_bigint_multiplications']+=2
            arrangements=factorial[h]//denominator;block=arrangements*symbols
            self.stats['setup_bigint_divisions']+=1;self.stats['setup_bigint_multiplications']+=1
            if not block:return
            start=self.size;self.size+=block;self.stats['setup_bigint_additions']+=1
            index=len(self.groups);self.lookup[counts]=index
            self.groups.append((counts,arrangements,symbols,start,self.size));self.ends.append(self.size)

        def recurse(w,remaining,weight,high,path):
            self.stats['histogram_DFS_nodes']+=1
            if self.stats['histogram_DFS_nodes']>max_candidates*(radius+2):
                raise CapacityLimit('histogram DFS limit; no partial domain returned')
            upper=a+(q-high)*radius
            if w>radius:
                if weight>=a:accept(path,weight,high)
                return
            for count in range(min(remaining,(upper-weight)//w)+1):
                recurse(w+1,remaining-count,weight+w*count,high,path+(count,))

        for high in range(q+1):
            if high and not self.alphabets[-1]:continue
            if radius==0:accept((),0,high)
            else:recurse(1,h-high,0,high,())
        self.stats['groups']=len(self.groups)
        self.stats['removal_cache_entries']=removal_weights.cache_info().currsize
        # This count includes each group's 5 integer fields, including counts,
        # and an independent end-prefix index. It is packed payload, not Python RAM.
        count_bytes=max(1,(h.bit_length()+7)//8)
        self.table_payload_bytes=sum(len(c)*count_bytes+sum(max(1,(v.bit_length()+7)//8) for v in (A,S,start,end,end))
            for c,A,S,start,end in self.groups)
        self.removal_cache_payload_upper=(radius+1)*count_bytes+max(1,((q*radius+1)+7)//8)
        self.removal_cache_payload_upper*=self.stats['removal_cache_entries']
        self.setup_payload_bytes=self.table_payload_bytes+self.removal_cache_payload_upper
        index_bytes=max(1,(len(self.groups).bit_length()+7)//8)
        self.rank_lookup_payload_bytes=len(self.groups)*((radius+2)*count_bytes+index_bytes)
        alphabet_index_bytes=max(1,(max(Q,h).bit_length()+7)//8)
        self.alphabet_payload_bytes=Q*(max(1,((Q-1).bit_length()+7)//8)+2*alphabet_index_bytes)
        self.persistent_payload_bytes=self.table_payload_bytes+self.rank_lookup_payload_bytes+self.alphabet_payload_bytes
        self.setup_work_payload_bytes=(self.removal_cache_payload_upper+
            sum(max(1,(x.bit_length()+7)//8) for x in factorial)+
            sum(max(1,(x.bit_length()+7)//8) for values in powers for x in values))
        self.setup_payload_bytes=self.persistent_payload_bytes+self.setup_work_payload_bytes
        self.setup_operand_max_bits=max(factorial[-1].bit_length(),max(v[-1].bit_length() for v in powers))+max(h.bit_length(),Q.bit_length())

    def _locate(self,rank,counters):
        lo=0;hi=len(self.ends)
        while lo<hi:
            mid=(lo+hi)//2
            if counters is not None:
                counters['rank_prefix_comparisons']+=1
                counters['rank_bigint_comparisons']+=1
            if rank<self.ends[mid]:hi=mid
            else:lo=mid+1
        return lo

    def unrank(self,rank,*,counters=None):
        if counters is not None:
            counters['rank_input_bound_comparisons']+=2
            counters['rank_bigint_comparisons']+=2
        if type(rank)!=int or not 0<=rank<self.size:raise ValueError('domain rank')
        index=self._locate(rank,counters);counts,A,S,start,_=self.groups[index]
        if counters is not None:
            counters['rank_bigint_subtractions']+=1;counters['rank_bigint_divmod']+=1
        arrangement,symbols=divmod(rank-start,S);counts=list(counts);labels=[]
        for remaining in range(self.h,0,-1):
            for category,count in enumerate(counts):
                if not count:continue
                if counters is not None:
                    counters['rank_bigint_multiplications']+=1;counters['rank_bigint_divisions']+=1
                block=A*count//remaining
                if counters is not None:
                    counters['rank_category_comparisons']+=1
                    counters['rank_bigint_comparisons']+=1
                if arrangement>=block:
                    if counters is not None:counters['rank_bigint_subtractions']+=1
                    arrangement-=block
                else:
                    labels.append(category);counts[category]-=1;A=block;break
            else:raise AssertionError('invalid multiset unranking')
        values=[0]*self.h
        for i in range(self.h-1,-1,-1):
            if counters is not None:counters['rank_bigint_divmod']+=1
            alphabet=self.alphabets[labels[i]];symbols,digit=divmod(symbols,len(alphabet))
            values[i]=alphabet[digit]
        if symbols or arrangement:raise AssertionError('unranking residual')
        return tuple(values)

    def rank(self,word):
        if len(word)!=self.h or any(type(x)!=int or x not in self.symbol_index for x in word):
            raise ValueError('invalid word')
        labels=[self.symbol_index[x][0] for x in word]
        counts=tuple(labels.count(i) for i in range(self.radius+2))
        index=self.lookup.get(counts)
        if index is None:raise ValueError('word outside union domain')
        _,A,S,start,_=self.groups[index];counts=list(counts);arrangement=0;symbols=0
        for remaining,category in zip(range(self.h,0,-1),labels):
            for before in range(category):arrangement+=A*counts[before]//remaining
            A=A*counts[category]//remaining;counts[category]-=1
        for x,category in zip(word,labels):
            symbols=symbols*len(self.alphabets[category])+self.symbol_index[x][1]
        return start+arrangement*S+symbols

    def generate(self,start=0,stop=None,*,max_records=100000,counters=None):
        stop=self.size if stop is None else stop
        if not 0<=start<=stop<=self.size:raise ValueError('rank interval')
        if stop-start>max_records:raise CapacityLimit('explicit generation limit')
        for rank in range(start,stop):yield self.unrank(rank,counters=counters)

    def sample_ranks(self,count,rng=None,*,max_records=100000,counters=None):
        if type(count)!=int or not 0<=count<=self.size:raise ValueError('subset size')
        if count>max_records:raise CapacityLimit('explicit sample limit')
        rng=rng if rng is not None else secrets.SystemRandom();seen=set();out=[]
        for j in range(self.size-count,self.size):
            if counters is not None:
                counters['sample_uniform_integer_calls']+=1;counters['sample_set_operations']+=2
            candidate=rng.randrange(j+1);chosen=j if candidate in seen else candidate
            seen.add(chosen);out.append(chosen)
        return tuple(out)

    def sample(self,count,rng=None,*,max_records=100000,counters=None):
        return tuple(self.unrank(r,counters=counters) for r in self.sample_ranks(count,rng,max_records=max_records,counters=counters))

    def generation_allowances(self,records,amortized=True):
        # Counts of explicitly defined high-level integer operations, not bit cost.
        h=self.h;categories=self.radius+2;G=len(self.groups)
        return {'sample_uniform_integer_calls':records if amortized else 0,
            'sample_set_operations':4*records if amortized else 0,
            'rank_prefix_comparisons':records*G.bit_length(),
            'rank_bigint_multiplications':records*h*categories,
            'rank_bigint_divisions':records*h*categories,
            'rank_bigint_subtractions':records*(1+h*categories),
            'rank_bigint_divmod':records*(h+1)}
