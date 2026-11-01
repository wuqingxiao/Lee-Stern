"""Exact counting identities; no protocol measurements or attack search."""
from collections import Counter,defaultdict
from itertools import combinations,product
from math import comb,factorial,log2,prod
from pathlib import Path
import csv,hashlib,json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
def lee(a,p):return min(a,p-a)
def bounded_words(C,p):
    dp={(0,0):1}
    for a,cap in enumerate(C):
        nxt=defaultdict(int)
        for (d,w),v in dp.items():
            for j in range(cap+1):nxt[d+j,w+j*lee(a,p)]+=v*comb(d+j,j)
        dp=dict(nxt)
    return dp
def subset_counts(classes):
    dp={(0,0):1}
    for weight,cap in classes:
        nxt=defaultdict(int)
        for (d,w),v in dp.items():
            for j in range(cap+1):nxt[d+j,w+j*weight]+=v*comb(cap,j)
        dp=dict(nxt)
    return dp
def orbit(C):return factorial(sum(C))//prod(factorial(c) for c in C)
cases=[(3,(1,1,1)),(3,(2,1,1)),(3,(3,0,0)),(5,(1,1,0,0,1)),(5,(1,2,0,0,1)),(7,(1,1,1,0,0,1,1))]
results=[]
for p,C in cases:
    n=sum(C);w=sum(lee(a,p)*c for a,c in enumerate(C));dp=bounded_words(C,p)
    brute=Counter()
    examined=0
    for d in range(n+1):
        for v in product(range(p),repeat=d):
            examined+=1;freq=Counter(v)
            if all(freq[a]<=C[a] for a in range(p)):brute[d,sum(lee(a,p) for a in v)]+=1
    assert all(dp.get(k,0)==brute.get(k,0) for k in set(dp)|set(brute))
    assert dp[n,w]==orbit(C)
    base=[a for a,c in enumerate(C) for _ in range(c)]
    proj=subset_counts([(lee(a,p),c) for a,c in enumerate(C)])
    brute_proj=Counter()
    for d in range(n+1):
        for I in combinations(range(n),d):brute_proj[d,sum(lee(base[i],p) for i in I)]+=1
    assert all(proj.get(k,0)==brute_proj.get(k,0) for k in set(proj)|set(brute_proj))
    T=[C[0]]+[C[i]+C[p-i] for i in range(1,(p+1)//2)]
    folded=subset_counts(list(enumerate(T)))
    assert proj==folded
    assert all(sum(v for (size,_),v in proj.items() if size==d)==comb(n,d) for d in range(n+1))
    results.append(dict(p=p,C=C,n=n,examined_words=examined,bounded_word_entries=len(dp),subset_entries=len(proj),normalization=dp[n,w],status='PASS'))

source=ROOT/'security/inputs/candidates.json'
rows=[]
for ix,c in enumerate(json.loads(source.read_text('utf-8')),1):
    C=c['C'];p=c['p'];n=c['n'];r=(p-1)//2
    T=[C[0]]+[C[a]+C[p-a] for a in range(1,r+1)]
    full=orbit(C);folded=factorial(n)//prod(factorial(m) for m in T)
    sign_factor=prod(comb(C[a]+C[p-a],C[a]) for a in range(1,r+1))
    amplitude=folded*2**(n-C[0])
    assert full==folded*sign_factor and amplitude>=full
    C_sha=hashlib.sha256(json.dumps(C,separators=(',',':')).encode()).hexdigest()
    rows.append(dict(candidate=f'C{ix}',recipe=c['primary_recipe'],C_sha256=C_sha,n=n,k=c['k'],p=p,
       weight=sum(lee(a,p)*v for a,v in enumerate(C)),zero_count=C[0],
       exact_SC=str(full),exact_folded=str(folded),exact_ST=str(amplitude),exact_fiber=str(sign_factor),
       log2_SC=log2(full),log2_folded=log2(folded),log2_ST=log2(amplitude),
       delta_folded=log2(sign_factor),delta_T_to_C=log2(amplitude)-log2(full)))
report=dict(date='2026-10-26',scope='finite exact combinatorial identities and inherited candidate cardinalities',
            cases=results,candidate_rows=rows,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
            counts_are_security_bits=False,new_performance_measurements=0,new_attack_runs=0,new_parameter_searches=0)
(HERE/'COUNTING_CHECKS.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n','utf-8')
fields=['candidate','recipe','C_sha256','n','k','p','weight','zero_count','log2_SC','log2_folded','delta_folded','delta_T_to_C']
with (HERE/'candidate_count_comparison.csv').open('w',encoding='utf-8',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');writer.writeheader();writer.writerows(rows)
print(json.dumps({'small_cases':len(results),'brute_words_examined':sum(x['examined_words'] for x in results),'candidate_rows':3,'status':'PASS'},ensure_ascii=False))
