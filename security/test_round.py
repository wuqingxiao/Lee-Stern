from collections import Counter
from fractions import Fraction as F
import itertools
import json
import math
from pathlib import Path
import random
import unittest

import final_model as model
from fixed_domain import FixedDomain
from execution import execute
from reference_kernel import Instance,prepare,mv,rank
from streamed import execute as oracle
from run_round import variants,read,HERE


class RoundTests(unittest.TestCase):
    def test_fixed_E_count_rank_and_complete_representations(self):
        for p,h,a,q,r,targetlow in ((7,3,1,1,1,(1,1)),(11,4,1,1,2,(1,1,0)),
            (11,5,1,2,2,(1,10,0)),(17,5,2,1,3,(1,1,2,0))):
            d=FixedDomain(p,h,a,q,r)
            lows=[x for x in itertools.product(d.alphabet,repeat=h-q) if sum(min(t,p-t) for t in x)==a]
            actual={free+x for free in itertools.product(range(p),repeat=q) for x in lows}
            self.assertEqual(len(actual),d.size)
            decoded={d.unrank(i) for i in range(d.size)};self.assertEqual(decoded,actual)
            for i in range(d.size):self.assertEqual(d.rank(d.unrank(i)),i)
            target=(0,)*q+targetlow
            full=sum(tuple((x-y)%p for x,y in zip(target,v)) in actual for v in actual)
            D=tuple(targetlow.count(i) for i in range(p));self.assertEqual(d.half_count(D),full)
            self.assertFalse(any(tuple(2*x%p for x in v)==target for v in actual))
        with self.assertRaises(ValueError):FixedDomain(11,4,1,1,3).half_count((1,2)+(0,)*9)

    def test_asymmetric_alpha_is_actual_ordered_partition_probability(self):
        C=(4,3)+(0,)*8+(1,);A=(0,2)+(0,)*9;B=(0,1)+(0,)*8+(1,)
        e=(0,1,1,0,1,10,0,0);hits=0
        for order in itertools.permutations(range(8)):
            w=tuple(e[i] for i in order)
            if w[0] or w[3]:continue
            if tuple(w[1:3].count(i) for i in range(11))!=A:continue
            if tuple(w[4:6].count(i) for i in range(11))!=B:continue
            hits+=1
        d=FixedDomain(11,3,1,1,1);cfg=dict(name='asym',p=11,n=8,k=2,C=C)
        t=dict(family='fixed_E',ell=4,a=1,q=1,radius=1,left=A,right=B)
        row=model.assess(cfg,t,d,6,1,4,'resident',2,2**30)
        self.assertEqual(F(hits,math.factorial(8)),F(1,140));self.assertEqual(F(row['alpha']),F(1,140))
        self.assertTrue(row['asymmetric']);self.assertGreater(F(row['success_lower']),F(1,2))

    def test_asymmetric_conditional_moments_all_rank_one_maps(self):
        p=5;base=((0,1),(1,0),(4,0),(2,1),(0,4));target=(1,1)
        left=[(base[0],base[1]),(base[1],base[0])]
        right=[(x,tuple((a-b)%p for a,b in zip(target,x))) for x in base]
        right=[pair for pair in right if pair[1] in base]
        whole=[(a[0]+b[0],frozenset(a),frozenset(b)) for a,b in itertools.product(left,right)]
        for L in (2,4,5):
            subsets=[frozenset(v) for v in itertools.combinations(base,L)]
            included=[[x for x,A,B in whole if A<=S and B<=T] for S,T in itertools.product(subsets,repeat=2)]
            total=total2=maps=0
            for G in itertools.product(range(p),repeat=4):
                if not any(G):continue
                maps+=1
                for vectors in included:
                    counts=Counter(sum(x*y for x,y in zip(G,v))%p for v in vectors)
                    total+=sum(counts.values());total2+=sum(x*x for x in counts.values())
            mean,second,_=model.old.offset_moments(p,4,1,5,L,2,4)
            self.assertEqual(F(total,maps*len(included)*p),mean)
            self.assertEqual(F(total2,maps*len(included)*p),second)

    def test_templates_cover_valid_larger_ranges_without_changing_C(self):
        config=read(HERE/'config.json');cfgs=read(HERE/'inputs/candidates.json');count=0
        for cfg in cfgs:
            for ell in config['fixed_E']['ells']:
                for r in config['fixed_E']['radii']:
                    for q in config['fixed_E']['q']:
                        for t in variants(cfg,ell,r,q,'fixed_E',config['fixed_E']['asymmetry_shifts']):
                            A=t['left'];B=t['right'];h=(cfg['k']+ell)//2
                            self.assertEqual(sum(A),h-q);self.assertEqual(sum(B),h-q)
                            self.assertEqual(sum(min(i,cfg['p']-i)*x for i,x in enumerate(A)),2*t['a'])
                            self.assertEqual(sum(min(i,cfg['p']-i)*x for i,x in enumerate(B)),2*t['a'])
                            self.assertTrue(all(x+y+(2*q if i==0 else 0)<=v for i,(x,y,v) in enumerate(zip(A,B,cfg['C']))));count+=1
        self.assertGreater(count,500)

    def test_repaired_packed_execution_on_actual_fixed_E_vectors(self):
        p=11;d=FixedDomain(p,3,1,1,1);e=(0,1,1,0,1,10,0,0);C=(4,3)+(0,)*8+(1,)
        H=((0,0,0,0,0,0,1,0),(0,0,0,0,0,0,0,1))+tuple(tuple(int(j==i) for j in range(8)) for i in (0,1,3,4))
        prepared=prepare(Instance(p,H,mv(H,e,p),C),range(6),1,1,1,1)
        ranks=tuple(range(6));words=tuple(d.unrank(i) for i in ranks)
        for mode in ('resident','tiled_partial'):
            for tau in (0,1,2):
                reference=oracle(prepared,words,words,2,offset=(tau,),all_solutions=True)
                result=execute(prepared,d,ranks,ranks,1,2,2,(tau,),mode,True)
                self.assertEqual(result['solutions'],reference['multiplicities']);self.assertEqual(result['final_active_areas'],{})
                self.assertLessEqual(result['observed_active_bytes'],result['reserved_bytes'])

    def test_fixed_E_unranking_and_stream_means_over_18720_maps(self):
        p=3;d=FixedDomain(p,2,1,1,1);cfg=dict(p=p,n=6,k=2);rng=random.Random(438)
        left=tuple(rng.sample(range(d.size),3));right=tuple(rng.sample(range(d.size),3))
        mem=model.layout(cfg,2,1,d.size,3,2,d,'tiled_partial',2,2**30)
        ops,setup,_=model.allowances(cfg,2,1,3,2,d,mem)
        vectors=[x for x in itertools.product(range(p),repeat=4) if any(x)]
        base=((0,0,0,0,1,0),(0,0,0,0,0,1));e=(1,1,1,1,0,0);total=Counter();n=0
        for low in itertools.product(vectors,repeat=2):
            if rank(low,p)!=2:continue
            H=base+tuple(x+(0,0) for x in low)
            prepared=prepare(Instance(p,H,mv(H,e,p),(2,4,0)),range(4),1,1,1,1)
            for tau in range(p):
                row=execute(prepared,d,left,right,1,2,2,(tau,),'tiled_partial',True)
                total.update(row['resources']);n+=1
        self.assertEqual(n,18720)
        for key,value in total.items():self.assertLessEqual(F(value,n),ops.get(key,F(0))+setup.get(key,F(0)),key)


if __name__=='__main__':
    r=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(RoundTests))
    (HERE/'TEST_RESULT.json').write_text(json.dumps(dict(passed=r.wasSuccessful(),tests=r.testsRun,
        matrix_offset_cases=18720,target_attack_executed=False),indent=2),encoding='utf-8')
    raise SystemExit(not r.wasSuccessful())
