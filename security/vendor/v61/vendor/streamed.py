"""Small explicit-list realization of the distinct-vector bounded-block schedule.

Prepared instances use the reference's normalization of source LAST u rows.
Counters are defined operation/reservation units, not CPU instructions or RAM.
"""
from collections import Counter
import itertools


def sort_rows(rows,key,c,name):
    a=list(rows);b=[None]*len(a);width=1
    while width<len(a):
        for start in range(0,len(a),2*width):
            mid=min(start+width,len(a));end=min(start+2*width,len(a));i=start;j=mid
            for out in range(start,end):
                take_left=j>=end
                if i<mid and j<end:c[name]+=1;take_left=key(a[i])<=key(a[j])
                if take_left:b[out]=a[i];i+=1
                else:b[out]=a[j];j+=1
        a,b=b,a;width*=2
    return a


def equal_range(rows,target,key,c,name):
    lo=0;hi=len(rows)
    while lo<hi:
        mid=(lo+hi)//2;c[name]+=1
        if key(rows[mid])<target:lo=mid+1
        else:hi=mid
    first=lo;hi=len(rows)
    while lo<hi:
        mid=(lo+hi)//2;c[name]+=1
        if target<key(rows[mid]):hi=mid
        else:lo=mid+1
    return range(first,lo)


def execute(prepared,left,right,block_records=1024,*,offset=None,all_solutions=False,max_records=100000):
    if block_records<1 or min(len(left),len(right))<1:raise ValueError('positive sizes')
    if prepared.instance.C is None:raise ValueError('costed schedule requires exact C')
    if max(len(left),len(right))>max_records:raise ValueError('explicit list limit')
    if len(set(left))!=len(left) or len(set(right))!=len(right):raise ValueError('distinct vector inputs required')
    p=prepared.instance.p;N=len(prepared.I);h=N//2;ell=len(prepared.B);u=prepared.u
    if offset is None:offset=(0,)*u
    offset=tuple(offset)
    if len(offset)!=u or any(type(x)!=int or not 0<=x<p for x in offset):raise ValueError('offset')
    second_target=tuple((x-y)%p for x,y in zip(prepared.target[:u],offset))
    n=prepared.instance.n;r=len(prepared.instance.H);j=len(prepared.J);c=Counter();leaves=[]
    for half,words in enumerate((left,right)):
        records=[]
        for word in words:
            if len(word)!=h:raise ValueError('half vector length')
            c['base_key_field_multiplications']+=h*ell;c['base_key_field_additions_upper']+=h*ell
            key=tuple(sum(row[half*h+i]*x for i,x in enumerate(word))%p for row in prepared.B)
            records.append((word,key))
        leaves.append(records)
    first=sort_rows(leaves[0],lambda rec:rec[1][:u],c,'base_sort_key_comparisons_upper')
    def merge(target,kind):
        if kind=='L2':c['L2_regeneration_passes_expected_upper']+=1
        for right,kr in leaves[1]:
            c['first_join_queries_expected_upper']+=1;c['first_join_target_field_subtractions_expected_upper']+=u
            wanted=tuple((x-y)%p for x,y in zip(target,kr[:u]))
            for index in equal_range(first,wanted,lambda rec:rec[1][:u],c,'first_join_key_comparisons_expected_upper'):
                left,kl=first[index]
                c['L1_records_expected_upper' if kind=='L1' else 'L2_generated_with_regeneration_expected_upper']+=1
                c['mid_vector_writes_expected_upper']+=N;c['mid_syndrome_additions_expected_upper']+=ell-u
                yield left+right,tuple((x+y)%p for x,y in zip(kl[u:],kr[u:]))
    stream=merge(offset,'L1');outputs=Counter()
    while True:
        block=list(itertools.islice(stream,block_records))
        if not block:break
        ordered=sort_rows(block,lambda rec:rec[1],c,'root_block_sort_comparisons_expected_upper')
        for y,ky in merge(second_target,'L2'):
            c['root_target_subtractions_expected_upper']+=ell-u
            wanted=tuple((x-z)%p for x,z in zip(prepared.target[u:],ky))
            for index in equal_range(ordered,wanted,lambda rec:rec[1],c,'root_boundary_comparisons_expected_upper'):
                c['raw_root_matches_expected_upper']+=1;c['root_sum_additions_expected_upper']+=N
                c['root_Lee_tests_expected_upper']+=N;c['root_vector_writes_expected_upper']+=3*n
                inside=tuple((x+z)%p for x,z in zip(ordered[index][0],y))
                if sum(min(x,p-x) for x in inside)!=prepared.v:continue
                c['completion_field_multiplications_expected_upper']+=j*N+r*n
                c['completion_field_additions_expected_upper']+=j*(N+1)+r*n
                c['histogram_increments_expected_upper']+=n;c['histogram_init_compare_expected_upper']+=2*p
                candidate=prepared.complete(inside)
                # A single histogram scan matches the O(n+p) fee model. Do
                # not call the slow reference's p separate word.count scans.
                histogram=[0]*p
                for value in candidate:histogram[value]+=1
                composition_ok=all(x==y for x,y in zip(histogram,prepared.instance.C))
                syndrome_ok=(tuple(sum(x*y for x,y in zip(row,candidate))%p
                    for row in prepared.instance.H)==prepared.instance.s) if composition_ok else False
                if syndrome_ok:
                    outputs[candidate]+=1
                    if not all_solutions:return {'solutions':tuple(sorted(outputs)),'multiplicities':dict(outputs),'operations':dict(c)}
    return {'solutions':tuple(sorted(outputs)),'multiplicities':dict(outputs),'operations':dict(c)}
