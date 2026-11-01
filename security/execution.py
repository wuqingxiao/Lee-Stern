"""Small packed-record oracle; cap reservations precede counted operations.

PGE preparation/domain preprocessing are separately bounded; Python heap/RSS
and a full target multi-proposal driver are not certified by this kernel.
"""
from collections import Counter
from contextlib import ExitStack
import final_model as model
from memory_contract import Arena,Records,Ledger,InstrumentedCounter,sort_records,equal_range


def prepare_with_reservation(instance,I,a,q,radius,u,ledger):
    """Small PGE oracle with its deterministic upper allowances prepaid.

    Charges are reservations, not a CPU instruction count. Failure still pays
    preparation; no base-list kernel is entered on a failed proposal.
    """
    from reference_kernel import prepare
    n=instance.n;m=len(instance.H);j=n-len(I);width=n+1+m
    for name,value in dict(field_multiplications=j*m*width,
        field_additions_subtractions=j*(m-1)*width,field_inversions=j,PGE_pivot_tests=j*m).items():
        ledger.charge(name,value)
    return prepare(instance,I,a,q,radius,u)


def execute(prepared,domain,left_ranks,right_ranks,u,block=8,tile=8,offset=None,
            mode='resident',all_solutions=False,max_records=10000,*,ledger=None):
    p=prepared.instance.p;N=len(prepared.I);h=N//2;ell=len(prepared.B);n=prepared.instance.n
    if not 1<=u<ell or prepared.u!=u or min(block,tile)<1:raise ValueError('geometry')
    L=len(left_ranks)
    if not 2<=L==len(right_ranks)<=max_records:raise ValueError('equal small subset sizes required')
    if len(set(left_ranks))!=L or len(set(right_ranks))!=L:raise ValueError('distinct ranks')
    if any(type(x)!=int or not 0<=x<domain.size for x in (*left_ranks,*right_ranks)):raise ValueError('rank range')
    offset=(0,)*u if offset is None else tuple(offset)
    if len(offset)!=u or any(type(x)!=int or not 0<=x<p for x in offset):raise ValueError('offset')
    if mode=='tiled_full' and (tuple(left_ranks)!=tuple(range(domain.size)) or tuple(right_ranks)!=tuple(range(domain.size))):
        raise ValueError('implicit full rank range required')
    cfg=dict(p=p,n=n,k=n-len(prepared.instance.H))
    mem=model.layout(cfg,ell,u,domain.size,L,block,domain,mode,tile,2**30)
    block=mem['root_block_records'];tile=mem['tile_records'];f=max(1,((p-1).bit_length()+7)//8)
    ledger=Ledger() if ledger is None else ledger;aux=InstrumentedCounter(ledger);stats=Counter();outputs=Counter()
    ledger.charge('allocation_zero_bytes',mem['reserved_bytes'])
    arena=Arena(mem);streams=[];rankwidth=max(1,(domain.size.bit_length()+7)//8)
    keypart=slice(h,h+u);rootkey=slice(N,N+ell-u)
    def answer():
        # active is intentionally a live mapping: context unwinding empties it
        # even when a result is returned at the first root.
        return dict(solutions=dict(outputs),operations=dict(stats),resources=dict(ledger.counts),unrank_aux=dict(aux),
            reserved_bytes=mem['reserved_bytes'],observed_active_bytes=arena.peak,
            peak_active_areas=arena.peak_names,lifetime_events=arena.events,
            final_active_areas=arena.active,process_RSS_measured=False)
    def rank_at(side,index):
        if mode!='tiled_partial':return (left_ranks,right_ranks)[side][index]
        ledger.charge('rank_array_byte_reads',rankwidth)
        view=arena.view('replayable_rank_arrays');pos=(side*L+index)*rankwidth
        return int.from_bytes(view[pos:pos+rankwidth],'little')
    def record(side,index):
        ledger.charge('streamed_leaf_vectors');stats['leaf_records_generated']+=1
        word=domain.unrank(rank_at(side,index),counters=aux)
        ledger.charge('field_multiplications',h*ell);ledger.charge('field_additions_subtractions',h*ell)
        key=tuple(sum(row[side*h+i]*x for i,x in enumerate(word))%p for row in prepared.B)
        return word+key
    with ExitStack() as static:
        dynamic={'base_left','base_right','base_sort_scratch','root_block','root_sort_scratch',
            'L1_tile','L2_tile','L1_sort_scratch','L2_sort_scratch','L1_right_record','L2_right_record'}
        for name in arena.regions:
            if name not in dynamic:static.enter_context(arena.lease(name))
        if mode=='tiled_partial':
            view=arena.view('replayable_rank_arrays')
            for side,ranks in enumerate((left_ranks,right_ranks)):
                for index,rank in enumerate(ranks):
                    ledger.charge('rank_array_byte_writes',rankwidth);pos=(side*L+index)*rankwidth
                    view[pos:pos+rankwidth]=rank.to_bytes(rankwidth,'little')
        if mode=='resident':
            static.enter_context(arena.lease('base_left'));static.enter_context(arena.lease('base_right'))
            left=Records(arena,'base_left',L,h+ell,f,mem['base_record_bytes'],ledger)
            right=Records(arena,'base_right',L,h+ell,f,mem['base_record_bytes'],ledger)
            for i in range(L):left.write(i,record(0,i));right.write(i,record(1,i))
            sort_records(left,L,keypart,'base_sort_scratch',ledger,stats,'base_sort_comparisons','base_sort_record_moves')
        def merge(target,kind):
            ledger.charge('join_pass_starts');stats['join_pass_starts']+=1
            if kind=='L2':stats['L2_regeneration_passes']+=1
            starts=(0,) if mode=='resident' else range(0,L,tile)
            with ExitStack() as right_stack:
                if mode!='resident':
                    right_stack.enter_context(arena.lease(kind+'_right_record'))
                    right_record=Records(arena,kind+'_right_record',1,h+ell,f,mem['base_record_bytes'],ledger)
                for start in starts:
                    with ExitStack() as tile_stack:
                        ledger.charge('base_tile_starts');stats['base_tile_starts']+=1
                        used=L if mode=='resident' else min(tile,L-start)
                        if mode=='resident':first=left
                        else:
                            tile_stack.enter_context(arena.lease(kind+'_tile'))
                            first=Records(arena,kind+'_tile',tile,h+ell,f,mem['base_record_bytes'],ledger)
                            for i in range(used):first.write(i,record(0,start+i))
                            sort_records(first,used,keypart,kind+'_sort_scratch',ledger,stats,'base_sort_comparisons','base_sort_record_moves')
                        for index in range(L):
                            if mode=='resident':rr=right.read(index)
                            else:right_record.write(0,record(1,index));rr=right_record.read(0)
                            ledger.charge('first_join_queries');stats['first_join_queries']+=1
                            ledger.charge('field_additions_subtractions',u)
                            wanted=tuple((x-y)%p for x,y in zip(target,rr[h:h+u]))
                            for pos in equal_range(first,used,wanted,keypart,ledger,stats,'first_join_key_comparisons'):
                                rec=first.read(pos);ledger.charge('first_join_records');stats[kind+'_records']+=1
                                ledger.charge('field_additions_subtractions',ell-u)
                                ledger.charge('vector_coordinate_writes',(ell-u) if kind=='L1' else N+2*(ell-u))
                                yield rec[:h]+rr[:h],tuple((x+y)%p for x,y in zip(rec[h+u:],rr[h+u:]))
        stream=merge(offset,'L1');streams.append(stream)
        second=tuple((x-y)%p for x,y in zip(prepared.target[:u],offset))
        ledger.charge('field_additions_subtractions',u)
        try:
            while True:
                with arena.lease('root_block'):
                    root=Records(arena,'root_block',block,N+ell-u,f,mem['mid_record_bytes'],ledger);used=0
                    while used<block:
                        try:word,key=next(stream)
                        except StopIteration:break
                        root.write(used,word+key);used+=1
                    if not used:break
                    stats['root_blocks']+=1
                    sort_records(root,used,rootkey,'root_sort_scratch',ledger,stats,'root_block_sort_comparisons','root_sort_record_moves')
                    other=merge(second,'L2');streams.append(other)
                    try:
                        for word,key in other:
                            stats['root_join_queries']+=1;ledger.charge('field_additions_subtractions',ell-u)
                            wanted=tuple((x-y)%p for x,y in zip(prepared.target[u:],key))
                            for index in equal_range(root,used,wanted,rootkey,ledger,stats,'root_key_comparisons'):
                                ledger.charge('root_matches');stats['raw_root_matches']+=1
                                ledger.charge('field_additions_subtractions',N);ledger.charge('vector_coordinate_writes',N)
                                inside=tuple((x+y)%p for x,y in zip(root.read(index)[:N],word))
                                ledger.charge('Lee_integer_additions',N);ledger.charge('field_coordinate_comparisons',N+1)
                                if sum(min(x,p-x) for x in inside)!=prepared.v:continue
                                j=len(prepared.J);m=len(prepared.instance.H)
                                ledger.charge('field_multiplications',j*N);ledger.charge('field_additions_subtractions',j*(N+1))
                                ledger.charge('vector_coordinate_writes',n);candidate=prepared.complete(inside);stats['completion_attempts']+=1
                                ledger.charge('histogram_counter_updates',p);hist=[0]*p
                                for x in candidate:ledger.charge('histogram_counter_updates');hist[x]+=1
                                ledger.charge('histogram_counter_comparisons',p)
                                if tuple(hist)!=prepared.instance.C:continue
                                ledger.charge('field_multiplications',m*n);ledger.charge('field_additions_subtractions',m*n)
                                ledger.charge('field_coordinate_comparisons',m)
                                syndrome=tuple(sum(x*y for x,y in zip(row,candidate))%p for row in prepared.instance.H)
                                if syndrome==prepared.instance.s:
                                    outputs[candidate]+=1
                                    if not all_solutions:return answer()
                    finally:other.close();streams.remove(other)
        finally:
            for generator in streams:generator.close()
    return answer()
