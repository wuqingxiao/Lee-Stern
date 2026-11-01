"""Disjoint named arenas and real small packed record buffers; no RSS claim."""
from contextlib import contextmanager
from collections import Counter


class ResourceCutoff(RuntimeError):
    def __init__(self,name,used,request,limit):
        self.details=dict(resource=name,used=used,request=request,limit=limit)
        super().__init__('cutoff BEFORE operation: '+name)


class Ledger:
    def __init__(self,caps=None):
        self.caps=None if caps is None else {k:int(v) for k,v in caps.items()}
        if self.caps is not None and any(v<0 for v in self.caps.values()):raise ValueError('negative cap')
        self.counts=Counter()
    def charge(self,name,amount=1):
        if type(amount)!=int or amount<0:raise ValueError('nonnegative integer charge')
        if not amount:return
        used=self.counts[name]
        if self.caps is not None:
            if name not in self.caps:raise ValueError('unbudgeted resource: '+name)
            if used+amount>self.caps[name]:raise ResourceCutoff(name,used,amount,self.caps[name])
        self.counts[name]+=amount


class InstrumentedCounter(Counter):
    """Reserve defined kernel primitives before the operation, retain auxiliaries."""
    def __init__(self,ledger):super().__init__();self.ledger=ledger
    def __setitem__(self,key,value):
        delta=value-self.get(key,0)
        mapping={'rank_bigint_subtractions':('rank_bigint_additions_subtractions',1),
            'rank_bigint_divmod':('rank_bigint_divisions',1),
            'key_comparisons':('rank_bigint_comparisons',1),
            'node_visits':('rank_tree_node_visits_upper',1),
            'nodes':('rank_tree_node_visits_upper',1),'rotations':('rank_tree_node_visits_upper',8),
            'uniform_integer_calls':('sample_uniform_integer_calls',1)}
        if key in mapping:
            name,factor=mapping[key];self.ledger.charge(name,delta*factor)
        elif key in ('rank_bigint_multiplications','rank_bigint_divisions','rank_bigint_comparisons'):
            self.ledger.charge(key,delta)
        super().__setitem__(key,value)


class Arena:
    def __init__(self,layout,*,max_backing_bytes=16*2**20):
        self.layout=layout;self.regions=layout['regions'];self.active={};self.peak=0;self.peak_names=[];self.events=[]
        intervals=sorted((v['offset'],v['offset']+v['bytes'],k) for k,v in self.regions.items() if v['bytes'])
        assert all(a[1]<=b[0] for a,b in zip(intervals,intervals[1:]))
        assert max((x[1] for x in intervals),default=0)==layout['reserved_bytes']
        if layout['reserved_bytes']>max_backing_bytes:raise ValueError('small packed oracle allocation limit')
        self.backing=bytearray(layout['reserved_bytes'])
    def note(self,event):
        used=sum(self.active.values())
        if used>self.peak:self.peak=used;self.peak_names=sorted(self.active)
        assert used<=self.layout['reserved_bytes']
        if len(self.events)<128:self.events.append(dict(event=event,active=sorted(self.active),bytes=used))
    @contextmanager
    def lease(self,name):
        if name in self.active:raise RuntimeError('arena double lease: '+name)
        self.active[name]=self.regions[name]['bytes'];self.note('acquire:'+name)
        try:yield
        finally:self.active.pop(name);self.note('release:'+name)
    def view(self,name):
        if name not in self.active:raise RuntimeError('use outside lifetime: '+name)
        r=self.regions[name];return memoryview(self.backing)[r['offset']:r['offset']+r['bytes']]


class Records:
    def __init__(self,arena,name,capacity,coordinates,fieldbytes,stride,ledger):
        self.arena=arena;self.name=name;self.capacity=capacity;self.coordinates=coordinates
        self.f=fieldbytes;self.stride=stride;self.ledger=ledger
        if capacity*stride>arena.regions[name]['bytes'] or coordinates*fieldbytes>stride:
            raise ValueError('packed capacity')
    def read(self,index):
        if not 0<=index<self.capacity:raise IndexError(index)
        b=self.arena.view(self.name);start=index*self.stride
        return tuple(int.from_bytes(b[start+i*self.f:start+(i+1)*self.f],'little') for i in range(self.coordinates))
    def key(self,index,part):
        if not 0<=index<self.capacity:raise IndexError(index)
        b=self.arena.view(self.name);start=index*self.stride
        return tuple(int.from_bytes(b[start+i*self.f:start+(i+1)*self.f],'little') for i in range(part.start,part.stop))
    def write(self,index,values):
        if not 0<=index<self.capacity or len(values)!=self.coordinates:raise ValueError('record shape')
        self.ledger.charge('vector_coordinate_writes',self.coordinates)
        b=self.arena.view(self.name);start=index*self.stride
        for i,x in enumerate(values):b[start+i*self.f:start+(i+1)*self.f]=int(x).to_bytes(self.f,'little')
    def copy_to(self,source,destination,index,move_unit):
        self.ledger.charge(move_unit)
        self.ledger.charge('sort_payload_bytes_read_upper',self.stride)
        self.ledger.charge('sort_payload_bytes_written_upper',self.stride)
        src=self.arena.view(self.name);dst=destination.arena.view(destination.name)
        dst[index*self.stride:(index+1)*self.stride]=src[source*self.stride:(source+1)*self.stride]


def compare(a,b,ledger,*,sort_read_bytes=0):
    if sort_read_bytes:ledger.charge('sort_payload_bytes_read_upper',sort_read_bytes)
    for x,y in zip(a,b):
        ledger.charge('field_coordinate_comparisons')
        if x<y:return -1
        if x>y:return 1
    return 0


def sort_records(main,count,key_slice,scratch_name,ledger,stats,stat,move_unit):
    with main.arena.lease(scratch_name):
        scratch=Records(main.arena,scratch_name,main.capacity,main.coordinates,main.f,main.stride,ledger)
        a=main;b=scratch;width=1
        while width<count:
            for start in range(0,count,2*width):
                mid=min(start+width,count);end=min(start+2*width,count);i=start;j=mid
                for target in range(start,end):
                    take=j>=end
                    if i<mid and j<end:
                        stats[stat]+=1
                        ledger.charge('sort_payload_bytes_read_upper',2*main.f*(key_slice.stop-key_slice.start))
                        left=a.key(i,key_slice);right=a.key(j,key_slice)
                        take=compare(left,right,ledger)<=0
                    if take:a.copy_to(i,b,target,move_unit);i+=1
                    else:a.copy_to(j,b,target,move_unit);j+=1
            a,b=b,a;width*=2
        if a is not main:
            for i in range(count):a.copy_to(i,main,i,move_unit)


def equal_range(main,count,target,key_slice,ledger,stats,stat):
    lo=0;hi=count
    while lo<hi:
        mid=(lo+hi)//2;stats[stat]+=1
        if compare(main.key(mid,key_slice),target,ledger)<0:lo=mid+1
        else:hi=mid
    first=lo;hi=count
    while lo<hi:
        mid=(lo+hi)//2;stats[stat]+=1
        if compare(target,main.key(mid,key_slice),ledger)<0:hi=mid
        else:lo=mid+1
    return range(first,lo)
