"""Exact batched fixed-bound sampling; same random-byte tape semantics."""
from batched_random import BufferedRandom

class FastRandom(BufferedRandom):
    def vector(self,stop,count):
        if type(stop)!=int or stop<1 or type(count)!=int or count<0:
            raise ValueError('invalid bound or count')
        if self.closed:raise ValueError('random buffer already discarded')
        if stop==1:return [0]*count
        bits=(stop-1).bit_length();size=(bits+7)//8;mask=(1<<bits)-1
        result=[]
        while len(result)<count:
            raw=self.take((count-len(result))*size)
            if size==1:
                result.extend(x for byte in raw if (x:=byte&mask)<stop)
            elif size==2:
                result.extend(x for i in range(0,len(raw),2)
                              if (x:=((raw[i]<<8)|raw[i+1])&mask)<stop)
            else:
                result.extend(x for i in range(0,len(raw),size)
                              if (x:=int.from_bytes(raw[i:i+size],'big')&mask)<stop)
        return result

    def shuffle(self,values):
        if self.closed:raise ValueError('random buffer already discarded')
        # Inline one/two-byte sampling and keep buffer state in local variables.
        buf=self.buffer;pos=self.position;source=self.source;blocksize=self.block_bytes
        try:
            for i in range(len(values)-1,0,-1):
                bits=i.bit_length();size=(bits+7)//8;mask=(1<<bits)-1
                while True:
                    while len(buf)-pos<size:
                        block=source(blocksize)
                        if type(block)!=bytes or len(block)!=blocksize:
                            raise ValueError('invalid random block')
                        buf=buf[pos:]+block;pos=0
                    if size==1:x=buf[pos]&mask
                    elif size==2:x=((buf[pos]<<8)|buf[pos+1])&mask
                    else:x=int.from_bytes(buf[pos:pos+size],'big')&mask
                    pos+=size
                    if x<=i:break
                values[i],values[x]=values[x],values[i]
        finally:self.buffer=buf;self.position=pos

class MaskRandom(FastRandom):
    shuffle=BufferedRandom.shuffle

