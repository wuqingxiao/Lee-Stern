"""Per-commitment OS randomness buffer with exact rejection sampling.

No deterministic PRNG, no modulo reduction of biased samples, no buffer reuse
between commitments. Not constant-time; this is a research implementation.
"""
import secrets


class BufferedRandom:
    def __init__(self,block_bytes=1024,source=None):
        if type(block_bytes)!=int or block_bytes<1: raise ValueError('positive buffer size required')
        self.block_bytes=block_bytes
        self.source=secrets.token_bytes if source is None else source
        self.buffer=b'';self.position=0;self.closed=False

    def take(self,count):
        if self.closed: raise ValueError('random buffer already discarded')
        if type(count)!=int or count<0: raise ValueError('invalid byte count')
        parts=[]
        while count:
            if self.position==len(self.buffer):
                self.buffer=self.source(self.block_bytes);self.position=0
                if type(self.buffer)!=bytes or len(self.buffer)!=self.block_bytes:
                    raise ValueError('random source returned an invalid block')
            size=min(count,len(self.buffer)-self.position)
            parts.append(self.buffer[self.position:self.position+size]);self.position+=size;count-=size
        return b''.join(parts)

    def randrange(self,stop):
        if type(stop)!=int or stop<1: raise ValueError('positive integer bound required')
        if self.closed: raise ValueError('random buffer already discarded')
        if stop==1: return 0
        bits=(stop-1).bit_length();size=(bits+7)//8;mask=(1<<bits)-1
        while True:
            # Common field/permutation bounds fit one or two bytes. Avoid a
            # list, join and int.from_bytes allocation for every coordinate.
            if size<=2:
                while len(self.buffer)-self.position<size:
                    block=self.source(self.block_bytes)
                    if type(block)!=bytes or len(block)!=self.block_bytes:
                        raise ValueError('random source returned an invalid block')
                    self.buffer=self.buffer[self.position:]+block;self.position=0
                i=self.position;self.position+=size
                raw=self.buffer[i] if size==1 else (self.buffer[i]<<8)|self.buffer[i+1]
                candidate=raw&mask
            else:
                candidate=int.from_bytes(self.take(size),'big')&mask
            if candidate<stop: return candidate

    def shuffle(self,values):
        for i in range(len(values)-1,0,-1):
            j=self.randrange(i+1);values[i],values[j]=values[j],values[i]

    def discard(self):
        self.buffer=b'';self.position=0;self.closed=True

