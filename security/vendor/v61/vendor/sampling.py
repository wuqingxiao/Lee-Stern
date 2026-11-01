"""Floyd subset sampler with an AVL set: worst-case logarithmic tree visits.

The Python nodes below check logic/counts. The memory specification reserves
fixed-width rank keys plus four uint64 fields per node, not Python object RAM.
"""
from collections import Counter


class Node:
    def __init__(self,key):self.key=key;self.left=None;self.right=None;self.height=1


def sample_ranks(size,count,rng,counters=None):
    if not 0<=count<=size:raise ValueError('sample size')
    c=counters if counters is not None else Counter()
    def height(node):return node.height if node else 0
    def refresh(node):node.height=1+max(height(node.left),height(node.right))
    def rotate_right(y):
        c['rotations']+=1;x=y.left;y.left=x.right;x.right=y;refresh(y);refresh(x);return x
    def rotate_left(x):
        c['rotations']+=1;y=x.right;x.right=y.left;y.left=x;refresh(x);refresh(y);return y
    def contains(node,key):
        while node:
            c['node_visits']+=1;c['key_comparisons']+=1
            if key==node.key:return True
            c['key_comparisons']+=1
            node=node.left if key<node.key else node.right
        return False
    def insert(node,key):
        if node is None:c['nodes']+=1;return Node(key)
        c['node_visits']+=1;c['key_comparisons']+=1
        if key<node.key:node.left=insert(node.left,key)
        else:node.right=insert(node.right,key)
        refresh(node);balance=height(node.left)-height(node.right)
        if balance>1:
            if height(node.left.left)<height(node.left.right):node.left=rotate_left(node.left)
            return rotate_right(node)
        if balance< -1:
            if height(node.right.right)<height(node.right.left):node.right=rotate_right(node.right)
            return rotate_left(node)
        return node
    root=None
    for j in range(size-count,size):
        c['uniform_integer_calls']+=1;candidate=rng.randrange(j+1)
        chosen=j if contains(root,candidate) else candidate
        root=insert(root,chosen)
        yield chosen
