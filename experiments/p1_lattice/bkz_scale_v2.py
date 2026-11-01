"""Small-scale primal uSVP experiment v2 (per-tour success check, parallel trials).

Instance (identical to bkz_small_scale.py): n_code = 2k, composition = frozen C
scaled to n_code by largest remainder, random systematic H = [I | A] over F_251,
s = e1 + A e2. Kannan embedding of dimension d = 2k + 1, target (e, 1).
Reduction: LLL, then progressive BKZ 2.0 (fplll C++, default pruning strategies)
with block sizes 10, 12, 14, ... and TOURS tours per block size. After every tour
all basis rows are checked for +-(e, 1). Records block size and tour of first
success, and saves the LLL profile (log2 GS norms) for simulator validation.
Usage: python bkz_scale_v2.py K1,K2,... REPS [TOURS]   (WORKERS env, default 6)
"""
import json, math, os, random, sys, time
from multiprocessing import Pool
from fpylll import BKZ, GSO, IntegerMatrix, LLL

P = 251
HERE = os.path.dirname(os.path.abspath(__file__))
BASE = {int(v): c for v, c in json.load(open(os.path.join(HERE, "instances.json")))[0]["counts_centered"].items()}

def scaled(n):
    tot = sum(BASE.values()); raw = {v: c * n / tot for v, c in BASE.items()}
    out = {v: int(x) for v, x in raw.items()}
    rem = n - sum(out.values())
    for v in sorted(raw, key=lambda v: raw[v] - out[v], reverse=True)[:rem]:
        out[v] += 1
    return out

def build(k, seed):
    rnd = random.Random(seed); n = 2 * k
    comp = scaled(n); e = [v for v, c in comp.items() for _ in range(c)]; rnd.shuffle(e)
    e1, e2 = e[:k], e[k:]
    A = [[rnd.randrange(P) for _ in range(k)] for _ in range(k)]
    s = [(e1[i] + sum(A[i][j] * e2[j] for j in range(k))) % P for i in range(k)]
    d = 2 * k + 1; B = IntegerMatrix(d, d)
    for i in range(k):
        B[i, i] = P
    for j in range(k):
        B[k + j, k + j] = 1
        for i in range(k):
            B[k + j, i] = (-A[i][j]) % P
    for i in range(k):
        B[2 * k, i] = s[i]
    B[2 * k, 2 * k] = 1
    return B, e1 + e2 + [1], math.sqrt(sum(x * x for x in e) / n)

def found(B, target):
    neg = [-x for x in target]
    for i in range(B.nrows):
        row = list(B[i])
        if row == target or row == neg:
            return i
    return None

def profile(B):
    M = GSO.Mat(B); M.update_gso()
    return [math.log2(M.get_r(i, i)) / 2 for i in range(B.nrows)]

def trial(args):
    k, seed, tours = args
    B, target, sigma = build(k, seed)
    t0 = time.time()
    B = LLL.reduction(B)
    res = dict(k=k, d=2 * k + 1, seed=seed, sigma=sigma, norm2=sum(x * x for x in target), lll_profile=profile(B))
    if found(B, target) is not None:
        res.update(beta=2, tour=0, time=time.time() - t0); return res
    for beta in range(10, 2 * k + 2, 2):
        par = BKZ.Param(block_size=beta, strategies=BKZ.DEFAULT_STRATEGY, max_loops=1,
                        flags=BKZ.MAX_LOOPS | BKZ.GH_BND)
        for t in range(1, tours + 1):
            B = BKZ.reduction(B, par)
            with open(os.path.join(HERE, "bkz_progress.log"), "a") as pf:
                pf.write("%s k=%d seed=%d beta=%d tour=%d t=%.0f\n" % (time.strftime("%H:%M:%S"), k, seed, beta, t, time.time() - t0))
            pos = found(B, target)
            if pos is not None:
                res.update(beta=beta, tour=t, row=pos, time=time.time() - t0); return res
    res.update(beta=None, time=time.time() - t0); return res

if __name__ == "__main__":
    ks = [int(x) for x in sys.argv[1].split(",")]
    reps = int(sys.argv[2]); tours = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    workers = int(os.environ.get("WORKERS", "6"))
    with open(os.path.join(HERE, "bkz_scale_v2.jsonl"), "a") as out, Pool(workers) as pool:
        for k in ks:
            tasks = [(k, 1000 * k + r, tours) for r in range(reps)]
            for res in pool.imap_unordered(trial, tasks):
                prof = res.pop("lll_profile")
                print(json.dumps(res), flush=True)
                res["lll_profile"] = prof; res["tours_per_beta"] = tours
                out.write(json.dumps(res) + "\n"); out.flush()