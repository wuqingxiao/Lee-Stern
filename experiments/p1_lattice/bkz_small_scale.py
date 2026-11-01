"""Small-scale primal uSVP experiment for FCSD with the frozen composition shape.

For each k: n_code = 2k, composition = frozen C scaled to n_code (largest remainder),
random systematic H = [I | A] over F_251, secret e with that exact composition,
s = e1 + A e2. Kannan embedding of the q-ary lattice, progressive BKZ 2.0 (fpylll).
Record the smallest beta at which (e, 1) appears as the first basis vector,
and compare with the 2016 prediction sqrt(beta)*sigma <= delta^(2beta-d-1) q^(m/d).
"""
import json, math, random, sys, time
from fpylll import IntegerMatrix, LLL, BKZ, GSO
from fpylll.algorithms.bkz2 import BKZReduction

P = 251
fz = json.load(open("instances.json"))[0]
base = {int(v): c for v, c in fz["counts_centered"].items()}

def scaled(n):
    tot = sum(base.values()); raw = {v: c * n / tot for v, c in base.items()}
    out = {v: int(x) for v, x in raw.items()}
    rem = n - sum(out.values())
    for v in sorted(raw, key=lambda v: raw[v] - out[v], reverse=True)[:rem]: out[v] += 1
    return out

def log_delta(b):
    return math.log(((b / (2 * math.pi * math.e)) * (math.pi * b) ** (1 / b)) ** (1 / (2 * (b - 1))))

def predict(k, sigma):
    m = k; best = None
    for mm in range(10, m + 1):
        d = k + mm + 1
        for b in range(40, d):
            if math.log(math.sqrt(b) * sigma) <= (2 * b - d - 1) * log_delta(b) + mm / d * math.log(P):
                best = b if best is None else min(best, b); break
    return best

def trial(k, seed):
    rnd = random.Random(seed); n = 2 * k
    comp = scaled(n); e = [v for v, c in comp.items() for _ in range(c)]; rnd.shuffle(e)
    e1, e2 = e[:k], e[k:]
    A = [[rnd.randrange(P) for _ in range(k)] for _ in range(k)]
    s = [(e1[i] + sum(A[i][j] * e2[j] for j in range(k))) % P for i in range(k)]
    # lattice of (x1, x2, t) with x1 + A x2 = t s (mod P); basis rows:
    d = 2 * k + 1; B = IntegerMatrix(d, d)
    for i in range(k): B[i, i] = P
    for j in range(k):
        B[k + j, k + j] = 1
        for i in range(k): B[k + j, i] = (-A[i][j]) % P
    for i in range(k): B[2 * k, i] = s[i]
    B[2 * k, 2 * k] = 1
    target = e1 + e2
    sigma = math.sqrt(sum(x * x for x in e) / n)
    LLL.reduction(B)
    M = GSO.Mat(B); M.update_gso(); bkz = BKZReduction(M)
    for b in range(10, d, 2):
        bkz(BKZ.EasyParam(b, max_loops=8))
        v = [B[0, i] for i in range(d)]
        if abs(v[-1]) == 1 and [x * v[-1] for x in v[:-1]] == target:
            return b, sigma
    return None, sigma

if __name__ == "__main__":
    ks = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else [30, 40, 50, 60]
    reps = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    rows = []
    for k in ks:
        got = []
        for r in range(reps):
            t0 = time.time(); b, sg = trial(k, 1000 * k + r); got.append(b)
            print(f"k={k} rep={r} beta_found={b} sigma={sg:.3f} t={time.time()-t0:.1f}s", flush=True)
        pred = predict(k, sg)
        ok = [x for x in got if x]
        rows.append(dict(k=k, d=2 * k + 1, beta_pred=pred, beta_found=got,
                         beta_median=sorted(ok)[len(ok) // 2] if ok else None))
        print(rows[-1], flush=True)
    json.dump(rows, open("bkz_small_scale.json", "w"), indent=1)

