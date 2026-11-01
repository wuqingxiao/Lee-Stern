import json, math
insts = json.load(open(r"F:\codex-pj\ac3\research\next_steps_128_20261101\p1_lattice\instances.json"))
def lmult(n, parts): return (math.lgamma(n + 1) - sum(math.lgamma(x + 1) for x in parts)) / math.log(2)
for inst in insts:
    n, k = inst["n"], inst["k"]; m = n - k
    cc = [c for c in inst["counts_centered"].values()]
    full = lmult(n, cc)
    half = lmult(m, [c * m / n for c in cc])
    print(f'{inst["name"]:36s} log2|S_C|={full:8.1f}  info-set enumeration ~ 2^{full - half:7.1f}')
