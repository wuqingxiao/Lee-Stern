"""P1: lattice-estimator on FCSD instances.

He = s with systematic H = [I_m | A] gives e1 + A e2 = s (mod p): LWE with
n = k secret coords (e2), m = n_code - k samples, error e1. Both parts are
coordinates of a vector with fixed full composition C, so per-coordinate
marginals are the empirical distribution of C (stddev, bounds, density exact).
support_size uses |S_C| restricted by length (multinomial of the marginal),
which is what exhaustive guessing must cover.
"""
import json, sys
from sage.all import RR, sqrt, factorial, log, ceil, oo
sys.path.insert(0, "/home/gql/lattice-estimator")
from estimator import LWE, ND, RC
from estimator.nd import NoiseDistribution

class FixedComposition(NoiseDistribution):
    def __init__(self, counts, n=None):
        self.counts = {int(v): int(c) for v, c in counts.items()}
        tot = sum(self.counts.values())
        mean = RR(sum(v * c for v, c in self.counts.items()) / tot)
        var = RR(sum(c * (v - mean) ** 2 for v, c in self.counts.items()) / tot)
        super().__init__(n=n, mean=mean, stddev=sqrt(var),
                         bounds=(min(self.counts), max(self.counts)),
                         _density=RR(1 - self.counts.get(0, 0) / tot))
        self._tot = tot
    def resize(self, new_n):
        return FixedComposition(self.counts, new_n)
    def __hash__(self):
        return hash(("FixedComposition", tuple(sorted(self.counts.items())), self.n))
    def support_size(self, fraction=1.0):
        n = len(self)
        # composition of an n-subvector scaled from C; multinomial count
        scaled = {v: RR(c * n / self._tot) for v, c in self.counts.items()}
        lg = log(factorial(n)) - sum(log(RR(x).gamma() * 1) if False else RR(x + 1).log_gamma() for x in scaled.values())
        return ceil(RR(fraction) * lg.exp())

def run(inst, mode):
    p, n, k = inst["p"], inst["n"], inst["k"]
    X = FixedComposition(inst["counts_centered"])
    P = LWE.Parameters(n=k, q=p, Xs=X, Xe=X, m=n - k, tag=inst["name"])
    print("=" * 70); print(P); print("stddev %.4f density %.4f" % (float(X.stddev), float(X._density)))
    if mode == "rough":
        return LWE.estimate.rough(P, jobs=8)
    return LWE.estimate(P, jobs=8, deny_list=("arora-gb", "bkw"))

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "rough"
    insts = json.load(open(sys.argv[2] if len(sys.argv) > 2 else "instances.json"))
    from estimator import schemes
    print("anchor Kyber512 rough:"); LWE.estimate.rough(schemes.Kyber512)
    for inst in insts:
        try:
            r = run(inst, mode)
        except Exception as ex:
            print("FAILED", inst["name"], repr(ex))
