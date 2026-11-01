import json, sys, traceback
sys.path.insert(0, "/home/gql/lattice-estimator")
from sage.all import log, RR
from estimator import LWE, ND
for name, n, k, sd in [("n1024", 1024, 512, 2.8441), ("wide680", 640, 320, 3.6433)]:
    P = LWE.Parameters(n=k, q=251, Xs=ND.DiscreteGaussian(sd), Xe=ND.DiscreteGaussian(sd), m=n - k, tag=name)
    try:
        b = LWE.coded_bkw(P); print(name, b)
    except Exception:
        traceback.print_exc(limit=6)