import json, sys
sys.path.insert(0, ".")
from run_estimator import FixedComposition, LWE
inst = json.load(open("instances.json"))[0]
X = FixedComposition(inst["counts_centered"])
P = LWE.Parameters(n=512, q=251, Xs=X, Xe=X, m=512, tag=inst["name"])
for name, f in [("dual", LWE.dual), ("dual_hybrid", LWE.dual_hybrid)]:
    try:
        print(name, f(P), flush=True)
    except Exception as ex:
        print(name, "FAILED", repr(ex), flush=True)