"""lattice-estimator primal_usvp prediction for the small experimental instances
(same composition scaling as bkz_scale_v2.py), core-SVP cost model ADPS16,
GSA and CN11 basis-shape models. Same model family as the n=640/1024 numbers."""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, "/home/gql/lattice-estimator"); sys.path.insert(0, HERE)
from run_estimator import FixedComposition, LWE, RC
import bkz_scale_v2 as bk
rows = []
for k in [40, 45, 50, 55, 60, 65, 70]:
    comp = {str(v): c for v, c in bk.scaled(2 * k).items() if c}
    X = FixedComposition(comp)
    P = LWE.Parameters(n=k, q=251, Xs=X, Xe=X, m=k)
    row = dict(k=k, d_full=2 * k + 1, sigma=float(X.stddev))
    for shape in ["gsa", "cn11"]:
        r = LWE.primal_usvp(P, red_cost_model=RC.ADPS16, red_shape_model=shape)
        row["beta_" + shape] = int(r["beta"]); row["d_" + shape] = int(r["d"])
    print(json.dumps(row), flush=True); rows.append(row)
json.dump(rows, open(os.path.join(HERE, "est_small.json"), "w"), indent=1)