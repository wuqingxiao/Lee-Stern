"""Bit-operation (RC.MATZOV) and core-SVP (RC.ADPS16) cost of BKZ-beta at the 2016 block size
and at the probabilistic-model block sizes (prob_calib.json), same lattice dimension d."""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, "/home/gql/lattice-estimator")
from estimator import RC
from sage.all import log
cal = json.load(open(os.path.join(HERE, "prob_calib.json")))
rows = []
for s in cal["scale"]:
    d = s["d"]; out = dict(name=s["name"], d=d)
    for tag, b in [("2016", s["beta_2016"]), ("t2_chi2_med", s["t2_chi2"]["beta_med"]),
                   ("t8_chi2_med", s["t8_chi2"]["beta_med"]), ("t8_chi2_q01", s["t8_chi2"]["beta_q01"])]:
        out[tag] = dict(beta=b, matzov=round(float(log(RC.MATZOV(b, d), 2)), 1),
                        adps16=round(float(log(RC.ADPS16(b, d), 2)), 1))
    print(json.dumps(out), flush=True); rows.append(out)
json.dump(rows, open(os.path.join(HERE, "matzov_shift.json"), "w"), indent=1)