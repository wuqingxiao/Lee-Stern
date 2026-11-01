import json, math
rows = json.load(open("arora_ge_dreg.json"))
def lb(a, b): return (math.lgamma(a + 1) - math.lgamma(b + 1) - math.lgamma(a - b + 1)) / math.log(2)
for r in rows:
    D, nv, s = r["dreg"], r["nv"], r["support"]
    lm = lb(nv + D, D)
    red = nv * math.log2(s)          # reduced monomials (each exponent < |supp|) are at most |supp|^nv
    r["log2_monomials"] = round(min(lm, red), 1)
    r["log2_cost_omega2"] = round(2 * r["log2_monomials"], 1)
    r["log2_cost_omega281"] = round(2.81 * r["log2_monomials"], 1)
    print(r["name"], "ps" if r["power_sums"] else "  ", "nv", nv, "dreg", D, "log2 mono", r["log2_monomials"], "cost w=2", r["log2_cost_omega2"])
json.dump(rows, open("arora_ge_dreg.json", "w"), indent=1)
