"""Roll the C-FIREWALL result files into one summary JSON."""
import json, glob, os, sys
R = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "/results"

def load(f): return json.load(open(os.path.join(R, f)))

groups = {
 "census_k1":   ["cfirewall_p3.json", "cfirewall_p4.json", "cfirewall_p5.json"],
 "sample_k1":   ["cfirewall_p6_sample.json", "cfirewall_p7_sample.json",
                 "cfirewall_p8_sample.json"],
 "kge2":        ["kfirewall_p4_k2.json", "kfirewall_p4_k3.json", "kfirewall_p5_k2.json",
                 "kfirewall_p6_k2.json", "kfirewall_p6_k3.json", "kfirewall_p7_k4.json"],
}
FIELDS = ["trial","abort","abort_nopath","abort_unamen","report","report_moved",
          "report_multiD","report_multiD_moved","checks","clean_trials",
          "violating_trials","violating_keys","violating_keys_allD","violating_keys_someD",
          "ctrl_trial","ctrl_report","ctrl_moved","ctrl_violating",
          "ctrl_violating_allD","ctrl_violating_someD","stmt_all","stmt_coherent",
          "stmt_rejected","mpdag","draw","draw_skipped_bigclass","ext_total"]

out = {"per_run": {}, "totals": {}}
for g, files in groups.items():
    tot = {k: 0 for k in FIELDS}
    for f in files:
        d = load(f); c = d["counts"]
        out["per_run"][f] = dict(p=d["p"], k=d.get("k", 1), mode=d["mode"],
                                 exhaustive=d.get("exhaustive"),
                                 elapsed_s=d["elapsed_s"],
                                 n_witnesses=d["n_witnesses"],
                                 counts={k: c.get(k, 0) for k in FIELDS if k in c})
        for k in FIELDS:
            tot[k] += c.get(k, 0)
    out["totals"][g] = tot

st = {}
for f in sorted(glob.glob(os.path.join(R, "structure_p*.json"))):
    d = json.load(open(f)); st[os.path.basename(f)] = dict(p=d["p"], mode=d.get("mode","census"),
                                                           counts=d["counts"])
out["structure"] = st

g = out["totals"]
allk1 = {k: g["census_k1"][k] + g["sample_k1"][k] for k in FIELDS}
out["grand_k1"] = allk1
out["grand_all_arms_trials"] = allk1["trial"] + g["kge2"]["trial"]
out["grand_all_arms_checks"] = allk1["checks"] + g["kge2"]["checks"]
out["violations_total"] = (allk1["violating_trials"] + g["kge2"]["violating_trials"])
n = out["grand_all_arms_trials"]
out["rule_of_three_upper_95"] = 3.0 / n
json.dump(out, open(os.path.join(R, "SUMMARY.json"), "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != "per_run" and k != "structure"}, indent=1))
