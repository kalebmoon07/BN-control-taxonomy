import json
import os
from automate_test import insts_ce, insts_practical, get_insts

from bntaxonomy.hierarchy import MultiInputSummary

cwd = os.path.dirname(os.path.abspath(__file__))

inst_groups = [
    result_dir.replace("instances", "results")
    for result_dir in [insts_practical, insts_ce]
]
hc = MultiInputSummary.from_inst_groups(inst_groups, "Hierarchy")
for a1, a2 in hc.ce_G.edges:
    print(
        f"{a1} -> {a2}:",
        [inst.name for inst in hc.ce_G.get_edge_data(a1, a2)["counterexamples"]],
    )
hc.save(f"{cwd}/results/_summary")

exp_names = hc.get_exp_names()
group_names = hc.get_exp_group_names()
with open(f"{cwd}/results/counterexample_group_list.json", "w") as f:
    json.dump(group_names, f, indent=4)

with open(f"{cwd}/results/counterexample_first_match.csv", "w") as f:
    f.write(hc.to_conflict_matrix_csv(full_ce=False))
with open(f"{cwd}/results/counterexample_full_match.csv", "w") as f:
    f.write(hc.to_conflict_matrix_csv(full_ce=True))