import os
from automate_test import insts_ce, insts_practical, insts_tonello, get_insts

from bntaxonomy.hierarchy import MultiInputSummary

cwd = os.path.dirname(os.path.abspath(__file__))

inst_groups = [result_dir.replace('instances','results') for result_dir in [insts_ce, insts_tonello, insts_practical]]
hc = MultiInputSummary.from_inst_groups(inst_groups, "Hierarchy")
for a1, a2 in hc.ce_G.edges:
    print(
        f"{a1} -> {a2}:",
        [inst.name for inst in hc.ce_G.get_edge_data(a1, a2)["counterexamples"]],
    )
hc.save(f"{cwd}/results/_summary")
print(hc.get_exp_names())
print(hc.to_conflict_matrix(full=True))
