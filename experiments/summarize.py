import os
from bntaxonomy.experiment import ExperimentHandler
from automate_test import insts_ce, insts_practical, insts_tonello, get_insts
import pprint

from bntaxonomy.hierarchy import MultiInputSummary

cwd = os.path.dirname(os.path.abspath(__file__))

inst_groups = [insts_ce, insts_tonello, insts_practical]
hc = MultiInputSummary.from_folders(
    [
        f"{inst_selected.replace('instances','results')}/{inst}"
        for inst_selected, inst in get_insts(inst_groups)
    ]
)
print(hc.get_exp_names())
print(hc.to_conflict_matrix())
for a1, a2 in hc.ce_G.edges:
    print(f"{a1} -> {a2}:", [inst.name for inst in hc.ce_G.get_edge_data(a1, a2)["counterexamples"]])
hc.save(f"{cwd}/results/_summary")
