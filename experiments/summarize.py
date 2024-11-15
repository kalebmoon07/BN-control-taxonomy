import os
from bntaxonomy.experiment import ExperimentHandler
from automate_test import insts_ce, insts_practical, insts_tonello, get_insts
import pprint

from bntaxonomy.hierarchy import Hierarchy

cwd = os.path.dirname(os.path.abspath(__file__))

inst_groups = [insts_practical, insts_tonello, insts_ce]
hc = Hierarchy.from_folders(
    [
        f"{inst_selected.replace('instances','results')}/{inst}"
        for inst_selected, inst in get_insts(inst_groups)
    ]
)
# for v1, v2 in [("PBN-perc","Caspo")]:
#     pprint.pprint(((v1,v2), [exp.name for exp in hc.ce_G.get_edge_data(v1,v2)["counterexamples"]]))
# hc = Hierarchy.from_folders([f"results_/{inst}" for inst in inst_list])
# for v1, v2 in [("Caspo","SM-brute-force")]:
#     pprint.pprint(((v1,v2), [exp.name for exp in hc.ce_G.get_edge_data(v1,v2)["counterexamples"]]))
hc.save(f"{cwd}/results/_summary")
