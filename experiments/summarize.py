import os
from bntaxonomy.experiment import ExperimentHandler
from automate_test import (
    insts_ce,
    insts_practical,
    insts_tonello,
    insts_all,
)
import pprint

from bntaxonomy.hierarchy import Hierarchy

cwd = os.path.dirname(os.path.abspath(__file__))

inst_list = insts_all
hc = Hierarchy.from_folders([f"{cwd}/results/{inst}" for inst in inst_list])
# for v1, v2 in [("SM-bf","Caspo")]:
#     pprint.pprint(((v1,v2), [exp.name for exp in hc.ce_G.get_edge_data(v1,v2)["counterexamples"]]))
# hc = Hierarchy.from_folders([f"results_/{inst}" for inst in inst_list])
# for v1, v2 in [("Caspo","SM-brute-force")]:
#     pprint.pprint(((v1,v2), [exp.name for exp in hc.ce_G.get_edge_data(v1,v2)["counterexamples"]]))
hc.save(f"{cwd}/results/_summary")
