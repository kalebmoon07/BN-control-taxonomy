import os
from bntaxonomy.dep.converters import bn_of_asynchronous_transition_graph
from colomoto.minibn import BooleanNetwork
import mpbn
import networkx as nx
import pydot
from bntaxonomy.experiment import ExperimentHandler


cwd = os.path.dirname(os.path.abspath(__file__))

insts_practical = f"{cwd}/instances/practical"
insts_ce = f"{cwd}/instances/counterexamples"


def get_insts(inst_groups: list[str]):
    for inst_selected in inst_groups:
        inst_list = sorted(os.listdir(inst_selected))
        for inst in inst_list:
            yield inst_selected, inst


if __name__ == "__main__":
    max_size = 2
    # inst_groups = [insts_ce, insts_practical]
    inst_groups = [insts_ce]

    for inst_selected, inst in get_insts(inst_groups):
        fpath = f"{inst_selected}/{inst}"
        bnet_fname = f"{fpath}/transition_formula.bnet"
        org_bnet = BooleanNetwork(data=bnet_fname)

        f = mpbn.MPBooleanNetwork(org_bnet)
        g = f.dynamics("synchronous")
        nx.nx_pydot.write_dot(g, f"{fpath}/STG_synchronous.dot")
        g = f.dynamics("asynchronous")
        nx.nx_pydot.write_dot(g, f"{fpath}/STG_asynchronous.dot")
