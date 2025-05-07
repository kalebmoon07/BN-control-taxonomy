import os
from colomoto.minibn import BooleanNetwork
import mpbn
import networkx as nx
from bntaxonomy.iface.pbn import make_pbn_primes_iface
from pyboolnet.attractors import (
    compute_attractors,
    compute_attractors_tarjan,
    compute_trap_spaces,
)
from bntaxonomy.utils.control import suppress_console_output

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
        print(f"\tRunning {inst}")
        fpath = f"{inst_selected}/{inst}"
        bnet_fname = f"{fpath}/transition_formula.bnet"
        org_bnet = BooleanNetwork(data=bnet_fname)
        primes = make_pbn_primes_iface(bnet_fname)

        f = mpbn.MPBooleanNetwork(org_bnet)
        for update in ["synchronous", "asynchronous"]:
            print(f"Update: {update}")
            g = f.dynamics(update)
            nx.nx_pydot.write_dot(g, f"{fpath}/STG_{update}.dot")
            with suppress_console_output():
                states, attrs = compute_attractors_tarjan(g)
            print("Fixed points: ", states)
            print("Cyclic attractors: ", attrs)

        print("Minimal trap spaces: ", compute_trap_spaces(primes))
