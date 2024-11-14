
from bntaxonomy.experiment import ExperimentRun
import networkx as nx
import subprocess
from itertools import product

class Hierarchy:
    def __init__(self, exp_list: list[ExperimentRun], name: str = "Hierarchy"):
        self.name = name
        self.exp_list = exp_list
        # self.G = graphviz.Digraph(self.name)
        self.G, self.ce_G = nx.DiGraph(), nx.DiGraph()
        nodes = set(node for exp in self.exp_list for node in exp.G.nodes)
        # nx.nx_pydot.from_pydot
        for v1, v2 in product(nodes, nodes):
            if v1 == v2:
                continue
            # missing_inst = [inst for inst, G in zip(inst_list, G_list) if not(G.has_node(v1) or G.has_node(v2))]
            # if missing_inst:
            #     print(f"Missing nodes: ({v1} -> {v2}) in {missing_inst}")
            # if all(
            #     nx.has_path(G, v1, v2)
            #     for exp in self.exp_list
            #     if (G:= exp.G).has_node(v1) and G.has_node(v2)
            # ):
            #     self.G.add_edge(v1, v2)
            ce_exp_list = [
                exp
                for exp in self.exp_list
                if (G := exp.G).has_node(v1)
                and G.has_node(v2)
                and (not nx.has_path(G, v1, v2))
            ]
            if ce_exp_list:
                self.ce_G.add_edge(v1, v2, counterexamples=ce_exp_list)
            else:
                self.G.add_edge(v1, v2)

    @staticmethod
    def from_folders(folders: list[str], name: str = "Hierarchy"):
        exp_list = [
            ExperimentRun.from_folder(folder, folder.split("/")[-1])
            for folder in folders
        ]
        return Hierarchy(exp_list, name)

    def save(self, fname: str):
        nx.nx_pydot.write_dot(self.G, f"{fname}.dot")
        tred_cmd = f"tred {fname}.dot | dot -T png > {fname}.png"
        process = subprocess.Popen(tred_cmd, shell=True)
        process.wait()
