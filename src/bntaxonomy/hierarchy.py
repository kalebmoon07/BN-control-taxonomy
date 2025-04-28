import json
import os
import networkx as nx
from itertools import combinations, product


from bntaxonomy.utils.control import CtrlResult, suppress_console_output
import bntaxonomy.utils.graph as graph_utils


class SingleInputSummary:
    def __init__(self, results: list[CtrlResult], name: str = "SingleInputSummary"):
        self.name = name
        self.results = results
        self.G = nx.DiGraph()
        for r1, r2 in combinations(self.results, 2):
            if r1.is_stronger_than(r2):
                self.G.add_edge(r2.name, r1.name)
            if r2.is_stronger_than(r1):
                self.G.add_edge(r1.name, r2.name)

    def save(self, fname: str):
        with suppress_console_output():
            dot_fname = f"{fname}.dot"
            tred_fname = f"{fname}_tred.dot"
            tred_scc_fname = f"{fname}_tred_scc.dot"
            graph_utils.write_dot(self.G, dot_fname)
            graph_utils.write_transitive_reduction(dot_fname, tred_fname)
            graph_utils.export_dot_png(tred_fname, f"{fname}.png")
            graph_utils.cluster_cycles(tred_fname, f"{fname}_tred_scc.dot")
            graph_utils.export_dot_png(f"{fname}_tred_scc.dot", f"{fname}_tred_scc.png")

    @staticmethod
    def from_folder(opath: str, name: str = ""):
        files = [fname for fname in os.listdir(opath) if fname.endswith(".json")]
        sol_list = [
            CtrlResult(fname[:-5], json.load(open(f"{opath}/{fname}")))
            for fname in files
        ]
        if not name:
            name = opath.split("/")[-1]
        return SingleInputSummary(sol_list, name)


class MultiInputSummary:
    def __init__(
        self, exp_list: list[SingleInputSummary], name: str = "MultiInputSummary"
    ):
        self.name = name
        self.exp_list = exp_list
        self.G, self.ce_G = nx.DiGraph(), nx.DiGraph()
        nodes = set(node for exp in self.exp_list for node in exp.G.nodes)
        self.G.add_nodes_from(nodes)
        self.ce_G.add_nodes_from(nodes)
        for v1, v2 in product(nodes, nodes):
            if v1 == v2:
                continue
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
            SingleInputSummary.from_folder(folder, folder.split("/")[-1])
            for folder in folders
        ]
        return MultiInputSummary(exp_list, name)

    def save(self, fname: str):
        with suppress_console_output():
            dot_fname = f"{fname}.dot"
            tred_fname = f"{fname}_tred.dot"
            graph_utils.write_dot(self.G, dot_fname)
            graph_utils.write_transitive_reduction(dot_fname, tred_fname)
            graph_utils.export_dot_png(tred_fname, f"{fname}.png")
            graph_utils.cluster_cycles(tred_fname, f"{fname}_tred_scc.dot")
            graph_utils.export_dot_png(f"{fname}_tred_scc.dot", f"{fname}_tred_scc.png")
