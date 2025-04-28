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

    def get_exp_names(self):
        return [exp.name for exp in self.exp_list]

    def to_conflict_matrix(self):
        # conflict_matrix = nx.to_numpy_array(self.ce_G, nodelist=self.ce_G.nodes())
        # return conflict_matrix
        column_str = """
                        \\newcolumntype{O}[2]{%
                            >{\\adjustbox{angle=#1,lap=\\width-(#2)}\\bgroup} l <{\\egroup}%
                        }
                        \\NewDocumentCommand{\\rot}{O{90} O{1em} m}{\\makebox[#2][l]{\\rotatebox{#1}{#3}}}%
                        \\setlength{\\cmidrulekern}{2.5pt}
                    """
        nodes = sorted(self.ce_G.nodes())  # fix consistent order
        n = len(nodes)
        matrix = [["" for _ in range(n)] for _ in range(n)]

        node_idx = {node: idx for idx, node in enumerate(nodes)}
        exp_names = self.get_exp_names()
        for src, dst, data in self.ce_G.edges(data=True):
            if "counterexamples" in data and data["counterexamples"]:
                example = data["counterexamples"][0].name  # assuming .name exists
                i = node_idx[src]
                j = node_idx[dst]
                matrix[i][j] = str(exp_names.index(example))

        # Build LaTeX tabular code
        latex = [column_str]
        latex.append("\\begin{tabular}{" + "|".join("c" * (n + 1)) + "}")
        latex.append("\\hline")

        # Header
        # header = [""] + list(nodes)
        header = [""] + [f"\\rot{{{node}}}" for node in nodes]
        latex.append(" & ".join(header) + " \\\\ \\hline")

        # Rows
        for i, row_node in enumerate(nodes):
            row = [row_node] + [matrix[i][j] for j in range(n)]
            latex.append(" & ".join(row) + " \\\\ \\hline")

        latex.append("\\end{tabular}")

        return "\n".join(latex)
