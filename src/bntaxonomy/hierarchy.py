from collections import defaultdict
import json
import os
import networkx as nx
from itertools import combinations, product
from colomoto.minibn import BooleanNetwork


from bntaxonomy.utils.control import CtrlResult, suppress_console_output
import bntaxonomy.utils.graph as graph_utils


class SingleInputSummary:
    def __init__(self, results: list[CtrlResult], name: str = "SingleInputSummary", bn: BooleanNetwork | None = None):
        self.name = name
        self.results = results
        self.G = nx.DiGraph()
        self.bn = bn

        for r1, r2 in combinations(self.results, 2):
            if r1.is_stronger_than(r2):
                self.G.add_edge(r2.name, r1.name)
            if r2.is_stronger_than(r1):
                self.G.add_edge(r1.name, r2.name)

    def save(self, fname: str):
        with suppress_console_output():
            dot_fname = f"{fname}.dot"
            tred_fname = f"{fname}_tred.dot"
            graph_utils.write_dot(self.G, dot_fname)
            graph_utils.write_transitive_reduction(dot_fname, tred_fname)
            graph_utils.export_dot_png(tred_fname, f"{fname}.png")
            graph_utils.cluster_cycles(tred_fname, tred_fname)
            graph_utils.export_dot_png(tred_fname, f"{fname}_tred.png")

    @staticmethod
    def from_folder(opath: str, name: str = "", bn: BooleanNetwork | None = None):
        files = [
            fname
            for fname in os.listdir(opath)
            if fname.endswith(".json") and not fname.endswith("_full.json")
        ]
        sol_list = [
            CtrlResult(fname[:-5], json.load(open(f"{opath}/{fname}")))
            for fname in files
        ]
        if not name:
            name = opath.split("/")[-1]
        return SingleInputSummary(sol_list, name, bn)


class MultiInputSummary:
    def __init__(
        self,
        exp_list: list[SingleInputSummary],
        name: str = "MultiInputSummary",
        exp_groups: dict[str, list[SingleInputSummary]] = dict(),
    ):
        self.name = name
        self.exp_list = exp_list
        self.exp_groups = exp_groups
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

    @staticmethod
    def from_inst_groups(groups: list[str], name: str = "Hierarchy"):
        exp_groups = defaultdict(list)
        exp_list = []
        for inst_selected in groups:
            group_name = os.path.basename(inst_selected)
            inst_bn_folder = inst_selected.replace("results", "instances")
            for exp_name in os.listdir(inst_selected):
                bn = BooleanNetwork.load(f"{inst_bn_folder}/{exp_name}/transition_formula.bnet")
                input_summary = SingleInputSummary.from_folder(
                    f"{inst_selected}/{exp_name}", exp_name, bn
                )
                exp_groups[group_name].append(input_summary)
                exp_list.append(input_summary)
        return MultiInputSummary(exp_list, name, exp_groups)

    def save(self, fname: str):
        with suppress_console_output():
            dot_fname = f"{fname}.dot"
            tred_fname = f"{fname}_tred.dot"
            graph_utils.write_dot(self.G, dot_fname)
            graph_utils.write_transitive_reduction(dot_fname, tred_fname)
            graph_utils.export_dot_png(tred_fname, f"{fname}.png")
            graph_utils.cluster_cycles(tred_fname, tred_fname)
            graph_utils.export_dot_png(tred_fname, f"{fname}_tred.png")

    def get_exp_names(self):
        return [exp.name for exp in self.exp_list]
    
    def get_exp_group_names(self):
        return {group: [exp.name for exp in exp_list] for group, exp_list in self.exp_groups.items()}

    def to_conflict_matrix(self, use_group_idx=True, full_ce=False):
        """
            Generates a LaTeX table representing the conflict matrix of the counterexamples.
            The table is formatted for use in a LaTeX document and includes rotation for headers.
            The function can be customized to use group indices or full counterexample names.
            Prints the LaTeX code for the table.
            Args:
                use_group_idx (bool): If True, uses group indices for the counterexamples; otherwise, uses the index of the counterexample in the list.
                full_ce (bool): If True, includes full counterexample names in the matrix; otherwise, only the first matched name is used.
        """
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

        group_idx_dict = {
            exp.name: chr(65 + idx) + str(exp_list.index(exp) + 1) ## A1, A2, B1, B2, etc.
            for idx, exp_list in enumerate(self.exp_groups.values())
            for exp in exp_list
        }
        
        def convert_to_str(value):
            if use_group_idx:
                return group_idx_dict.get(value, value)
            else:
                return str(exp_names.index(value))
        
        for src, dst, data in self.ce_G.edges(data=True):
            i = node_idx[src]
            j = node_idx[dst]
            if full_ce:
                names = sorted([convert_to_str(exp.name) for exp in data["counterexamples"]])
                matrix[i][j] = f"{','.join(names)}"
            elif "counterexamples" in data and data["counterexamples"]:
                example = data["counterexamples"][0].name  # assuming .name exists
                matrix[i][j] = convert_to_str(example)


        # Build LaTeX tabular code
        nodes = [node.replace('_','\\_') for node in sorted(self.ce_G.nodes())]  # fix consistent order

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

    def to_conflict_matrix_csv(self, use_group_idx=True, full_ce=False):
        # Get sorted list of nodes for consistent ordering
        nodes = sorted(self.ce_G.nodes())
        
        # Build header
        matrix = [",".join([""] + nodes)]
        
        exp_names = self.get_exp_names()
        group_idx_dict = {
            exp.name: chr(65 + idx) + str(exp_list.index(exp) + 1) ## A1, A2, B1, B2, etc.
            for idx, exp_list in enumerate(self.exp_groups.values())
            for exp in exp_list
        }
        def convert_to_str(value):
            if use_group_idx:
                return group_idx_dict.get(value, value)
            else:
                return str(exp_names.index(value))

        # Populate the matrix
        for src in nodes:
            row = [src]
            for dst in nodes:
                if src == dst:
                    row.append("")  # No self-loop
                elif self.ce_G.has_edge(src, dst):
                    ce_list = self.ce_G.edges[src, dst]["counterexamples"]
                    if ce_list:
                        if full_ce:
                            names = sorted([convert_to_str(exp.name) for exp in ce_list])
                            row.append('"'+ ",".join(names)+'"')
                        else:
                            row.append(convert_to_str(ce_list[0].name))  # Use the name of the first counterexample
                    else:
                        row.append("")
                else:
                    row.append("")
            matrix.append(",".join(row))
        return "\n".join(matrix)