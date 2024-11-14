import json
import os
import mpbn
from colomoto.minibn import BooleanNetwork
from itertools import product
import networkx as nx
from bntaxonomy.dep.converters import bn_of_asynchronous_transition_graph

l1 = 0.8
l2 = 2.0
l3 = 4.0
topology = {
    2: {"11": (l1, l1), "01": (-l1, l1), "00": (-l1, -l1), "10": (l1, -l1)},
    3: {
        "110": (l1, l1),
        "010": (-l1, l1),
        "000": (-l1, -l1),
        "100": (l1, -l1),
        "111": (l2, l2),
        "011": (-l2, l2),
        "001": (-l2, -l2),
        "101": (l2, -l2),
    },
}
drawing_options = {
    "node_size": 1000,
    "node_color": "none",
    "edgecolors": "k",
    "linewidths": 2.0,
}


class InstanceGen:
    def __init__(self, inst_name, size: int):
        self.inst_name = inst_name
        self.size = size
        self.states = list("".join(values) for values in product("01", repeat=size))
        self.G = nx.DiGraph()
        self.G.add_nodes_from(self.states)

    def add_edge(self, source, target):
        if sum(1 for a, b in zip(source, target) if a != b) != 1:
            raise ValueError(f"Invalid edge: {source} -> {target}")
        self.G.add_edge(source, target)

    def show_stg(self):
        nx.draw(self.G, pos=topology[self.size], with_labels=True, **drawing_options)

    def to_bnet(self):
        return bn_of_asynchronous_transition_graph(
            self.G, [f"x{i}" for i in range(1, 1 + self.size)]
        )

    def save(self, fpath: str, inputs=dict(), target=dict()):
        os.makedirs(fpath, exist_ok=True)
        self.to_bnet().save(f"{fpath}/transition_formula.bnet")
        with open(f"{fpath}/setting.json", "w") as _f:
            json.dump({"inputs": inputs, "target": target}, _f)

        # self.G.save(path)
