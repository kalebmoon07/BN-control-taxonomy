import re
import subprocess
import networkx as nx
import pydot


def write_dot(G: nx.DiGraph, output_file: str):
    nx.nx_pydot.write_dot(G, output_file)
    clean_and_sort_dot(output_file, output_file)


def write_transitive_reduction(input_file: str, output_file: str):
    subprocess.run(f"tred {input_file} | dot -Tdot -o {output_file}", shell=True)
    clean_and_sort_dot(output_file, output_file)


def export_dot_png(input_file: str, output_file: str):
    subprocess.run(f"dot -Tpng {input_file} -o {output_file}", shell=True)


def clean_and_sort_dot(input_file: str, output_file: str):
    graph = pydot.graph_from_dot_file(input_file)[0]
    new_graph = pydot.Dot(
        graph_name="G",
        graph_type=graph.get_type(),
        strict=False,
    )

    # Canonical graph attrs
    gattrs = dict(graph.get_attributes())
    gattrs.pop("bb", None)
    for k in sorted(gattrs):
        new_graph.set(k, gattrs[k])

    # Nodes: sort + rebuild fresh with sorted attrs
    for n in sorted(graph.get_nodes(), key=lambda x: x.get_name()):
        nname = n.get_name()
        if nname in ["graph", "node"]:
            continue
        attrs = dict(n.get_attributes())
        for a in ["pos", "height", "width"]:
            attrs.pop(a, None)
        attrs = {k: attrs[k] for k in sorted(attrs)}
        new_graph.add_node(pydot.Node(nname, **attrs))

    # Edges: sort + rebuild fresh with sorted attrs
    for e in sorted(
        graph.get_edges(), key=lambda e: (e.get_source(), e.get_destination())
    ):
        attrs = dict(e.get_attributes())
        attrs.pop("pos", None)
        attrs = {k: attrs[k] for k in sorted(attrs)}
        new_graph.add_edge(pydot.Edge(e.get_source(), e.get_destination(), **attrs))

    # Write baseline DOT
    dot_str = new_graph.to_string()
    dot_str = dot_str.replace("\r\n", "\n").replace("\r", "\n")

    # Write with forced LF
    with open(output_file, "w", newline="\n", encoding="utf-8") as f:
        f.write(dot_str)



def cluster_cycles(input_dot: str, output_dot: str):
    graphs = pydot.graph_from_dot_file(input_dot)
    graph = graphs[0]

    G_nx = nx.DiGraph()
    for edge in graph.get_edges():
        G_nx.add_edge(edge.get_source(), edge.get_destination())

    # Find SCCs
    sccs = list(nx.strongly_connected_components(G_nx))
    sccs = sorted(
        [sorted(scc) for scc in sccs],  # sort nodes within each SCC
        key=lambda x: x[0],  # sort SCCs by first element
    )

    node_to_cluster = {}
    clusters = {}

    cluster_id = 0
    for scc in sccs:
        if len(scc) > 1:
            cluster_name = f"cluster_{cluster_id}"
            for node in scc:
                node_to_cluster[node] = cluster_name
            clusters[cluster_name] = scc
            cluster_id += 1

    # New graph
    new_graph = pydot.Dot(graph_type="digraph", strict=True)
    new_graph.set_graph_defaults(fontsize="10", fontname="Verdana", compound="true")
    new_graph.set_node_defaults(shape="record", fontsize="10", fontname="Verdana")
    new_graph.set_graph_defaults(nodesep=0.3, ranksep=0.8)

    # Create clusters
    for cluster_name, nodes in clusters.items():
        subg = pydot.Subgraph(graph_name=cluster_name)
        subg.set_label(f"cluster_{cluster_name.split('_')[-1]}")
        subg.set_graph_defaults(
            style="dashed",
            color="black",
        )

        for node_name in nodes:
            node = graph.get_node(node_name)[0]
            subg.add_node(node)

        new_graph.add_subgraph(subg)

    # Add standalone nodes
    clustered_nodes = {node for nodes in clusters.values() for node in nodes}
    for node in graph.get_nodes():
        node_name = node.get_name()
        if node_name not in clustered_nodes and node_name != "node":
            new_graph.add_node(node)

    # Add all edges, but annotate ltail/lhead if source/target inside a cluster
    for edge in graph.get_edges():
        src = edge.get_source()
        dst = edge.get_destination()

        attributes = {}
        if src in node_to_cluster:
            attributes["ltail"] = node_to_cluster[src]
        if dst in node_to_cluster:
            attributes["lhead"] = node_to_cluster[dst]

        if (
            (src in node_to_cluster)
            & (dst in node_to_cluster)
            & (node_to_cluster.get(src) == node_to_cluster.get(dst))
        ):
            continue
        new_edge = pydot.Edge(src, dst, **attributes)
        new_graph.add_edge(new_edge)

    # Save
    dot_str = new_graph.to_string()
    dot_str = dot_str.replace("\r\n", "\n").replace("\r", "\n")
    # remove empty lines
    lines = [ln for ln in dot_str.split("\n") if ln.strip() != ""]
    dot_str = "\n".join(lines) + "\n"
    with open(output_dot, "w", newline="\n", encoding="utf-8") as f:
        f.write(dot_str)
