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
    graphs = pydot.graph_from_dot_file(input_file)
    graph: pydot.Dot = graphs[0]
    new_graph = pydot.Dot(
        graph_name=graph.get_name(),
        graph_type=graph.get_type(),
        strict=graph.get_strict(False),  # careful: call get_strict
    )

    # Sort nodes by name
    seq: int = 1
    node_list = sorted(graph.get_nodes(), key=lambda x: x.get_name())
    for node in node_list:
        # Clean attributes
        node: pydot.Node = node
        if node.get_name() in ["graph", "node"]:
            continue
        node_attrs = node.get_attributes()
        for attr in ["pos", "height", "width"]:
            node_attrs.pop(attr, None)
        nodeG: pydot.Graph = node.obj_dict["parent_graph"]
        nodeG.set_sequence(seq)
        seq += 1
        new_graph.add_node(node)

    # Sort edges by (source, destination)
    edge_list = sorted(
        graph.get_edges(), key=lambda e: (e.get_source(), e.get_destination())
    )
    for edge in edge_list:
        edge: pydot.Edge = edge
        edge_attrs = edge.get_attributes()
        edge_attrs.pop("pos", None)
        edgeG: pydot.Graph = edge.obj_dict["parent_graph"]
        edgeG.set_sequence(seq)
        seq += 1
        new_graph.add_edge(edge)

    # Clean top-level graph attributes
    graph_attrs = new_graph.get_attributes()
    graph_attrs.pop("bb", None)

    # Save to output
    new_graph.write(output_file)


def cluster_cycles(input_dot: str, output_dot: str):
    graphs = pydot.graph_from_dot_file(input_dot)
    graph = graphs[0]

    G_nx = nx.DiGraph()
    for edge in graph.get_edges():
        G_nx.add_edge(edge.get_source(), edge.get_destination())

    # Find SCCs
    sccs = list(nx.strongly_connected_components(G_nx))
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
    new_graph.write(output_dot)
