from __future__ import annotations

from ._imports import nx


def _group_nodes_by_color(labels: dict) -> dict:
    """Invert the labels dict, grouping node IDs by their WL color integer."""
    color2nodes: dict = {}
    for v, c in labels.items():
        color2nodes.setdefault(c, []).append(v)
    return color2nodes


def build_flip_graph_from_labels(G: nx.Graph, labels: dict) -> tuple[nx.Graph, dict]:
    """Construct the flip graph from the original graph and WL node labels."""
    color2nodes = _group_nodes_by_color(labels)
    colors = list(color2nodes.keys())

    nodes = list(G.nodes())
    info = {"within_copy": 0, "within_flip": 0, "between_copy": 0, "between_flip": 0}
    edges: list[tuple] = []

    for c in colors:
        Ci = color2nodes[c]
        n = len(Ci)
        if n < 2:
            continue

        S = set(Ci)
        m = sum(1 for u in Ci for v in G.neighbors(u) if v in S and u < v)
        M = n * (n - 1) // 2

        if m > M / 2:
            for i in range(n):
                for j in range(i + 1, n):
                    if not G.has_edge(Ci[i], Ci[j]):
                        edges.append((Ci[i], Ci[j]))
            info["within_flip"] += 1
        else:
            for u in Ci:
                for v in G.neighbors(u):
                    if v in S and u < v:
                        edges.append((u, v))
            info["within_copy"] += 1

    for ci in range(len(colors)):
        for cj in range(ci + 1, len(colors)):
            Ci = color2nodes[colors[ci]]
            Cj = color2nodes[colors[cj]]
            ni, nj = len(Ci), len(Cj)
            if ni == 0 or nj == 0:
                info["between_copy"] += 1
                continue
            if ni > nj:
                Ci, Cj = Cj, Ci
                ni, nj = nj, ni

            Sj = set(Cj)
            m = sum(1 for u in Ci for v in G.neighbors(u) if v in Sj)
            M = ni * nj

            if m > M / 2:
                for u in Ci:
                    for v in Cj:
                        if not G.has_edge(u, v):
                            edges.append((u, v))
                info["between_flip"] += 1
            else:
                for u in Ci:
                    for v in G.neighbors(u):
                        if v in Sj:
                            edges.append((u, v))
                info["between_copy"] += 1

    F = nx.Graph()
    F.add_nodes_from(nodes)
    F.add_edges_from(edges)
    return F, info
