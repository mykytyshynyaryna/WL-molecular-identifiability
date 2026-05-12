from __future__ import annotations

from collections import defaultdict
from collections.abc import Hashable
from typing import Any

import matplotlib.pyplot as plt
import networkx as nx

from .wl import (
    compute_label_histogram,
    compute_wl_coloring,
    compute_wl_stabilization_steps,
    normalize_node_labels,
)


def _group_nodes_by_label(
    G: nx.Graph,
    labels: Any,
) -> dict[Hashable, list[Hashable]]:
    """Return a dict mapping each label to the sorted list of nodes carrying that label."""
    labels = normalize_node_labels(G, labels)

    classes: dict[Any, list[Any]] = defaultdict(list)
    for v in G.nodes():
        classes[labels[v]].append(v)

    grouped: dict[Hashable, list[Hashable]] = {}
    for lab, nodes in classes.items():
        grouped[lab] = sorted(nodes, key=lambda x: repr(x))

    return grouped


def _sort_color_class_keys(classes: dict[Hashable, list[Hashable]]) -> list[Hashable]:
    """Return the keys of classes in a stable, type-then-repr order for deterministic display."""
    return sorted(
        classes.keys(),
        key=lambda x: (str(type(x)), repr(x)),
    )


def _build_color_palette(unique_labels: list[Hashable]) -> dict[Hashable, int]:
    """Map each unique label to a distinct integer index for use as a matplotlib color index."""
    unique_sorted = sorted(set(unique_labels), key=repr)
    return {lab: i for i, lab in enumerate(unique_sorted)}


def _map_labels_to_colors(
    labels: dict[Hashable, Hashable],
    color_map: dict[Hashable, int],
    nodelist: list[Hashable],
) -> list[int]:
    """Translate each node's label to its palette integer, ordered according to nodelist."""
    return [color_map[labels[v]] for v in nodelist]


def draw_graph_wl_coloring(
    G: nx.Graph,
    labels: Any,
    *,
    layout: str = "spring",
    seed: int = 42,
    title: str | None = None,
    show_node_ids: bool = True,
    show_wl_labels: bool = False,
    pos: dict[Any, Any] | None = None,
    node_size: int = 700,
    font_size: int = 9,
) -> tuple[dict[Any, Any], dict[Hashable, int]]:
    """
    Draw a graph with nodes colored by their WL labels.
    """
    labels = normalize_node_labels(G, labels)
    nodelist = list(G.nodes())

    color_map = _build_color_palette(list(labels.values()))
    node_colors = _map_labels_to_colors(labels, color_map, nodelist)

    if pos is None:
        if layout == "spring":
            pos = nx.spring_layout(G, seed=seed)
        elif layout == "kamada_kawai":
            pos = nx.kamada_kawai_layout(G)
        elif layout == "circular":
            pos = nx.circular_layout(G)
        else:
            raise ValueError("layout must be 'spring', 'kamada_kawai', or 'circular'")

    plt.figure(figsize=(7, 5))
    if title:
        plt.title(title)

    nx.draw(
        G,
        pos,
        nodelist=nodelist,
        node_color=node_colors,
        with_labels=show_node_ids,
        node_size=node_size,
        font_size=font_size,
    )

    if show_wl_labels:
        wl_text = {v: str(labels[v])[:8] for v in nodelist}
        shifted_pos = {v: (xy[0], xy[1] - 0.08) for v, xy in pos.items()}
        nx.draw_networkx_labels(
            G,
            shifted_pos,
            labels=wl_text,
            font_size=max(font_size - 1, 6),
        )

    plt.axis("off")
    plt.tight_layout()
    plt.show()

    return pos, color_map


def inspect_wl_behavior(
    G: nx.Graph,
    *,
    fixed_wl_steps: int | None = None,
    atom_attr: str = "atomic_num",
    draw: bool = True,
) -> dict[str, Any]:
    """
    Diagnostic helper for comparing topological WL vs atom-aware WL.
    """
    if G.number_of_nodes() == 0:
        raise ValueError("Graph is empty")

    if fixed_wl_steps is None:
        fixed_wl_steps = compute_wl_stabilization_steps(G, label_attr=None)

    wl_top = compute_wl_coloring(
        G,
        label_attr=None,
        max_iter=fixed_wl_steps,
        store_history=False,
        mode="fixed",
    )

    wl_atom = compute_wl_coloring(
        G,
        label_attr=atom_attr,
        max_iter=fixed_wl_steps,
        store_history=False,
        mode="fixed",
    )

    top_hist = compute_label_histogram(G, wl_top)
    atom_hist = compute_label_histogram(G, wl_atom)

    top_classes = _group_nodes_by_label(G, wl_top)
    atom_classes = _group_nodes_by_label(G, wl_atom)

    result: dict[str, Any] = {
        "fixed_wl_steps": fixed_wl_steps,
        "stabilization_top": compute_wl_stabilization_steps(G, label_attr=None),
        "stabilization_atom": compute_wl_stabilization_steps(G, label_attr=atom_attr),
        "topological": {
            "summary": top_hist,
            "classes": {k: top_classes[k] for k in _sort_color_class_keys(top_classes)},
        },
        "atomic": {
            "summary": atom_hist,
            "classes": {k: atom_classes[k] for k in _sort_color_class_keys(atom_classes)},
        },
        "color_ratio_atom_to_top": (
            atom_hist["n_classes"] / top_hist["n_classes"] if top_hist["n_classes"] > 0 else None
        ),
    }

    print("=== WL INSPECTION ===")
    print(f"Fixed WL steps: {result['fixed_wl_steps']}")
    print(f"Top stabilization: {result['stabilization_top']}")
    print(f"Atom stabilization: {result['stabilization_atom']}")
    print()
    print("Topological WL summary:")
    print(top_hist)
    print()
    print("Atom-aware WL summary:")
    print(atom_hist)

    if draw:
        draw_graph_wl_coloring(
            G,
            wl_top,
            title="Topological WL coloring",
            show_node_ids=True,
            show_wl_labels=False,
        )
        draw_graph_wl_coloring(
            G,
            wl_atom,
            title=f"Atom-aware WL coloring ({atom_attr})",
            show_node_ids=True,
            show_wl_labels=False,
        )

    return result
