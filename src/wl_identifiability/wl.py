from __future__ import annotations

from collections import Counter
from collections.abc import Hashable, Iterable
from typing import Any

import networkx as nx

LabelSource = str | Iterable[str] | None


def _normalize_label_value(x: object) -> Any:
    """Convert an arbitrary node-attribute value to a hashable, comparable form."""
    if x is None:
        return "__NONE__"
    if isinstance(x, (str, int, float, bool)):
        return x
    if isinstance(x, tuple):
        return tuple(_normalize_label_value(i) for i in x)
    if isinstance(x, list):
        return tuple(_normalize_label_value(i) for i in x)
    if isinstance(x, set):
        return tuple(sorted(_normalize_label_value(i) for i in x))
    if isinstance(x, dict):
        return tuple(sorted((k, _normalize_label_value(v)) for k, v in x.items()))
    return repr(x)


def _resolve_label_attribute(label_attr: LabelSource) -> str | list[str] | None:
    """Normalise label_attr to None, a single string, or a list of strings."""
    if label_attr is None:
        return None
    if isinstance(label_attr, str):
        s = label_attr.strip()
        if "," in s:
            parts = [p.strip() for p in s.split(",") if p.strip()]
            return parts
        return s
    return list(label_attr)


def _initialize_node_labels(
    G: nx.Graph,
    label_attr: LabelSource = None,
    missing: str = "__MISSING__",
    default: str = "__UNLABELED__",
) -> dict[Hashable, int]:
    """Read node attributes from G and map them to compact integer labels at iteration 0."""
    label_attr = _resolve_label_attribute(label_attr)

    raw: dict[Any, Any] = {}
    for v in G.nodes:
        if label_attr is None:
            raw[v] = default
        elif isinstance(label_attr, str):
            raw[v] = G.nodes[v].get(label_attr, missing)
        else:
            raw[v] = tuple(G.nodes[v].get(a, missing) for a in label_attr)

        raw[v] = _normalize_label_value(raw[v])

    uniq = sorted(set(raw.values()), key=lambda z: repr(z))
    mp: dict[Any, int] = {val: i for i, val in enumerate(uniq)}
    return {v: mp[raw[v]] for v in G.nodes}


def _refine_wl_labels_once(G: nx.Graph, labels: dict[Hashable, int]) -> dict[Hashable, int]:
    """Perform one WL refinement step: hash each node's label together with its sorted neighbour labels."""
    sig: dict[Any, tuple[int, tuple[int, ...]]] = {}
    for v in G.nodes:
        neigh = sorted(labels[u] for u in G.neighbors(v))
        sig[v] = (labels[v], tuple(neigh))

    uniq = sorted(set(sig.values()), key=lambda z: repr(z))
    mp: dict[Any, int] = {val: i for i, val in enumerate(uniq)}
    return {v: mp[sig[v]] for v in G.nodes}


def compute_wl_coloring(
    G: nx.Graph,
    label_attr: LabelSource = None,
    max_iter: int | None = None,
    store_history: bool = True,
    mode: str = "bounded",
) -> dict[str, Any]:
    """
    Compute Weisfeiler-Lehman node coloring for the given graph.
    """
    if G.number_of_nodes() == 0:
        return {
            "converged": True,
            "converge_iter": 0,
            "labels": {},
            "history": [],
            "iterations": 0,
            "iter_cap": 0,
        }

    n = G.number_of_nodes()

    if mode == "fixed":
        if max_iter is None:
            raise ValueError("mode='fixed' requires max_iter (K).")
        iter_cap = int(max_iter)
        stop_on_converge = False
    elif mode in ("until_converge", "bounded"):
        iter_cap = n if max_iter is None else int(max_iter)
        stop_on_converge = True
    else:
        raise ValueError(f"Unknown mode: {mode}")

    labels: dict[Hashable, int] = _initialize_node_labels(G, label_attr=label_attr)

    history: list[dict[Hashable, int]] = []
    if store_history:
        history.append(labels.copy())

    converged = False
    converge_iter: int | None = None
    iterations = 0

    for t in range(1, iter_cap + 1):
        new_labels = _refine_wl_labels_once(G, labels)
        iterations = t

        if store_history:
            history.append(new_labels.copy())

        if (not converged) and (len(set(new_labels.values())) == len(set(labels.values()))):
            converged = True
            converge_iter = t
            if stop_on_converge:
                labels = new_labels
                break

        labels = new_labels

    return {
        "converged": converged,
        "converge_iter": converge_iter,
        "labels": labels,
        "history": history,
        "iterations": iterations,
        "iter_cap": iter_cap,
    }


def compute_fixed_wl_steps_from_topology(G: nx.Graph, cap: int | None = None) -> int:
    """Compute K from topology (no label_attr), then run exactly K steps."""
    if G.number_of_nodes() == 0:
        return 0

    max_steps = G.number_of_nodes() if cap is None else min(cap, G.number_of_nodes())
    out = compute_wl_coloring(
        G,
        label_attr=None,
        mode="until_converge",
        max_iter=max_steps,
        store_history=False,
    )
    K = out["converge_iter"] if out["converge_iter"] is not None else out["iterations"]
    return int(K)


def compute_wl_with_fixed_steps(G: nx.Graph, K: int, label_attr: LabelSource = None) -> dict[str, Any]:
    """
    Run WL coloring for a fixed number of iterations.
    """
    return compute_wl_coloring(G, label_attr=label_attr, mode="fixed", max_iter=K, store_history=False)


def compute_wl_stabilization_steps(
    G: nx.Graph,
    *,
    label_attr: LabelSource = None,
    cap: int | None = None,
) -> int:
    """
    Compute how many WL iterations are needed until stabilization.
    If WL does not converge before the cap, returns the number of executed iterations.
    """
    out = compute_wl_coloring(
        G,
        label_attr=label_attr,
        max_iter=cap,
        store_history=False,
        mode="until_converge",
    )
    return int(out["converge_iter"]) if out["converged"] else int(out["iterations"])


def normalize_node_labels(
    G: nx.Graph,
    labels: Any,
) -> dict[Hashable, Hashable]:
    """Normalise labels input: accepts a plain dict, a WL result dict, or a history list."""
    if isinstance(labels, dict):
        if "labels" in labels and isinstance(labels["labels"], dict):
            labels = labels["labels"]
        elif "history" in labels and isinstance(labels["history"], list):
            if not labels["history"]:
                raise ValueError("WL history is empty")
            labels = labels["history"][-1]

    elif isinstance(labels, list):
        if not labels:
            raise ValueError("labels history is empty")
        labels = labels[-1]

    if not isinstance(labels, dict):
        raise TypeError("labels must be a dict, history list, or WL result dict")

    missing = [v for v in G.nodes() if v not in labels]
    if missing:
        preview = missing[:5]
        suffix = "..." if len(missing) > 5 else ""
        raise ValueError(f"Missing labels for nodes: {preview}{suffix}")

    extra = [v for v in labels if v not in G]
    if extra:
        preview = extra[:5]
        suffix = "..." if len(extra) > 5 else ""
        raise ValueError(f"Labels contain nodes not in graph: {preview}{suffix}")

    return {v: labels[v] for v in G.nodes()}


def compute_label_histogram(
    G: nx.Graph,
    labels: Any,
    *,
    sort_by: str = "size_desc",
) -> dict[str, Any]:
    """
    Compute a summary of WL color classes.
    """
    labels = normalize_node_labels(G, labels)
    hist = Counter(labels.values())

    if sort_by == "size_desc":
        items = sorted(hist.items(), key=lambda x: (-x[1], repr(x[0])))
    elif sort_by == "label":
        items = sorted(hist.items(), key=lambda x: repr(x[0]))
    else:
        raise ValueError("sort_by must be 'size_desc' or 'label'")

    hist_sorted = dict(items)

    return {
        "n_nodes": G.number_of_nodes(),
        "n_classes": len(hist),
        "largest_class_size": max(hist.values()) if hist else 0,
        "singleton_classes": sum(1 for c in hist.values() if c == 1),
        "histogram": hist_sorted,
    }
