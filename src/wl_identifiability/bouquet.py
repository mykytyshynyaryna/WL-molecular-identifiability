from __future__ import annotations

from ._imports import Hashable, nx

from .wl import compute_wl_coloring



def _find_unique_induced_c5_cycle(G: nx.Graph) -> list | None:
    """Return the unique induced C5 node list via cycle_basis, or None if no such cycle exists."""
    cycles = nx.cycle_basis(G)
    if len(cycles) != 1:
        return None
    cyc = cycles[0]
    if len(cyc) != 5:
        return None
    C = set(cyc)
    H = G.subgraph(C)
    if H.number_of_edges() != 5:
        return None
    if any(H.degree(v) != 2 for v in H.nodes()):
        return None
    return list(cyc)


def _extract_cycle_edges(G: nx.Graph, cycle_nodes: list) -> list[tuple]:
    """Return all edges of G whose both endpoints are in cycle_nodes."""
    edges = []
    for u in cycle_nodes:
        for v in cycle_nodes:
            if u >= v:
                continue
            if G.has_edge(u, v):
                edges.append((u, v))
    return edges


def _build_graph_without_edges(G: nx.Graph, edges: list[tuple]) -> nx.Graph:
    """Return a new graph with the same nodes/attributes as G, minus the given edges."""
    excluded = {(min(u, v), max(u, v)) for u, v in edges}
    H = nx.create_empty_copy(G)
    H.add_edges_from(
        (u, v, d) for u, v, d in G.edges(data=True)
        if (min(u, v), max(u, v)) not in excluded
    )
    return H

def _compute_canonical_cycle_order(G: nx.Graph, cycle_nodes: list) -> list:
    """Return cycle_nodes in a canonical traversal order (smallest node first, lexicographically earlier direction)."""
    n = len(cycle_nodes)
    C = set(cycle_nodes)
    H = G.subgraph(C)

    a = min(cycle_nodes)
    nbrs = sorted(H.neighbors(a))
    if len(nbrs) != 2:
        return sorted(cycle_nodes)

    def _build_order(start_next: object) -> list | None:
        order = [a, start_next]
        while len(order) < n:
            prev, cur = order[-2], order[-1]
            nxts = [x for x in H.neighbors(cur) if x != prev]
            if len(nxts) != 1:
                return None
            order.append(nxts[0])
        return order

    o1, o2 = _build_order(nbrs[0]), _build_order(nbrs[1])
    valid = [o for o in (o1, o2) if o is not None and len(o) == n]
    if len(valid) == 2:
        return valid[0] if tuple(valid[0]) <= tuple(valid[1]) else valid[1]
    return valid[0] if valid else sorted(cycle_nodes)

def _is_valid_tree_structure(G: nx.Graph) -> bool:
    """Return True if G is a connected tree."""
    n = G.number_of_nodes()
    if n == 0:
        return False
    if G.number_of_edges() != n - 1:
        return False
    return nx.is_connected(G)


def _split_bouquet_component_into_petals(G: nx.Graph):
    """
    Split G into a canonical C₅ order and five rooted petal subgraphs.

    Returns (True, cyc_ord, petals) on success, or (False, None, None) if G
    does not have exactly one induced C₅ with five valid tree petals.
    """
    if G.number_of_nodes() == 0:
        return False, None, None

    cyc = _find_unique_induced_c5_cycle(G)
    if cyc is None:
        return False, None, None

    cyc_ord = _compute_canonical_cycle_order(G, cyc)
    C = set(cyc_ord)

    cycle_edges = _extract_cycle_edges(G, cyc_ord)
    H = _build_graph_without_edges(G, cycle_edges)

    petals = {}
    for comp_nodes in nx.connected_components(H):
        comp = set(comp_nodes)
        inter = list(comp.intersection(C))
        if len(inter) != 1:
            return False, None, None
        root = inter[0]
        T = H.subgraph(comp).copy()
        if not _is_valid_tree_structure(T):
            return False, None, None
        petals[root] = T

    if len(petals) != 5 or len(set(petals.keys())) != 5:
        return False, None, None

    return True, cyc_ord, petals


def _are_all_petals_isomorphic(petals: dict, labels: dict | None = None) -> bool:
    """
    Return True iff all five petal trees are mutually isomorphic (Definition 9, Kiefer).

    If labels are provided, isomorphism is tested as exact vertex-colored rooted-tree
    signature comparison; otherwise uses the uncolored AHU canonical signature.
    Both paths are exact (no WL approximation).
    """
    petal_items = list(petals.items())
    if len(petal_items) < 2:
        return True

    if labels is not None:
        sigs = [str(_rooted_colored_signature(T, root, labels)) for root, T in petal_items]
    else:
        sigs = [str(rooted_tree_signature(T, root)) for root, T in petal_items]
    return len(set(sigs)) == 1




def _sig_from_bouquet_result(result: dict) -> str | None:
    """
    Extract the canonical AHU petal signature from an already-validated bouquet result dict.

    All five petal signatures are equal after acceptance, so ``tree_sigs[0]`` is the
    canonical certificate. Returns None only if the result carries no tree_signatures.
    """
    tree_sigs = result.get("tree_signatures") or []
    return tree_sigs[0] if tree_sigs else None


def _compute_bouquet_signature(G: nx.Graph, labels: dict | None = None) -> str | None:
    """
    Return a canonical AHU isomorphism certificate for bouquet G, or None if G is not a bouquet.

    Since all five petals are mutually isomorphic by the bouquet definition, the AHU
    signature of any single petal uniquely identifies the bouquet up to isomorphism —
    signature equality replaces full graph isomorphism checks entirely.
    """
    result = is_bouquet_component(G, method="optimized", labels=labels)
    if not result["is_bouquet"]:
        return None
    return _sig_from_bouquet_result(result)


def _not_bouquet_result(method: str, reason: str) -> dict:
    """Return a structured negative-result dict with all fields set to None."""
    return {
        "is_bouquet": False,
        "method": method,
        "reason": reason,
        "cycle_nodes": None,
        "cycle_length": None,
        "n_rooted_trees": None,
        "tree_signatures": None,
    }


def _add_root_marker(T: nx.Graph, root: Hashable, labels: dict | None) -> nx.Graph:
    """Return a copy of T with ``__is_root`` and (if labels given) ``__color`` node attributes."""
    H = T.copy()
    for v in H.nodes():
        H.nodes[v]["__is_root"] = v == root
        if labels is not None:
            H.nodes[v]["__color"] = labels.get(v, "__missing__")
    return H


def _extract_petals(G: nx.Graph, cyc_ord: list) -> tuple[bool, dict | None]:
    """
    Remove cycle edges from G and return five rooted petal subgraphs keyed by root.

    Returns (True, petals_dict) on success, or (False, None) if the components
    do not form exactly five valid trees each attached through one cycle vertex.
    """
    C = set(cyc_ord)
    cycle_edges = _extract_cycle_edges(G, cyc_ord)
    G_no_cycle = _build_graph_without_edges(G, cycle_edges)

    petals: dict[Hashable, nx.Graph] = {}
    for comp_nodes in nx.connected_components(G_no_cycle):
        comp = set(comp_nodes)
        inter = list(comp.intersection(C))
        if len(inter) != 1:
            return False, None
        root = inter[0]
        T = G_no_cycle.subgraph(comp).copy()
        if not _is_valid_tree_structure(T):
            return False, None
        petals[root] = T

    if len(petals) != 5:
        return False, None

    return True, petals




def rooted_tree_signature(T: nx.Graph, root: Hashable) -> tuple:
    """
    Return an exact canonical nested-tuple signature for rooted tree (T, root).

    Two rooted trees are isomorphic iff their signatures are equal; no WL
    approximation is involved.
    """
    def _sig(v: Hashable, parent: Hashable | None) -> tuple:
        return tuple(sorted(_sig(c, v) for c in T.neighbors(v) if c != parent))

    return _sig(root, None)


def _rooted_colored_signature(T: nx.Graph, root: Hashable, labels: dict) -> tuple:
    """Return a canonical AHU signature for a rooted tree with per-node colors from labels."""
    def _sig(v: Hashable, parent: Hashable | None) -> tuple:
        color = labels.get(v, None)
        return (color, tuple(sorted(_sig(c, v) for c in T.neighbors(v) if c != parent)))

    return _sig(root, None)



def _find_c5_nodes_adj(adj: dict) -> list | None:
    """Leaf-stripping C₅ finder on a plain adjacency dict — no NX view overhead."""
    degree = {v: len(nbs) for v, nbs in adj.items()}
    queue = [v for v, d in degree.items() if d == 1]
    removed: set = set()
    while queue:
        v = queue.pop()
        removed.add(v)
        for nb in adj[v]:
            if nb not in removed:
                degree[nb] -= 1
                if degree[nb] == 1:
                    queue.append(nb)
    cycle_nodes = [v for v in adj if v not in removed]
    if len(cycle_nodes) != 5:
        return None
    if any(degree[v] != 2 for v in cycle_nodes):
        return None
    return cycle_nodes


def _canonical_cycle_order_adj(adj: dict, cycle_nodes: list) -> list:
    """Canonical C₅ traversal order using plain adjacency dict (no subgraph view)."""
    C = set(cycle_nodes)
    a = min(cycle_nodes)
    nbrs = sorted(nb for nb in adj[a] if nb in C)
    if len(nbrs) != 2:
        return sorted(cycle_nodes)

    def _build_order(start_next: object) -> list | None:
        order = [a, start_next]
        while len(order) < 5:
            prev, cur = order[-2], order[-1]
            nxts = [x for x in adj[cur] if x != prev and x in C]
            if len(nxts) != 1:
                return None
            order.append(nxts[0])
        return order

    o1, o2 = _build_order(nbrs[0]), _build_order(nbrs[1])
    valid = [o for o in (o1, o2) if o and len(o) == 5]
    if len(valid) == 2:
        return valid[0] if tuple(valid[0]) <= tuple(valid[1]) else valid[1]
    return valid[0] if valid else sorted(cycle_nodes)


def _petal_signatures_adj(
    adj: dict,
    cycle_nodes: list,
    labels: dict | None,
) -> tuple[bool, list | None]:
    """BFS petal discovery and AHU signatures on a plain adjacency dict."""
    C: set = set(cycle_nodes)
    cycle_edge_set: set = {
        (min(u, v), max(u, v))
        for u in C
        for v in adj[u]
        if v in C
    }

    visited: set = set()
    comps: list = []
    for start in adj:
        if start in visited:
            continue
        comp: list = []
        queue = [start]
        visited.add(start)
        qi = 0
        while qi < len(queue):
            v = queue[qi]
            qi += 1
            comp.append(v)
            for nb in adj[v]:
                if nb not in visited and (min(v, nb), max(v, nb)) not in cycle_edge_set:
                    visited.add(nb)
                    queue.append(nb)
        comps.append(comp)

    if len(comps) != 5:
        return False, None

    petal_pairs: list = []
    for comp in comps:
        roots = [n for n in comp if n in C]
        if len(roots) != 1:
            return False, None
        petal_pairs.append((roots[0], set(comp)))

    if labels is not None:
        def _sig(v: Hashable, parent: Hashable | None, pset: set) -> tuple:
            color = labels.get(v, None)
            return (color, tuple(sorted(
                _sig(c, v, pset) for c in adj[v] if c != parent and c in pset
            )))
    else:
        def _sig(v: Hashable, parent: Hashable | None, pset: set) -> tuple:
            return tuple(sorted(
                _sig(c, v, pset) for c in adj[v] if c != parent and c in pset
            ))

    return True, [str(_sig(r, None, ns)) for r, ns in petal_pairs]


def _check_bouquet_adj(
    adj: dict,
    n: int,
    m: int,
    labels: dict | None,
) -> dict:
    """Core bouquet check on a plain adjacency dict with pre-computed n and m."""
    method = "optimized"
    if m == n - 1:
        return _not_bouquet_result(method, "is_tree")
    if m != n:
        return _not_bouquet_result(method, "not_one_cycle_component")

    cyc = _find_c5_nodes_adj(adj)
    if cyc is None:
        return _not_bouquet_result(method, "no_unique_induced_c5")

    ok, str_sigs = _petal_signatures_adj(adj, cyc, labels)
    if not ok:
        return _not_bouquet_result(method, "invalid_petal_structure")

    if len(set(str_sigs)) != 1:
        return {
            "is_bouquet": False,
            "method": method,
            "reason": "petals_not_isomorphic",
            "cycle_nodes": cyc,
            "cycle_length": 5,
            "n_rooted_trees": 5,
            "tree_signatures": str_sigs,
        }

    return {
        "is_bouquet": True,
        "method": method,
        "reason": "ok",
        "cycle_nodes": cyc,
        "cycle_length": 5,
        "n_rooted_trees": 5,
        "tree_signatures": str_sigs,
    }



def _check_bouquet_component_baseline(
    G: nx.Graph,
    labels: dict | None = None,
) -> dict:
    """
    Reference bouquet checker using nx.cycle_basis and nx.is_isomorphic.

    Verifies connectivity, a unique induced C₅, five valid petal trees, and
    pairwise isomorphism of those petals (with color matching when labels given).
    """
    method = "baseline"

    if not nx.is_connected(G):
        return _not_bouquet_result(method, "not_connected")

    if nx.is_tree(G):
        return _not_bouquet_result(method, "is_tree")

    cycles = nx.cycle_basis(G)
    if len(cycles) != 1:
        return _not_bouquet_result(method, "not_exactly_one_cycle")

    cyc = cycles[0]
    if len(cyc) != 5:
        return _not_bouquet_result(method, "cycle_length_not_5")

    C = set(cyc)
    H_cyc = G.subgraph(C)
    if H_cyc.number_of_edges() != 5 or any(
        H_cyc.degree(v) != 2 for v in H_cyc.nodes()
    ):
        return _not_bouquet_result(method, "cycle_not_induced_c5")

    cyc_ord = _compute_canonical_cycle_order(G, cyc)
    ok, petals = _extract_petals(G, cyc_ord)
    if not ok:
        return _not_bouquet_result(method, "invalid_petal_structure")

    petal_list = list(petals.items())
    ref_root, ref_T = petal_list[0]

    if labels is not None:
        node_match = lambda a, b: (
            a["__is_root"] == b["__is_root"] and a["__color"] == b["__color"]
        )
    else:
        node_match = lambda a, b: a["__is_root"] == b["__is_root"]

    ref_marked = _add_root_marker(ref_T, ref_root, labels)
    for root, T in petal_list[1:]:
        if not nx.is_isomorphic(
            _add_root_marker(T, root, labels), ref_marked, node_match=node_match
        ):
            if labels is not None:
                sigs = [str(_rooted_colored_signature(t, r, labels)) for r, t in petal_list]
            else:
                sigs = [str(rooted_tree_signature(t, r)) for r, t in petal_list]
            return {
                "is_bouquet": False,
                "method": method,
                "reason": "petals_not_isomorphic",
                "cycle_nodes": cyc_ord,
                "cycle_length": 5,
                "n_rooted_trees": len(petals),
                "tree_signatures": sigs,
            }

    if labels is not None:
        sigs = [str(_rooted_colored_signature(t, r, labels)) for r, t in petal_list]
    else:
        sigs = [str(rooted_tree_signature(t, r)) for r, t in petal_list]
    return {
        "is_bouquet": True,
        "method": method,
        "reason": "ok",
        "cycle_nodes": cyc_ord,
        "cycle_length": 5,
        "n_rooted_trees": 5,
        "tree_signatures": sigs,
    }




def _check_bouquet_component_optimized(
    G: nx.Graph,
    labels: dict | None = None,
) -> dict:
    """
    Fast bouquet checker using edge-count pre-filter, leaf-stripping, and AHU signatures.

    Materialises the adjacency dict once from G (avoiding repeated subgraph-view
    coreviews overhead), then delegates entirely to the pure-dict helpers.
    """
    adj = {v: list(G.adj[v]) for v in G.nodes()}
    n = len(adj)
    m = sum(len(nbs) for nbs in adj.values()) // 2
    return _check_bouquet_adj(adj, n, m, labels)




def is_bouquet_component(
    G: nx.Graph,
    method: str = "optimized",
    labels: dict | None = None,
) -> dict:
    """
    Check whether connected component G is a valid bouquet (Definition 9, Kiefer).

    A bouquet is five mutually isomorphic rooted trees whose roots form a C₅.
    Pass labels (C²-partition colors) to enforce vertex-colored isomorphism.

    Returns a dict with keys: is_bouquet, method, reason, cycle_nodes,
    cycle_length, n_rooted_trees, tree_signatures.
    """
    if method == "baseline":
        return _check_bouquet_component_baseline(G, labels=labels)
    if method == "optimized":
        return _check_bouquet_component_optimized(G, labels=labels)
    raise ValueError(
        f"Unknown method {method!r}. Supported values: 'baseline', 'optimized'."
    )



_COMPONENT_REASON_MAP: dict[str, str] = {
    "is_tree":                  "component_not_tree_or_bouquet",
    "not_one_cycle_component":  "multiple_cycles_in_component",
    "no_unique_induced_c5":     "cycle_not_5",
    "invalid_petal_structure":  "invalid_petal_structure",
    "petals_not_isomorphic":    "petals_not_isomorphic",
    "not_connected":            "disconnected_invalid_structure",
    "not_exactly_one_cycle":    "multiple_cycles_in_component",
    "cycle_length_not_5":       "cycle_not_5",
    "cycle_not_induced_c5":     "not_induced_c5",
}


def analyze_bouquet_forest_structure(
    F: nx.Graph,
    labels: dict | None = None,
    method: str = "optimized",
) -> dict:
    """
    Test whether flip graph F is a bouquet forest (Theorem 17, Kiefer).

    Returns True iff every non-tree component is a valid bouquet and all
    bouquets are pairwise non-isomorphic — the necessary and sufficient
    condition for C²-identifiability of the original graph.

    Returns a dict with keys:
      is_bouquet_forest    — True iff F satisfies Theorem 17 (all non-tree
                             components are pairwise non-isomorphic bouquets).
      has_bouquet_component — True iff at least one connected component of F
                             is a valid bouquet, regardless of whether the
                             whole graph is a bouquet forest.  This field is
                             independent of is_bouquet_forest.
      bouquets             — list of canonical petal signatures for each
                             bouquet component found (may be partial when
                             invalid components are also present).
      non_identifiable     — True iff two bouquets share the same signature.
      reason               — short code describing why is_bouquet_forest=False,
                             or "ok" when True.

    All components are always visited so that has_bouquet_component reflects
    the full flip graph, even when an invalid component is encountered first.

    Reason codes for is_bouquet_forest=False:
      cycle_not_5                  — unicyclic component whose cycle length ≠ 5
      multiple_cycles_in_component — component has more than one cycle
      invalid_petal_structure      — petals do not form five valid rooted trees
      petals_not_isomorphic        — petal trees have different structures
      not_induced_c5               — five-cycle found but it is not induced
      disconnected_invalid_structure — baseline: component is not connected
      component_not_tree_or_bouquet  — fallback for unrecognised structural issue
      duplicate_bouquet_signature  — two bouquets share the same isomorphism class
    """
    bouquets: list[str] = []
    has_bouquet_component: bool = False
    invalid_reason: str | None = None

    for comp_nodes in nx.connected_components(F):
        n_comp = len(comp_nodes)
        m_comp = sum(len(F[v]) for v in comp_nodes) // 2
        if m_comp == n_comp - 1:
            continue

        if method == "optimized":
            adj = {v: list(F[v]) for v in comp_nodes}
            result = _check_bouquet_adj(adj, n_comp, m_comp, labels)
        else:
            H = F.subgraph(comp_nodes)
            result = _check_bouquet_component_baseline(H, labels=labels)

        if not result["is_bouquet"]:
            if invalid_reason is None:
                internal = result.get("reason", "unknown_error")
                invalid_reason = _COMPONENT_REASON_MAP.get(
                    internal, "component_not_tree_or_bouquet"
                )
            continue

        sig = _sig_from_bouquet_result(result)
        if sig is None:
            if invalid_reason is None:
                invalid_reason = "component_not_tree_or_bouquet"
            continue

        bouquets.append(sig)
        has_bouquet_component = True

    if invalid_reason is not None:
        return {
            "is_bouquet_forest": False,
            "has_bouquet_component": has_bouquet_component,
            "bouquets": bouquets,
            "non_identifiable": False,
            "reason": invalid_reason,
        }

    seen: set = set()
    non_identifiable = False
    for s in bouquets:
        if s in seen:
            non_identifiable = True
            break
        seen.add(s)

    if non_identifiable:
        return {
            "is_bouquet_forest": False,
            "has_bouquet_component": True,
            "bouquets": bouquets,
            "non_identifiable": True,
            "reason": "duplicate_bouquet_signature",
        }

    return {
        "is_bouquet_forest": True,
        "has_bouquet_component": has_bouquet_component,
        "bouquets": bouquets,
        "non_identifiable": False,
        "reason": "ok",
    }


def check_bouquet_component(
    G: nx.Graph,
    labels: dict | None = None,
) -> tuple[bool, dict | None]:
    """
    Legacy public API — prefer ``is_bouquet_component`` for new code.

    Returns (True, info_dict) if G is a valid bouquet, else (False, None).
    info_dict keys: cycle, petals, isomorphic_petals.
    """
    if not nx.is_connected(G):
        return False, None

    ok, cycle, petals = _split_bouquet_component_into_petals(G)
    if not ok:
        return False, None

    if not _are_all_petals_isomorphic(petals, labels=labels):
        return False, None

    return True, {"cycle": cycle, "petals": petals, "isomorphic_petals": True}


def compare_bouquet_methods(G: nx.Graph, labels: dict | None = None) -> dict:
    """Run both checkers on G and return their results for side-by-side comparison."""
    return {
        "baseline": _check_bouquet_component_baseline(G, labels=labels),
        "optimized": _check_bouquet_component_optimized(G, labels=labels),
    }
