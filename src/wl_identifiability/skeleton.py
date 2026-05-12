"""
Skeleton graph S_G construction and Lemma 16 validation.

Implements the theoretical objects from Kiefer's paper that constrain the
structure of C²-identified graphs:

  Notation 15  — biregular relations between partition classes (□, ≐, ≪)
  Section 4    — skeleton graph S_G over C²-partition classes
  Lemma 16     — three structural path conditions on S_G

Theorem 17 (Kiefer): G is identified by C² iff its flip is a bouquet forest.
The skeleton and Lemma 16 conditions are intermediate steps in the proof of
the necessary direction (identified → bouquet forest structure).
"""

from __future__ import annotations

from typing import Any

import networkx as nx

REL_EMPTY: str = "empty"

REL_MATCHING: str = "matching"

REL_FAN: str = "fan"

REL_IRREGULAR: str = "irregular"


def classify_between_class_relation(
    F: nx.Graph,
    P: list[Any],
    Q: list[Any],
) -> tuple[str, int, int]:
    """
    Classify the biregular edge relation between partition classes P and Q
    in flip graph F, per Notation 15 (Kiefer).

    Parameters
    ----------
    F : flip graph (nx.Graph)
    P : list of nodes belonging to one partition class
    Q : list of nodes belonging to another partition class

    Returns
    -------
    (relation, k, l) where
      k = uniform degree of every P-vertex toward Q  (-1 if irregular)
      l = uniform degree of every Q-vertex toward P  (-1 if irregular)

    Relations
    ---------
    REL_EMPTY    : k=0, l=0
    REL_MATCHING : k=1, l=1      (P ≐ Q)
    REL_FAN      : k≥2, l=1      (P ≪ Q)  or  k=1, l≥2  (P ≫ Q = Q ≪ P)
    REL_IRREGULAR: non-biregular, or k≥2 and l≥2, etc.

    Lemma 14 (Kiefer): C² identifies a flipped biregular graph iff k≤1 or l≤1.
    REL_IRREGULAR is returned whenever Lemma 14 would be violated.
    """
    if not P or not Q:
        return REL_EMPTY, 0, 0

    Qset = set(Q)
    Pset = set(P)

    k_vals = [sum(1 for v in F.neighbors(u) if v in Qset) if F.has_node(u) else 0 for u in P]
    l_vals = [sum(1 for v in F.neighbors(u) if v in Pset) if F.has_node(u) else 0 for u in Q]

    if len(set(k_vals)) != 1 or len(set(l_vals)) != 1:
        return REL_IRREGULAR, -1, -1

    k, l_deg = k_vals[0], l_vals[0]

    if k == 0 and l_deg == 0:
        return REL_EMPTY, 0, 0
    if k == 1 and l_deg == 1:
        return REL_MATCHING, 1, 1
    if (k >= 2 and l_deg == 1) or (k == 1 and l_deg >= 2):
        return REL_FAN, k, l_deg

    return REL_IRREGULAR, k, l_deg


def classify_within_class_structure(F: nx.Graph, P: list[Any]) -> str:
    """
    Classify the induced subgraph on partition class P in flip graph F.

    Returns one of:
      'empty'    — no edges (isolated vertices)
      'matching' — perfect matching on all nodes of P
      'cycle5'   — exactly one induced C₅ covering all 5 nodes of P
      'other'    — any other structure

    Lemma 13 (Kiefer): for a C²-identified flipped graph, the only valid
    within-class induced structures are 'empty', 'matching', and 'cycle5'.
    'other' signals a Lemma 13 violation.
    """
    if len(P) <= 1:
        return "empty"

    H = F.subgraph(P).copy()
    H.add_nodes_from(P)
    n = H.number_of_nodes()
    m = H.number_of_edges()

    if m == 0:
        return "empty"

    degrees = [d for _, d in H.degree()]

    if all(d == 1 for d in degrees) and n % 2 == 0 and m == n // 2:
        return "matching"

    active = [v for v in H.nodes() if H.degree(v) > 0]
    if len(active) == 5:
        H5 = H.subgraph(active)
        if H5.number_of_edges() == 5 and all(d == 2 for _, d in H5.degree()) and nx.is_connected(H5):
            return "cycle5"

    return "other"


def build_skeleton(
    F: nx.Graph,
    labels: dict[Any, Any],
) -> tuple[nx.DiGraph, dict[tuple[Any, Any], tuple[str, int, int]], dict[Any, str]]:
    """
    Build skeleton graph S_G (Kiefer, Section 4) from flip graph F and
    the C²-partition assignment ``labels``.

    The skeleton S_G abstracts the inter-class connectivity of the flip:
      - Vertices : one per C²-partition class (keyed by its color label)
      - Directed arcs for REL_FAN (≪): P → Q when P ≪ Q (P fans into Q)
      - Antiparallel arcs for REL_MATCHING (≐): P ↔ Q stored as two arcs
      - No arc for REL_EMPTY (□) or REL_IRREGULAR

    Node attributes on skeleton nodes:
      within_structure : 'empty' | 'matching' | 'cycle5' | 'other'
      class_size       : number of original nodes in this class

    Edge attributes on skeleton arcs:
      relation : REL_MATCHING or REL_FAN
      k, l     : degree values from the biregular structure

    Parameters
    ----------
    F      : flip graph (nx.Graph), output of build_flip_graph_from_labels
    labels : dict {node → color}, the C²-partition assignment used to build F

    Returns
    -------
    skeleton              : nx.DiGraph over color labels
    class_relations       : dict {(ci, cj) → (relation, k, l)} for all
                            unordered class pairs (ci < cj in insertion order)
    within_class_structures : dict {ci → structure_label}
    """
    color2nodes: dict[Any, list[Any]] = {}
    for v, c in labels.items():
        color2nodes.setdefault(c, []).append(v)

    colors = list(color2nodes.keys())

    within_class_structures: dict[Any, str] = {c: classify_within_class_structure(F, color2nodes[c]) for c in colors}

    S: nx.DiGraph = nx.DiGraph()
    for c in colors:
        S.add_node(
            c,
            within_structure=within_class_structures[c],
            class_size=len(color2nodes[c]),
        )

    class_relations: dict[tuple[Any, Any], tuple[str, int, int]] = {}

    for i in range(len(colors)):
        ci = colors[i]
        Pi = color2nodes[ci]
        for j in range(i + 1, len(colors)):
            cj = colors[j]
            Pj = color2nodes[cj]

            rel, k, l_deg = classify_between_class_relation(F, Pi, Pj)
            class_relations[(ci, cj)] = (rel, k, l_deg)

            if rel == REL_MATCHING:
                S.add_edge(ci, cj, relation=REL_MATCHING, k=k, l=l_deg)
                S.add_edge(cj, ci, relation=REL_MATCHING, k=l_deg, l=k)
            elif rel == REL_FAN:
                if k >= 2 and l_deg == 1:
                    S.add_edge(ci, cj, relation=REL_FAN, k=k, l=l_deg)
                else:
                    S.add_edge(cj, ci, relation=REL_FAN, k=l_deg, l=k)

    return S, class_relations, within_class_structures


def check_lemma16_conditions(
    skeleton: nx.DiGraph,
    within_class_structures: dict[Any, str],
) -> dict[str, Any]:
    """
    Check the three structural conditions of Lemma 16 (Kiefer) on S_G.

    These conditions are necessary for C²-identification of the underlying
    graph. Their combined effect forces S_G to be a forest with at most one
    'exception class' (matching or 5-cycle) per connected component, and
    prevents opposing fan sources from appearing on the same path.

    Conditions
    ----------
    1. No undirected path P_0, ..., P_t in S_G where
       P_0 ≪ P_1  (REL_FAN arc P_0 → P_1) and
       P_{t-1} ≫ P_t  (REL_FAN arc P_t → P_{t-1}).
       Prevents two opposing fan sources on the same path.

    2. No undirected path P_0, ..., P_t where
       P_0 ≪ P_1  (REL_FAN arc P_0 → P_1) and
       P_t has within-class structure 'matching' or 'cycle5'.
       Prevents fan sources from reaching exception classes.

    3. Each connected component of S_G contains at most one class with
       within-class structure 'matching' or 'cycle5'.
       Allows at most one 'exception' per component.

    Note: conditions 1 and 2 are checked via shortest paths in the undirected
    skeleton. For acyclic skeletons (forests) this is exact — there is only
    one path between any two nodes. For skeletons with cycles the check may
    be incomplete; an acyclicity violation is reported separately.

    Parameters
    ----------
    skeleton               : nx.DiGraph, output of build_skeleton
    within_class_structures : dict {color → structure_label}

    Returns
    -------
    dict with keys:
      skeleton_is_forest   : bool
      condition1_satisfied : bool
      condition2_satisfied : bool
      condition3_satisfied : bool
      all_satisfied        : bool
      violations           : list[str]  — human-readable descriptions
    """
    violations: list[str] = []

    undirected = skeleton.to_undirected()

    skeleton_is_forest = nx.is_forest(undirected)
    if not skeleton_is_forest:
        violations.append(
            "Skeleton S_G contains cycles; a forest is required for C²-identification. "
            "Conditions 1 and 2 below are checked via shortest paths only."
        )

    fan_targets: dict[Any, set[Any]] = {}
    for u, v, data in skeleton.edges(data=True):
        if data.get("relation") == REL_FAN:
            fan_targets.setdefault(u, set()).add(v)

    exception_nodes: set[Any] = {c for c, s in within_class_structures.items() if s in ("matching", "cycle5")}

    cond1 = True
    fan_source_list = list(fan_targets.keys())
    for i in range(len(fan_source_list)):
        u = fan_source_list[i]
        for j in range(i + 1, len(fan_source_list)):
            v = fan_source_list[j]
            if not nx.has_path(undirected, u, v):
                continue
            path = nx.shortest_path(undirected, u, v)
            if len(path) < 2:
                continue
            first_is_fan = path[1] in fan_targets.get(u, set())
            last_is_fan_back = path[-2] in fan_targets.get(v, set())
            if first_is_fan and last_is_fan_back:
                cond1 = False
                violations.append(
                    f"Condition 1 violated: opposing fan sources {u} (→{path[1]}) and {v} (→{path[-2]}) on path {path}"
                )

    cond2 = True
    for u, targets in fan_targets.items():
        for ex in exception_nodes:
            if u == ex:
                continue
            if not nx.has_path(undirected, u, ex):
                continue
            path = nx.shortest_path(undirected, u, ex)
            if len(path) < 2:
                continue
            if path[1] in targets:
                cond2 = False
                violations.append(
                    f"Condition 2 violated: fan source {u} (→{path[1]}) "
                    f"reaches exception class {ex} "
                    f"(structure='{within_class_structures.get(ex)}') "
                    f"via path {path}"
                )

    cond3 = True
    for comp in nx.connected_components(undirected):
        ex_in_comp = [c for c in comp if c in exception_nodes]
        if len(ex_in_comp) > 1:
            cond3 = False
            violations.append(
                f"Condition 3 violated: {len(ex_in_comp)} exception classes in one component: {ex_in_comp}"
            )

    all_ok = skeleton_is_forest and cond1 and cond2 and cond3

    return {
        "skeleton_is_forest": skeleton_is_forest,
        "condition1_satisfied": cond1,
        "condition2_satisfied": cond2,
        "condition3_satisfied": cond3,
        "all_satisfied": all_ok,
        "violations": violations,
    }
