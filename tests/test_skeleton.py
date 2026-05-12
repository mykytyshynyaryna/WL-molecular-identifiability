"""Tests for skeleton.py — Notation 15 relations, skeleton S_G, and Lemma 16.

Covers:
  classify_between_class_relation  — Notation 15 biregular classifier
  classify_within_class_structure  — Lemma 13 within-class classifier
  build_skeleton                   — skeleton S_G construction
  check_lemma16_conditions         — three path conditions of Lemma 16
"""

import networkx as nx

from wl_identifiability.skeleton import (
    REL_EMPTY,
    REL_FAN,
    REL_IRREGULAR,
    REL_MATCHING,
    build_skeleton,
    check_lemma16_conditions,
    classify_between_class_relation,
    classify_within_class_structure,
)


class TestClassifyBetweenClassRelation:
    """
    Tests for the biregular relation classifier between two partition classes.
    Based on Notation 15 (Kiefer): □ (empty), ≐ (matching), ≪ (fan).
    """

    def test_empty_no_edges(self):
        """No edges between classes → REL_EMPTY."""
        F = nx.Graph()
        F.add_nodes_from([0, 1, 2, 3])
        rel, _k, _l_deg = classify_between_class_relation(F, [0, 1], [2, 3])
        assert rel == REL_EMPTY
        assert _k == 0 and _l_deg == 0

    def test_empty_isolated_nodes(self):
        """Explicit isolated node groups → REL_EMPTY."""
        F = nx.Graph()
        F.add_edges_from([(0, 1)])
        rel, _k, _l_deg = classify_between_class_relation(F, [0, 1], [2, 3])
        assert rel == REL_EMPTY

    def test_perfect_matching(self):
        """One-to-one matching between classes → REL_MATCHING, k=l=1."""
        F = nx.Graph()
        F.add_edges_from([(0, 2), (1, 3)])
        rel, k, l_deg = classify_between_class_relation(F, [0, 1], [2, 3])
        assert rel == REL_MATCHING
        assert k == 1 and l_deg == 1

    def test_fan_p_to_q(self):
        """P-vertices each have 2 neighbours in Q; Q-vertices each have 1 in P → P ≪ Q."""
        P = [0, 1, 2]
        Q = [3, 4, 5, 6, 7, 8]
        F = nx.Graph()
        F.add_edges_from([(0, 3), (0, 4), (1, 5), (1, 6), (2, 7), (2, 8)])
        rel, k, l_deg = classify_between_class_relation(F, P, Q)
        assert rel == REL_FAN
        assert k == 2 and l_deg == 1

    def test_fan_q_to_p(self):
        """Q-vertices each have 2 neighbours in P; P-vertices each have 1 → P ≫ Q (Q ≪ P)."""
        P = [0, 1, 2, 3, 4, 5]
        Q = [6, 7, 8]
        F = nx.Graph()
        F.add_edges_from([(0, 6), (1, 6), (2, 7), (3, 7), (4, 8), (5, 8)])
        rel, k, l_deg = classify_between_class_relation(F, P, Q)
        assert rel == REL_FAN
        assert k == 1 and l_deg == 2

    def test_irregular_non_biregular(self):
        """Non-uniform degrees → REL_IRREGULAR."""
        F = nx.Graph()
        F.add_edges_from([(0, 2), (0, 3), (1, 4)])
        rel, k, l_deg = classify_between_class_relation(F, [0, 1], [2, 3, 4])
        assert rel == REL_IRREGULAR
        assert k == -1 and l_deg == -1

    def test_irregular_both_high_degree(self):
        """k=2, l=2 → does not match any Notation-15 case → REL_IRREGULAR."""
        F = nx.complete_bipartite_graph(2, 2)
        P = [0, 1]
        Q = [2, 3]
        rel, _k, _l_deg = classify_between_class_relation(F, P, Q)
        assert rel == REL_IRREGULAR

    def test_empty_classes(self):
        """Empty class lists → REL_EMPTY."""
        F = nx.Graph()
        rel, _k, _l_deg = classify_between_class_relation(F, [], [0, 1])
        assert rel == REL_EMPTY


class TestClassifyWithinClassStructure:
    """
    Tests for the within-class structure classifier.
    Lemma 13 (Kiefer): valid flip within-class structures are empty, matching, cycle5.
    """

    def test_empty_no_edges(self):
        F = nx.Graph()
        F.add_nodes_from([0, 1, 2])
        assert classify_within_class_structure(F, [0, 1, 2]) == "empty"

    def test_empty_single_node(self):
        F = nx.Graph()
        F.add_node(0)
        assert classify_within_class_structure(F, [0]) == "empty"

    def test_perfect_matching_4_nodes(self):
        """Two disjoint edges on 4 nodes → perfect matching."""
        F = nx.Graph()
        F.add_edges_from([(0, 1), (2, 3)])
        assert classify_within_class_structure(F, [0, 1, 2, 3]) == "matching"

    def test_perfect_matching_6_nodes(self):
        """Three disjoint edges on 6 nodes → perfect matching."""
        F = nx.Graph()
        F.add_edges_from([(0, 1), (2, 3), (4, 5)])
        assert classify_within_class_structure(F, [0, 1, 2, 3, 4, 5]) == "matching"

    def test_cycle5(self):
        F = nx.cycle_graph(5)
        assert classify_within_class_structure(F, list(range(5))) == "cycle5"

    def test_other_triangle(self):
        F = nx.cycle_graph(3)
        assert classify_within_class_structure(F, [0, 1, 2]) == "other"

    def test_other_cycle6(self):
        F = nx.cycle_graph(6)
        assert classify_within_class_structure(F, list(range(6))) == "other"

    def test_other_odd_non_matching(self):
        """3 nodes, one edge → cannot be a perfect matching."""
        F = nx.Graph()
        F.add_edges_from([(0, 1)])
        assert classify_within_class_structure(F, [0, 1, 2]) == "other"

    def test_within_class_uses_only_p_nodes(self):
        """
        Edges outside P must not affect the classification of P.
        """
        F = nx.cycle_graph(5)
        F.add_edge(0, 5)
        assert classify_within_class_structure(F, list(range(5))) == "cycle5"


def _make_two_class_matching_flip() -> tuple[nx.Graph, dict]:
    """
    Flip graph: P=[0,1], Q=[2,3], perfect matching 0-2 and 1-3.
    Expected: REL_MATCHING between classes 0 and 1, no within-class edges.
    """
    F = nx.Graph()
    F.add_edges_from([(0, 2), (1, 3)])
    labels = {0: "A", 1: "A", 2: "B", 3: "B"}
    return F, labels


def _make_fan_flip() -> tuple[nx.Graph, dict]:
    """
    Flip graph: P=[0,1,2] (class 'A'), Q=[3,4,5,6,7,8] (class 'B').
    Each A-node fans to 2 B-nodes (k=2), each B-node hit by 1 A-node (l=1).
    Expected: REL_FAN with arc A→B.
    """
    F = nx.Graph()
    F.add_edges_from([(0, 3), (0, 4), (1, 5), (1, 6), (2, 7), (2, 8)])
    labels = {0: "A", 1: "A", 2: "A", 3: "B", 4: "B", 5: "B", 6: "B", 7: "B", 8: "B"}
    return F, labels


class TestBuildSkeleton:
    def test_two_classes_matching_has_arc_both_ways(self):
        """REL_MATCHING → two antiparallel arcs in the directed skeleton."""
        F, labels = _make_two_class_matching_flip()
        S, _relations, _within = build_skeleton(F, labels)
        assert S.has_edge("A", "B")
        assert S.has_edge("B", "A")
        assert S["A"]["B"]["relation"] == REL_MATCHING
        assert S["B"]["A"]["relation"] == REL_MATCHING

    def test_two_classes_matching_relation_stored(self):
        F, labels = _make_two_class_matching_flip()
        _, relations, _ = build_skeleton(F, labels)
        rel, k, l_deg = relations[("A", "B")]
        assert rel == REL_MATCHING
        assert k == 1 and l_deg == 1

    def test_two_classes_empty_within_structures(self):
        """No within-class edges → both classes report 'empty'."""
        F, labels = _make_two_class_matching_flip()
        _, _, within = build_skeleton(F, labels)
        assert within["A"] == "empty"
        assert within["B"] == "empty"

    def test_fan_arc_direction(self):
        """P ≪ Q (k=2, l=1) → directed arc P → Q, no arc Q → P."""
        F, labels = _make_fan_flip()
        S, _, _ = build_skeleton(F, labels)
        assert S.has_edge("A", "B")
        assert not S.has_edge("B", "A")
        assert S["A"]["B"]["relation"] == REL_FAN

    def test_fan_arc_attributes(self):
        F, labels = _make_fan_flip()
        S, _, _ = build_skeleton(F, labels)
        assert S["A"]["B"]["k"] == 2
        assert S["A"]["B"]["l"] == 1

    def test_irregular_relation_no_arc(self):
        """REL_IRREGULAR between classes → no arc in skeleton."""
        F = nx.complete_bipartite_graph(2, 2)
        labels = {0: "A", 1: "A", 2: "B", 3: "B"}
        S, relations, _ = build_skeleton(F, labels)
        assert not S.has_edge("A", "B")
        assert not S.has_edge("B", "A")
        rel, _, _ = relations[("A", "B")]
        assert rel == REL_IRREGULAR

    def test_node_attributes_set(self):
        F, labels = _make_two_class_matching_flip()
        S, _, within = build_skeleton(F, labels)
        for c in ["A", "B"]:
            assert "within_structure" in S.nodes[c]
            assert "class_size" in S.nodes[c]
            assert S.nodes[c]["within_structure"] == within[c]

    def test_within_class_cycle5_detected(self):
        """A class with a C5 internal structure is classified as 'cycle5'."""
        F = nx.cycle_graph(5)
        labels = dict.fromkeys(range(5), "A")
        _, _, within = build_skeleton(F, labels)
        assert within["A"] == "cycle5"

    def test_three_isolated_classes_no_edges(self):
        F = nx.Graph()
        F.add_nodes_from([0, 1, 2, 3, 4, 5])
        labels = {0: "A", 1: "A", 2: "B", 3: "B", 4: "C", 5: "C"}
        S, _relations, _ = build_skeleton(F, labels)
        assert S.number_of_edges() == 0
        assert S.number_of_nodes() == 3


def _skeleton_from_edges(
    fan_edges: list[tuple],
    matching_edges: list[tuple],
    within: dict,
) -> tuple[nx.DiGraph, dict]:
    """
    Build a skeleton DiGraph directly from edge specifications.

    fan_edges     : list of (u, v) — directed REL_FAN arcs u → v
    matching_edges: list of (u, v) — undirected REL_MATCHING (stored as u↔v)
    within        : dict {node → structure_label}
    """
    S = nx.DiGraph()
    all_nodes = set(within.keys())
    for u, v in fan_edges + matching_edges:
        all_nodes.update([u, v])
    for n in all_nodes:
        S.add_node(n, within_structure=within.get(n, "empty"))
    for u, v in fan_edges:
        S.add_edge(u, v, relation=REL_FAN)
    for u, v in matching_edges:
        S.add_edge(u, v, relation=REL_MATCHING)
        S.add_edge(v, u, relation=REL_MATCHING)
    return S, within


class TestCheckLemma16Conditions:
    def test_empty_skeleton_all_satisfied(self):
        """No edges, no exceptions → all conditions trivially hold."""
        S = nx.DiGraph()
        S.add_nodes_from(["A", "B", "C"])
        within = {"A": "empty", "B": "empty", "C": "empty"}
        result = check_lemma16_conditions(S, within)
        assert result["all_satisfied"] is True
        assert result["violations"] == []

    def test_single_fan_no_exceptions_all_satisfied(self):
        """One fan arc, no exception classes → conditions 1 and 2 trivially hold."""
        S, within = _skeleton_from_edges(
            fan_edges=[("A", "B")],
            matching_edges=[],
            within={"A": "empty", "B": "empty"},
        )
        result = check_lemma16_conditions(S, within)
        assert result["condition1_satisfied"] is True
        assert result["condition2_satisfied"] is True
        assert result["condition3_satisfied"] is True

    def test_single_matching_all_satisfied(self):
        S, within = _skeleton_from_edges(
            fan_edges=[],
            matching_edges=[("A", "B")],
            within={"A": "empty", "B": "empty"},
        )
        result = check_lemma16_conditions(S, within)
        assert result["all_satisfied"] is True

    def test_condition1_opposing_fans_on_path(self):
        """
        Skeleton: A → B ← C (both A and C are fan sources pointing into B).
        Undirected path A-B-C: starts with fan A→B, ends with fan C→B (backward).
        → Condition 1 violated.
        """
        S, within = _skeleton_from_edges(
            fan_edges=[("A", "B"), ("C", "B")],
            matching_edges=[],
            within={"A": "empty", "B": "empty", "C": "empty"},
        )
        result = check_lemma16_conditions(S, within)
        assert result["condition1_satisfied"] is False
        assert any("Condition 1" in v for v in result["violations"])

    def test_condition1_single_fan_source_satisfied(self):
        """A → B → C (only A is a fan source) → no opposing pair → satisfied."""
        S, within = _skeleton_from_edges(
            fan_edges=[("A", "B"), ("A", "C")],
            matching_edges=[],
            within={"A": "empty", "B": "empty", "C": "empty"},
        )
        result = check_lemma16_conditions(S, within)
        assert result["condition1_satisfied"] is True

    def test_condition1_disconnected_fans_satisfied(self):
        """Two fan sources in disconnected components → no shared path → satisfied."""
        S, within = _skeleton_from_edges(
            fan_edges=[("A", "B"), ("C", "D")],
            matching_edges=[],
            within={"A": "empty", "B": "empty", "C": "empty", "D": "empty"},
        )
        result = check_lemma16_conditions(S, within)
        assert result["condition1_satisfied"] is True

    def test_condition2_fan_reaches_cycle5(self):
        """
        Fan source A → B, and B has a cycle5 structure.
        Path A-B: starts with fan from A, ends at exception B.
        → Condition 2 violated.
        """
        S, within = _skeleton_from_edges(
            fan_edges=[("A", "B")],
            matching_edges=[],
            within={"A": "empty", "B": "cycle5"},
        )
        result = check_lemma16_conditions(S, within)
        assert result["condition2_satisfied"] is False
        assert any("Condition 2" in v for v in result["violations"])

    def test_condition2_fan_reaches_matching(self):
        """Fan source A → B, B has a matching structure → Condition 2 violated."""
        S, within = _skeleton_from_edges(
            fan_edges=[("A", "B")],
            matching_edges=[],
            within={"A": "empty", "B": "matching"},
        )
        result = check_lemma16_conditions(S, within)
        assert result["condition2_satisfied"] is False

    def test_condition2_matching_edge_to_exception_satisfied(self):
        """
        A ≐ B (matching arc, not fan), and B has a cycle5.
        The path from A to B does NOT start with a fan arc from A.
        → Condition 2 satisfied.
        """
        S, within = _skeleton_from_edges(
            fan_edges=[],
            matching_edges=[("A", "B")],
            within={"A": "empty", "B": "cycle5"},
        )
        result = check_lemma16_conditions(S, within)
        assert result["condition2_satisfied"] is True

    def test_condition2_fan_not_reaching_exception_satisfied(self):
        """
        Fan source A → B, exception at C, A and C are in different components.
        → Condition 2 satisfied (A cannot reach C).
        """
        S, within = _skeleton_from_edges(
            fan_edges=[("A", "B")],
            matching_edges=[],
            within={"A": "empty", "B": "empty", "C": "cycle5"},
        )
        result = check_lemma16_conditions(S, within)
        assert result["condition2_satisfied"] is True

    def test_condition3_two_exceptions_in_component(self):
        """
        Two cycle5 classes X and Y connected by a matching arc.
        → Condition 3 violated (two exceptions in one component).
        """
        S, within = _skeleton_from_edges(
            fan_edges=[],
            matching_edges=[("X", "Y")],
            within={"X": "cycle5", "Y": "cycle5"},
        )
        result = check_lemma16_conditions(S, within)
        assert result["condition3_satisfied"] is False
        assert any("Condition 3" in v for v in result["violations"])

    def test_condition3_two_exceptions_separate_components_satisfied(self):
        """Two cycle5 classes in separate disconnected components → satisfied."""
        S, within = _skeleton_from_edges(
            fan_edges=[],
            matching_edges=[],
            within={"X": "cycle5", "Y": "cycle5"},
        )
        result = check_lemma16_conditions(S, within)
        assert result["condition3_satisfied"] is True

    def test_condition3_one_exception_per_component_satisfied(self):
        """Two components, each with one exception → satisfied."""
        S, within = _skeleton_from_edges(
            fan_edges=[],
            matching_edges=[("A", "B"), ("C", "D")],
            within={"A": "cycle5", "B": "empty", "C": "empty", "D": "matching"},
        )
        result = check_lemma16_conditions(S, within)
        assert result["condition3_satisfied"] is True

    def test_condition3_matching_and_cycle5_same_component_violated(self):
        """One component with both a matching class and a cycle5 class → violated."""
        S, within = _skeleton_from_edges(
            fan_edges=[],
            matching_edges=[("A", "B"), ("B", "C")],
            within={"A": "matching", "B": "empty", "C": "cycle5"},
        )
        result = check_lemma16_conditions(S, within)
        assert result["condition3_satisfied"] is False

    def test_non_forest_flagged(self):
        """A skeleton with a cycle is not a forest — this is reported."""
        S = nx.DiGraph()
        S.add_edge("A", "B", relation=REL_MATCHING)
        S.add_edge("B", "A", relation=REL_MATCHING)
        S.add_edge("B", "C", relation=REL_MATCHING)
        S.add_edge("C", "B", relation=REL_MATCHING)
        S.add_edge("C", "A", relation=REL_MATCHING)
        S.add_edge("A", "C", relation=REL_MATCHING)
        within = {"A": "empty", "B": "empty", "C": "empty"}
        result = check_lemma16_conditions(S, within)
        assert result["skeleton_is_forest"] is False
        assert result["all_satisfied"] is False
        assert len(result["violations"]) >= 1

    def test_forest_skeleton_flagged_true(self):
        """A chain A-B-C (matching arcs) is a forest."""
        S, within = _skeleton_from_edges(
            fan_edges=[],
            matching_edges=[("A", "B"), ("B", "C")],
            within={"A": "empty", "B": "empty", "C": "empty"},
        )
        result = check_lemma16_conditions(S, within)
        assert result["skeleton_is_forest"] is True

    def test_flip_with_fan_and_exception_violates_condition2(self):
        """
        Construct a flip graph where class A fans into class B,
        and class B has a C5 internal structure.
        After building skeleton: condition 2 should be violated.
        """
        A_nodes = list(range(5))
        B_nodes = list(range(5, 15))
        F = nx.Graph()
        F.add_nodes_from(A_nodes + B_nodes)
        for i, a in enumerate(A_nodes):
            F.add_edge(a, B_nodes[2 * i])
            F.add_edge(a, B_nodes[2 * i + 1])
        B_cycle = B_nodes[:5]
        for i in range(5):
            F.add_edge(B_cycle[i], B_cycle[(i + 1) % 5])

        labels = dict.fromkeys(A_nodes, "A")
        labels.update(dict.fromkeys(B_nodes, "B"))

        S, _, within = build_skeleton(F, labels)
        result = check_lemma16_conditions(S, within)

        assert within["B"] == "cycle5"
        assert S.has_edge("A", "B")
        assert S["A"]["B"]["relation"] == REL_FAN
        assert result["condition2_satisfied"] is False
