"""Tests for wl.py — Weisfeiler-Lehman coloring."""

from collections import Counter

import networkx as nx
import pytest

from wl_identifiability.wl import (
    _initialize_node_labels,
    _refine_wl_labels_once,
    compute_fixed_wl_steps_from_topology,
    compute_wl_coloring,
    compute_wl_stabilization_steps,
    compute_wl_with_fixed_steps,
)


class TestInitializeNodeLabels:
    def test_unlabeled_assigns_same_label_to_all(self):
        G = nx.path_graph(4)
        labels = _initialize_node_labels(G)
        assert set(labels.keys()) == set(G.nodes())
        assert len(set(labels.values())) == 1

    def test_single_attr_distinguishes_nodes(self):
        G = nx.Graph()
        G.add_node(0, color="red")
        G.add_node(1, color="blue")
        G.add_edge(0, 1)
        labels = _initialize_node_labels(G, label_attr="color")
        assert labels[0] != labels[1]

    def test_missing_attr_uses_sentinel(self):
        G = nx.Graph()
        G.add_node(0, color="red")
        G.add_node(1)
        labels = _initialize_node_labels(G, label_attr="color", missing="__MISSING__")
        assert labels[0] != labels[1]


class TestRefineWlLabelsOnce:
    def test_symmetric_graph_keeps_same_number_of_classes(self):
        G = nx.cycle_graph(6)
        labels = _initialize_node_labels(G)
        new_labels = _refine_wl_labels_once(G, labels)
        assert len(set(new_labels.values())) == 1

    def test_path_endpoints_distinguished_from_interior(self):
        G = nx.path_graph(5)
        labels = _initialize_node_labels(G)
        refined = _refine_wl_labels_once(G, labels)
        assert refined[0] == refined[4]
        assert refined[1] == refined[3]
        assert refined[0] != refined[2]


class TestComputeWlColoring:
    def test_empty_graph_returns_converged(self):
        result = compute_wl_coloring(nx.Graph())
        assert result["converged"] is True
        assert result["labels"] == {}

    def test_single_node_converges_immediately(self):
        G = nx.Graph()
        G.add_node(0)
        result = compute_wl_coloring(G)
        assert result["converged"] is True
        assert 0 in result["labels"]

    def test_complete_graph_all_same_label(self):
        result = compute_wl_coloring(nx.complete_graph(5))
        assert len(set(result["labels"].values())) == 1

    def test_path_graph_converges(self):
        result = compute_wl_coloring(nx.path_graph(6))
        assert result["converged"] is True
        assert result["converge_iter"] is not None

    def test_fixed_mode_runs_exact_iterations(self):
        result = compute_wl_coloring(nx.path_graph(10), mode="fixed", max_iter=3)
        assert result["iterations"] == 3

    def test_isomorphic_graphs_same_histogram(self):
        r1 = compute_wl_coloring(nx.cycle_graph(6))
        r2 = compute_wl_coloring(nx.cycle_graph(6))
        assert Counter(r1["labels"].values()) == Counter(r2["labels"].values())

    def test_non_isomorphic_graphs_different_histograms(self):
        r1 = compute_wl_coloring(nx.path_graph(5))
        r2 = compute_wl_coloring(nx.cycle_graph(5))
        assert Counter(r1["labels"].values()) != Counter(r2["labels"].values())

    @pytest.mark.parametrize(
        "store_history,expect_empty",
        [
            pytest.param(True, False, id="history_stored"),
            pytest.param(False, True, id="history_suppressed"),
        ],
    )
    def test_history_flag(self, store_history, expect_empty):
        result = compute_wl_coloring(nx.path_graph(4), store_history=store_history)
        if expect_empty:
            assert result["history"] == []
        else:
            assert len(result["history"]) >= 1

    @pytest.mark.parametrize(
        "kwargs",
        [
            pytest.param({"mode": "fixed", "max_iter": None}, id="fixed_without_max_iter"),
            pytest.param({"mode": "bad_mode", "max_iter": 5}, id="unknown_mode"),
        ],
    )
    def test_raises_value_error(self, kwargs):
        with pytest.raises(ValueError):
            compute_wl_coloring(nx.path_graph(4), **kwargs)


class TestComputeFixedWlStepsFromTopology:
    def test_empty_graph_returns_zero(self):
        assert compute_fixed_wl_steps_from_topology(nx.Graph()) == 0

    def test_single_node_returns_non_negative(self):
        G = nx.Graph()
        G.add_node(0)
        assert compute_fixed_wl_steps_from_topology(G) >= 0

    def test_returns_int(self):
        assert isinstance(compute_fixed_wl_steps_from_topology(nx.path_graph(8)), int)

    def test_cap_is_respected(self):
        assert compute_fixed_wl_steps_from_topology(nx.path_graph(20), cap=3) <= 3


class TestComputeWlStabilizationSteps:
    def test_complete_graph_stabilizes_in_one_step(self):
        assert compute_wl_stabilization_steps(nx.complete_graph(5)) == 1

    def test_path_graph_returns_positive_int(self):
        steps = compute_wl_stabilization_steps(nx.path_graph(6))
        assert isinstance(steps, int)
        assert steps >= 1


class TestComputeWlWithFixedSteps:
    @pytest.mark.parametrize("K", [0, 1, 4])
    def test_runs_exactly_k_steps(self, K):
        result = compute_wl_with_fixed_steps(nx.path_graph(8), K=K)
        assert result["iterations"] == K
