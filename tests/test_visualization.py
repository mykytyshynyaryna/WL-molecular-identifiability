"""Tests for WL output processing (wl.py) and visualization helpers (visualization.py)."""

import networkx as nx
import pytest

from wl_identifiability.visualization import (
    _build_color_palette,
    _group_nodes_by_label,
    _sort_color_class_keys,
)
from wl_identifiability.wl import (
    compute_label_histogram,
    compute_wl_coloring,
    normalize_node_labels,
)


def _wl_result(G, **kwargs):
    return compute_wl_coloring(G, **kwargs)


class TestNormalizeNodeLabels:
    def test_accepts_plain_dict(self):
        G = nx.path_graph(3)
        labels = {0: "A", 1: "B", 2: "A"}
        assert normalize_node_labels(G, labels) == labels

    def test_accepts_wl_result_dict(self):
        G = nx.path_graph(4)
        result = normalize_node_labels(G, _wl_result(G))
        assert set(result.keys()) == set(G.nodes())

    def test_accepts_history_list(self):
        G = nx.path_graph(4)
        wl = _wl_result(G, store_history=True)
        result = normalize_node_labels(G, wl["history"])
        assert set(result.keys()) == set(G.nodes())

    @pytest.mark.parametrize(
        "labels,exc_type,match",
        [
            pytest.param({0: "A", 1: "B"}, ValueError, "Missing labels", id="missing_node"),
            pytest.param({0: "A", 1: "B", 99: "X"}, ValueError, "not in graph", id="extra_node"),
            pytest.param("bad_input", TypeError, None, id="wrong_type"),
            pytest.param([], ValueError, None, id="empty_history"),
        ],
    )
    def test_raises(self, labels, exc_type, match):
        G = nx.path_graph(3) if isinstance(labels, dict) and len(labels) == 2 else nx.path_graph(2)
        with pytest.raises(exc_type, match=match):
            normalize_node_labels(G, labels)


class TestComputeLabelHistogram:
    def setup_method(self):
        self.G = nx.path_graph(5)
        self.wl = _wl_result(self.G)

    @pytest.mark.parametrize(
        "key",
        [
            "n_nodes",
            "n_classes",
            "largest_class_size",
            "singleton_classes",
            "histogram",
        ],
    )
    def test_result_has_key(self, key):
        assert key in compute_label_histogram(self.G, self.wl)

    def test_n_nodes_correct(self):
        assert compute_label_histogram(self.G, self.wl)["n_nodes"] == 5

    def test_n_classes_positive(self):
        assert compute_label_histogram(self.G, self.wl)["n_classes"] >= 1

    def test_largest_class_size_leq_n_nodes(self):
        assert compute_label_histogram(self.G, self.wl)["largest_class_size"] <= 5

    def test_complete_graph_one_class(self):
        G = nx.complete_graph(4)
        result = compute_label_histogram(G, _wl_result(G))
        assert result["n_classes"] == 1
        assert result["largest_class_size"] == 4

    def test_invalid_sort_by_raises(self):
        with pytest.raises(ValueError):
            compute_label_histogram(self.G, self.wl, sort_by="nonsense")

    def test_sort_by_label(self):
        assert "histogram" in compute_label_histogram(self.G, self.wl, sort_by="label")


class TestGroupNodesByLabel:
    def test_returns_dict(self):
        G = nx.path_graph(4)
        assert isinstance(_group_nodes_by_label(G, _wl_result(G)), dict)

    def test_all_nodes_covered(self):
        G = nx.path_graph(6)
        groups = _group_nodes_by_label(G, _wl_result(G))
        assert {n for nodes in groups.values() for n in nodes} == set(G.nodes())

    def test_complete_graph_one_group(self):
        G = nx.complete_graph(5)
        groups = _group_nodes_by_label(G, _wl_result(G))
        assert len(groups) == 1
        assert len(next(iter(groups.values()))) == 5


class TestSortColorClassKeys:
    def test_returns_list_covering_all_keys(self):
        classes = {2: [0], 0: [1], 1: [2]}
        keys = _sort_color_class_keys(classes)
        assert isinstance(keys, list)
        assert set(keys) == {0, 1, 2}

    def test_deterministic(self):
        classes = {3: [], 1: [], 2: []}
        assert _sort_color_class_keys(classes) == _sort_color_class_keys(classes)


class TestBuildColorPalette:
    def test_unique_labels_get_unique_indices(self):
        assert len(set(_build_color_palette(["A", "B", "C"]).values())) == 3

    def test_duplicate_labels_collapsed(self):
        assert len(_build_color_palette(["A", "A", "B"])) == 2

    def test_indices_start_at_zero(self):
        assert min(_build_color_palette(["X", "Y"]).values()) == 0
