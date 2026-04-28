"""Tests for flip_graph.py — flip-graph construction from WL labels."""

import networkx as nx

from wl_identifiability.flip_graph import (
    _group_nodes_by_color,
    build_flip_graph_from_labels,
)


class TestGroupNodesByColor:
    def test_basic_grouping(self):
        labels = {0: "A", 1: "B", 2: "A", 3: "C"}
        groups = _group_nodes_by_color(labels)
        assert set(groups["A"]) == {0, 2}
        assert groups["B"] == [1]
        assert groups["C"] == [3]

    def test_all_same_color(self):
        labels = {0: 1, 1: 1, 2: 1}
        groups = _group_nodes_by_color(labels)
        assert set(groups[1]) == {0, 1, 2}

    def test_all_different_colors(self):
        labels = {0: 0, 1: 1, 2: 2}
        groups = _group_nodes_by_color(labels)
        assert len(groups) == 3


class TestBuildFlipGraphFromLabels:
    def _path4_labels(self):
        G = nx.path_graph(4)
        labels = {0: 0, 1: 1, 2: 1, 3: 0}
        return G, labels

    def test_returns_graph_and_info(self):
        G, labels = self._path4_labels()
        F, info = build_flip_graph_from_labels(G, labels)
        assert isinstance(F, nx.Graph)
        assert isinstance(info, dict)

    def test_flip_graph_has_same_nodes(self):
        G, labels = self._path4_labels()
        F, _ = build_flip_graph_from_labels(G, labels)
        assert set(F.nodes()) == set(G.nodes())

    def test_info_has_expected_keys(self):
        G, labels = self._path4_labels()
        _, info = build_flip_graph_from_labels(G, labels)
        for key in ("within_copy", "within_flip", "between_copy", "between_flip"):
            assert key in info

    def test_info_counts_are_non_negative(self):
        G, labels = self._path4_labels()
        _, info = build_flip_graph_from_labels(G, labels)
        for v in info.values():
            assert v >= 0

    def test_single_color_class_no_between(self):
        G = nx.complete_graph(4)
        labels = {0: 0, 1: 0, 2: 0, 3: 0}
        _, info = build_flip_graph_from_labels(G, labels)
        assert info["between_copy"] == 0
        assert info["between_flip"] == 0

    def test_all_different_labels_no_within(self):
        G = nx.path_graph(4)
        labels = {0: 0, 1: 1, 2: 2, 3: 3}
        _, info = build_flip_graph_from_labels(G, labels)
        assert info["within_copy"] == 0
        assert info["within_flip"] == 0

    def test_dense_graph_uses_flip(self):
        G = nx.complete_graph(4)
        labels = {0: 0, 1: 0, 2: 1, 3: 1}
        _, info = build_flip_graph_from_labels(G, labels)
        assert info["between_flip"] >= 1
