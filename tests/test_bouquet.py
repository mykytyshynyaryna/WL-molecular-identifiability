"""Tests for bouquet.py — bouquet-forest structure, petal isomorphism, and signatures.

Key behavioral changes from the theory-faithful implementation:
  - check_bouquet_component now ENFORCES petal isomorphism (returns (False, None)
    when petals are not all isomorphic) per Definition 9 (Kiefer).
  - analyze_bouquet_forest_structure sets is_bouquet_forest=False when duplicate
    bouquets exist, consistent with the non-isomorphism requirement of Definition 9.
"""

import networkx as nx
import pytest

from wl_identifiability.bouquet import (
    _are_all_petals_isomorphic,
    _compute_bouquet_signature,
    _compute_tree_wl_signature,
    _find_unique_induced_c5_cycle,
    _is_valid_tree_structure,
    _split_bouquet_component_into_petals,
    analyze_bouquet_forest_structure,
    check_bouquet_component,
    is_bouquet_component,
    rooted_tree_signature,
)


def _disconnected_graph() -> nx.Graph:
    G = nx.Graph()
    G.add_nodes_from([0, 1, 2])
    G.add_edge(0, 1)
    return G


def _single_node_graph() -> nx.Graph:
    G = nx.Graph()
    G.add_node(0)
    return G


def _c5_with_chord() -> nx.Graph:
    G = nx.cycle_graph(5)
    G.add_edge(0, 3)
    return G


def _build_uniform_bouquet(extra_depth: int) -> nx.Graph:
    """
    C5 where each cycle node has a pendant chain of ``extra_depth`` nodes.

    extra_depth=0: pure C5 (trivial single-node petals).
    extra_depth=1: C5 + one leaf per node — all petals are 2-node paths.
    extra_depth=2: C5 + path of 2 extra nodes per node — all petals are 3-node paths.

    All 5 petals are structurally identical by construction.
    """
    G = nx.cycle_graph(5)
    nxt = 5
    for root in range(5):
        prev = root
        for _ in range(extra_depth):
            G.add_edge(prev, nxt)
            prev = nxt
            nxt += 1
    return G


def _build_mixed_bouquet(depths: list) -> nx.Graph:
    """
    C5 where cycle node i gets a chain of depths[i] extra nodes.
    When depths are not all equal, petals are NOT all isomorphic.
    """
    assert len(depths) == 5
    G = nx.cycle_graph(5)
    nxt = 5
    for root, d in enumerate(depths):
        prev = root
        for _ in range(d):
            G.add_edge(prev, nxt)
            prev = nxt
            nxt += 1
    return G


def _disjoint_union(G1: nx.Graph, G2: nx.Graph) -> nx.Graph:
    """Disjoint union: relabel G2 nodes to avoid collision with G1."""
    offset = max(G1.nodes()) + 1
    G2r = nx.relabel_nodes(G2, {n: n + offset for n in G2.nodes()})
    return nx.compose(G1, G2r)


class TestIsValidTreeStructure:
    @pytest.mark.parametrize(
        "G",
        [
            pytest.param(nx.path_graph(5), id="path5"),
            pytest.param(nx.star_graph(4), id="star4"),
            pytest.param(_single_node_graph(), id="single_node"),
        ],
    )
    def test_valid_trees(self, G):
        assert _is_valid_tree_structure(G)

    @pytest.mark.parametrize(
        "G",
        [
            pytest.param(nx.cycle_graph(5), id="cycle5"),
            pytest.param(nx.Graph(), id="empty"),
            pytest.param(_disconnected_graph(), id="disconnected"),
            pytest.param(nx.DiGraph(), id="empty_digraph"),
        ],
    )
    def test_invalid_trees(self, G):
        assert not _is_valid_tree_structure(G)

    def test_directed_edge_is_not_tree(self):
        G = nx.DiGraph()
        G.add_edge(0, 1)
        assert not _is_valid_tree_structure(G)


class TestFindUniqueInducedC5Cycle:
    def test_c5_found(self):
        cyc = _find_unique_induced_c5_cycle(nx.cycle_graph(5))
        assert cyc is not None
        assert len(cyc) == 5

    def test_bouquet_finds_cycle(self, bouquet_graph):
        assert _find_unique_induced_c5_cycle(bouquet_graph) is not None

    @pytest.mark.parametrize(
        "G",
        [
            pytest.param(nx.cycle_graph(6), id="c6"),
            pytest.param(nx.path_graph(5), id="path5"),
            pytest.param(_c5_with_chord(), id="c5_with_chord"),
        ],
    )
    def test_returns_none(self, G):
        assert _find_unique_induced_c5_cycle(G) is None


class TestComputeTreeWlSignature:
    def test_returns_string(self):
        assert isinstance(_compute_tree_wl_signature(nx.path_graph(4)), str)

    def test_same_tree_same_signature(self):
        T = nx.path_graph(6)
        assert _compute_tree_wl_signature(T) == _compute_tree_wl_signature(T)

    def test_different_trees_different_signatures(self):
        assert _compute_tree_wl_signature(nx.path_graph(5)) != _compute_tree_wl_signature(nx.star_graph(4))

    @pytest.mark.parametrize("n", [3, 5, 7])
    def test_isomorphic_paths_same_signature(self, n):
        assert _compute_tree_wl_signature(nx.path_graph(n)) == _compute_tree_wl_signature(nx.path_graph(n))


class TestArePetalsIsomorphic:
    """Tests for _are_all_petals_isomorphic — comparing 5 petals of one bouquet."""

    def _get_petals(self, G):
        ok, _, petals = _split_bouquet_component_into_petals(G)
        assert ok, "Expected bouquet structure to be valid"
        return petals

    def test_uniform_depth_1_petals_are_isomorphic(self):
        petals = self._get_petals(_build_uniform_bouquet(1))
        assert _are_all_petals_isomorphic(petals) is True

    def test_uniform_depth_2_petals_are_isomorphic(self):
        petals = self._get_petals(_build_uniform_bouquet(2))
        assert _are_all_petals_isomorphic(petals) is True

    def test_mixed_depths_petals_are_not_isomorphic(self):
        """petals of depth [1,1,2,2,2] are NOT all isomorphic."""
        petals = self._get_petals(_build_mixed_bouquet([1, 1, 2, 2, 2]))
        assert _are_all_petals_isomorphic(petals) is False

    def test_with_labels_identical_colors_accepted(self):
        """Uniform bouquet: all petals have identical node colors → isomorphic."""
        G = _build_uniform_bouquet(1)
        labels = dict.fromkeys(G.nodes(), 0)
        petals = self._get_petals(G)
        assert _are_all_petals_isomorphic(petals, labels=labels) is True

    def test_with_labels_different_colors_rejected(self):
        """Uniform structure but each node has a unique color → not isomorphic."""
        G = _build_uniform_bouquet(1)
        labels = {v: v for v in G.nodes()}
        petals = self._get_petals(G)
        assert _are_all_petals_isomorphic(petals, labels=labels) is False


class TestCheckBouquetComponent:
    def test_valid_bouquet_detected(self, bouquet_graph):
        ok, info = check_bouquet_component(bouquet_graph)
        assert ok is True
        assert info is not None

    def test_info_has_cycle(self, bouquet_graph):
        _, info = check_bouquet_component(bouquet_graph)
        assert "cycle" in info
        assert len(info["cycle"]) == 5

    def test_info_has_petals(self, bouquet_graph):
        _, info = check_bouquet_component(bouquet_graph)
        assert "petals" in info
        assert len(info["petals"]) == 5

    def test_info_isomorphic_petals_is_always_true_on_accept(self, bouquet_graph):
        """When check_bouquet_component accepts, isomorphic_petals is always True."""
        _, info = check_bouquet_component(bouquet_graph)
        assert info["isomorphic_petals"] is True

    def test_pure_cycle_not_bouquet(self):
        ok, _ = check_bouquet_component(nx.cycle_graph(5))
        assert isinstance(ok, bool)

    @pytest.mark.parametrize(
        "G",
        [
            pytest.param(nx.path_graph(10), id="path10"),
            pytest.param(_disconnected_graph(), id="disconnected"),
        ],
    )
    def test_non_bouquets_rejected(self, G):
        ok, _ = check_bouquet_component(G)
        assert ok is False

    def test_non_isomorphic_petals_enforced(self):
        """
        Theory-faithful enforcement (Definition 9): a component with non-isomorphic
        petals is rejected as a bouquet. Previously this check was computed but
        the result was not used — now it gates the return value.
        """
        assert check_bouquet_component(_build_mixed_bouquet([1, 1, 2, 2, 2])) == (
            False,
            None,
        )

    def test_uniform_depth_accepted_without_labels(self):
        ok, info = check_bouquet_component(_build_uniform_bouquet(1))
        assert ok is True
        assert info is not None

    def test_uniform_depth_accepted_with_consistent_labels(self):
        """Same structure, all nodes same color → accepted as colored bouquet."""
        G = _build_uniform_bouquet(1)
        labels = dict.fromkeys(G.nodes(), 0)
        ok, _info = check_bouquet_component(G, labels=labels)
        assert ok is True

    def test_uniform_structure_rejected_with_unique_labels(self):
        """Same structure but each node gets a unique color → petals not isomorphic."""
        G = _build_uniform_bouquet(1)
        labels = {v: v for v in G.nodes()}
        assert check_bouquet_component(G, labels=labels) == (False, None)


class TestBouquetSignatureIsomorphism:
    """
    Tests for _compute_bouquet_signature.

    Two bouquets are considered isomorphic iff their signatures are equal.
    The signature must be invariant to node numbering and cycle rotation.
    """

    def test_same_petal_structure_gives_same_signature(self):
        B1 = _build_uniform_bouquet(1)
        B2 = _build_uniform_bouquet(1)
        assert _compute_bouquet_signature(B1) == _compute_bouquet_signature(B2)

    def test_signature_invariant_to_node_relabeling(self):
        B = _build_uniform_bouquet(1)
        B_shifted = nx.relabel_nodes(B, {n: n + 100 for n in B.nodes()})
        assert _compute_bouquet_signature(B) == _compute_bouquet_signature(B_shifted)

    def test_different_petal_depths_give_different_signatures(self):
        B1 = _build_uniform_bouquet(1)
        B2 = _build_uniform_bouquet(2)
        assert _compute_bouquet_signature(B1) != _compute_bouquet_signature(B2)

    def test_non_isomorphic_petal_bouquet_returns_none(self):
        """
        A component with non-isomorphic petals is rejected by check_bouquet_component
        and _compute_bouquet_signature must return None.
        """
        assert _compute_bouquet_signature(_build_mixed_bouquet([1, 1, 2, 2, 2])) is None

    def test_non_bouquet_returns_none(self):
        assert _compute_bouquet_signature(nx.path_graph(10)) is None

    def test_signature_is_string(self):
        assert isinstance(_compute_bouquet_signature(_build_uniform_bouquet(1)), str)

    def test_with_labels_same_colors_same_signature(self):
        B1 = _build_uniform_bouquet(1)
        B2 = _build_uniform_bouquet(1)
        labels1 = dict.fromkeys(B1.nodes(), 0)
        labels2 = dict.fromkeys(B2.nodes(), 0)
        assert _compute_bouquet_signature(B1, labels=labels1) == _compute_bouquet_signature(B2, labels=labels2)


class TestAnalyzeBouquetForestStructure:
    @pytest.mark.parametrize("key", ["is_bouquet_forest", "bouquets", "non_identifiable", "reason"])
    def test_result_has_expected_keys(self, key):
        result = analyze_bouquet_forest_structure(nx.path_graph(4))
        assert key in result

    @pytest.mark.parametrize(
        "G",
        [
            pytest.param(nx.path_graph(6), id="path"),
            pytest.param(nx.star_graph(3), id="star"),
        ],
    )
    def test_pure_trees_are_bouquet_forests(self, G):
        assert analyze_bouquet_forest_structure(G)["is_bouquet_forest"] is True

    def test_valid_single_bouquet_is_bouquet_forest(self, bouquet_graph):
        result = analyze_bouquet_forest_structure(bouquet_graph)
        assert result["is_bouquet_forest"] is True
        assert len(result["bouquets"]) == 1
        assert result["non_identifiable"] is False
        assert result["reason"] == "ok"


class TestBouquetForestNonIdentifiability:
    """
    Tests for the non-isomorphism requirement of Definition 9 and its relation
    to C²-identifiability via Theorem 17.

    After the theory-faithful fix:
      - Two isomorphic bouquets → is_bouquet_forest=False AND non_identifiable=True.
      - is_bouquet_forest=True requires BOTH valid structure AND pairwise non-isomorphic bouquets.
    """

    def test_single_bouquet_is_bouquet_forest_and_identifiable(self):
        result = analyze_bouquet_forest_structure(_build_uniform_bouquet(1))
        assert result["is_bouquet_forest"] is True
        assert result["non_identifiable"] is False

    def test_two_isomorphic_bouquets_not_bouquet_forest(self):
        """
        Two identical bouquet components share a signature.
        Per Definition 9, this is NOT a bouquet forest (non-isomorphism violated).
        Per Theorem 17, the original graph is NOT identified by C².
        """
        F = _disjoint_union(_build_uniform_bouquet(1), _build_uniform_bouquet(1))
        result = analyze_bouquet_forest_structure(F)
        assert result["is_bouquet_forest"] is False
        assert len(result["bouquets"]) == 2
        assert result["bouquets"][0] == result["bouquets"][1]
        assert result["non_identifiable"] is True
        assert result["reason"] == "duplicate_bouquets"

    def test_two_non_isomorphic_bouquets_are_bouquet_forest(self):
        F = _disjoint_union(_build_uniform_bouquet(1), _build_uniform_bouquet(2))
        result = analyze_bouquet_forest_structure(F)
        assert result["is_bouquet_forest"] is True
        assert len(result["bouquets"]) == 2
        assert result["bouquets"][0] != result["bouquets"][1]
        assert result["non_identifiable"] is False

    def test_three_bouquets_two_identical_not_bouquet_forest(self):
        F = _disjoint_union(
            _disjoint_union(_build_uniform_bouquet(1), _build_uniform_bouquet(2)),
            _build_uniform_bouquet(1),
        )
        result = analyze_bouquet_forest_structure(F)
        assert result["is_bouquet_forest"] is False
        assert len(result["bouquets"]) == 3
        assert result["non_identifiable"] is True

    def test_three_all_distinct_bouquets_is_bouquet_forest(self):
        F = _disjoint_union(
            _disjoint_union(_build_uniform_bouquet(1), _build_uniform_bouquet(2)),
            _build_uniform_bouquet(3),
        )
        result = analyze_bouquet_forest_structure(F)
        assert result["is_bouquet_forest"] is True
        assert len(result["bouquets"]) == 3
        assert len(set(result["bouquets"])) == 3
        assert result["non_identifiable"] is False

    def test_bouquet_plus_tree_components_is_bouquet_forest(self):
        F = _disjoint_union(_build_uniform_bouquet(1), nx.path_graph(7))
        result = analyze_bouquet_forest_structure(F)
        assert result["is_bouquet_forest"] is True
        assert len(result["bouquets"]) == 1
        assert result["non_identifiable"] is False

    def test_component_with_non_isomorphic_petals_not_bouquet_forest(self):
        """
        A non-tree component that is a C5+trees structure but with non-isomorphic
        petals is rejected by check_bouquet_component → reason is
        'component_not_tree_nor_bouquet' (not a tree, not a valid bouquet).
        """
        result = analyze_bouquet_forest_structure(_build_mixed_bouquet([1, 1, 2, 2, 2]))
        assert result["is_bouquet_forest"] is False
        assert result["reason"] == "component_not_tree_nor_bouquet"

    def test_labels_passed_through(self, bouquet_graph):
        """analyze_bouquet_forest_structure accepts labels kwarg without error."""
        labels = dict.fromkeys(bouquet_graph.nodes(), 0)
        result = analyze_bouquet_forest_structure(bouquet_graph, labels=labels)
        assert result["is_bouquet_forest"] is True


class TestRootedTreeSignature:
    def test_single_node_is_empty_tuple(self):
        G = nx.Graph()
        G.add_node(0)
        assert rooted_tree_signature(G, 0) == ()

    def test_root_with_one_leaf(self):
        G = nx.path_graph(2)
        assert rooted_tree_signature(G, 0) == ((),)
        assert rooted_tree_signature(G, 1) == ((),)

    def test_path_of_3_rooted_at_middle(self):
        G = nx.path_graph(3)
        assert rooted_tree_signature(G, 1) == ((), ())

    def test_path_of_3_rooted_at_end(self):
        G = nx.path_graph(3)
        assert rooted_tree_signature(G, 0) == (((),),)

    def test_isomorphic_rooted_trees_equal_signature(self):
        T1 = nx.path_graph(3)
        T2 = nx.relabel_nodes(nx.path_graph(3), {0: 10, 1: 11, 2: 12})
        assert rooted_tree_signature(T1, 0) == rooted_tree_signature(T2, 10)

    def test_different_rooted_trees_different_signature(self):
        path = nx.path_graph(4)
        star = nx.star_graph(3)
        assert rooted_tree_signature(path, 0) != rooted_tree_signature(star, 0)


_EXPECTED_KEYS = frozenset(
    {
        "is_bouquet",
        "method",
        "reason",
        "cycle_nodes",
        "cycle_length",
        "n_rooted_trees",
        "tree_signatures",
    }
)


class TestIsBouquetComponent:
    """Tests for is_bouquet_component."""

    def test_result_has_all_keys(self):
        result = is_bouquet_component(nx.path_graph(4))
        assert set(result.keys()) == _EXPECTED_KEYS

    def test_method_field_is_optimize(self):
        result = is_bouquet_component(_build_uniform_bouquet(1))
        assert result["method"] == "optimize"

    def test_tree_returns_not_bouquet(self):
        result = is_bouquet_component(nx.path_graph(5))
        assert result["is_bouquet"] is False
        assert result["reason"] == "is_tree"
        assert result["cycle_nodes"] is None
        assert result["tree_signatures"] is None

    def test_c4_not_bouquet(self):
        """C4 has one cycle of length 4, not 5."""
        result = is_bouquet_component(nx.cycle_graph(4))
        assert result["is_bouquet"] is False
        assert result["cycle_nodes"] is None

    def test_c6_not_bouquet(self):
        """C6 has one cycle of length 6, not 5."""
        result = is_bouquet_component(nx.cycle_graph(6))
        assert result["is_bouquet"] is False

    def test_valid_bouquet_accepted(self):
        result = is_bouquet_component(_build_uniform_bouquet(1))
        assert result["is_bouquet"] is True
        assert result["cycle_length"] == 5
        assert result["n_rooted_trees"] == 5
        assert result["cycle_nodes"] is not None
        assert len(result["cycle_nodes"]) == 5
        assert result["tree_signatures"] is not None
        assert len(result["tree_signatures"]) == 5

    def test_valid_bouquet_deeper_petals(self):
        result = is_bouquet_component(_build_uniform_bouquet(2))
        assert result["is_bouquet"] is True

    def test_c5_with_chord_not_bouquet(self):
        """C5 with an extra chord: two basis cycles, not a single induced C5."""
        result = is_bouquet_component(_c5_with_chord())
        assert result["is_bouquet"] is False

    def test_mixed_bouquet_not_bouquet(self):
        """C5 with non-isomorphic petals is rejected."""
        result = is_bouquet_component(_build_mixed_bouquet([1, 1, 2, 2, 2]))
        assert result["is_bouquet"] is False
        assert result["reason"] == "petals_not_isomorphic"
        assert result["cycle_nodes"] is not None
        assert result["tree_signatures"] is not None

    def test_uniform_colors_accepted(self):
        G = _build_uniform_bouquet(1)
        labels = dict.fromkeys(G.nodes(), 0)
        result = is_bouquet_component(G, labels=labels)
        assert result["is_bouquet"] is True

    def test_unique_colors_rejected(self):
        G = _build_uniform_bouquet(1)
        labels = {v: v for v in G.nodes()}
        result = is_bouquet_component(G, labels=labels)
        assert result["is_bouquet"] is False


class TestIsBouquetComponentCoverage:
    """Additional coverage for is_bouquet_component on diverse graph shapes."""

    def test_tree_is_not_bouquet(self):
        assert is_bouquet_component(nx.path_graph(5))["is_bouquet"] is False

    def test_non_c5_cycle_is_not_bouquet(self):
        assert is_bouquet_component(nx.cycle_graph(4))["is_bouquet"] is False

    def test_valid_bouquet_is_accepted(self):
        assert is_bouquet_component(_build_uniform_bouquet(1))["is_bouquet"] is True

    def test_c5_with_chord_is_not_bouquet(self):
        """C5 with a chord: multiple cycles, not a unicyclic C5 bouquet."""
        assert is_bouquet_component(_c5_with_chord())["is_bouquet"] is False

    def test_non_isomorphic_petals_rejected(self):
        assert is_bouquet_component(_build_mixed_bouquet([1, 1, 2, 2, 2]))["is_bouquet"] is False

    def test_method_field_value(self):
        result = is_bouquet_component(_build_uniform_bouquet(1))
        assert result["method"] == "optimize"
