"""Tests for experiments.py — single-molecule analysis and pipeline orchestration."""

import networkx as nx
import polars as pl
import pytest

from wl_identifiability.experiments import (
    estimate_global_wl_steps,
    _compute_wl_summary,
    _compute_flip_graph_summary,
    _compute_color_ratio,
    analyze_single_molecule,
    run_molecule_analysis_pipeline,
)





def _make_wl_result(n_labels=3, converged=True, iterations=2):
    labels = {i: i % n_labels for i in range(6)}
    return {"labels": labels, "converged": converged, "iterations": iterations}




class TestEstimateGlobalWlSteps:
    @pytest.mark.parametrize("key", ["K_max", "K_p95", "all_iters"])
    def test_result_has_key(self, key):
        assert key in estimate_global_wl_steps([nx.path_graph(4)], cap=10)

    def test_empty_graphs_list_returns_zeros(self):
        result = estimate_global_wl_steps([], cap=10)
        assert result["K_max"] == 0
        assert result["K_p95"] == 0
        assert result["all_iters"] == []

    def test_k_max_is_int(self):
        assert isinstance(
            estimate_global_wl_steps([nx.path_graph(6)], cap=20)["K_max"], int
        )

    def test_cap_respected(self):
        assert estimate_global_wl_steps([nx.path_graph(30)], cap=3)["K_max"] <= 3

    def test_k_p95_lte_k_max(self):
        result = estimate_global_wl_steps(
            [nx.path_graph(i) for i in range(3, 10)], cap=20
        )
        assert result["K_p95"] <= result["K_max"]




class TestComputeWlSummary:
    @pytest.mark.parametrize("prefix", ["top", "atom", "x"])
    def test_returns_prefixed_keys(self, prefix):
        summary = _compute_wl_summary(_make_wl_result(), prefix)
        assert f"wl_converged_{prefix}" in summary
        assert f"wl_iters_{prefix}" in summary
        assert f"n_colors_{prefix}" in summary

    def test_n_colors_counts_unique_labels(self):
        wl = {"labels": {0: 0, 1: 1, 2: 0, 3: 2}, "converged": True, "iterations": 1}
        assert _compute_wl_summary(wl, "x")["n_colors_x"] == 3

    def test_empty_labels_gives_zero_colors(self):
        wl = {"labels": {}, "converged": True, "iterations": 0}
        assert _compute_wl_summary(wl, "top")["n_colors_top"] == 0




class TestComputeFlipGraphSummary:
    @pytest.mark.parametrize(
        "key",
        [
            "within_copy",
            "within_flip",
            "between_copy",
            "between_flip",
            "n_flipped_edges",
        ],
    )
    def test_result_has_key(self, key):
        info = {
            "within_copy": 1,
            "within_flip": 2,
            "between_copy": 3,
            "between_flip": 4,
        }
        assert key in _compute_flip_graph_summary(info)

    def test_n_flipped_edges_is_sum_of_flips(self):
        info = {
            "within_copy": 0,
            "within_flip": 3,
            "between_copy": 0,
            "between_flip": 2,
        }
        assert _compute_flip_graph_summary(info)["n_flipped_edges"] == 5

    def test_missing_keys_default_to_zero(self):
        assert _compute_flip_graph_summary({})["n_flipped_edges"] == 0




class TestComputeColorRatio:
    @pytest.mark.parametrize(
        "top_labels,atom_labels,expected",
        [
            pytest.param(
                {0: 0, 1: 0, 2: 0}, {0: 0, 1: 1, 2: 2}, pytest.approx(3.0), id="ratio_3"
            ),
            pytest.param({0: 0, 1: 1}, {0: 0, 1: 1}, pytest.approx(1.0), id="equal"),
            pytest.param({}, {0: 0}, None, id="zero_top"),
        ],
    )
    def test_color_ratio(self, top_labels, atom_labels, expected):
        assert (
            _compute_color_ratio({"labels": top_labels}, {"labels": atom_labels})
            == expected
        )




class TestAnalyzeSingleMolecule:
    def test_valid_smiles_returns_ok(self, valid_smiles):
        result = analyze_single_molecule(valid_smiles, "zinc_001", fixed_wl_steps=3)
        assert result["ok"] is True
        assert result["stage"] == "done"

    def test_invalid_smiles_returns_parse_error(self, invalid_smiles):
        result = analyze_single_molecule(invalid_smiles, "zinc_bad", fixed_wl_steps=3)
        assert result["ok"] is False
        assert result["stage"] == "parse"
        assert result["reason"] == "rdkit_parse_failed"

    def test_parse_failure_warns_to_stderr(self, invalid_smiles, capsys):
        analyze_single_molecule(invalid_smiles, "zinc_bad", fixed_wl_steps=3)
        captured = capsys.readouterr()
        assert "WARN" in captured.err
        assert "zinc_bad" in captured.err

    def test_zero_wl_steps_does_not_crash(self, valid_smiles):
        assert "ok" in analyze_single_molecule(valid_smiles, "z1", fixed_wl_steps=0)

    def test_larger_molecule(self):
        result = analyze_single_molecule(
            "CC(=O)Oc1ccccc1C(=O)O", "zinc_asp", fixed_wl_steps=5
        )
        assert result["ok"] is True

    @pytest.mark.parametrize(
        "smiles,zinc_id,expected_id,expected_smiles",
        [
            pytest.param("CCO", "ZINC123", "ZINC123", "CCO", id="ethanol"),
        ],
    )
    def test_result_identity_fields(
        self, smiles, zinc_id, expected_id, expected_smiles
    ):
        result = analyze_single_molecule(smiles, zinc_id, fixed_wl_steps=3)
        assert result["zinc_id"] == expected_id
        assert result["smiles"] == expected_smiles

    @pytest.mark.parametrize(
        "field",
        [
            pytest.param("n_colors_top", id="wl_top_colors"),
            pytest.param("n_colors_atom", id="wl_atom_colors"),
            pytest.param("wl_converged_top", id="wl_converged"),
            pytest.param("n_flipped_edges", id="flip_edges"),
            pytest.param("within_copy", id="flip_within_copy"),
            pytest.param("top_bouquet_forest_verdict", id="top_bouquet_forest_verdict"),
            pytest.param("atom_bouquet_forest_verdict", id="atom_bouquet_forest_verdict"),
            pytest.param("top_bouquet_forest_reason", id="top_bouquet_forest_reason"),
            pytest.param("atom_bouquet_forest_reason", id="atom_bouquet_forest_reason"),
            pytest.param("top_has_bouquet_component", id="top_has_bouquet_component"),
            pytest.param("atom_has_bouquet_component", id="atom_has_bouquet_component"),
        ],
    )
    def test_valid_result_has_field(self, field, valid_smiles):
        result = analyze_single_molecule(valid_smiles, "z1", fixed_wl_steps=3)
        assert field in result




class TestRunMoleculeAnalysisPipeline:
    def _df(self, rows):
        return pl.DataFrame(rows)

    def test_returns_dataframe_and_bad_count(self, valid_smiles):
        df = self._df([{"smiles": valid_smiles, "zinc_id": "z1"}])
        result_df, bad = run_molecule_analysis_pipeline(df, fixed_wl_steps=3)
        assert isinstance(result_df, pl.DataFrame)
        assert isinstance(bad, int)

    def test_bad_count_zero_for_valid_molecules(self, valid_smiles):
        df = self._df(
            [
                {"smiles": valid_smiles, "zinc_id": "z1"},
                {"smiles": "C", "zinc_id": "z2"},
            ]
        )
        _, bad = run_molecule_analysis_pipeline(df, fixed_wl_steps=3)
        assert bad == 0

    def test_bad_count_increments_for_invalid_smiles(self, valid_smiles, invalid_smiles):
        df = self._df(
            [
                {"smiles": valid_smiles, "zinc_id": "z1"},
                {"smiles": invalid_smiles, "zinc_id": "z2"},
            ]
        )
        _, bad = run_molecule_analysis_pipeline(df, fixed_wl_steps=3)
        assert bad == 1

    @pytest.mark.parametrize("n_rows", [1, 2, 3])
    def test_result_df_has_one_row_per_input(self, n_rows, valid_smiles):
        df = self._df(
            [{"smiles": valid_smiles, "zinc_id": f"z{i}"} for i in range(n_rows)]
        )
        result_df, _ = run_molecule_analysis_pipeline(df, fixed_wl_steps=3)
        assert len(result_df) == n_rows

    def test_result_df_contains_ok_column(self, valid_smiles):
        df = self._df([{"smiles": valid_smiles, "zinc_id": "z1"}])
        result_df, _ = run_molecule_analysis_pipeline(df, fixed_wl_steps=3)
        assert "ok" in result_df.columns
