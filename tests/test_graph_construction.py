"""Tests for graph_construction.py — RDKit molecule → NetworkX graph."""

import networkx as nx
import pytest
from rdkit import Chem

from wl_identifiability.graph_construction import convert_rdkit_molecule_to_nx_graph


def mol(smiles: str) -> Chem.Mol:
    m = Chem.MolFromSmiles(smiles)
    assert m is not None, f"Invalid SMILES: {smiles}"
    return m


class TestConvertRdkitMoleculeToNxGraph:
    def test_returns_nx_graph(self):
        G = convert_rdkit_molecule_to_nx_graph(mol("CC"))
        assert isinstance(G, nx.Graph)

    def test_graph_is_undirected(self):
        G = convert_rdkit_molecule_to_nx_graph(mol("CCC"))
        assert not G.is_directed()

    def test_propane_is_connected(self):
        G = convert_rdkit_molecule_to_nx_graph(mol("CCC"))
        assert nx.is_connected(G)

    @pytest.mark.parametrize(
        "smiles,n_nodes,n_edges",
        [
            pytest.param("CC", 2, 1, id="ethane"),
            pytest.param("c1ccccc1", 6, 6, id="benzene"),
            pytest.param("CCC", 3, 2, id="propane"),
        ],
    )
    def test_node_and_edge_counts(self, smiles, n_nodes, n_edges):
        G = convert_rdkit_molecule_to_nx_graph(mol(smiles))
        assert G.number_of_nodes() == n_nodes
        assert G.number_of_edges() == n_edges

    @pytest.mark.parametrize("attr", ["atomic_num", "is_aromatic", "formal_charge"])
    def test_node_has_attribute(self, attr):
        G = convert_rdkit_molecule_to_nx_graph(mol("CO"))
        for _, data in G.nodes(data=True):
            assert attr in data

    @pytest.mark.parametrize("attr", ["bond_type", "is_aromatic"])
    def test_edge_has_attribute(self, attr):
        G = convert_rdkit_molecule_to_nx_graph(mol("CC"))
        for _, _, data in G.edges(data=True):
            assert attr in data

    @pytest.mark.parametrize(
        "smiles,expected_num",
        [
            pytest.param("CC", 6, id="carbon"),
            pytest.param("O", 8, id="oxygen"),
        ],
    )
    def test_atomic_num(self, smiles, expected_num):
        G = convert_rdkit_molecule_to_nx_graph(mol(smiles))
        for _, data in G.nodes(data=True):
            assert data["atomic_num"] == expected_num
