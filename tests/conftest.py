"""
Shared fixtures and constants for the wl_identifiability test suite.

Module-level constants (VALID_SMILES, INVALID_SMILES) are intentionally kept
alongside their fixture counterparts so that @pytest.mark.parametrize
decorators can reference them directly without needing an import.
"""

import pytest
import networkx as nx


VALID_SMILES = "CCO"
INVALID_SMILES = "not_smiles"


@pytest.fixture
def valid_smiles() -> str:
    """SMILES string that always parses successfully (ethanol)."""
    return VALID_SMILES


@pytest.fixture
def invalid_smiles() -> str:
    """SMILES string that always fails RDKit parsing."""
    return INVALID_SMILES




@pytest.fixture
def bouquet_graph() -> nx.Graph:
    """
    Standard test bouquet: C₅ (nodes 0–4) with one pendant leaf per cycle
    node (nodes 5–9). All five petals are isomorphic single-edge trees.
    """
    G = nx.cycle_graph(5)
    for i in range(5):
        G.add_edge(i, 5 + i)
    return G
