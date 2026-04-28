from __future__ import annotations

from ._imports import np, nx, Chem


def convert_rdkit_molecule_to_nx_graph(mol: Chem.Mol) -> nx.Graph:
    """
    Convert an RDKit molecule into a NetworkX graph using adjacency matrix.
    """
    n = mol.GetNumAtoms()

    A = np.zeros((n, n), dtype=np.uint8)
    edge_attrs = {}

    for bond in mol.GetBonds():
        u = bond.GetBeginAtomIdx()
        v = bond.GetEndAtomIdx()

        A[u, v] = 1
        A[v, u] = 1

        edge_attrs[(u, v)] = {
            "bond_type": str(bond.GetBondType()),
            "is_aromatic": bond.GetIsAromatic(),
        }
        edge_attrs[(v, u)] = edge_attrs[(u, v)]

    G = nx.from_numpy_array(A)

    for atom in mol.GetAtoms():
        i = atom.GetIdx()
        G.nodes[i].update(
            {
                "atomic_num": atom.GetAtomicNum(),
                "is_aromatic": atom.GetIsAromatic(),
                "formal_charge": atom.GetFormalCharge(),
            }
        )

    for u, v in G.edges():
        if (u, v) in edge_attrs:
            G[u][v].update(edge_attrs[(u, v)])

    return G
