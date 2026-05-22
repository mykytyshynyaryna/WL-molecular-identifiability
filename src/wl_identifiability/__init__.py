"""
wl_identifiability — 1-WL identifiability analysis for molecular graphs.

Public API:
  Quick checks (simplest entry points)
    is_smi_identifiable   — SMILES string  → bool
    is_mol_identifiable   — RDKit Mol      → bool
    is_graph_identifiable — NetworkX Graph → bool

  WL coloring
    compute_wl_coloring
    compute_fixed_wl_steps_from_topology
    compute_wl_with_fixed_steps
    compute_wl_stabilization_steps
    compute_label_histogram

  Graph construction
    convert_rdkit_molecule_to_nx_graph

  Flip graph
    build_flip_graph_from_labels

  Skeleton (Kiefer, Notation 15 / Section 4 / Lemma 16)
    classify_between_class_relation
    classify_within_class_structure
    build_skeleton
    check_lemma16_conditions
    REL_EMPTY, REL_MATCHING, REL_FAN, REL_IRREGULAR

  Bouquet detection
    check_bouquet_component
    analyze_bouquet_forest_structure

  Visualization
    draw_graph_wl_coloring
    inspect_wl_behavior

  Pipeline
    estimate_global_wl_steps
    estimate_fixed_wl_steps_from_dataframe
    analyze_single_molecule
    run_molecule_analysis_pipeline
"""

from .bouquet import (
    analyze_bouquet_forest_structure,
    check_bouquet_component,
    is_bouquet_component,
    rooted_tree_signature,
)
from .experiments import (
    analyze_single_molecule,
    estimate_fixed_wl_steps_from_dataframe,
    estimate_global_wl_steps,
    is_graph_identifiable,
    is_mol_identifiable,
    is_smi_identifiable,
    run_molecule_analysis_pipeline,
)
from .flip_graph import build_flip_graph_from_labels
from .graph_construction import convert_rdkit_molecule_to_nx_graph
from .skeleton import (
    REL_EMPTY,
    REL_FAN,
    REL_IRREGULAR,
    REL_MATCHING,
    build_skeleton,
    check_lemma16_conditions,
    classify_between_class_relation,
    classify_within_class_structure,
)
from .visualization import draw_graph_wl_coloring, inspect_wl_behavior
from .wl import (
    compute_fixed_wl_steps_from_topology,
    compute_label_histogram,
    compute_wl_coloring,
    compute_wl_stabilization_steps,
    compute_wl_with_fixed_steps,
    normalize_node_labels,
)

__all__ = [
    # Quick checks — simplest entry points
    "is_smi_identifiable",
    "is_mol_identifiable",
    "is_graph_identifiable",
    # Lower-level API
    "REL_EMPTY",
    "REL_FAN",
    "REL_IRREGULAR",
    "REL_MATCHING",
    "analyze_bouquet_forest_structure",
    "analyze_single_molecule",
    "build_flip_graph_from_labels",
    "build_skeleton",
    "check_bouquet_component",
    "check_lemma16_conditions",
    "classify_between_class_relation",
    "classify_within_class_structure",
    "compute_fixed_wl_steps_from_topology",
    "compute_label_histogram",
    "compute_wl_coloring",
    "compute_wl_stabilization_steps",
    "compute_wl_with_fixed_steps",
    "convert_rdkit_molecule_to_nx_graph",
    "draw_graph_wl_coloring",
    "estimate_fixed_wl_steps_from_dataframe",
    "estimate_global_wl_steps",
    "inspect_wl_behavior",
    "is_bouquet_component",
    "normalize_node_labels",
    "rooted_tree_signature",
    "run_molecule_analysis_pipeline",
]
