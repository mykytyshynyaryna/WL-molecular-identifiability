"""
wl_identifiability — 1-WL identifiability analysis for molecular graphs.

Public API:
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

from .wl import (
    compute_wl_coloring,
    compute_fixed_wl_steps_from_topology,
    compute_wl_with_fixed_steps,
    compute_wl_stabilization_steps,
    compute_label_histogram,
    normalize_node_labels,
)
from .graph_construction import convert_rdkit_molecule_to_nx_graph
from .flip_graph import build_flip_graph_from_labels
from .skeleton import (
    classify_between_class_relation,
    classify_within_class_structure,
    build_skeleton,
    check_lemma16_conditions,
    REL_EMPTY,
    REL_MATCHING,
    REL_FAN,
    REL_IRREGULAR,
)
from .bouquet import (
    check_bouquet_component,
    analyze_bouquet_forest_structure,
    is_bouquet_component,
    rooted_tree_signature,
)
from .visualization import draw_graph_wl_coloring, inspect_wl_behavior
from .experiments import (
    estimate_global_wl_steps,
    estimate_fixed_wl_steps_from_dataframe,
    analyze_single_molecule,
    run_molecule_analysis_pipeline,
)

__all__ = [
    "compute_wl_coloring",
    "compute_fixed_wl_steps_from_topology",
    "compute_wl_with_fixed_steps",
    "compute_wl_stabilization_steps",
    "convert_rdkit_molecule_to_nx_graph",
    "build_flip_graph_from_labels",
    "classify_between_class_relation",
    "classify_within_class_structure",
    "build_skeleton",
    "check_lemma16_conditions",
    "REL_EMPTY",
    "REL_MATCHING",
    "REL_FAN",
    "REL_IRREGULAR",
    "check_bouquet_component",
    "analyze_bouquet_forest_structure",
    "is_bouquet_component",
    "rooted_tree_signature",
    "compute_label_histogram",
    "normalize_node_labels",
    "draw_graph_wl_coloring",
    "inspect_wl_behavior",
    "estimate_global_wl_steps",
    "estimate_fixed_wl_steps_from_dataframe",
    "analyze_single_molecule",
    "run_molecule_analysis_pipeline",
]
