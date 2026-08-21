"""
compare.py

Computes overlap between two disease fingerprints (see build_fingerprint.py).

Each fingerprint has four sets: targets, pathways, genetics, adverse_events.
This module turns those into four Jaccard-similarity scores plus the actual
shared items (for the evidence record in evidence.py). The fifth scoring
dimension, ot_indirect_score, is pair-level rather than fingerprint-based and
is computed separately by fetch_ot_indirect.py, then merged into these scores
by pipeline.py before compute_adjacency_score() combines all five.

See docs/adjacency_definition.md for the full scoring model.
"""

from typing import Dict, Set

REQUIRED_FINGERPRINT_KEYS = {"targets", "pathways", "genetics", "adverse_events"}


def jaccard_similarity(set_a: Set[str], set_b: Set[str]) -> float:
    """
    |A ∩ B| / |A ∪ B|. Returns 0.0 when both sets are empty (no signal,
    not a division-by-zero error).
    """
    union = set_a | set_b
    if not union:
        return 0.0
    return len(set_a & set_b) / len(union)


def _validate_fingerprint(fp: dict, label: str) -> None:
    missing = REQUIRED_FINGERPRINT_KEYS - set(fp.keys())
    if missing:
        raise ValueError(f"{label} is missing fingerprint keys: {missing}")


def compare_fingerprints(fp_a: dict, fp_b: dict) -> Dict[str, float]:
    """
    Compare two fingerprints and return the four fingerprint-based overlap
    scores: target_score, pathway_score, genetics_score, ae_score.

    Raises ValueError if either fingerprint is missing a required key.
    """
    _validate_fingerprint(fp_a, "fp_a")
    _validate_fingerprint(fp_b, "fp_b")

    return {
        "target_score": jaccard_similarity(fp_a["targets"], fp_b["targets"]),
        "pathway_score": jaccard_similarity(fp_a["pathways"], fp_b["pathways"]),
        "genetics_score": jaccard_similarity(fp_a["genetics"], fp_b["genetics"]),
        "ae_score": jaccard_similarity(fp_a["adverse_events"], fp_b["adverse_events"]),
    }


def get_shared_evidence(fp_a: dict, fp_b: dict) -> Dict[str, Set[str]]:
    """
    Return the actual shared items behind each overlap score (not just the
    counts), so every adjacency score traces to specific evidence rows.
    """
    _validate_fingerprint(fp_a, "fp_a")
    _validate_fingerprint(fp_b, "fp_b")

    return {
        "shared_targets": fp_a["targets"] & fp_b["targets"],
        "shared_pathways": fp_a["pathways"] & fp_b["pathways"],
        "shared_genetics": fp_a["genetics"] & fp_b["genetics"],
        "shared_adverse_events": fp_a["adverse_events"] & fp_b["adverse_events"],
    }
