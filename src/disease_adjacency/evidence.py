"""
evidence.py

Formats the shared-evidence rows (from compare.get_shared_evidence) plus
the per-dimension overlap scores into a clean, audit-friendly record
attached to each adjacency score.

Every adjacency score in the output must trace to the shared-biology rows
in this record. A score with no attached evidence rows is invalid output
under the project's traceability requirement.
"""

from typing import Dict, Set


def build_evidence_record(
    disease_a: str,
    disease_b: str,
    shared_evidence: Dict[str, Set[str]],
    overlap_scores: Dict[str, float],
    adjacency_score: float,
) -> dict:
    """
    Combine a disease pair, their shared evidence sets, per-dimension
    overlap scores, and the final combined score into one inspectable record.

    Parameters
    ----------
    disease_a, disease_b : display names of the two diseases
    shared_evidence      : from compare.get_shared_evidence() -- the actual
                           shared items per dimension
    overlap_scores       : the five per-dimension scores before weighting
                           (target_score, pathway_score, genetics_score,
                            ot_indirect_score, ae_score)
    adjacency_score      : the final weighted combined score

    Returns a dict that is written directly to output/adjacency_scores.json.
    """
    return {
        "disease_a": disease_a,
        "disease_b": disease_b,
        "adjacency_score": adjacency_score,
        # Per-dimension breakdown so every score is decomposable
        "components": {
            "target_score":       round(overlap_scores.get("target_score",       0.0), 4),
            "pathway_score":      round(overlap_scores.get("pathway_score",      0.0), 4),
            "genetics_score":     round(overlap_scores.get("genetics_score",     0.0), 4),
            "ot_indirect_score":  round(overlap_scores.get("ot_indirect_score",  0.0), 4),
            "ae_score":           round(overlap_scores.get("ae_score",           0.0), 4),
        },
        # Shared-biology rows that justify each component score
        "shared_targets":        sorted(shared_evidence.get("shared_targets",        set())),
        "shared_pathways":       sorted(shared_evidence.get("shared_pathways",       set())),
        "shared_genetics":       sorted(shared_evidence.get("shared_genetics",       set())),
        "shared_adverse_events": sorted(shared_evidence.get("shared_adverse_events", set())),
    }
