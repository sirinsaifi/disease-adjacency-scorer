"""
score.py

Combines five overlap scores into a single adjacency score using a
weighted formula.

Dimensions (all in [0, 1]):
  target_score      -- Jaccard similarity on associated gene/protein sets
  pathway_score     -- Jaccard similarity on Reactome pathway sets
  genetics_score    -- Weighted overlap on GWAS/genetics evidence
  ot_indirect_score -- EFO ontology-proximity score (from fetch_ot_indirect)
  ae_score          -- Jaccard similarity on FAERS adverse-event term sets

Weights live in config/settings.yaml (scoring.weights). Pass an explicit
`weights` dict to override (e.g. in tests or sensitivity analysis).
See docs/adjacency_definition.md for the justification behind the defaults.
"""

from typing import Dict, Optional

from disease_adjacency.config import load_config

EXPECTED_WEIGHT_KEYS = {
    "target_score",
    "pathway_score",
    "genetics_score",
    "ot_indirect_score",
    "ae_score",
}


def get_configured_weights() -> Dict[str, float]:
    """Read the current scoring weights from config/settings.yaml."""
    return load_config()["scoring"]["weights"]


def validate_weights(weights: Dict[str, float]) -> None:
    """
    Weights must contain exactly the five expected keys and sum to 1.0
    (within floating-point tolerance). Raises ValueError if invalid.
    """
    if set(weights.keys()) != EXPECTED_WEIGHT_KEYS:
        raise ValueError(
            f"Weights must have exactly these keys: {EXPECTED_WEIGHT_KEYS}. "
            f"Got: {set(weights.keys())}"
        )

    total = sum(weights.values())
    if not (0.999 <= total <= 1.001):
        raise ValueError(f"Weights must sum to 1.0, got {total:.4f}")


def compute_adjacency_score(
    overlap_scores: Dict[str, float],
    weights: Optional[Dict[str, float]] = None,
) -> float:
    """
    Combine five overlap scores into one adjacency score using a weighted sum.

    overlap_scores: dict with keys:
        target_score, pathway_score, genetics_score,
        ot_indirect_score, ae_score   (each between 0 and 1)
    weights: optional override; defaults to config/settings.yaml values.

    Returns a float in [0, 1] rounded to 4 decimal places.
    """
    if weights is None:
        weights = get_configured_weights()

    validate_weights(weights)

    if not EXPECTED_WEIGHT_KEYS.issubset(overlap_scores.keys()):
        missing = EXPECTED_WEIGHT_KEYS - set(overlap_scores.keys())
        raise ValueError(f"overlap_scores is missing keys: {missing}")

    for key in EXPECTED_WEIGHT_KEYS:
        value = overlap_scores[key]
        if not (0.0 <= value <= 1.0):
            raise ValueError(f"{key} must be between 0 and 1, got {value}")

    score = sum(overlap_scores[key] * weights[key] for key in EXPECTED_WEIGHT_KEYS)

    # Guard against floating-point drift outside [0, 1]
    return round(max(0.0, min(1.0, score)), 4)
