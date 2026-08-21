"""
pipeline.py

Orchestrates the full adjacency scoring pipeline:
  resolve disease names
    → fetch OT associations + AE fingerprint
    → build fingerprint (targets / pathways / genetics / adverse_events)
    → fetch OT indirect (ontology-proximity) per pair
    → compare fingerprints → combine 5 scores → build evidence record
    → save output

Disease names come from config/diseases.csv (one anchor + N candidates).
Each name is resolved via the OLS API in resolve_disease.py. Adding a new
disease means adding a row to that CSV -- no EFO ID lookup required by hand.

If a name doesn't resolve with high confidence the pipeline stops and
reports candidates for human review rather than guessing.

Run with: poetry run python src/disease_adjacency/pipeline.py
"""

import csv
import json
import os

from disease_adjacency.config import load_config
from disease_adjacency.resolve_disease import (
    resolve_with_cache,
    suggest_candidates_for_confirmation,
)
from disease_adjacency.fetch_opentargets import fetch_disease_associations
from disease_adjacency.fetch_adverse_events import fetch_ae_fingerprint
from disease_adjacency.fetch_ot_indirect import fetch_ot_indirect_score
from disease_adjacency.build_fingerprint import build_fingerprint
from disease_adjacency.compare import compare_fingerprints, get_shared_evidence
from disease_adjacency.score import compute_adjacency_score
from disease_adjacency.evidence import build_evidence_record


class DiseaseResolutionError(Exception):
    """Raised when a disease name can't be auto-resolved with confidence."""
    def __init__(self, name, candidates):
        self.name = name
        self.candidates = candidates
        candidate_lines = "\n".join(
            f"  - {c.label} ({c.efo_id}) confidence={c.confidence} [{c.match_type}]"
            for c in candidates
        )
        super().__init__(
            f"Could not confidently resolve disease name '{name}'.\n"
            f"Top candidates:\n{candidate_lines}\n"
            f"Confirm the correct one and add it via resolve_disease.confirm_and_cache(), "
            f"or correct the spelling in config/diseases.csv."
        )


def load_disease_list() -> tuple:
    """
    Reads config/diseases.csv. Returns (anchor_name, [candidate_names]).
    """
    disease_list_path = load_config()["paths"]["disease_list"]

    anchor_name = None
    candidates = []

    with open(disease_list_path, "r") as f:
        for row in csv.DictReader(f):
            role = row["role"].strip().lower()
            name = row["disease_name"].strip()
            if role == "anchor":
                anchor_name = name
            elif role == "candidate":
                candidates.append(name)

    if not anchor_name:
        raise ValueError(f"No anchor disease found in {disease_list_path}")

    return anchor_name, candidates


def resolve_or_raise(name: str) -> str:
    """
    Resolve a disease name to an EFO ID. Raises DiseaseResolutionError with
    ranked candidates if it can't be resolved with confidence.
    """
    match = resolve_with_cache(name)
    if match is None:
        candidates = suggest_candidates_for_confirmation(name)
        raise DiseaseResolutionError(name, candidates)
    return match.efo_id


def build_disease_fingerprint(name: str) -> tuple:
    """
    Resolve a disease name and build its biology fingerprint (targets,
    pathways, genetics, adverse_events). Returns (efo_id, fingerprint).
    Used for both the anchor and each candidate -- the anchor's fingerprint
    is built once and reused across every pair.
    """
    efo_id = resolve_or_raise(name)
    raw = fetch_disease_associations(efo_id)
    ae = fetch_ae_fingerprint(efo_id)
    fingerprint = build_fingerprint(raw, adverse_events=ae)
    return efo_id, fingerprint


def score_pair(anchor_name: str, anchor_efo_id: str, anchor_fp: dict, candidate_name: str) -> dict:
    """
    Score one candidate disease against an already-built anchor fingerprint.
    Resolves the candidate, fetches/builds its fingerprint, computes all five
    overlap scores, and returns a full evidence record (same shape as one
    entry in output/adjacency_scores.json). Raises DiseaseResolutionError if
    the candidate name can't be confidently resolved.
    """
    candidate_efo_id, candidate_fp = build_disease_fingerprint(candidate_name)

    # Compute the four fingerprint-based Jaccard scores
    overlap_scores = compare_fingerprints(anchor_fp, candidate_fp)

    # Add the fifth score: EFO ontology-proximity (pair-level, not fingerprint-based)
    overlap_scores["ot_indirect_score"] = fetch_ot_indirect_score(anchor_efo_id, candidate_efo_id)

    # Combine into final adjacency score
    adjacency_score = compute_adjacency_score(overlap_scores)

    # Build inspectable evidence record
    shared_evidence = get_shared_evidence(anchor_fp, candidate_fp)
    return build_evidence_record(
        disease_a=anchor_name,
        disease_b=candidate_name,
        shared_evidence=shared_evidence,
        overlap_scores=overlap_scores,
        adjacency_score=adjacency_score,
    )


def run_pipeline() -> list:
    anchor_name, candidate_names = load_disease_list()

    # --- Resolve and fetch anchor ---
    anchor_efo_id, anchor_fp = build_disease_fingerprint(anchor_name)
    print(f"Anchor: {anchor_name} ({anchor_efo_id})")

    print(f"  Anchor fingerprint: {len(anchor_fp['targets'])} targets, "
          f"{len(anchor_fp['pathways'])} pathways, "
          f"{len(anchor_fp['genetics'])} genetics, "
          f"{len(anchor_fp['adverse_events'])} AE terms")

    results = []

    for name in candidate_names:
        print(f"\nScoring: {name}")
        record = score_pair(anchor_name, anchor_efo_id, anchor_fp, name)
        results.append(record)

        c = record["components"]
        print(f"  adj_score={record['adjacency_score']} | "
              f"target={c['target_score']} | "
              f"pathway={c['pathway_score']} | "
              f"genetics={c['genetics_score']} | "
              f"ot_indirect={c['ot_indirect_score']} | "
              f"ae={c['ae_score']}")

    results.sort(key=lambda r: r["adjacency_score"], reverse=True)
    return results


def save_results(results: list) -> None:
    output_path = load_config()["paths"]["output_scores"]
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(results, f, indent=2)


if __name__ == "__main__":
    try:
        results = run_pipeline()
        save_results(results)
        print(f"\nDone. Scored {len(results)} disease pairs.")
        print(f"Output: {load_config()['paths']['output_scores']}")
    except DiseaseResolutionError as e:
        print(f"\nPipeline stopped -- needs human input:\n{e}")
