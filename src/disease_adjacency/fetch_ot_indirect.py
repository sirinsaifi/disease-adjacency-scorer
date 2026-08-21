"""
fetch_ot_indirect.py

Computes the "OT indirect" / ontology-proximity score (dimension 4): how
close two diseases are in the EFO hierarchy, queried directly from the
EMBL-EBI Ontology Lookup Service (OLS4) API rather than re-querying Open
Targets. See docs/adjacency_definition.md for the full rationale.

Score bands:
  1.0  = same disease (identical EFO IDs)
  0.67 = direct parent/child relationship in EFO
  0.49 = share at least one EFO ancestor (siblings/cousins)
  0.0  = no shared ancestors, or the OLS4 lookup failed (conservative)

OLS4 terms are addressed by IRI, not by short-form ID, and the endpoint
requires the IRI to be percent-encoded twice. Most EFO-native terms live
under http://www.ebi.ac.uk/efo/, while terms EFO imports from other OBO
ontologies (MONDO, HP, Orphanet, ...) live under the shared OBO PURL
namespace http://purl.obolibrary.org/obo/ -- this is a stable, documented
OLS/OBO convention, not a guess specific to this dataset.

Every pair's evidence (shared-ancestor count, direct-relationship flag,
score) is cached to data/raw/ot_indirect_{anchor}_{candidate}.json so the
score is traceable and repeated runs don't re-hit OLS4.
"""

import json
import os
import time
import urllib.parse
from typing import Set

import requests

from disease_adjacency.config import load_config

# Not tunable policy knobs -- these define the scoring scale itself (identical
# disease = 1.0, no relationship found = 0.0), so they stay as code constants
# rather than config values. direct_relationship_score and
# shared_ancestor_score_ceiling (the tunable bands in between) live in
# config/settings.yaml under "ols".
SAME_DISEASE_SCORE = 1.0
NO_RELATIONSHIP_SCORE = 0.0


def _efo_iri(efo_id: str) -> str:
    if efo_id.startswith("EFO_"):
        return f"http://www.ebi.ac.uk/efo/{efo_id}"
    return f"http://purl.obolibrary.org/obo/{efo_id}"


def _double_encode(iri: str) -> str:
    # OLS4's terms-by-IRI endpoint requires the IRI to be percent-encoded twice.
    return urllib.parse.quote(urllib.parse.quote(iri, safe=""), safe="")


def _ols_related_ids(
    efo_id: str,
    relation: str,
    term_base_url: str,
    max_retries: int,
    retry_backoff_seconds: int,
    page_size: int,
) -> Set[str]:
    """
    Fetch related term short-forms (e.g. "MONDO_0005138") from OLS4.
    relation: "parents" or "ancestors".
    Returns an empty set if the term has none, or if the lookup fails after
    retries -- callers treat "no relationship found" the same as "none exist".
    """
    encoded_iri = _double_encode(_efo_iri(efo_id))
    url = f"{term_base_url}/{encoded_iri}/{relation}"

    last_error = None
    for attempt in range(1, max_retries + 1):
        try:
            response = requests.get(url, params={"size": page_size}, timeout=20)
            if response.status_code == 404:
                return set()  # term has no parents/ancestors recorded (e.g. a root term)
            response.raise_for_status()
            terms = response.json().get("_embedded", {}).get("terms", []) or []
            return {t["short_form"] for t in terms if t.get("short_form")}
        except requests.RequestException as e:
            last_error = e
            if attempt < max_retries:
                time.sleep(attempt * retry_backoff_seconds)

    print(f"Warning: OLS4 {relation} lookup failed for {efo_id}: {last_error}")
    return set()


def _save_evidence(cache_path: str, evidence: dict) -> None:
    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    with open(cache_path, "w") as f:
        json.dump(evidence, f, indent=2)


def fetch_ot_indirect_score(anchor_efo_id: str, efo_id: str, use_cache: bool = True) -> float:
    """
    Compute the ontology-proximity score between two EFO-resolved diseases.
    Checks the local per-pair cache first if use_cache is True.
    """
    config = load_config()
    raw_data_dir = config["open_targets"]["raw_data_dir"]
    ols_config = config["ols"]
    cache_path = os.path.join(raw_data_dir, f"ot_indirect_{anchor_efo_id}_{efo_id}.json")

    if use_cache and os.path.exists(cache_path):
        with open(cache_path, "r") as f:
            return json.load(f)["score"]

    if anchor_efo_id == efo_id:
        evidence = {
            "anchor_efo_id": anchor_efo_id,
            "candidate_efo_id": efo_id,
            "direct_parent_child": False,
            "shared_ancestor_count": None,
            "score": SAME_DISEASE_SCORE,
        }
        _save_evidence(cache_path, evidence)
        return SAME_DISEASE_SCORE

    term_base_url = ols_config["term_base_url"]
    max_retries = ols_config["max_retries"]
    retry_backoff_seconds = ols_config["retry_backoff_seconds"]
    page_size = ols_config["ancestor_page_size"]

    anchor_parents = _ols_related_ids(anchor_efo_id, "parents", term_base_url, max_retries, retry_backoff_seconds, page_size)
    candidate_parents = _ols_related_ids(efo_id, "parents", term_base_url, max_retries, retry_backoff_seconds, page_size)
    direct_parent_child = efo_id in anchor_parents or anchor_efo_id in candidate_parents

    if direct_parent_child:
        score = ols_config["direct_relationship_score"]
        shared_ancestor_count = None
    else:
        anchor_ancestors = _ols_related_ids(anchor_efo_id, "ancestors", term_base_url, max_retries, retry_backoff_seconds, page_size)
        candidate_ancestors = _ols_related_ids(efo_id, "ancestors", term_base_url, max_retries, retry_backoff_seconds, page_size)
        shared_ancestor_count = len(anchor_ancestors & candidate_ancestors)
        score = ols_config["shared_ancestor_score_ceiling"] if shared_ancestor_count > 0 else NO_RELATIONSHIP_SCORE

    evidence = {
        "anchor_efo_id": anchor_efo_id,
        "candidate_efo_id": efo_id,
        "direct_parent_child": direct_parent_child,
        "shared_ancestor_count": shared_ancestor_count,
        "score": score,
    }
    _save_evidence(cache_path, evidence)
    return score
