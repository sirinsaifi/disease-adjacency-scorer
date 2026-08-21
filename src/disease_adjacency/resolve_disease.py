"""
resolve_disease.py

Resolves a free-text disease name (possibly vague or misspelled) to an EFO
ID, using the EMBL-EBI Ontology Lookup Service (OLS) API.

API: https://www.ebi.ac.uk/ols4/api/search?q=<name>&ontology=efo&rows=<n>
Response shape (confirmed against the live schema):
{
  "response": {
    "docs": [
      {"short_form": "EFO_0003060", "label": "non-small cell lung carcinoma",
       "synonym": ["NSCLC", "non-small cell lung cancer"], "iri": "..."},
      ...
    ]
  }
}

This is deliberately NOT a silent auto-guesser. A wrong disease match
silently poisons every downstream score, so:
- Exact label or exact synonym matches are treated as high confidence.
- Everything else gets a similarity score (0-1) via difflib, and is
  returned as a *candidate to confirm*, not auto-accepted, unless it
  clears AUTO_ACCEPT_THRESHOLD.
- No results at all is reported explicitly, never silently swallowed.

Resolved names are cached to data/efo/resolved_diseases.csv so repeated
pipeline runs don't re-hit the OLS API for the same name.
"""

import csv
import difflib
import os
import time
from dataclasses import dataclass
from typing import List, Optional

import requests

from disease_adjacency.config import load_config

# A match at or above this similarity is auto-accepted without confirmation.
# Exact label/synonym matches always score 1.0 and clear this automatically.
# Sourced from config/settings.yaml (ols.auto_accept_threshold) at import time;
# kept as a module attribute since DiseaseMatch.is_auto_acceptable() and
# existing callers/tests reference it directly.
AUTO_ACCEPT_THRESHOLD = load_config()["ols"]["auto_accept_threshold"]


@dataclass
class DiseaseMatch:
    query: str          # what the user typed
    efo_id: str         # e.g. "EFO_0003060"
    label: str          # official EFO label, e.g. "non-small cell lung carcinoma"
    confidence: float   # 0.0 - 1.0
    match_type: str     # "exact_label" | "exact_synonym" | "fuzzy"

    def is_auto_acceptable(self) -> bool:
        return self.confidence >= AUTO_ACCEPT_THRESHOLD


def _score_candidate(query: str, label: str, synonyms: List[str]) -> tuple:
    """
    Returns (confidence, match_type) for a single OLS candidate against
    the user's query string.
    """
    query_lower = query.strip().lower()
    label_lower = label.strip().lower()

    if query_lower == label_lower:
        return 1.0, "exact_label"

    for syn in synonyms:
        if query_lower == syn.strip().lower():
            return 1.0, "exact_synonym"

    # Fuzzy fallback: best similarity against the label and each synonym
    best_ratio = difflib.SequenceMatcher(None, query_lower, label_lower).ratio()
    for syn in synonyms:
        ratio = difflib.SequenceMatcher(None, query_lower, syn.strip().lower()).ratio()
        best_ratio = max(best_ratio, ratio)

    return round(best_ratio, 4), "fuzzy"


def _query_ols(name: str, search_url: str, max_retries: int, retry_backoff_seconds: int, rows: int) -> list:
    """Query the live OLS API for candidate EFO terms matching a name."""
    params = {"q": name, "ontology": "efo", "rows": rows}

    last_error = None
    for attempt in range(1, max_retries + 1):
        try:
            response = requests.get(search_url, params=params, timeout=20)
            response.raise_for_status()
            docs = response.json().get("response", {}).get("docs", [])
            return docs
        except requests.RequestException as e:
            last_error = e
            if attempt < max_retries:
                time.sleep(attempt * retry_backoff_seconds)

    raise RuntimeError(f"OLS API request failed for '{name}' after {max_retries} attempts: {last_error}")


def resolve_disease_name(name: str) -> List[DiseaseMatch]:
    """
    Resolve a free-text disease name to ranked EFO candidates.

    Returns a list of DiseaseMatch, sorted by confidence descending.
    Returns an empty list if OLS has no matches at all -- this is reported
    explicitly by the caller, never silently treated as "no disease exists".
    """
    ols_config = load_config()["ols"]
    docs = _query_ols(
        name,
        search_url=ols_config["search_url"],
        max_retries=ols_config["max_retries"],
        retry_backoff_seconds=ols_config["retry_backoff_seconds"],
        rows=ols_config["search_result_rows"],
    )

    matches = []
    for doc in docs:
        efo_id = doc.get("short_form")
        label = doc.get("label")
        synonyms = doc.get("synonym", []) or []

        if not efo_id or not label:
            continue  # skip malformed entries rather than crash

        confidence, match_type = _score_candidate(name, label, synonyms)
        matches.append(DiseaseMatch(
            query=name,
            efo_id=efo_id,
            label=label,
            confidence=confidence,
            match_type=match_type,
        ))

    matches.sort(key=lambda m: m.confidence, reverse=True)
    return matches


def resolve_with_cache(name: str, use_cache: bool = True) -> Optional[DiseaseMatch]:
    """
    Resolve a disease name, checking the local cache first. Only returns a
    match if it is auto-acceptable (confidence >= AUTO_ACCEPT_THRESHOLD).
    If the best match is below threshold, returns None so the caller can
    surface the candidates for human confirmation instead of guessing.
    """
    if use_cache:
        cached = _read_cache().get(name.strip().lower())
        if cached:
            return cached

    candidates = resolve_disease_name(name)

    if not candidates:
        return None

    best = candidates[0]
    if not best.is_auto_acceptable():
        return None  # caller must confirm -- see suggest_candidates_for_confirmation

    _append_to_cache(best)
    return best


def suggest_candidates_for_confirmation(name: str) -> List[DiseaseMatch]:
    """
    Returns the top-N ranked candidates for a name that did not auto-resolve,
    so a human can pick the correct one instead of the system guessing.
    """
    top_n = load_config()["ols"]["confirmation_candidates"]
    return resolve_disease_name(name)[:top_n]


def confirm_and_cache(match: DiseaseMatch) -> None:
    """Call this once a human has confirmed a low-confidence match is correct."""
    _append_to_cache(match)


def _read_cache() -> dict:
    resolved_cache_path = load_config()["ols"]["resolved_cache_path"]
    if not os.path.exists(resolved_cache_path):
        return {}

    cache = {}
    with open(resolved_cache_path, "r", newline="") as f:
        for row in csv.DictReader(f):
            cache[row["query"].strip().lower()] = DiseaseMatch(
                query=row["query"],
                efo_id=row["efo_id"],
                label=row["label"],
                confidence=float(row["confidence"]),
                match_type=row["match_type"],
            )
    return cache


def _append_to_cache(match: DiseaseMatch) -> None:
    resolved_cache_path = load_config()["ols"]["resolved_cache_path"]
    os.makedirs(os.path.dirname(resolved_cache_path), exist_ok=True)
    file_exists = os.path.exists(resolved_cache_path)

    with open(resolved_cache_path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["query", "efo_id", "label", "confidence", "match_type"])
        if not file_exists:
            writer.writeheader()
        writer.writerow({
            "query": match.query,
            "efo_id": match.efo_id,
            "label": match.label,
            "confidence": match.confidence,
            "match_type": match.match_type,
        })
