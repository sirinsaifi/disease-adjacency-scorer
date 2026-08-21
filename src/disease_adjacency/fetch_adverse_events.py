"""
fetch_adverse_events.py

Builds a disease's adverse-event (AE) fingerprint: the union of MedDRA
adverse-event terms reported for the drugs used to treat it.

Steps (see docs/adjacency_definition.md, dimension 5):
  1. Pull the disease's top known drugs (approved + clinical candidates)
     from Open Targets (`drugAndClinicalCandidates`).
  2. For each drug, query openFDA FAERS for its top MedDRA adverse-event
     terms (`patient.reaction.reactionmeddrapt.exact`, ranked by count).
  3. Union all terms into one AE fingerprint for the disease.

A drug with no FAERS records (openFDA returns 404 for "no matches") just
contributes an empty set rather than failing the whole fingerprint -- FAERS
coverage is inherently patchy and one missing drug isn't a pipeline error.

Every disease's AE fingerprint is cached to data/raw/ae_{efo_id}.json so
repeated pipeline runs don't re-hit Open Targets/openFDA.
"""

import json
import os
import time
from typing import List, Set

import requests

from disease_adjacency.config import load_config

KNOWN_DRUGS_QUERY = """
query DiseaseKnownDrugs($efoId: String!) {
  disease(efoId: $efoId) {
    drugAndClinicalCandidates {
      rows {
        drug {
          name
        }
      }
    }
  }
}
"""


def _fetch_known_drug_names(
    efo_id: str, top_n: int, api_url: str, max_retries: int, retry_backoff_seconds: int
) -> List[str]:
    """
    Fetch the disease's known drugs from Open Targets and return up to
    `top_n` distinct drug names, in the order Open Targets returns them.
    The API has no server-side limit/pagination for this field, so the
    slicing happens client-side.
    """
    last_error = None
    for attempt in range(1, max_retries + 1):
        try:
            response = requests.post(
                api_url,
                json={"query": KNOWN_DRUGS_QUERY, "variables": {"efoId": efo_id}},
                timeout=30,
            )
            response.raise_for_status()
            payload = response.json()

            if "errors" in payload:
                raise RuntimeError(f"Open Targets API returned errors for {efo_id}: {payload['errors']}")

            disease = payload.get("data", {}).get("disease")
            if disease is None:
                raise ValueError(f"No disease found for EFO ID {efo_id} -- check the ID is correct.")

            rows = (disease.get("drugAndClinicalCandidates") or {}).get("rows") or []
            names = []
            seen = set()
            for row in rows:
                name = (row.get("drug") or {}).get("name")
                if name and name not in seen:
                    seen.add(name)
                    names.append(name)
                if len(names) >= top_n:
                    break
            return names
        except (requests.RequestException, RuntimeError) as e:
            last_error = e
            if attempt < max_retries:
                time.sleep(attempt * retry_backoff_seconds)

    raise RuntimeError(
        f"Failed to fetch known drugs for {efo_id} after {max_retries} attempts: {last_error}"
    )


def _fetch_faers_terms(
    drug_name: str, top_n: int, api_url: str, max_retries: int, retry_backoff_seconds: int
) -> Set[str]:
    """
    Query openFDA FAERS for a drug's top MedDRA adverse-event terms.
    Returns an empty set (not an error) if openFDA has no records for the
    drug, or if the lookup ultimately fails after retries -- FAERS coverage
    is patchy and one drug's failure shouldn't sink the whole fingerprint.
    """
    params = {
        "search": f'patient.drug.medicinalproduct:"{drug_name}"',
        "count": "patient.reaction.reactionmeddrapt.exact",
        "limit": top_n,
    }

    last_error = None
    for attempt in range(1, max_retries + 1):
        try:
            response = requests.get(api_url, params=params, timeout=20)
            if response.status_code == 404:
                return set()  # no FAERS reports for this drug
            response.raise_for_status()
            results = response.json().get("results", []) or []
            return {row["term"] for row in results if row.get("term")}
        except requests.RequestException as e:
            last_error = e
            if attempt < max_retries:
                time.sleep(attempt * retry_backoff_seconds)

    print(f"Warning: FAERS lookup failed for '{drug_name}' after {max_retries} attempts: {last_error}")
    return set()


def fetch_ae_fingerprint(efo_id: str, use_cache: bool = True) -> Set[str]:
    """
    Build the AE fingerprint (union of MedDRA terms) for a disease.
    Checks the local cache first if use_cache is True.
    """
    config = load_config()
    ot_config = config["open_targets"]
    fda_config = config["openfda"]

    raw_data_dir = ot_config["raw_data_dir"]
    cache_path = os.path.join(raw_data_dir, f"ae_{efo_id}.json")

    if use_cache and os.path.exists(cache_path):
        with open(cache_path, "r") as f:
            return set(json.load(f))

    # Drug count per disease is controlled by openfda.top_drugs_per_disease
    # in settings.yaml, not by anything in the open_targets section -- the
    # drugs themselves come from Open Targets, but the *count* is an
    # openFDA/FAERS-side concern, since that's what top_ae_count pairs with.
    drug_names = _fetch_known_drug_names(
        efo_id,
        top_n=fda_config["top_drugs_per_disease"],
        api_url=ot_config["api_url"],
        max_retries=ot_config["max_retries"],
        retry_backoff_seconds=ot_config["retry_backoff_seconds"],
    )

    ae_terms: Set[str] = set()
    for drug_name in drug_names:
        ae_terms |= _fetch_faers_terms(
            drug_name,
            top_n=fda_config["top_ae_count"],
            api_url=fda_config["api_url"],
            max_retries=ot_config["max_retries"],
            retry_backoff_seconds=fda_config["retry_backoff_seconds"],
        )

    os.makedirs(raw_data_dir, exist_ok=True)
    with open(cache_path, "w") as f:
        json.dump(sorted(ae_terms), f, indent=2)

    return ae_terms
