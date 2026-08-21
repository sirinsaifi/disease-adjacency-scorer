"""
fetch_opentargets.py

Pulls target, pathway, and genetics association data for a given disease
from the Open Targets Platform GraphQL API.

API endpoint: https://api.platform.opentargets.org/api/v4/graphql
Docs: https://platform-docs.opentargets.org
Interactive schema browser: https://api.platform.opentargets.org/api/v4/graphql/browser

Every response is cached to data/raw/{efo_id}.json so repeated pipeline runs
don't re-hit the live API.
"""

import json
import os
import time

import requests

from disease_adjacency.config import load_config

# Open Targets datatype ids used to split evidence into genetics vs pathway.
# Full list of datatype ids: genetic_association, somatic_mutation, known_drug,
# affected_pathway, literature, rna_expression, animal_model.
# (The genetics/pathway datatype ids actually used are read from
# config/settings.yaml's fingerprint section by build_fingerprint.py.)

# The real GraphQL query. For each disease we pull its top associated targets
# (sorted by overall association score), the per-datatype score breakdown
# (used to isolate the genetics and pathway signal), and each target's
# Reactome pathway membership (used to build the pathway-overlap set).
ASSOCIATED_TARGETS_QUERY = """
query DiseaseAssociatedTargets($efoId: String!, $size: Int!) {
  disease(efoId: $efoId) {
    id
    name
    associatedTargets(page: { index: 0, size: $size }) {
      count
      rows {
        score
        datatypeScores {
          id
          score
        }
        target {
          id
          approvedSymbol
          pathways {
            pathway
          }
        }
      }
    }
  }
}
"""


def _run_query(efo_id: str, size: int, api_url: str) -> dict:
    """Execute the GraphQL query against the live Open Targets API."""
    response = requests.post(
        api_url,
        json={
            "query": ASSOCIATED_TARGETS_QUERY,
            "variables": {"efoId": efo_id, "size": size},
        },
        timeout=30,
    )
    response.raise_for_status()
    payload = response.json()

    if "errors" in payload:
        raise RuntimeError(f"Open Targets API returned errors for {efo_id}: {payload['errors']}")

    if payload.get("data", {}).get("disease") is None:
        raise ValueError(f"No disease found for EFO ID {efo_id} -- check the ID is correct.")

    return payload["data"]["disease"]


def fetch_disease_associations(efo_id: str, use_cache: bool = True) -> dict:
    """
    Fetch target/pathway/genetics associations for a disease by its EFO ID.
    Checks local cache first if use_cache is True. Retries on transient
    network/API failures with simple linear backoff.

    All tunable values (API URL, fetch size, retry count, cache dir) come
    from config/settings.yaml, not hardcoded constants.
    """
    config = load_config()["open_targets"]
    raw_data_dir = config["raw_data_dir"]
    cache_path = os.path.join(raw_data_dir, f"{efo_id}.json")

    if use_cache and os.path.exists(cache_path):
        with open(cache_path, "r") as f:
            return json.load(f)

    last_error = None
    for attempt in range(1, config["max_retries"] + 1):
        try:
            data = _run_query(efo_id, size=config["fetch_size"], api_url=config["api_url"])
            save_to_cache(efo_id, data)
            return data
        except (requests.RequestException, RuntimeError) as e:
            last_error = e
            if attempt < config["max_retries"]:
                time.sleep(attempt * config["retry_backoff_seconds"])  # linear backoff

    raise RuntimeError(
        f"Failed to fetch data for {efo_id} after {config['max_retries']} attempts: {last_error}"
    )


def save_to_cache(efo_id: str, data: dict) -> None:
    raw_data_dir = load_config()["open_targets"]["raw_data_dir"]
    os.makedirs(raw_data_dir, exist_ok=True)
    cache_path = os.path.join(raw_data_dir, f"{efo_id}.json")
    with open(cache_path, "w") as f:
        json.dump(data, f, indent=2)
