"""
Tests for fetch_adverse_events.py -- building a disease's AE fingerprint
from Open Targets known drugs + openFDA FAERS. Uses mocked HTTP responses,
so no live network call is made during tests.

Run with: poetry run pytest tests/test_fetch_adverse_events.py -v
"""

import json

from disease_adjacency.fetch_adverse_events import fetch_ae_fingerprint


class FakeResponse:
    """Lightweight requests.Response stand-in, same pattern as
    test_resolve_disease.py's make_ols_response()."""

    def __init__(self, json_data, status_code=200):
        self._json_data = json_data
        self.status_code = status_code

    def raise_for_status(self):
        pass  # tests only ever exercise the 200 and 404 paths, never a raise here

    def json(self):
        return self._json_data


def known_drugs_response(drug_names):
    return FakeResponse({
        "data": {
            "disease": {
                "drugAndClinicalCandidates": {
                    "rows": [{"drug": {"name": name}} for name in drug_names]
                }
            }
        }
    })


def faers_response(terms, status_code=200):
    return FakeResponse({"results": [{"term": t} for t in terms]}, status_code=status_code)


def make_config(tmp_path):
    """
    Fake settings.yaml shape with raw_data_dir pointed at a pytest tmp_path,
    so tests never read or write the real project's data/raw/ cache.
    max_retries=1 keeps retry-exhaustion paths from sleeping for real.
    """
    return {
        "open_targets": {
            "raw_data_dir": str(tmp_path),
            "api_url": "https://fake-opentargets.test/graphql",
            "max_retries": 1,
            "retry_backoff_seconds": 0,
        },
        "openfda": {
            "api_url": "https://fake-openfda.test/event.json",
            "top_drugs_per_disease": 5,
            "top_ae_count": 10,
            "retry_backoff_seconds": 0,
        },
    }


class TestFetchAeFingerprint:

    def test_fetch_ae_fingerprint_returns_set(self, mocker, tmp_path):
        mocker.patch("disease_adjacency.fetch_adverse_events.load_config", return_value=make_config(tmp_path))
        mocker.patch(
            "disease_adjacency.fetch_adverse_events.requests.post",
            return_value=known_drugs_response(["DRUG_A", "DRUG_B"]),
        )
        mocker.patch(
            "disease_adjacency.fetch_adverse_events.requests.get",
            return_value=faers_response(["NAUSEA", "FATIGUE", "HEADACHE"]),
        )

        result = fetch_ae_fingerprint("EFO_TEST_0001", use_cache=False)

        assert isinstance(result, set)
        assert {"NAUSEA", "FATIGUE", "HEADACHE"}.issubset(result)

    def test_fetch_ae_fingerprint_deduplicates_terms(self, mocker, tmp_path):
        mocker.patch("disease_adjacency.fetch_adverse_events.load_config", return_value=make_config(tmp_path))
        mocker.patch(
            "disease_adjacency.fetch_adverse_events.requests.post",
            return_value=known_drugs_response(["DRUG_A", "DRUG_B"]),
        )
        # DRUG_A and DRUG_B share "FATIGUE" -- it must appear only once in the result set.
        mocker.patch(
            "disease_adjacency.fetch_adverse_events.requests.get",
            side_effect=[
                faers_response(["PAIN", "FATIGUE"]),
                faers_response(["FATIGUE", "NAUSEA"]),
            ],
        )

        result = fetch_ae_fingerprint("EFO_TEST_0002", use_cache=False)

        assert result == {"PAIN", "FATIGUE", "NAUSEA"}

    def test_fetch_ae_fingerprint_handles_404_gracefully(self, mocker, tmp_path):
        mocker.patch("disease_adjacency.fetch_adverse_events.load_config", return_value=make_config(tmp_path))
        mocker.patch(
            "disease_adjacency.fetch_adverse_events.requests.post",
            return_value=known_drugs_response(["DRUG_WITH_NO_FAERS_RECORDS"]),
        )
        mocker.patch(
            "disease_adjacency.fetch_adverse_events.requests.get",
            return_value=FakeResponse({}, status_code=404),
        )

        result = fetch_ae_fingerprint("EFO_TEST_0003", use_cache=False)

        assert isinstance(result, set)
        assert result == set()

    def test_fetch_ae_fingerprint_uses_cache(self, mocker, tmp_path):
        efo_id = "EFO_TEST_0004"
        cache_path = tmp_path / f"ae_{efo_id}.json"
        cache_path.write_text(json.dumps(["CACHED_TERM_A", "CACHED_TERM_B"]))

        mocker.patch("disease_adjacency.fetch_adverse_events.load_config", return_value=make_config(tmp_path))
        mock_get = mocker.patch("disease_adjacency.fetch_adverse_events.requests.get")
        mock_post = mocker.patch("disease_adjacency.fetch_adverse_events.requests.post")

        result = fetch_ae_fingerprint(efo_id, use_cache=True)

        assert result == {"CACHED_TERM_A", "CACHED_TERM_B"}
        mock_get.assert_not_called()
        mock_post.assert_not_called()
