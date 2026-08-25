"""
Tests for fetch_ot_indirect.py -- the EFO ontology-proximity score.
Uses mocked HTTP responses / mocked internals, so no live network call is
made during tests.

Note: there is no _fetch_ancestors() function in fetch_ot_indirect.py --
the real internal helper is _ols_related_ids(efo_id, relation, ...), used
for both "parents" and "ancestors" lookups via its `relation` argument.
Tests that need to control what OLS4 "returns" mock that function directly.

Run with: poetry run pytest tests/test_fetch_ot_indirect.py -v
"""

import json

import pytest
import requests

from disease_adjacency.fetch_ot_indirect import fetch_ot_indirect_score


class FakeResponse:
    def __init__(self, json_data, status_code=200):
        self._json_data = json_data
        self.status_code = status_code

    def raise_for_status(self):
        pass

    def json(self):
        return self._json_data


def empty_terms_response():
    return FakeResponse({"_embedded": {"terms": []}}, status_code=200)


def make_config(tmp_path):
    """
    Fake settings.yaml shape with raw_data_dir pointed at a pytest tmp_path,
    so tests never read or write the real project's data/raw/ cache.
    """
    return {
        "open_targets": {"raw_data_dir": str(tmp_path)},
        "ols": {
            "term_base_url": "https://fake-ols4.test/terms",
            "max_retries": 1,
            "retry_backoff_seconds": 0,
            "ancestor_page_size": 200,
            "direct_relationship_score": 0.67,
            "shared_ancestor_score_ceiling": 0.49,
        },
    }


class TestFetchOtIndirectScore:

    def test_same_disease_returns_1(self, mocker, tmp_path):
        mocker.patch("disease_adjacency.fetch_ot_indirect.load_config", return_value=make_config(tmp_path))
        mock_get = mocker.patch("disease_adjacency.fetch_ot_indirect.requests.get")

        result = fetch_ot_indirect_score("EFO_0000001", "EFO_0000001", use_cache=False)

        assert result == 1.0
        mock_get.assert_not_called()

    def test_direct_relationship_score(self, mocker, tmp_path):
        config = make_config(tmp_path)
        mocker.patch("disease_adjacency.fetch_ot_indirect.load_config", return_value=config)

        # The "direct" tier is decided from the *parents* relation, not
        # ancestors -- the candidate is the anchor's direct parent here.
        def fake_related_ids(efo_id, relation, *args, **kwargs):
            if relation == "parents" and efo_id == "EFO_ANCHOR":
                return {"EFO_CANDIDATE"}
            return set()

        mocker.patch("disease_adjacency.fetch_ot_indirect._ols_related_ids", side_effect=fake_related_ids)

        result = fetch_ot_indirect_score("EFO_ANCHOR", "EFO_CANDIDATE", use_cache=False)

        assert result == config["ols"]["direct_relationship_score"]

    def test_shared_ancestors_produces_score_between_0_and_1(self, mocker, tmp_path):
        mocker.patch("disease_adjacency.fetch_ot_indirect.load_config", return_value=make_config(tmp_path))

        def fake_related_ids(efo_id, relation, *args, **kwargs):
            if relation == "parents":
                return set()  # no direct relationship either way
            # ancestors -- overlapping sets
            return {"EFO_SHARED", "EFO_ANCHOR_ONLY"} if efo_id == "EFO_ANCHOR" else {"EFO_SHARED"}

        mocker.patch("disease_adjacency.fetch_ot_indirect._ols_related_ids", side_effect=fake_related_ids)

        result = fetch_ot_indirect_score("EFO_ANCHOR", "EFO_CANDIDATE", use_cache=False)

        assert 0.0 < result < 1.0

    def test_no_shared_ancestors_returns_0(self, mocker, tmp_path):
        mocker.patch("disease_adjacency.fetch_ot_indirect.load_config", return_value=make_config(tmp_path))

        def fake_related_ids(efo_id, relation, *args, **kwargs):
            if relation == "parents":
                return set()
            return {"EFO_ANCHOR_ANCESTOR"} if efo_id == "EFO_ANCHOR" else {"EFO_CANDIDATE_ANCESTOR"}

        mocker.patch("disease_adjacency.fetch_ot_indirect._ols_related_ids", side_effect=fake_related_ids)

        result = fetch_ot_indirect_score("EFO_ANCHOR", "EFO_CANDIDATE", use_cache=False)

        assert result == 0.0

    def test_api_failure_returns_0_gracefully(self, mocker, tmp_path):
        mocker.patch("disease_adjacency.fetch_ot_indirect.load_config", return_value=make_config(tmp_path))
        mocker.patch(
            "disease_adjacency.fetch_ot_indirect.requests.get",
            side_effect=requests.RequestException("network down"),
        )

        result = fetch_ot_indirect_score("EFO_ANCHOR", "EFO_CANDIDATE", use_cache=False)

        assert result == 0.0

    def test_result_is_cached(self, mocker, tmp_path):
        config = make_config(tmp_path)
        mocker.patch("disease_adjacency.fetch_ot_indirect.load_config", return_value=config)
        mock_get = mocker.patch(
            "disease_adjacency.fetch_ot_indirect.requests.get",
            return_value=empty_terms_response(),
        )

        anchor, candidate = "EFO_ANCHOR", "EFO_CANDIDATE"
        first = fetch_ot_indirect_score(anchor, candidate, use_cache=True)

        cache_path = tmp_path / f"ot_indirect_{anchor}_{candidate}.json"
        assert cache_path.exists()
        with open(cache_path) as f:
            assert json.load(f)["score"] == first

        assert mock_get.call_count > 0  # the first run genuinely hit the (mocked) network
        calls_after_first_run = mock_get.call_count

        second = fetch_ot_indirect_score(anchor, candidate, use_cache=True)

        assert second == first
        assert mock_get.call_count == calls_after_first_run  # no new calls on the cached run
