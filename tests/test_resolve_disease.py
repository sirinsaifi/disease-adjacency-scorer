"""
Tests for resolve_disease.py — resolving free-text disease names (including
vague or misspelled ones) to EFO IDs, without ever silently guessing wrong.

Uses pytest-mock to fake the OLS API response, so no live network call is
made during tests.

Run with: poetry run pytest tests/test_resolve_disease.py -v
"""

import pytest

from disease_adjacency.resolve_disease import (
    resolve_disease_name,
    _score_candidate,
    AUTO_ACCEPT_THRESHOLD,
)


def make_ols_response(docs):
    """Build a fake requests.Response-like object for mocking."""
    class FakeResponse:
        def raise_for_status(self):
            pass

        def json(self):
            return {"response": {"docs": docs}}

    return FakeResponse()


class TestScoreCandidate:

    def test_exact_label_match_scores_one(self):
        confidence, match_type = _score_candidate(
            "non-small cell lung carcinoma", "non-small cell lung carcinoma", []
        )
        assert confidence == 1.0
        assert match_type == "exact_label"

    def test_exact_label_match_is_case_insensitive(self):
        confidence, match_type = _score_candidate(
            "NON-SMALL CELL LUNG CARCINOMA", "non-small cell lung carcinoma", []
        )
        assert confidence == 1.0
        assert match_type == "exact_label"

    def test_exact_synonym_match_scores_one(self):
        confidence, match_type = _score_candidate(
            "NSCLC", "non-small cell lung carcinoma", ["NSCLC", "non-small cell lung cancer"]
        )
        assert confidence == 1.0
        assert match_type == "exact_synonym"

    def test_misspelled_name_gets_fuzzy_score_below_one(self):
        # "carcinoma" -> "carcinom" (missing a letter), should still score high but not 1.0
        confidence, match_type = _score_candidate(
            "non-small cell lung carcinom", "non-small cell lung carcinoma", []
        )
        assert match_type == "fuzzy"
        assert 0.9 < confidence < 1.0

    def test_completely_unrelated_name_scores_low(self):
        confidence, match_type = _score_candidate(
            "diabetes", "non-small cell lung carcinoma", []
        )
        assert match_type == "fuzzy"
        assert confidence < 0.5


class TestResolveDiseaseName:

    def test_exact_match_returns_top_candidate_with_full_confidence(self, mocker):
        mock_docs = [
            {"short_form": "EFO_0003060", "label": "non-small cell lung carcinoma", "synonym": ["NSCLC"]},
        ]
        mocker.patch(
            "disease_adjacency.resolve_disease.requests.get",
            return_value=make_ols_response(mock_docs),
        )

        results = resolve_disease_name("non-small cell lung carcinoma")
        assert len(results) == 1
        assert results[0].efo_id == "EFO_0003060"
        assert results[0].confidence == 1.0
        assert results[0].is_auto_acceptable()

    def test_misspelled_name_still_finds_best_candidate_but_may_need_confirmation(self, mocker):
        mock_docs = [
            {"short_form": "EFO_0003060", "label": "non-small cell lung carcinoma", "synonym": []},
        ]
        mocker.patch(
            "disease_adjacency.resolve_disease.requests.get",
            return_value=make_ols_response(mock_docs),
        )

        results = resolve_disease_name("non small cell lung carcinom")  # missing hyphen + typo
        assert len(results) == 1
        assert results[0].efo_id == "EFO_0003060"
        # Should score reasonably high but we don't assert auto-acceptance here --
        # that's a judgment call the AUTO_ACCEPT_THRESHOLD makes explicit.

    def test_ambiguous_name_returns_multiple_ranked_candidates(self, mocker):
        mock_docs = [
            {"short_form": "EFO_0000702", "label": "small cell lung carcinoma", "synonym": []},
            {"short_form": "EFO_0003060", "label": "non-small cell lung carcinoma", "synonym": []},
        ]
        mocker.patch(
            "disease_adjacency.resolve_disease.requests.get",
            return_value=make_ols_response(mock_docs),
        )

        results = resolve_disease_name("lung carcinoma")
        assert len(results) == 2
        # Results must be sorted by confidence descending
        assert results[0].confidence >= results[1].confidence

    def test_no_results_returns_empty_list_not_a_crash(self, mocker):
        mocker.patch(
            "disease_adjacency.resolve_disease.requests.get",
            return_value=make_ols_response([]),
        )

        results = resolve_disease_name("completely made up disease xyz123")
        assert results == []

    def test_malformed_doc_missing_short_form_is_skipped(self, mocker):
        mock_docs = [
            {"label": "no id here", "synonym": []},  # missing short_form
            {"short_form": "EFO_0003060", "label": "non-small cell lung carcinoma", "synonym": []},
        ]
        mocker.patch(
            "disease_adjacency.resolve_disease.requests.get",
            return_value=make_ols_response(mock_docs),
        )

        results = resolve_disease_name("non-small cell lung carcinoma")
        assert len(results) == 1
        assert results[0].efo_id == "EFO_0003060"


class TestAutoAcceptThreshold:

    def test_threshold_is_below_one_to_allow_near_exact_matches(self):
        # Sanity check the configured threshold is sensible (not 1.0, which
        # would make the "misspelled name" feature pointless)
        assert 0.0 < AUTO_ACCEPT_THRESHOLD < 1.0
