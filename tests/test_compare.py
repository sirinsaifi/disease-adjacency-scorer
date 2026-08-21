"""
Tests for compare.py — overlap calculations between disease fingerprints.
Run with: poetry run pytest tests/test_compare.py -v
"""

import pytest
from disease_adjacency.compare import (
    jaccard_similarity,
    compare_fingerprints,
    get_shared_evidence,
)


class TestJaccardSimilarity:

    def test_partial_overlap(self):
        set_a = {"EGFR", "KRAS", "ALK", "TP53", "MET"}
        set_b = {"TP53", "RB1", "MYC", "KRAS"}
        # shared = {KRAS, TP53} = 2, union = 7 -> 2/7
        result = jaccard_similarity(set_a, set_b)
        assert result == pytest.approx(2 / 7, rel=1e-3)

    def test_identical_sets_give_score_of_one(self):
        set_a = {"EGFR", "KRAS"}
        set_b = {"EGFR", "KRAS"}
        assert jaccard_similarity(set_a, set_b) == 1.0

    def test_no_overlap_gives_zero(self):
        set_a = {"EGFR", "KRAS"}
        set_b = {"BRCA1", "BRCA2"}
        assert jaccard_similarity(set_a, set_b) == 0.0

    def test_both_empty_sets_give_zero_not_error(self):
        # Edge case: no data for either disease should not crash or divide by zero
        assert jaccard_similarity(set(), set()) == 0.0

    def test_one_empty_set_gives_zero(self):
        set_a = {"EGFR"}
        set_b = set()
        assert jaccard_similarity(set_a, set_b) == 0.0


class TestCompareFingerprints:

    def test_returns_all_four_scores(self):
        fp_a = {
            "targets": {"EGFR", "KRAS"},
            "pathways": {"MAPK signaling"},
            "genetics": {"rs123"},
            "adverse_events": {"pneumonitis", "fatigue"},
        }
        fp_b = {
            "targets": {"KRAS", "TP53"},
            "pathways": {"MAPK signaling", "PI3K-AKT signaling"},
            "genetics": {"rs456"},
            "adverse_events": {"fatigue"},
        }
        result = compare_fingerprints(fp_a, fp_b)
        assert set(result.keys()) == {"target_score", "pathway_score", "genetics_score", "ae_score"}
        assert result["target_score"] == pytest.approx(1 / 3, rel=1e-3)
        assert result["pathway_score"] == pytest.approx(1 / 2, rel=1e-3)
        assert result["genetics_score"] == 0.0
        assert result["ae_score"] == pytest.approx(1 / 2, rel=1e-3)

    def test_missing_key_raises_error(self):
        fp_a = {"targets": {"EGFR"}, "pathways": set(), "genetics": set()}  # missing "adverse_events"
        fp_b = {
            "targets": {"EGFR"},
            "pathways": set(),
            "genetics": set(),
            "adverse_events": set(),
        }
        with pytest.raises(ValueError):
            compare_fingerprints(fp_a, fp_b)


class TestGetSharedEvidence:

    def test_returns_actual_shared_items_not_just_counts(self):
        fp_a = {
            "targets": {"EGFR", "KRAS"},
            "pathways": {"MAPK signaling"},
            "genetics": {"rs123"},
            "adverse_events": {"pneumonitis", "fatigue"},
        }
        fp_b = {
            "targets": {"KRAS", "TP53"},
            "pathways": {"MAPK signaling"},
            "genetics": {"rs999"},
            "adverse_events": {"fatigue"},
        }
        evidence = get_shared_evidence(fp_a, fp_b)
        assert evidence["shared_targets"] == {"KRAS"}
        assert evidence["shared_pathways"] == {"MAPK signaling"}
        assert evidence["shared_genetics"] == set()
        assert evidence["shared_adverse_events"] == {"fatigue"}
