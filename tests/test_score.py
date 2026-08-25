"""
Tests for score.py — the weighted adjacency scoring formula.
Run with: poetry run pytest tests/test_score.py -v
"""

import pytest
from disease_adjacency.score import (
    compute_adjacency_score,
    validate_weights,
    get_configured_weights,
)


class TestValidateWeights:

    def test_configured_weights_sum_to_one(self):
        # Protects against a typo in config/settings.yaml silently breaking the formula
        weights = get_configured_weights()
        validate_weights(weights)  # should not raise

    def test_weights_not_summing_to_one_raises_error(self):
        bad_weights = {
            "target_score": 0.5,
            "pathway_score": 0.5,
            "genetics_score": 0.5,
            "ot_indirect_score": 0.5,
            "ae_score": 0.5,
        }
        with pytest.raises(ValueError):
            validate_weights(bad_weights)

    def test_missing_key_raises_error(self):
        bad_weights = {
            "target_score": 0.5,
            "pathway_score": 0.5,
            "genetics_score": 0.5,
            "ot_indirect_score": 0.5,
            # missing ae_score
        }
        with pytest.raises(ValueError):
            validate_weights(bad_weights)

    def test_all_five_correct_keys_does_not_raise(self):
        weights = {
            "target_score": 0.30,
            "pathway_score": 0.20,
            "genetics_score": 0.20,
            "ot_indirect_score": 0.15,
            "ae_score": 0.15,
        }
        validate_weights(weights)  # should not raise

    def test_only_three_keys_raises_error(self):
        weights = {
            "target_score": 0.34,
            "pathway_score": 0.33,
            "genetics_score": 0.33,
        }
        with pytest.raises(ValueError):
            validate_weights(weights)

    def test_sum_within_tolerance_band_does_not_raise(self):
        # validate_weights accepts 0.999 <= sum <= 1.001 (inclusive) as a
        # floating-point tolerance band, not a strict boundary -- both ends
        # of that band are valid, not rejected.
        weights_high_end = {
            "target_score": 0.301, "pathway_score": 0.20,
            "genetics_score": 0.20, "ot_indirect_score": 0.15, "ae_score": 0.15,
        }  # sums to 1.001
        weights_low_end = {
            "target_score": 0.299, "pathway_score": 0.20,
            "genetics_score": 0.20, "ot_indirect_score": 0.15, "ae_score": 0.15,
        }  # sums to 0.999
        validate_weights(weights_high_end)  # should not raise
        validate_weights(weights_low_end)  # should not raise

    def test_sum_outside_tolerance_band_raises_error(self):
        weights_too_high = {
            "target_score": 0.302, "pathway_score": 0.20,
            "genetics_score": 0.20, "ot_indirect_score": 0.15, "ae_score": 0.15,
        }  # sums to 1.002
        weights_too_low = {
            "target_score": 0.298, "pathway_score": 0.20,
            "genetics_score": 0.20, "ot_indirect_score": 0.15, "ae_score": 0.15,
        }  # sums to 0.998
        with pytest.raises(ValueError):
            validate_weights(weights_too_high)
        with pytest.raises(ValueError):
            validate_weights(weights_too_low)


class TestComputeAdjacencyScore:

    def test_known_example_matches_expected_value(self):
        # Worked example -- weights passed explicitly so this test's result
        # doesn't depend on whatever happens to be configured in
        # config/settings.yaml (which changes as scoring.weights is tuned).
        overlap_scores = {
            "target_score": 0.29,
            "pathway_score": 0.55,
            "genetics_score": 0.15,
            "ot_indirect_score": 0.67,
            "ae_score": 0.40,
        }
        weights = {
            "target_score": 0.30,
            "pathway_score": 0.20,
            "genetics_score": 0.20,
            "ot_indirect_score": 0.15,
            "ae_score": 0.15,
        }
        result = compute_adjacency_score(overlap_scores, weights=weights)
        # (0.30*0.29) + (0.20*0.55) + (0.20*0.15) + (0.15*0.67) + (0.15*0.40) = 0.3875
        assert result == pytest.approx(0.3875, rel=1e-2)

    def test_perfect_overlap_gives_score_of_one(self):
        overlap_scores = {
            "target_score": 1.0,
            "pathway_score": 1.0,
            "genetics_score": 1.0,
            "ot_indirect_score": 1.0,
            "ae_score": 1.0,
        }
        assert compute_adjacency_score(overlap_scores) == 1.0

    def test_zero_overlap_gives_score_of_zero(self):
        overlap_scores = {
            "target_score": 0.0,
            "pathway_score": 0.0,
            "genetics_score": 0.0,
            "ot_indirect_score": 0.0,
            "ae_score": 0.0,
        }
        assert compute_adjacency_score(overlap_scores) == 0.0

    def test_score_always_within_valid_range(self):
        overlap_scores = {
            "target_score": 0.8,
            "pathway_score": 0.9,
            "genetics_score": 0.7,
            "ot_indirect_score": 0.67,
            "ae_score": 0.6,
        }
        result = compute_adjacency_score(overlap_scores)
        assert 0.0 <= result <= 1.0

    def test_custom_weights_are_respected(self):
        overlap_scores = {
            "target_score": 1.0,
            "pathway_score": 0.0,
            "genetics_score": 0.0,
            "ot_indirect_score": 0.0,
            "ae_score": 0.0,
        }
        custom_weights = {
            "target_score": 1.0,
            "pathway_score": 0.0,
            "genetics_score": 0.0,
            "ot_indirect_score": 0.0,
            "ae_score": 0.0,
        }
        result = compute_adjacency_score(overlap_scores, weights=custom_weights)
        assert result == 1.0

    def test_out_of_range_overlap_score_raises_error(self):
        overlap_scores = {
            "target_score": 1.5,  # invalid — overlap scores must be 0-1
            "pathway_score": 0.5,
            "genetics_score": 0.5,
            "ot_indirect_score": 0.5,
            "ae_score": 0.5,
        }
        with pytest.raises(ValueError):
            compute_adjacency_score(overlap_scores)

    def test_missing_overlap_key_raises_error(self):
        overlap_scores = {
            "target_score": 0.5,
            "pathway_score": 0.5,
            "genetics_score": 0.5,
            "ot_indirect_score": 0.5,
            # missing ae_score
        }
        with pytest.raises(ValueError):
            compute_adjacency_score(overlap_scores)
