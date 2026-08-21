"""
Tests for build_fingerprint.py — parsing a raw Open Targets response into a
clean biology fingerprint. Uses realistic mocked responses shaped exactly
like the real API (see fetch_opentargets.py's ASSOCIATED_TARGETS_QUERY),
so no live network call is needed.

Run with: poetry run pytest tests/test_build_fingerprint.py -v
"""

from disease_adjacency.build_fingerprint import build_fingerprint


def make_row(symbol, overall_score, genetics_score=0.0, pathways=None):
    return {
        "score": overall_score,
        "datatypeScores": [
            {"id": "genetic_association", "score": genetics_score},
        ],
        "target": {
            "id": f"ENSG_{symbol}",
            "approvedSymbol": symbol,
            "pathways": [{"pathway": p} for p in (pathways or [])],
        },
    }


class TestBuildFingerprint:

    def test_extracts_targets_above_threshold(self):
        raw = {
            "associatedTargets": {
                "rows": [
                    make_row("EGFR", 0.8, genetics_score=0.5, pathways=["Signaling by EGFR"]),
                    make_row("NOISE", 0.01),  # below MIN_TARGET_SCORE, should be excluded
                ]
            }
        }
        fp = build_fingerprint(raw)
        assert fp["targets"] == {"EGFR"}
        assert "NOISE" not in fp["targets"]

    def test_extracts_pathways_only_for_included_targets(self):
        raw = {
            "associatedTargets": {
                "rows": [
                    make_row("EGFR", 0.8, pathways=["Signaling by EGFR", "MAPK signaling"]),
                ]
            }
        }
        fp = build_fingerprint(raw)
        assert fp["pathways"] == {"Signaling by EGFR", "MAPK signaling"}

    def test_genetics_uses_genetic_association_datatype_score(self):
        raw = {
            "associatedTargets": {
                "rows": [
                    make_row("EGFR", 0.8, genetics_score=0.5),   # above MIN_GENETICS_SCORE
                    make_row("KRAS", 0.5, genetics_score=0.01),  # below MIN_GENETICS_SCORE
                ]
            }
        }
        fp = build_fingerprint(raw)
        assert fp["genetics"] == {"EGFR"}

    def test_handles_empty_response_without_crashing(self):
        raw = {"associatedTargets": {"rows": []}}
        fp = build_fingerprint(raw)
        assert fp == {"targets": set(), "pathways": set(), "genetics": set(), "adverse_events": set()}

    def test_handles_missing_associated_targets_key(self):
        # Defensive: a malformed response shouldn't crash the whole pipeline
        raw = {}
        fp = build_fingerprint(raw)
        assert fp == {"targets": set(), "pathways": set(), "genetics": set(), "adverse_events": set()}

    def test_skips_rows_with_no_approved_symbol(self):
        raw = {
            "associatedTargets": {
                "rows": [
                    {"score": 0.9, "datatypeScores": [], "target": {"id": "ENSG1", "approvedSymbol": None}},
                    make_row("EGFR", 0.8),
                ]
            }
        }
        fp = build_fingerprint(raw)
        assert fp["targets"] == {"EGFR"}
