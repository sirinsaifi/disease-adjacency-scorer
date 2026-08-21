"""
validate.py

Compares the pipeline's ranked adjacency output against the hand-built
gold set, computing a Spearman rank correlation and printing a full
match/mismatch table.

Run with: poetry run python validation/validate.py
"""

import csv
import json
import sys

from disease_adjacency.config import load_config

# This script prints Unicode markers (✓ ✗ △). Windows consoles default to a
# codepage (e.g. cp1252) that can't encode them, which crashes print() --
# force UTF-8 stdout so `poetry run python validation/validate.py` works
# out of the box on any platform, matching the documented run instructions.
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

try:
    from scipy.stats import spearmanr
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False


def load_gold_set(path: str = None) -> dict:
    """Returns {disease_name: expected_rank}."""
    path = path or load_config()["paths"]["gold_set"]
    gold = {}
    with open(path, "r") as f:
        for row in csv.DictReader(f):
            gold[row["disease_name"]] = int(row["expected_rank"])
    return gold


def load_pipeline_results(path: str = None) -> list:
    path = path or load_config()["paths"]["output_scores"]
    with open(path, "r") as f:
        return json.load(f)


def compare_rankings(gold: dict, results: list) -> None:
    validation_config = load_config()["validation"]
    target_rho = validation_config["target_rho"]
    moderate_rho = validation_config["moderate_rho"]

    print(f"\n{'Disease':<45} {'Gold':>6} {'Actual':>7} {'Score':>8}  Match?")
    print("-" * 80)

    gold_ranks   = []
    actual_ranks = []

    for actual_rank, record in enumerate(results, start=1):
        name = record["disease_b"]
        expected_rank = gold.get(name)

        if expected_rank is None:
            match_label = "not in gold"
        else:
            match_label = "✓" if expected_rank == actual_rank else "✗"
            gold_ranks.append(expected_rank)
            actual_ranks.append(actual_rank)

        print(
            f"{name:<45} {str(expected_rank or '?'):>6} {actual_rank:>7} "
            f"{record['adjacency_score']:>8.4f}  {match_label}"
        )

    print("-" * 80)

    if SCIPY_AVAILABLE and len(gold_ranks) >= 3:
        rho, p_value = spearmanr(gold_ranks, actual_ranks)
        print(f"\nSpearman ρ = {rho:.4f}  (p = {p_value:.4f})")
        if rho >= target_rho:
            print(f"✓ Ranking correlation is strong (ρ ≥ {target_rho:.2f})")
        elif rho >= moderate_rho:
            print(f"△ Ranking correlation is moderate (ρ ≥ {moderate_rho:.2f}) — consider weight tuning")
        else:
            print(f"✗ Ranking correlation is weak (ρ < {moderate_rho:.2f}) — review weights and fingerprint thresholds")
    elif not SCIPY_AVAILABLE:
        print("\nInstall scipy for Spearman correlation: pip install scipy")
    else:
        print("\nNeed at least 3 matched pairs for Spearman correlation.")


def component_summary(results: list) -> None:
    """Print a breakdown of which dimension is driving each score."""
    print(f"\n{'Disease':<45} {'target':>7} {'pathway':>8} {'genetics':>9} {'ot_ind':>7} {'ae':>6}")
    print("-" * 90)
    for record in results:
        c = record.get("components", {})
        print(
            f"{record['disease_b']:<45} "
            f"{c.get('target_score', 0):.3f}   "
            f"{c.get('pathway_score', 0):.3f}    "
            f"{c.get('genetics_score', 0):.3f}      "
            f"{c.get('ot_indirect_score', 0):.3f}   "
            f"{c.get('ae_score', 0):.3f}"
        )


if __name__ == "__main__":
    gold    = load_gold_set()
    results = load_pipeline_results()

    compare_rankings(gold, results)
    component_summary(results)
