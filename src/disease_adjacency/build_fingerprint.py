"""
build_fingerprint.py

Turns a raw Open Targets API response (from fetch_opentargets.py) plus an
optional adverse-event term set (from fetch_adverse_events.py) into a clean
"biology fingerprint": a dict of four sets ready for comparison in compare.py.

    {
        "targets":         set of approved gene symbols above min_target_score
        "pathways":        set of Reactome pathway names from those targets
        "genetics":        set of gene symbols with genetics evidence above threshold
        "adverse_events":  set of MedDRA AE terms (passed in, not derived from OT)
    }

The adverse_events set comes from fetch_adverse_events.fetch_ae_fingerprint()
and is passed in as an argument so this function stays testable without
hitting the openFDA API.

Thresholds come from config/settings.yaml (fingerprint section).
"""

from typing import Dict, Optional, Set

from disease_adjacency.config import load_config


def build_fingerprint(
    raw_response: dict,
    adverse_events: Optional[Set[str]] = None,
) -> Dict[str, Set[str]]:
    """
    Extract targets, pathways, genetics, and adverse-event sets from a raw
    Open Targets disease association response.

    raw_response : the `disease` object returned by fetch_opentargets.py
    adverse_events : AE term set from fetch_adverse_events.fetch_ae_fingerprint()
                     (pass None or omit to get an empty set -- pipeline always
                      passes the real set)

    Returns a fingerprint dict with four sets.
    """
    config = load_config()["fingerprint"]
    min_target_score  = config["min_target_score"]
    min_genetics_score = config["min_genetics_score"]
    genetics_datatype_id = config["genetics_datatype_id"]

    rows = (raw_response.get("associatedTargets") or {}).get("rows") or []

    targets:  Set[str] = set()
    pathways: Set[str] = set()
    genetics: Set[str] = set()

    for row in rows:
        overall_score = row.get("score", 0.0)
        target_info   = row.get("target") or {}
        symbol        = target_info.get("approvedSymbol")

        if not symbol:
            continue  # skip malformed rows

        if overall_score >= min_target_score:
            targets.add(symbol)

            for pathway in target_info.get("pathways", []) or []:
                name = pathway.get("pathway")
                if name:
                    pathways.add(name)

        datatype_scores = row.get("datatypeScores") or []
        genetics_score = next(
            (d["score"] for d in datatype_scores if d.get("id") == genetics_datatype_id),
            0.0,
        )
        if genetics_score >= min_genetics_score:
            genetics.add(symbol)

    return {
        "targets":        targets,
        "pathways":       pathways,
        "genetics":       genetics,
        "adverse_events": adverse_events if adverse_events is not None else set(),
    }
