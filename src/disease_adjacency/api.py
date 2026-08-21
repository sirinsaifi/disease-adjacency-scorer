"""
api.py

A thin HTTP layer over the scoring pipeline, so the interactive frontend
(output/index.html) can run live adjacency queries instead of only reading
the pre-computed output/adjacency_scores.json.

Every endpoint reuses the exact same functions pipeline.py's CLI run uses
(resolve_or_raise, build_disease_fingerprint, score_pair, ...) -- this file
adds no scoring logic of its own, only request/response plumbing. Disease
fetches still go through each module's existing on-disk cache (data/raw/,
data/efo/), so scoring a pair that's already been run costs no live API
calls; a genuinely new disease name still hits Open Targets/OLS4/openFDA
live, same as running the CLI pipeline would.

Run with (from the project root):
    poetry run uvicorn disease_adjacency.api:app --reload --port 8000

Note: NOT `uvicorn src.disease_adjacency.api:app` -- every other module in
this package imports as `from disease_adjacency.xxx import yyy` (absolute,
no "src." prefix), matching how `poetry install` registers the package
(see pyproject.toml's `packages = [{ include = "disease_adjacency", from
= "src" }]`). Loading this file as `src.disease_adjacency.api` would still
import fine on its own, but its internal `from disease_adjacency...`
imports would then fail to resolve (or, with PYTHONPATH=src set, load a
second, separate copy of every module under a different qualified name)
since "src" was never on sys.path for that project.
"""

import csv
import json
import os

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

from disease_adjacency.config import load_config
from disease_adjacency.pipeline import (
    DiseaseResolutionError,
    build_disease_fingerprint,
    score_pair,
)

app = FastAPI(title="Disease Adjacency Scorer API")

# Local dev tool only -- the frontend is served from a different localhost
# port (python3 -m http.server), so the browser treats it as cross-origin.
# Wide open is fine here since nothing in this API is exposed beyond
# localhost and it reads/writes only this project's own local files.
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)


class ScoreRequest(BaseModel):
    anchor: str
    candidates: list[str]


def _resolution_error_body(e: DiseaseResolutionError) -> dict:
    return {
        "error": f"Could not confidently resolve disease name '{e.name}'.",
        "disease": e.name,
        "candidates": [
            {
                "label": c.label,
                "efo_id": c.efo_id,
                "confidence": c.confidence,
                "match_type": c.match_type,
            }
            for c in e.candidates
        ],
    }


@app.get("/diseases")
def get_diseases():
    """
    Returns the anchor and candidate disease names from config/diseases.csv,
    so the frontend can pre-fill the anchor field and populate the
    candidate multi-select.
    """
    disease_list_path = load_config()["paths"]["disease_list"]

    anchor = None
    candidates = []
    with open(disease_list_path, "r") as f:
        for row in csv.DictReader(f):
            role = row["role"].strip().lower()
            name = row["disease_name"].strip()
            if role == "anchor":
                anchor = name
            elif role == "candidate":
                candidates.append(name)

    return {"anchor": anchor, "candidates": candidates}


@app.post("/score")
def post_score(req: ScoreRequest):
    """
    Scores each of req.candidates against req.anchor using the same
    scoring path as the CLI pipeline. Candidates that fail to resolve are
    reported in "errors" rather than failing the whole request -- one bad
    disease name shouldn't block scoring the others.
    """
    try:
        anchor_efo_id, anchor_fp = build_disease_fingerprint(req.anchor)
    except DiseaseResolutionError as e:
        raise HTTPException(status_code=422, detail=_resolution_error_body(e))

    results = []
    errors = []
    for name in req.candidates:
        try:
            results.append(score_pair(req.anchor, anchor_efo_id, anchor_fp, name))
        except DiseaseResolutionError as e:
            errors.append(_resolution_error_body(e))

    results.sort(key=lambda r: r["adjacency_score"], reverse=True)
    return {"anchor": req.anchor, "anchor_efo_id": anchor_efo_id, "results": results, "errors": errors}


@app.get("/results")
def get_results():
    """Returns the full pre-computed output/adjacency_scores.json as-is."""
    output_path = load_config()["paths"]["output_scores"]
    if not os.path.exists(output_path):
        raise HTTPException(
            status_code=404,
            detail=f"No precomputed results at {output_path} -- run the pipeline first.",
        )
    with open(output_path, "r") as f:
        return json.load(f)


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8000)
