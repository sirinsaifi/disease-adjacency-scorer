# Disease Adjacency Scorer

Scores mechanistic adjacency between disease indication pairs, anchored on 2L non-small cell lung carcinoma (NSCLC).

## What it does

Scores mechanistic adjacency between disease indication pairs using five biological evidence dimensions sourced from Open Targets, EFO/OLS4, and openFDA FAERS. Every score decomposes into an inspectable evidence record — the specific shared genes, pathways, genetic loci, and adverse-event terms that produced it — not a black-box number. Validated against a hand-built gold set of 12 candidate diseases anchored on 2L NSCLC.

## Evidence dimensions

| Dimension | Source | What it captures | Weight |
|---|---|---|---|
| target_score | Open Targets | Shared associated genes | 0.55 |
| genetics_score | Open Targets GWAS | Shared genetic loci | 0.25 |
| ae_score | openFDA FAERS | Shared adverse-event profile | 0.15 |
| pathway_score | Open Targets / Reactome | Shared biological pathways | 0.04 |
| ot_indirect_score | EFO via OLS4 | Ontology proximity | 0.01 |

## Quickstart

```
poetry install
poetry run python src/disease_adjacency/pipeline.py
poetry run python validation/validate.py
```

## Interactive frontend

Terminal 1 (API):
```
poetry run uvicorn disease_adjacency.api:app --reload --port 8000
```

Terminal 2 (frontend):
```
python3 -m http.server 8080 --directory output
```

Open: http://localhost:8080/index.html

## Validation result

Spearman ρ = 0.6154 across 12 candidate diseases vs. the hand-built gold set.

Target was ρ ≥ 0.70 — see `docs/adjacency_definition.md` for the documented limitations explaining the gap.

## Project structure

```
src/disease_adjacency/
  pipeline.py             CLI entry point -- orchestrates the full scoring run
  api.py                  FastAPI layer for the interactive frontend
  resolve_disease.py      Free-text disease name -> EFO ID, via OLS4
  fetch_opentargets.py    Target/pathway/genetics associations from Open Targets
  fetch_adverse_events.py Known-drug AE fingerprint from Open Targets + openFDA FAERS
  fetch_ot_indirect.py    EFO ontology-proximity score via OLS4
  build_fingerprint.py    Raw API responses -> {targets, pathways, genetics, adverse_events}
  compare.py              Jaccard similarity across the four fingerprint dimensions
  score.py                Combines the five overlap scores into one weighted adjacency score
  evidence.py             Builds the inspectable per-pair evidence record
  config.py               Loads config/settings.yaml

config/
  settings.yaml           All tunable values -- API URLs, thresholds, weights
  diseases.csv            Anchor + candidate disease list

validation/
  gold_set.csv            Hand-ranked 12 diseases for Spearman evaluation
  validate.py             Compares pipeline output to the gold set

docs/
  adjacency_definition.md Scoring model definition, validation results, known limitations

output/
  index.html              Interactive frontend (live queries + baseline results)
  adjacency_scores.json   Pre-computed 12-disease baseline run
  run_summary.md          Static run summary

tests/
  test_*.py               Unit tests, mocked HTTP -- no live API calls
```

## Configuration

All tunable values live in `config/settings.yaml`.

To add diseases: edit `config/diseases.csv`.
To adjust weights: edit `scoring.weights` in `settings.yaml` (must sum to 1.0).

## Running tests

```
poetry run pytest
```

Expected: all tests pass (no count — it will grow over time).
