# Disease Adjacency Scorer — Run Summary

## What this is

A scored, evidence-linked adjacency system for disease pairs, anchored on
2L NSCLC (non-small cell lung carcinoma). It combines five biological
evidence dimensions — target overlap, pathway overlap, genetics overlap,
EFO ontology proximity, and adverse-event overlap — sourced from Open
Targets, EFO/OLS4, and openFDA FAERS into a single weighted adjacency
score per disease pair. Every score traces back to the specific shared
genes, pathways, genetic loci, and adverse-event terms that produced it —
there are no black-box numbers.

## How to run

1. Install: `poetry install`
2. Add diseases to `config/diseases.csv` (`role`: `anchor` or `candidate`)
3. Run pipeline: `poetry run python src/disease_adjacency/pipeline.py`
4. Run validation: `poetry run python validation/validate.py`
5. Output: `output/adjacency_scores.json`

## Final configuration

| Weight | Value |
|---|---|
| target_score | 0.55 |
| genetics_score | 0.25 |
| ae_score | 0.15 |
| pathway_score | 0.04 |
| ot_indirect_score | 0.01 |

`fingerprint.min_target_score`: **0.1**

Chosen via sensitivity analysis against the original weights
(0.30/0.20/0.20/0.15/0.15) and one alternative set — see
`docs/adjacency_definition.md` for the full rationale.

## Validation result

**Spearman ρ = 0.6154** across 12 candidate diseases vs. the hand-built gold set.

| Actual rank | Disease | Gold rank | adj_score | Match |
|---|---|---|---|---|
| 1 | lung adenocarcinoma | 2 | 0.5531 | ✗ |
| 2 | squamous cell lung carcinoma | 3 | 0.4434 | ✗ |
| 3 | head and neck squamous cell carcinoma | 5 | 0.3305 | ✗ |
| 4 | urinary bladder carcinoma | 9 | 0.3288 | ✗ |
| 5 | small cell lung carcinoma | 1 | 0.3154 | ✗ |
| 6 | pancreatic adenocarcinoma | 10 | 0.3037 | ✗ |
| 7 | mesothelioma | 4 | 0.3028 | ✗ |
| 8 | ovarian carcinoma | 8 | 0.3002 | ✓ |
| 9 | thyroid gland carcinoma | 11 | 0.2951 | ✗ |
| 10 | breast carcinoma | 7 | 0.2670 | ✗ |
| 11 | colorectal carcinoma | 6 | 0.2313 | ✗ |
| 12 | chronic obstructive pulmonary disease | 12 | 0.0687 | ✓ |

Two exact gold matches: **COPD** (rank 12, correctly the least mechanistically
adjacent — inflammatory rather than oncogenic) and **ovarian carcinoma** (rank 8).

## Output format

Example record (`lung adenocarcinoma`, ranked 1st) from
`output/adjacency_scores.json`. The `shared_*` arrays are large (up to
776 entries) — shown here truncated with counts noted; the real file has
every item.

```json
{
  "disease_a": "non-small cell lung carcinoma",
  "disease_b": "lung adenocarcinoma",
  "adjacency_score": 0.5531,
  "components": {
    "target_score": 0.5444,
    "pathway_score": 0.7951,
    "genetics_score": 0.6786,
    "ot_indirect_score": 0.67,
    "ae_score": 0.3033
  },
  "shared_targets": ["ACVR1B", "AFF3", "ALK", "APC", "AR", "..."],
  "shared_pathways": [
    "ADORA2B mediated anti-inflammatory cytokines production",
    "AKT phosphorylates targets in the cytosol",
    "ALK mutants bind TKIs",
    "..."
  ],
  "shared_genetics": ["ACVR1B", "ATM", "BAZ1A", "BRCA2", "CDKN2A", "..."],
  "shared_adverse_events": [
    "ABDOMINAL DISCOMFORT",
    "ABDOMINAL DISTENSION",
    "ABDOMINAL PAIN",
    "..."
  ]
}
```

| Field | Count for this record |
|---|---|
| shared_targets | 141 |
| shared_pathways | 776 |
| shared_genetics | 19 |
| shared_adverse_events | 91 |

## Known limitations

- **OT indirect score saturation**: `ot_indirect_score` is 0.49 for 10 of 12
  candidates, leaving it with almost no ranking power — full explanation in
  `docs/adjacency_definition.md`.
- **Pathway score saturation**: `pathway_score` sits in a compressed
  0.65–0.80 band from broad, genuinely-shared oncogenic Reactome pathways —
  full explanation in `docs/adjacency_definition.md`.
- **SCLC genetics sparsity**: small cell lung carcinoma's `genetics_score`
  (0.042) likely understates real shared TP53/RB1 biology due to incomplete
  Open Targets GWAS curation — full explanation in `docs/adjacency_definition.md`.

## Data sources

| Source | Provides | URL | Access method |
|---|---|---|---|
| Open Targets Platform | Target associations, pathway membership, genetics evidence, known drugs | https://platform.opentargets.org | GraphQL API |
| EFO via OLS4 (EMBL-EBI Ontology Lookup Service) | Disease name resolution, EFO ontology hierarchy (parents/ancestors) | https://www.ebi.ac.uk/ols4 | REST API |
| openFDA FAERS | Post-market adverse-event reports per drug (MedDRA terms) | https://open.fda.gov/apis/drug/event | REST API |
