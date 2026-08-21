# Defining "adjacent": how this project scores disease similarity

## What adjacency means here

Two indications are "adjacent" to the degree they share the same underlying
biology — the same disease-driving genes, the same pathways, the same
genetic risk factors, the same ontological position in the disease hierarchy,
and similar adverse-event profiles from the drugs used to treat them.

This is distinct from clinical or symptomatic similarity: two diseases can
look unrelated on the surface but share a mechanistic driver that makes a
drug repurposing hypothesis credible. The adjacency score is designed to
surface exactly that kind of hidden biological kinship.

## The five evidence dimensions

| # | Signal | Source | What it captures |
|---|---|---|---|
| 1 | Target overlap | Open Targets target associations | Same gene/protein drives both diseases |
| 2 | Pathway overlap | Open Targets / Reactome | Different genes, same biological pathway |
| 3 | Genetics overlap | Open Targets genetics evidence | Same genetic variants/loci implicated in both |
| 4 | OT indirect / ontology proximity | EFO hierarchy via OLS4 API | How close the two diseases are in the EFO disease tree |
| 5 | Adverse-event overlap | openFDA FAERS + OT knownDrugs | Drugs used in both diseases share a side-effect profile |

Dimensions 1–3 are computed as Jaccard similarity on sets derived from Open
Targets fingerprints. Dimension 4 is a pair-level ontology-distance score
(not Jaccard). Dimension 5 is Jaccard over MedDRA adverse-event term sets
aggregated from FAERS across each disease's known drugs.

## Jaccard similarity

All set-based dimensions use:

```
jaccard(A, B) = |A ∩ B| / |A ∪ B|
```

Returns 0.0 when both sets are empty (no signal, not an error).

## Ontology-proximity score (dimension 4)

```
ot_indirect_score(A, B) = f(shared EFO ancestors, direct parent/child relationship)
```

Where:
- 1.0 = same disease (identical EFO IDs)
- ~0.67 = direct parent/child in EFO (e.g. lung adenocarcinoma ↔ NSCLC)
- ~0.49 = siblings or cousins in EFO (share several ancestors)
- 0.0 = no shared ancestors within 6 hops (or lookup failed — conservative)

This operationalises the Open Targets "indirect associations" concept: OT
propagates target-disease evidence up the EFO hierarchy, so diseases close
in the tree naturally share indirect evidence. We measure that proximity
directly from EFO rather than re-querying OT.

## Adverse-event overlap (dimension 5)

For each disease:
1. Pull its top 15 known drugs from Open Targets (approved and clinical-trial).
2. For each drug, query openFDA FAERS for its top 50 MedDRA adverse-event terms.
3. Take the union of all AE terms as the disease's AE fingerprint.

Jaccard similarity between two diseases' AE fingerprints is the ae_score.

Rationale: diseases treated by drugs with similar side-effect profiles share
underlying biology — the adverse-event space reflects pharmacological
mechanism, not just clinical practice. This adds a phenotypic layer that
pure target/pathway overlap can miss (e.g. two diseases with few shared
targets but the same immune-checkpoint biology will share AE terms like
pneumonitis, colitis, fatigue).

## The combination formula

```
adjacency_score = (0.55 × target_score)
                + (0.04 × pathway_score)
                + (0.25 × genetics_score)
                + (0.01 × ot_indirect_score)
                + (0.15 × ae_score)
```

All weights are set in config/settings.yaml (scoring.weights) and must sum
to 1.0. The formula is a weighted sum; each dimension is bounded [0, 1].

This is the final, empirically-validated weighting ("Set A"). It differs
substantially from the originally-hypothesized weights below -- see
"Why these weights changed" for what the sensitivity analysis found.

## Why these weights (original hypothesis)

The initial weighting, before validation against the gold set, was:

```
target_score: 0.30, pathway_score: 0.20, genetics_score: 0.20,
ot_indirect_score: 0.15, ae_score: 0.15
```

- **Target overlap (0.30)** — the most direct mechanistic evidence. A shared
  gene/protein is a clear, testable repurposing rationale.

- **Pathway overlap (0.20)** — one step more indirect. Shared pathway
  membership means different genes driving the same biological programme,
  which is still a strong mechanistic connection.

- **Genetics overlap (0.20)** — GWAS and genetic-association evidence is an
  independent data type. A genetic signal for both diseases in the same gene
  is strong corroborating evidence for a causal, not correlative, link.

- **OT indirect / ontology proximity (0.15)** — the EFO hierarchy encodes
  curated biomedical knowledge about disease classification. Proximity in
  that ontology is a real prior for shared biology.

- **Adverse-event overlap (0.15)** — adds a phenotypic-pharmacological layer.
  Most valuable when the first three dimensions produce weak signal for
  diseases that are clinically similar but genetically divergent.

## Why these weights changed

Validating the original weights against the hand-built gold set (12
candidates ranked against 2L NSCLC) gave Spearman ρ = 0.5455 -- moderate,
well short of the ρ ≥ 0.70 target. Inspecting the per-dimension component
breakdown across all 12 candidates showed two of the five dimensions had
collapsed to near-constants and were doing no ranking work at all:

- `ot_indirect_score` was 0.49 for 10 of the 12 candidates (everything
  except lung adenocarcinoma's direct-parent-child 0.67) -- nearly every
  candidate sits at the same EFO "cousin" distance from NSCLC, so this
  dimension can't discriminate between them.
- `pathway_score` sat in a narrow 0.65–0.80 band for almost every candidate.
  Solid tumors genuinely share core oncogenic pathway machinery (PI3K/AKT,
  RAS/MAPK, cell-cycle checkpoints) at real, high rates, which compresses
  this dimension's usable range for this candidate set.

Two alternative weight sets were tested (both stripping weight from the
two saturated dimensions down to a nominal 0.04/0.01 floor -- kept nonzero,
not deleted, since the dimensions are still real evidence and should count
for the rare pair where they do differ):

| Weight set | target | pathway | genetics | ot_indirect | ae | ρ |
|---|---|---|---|---|---|---|
| Original | 0.30 | 0.20 | 0.20 | 0.15 | 0.15 | 0.5455 |
| Set A | 0.55 | 0.04 | 0.25 | 0.01 | 0.15 | **0.6154** |
| Set B | 0.45 | 0.04 | 0.30 | 0.01 | 0.20 | 0.4965 |

Set A (target-dominant) won: pushing more weight onto genetics and ae_score
than target_score (Set B) made things worse, since ae_score is also only
weakly discriminating (0.25–0.30 band across nearly all candidates) --
amplifying its weight amplified noise, not signal. Set A is the final
configuration. A `fingerprint.min_target_score` threshold sweep (0.1 → 0.3
→ 0.4 → 0.5) was also tried to fix the pathway saturation at its source, but
made the ranking monotonically worse at every step (ρ: 0.6154 → 0.5385 →
0.3217 → 0.2657) -- see "Known Limitations" for why, and why the threshold
was reverted to its original 0.1.

These weights are still not a fixed law -- if new candidate diseases are
added and the ranking degrades, re-run the sensitivity analysis rather than
assuming Set A generalizes.

## Traceability requirement

Every adjacency score must resolve to specific rows:
- `shared_targets` → the gene symbols that produced the target_score
- `shared_pathways` → the Reactome pathway names that produced the pathway_score
- `shared_genetics` → the gene symbols with genetics evidence in both diseases
- `shared_adverse_events` → the MedDRA terms shared across both AE profiles

The `ot_indirect_score` traces to: the shared EFO ancestor count and whether
a direct parent/child relationship exists, stored in `data/raw/ot_indirect_*.json`.

A score with no attached evidence rows is treated as invalid output.

## Validation approach

The ranked output for the anchor indication (2L NSCLC) against its 12
candidate neighbours is checked against the hand-built gold set
(validation/gold_set.csv) using Spearman rank correlation (scipy.stats.spearmanr).

Target: ρ ≥ 0.70. If the formula's ranking diverges meaningfully from
reviewer judgement, weights are adjusted and re-validated before being
treated as final.

## Validation results (final)

**Spearman ρ = 0.6154** (p = 0.0332) across all 12 candidates, using the
Set A weights and the original min_target_score = 0.1. This is moderate
correlation, below the ρ ≥ 0.70 target -- see "Known Limitations" for why
this is the practical ceiling found, not a tuning shortfall.

Full ranked comparison:

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

Two exact rank matches: COPD (correctly ranked least-adjacent, consistent
with it being inflammatory rather than oncogenic) and ovarian carcinoma.
The largest misses: small cell lung carcinoma (gold #1, actual #5) and
urinary bladder carcinoma (gold #9, actual #4) -- see below.

## Known limitations

- **OT indirect score saturation**: `ot_indirect_score` is 0.49 for 10 of
  the 12 candidates (everything except lung adenocarcinoma, which has a
  direct EFO parent/child relationship at 0.67) -- nearly every candidate
  sits at the same "EFO cousin" ontology distance from NSCLC, so this
  dimension carries almost no ranking information for this candidate set.

- **Pathway score saturation**: `pathway_score` sits in a compressed
  0.65–0.80 band for nearly every candidate (COPD's 0.325 is the only real
  outlier), because solid tumors genuinely share core oncogenic pathway
  machinery (PI3K/AKT, RAS/MAPK, cell-cycle checkpoints) at real, high
  rates -- raising `fingerprint.min_target_score` to filter this out was
  tried (0.1 → 0.5) and made the ranking monotonically worse at every step,
  because Open Targets' evidence density varies enormously by disease
  (NSCLC has 200/200 fetched targets scoring ≥0.43; small cell lung
  carcinoma and mesothelioma have far fewer high-scoring targets), so a
  global score threshold penalizes less-studied diseases regardless of
  true biological similarity.

- **SCLC genetics sparsity**: small cell lung carcinoma's `genetics_score`
  is 0.042 despite the gold set's rationale citing shared TP53/RB1 genetic
  drivers with NSCLC -- likely reflecting incomplete GWAS/genetic-association
  curation for SCLC in Open Targets rather than a true absence of shared
  genetic evidence, and it's the main reason SCLC under-ranks (gold #1,
  actual #5).
