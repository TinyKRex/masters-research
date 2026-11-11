# Metagenomic enzyme screening with an iterative learning loop

*M.S. thesis · KAIST Graduate School of Engineering Biology · KRIBB Synthetic Biology Research
Center · 2024–2027*

**In one line.** Ranked 5.3 million environmental genes with a protein language model, tested
855 of them on an automated biofoundry, and used the results — failures included — to train a
selector that concentrates working enzymes **8.6× more densely** than random into the top 5% of
its ranking.

---

## Part 1 — for a reader outside the field

### The problem

Soil and seawater are full of microbes nobody can grow in a lab. Their genes can still be read,
and among them are enzymes — catalysts that let industry replace steps needing heat, solvents or
heavy metals. One environmental sample yields millions of candidate genes.

Testing one candidate costs money and days of work. Testing millions is impossible. **The
shortlist is the whole problem.**

### How shortlists are made today, and why it limits you

By sequence similarity: find genes that resemble enzymes we already know. Cheap, reliable, and
by construction it returns more of what we already have. The genuinely unfamiliar enzymes — the
commercially interesting ones — are exactly what it misses.

### What I changed

**1 · Rank by meaning, not by spelling.**
A protein language model turns each gene into a vector encoding structure and chemistry, not
just letters. Two enzymes doing the same job can land close together there even when their
sequences barely resemble each other — so the ranking reaches candidates similarity search
cannot return.

**2 · Measure at scale, under identical conditions.**
855 candidates were built as DNA and expressed in a cell-free system — protein made in a tube
rather than inside a living cell — then assayed on an automated biofoundry. Automation here is
not about speed. It is so that a round run in March is comparable to one run in November.

**3 · Feed the failures back.**
The candidates that did *not* work are the valuable part. Similarity can tell you what looks
promising; only an experiment tells you what looks promising **and fails**. Those confirmed
negatives are what make a learned ranking possible at all.

### Where it stands, stated honestly

On Round 1 data the trained model reaches **8.6× enrichment** in the top 5%, against 7.0× for
the standard similarity method.

But the confidence interval on that improvement **includes zero** (−2.6 to +7.4). Scored
retrospectively, I cannot claim the model beats the baseline. So Round 2 settles it the only way
that counts: **190 candidates chosen by the model and 190 by similarity alone, no overlap, run
side by side in the same experiment.** The predicted outcome was recorded in writing before the
experiment began.

### Why the honest version is the interesting one

A great deal of machine learning for biology reports how well a model ranks a dataset and stops
there. The number that decides whether a discovery programme is affordable is the hit rate in
the *next* experiment — and almost nobody reports it, because it means running the experiment
and publishing the result when it disappoints. This project is built to report it either way.

---

## Part 2 — technical

### Pipeline

```
16 metagenome samples
  └─ ORF prediction ──────────────────────────────► ~5.3 M ORFs
       └─ ProtT5 embeddings (1024-d)
            └─ cosine similarity to the mean embedding of 98 characterised lipases
                 └─ structural filter: ESMFold pLDDT ≥ 70, Foldseek fold match,
                    FoldDisco catalytic-triad geometry
                      └─ Round 1: 855 ORFs across 9 libraries
                           └─ tailing-PCR cassettes → in-house GroEL/ES S30 cell-free
                              expression → pNP-ester assay at 405 nm
                                └─ labels ──► random forest selector
                                     └─ Round 2: 190 model vs 190 similarity
```

### Round 1 outcome

| | |
|---|---|
| ORFs ordered | 855 |
| Measured after construction and QC losses | 703 |
| Unique sequences after merging duplicates | 637 |
| Confirmed active | 17 (2.67%) |
| Distinct 30%-identity clusters among the actives | 17 of 17 |

Every active sequence falls in a different cluster, so the positives are not one family
rediscovered seventeen times.

### Selector and evaluation

Random forest over 85 features — structural annotations, sequence descriptors and embedding
similarity scores — evaluated by **StratifiedGroupKFold** in which whole 30%-identity clusters
are held out, repeated 30×, scored on pooled out-of-fold predictions.

The four label-dependent embedding scalars and the near-constant column filter are both computed
**inside each fold**. Either one computed outside leaks test labels.

| Ranking criterion | EF@5% | EF@10% |
|---|---|---|
| **Random forest** | **8.63** | **5.87** |
| Embedding cosine alone | 7.03 | 4.68 |
| Structural similarity alone | 5.85 | 5.27 |
| Catalytic triad present (rule) | 2.50 | — |
| Predicted fold confidence (negative control) | 2.34 | 1.17 |

### Caveats, as stated in the thesis

- Paired cluster-level bootstrap: the model's advantage over cosine is **+2.0 EF points at 5%,
  95% CI −2.6 to +7.4**. Not established retrospectively. Hence the prospective trial.
- Permutation test confirms the ranking is non-random (*p* = 0.003).
- A *random* rather than cluster-aware split scores **lower**, which argues against the result
  being homology leakage.
- Predicted fold confidence alone reaches 2.34×. Part of what the assay measures is therefore
  whether a protein expresses and folds in this particular cell-free system, not only whether it
  is catalytically correct. That is a property of activity-based screening generally, and it is
  reported rather than hidden.

### Prospective prediction, recorded 2026-11-11

Written down before Round 2 was run. Not a formal preregistration.

- Model arm (190 candidates): **30–53 confirmed actives**, central estimate 38
- Similarity arm (190 candidates): **≥ 5 confirmed actives**

Both assume the background activity rate inside the applicability domain stays near 2.67%, which
may be optimistic.

---

## Contents

See [`code/`](code/) for the run order and the input schema. Metagenome assemblies, ORF
sequences, candidate lists and hit sequences are not included — they are part of ongoing work at
KRIBB.

**Stack.** ProtT5, ESM-C · ESMFold · Foldseek · FoldDisco · MMseqs2 · scikit-learn ·
Echo 525 acoustic dispenser and automated liquid handlers
