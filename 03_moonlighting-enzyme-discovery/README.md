# Moonlighting enzyme discovery — MetaForge / IDEA-B

*IDEA-B Synthetic Biology Competition, MSIT / Korea Synthetic Biology Council · Mar.–Oct. 2026*
*Student team of 4, with an advising researcher and professor*
*My role: candidate selection and validation design*

**Finalist** — one of 15 graduate-division teams from 104 applicants across 37 institutions, in
Korea's first national synthetic biology competition. Carried out as a commissioned research
project at KRIBB.

**In one line.** Cut 14 candidate reaction systems to the 7 that were actually measurable, which
is what let a 4-person team run **~10,000 measurements in 3 days** inside an 8-month budget.

---

## The problem

A moonlighting protein does two unrelated jobs. The same molecule that catalyses a metabolic
reaction may also bind DNA or signal to another cell.

Biology has known this for decades. Nearly every computational tool for enzyme discovery
nonetheless assumes **one protein, one function**, because that assumption makes the data tidy.

## Why it costs money

Multifunctional proteins are heavily over-represented among proteins linked to **drug side
effects**: a drug designed against one activity can disrupt a second activity nobody modelled.
Finding and flagging these proteins is directly useful to anyone developing a drug or an
engineered strain — and single-function models are structurally incapable of it.

## The approach

Rather than force each protein into one category, design an embedding space where a protein with
two functions sits meaningfully **between both**. In practice: metric learning with pairs
defined by overlap between *sets* of reaction classes, so a two-function protein is never split
into two training examples.

Detection is by **cell-free biosensor**, which is what sets the real constraint — the design
below is dictated by what the sensor can actually read.

## My contribution

### 1 · Deciding what to discard

The experimental budget covered eight months. We began with 14 enzyme–substrate systems from
BRENDA, UniProt and the literature, and ran a control experiment on each:

1. Does the biosensor cross-react with the substrate itself?
2. Is it sensitive enough at the concentrations we can reach?

**Half the systems failed one of the two and were dropped before any screening began.**

Cutting half the plan in month two is the only reason the other half finished. It is the part of
this project I would point to in an interview.

### 2 · Designing the screen the survivors made possible

The surviving 7 systems let the team screen **380 metagenomic candidates × 7 substrates ×
4 reaction conditions — roughly 10,000 measurements — in three days** on the biofoundry.

### 3 · Selecting the candidates

The 380 candidates were selected in embedding space, using the same approach as the thesis
project: rank by distance to a reference set rather than by sequence alignment.

---

`code/` holds the database analysis, the embedding model and the target selection pipeline.

**Stack.** ESM-C · metric learning (set-overlap pairs) · BRENDA, UniProt, MoonProtDB ·
cell-free biosensor screening · biofoundry automation
