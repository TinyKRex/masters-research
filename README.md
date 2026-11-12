# Finding useful enzymes with fewer experiments

**Yunkyeong Kang** — M.S. Engineering Biology, KAIST · Synthetic Biology Research Center, KRIBB
[tinykrex.github.io](https://tinykrex.github.io/) · 2024–2027
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24149119.svg)](https://doi.org/10.5281/zenodo.24149119)

Four projects from my master's work. Each folder has a plain-language summary, the core code,
and the numbers — including the ones that did not come out the way I wanted.

---

## The problem, in business terms

Enzymes replace chemical process steps that need heat, solvents or heavy metals. Finding a new
one is a search problem with an expensive test: a single environmental DNA sample yields
**millions of candidate genes**, and each candidate costs money and days of lab work to verify.

Nobody runs out of candidates. They run out of budget for testing.

So the metric that decides whether a discovery programme is viable is the **hit rate** — how many
of the candidates you test actually work. Going from 1 hit per 100 tested to 3 per 100 cuts the
cost of a confirmed enzyme by two thirds. Everything below is in service of moving that number.

## What I built

| | |
|---|---|
| **Scale** | Ranked 5.3 M metagenomic genes using a protein language model, then tested 855 of them in the lab |
| **Platform** | Automated cell-free expression and assay on a biofoundry; in-house reagent at **2.5× the yield** of the commercial kit |
| **Result** | A ranking model that concentrates working enzymes **8.6× more densely** than random into the top 5% of its list |
| **Validation** | A 380-candidate prospective trial comparing the model head-to-head against the standard method, with the predicted outcome recorded before the experiment |

The part I would defend in an interview is the last row. A model scored on its own training data
is a claim about the data; a hit rate in the next experiment is a claim about the world. Most
published work in this area reports the first and calls it the second.

---

## Projects

### [01 · Metagenomic enzyme screening with an iterative learning loop](01_metagenomic-enzyme-screening/)
*M.S. thesis · Python, scikit-learn, protein language models, biofoundry automation*

Built an end-to-end pipeline from 5.3 M environmental genes to measured enzyme activity, and
closed the loop by training the next round's selector on the previous round's results —
**including the failures**, which similarity search cannot provide and which are what make a
learned ranking possible.

- **8.6× enrichment** at the top 5% of the ranking, against 7.0× for the standard
  similarity method and 2.3× for a negative control
- Validated with cluster-aware cross-validation and a cluster-level bootstrap, so the
  confidence interval reflects how correlated biological sequences really are
- Reported honestly: the advantage over the standard method has a 95% CI of **−2.6 to +7.4**.
  That is why Round 2 is a prospective trial and not a claim

### [02 · Enzyme function prediction from sequence](02_enzyme-ec-classification/)
*Groundwork · PyTorch, CNN / LSTM / GRU / transformer, ESM embeddings*

Tested whether a pretrained protein model's internal representation encodes **function**, not
just sequence resemblance — the assumption the whole screening strategy rests on. A small head
on frozen embeddings matched or beat models trained end-to-end on sequence, at a fraction of the
compute.

Related result (KMB 2025, collaborative): embeddings normally discarded during pooling carry
whole-protein information, improving bacterial essential-gene prediction by **up to 17%**.

### [03 · Moonlighting enzyme discovery — MetaForge / IDEA-B](03_moonlighting-enzyme-discovery/)
*Competition project, funded at KRIBB · team of 6, my role: candidate selection and validation design*

Some proteins do two unrelated jobs, and they are over-represented among proteins implicated in
**drug side effects** — a commercial risk that single-function models cannot see.

My contribution was scope control. Of 14 enzyme–substrate systems, control experiments showed
**half were unmeasurable** with our biosensor. Cutting them before screening started is what let
the team run **380 candidates × 7 substrates × 4 conditions — ~10,000 measurements — in 3 days**
inside an 8-month budget.

### [04 · Biofoundry automation and assay standardisation](04_biofoundry-automation/)
*Infrastructure · acoustic liquid handling, cell-free protein synthesis, plate-level QC*

Made rounds run months apart comparable to each other, which is the precondition for any
iterative strategy meaning anything.

- In-house cell-free extract at **~2.5× the expression yield** of the commercial kit, and
  detecting activity on three substrates the kit had scored as failures
- Automated 96/384-well expression and assay with an expression control on every plate
- Hit calling from each plate's own background rather than a fixed threshold

---

## What this says about how I work

- **I pick the metric that costs money**, not the one that is easy to improve. Hit rate per
  experiment, not model accuracy.
- **I design the test that could prove me wrong.** The Round 2 trial is built so that a null
  result is publishable and interpretable.
- **I cut scope early.** Half the IDEA-B reaction systems were dropped in month two, which is
  the only reason the rest finished.
- **I report the confidence interval.** Including when it contains zero.

## Repository layout

```
<project>/README.md     plain-language summary, then the technical detail
<project>/code/         the scripts behind the result, paths parameterised
<project>/figures/
```

Sequence data, candidate lists and unpublished hit sequences are **not included** — they are part
of ongoing work at KRIBB. Each script documents the input schema it expects, so the methods run
on other data.

## Citing

> Kang, Y. (2027). *Enhancing metagenomic enzyme screening efficiency through an embedding-based
> iterative learning strategy* [M.S. thesis]. Korea Advanced Institute of Science and Technology.

This repository is archived at **[10.5281/zenodo.24149119](https://doi.org/10.5281/zenodo.24149119)** — a concept DOI that always resolves to the
latest version.

## Licence

Code: MIT. Text and figures: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
Contact via [my site](https://tinykrex.github.io/).
