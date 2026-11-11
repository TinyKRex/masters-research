# Enzyme function prediction from sequence

*Groundwork for the thesis · 2024–2025*

**In one line.** Established that a pretrained protein model's internal representation encodes
**function**, not just sequence resemblance — the assumption the entire screening strategy rests
on — by showing a small head on frozen embeddings matches end-to-end sequence models at a
fraction of the compute.

---

## The problem

Every enzyme has a four-part code describing the reaction it catalyses. For most proteins ever
sequenced, nobody has measured it. Can it be predicted from the amino-acid sequence alone?

## Why I needed to know

Before building a screening strategy on top of a protein language model, I needed evidence that
the model's representation carries functional information rather than only sequence similarity.
If it does, a small classifier on top of a **frozen** model should recover reaction classes it
was never explicitly taught.

That is a falsifiable precondition, and it is why this came first.

## What I did

- Assembled a labelled dataset from UniProt/SwissProt (curated) and TrEMBL (automatic), with the
  reaction hierarchy flattened to the level where class sizes stayed usable
- Trained CNN, LSTM, GRU and transformer models directly on sequence
- Compared those against a small head on frozen ESM-3 / ESM-C embeddings

## Result

The embedding-based models reached comparable or better accuracy at a fraction of the training
cost. **The expensive part of understanding a protein has already been paid for by the
pretrained model** — which is the finding the thesis builds on.

## Related result

Collaborative study, **KMB 2025**: single-token embeddings — normally discarded when a protein
is summarised by averaging over its residues — carry global information about the whole protein,
improving bacterial essential-gene prediction by **up to 17%** over average pooling.

---

`code/` holds data collection, preprocessing, embedding generation and training scripts.
Trained weights and full datasets are not included; both are reconstructible from the scripts
and the UniProt releases they name.

**Stack.** PyTorch · ESM-3 / ESM-C · UniProt / SwissProt / TrEMBL
