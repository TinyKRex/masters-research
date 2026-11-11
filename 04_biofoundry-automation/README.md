# Biofoundry automation and assay standardisation

*Infrastructure behind the screening projects · KRIBB Synthetic Biology Research Center*

**In one line.** Replaced a commercial reagent with an in-house preparation at **~2.5× the
expression yield**, and standardised the screening workflow so that rounds run months apart are
comparable to each other.

---

## Why this exists

An iterative strategy only means something if each round is measured the same way as the last.
If Round 2 looks better because the reagents were fresher, nothing has been learned. Most of
this work is removing those excuses before they can be made.

## What was built

**Library construction.** Linear expression cassettes built by tailing PCR, carrying the
promoter and terminator in the primer tails. This replaced a longer route through purification,
quantification and assembly — fewer steps, fewer places to lose material, and a per-candidate
cost low enough to order hundreds at a time. A gradient PCR resolved the melting-temperature gap
between the two anchors.

**Cell-free extract, made in-house.** Proteins are expressed in a tube rather than in a living
cell, which removes the need to keep hundreds of candidate strains alive. I prepared S30 extract
from a chaperone-overexpressing (GroEL/GroES) strain instead of buying a kit:

- **~2.5× the expression yield** by GFP measurement
- Removed batch-to-batch variability from a supplier we do not control
- Detected activity on three substrates in the positive control **that the kit had missed** —
  meaning the kit had been scoring working enzymes as failures

**Automation.** An Echo 525 acoustic dispenser and automated liquid handlers run the 96- and
384-well expression, dilution and assay steps, with an sfGFP expression control on every plate
so plate-to-plate variation is measured rather than assumed away.

**Hit calling without a fixed threshold.** Each plate's own background defines a robust Z-score,
and a candidate is called only if it clears both Z ≥ 3 and that plate's positive control. The
threshold moves with the plate instead of being chosen once and defended forever.

## Related tool

`echoprep` — turns a plate layout into the pick-lists the dispenser and liquid handlers consume,
which is the step where manual transcription errors used to enter.

---

`code/` holds the protocol scripts, plate layout generation and assay analysis.

**Stack.** Echo 525 acoustic dispenser · automated liquid handlers · cell-free protein synthesis
(S30) · pNP-ester activity assays · Python
