# Open items

Last updated 2026-11-12.

## Before the thesis is submitted

- [ ] **Set repository topics.** The API still reports none. About ⚙️ → Topics:
      `protein-language-models`, `metagenomics`, `enzyme-engineering`,
      `cell-free-protein-synthesis`, `biofoundry`, `machine-learning`
- [ ] **Fill `code/` for projects 02, 03 and 04.** Only `01_` is populated. Strip notebook
      outputs first: `pip install nbstripout && find . -name '*.ipynb' -exec nbstripout {} +`
- [ ] **Add figures.** Every `figures/` directory is empty.
- [ ] **Release `v1.1.0`** once the above is done, and update the DOI in the thesis
      bibliography. The `v1.0.0` Zenodo record carries an earlier form of the thesis title;
      a new release corrects it. The concept DOI does not change.

## After Round 2 (experiment 2026-11-18 – 11-22)

- [ ] Replace the prospective prediction in `01_metagenomic-enzyme-screening/README.md` with
      the measured result — whichever way it comes out. The prediction was recorded on
      2026-11-11 precisely so that it could be checked.

## Standing rules

- No sequence data, candidate lists or hit sequences in this repository.
- Work through `PUBLICATION_CHECKLIST.md` before any release.
