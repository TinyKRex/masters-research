# Checklist before making this repository public

Work through this once. Everything on it is easier to do before the first push than after —
a public commit is retrievable from the git history and from third-party mirrors even if the
file is deleted in a later commit.

## 1. Permission

- [ ] **Advisor (Prof. Haseong Kim) has approved public release.** The underlying data belongs
      to an ongoing KRIBB project and the thesis has not been defended. This is not optional.
- [ ] Checked whether any result is intended for a patent filing. **Public disclosure of a hit
      sequence can destroy its novelty**, which in most jurisdictions cannot be undone.
- [ ] Checked whether a manuscript in preparation has a journal embargo that this would breach.
- [ ] Confirmed which code is lab-owned (`sblabkribb/*`) versus personally authored, and that the
      lab is content for the lab-owned portion to carry an MIT licence.

## 2. What must not be in the repository

- [ ] No metagenome assemblies, ORF sequences, or FASTA files of any kind
- [ ] No candidate lists with sequences attached — including the Round 1 855 and the Round 2 380
- [ ] No sequences of confirmed hits
- [ ] No raw plate readings that could be joined back to sequences
- [ ] No internal KRIBB documents, costings, supplier pricing, or unpublished protocols the lab
      treats as confidential
- [ ] No API keys, tokens, `.env` files, certificates (`huggingface.crt` and similar), or
      absolute paths containing personal directory names
- [ ] No email addresses — contact routes through the CV site only
- [ ] No photographs of lab members or of equipment with visible asset or serial numbers

Run before the first push:

```bash
cd ~/github/masters-research
grep -rinE '(api[_-]?key|secret|token|password|BEGIN [A-Z ]*PRIVATE KEY)' . --exclude-dir=.git
grep -rnE '[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}' . --exclude-dir=.git
grep -rn '/home/\|/Users/\|C:\\\\' . --exclude-dir=.git
grep -rlE '^>[A-Za-z0-9]' . --exclude-dir=.git          # stray FASTA
find . -path ./.git -prune -o -size +5M -print          # anything large enough to be data
```

Notebooks keep their outputs. Strip them:

```bash
pip install nbstripout && find . -name '*.ipynb' -exec nbstripout {} +
```

## 3. Make the history clean from the start

- [ ] This is a **new repository**, not a copy of `YK_project` or `enzyme_classification` with
      files deleted. Copy in the files you want; do not import the old history.
- [ ] `git log` shows only commits you intend to publish.
- [ ] `git config user.email` is a GitHub `noreply` address or a public address, not the
      institutional one — commit metadata is permanently public.

```bash
git -C ~/github/masters-research config user.email "<your-github-noreply>@users.noreply.github.com"
```

## 4. After publishing

- [ ] Repository description and topics set (`protein-language-models`, `metagenomics`,
      `enzyme-engineering`, `cell-free-protein-synthesis`, `biofoundry`)
- [ ] Zenodo linked to the repository, a `v1.0` release tagged, DOI minted
- [ ] DOI written into `CITATION.cff`, the top-level `README.md`, and the thesis bibliography
- [ ] Linked from the CV site (`projects.md` and `research.md`)
- [ ] Re-read the top-level README as a stranger: does someone outside the field learn what the
      work was for within thirty seconds?
