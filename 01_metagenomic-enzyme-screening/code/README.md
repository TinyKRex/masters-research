# Code

Run order:

```bash
python build_features.py   --dataset … --clusters … --annotations … \
                           --emb … --ref-emb … --pool-fasta …
python train_rf.py         --features features.parquet --emb emb.npy \
                           --imbalance balanced_subsample --tag bs
python evaluate_baselines.py     --features features.parquet
python compare_rf_vs_baseline.py --features features.parquet --oof oof30.npy
python select_cosine.py    …    # similarity arm, 190 candidates
python select_rf.py        …    # model arm, 190 candidates
```

Every path is a command-line argument; nothing is hard-coded.

## Inputs

| File | Contents |
|---|---|
| `dataset.csv` | one row per unique sequence: `labels`, `aa`, `y` (0/1 activity), `n_copy` |
| `clusters.csv` | `cluster` — MMseqs2 assignment at 30% identity, row-aligned to `dataset.csv` |
| `annotations.parquet` | structural and sequence features keyed by `candidate_label` |
| `*.h5` | ProtT5 embeddings, indexed by `uniq_idx` |
| `ref_prott5.h5` | embeddings of the 98 characterised reference lipases |
| `pool.fasta` | the screening pool, headers carrying `uniq_idx` |
| `orf_index.parquet` | `uniq_idx` → source metagenome, contig, length |
| `pool_features.parquet` | the same feature columns as `annotations.parquet`, for the pool |
| `oof30.npy` | `(repeats, n_rows)` out-of-fold scores written by `train_rf.py` |

Sequence data is not included in this repository. The schema above is what the
scripts need; they run on any dataset in that shape.

## External helpers

`emb_io` and `design_primers` belong to the lab's shared pipeline and are not
redistributed here. Their interfaces:

```python
load_emb(pattern)       -> (idx: int64[n], emb: float32[n, 1024])  # idx sorted
read_fasta_map(path)    -> DataFrame[uniq_idx, aa]
load_cds(dir, groups)   -> {group: {seq_id: (nt_sequence, description)}}
pick(nt)                -> (primer, tm, warning)                   # 5' anchor
rc(nt)                  -> reverse complement
```

## Requirements

Python 3.10+, `numpy`, `pandas`, `pyarrow`, `scikit-learn`, `scipy`, `h5py`.

A full 30-repeat run of `train_rf.py` takes about 20 minutes on 16 cores.
