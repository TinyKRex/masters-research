"""Assemble the training feature table.

Structural and sequence features come from the annotation pipeline; raw
embeddings are kept out of the table. Only `cos_lipase` is computed here,
because it is the one embedding feature that does not use labels and can
therefore be computed outside the cross-validation fold. The label-dependent
scalars are built inside the fold by train_rf.py.

Excluded on purpose:
  nucleotide composition   computed differently in the screening pool, so the
                           distributions do not line up
  docking scores           100% missing in the pool

    python build_features.py --dataset dataset.csv --clusters clusters.csv \
        --annotations annotations.parquet --emb 'prott5_*.h5' \
        --ref-emb ref_prott5.h5 --pool-fasta pool.fasta --out-dir .
"""
import argparse
import os

import numpy as np
import pandas as pd

from emb_io import load_emb, read_fasta_map

DROP = ['gc', 'gc3', 'gc_head45', 'rare_codon_frac'] + [
    f'dock_{s}_{m}' for s in ('pNPB', 'pNPO', 'pNPP')
    for m in ('best_productive_affinity', 'n_productive')]
META = {'candidate_label', 'group', 'grade', 'weight'}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--dataset', required=True, help='one row per unique sequence, with y')
    p.add_argument('--clusters', required=True, help='MMseqs2 30% identity assignment')
    p.add_argument('--annotations', required=True, help='structure/sequence feature table')
    p.add_argument('--emb', required=True)
    p.add_argument('--ref-emb', required=True, help='embeddings of characterised lipases')
    p.add_argument('--pool-fasta', required=True)
    p.add_argument('--out-dir', default='.')
    a = p.parse_args()

    ds = pd.read_csv(a.dataset)
    ds['cluster'] = pd.read_csv(a.clusters).cluster
    ds['rep'] = ds.labels.str.split(';').str[0]

    f = pd.read_parquet(a.annotations)
    cols = [c for c in f.columns
            if c not in META | set(DROP) and not c.startswith(('prott5_', 'esmc_'))]
    f = f[['candidate_label'] + cols].drop_duplicates('candidate_label')

    ds = ds.drop(columns=[c for c in ds.columns if c in set(cols)])
    d = ds.merge(f.rename(columns={'candidate_label': 'rep'}), on='rep', how='left')
    miss = d[cols].isna().all(1).sum()
    if miss:
        print(f'warning: {miss} sequences have no structural features')

    idx, emb = load_emb(a.emb)
    _, ref = load_emb(a.ref_emb)
    u = read_fasta_map(a.pool_fasta).set_index('aa').uniq_idx.reindex(d.aa)
    if u.isna().any():
        raise SystemExit(f'{int(u.isna().sum())} sequences are not in the pool fasta')
    X = emb[np.searchsorted(idx, u.to_numpy().astype('int64'))]

    c = ref.mean(0)
    d['cos_lipase'] = (X @ c) / (np.linalg.norm(X, axis=1) * np.linalg.norm(c))

    os.makedirs(a.out_dir, exist_ok=True)
    out = os.path.join(a.out_dir, 'features.parquet')
    d.drop(columns=['aa']).to_parquet(out, index=False)
    np.save(os.path.join(a.out_dir, 'emb.npy'), X)

    num = [c_ for c_ in cols if pd.api.types.is_numeric_dtype(d[c_])]
    dead = [c_ for c_ in num
            if d[c_].value_counts(normalize=True, dropna=False).iloc[0] > 0.99]
    print(f'{len(d)} rows / {len(num)} tabular features (+ cos_lipase)')
    print(f'{len(dead)} are near-constant here and will be dropped inside each fold')
    print(f'-> {out}, {os.path.join(a.out_dir, "emb.npy")} {X.shape}')


if __name__ == '__main__':
    main()
