"""Paired comparison of the trained ranker against the similarity baselines.

Resampling is done over 30%-identity clusters, not over rows: sequences in the
same cluster are not independent observations, and bootstrapping rows would
produce a confidence interval that is too narrow to mean anything.
"""
import argparse

import numpy as np
import pandas as pd
from scipy.stats import rankdata

p = argparse.ArgumentParser()
p.add_argument('--features', required=True)
p.add_argument('--oof', required=True, help='npy of shape (repeats, n_rows)')
p.add_argument('--boot', type=int, default=2000)
a = p.parse_args()

O = np.load(a.oof)
f = pd.read_parquet(a.features)
y = f.y.to_numpy().astype(int)
cos = f.cos_lipase.to_numpy()
tm = np.nan_to_num(f.ref_tm.to_numpy(), nan=-1.0)
cl = f.cluster.to_numpy()


def ef(score, yy, frac):
    k = max(1, int(round(frac * len(yy))))
    return yy[np.argsort(-score)[:k]].mean() / yy.mean()


print(f'n={len(y)}  P={y.sum()}  pi={y.mean():.4f}\n')
for frac in (0.05, 0.10):
    e = np.array([ef(O[r], y, frac) for r in range(O.shape[0])])
    ec, et = ef(cos, y, frac), ef(tm, y, frac)
    print(f'EF@{int(frac*100)}%  RF {e.mean():.2f} (min {e.min():.2f}, max {e.max():.2f})  '
          f'cosine {ec:.2f}  structure {et:.2f}  |  '
          f'RF>cosine {int((e > ec).sum())}/{O.shape[0]}, '
          f'RF>structure {int((e > et).sum())}/{O.shape[0]}')

rfa = np.mean([rankdata(O[r]) for r in range(O.shape[0])], axis=0)
print(f'\nrank-averaged RF  EF@5% {ef(rfa, y, .05):.2f}  EF@10% {ef(rfa, y, .10):.2f}')

rng = np.random.default_rng(0)
groups = [np.where(cl == c)[0] for c in np.unique(cl)]
for frac in (0.05, 0.10):
    D = []
    for _ in range(a.boot):
        idx = np.concatenate([groups[i] for i in
                              rng.integers(0, len(groups), len(groups))])
        yy = y[idx]
        if yy.sum() < 1:
            continue
        D.append(ef(rfa[idx], yy, frac) - ef(cos[idx], yy, frac))
    D = np.array(D)
    lo, hi = np.percentile(D, [2.5, 97.5])
    print(f'dEF@{int(frac*100)}% (RF - cosine)  mean {D.mean():.2f}  '
          f'95% CI {lo:.2f} to {hi:.2f}  P(RF>cosine) {(D > 0).mean():.3f}')
