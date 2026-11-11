"""Train and evaluate the Round 2 candidate ranker.

Repeated StratifiedGroupKFold cross-validation, pooled out-of-fold scoring.
Groups are MMseqs2 30%-identity clusters, so no near-duplicate of a test
sequence is ever in the training fold.

Two things are computed inside the fold, because both would otherwise leak
the test labels:
  1. the four embedding scalars (they are defined relative to known positives)
  2. the near-constant column filter (judged on training rows only)

    python train_rf.py --features features.parquet --emb emb.npy
    python train_rf.py ... --imbalance balanced_subsample --tag bs
    python train_rf.py ... --random-split --tag random   # homology-leak check
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold

# Early-retrieval metrics are defined on fractions, not counts, so they stay
# comparable across evaluation sets of different size (EF@1%, EF@5% convention).
FRACS = [0.05, 0.10]
NEAR_CONST = 0.99
META = {'labels', 'groups', 'grade', 'y', 'n_copy', 'cluster', 'rep', 'aa'}


def precision_at_k(y, s, k):
    k = min(k, len(y))
    return float(y[np.argsort(-s)[:k]].mean())


def metrics(y, s):
    pi = y.mean()
    out = dict(ap=average_precision_score(y, s), roc=roc_auc_score(y, s))
    for f in FRACS:
        k = max(1, int(round(f * len(y))))
        pk = precision_at_k(y, s, k)
        out[f'p_at_{int(f * 100)}'] = pk
        out[f'ef_at_{int(f * 100)}'] = pk / pi
    return out


def fold_features(X, tr, te, y):
    """Embedding scalars built from training-fold labels only."""
    Xn = X / np.linalg.norm(X, axis=1, keepdims=True)
    pos, neg = Xn[tr][y[tr] == 1], Xn[tr][y[tr] == 0]
    cp = pos.mean(0); cp /= np.linalg.norm(cp)
    cn = neg.mean(0); cn /= np.linalg.norm(cn)
    out = {}
    for name, idx in (('tr', tr), ('te', te)):
        A = Xn[idx]
        S = A @ pos.T
        k5 = min(5, S.shape[1])
        out[name] = np.column_stack([
            A @ cp,                                # cos_pos
            S.max(1),                              # knn1_pos
            np.sort(S, axis=1)[:, -k5:].mean(1),   # knn5_pos
            (A @ cp) - (A @ cn),                   # cos_margin
        ])
    return out['tr'], out['te']


def live_cols(F, tr):
    """Drop columns that are near-constant on the training rows."""
    keep = []
    for j in range(F.shape[1]):
        v = F[tr, j]
        v = v[~np.isnan(v)]
        if len(v) and (np.unique(v, return_counts=True)[1].max() / len(v)) <= NEAR_CONST:
            keep.append(j)
    return np.array(keep)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--features', required=True, help='parquet from build_features.py')
    p.add_argument('--emb', required=True, help='npy, row order matches --features')
    p.add_argument('--out-dir', default='.')
    p.add_argument('--repeats', type=int, default=30)
    p.add_argument('--folds', type=int, default=5)
    p.add_argument('--imbalance', default='none', choices=['none', 'balanced_subsample'])
    p.add_argument('--max-features', default='sqrt')
    p.add_argument('--no-emb', action='store_true')
    p.add_argument('--no-cos-lipase', action='store_true')
    p.add_argument('--random-split', action='store_true',
                   help='ignore clusters; measures how large a homology leak would be')
    p.add_argument('--tag', default='main')
    a = p.parse_args()

    d = pd.read_parquet(a.features)
    X = np.load(a.emb)
    y = d.y.to_numpy()
    grp = d.cluster.to_numpy()

    base = [c for c in d.columns if c not in META and pd.api.types.is_numeric_dtype(d[c])]
    if a.no_cos_lipase and 'cos_lipase' in base:
        base.remove('cos_lipase')
    F = np.nan_to_num(d[base].to_numpy(dtype='float64'), nan=-1.0)
    emb_names = ['cos_pos', 'knn1_pos', 'knn5_pos', 'cos_margin']

    print(f'{len(y)} rows / {int(y.sum())} positive ({y.mean():.2%}) / '
          f'{len(set(grp))} clusters')
    print(f'{len(base)} tabular features + {0 if a.no_emb else 4} in-fold embedding scalars')

    rows, oof_last, imp = [], None, []
    for rep in range(a.repeats):
        cv = (StratifiedKFold(a.folds, shuffle=True, random_state=rep) if a.random_split
              else StratifiedGroupKFold(a.folds, shuffle=True, random_state=rep))
        oof = np.full(len(y), np.nan)
        for tr, te in cv.split(F, y, None if a.random_split else grp):
            cols = live_cols(F, tr)
            Ztr, Zte = F[np.ix_(tr, cols)], F[np.ix_(te, cols)]
            names = [base[j] for j in cols]
            if not a.no_emb:
                Etr, Ete = fold_features(X, tr, te, y)
                Ztr, Zte = np.hstack([Ztr, Etr]), np.hstack([Zte, Ete])
                names += emb_names
            kw = dict(n_estimators=2000, min_samples_leaf=2,
                      max_features=a.max_features, n_jobs=-1, random_state=rep)
            if a.imbalance == 'balanced_subsample':
                kw['class_weight'] = 'balanced_subsample'
            m = RandomForestClassifier(**kw).fit(Ztr, y[tr])
            oof[te] = m.predict_proba(Zte)[:, 1]
            if rep == 0:
                imp.append(pd.Series(m.feature_importances_, index=names))
        r = metrics(y, oof)
        r['repeat'] = rep
        rows.append(r)
        oof_last = oof
        print(f'\rrepeat {rep + 1}/{a.repeats}  AP={r["ap"]:.3f} '
              f'EF@5%={r["ef_at_5"]:.1f} EF@10%={r["ef_at_10"]:.1f}', end='', flush=True)
    print()

    res = pd.DataFrame(rows)
    os.makedirs(a.out_dir, exist_ok=True)
    res.to_csv(os.path.join(a.out_dir, f'cv_{a.tag}.csv'), index=False)
    pd.DataFrame({'y': y, 'score': oof_last, 'labels': d.labels, 'cluster': grp}).to_csv(
        os.path.join(a.out_dir, f'oof_{a.tag}.csv'), index=False)
    pd.concat(imp, axis=1).mean(1).sort_values(ascending=False).to_csv(
        os.path.join(a.out_dir, f'mdi_{a.tag}.csv'))

    labs = [('ap', 'PR-AUC'), ('roc', 'ROC-AUC')]
    for f in FRACS:
        labs += [(f'p_at_{int(f*100)}', f'P@{int(f*100)}%'),
                 (f'ef_at_{int(f*100)}', f'EF@{int(f*100)}%')]
    print(f'\n=== {a.tag} ({a.repeats} repeats, '
          f'{"random" if a.random_split else "grouped"} {a.folds}-fold, '
          f'imbalance={a.imbalance}) ===')
    print(f'{"metric":<14}{"mean":>8}{"SD":>8}{"2.5%":>8}{"97.5%":>8}')
    for m_, lab in labs:
        v = res[m_]
        print(f'{lab:<14}{v.mean():>8.3f}{v.std():>8.3f}'
              f'{v.quantile(.025):>8.3f}{v.quantile(.975):>8.3f}')
    print(f'\nbase rate {y.mean():.4f}  ->  random PR-AUC = {y.mean():.4f}, EF = 1.00')
    json.dump({'tag': a.tag, **{k: float(res[k].mean()) for k, _ in labs}},
              open(os.path.join(a.out_dir, f'cv_{a.tag}.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
