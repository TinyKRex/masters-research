"""Model arm: 190 Round 2 candidates from the trained ranker.

Applicability domain, applied before the model is allowed to score anything:

  cosine rank <= max rank of any training row
      The training set was itself selected by cosine, so it occupies a bounded
      band of that axis. Outside it the model extrapolates, and because its
      score correlates negatively with cosine, widening the band makes it pick
      the least lipase-like sequences available. The cut is read off the
      training rows, not chosen.

  pLDDT >= 70
      The lowest pLDDT in the training set is 71.70. Below that the model has
      seen nothing. This is the edge of the training region, not a quality bar.

Excluded: sequences identical to Round 1 ORFs, and the 190 already taken by the
similarity arm. The overlap between the two arms is reported before removal —
it measures how much the two methods agree, which is part of the result.

    python select_rf.py --features features.parquet --emb emb.npy \
        --pool-table pool_features.parquet --dataset dataset.csv \
        --cds-dir CDS/ --emb-h5 'prott5_*.h5' --ref-emb ref_prott5.h5 \
        --orf-index orf_index.parquet --pool-fasta pool.fasta \
        --round1-aa round1.txt --cosine-arm round2_cosine190.csv
"""
import argparse
import os

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier

import train_rf as T
from design_primers import load_cds, pick, rc
from emb_io import load_emb, read_fasta_map

MIN_TM, MAX_TM = 55.0, 63.0
MIN_PLDDT = 70.0
N_ARM = 190
SCAN = 1200

# With 17 positives against 620 negatives, an unweighted forest grows many trees
# that never see a positive. Paired on the same split seeds, balanced_subsample
# wins in all 30 repeats: PR-AUC 0.127 -> 0.207, EF@5% 6.01 -> 8.63.
RF_KW = dict(n_estimators=2000, min_samples_leaf=2, max_features='sqrt',
             class_weight='balanced_subsample', n_jobs=-1)


def emb_scalars(Xq, Xtr, ytr, ref):
    """Same definition as train_rf.fold_features, fitted on all training rows."""
    Qn = Xq / np.linalg.norm(Xq, axis=1, keepdims=True)
    Tn = Xtr / np.linalg.norm(Xtr, axis=1, keepdims=True)
    pos, neg = Tn[ytr == 1], Tn[ytr == 0]
    cp = pos.mean(0); cp /= np.linalg.norm(cp)
    cn = neg.mean(0); cn /= np.linalg.norm(cn)
    cl = ref.mean(0)
    S = Qn @ pos.T
    k5 = min(5, S.shape[1])
    return np.column_stack([
        (Xq @ cl) / (np.linalg.norm(Xq, axis=1) * np.linalg.norm(cl)),
        Qn @ cp, S.max(1), np.sort(S, axis=1)[:, -k5:].mean(1), (Qn @ cp) - (Qn @ cn),
    ])


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--features', required=True)
    p.add_argument('--emb', required=True)
    p.add_argument('--pool-table', required=True)
    p.add_argument('--dataset', required=True)
    p.add_argument('--cds-dir', required=True)
    p.add_argument('--emb-h5', required=True)
    p.add_argument('--ref-emb', required=True)
    p.add_argument('--orf-index', required=True)
    p.add_argument('--pool-fasta', required=True)
    p.add_argument('--round1-aa', required=True)
    p.add_argument('--cosine-arm')
    p.add_argument('--out-dir', default='.')
    p.add_argument('--n', type=int, default=N_ARM)
    p.add_argument('--cos-rank-cut', type=int, default=None)
    p.add_argument('--allow-cosine-overlap', action='store_true')
    a = p.parse_args()

    # ---- fit ----
    d = pd.read_parquet(a.features)
    Xtr = np.load(a.emb)
    y = d.y.to_numpy()
    # cos_lipase is recomputed by emb_scalars(); keeping the column as well would
    # break the lookup, since the pool table does not carry it.
    base = [c for c in d.columns if c not in T.META and c != 'cos_lipase'
            and pd.api.types.is_numeric_dtype(d[c])]
    Ftr = np.nan_to_num(d[base].to_numpy(dtype='float64'), nan=-1.0)
    cols = T.live_cols(Ftr, np.arange(len(y)))
    names = [base[j] for j in cols]
    print(f'training on {len(y)} rows / {int(y.sum())} positive / '
          f'{len(cols)} features + 5 embedding scalars')

    idx, emb = load_emb(a.emb_h5)
    _, ref = load_emb(a.ref_emb)
    m = RandomForestClassifier(**RF_KW, random_state=0).fit(
        np.hstack([Ftr[:, cols], emb_scalars(Xtr, Xtr, y, ref)]), y)

    # ---- applicability domain ----
    pt = pd.read_parquet(a.pool_table)
    n0 = len(pt)
    pt = pt[pt.plddt_mean >= MIN_PLDDT].reset_index(drop=True)
    print(f'pLDDT >= {MIN_PLDDT:g}: {len(pt):,}/{n0:,} kept')

    fa = read_fasta_map(a.pool_fasta)
    c = np.vstack([ref, Xtr[y == 1]]).mean(0)
    cos_all = (emb @ c) / (np.linalg.norm(emb, axis=1) * np.linalg.norm(c))
    cos_rank = pd.Series(
        np.argsort(np.argsort(-cos_all, kind='stable'), kind='stable') + 1, index=idx)
    tr_rank = cos_rank.reindex(
        fa[fa.aa.isin(set(pd.read_csv(a.dataset).aa))].uniq_idx).dropna()
    cut = int(tr_rank.max()) if a.cos_rank_cut is None else a.cos_rank_cut
    n0 = len(pt)
    pt = pt[pt.uniq_idx.map(cos_rank) <= cut].reset_index(drop=True)
    print(f'cosine rank <= {cut:,} (max over {len(tr_rank)} training rows): '
          f'{len(pt):,}/{n0:,} kept')

    # ---- score ----
    missing = [c_ for c_ in names if c_ not in pt.columns]
    if missing:
        raise SystemExit(f'{len(missing)} feature columns missing from the pool table: '
                         f'{missing[:10]}')
    Fp = np.nan_to_num(pt[names].to_numpy(dtype='float64'), nan=-1.0)
    pu = pt.uniq_idx.to_numpy()
    Xp = emb[np.searchsorted(idx, pu)]
    score = m.predict_proba(np.hstack([Fp, emb_scalars(Xp, Xtr, y, ref)]))[:, 1]

    oi = pd.read_parquet(a.orf_index)
    first = oi.drop_duplicates('uniq_idx')[['uniq_idx', 'group', 'seq_id', 'aa_len']]
    groups = oi.groupby('uniq_idx').group.apply(lambda s: ';'.join(sorted(set(s))))
    rank = (pd.DataFrame({'uniq_idx': pu, 'rf_score': score})
            .sort_values('rf_score', ascending=False).reset_index(drop=True))
    rank['rank_all'] = np.arange(1, len(rank) + 1)
    rank = rank.merge(first, on='uniq_idx', how='left')
    rank['groups'] = rank.uniq_idx.map(groups)

    os.makedirs(a.out_dir, exist_ok=True)
    rank.to_csv(os.path.join(a.out_dir, 'rf_pool_scored.csv'), index=False)
    print(f'scored {len(rank):,} — {rank.rf_score.min():.4f} to {rank.rf_score.max():.4f}')

    # ---- exclusions ----
    r1 = set(open(a.round1_aa).read().split())
    used = set(fa[fa.aa.isin(r1)].uniq_idx)
    print(f'identical to a Round 1 ORF: {int(rank.uniq_idx.isin(used).sum())}')

    cos_idx = set()
    if a.cosine_arm and os.path.exists(a.cosine_arm):
        cos_idx = set(pd.read_csv(a.cosine_arm).uniq_idx)
        ov = int(rank.head(a.n).uniq_idx.isin(cos_idx).sum())
        print(f'** overlap between the top {a.n} model picks and the similarity arm: '
              f'{ov} ({ov / a.n * 100:.1f}%) **')

    drop = used | (set() if a.allow_cosine_overlap else cos_idx)
    cand = rank[~rank.uniq_idx.isin(drop)].head(SCAN).reset_index(drop=True)

    # ---- primers ----
    cds = load_cds(a.cds_dir, cand.group.dropna().unique())
    rows = []
    for _, r in cand.iterrows():
        rec = cds.get(r.group, {}).get(r.seq_id)
        if rec is None:
            continue
        nt = 'ATG' + rec[0][3:]
        fw, tmf, warn_f = pick(nt)
        rv, tmr, warn_r = pick(rc(nt))
        rows.append(dict(uniq_idx=r.uniq_idx, rf_score=r.rf_score, rank_all=r.rank_all,
                         group=r.group, groups=r.groups, sequence_id=r.seq_id,
                         aa_len=r.aa_len, amplicon_bp=len(nt), start_codon=rec[0][:3],
                         primer_fw=fw, primer_rv=rv,
                         tm_fw=round(tmf, 2), tm_rv=round(tmr, 2),
                         warn_fw=warn_f, warn_rv=warn_r,
                         tm_pass=MIN_TM <= min(tmf, tmr) and max(tmf, tmr) < MAX_TM,
                         whole_sequence=nt))
    dd = pd.DataFrame(rows)
    print(f'designed {len(dd):,} / passed Tm {int(dd.tm_pass.sum()):,} '
          f'({dd.tm_pass.mean() * 100:.1f}%)')

    out = dd[dd.tm_pass].head(a.n).drop(columns='tm_pass').reset_index(drop=True)
    out.insert(0, 'rank_arm', np.arange(1, len(out) + 1))
    out['arm'] = 'rf'
    if len(out) < a.n:
        print(f'warning: only {len(out)} selected; raise SCAN')

    f = os.path.join(a.out_dir, 'round2_rf190.csv')
    out.to_csv(f, index=False)
    print(f'\n{len(out)} selected -> {f}')
    print(f'  rf_score {out.rf_score.min():.4f} to {out.rf_score.max():.4f}')
    print(f'  pool rank {int(out.rank_all.min()):,} to {int(out.rank_all.max()):,}')


if __name__ == '__main__':
    main()
