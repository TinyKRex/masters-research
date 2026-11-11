"""Similarity-only arm: 190 Round 2 candidates by cosine distance alone.

Reference point   mean embedding of 98 characterised bacterial lipases plus the
                  17 confirmed Round 1 hits
Pool              ORFs passing the 250-410 aa length filter
Excluded          sequences identical to any Round 1 ORF — already synthesised

Primer Tm is constrained to 55-63 C. The upper bound is measured, not assumed:
in Round 1, primers at Tm >= 63 C failed at 19.6% against 6.4% below it
(Fisher p = 0.0036). The lower bound is the annealing temperature. Round 1 never
tested primers below 54.1 C, so the observed trend is not extrapolated past it.
Length, GC and hairpin showed no effect in the same data and are not constrained.

    python select_cosine.py --cds-dir CDS/ --emb 'prott5_*.h5' \
        --ref-emb ref_prott5.h5 --orf-index orf_index.parquet \
        --pool-fasta pool.fasta --dataset dataset.csv --round1-aa round1.txt
"""
import argparse
import os

import numpy as np
import pandas as pd

from design_primers import load_cds, pick, rc
from emb_io import load_emb, read_fasta_map

MIN_TM, MAX_TM = 55.0, 63.0
N_ARM = 190
SCAN = 1200          # candidates to design primers for; about half pass Tm


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--cds-dir', required=True)
    p.add_argument('--emb', required=True)
    p.add_argument('--ref-emb', required=True)
    p.add_argument('--orf-index', required=True)
    p.add_argument('--pool-fasta', required=True)
    p.add_argument('--dataset', required=True)
    p.add_argument('--round1-aa', required=True, help='one Round 1 sequence per line')
    p.add_argument('--out-dir', default='.')
    p.add_argument('--n', type=int, default=N_ARM)
    a = p.parse_args()

    idx, emb = load_emb(a.emb)
    _, ref = load_emb(a.ref_emb)
    fa = read_fasta_map(a.pool_fasta)
    ds = pd.read_csv(a.dataset)

    pos = fa.merge(ds[ds.y == 1][['aa']], on='aa')
    hit_emb = emb[np.searchsorted(idx, pos.uniq_idx.values)]
    print(f'reference = {len(ref)} known lipases + {len(hit_emb)} confirmed Round 1 hits')

    c = np.vstack([ref, hit_emb]).mean(0)
    cos = (emb @ c) / (np.linalg.norm(emb, axis=1) * np.linalg.norm(c))

    oi = pd.read_parquet(a.orf_index)
    first = oi.drop_duplicates('uniq_idx')[['uniq_idx', 'group', 'seq_id', 'aa_len']]
    groups = oi.groupby('uniq_idx').group.apply(lambda s: ';'.join(sorted(set(s))))

    rank = pd.DataFrame({'uniq_idx': idx, 'cosine': cos}).sort_values(
        'cosine', ascending=False)
    rank['rank_all'] = np.arange(1, len(rank) + 1)
    rank = rank.merge(first, on='uniq_idx', how='left')
    rank['groups'] = rank.uniq_idx.map(groups)

    os.makedirs(a.out_dir, exist_ok=True)
    np.save(os.path.join(a.out_dir, 'centroid_prott5.npy'), c)
    rank.to_csv(os.path.join(a.out_dir, 'cosine_ranked.csv.gz'), index=False)

    r1 = set(open(a.round1_aa).read().split())
    used = set(fa[fa.aa.isin(r1)].uniq_idx)
    print(f'excluded {len(used)} sequences identical to Round 1 ORFs')
    cand = rank[~rank.uniq_idx.isin(used)].head(SCAN).reset_index(drop=True)

    cds = load_cds(a.cds_dir, cand.group.dropna().unique())
    rows = []
    for _, r in cand.iterrows():
        rec = cds.get(r.group, {}).get(r.seq_id)
        if rec is None:
            continue
        nt = 'ATG' + rec[0][3:]                   # uniform start codon
        fw, tmf, warn_f = pick(nt)
        rv, tmr, warn_r = pick(rc(nt))
        rows.append(dict(uniq_idx=r.uniq_idx, cosine=r.cosine, rank_all=r.rank_all,
                         group=r.group, groups=r.groups, sequence_id=r.seq_id,
                         aa_len=r.aa_len, amplicon_bp=len(nt), start_codon=rec[0][:3],
                         primer_fw=fw, primer_rv=rv,
                         tm_fw=round(tmf, 2), tm_rv=round(tmr, 2),
                         warn_fw=warn_f, warn_rv=warn_r,
                         tm_pass=MIN_TM <= min(tmf, tmr) and max(tmf, tmr) < MAX_TM,
                         whole_sequence=nt))
    d = pd.DataFrame(rows)
    print(f'designed {len(d):,} / passed Tm {MIN_TM:g}-{MAX_TM:g} C '
          f'{int(d.tm_pass.sum()):,} ({d.tm_pass.mean() * 100:.1f}%)')

    out = d[d.tm_pass].head(a.n).drop(columns='tm_pass').reset_index(drop=True)
    out.insert(0, 'rank_arm', np.arange(1, len(out) + 1))
    out['arm'] = 'cosine'
    if len(out) < a.n:
        print(f'warning: only {len(out)} selected; raise SCAN')

    f = os.path.join(a.out_dir, 'round2_cosine190.csv')
    out.to_csv(f, index=False)
    print(f'\n{len(out)} selected -> {f}')
    print(f'  cosine {out.cosine.min():.4f} to {out.cosine.max():.4f}')
    print(f'  overall rank {int(out.rank_all.min()):,} to {int(out.rank_all.max()):,}')


if __name__ == '__main__':
    main()
