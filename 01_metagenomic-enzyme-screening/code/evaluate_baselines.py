"""Enrichment of every baseline, on the same metric as the trained ranker.

The baselines do not learn, so they need no cross-validation and are applied
directly to all rows. That favours them over the model, which is scored
out-of-fold. The comparison is deliberately set up that way.
"""
import argparse

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

p = argparse.ArgumentParser()
p.add_argument('--features', required=True)
p.add_argument('--rf-ef5', type=float, help='EF@5% of the trained model, for the last row')
p.add_argument('--rf-ef10', type=float)
a = p.parse_args()

d = pd.read_parquet(a.features)
y = d.y.to_numpy().astype(int)
n, P = len(y), y.sum()
pi = P / n
rng = np.random.default_rng(0)
KS = {f: max(1, int(round(f * n))) for f in (0.05, 0.10)}

print(f'n={n}  P={P}  pi={pi:.4f}   k: ' +
      ', '.join(f'{int(f*100)}%={k}' for f, k in KS.items()) + '\n')

HYD = [c for c in d.columns if c.startswith('pfam_top_') and any(t in c for t in (
    'Abhydrolase', 'Esterase', 'Lipase', 'BD-FAE', 'DLH', 'Hydrolase',
    'Peptidase_S15', 'AXE1', 'Palm_thioest', 'Chlorophyllase'))]

rules = {
    '(1) catalytic triad only': d.catalytic_triad_ok.to_numpy() > 0,
    '(2) triad + GxSxG': (d.catalytic_triad_ok.to_numpy() > 0) & (d.GxSxG.to_numpy() > 0),
    '(3) alpha/beta hydrolase fold': d.alpha_beta_hydrolase.to_numpy() > 0,
    '(4) Pfam hydrolase family': d[HYD].to_numpy().sum(1) > 0,
    '(5) any lipase motif': (d[['GxSxG', 'GDSL', 'HGG', 'HSL_pattern',
                                'SxxK']].to_numpy() > 0).any(1),
}

print(f'{"rule":<40}{"sel":>6}{"hit":>5}{"prec":>8}{"EF":>7}{"Fisher p":>11}')
print('-' * 77)
for name, m in rules.items():
    sel, hit = int(m.sum()), int(y[m].sum())
    prec = hit / sel if sel else np.nan
    _, pv = fisher_exact([[hit, sel - hit], [P - hit, (n - sel) - (P - hit)]],
                         alternative='greater')
    print(f'{name:<40}{sel:>6}{hit:>5}{prec:>8.3f}{prec / pi if sel else np.nan:>7.2f}{pv:>11.4f}')


def ef_at(score, ties_random=True, reps=1000):
    out = {}
    for f, k in KS.items():
        if ties_random:   # binary rules: break ties randomly and average
            pr = float(np.mean([y[np.lexsort((rng.random(n), -score))[:k]].mean()
                                for _ in range(reps)]))
        else:
            pr = y[np.argsort(-score)[:k]].mean()
        out[f] = (pr, pr / pi)
    return out


print(f'\n{"ranking":<40}{"P@5%":>8}{"EF@5%":>8}{"P@10%":>9}{"EF@10%":>9}')
print('-' * 77)
scores = {
    'cosine to lipase reference': d.cos_lipase.to_numpy(),
    'structure: TM-score to reference': np.nan_to_num(d.ref_tm.to_numpy(), nan=-1),
    'pLDDT (negative control)': np.nan_to_num(d.plddt_mean.to_numpy(), nan=-1),
}
for name, s in scores.items():
    r = ef_at(s, ties_random=False)
    print(f'{name:<40}{r[0.05][0]:>8.3f}{r[0.05][1]:>8.2f}{r[0.10][0]:>9.3f}{r[0.10][1]:>9.2f}')

for name, m in rules.items():
    r = ef_at(m.astype(float))
    print(f'{name + " [binary]":<40}{r[0.05][0]:>8.3f}{r[0.05][1]:>8.2f}'
          f'{r[0.10][0]:>9.3f}{r[0.10][1]:>9.2f}')

gate = (d.catalytic_triad_ok.to_numpy() > 0).astype(float) * 2 + d.cos_lipase.to_numpy()
r = ef_at(gate, ties_random=False)
print(f'{"triad gate then cosine":<40}{r[0.05][0]:>8.3f}{r[0.05][1]:>8.2f}'
      f'{r[0.10][0]:>9.3f}{r[0.10][1]:>9.2f}')

if a.rf_ef5:
    print(f'{"random forest (pooled OOF)":<40}{"":>8}{a.rf_ef5:>8.2f}{"":>9}{a.rf_ef10:>9.2f}')
