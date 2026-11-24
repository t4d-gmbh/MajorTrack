"""Validation of the truth-free parameter selection (tuning/majortrack_tune.py)
on the synthetic scenarios.

For each scenario, the interaction snapshots are treated as raw data: the
scan / select steps of majortrack_tune.py choose the window, the history and
the absorption option without using the ground truth. The ground truth is then
used to evaluate every scanned configuration (global NMI with the reference
DCs, i.e. MajorTrack on the true partitions) and to locate the selected
configuration in this ranking, together with a naive configuration (w=1,
h=1) and the best configuration in hindsight.

Outputs:
  data/<scenario>/tuning/            scan.csv, scan_labels.npz, selection.json,
                                     selection.pdf, validation.csv
  data/tuning_validation.csv         one row per scenario
  data/tuning_validation.pdf         ranking of the selected configurations
  data/tuning_heatmaps.pdf           reference / naive / selected / best (w<=10)
                                     heatmaps

Usage (from synthetic_simulations/):
  python src/validate_tuning.py [substring ...]
"""
import glob
import json
import os
import sys

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use('Agg')
from matplotlib import pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', '..', 'tuning'))
import majortrack_tune as mtt  # noqa: E402
from run_majortrack import DIR_DATA, read_scenario  # noqa: E402
from window_experiment import accuracy, event_detection, event_signatures  # noqa: E402

SHOWCASE = ['problem_merge_N-200_K-10_pin-0.5_pbet-0.05',
            'problem_merge-split_N-200_K-10_pin-0.5_pbet-0.05',
            'problem_merge-split-drift_N-200_K-10_pin-0.2_pbet-0.05']


def short(label):
    return label.split('_N')[0][len('problem_'):] + r', $p_{in}$=' + label.split('pin-')[1][:3]


def validate(folder):
    label, mats, truth = read_scenario(folder)
    out = os.path.join(folder, 'tuning')
    mats_sym = [((m + m.T) > 0).astype(float) for m in mats]
    df, labels = mtt.scan(mats_sym, out=out, verbose=False)
    sel = mtt.select(df)
    with open(os.path.join(out, 'selection.json'), 'w') as f:
        json.dump(sel, f, indent=2)
    mtt.selection_figure(df, sel, os.path.join(out, 'selection.pdf'))
    # evaluation with the ground truth
    valid = np.array([len(np.unique(x)) > 1 for x in truth])
    ref = mtt.run_majortrack(truth, 1, True)
    true_sig = event_signatures(truth, list(range(len(truth))))
    acc = []
    for key, dc in labels.items():
        w, h, a = (int(x[1:]) for x in key.split('_'))
        step, glob_, glob_ref = accuracy(dc, truth, w, valid, ref)
        f1 = event_detection(dc, w, true_sig)[2]
        acc.append(dict(window=w, history=h, absorption=bool(a), key=key,
                        mean_step_NMI=step, global_NMI_ref=glob_ref, event_F1=f1))
    acc = pd.DataFrame(acc)
    val = df.merge(acc, on=['window', 'history', 'absorption'])
    val.to_csv(os.path.join(out, 'validation.csv'), index=False)
    key_sel = f'w{sel["window"]}_h{sel["history"]}_a1'
    best = acc.loc[acc.global_NMI_ref.idxmax()]
    a_i = acc.set_index('key')
    s_sel = a_i.loc[key_sel, 'global_NMI_ref']
    s_naive = a_i.loc['w1_h1_a1', 'global_NMI_ref']
    # best configuration among windows not longer than the planted regimes
    # (longer windows merge consecutive events and lose temporal resolution)
    short_w = acc[acc.window <= 10]
    best10 = short_w.loc[short_w.global_NMI_ref.idxmax()]
    row = dict(scenario=label, sel_window=sel['window'], sel_history=sel['history'],
               window_min_connections=sel['window_min_connections'],
               NMI_selected=s_sel, NMI_naive=s_naive, NMI_best=best.global_NMI_ref,
               best_key=best.key, NMI_best_w10=best10.global_NMI_ref, best_w10_key=best10.key,
               F1_selected=a_i.loc[key_sel, 'event_F1'], F1_naive=a_i.loc['w1_h1_a1', 'event_F1'],
               F1_best=best.event_F1, F1_max=acc.event_F1.max(),
               stepNMI_selected=a_i.loc[key_sel, 'mean_step_NMI'],
               stepNMI_naive=a_i.loc['w1_h1_a1', 'mean_step_NMI'],
               stepNMI_best=best.mean_step_NMI, n_configs=len(acc),
               rank_selected=int((acc.global_NMI_ref > s_sel).sum()) + 1,
               percentile_selected=100 * (acc.global_NMI_ref <= s_sel).mean(),
               percentile_selected_w10=100 * (short_w.global_NMI_ref <= s_sel).mean())
    return row, labels, truth, ref, key_sel, best10.key


def heatmap_figure(cases, path):
    fig, axes = plt.subplots(len(cases), 4, figsize=(14, 3.2 * len(cases)), squeeze=False)
    for i, (label, labels, truth, ref, key_sel, key_best) in enumerate(cases):
        order = np.lexsort(truth[::-1])
        mtt.membership_heatmap(axes[i, 0], ref, order, 'reference DCs')
        for j, (key, name) in enumerate([('w1_h1_a1', 'naive'), (key_sel, 'selected'),
                                         (key_best, 'best in hindsight ($w\\leq10$)')], start=1):
            w, h, _ = (int(x[1:]) for x in key.split('_'))
            mtt.membership_heatmap(axes[i, j], labels[key], order,
                                   f'{name}: $w$={w}, $h$={h}', w)
            axes[i, j].set_xlim(-0.5, truth.shape[0] - 0.5)
        axes[i, 0].set_ylabel(short(label) + '\nelement')
        for ax in axes[i]:
            for b in range(10, truth.shape[0], 10):
                ax.axvline(b - 0.5, color='k', lw=0.4, ls=':')
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def ranking_figure(folders, rows, path):
    fig, ax = plt.subplots(figsize=(11, 3.8))
    for i, (f, r) in enumerate(zip(folders, rows)):
        v = pd.read_csv(os.path.join(f, 'tuning', 'validation.csv'))
        jitter = np.random.default_rng(0).uniform(-0.18, 0.18, len(v))
        ax.scatter(i + jitter, v.global_NMI_ref, s=6, color='lightgrey', zorder=1)
        ax.scatter(i, r['NMI_naive'], marker='v', color='#d95f02', s=40, zorder=2,
                   label='naive ($w$=1, $h$=1)' if i == 0 else None)
        ax.scatter(i, r['NMI_selected'], marker='o', color='#1b9e77', s=45, zorder=3,
                   label='selected (truth-free)' if i == 0 else None)
        ax.scatter(i, r['NMI_best'], marker='*', color='k', s=55, zorder=2,
                   label='best in hindsight' if i == 0 else None)
        ax.scatter(i, r['NMI_best_w10'], marker='*', facecolor='none', edgecolor='k', s=55,
                   zorder=2, label=r'best in hindsight, $w\leq10$' if i == 0 else None)
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels([short(r['scenario']) for r in rows], rotation=30, ha='right', fontsize=8)
    ax.set_ylabel('global NMI (reference DCs)')
    ax.legend(fontsize=8, loc='lower right')
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


if __name__ == '__main__':
    folders = sorted(f for f in glob.glob(os.path.join(DIR_DATA, 'problem_*')) if os.path.isdir(f))
    selected = sys.argv[1:]
    folders = [f for f in folders if not selected or any(s in os.path.basename(f) for s in selected)]
    rows, cases = [], {}
    for f in folders:
        print('validating', os.path.basename(f), flush=True)
        row, labels, truth, ref, key_sel, key_best = validate(f)
        rows.append(row)
        cases[row['scenario']] = (row['scenario'], labels, truth, ref, key_sel, key_best)
    pd.DataFrame(rows).to_csv(os.path.join(DIR_DATA, 'tuning_validation.csv'), index=False)
    ranking_figure(folders, rows, os.path.join(DIR_DATA, 'tuning_validation.pdf'))
    show = [cases[s] for s in SHOWCASE if s in cases]
    if show:
        heatmap_figure(show, os.path.join(DIR_DATA, 'tuning_heatmaps.pdf'))
    print(pd.DataFrame(rows).round(3).to_string())
