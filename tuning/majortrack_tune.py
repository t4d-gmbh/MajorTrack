"""Truth-free selection of the window size, the history parameter and the
absorption option of MajorTrack for raw temporal interaction data.

The procedure follows the synthetic benchmark described in
DCD_manuscript/synthetic_experiment/synthetic_experiment.tex:

  scan    For every window size w (snapshots aggregated per window), detect
          static communities in each window (Leiden, modularity; repeated with
          several seeds, diagnostics averaged over the replicates) and run
          MajorTrack for every history h (in windows) and absorption on/off.
          Store the DC labels and truth-free diagnostics:
            - number of DCs,
            - total consistency C_tot (mean over DCs of the mean Jaccard
              similarity of their members between consecutive windows),
            - connections between DCs per transition (merges, splits, shifts).
  select  Choose w as the window that maximises C_tot (averaged over the
          histories h <= 5 windows) and report the window that minimises the
          connections per transition as a cross-check; then choose h as the
          smallest history whose C_tot is within `tol` of the maximum at the
          selected window. Absorption is kept on (it removes transient DCs
          without affecting accuracy in the benchmark).
  run     Run MajorTrack with the selected (or user-given) parameters and write
          the DC membership of every element and a membership heatmap.

Input formats (--input):
  *.csv   edge list with columns time,source,target[,weight]; times are
          integer snapshot indices (snapshots without edges are kept as empty
          snapshots) or other sortable values (one snapshot per distinct
          value); node identifiers are arbitrary.
  *.npy   array of shape (T, N, N) with one interaction matrix per snapshot.
  *.RDS   list of T interaction matrices (as produced by synthetic_simulations).

Elements without any interaction in a window are considered absent from that
window (label -1 in the outputs). Window sizes that produce a window without
any interaction are skipped, since MajorTrack requires groups at every time
point.

Usage:
  python majortrack_tune.py scan   --input data.csv --out results/
  python majortrack_tune.py select --out results/
  python majortrack_tune.py run    --input data.csv --out results/ [--window W --history H --no-absorption]
"""
import argparse
import json
import os
import warnings

import igraph as ig
import leidenalg
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use('Agg')
warnings.filterwarnings('ignore', category=SyntaxWarning)
from matplotlib import pyplot as plt  # noqa: E402
from majortrack import MajorTrack  # noqa: E402

DEFAULT_WINDOWS = [1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 25]
DEFAULT_HISTORIES = [1, 2, 3, 5, 10]
MIN_WINDOWS = 4          # a window size is scanned only if it yields >= 4 windows
SEED = 12345


# --------------------------------------------------------------------------- io
def load_snapshots(path):
    """Return (list of symmetric binary/weighted N x N arrays, node names)."""
    ext = os.path.splitext(path)[1].lower()
    if ext == '.csv':
        e = pd.read_csv(path)
        if 'weight' not in e:
            e['weight'] = 1.0
        nodes = np.array(sorted(set(e.source) | set(e.target), key=str))
        idx = {v: i for i, v in enumerate(nodes)}
        times = sorted(e.time.unique())
        if pd.api.types.is_integer_dtype(e.time):
            # integer snapshots: keep empty snapshots so that windows stay aligned
            times = list(range(int(min(times)), int(max(times)) + 1))
        mats = []
        for t in times:
            m = np.zeros((len(nodes), len(nodes)))
            s = e[e.time == t]
            for a, b, w in zip(s.source, s.target, s.weight):
                m[idx[a], idx[b]] += w
            mats.append(m)
    elif ext == '.npy':
        mats = list(np.load(path).astype(float))
        nodes = np.arange(mats[0].shape[0])
    elif ext == '.rds':
        import rdata
        mats = [np.asarray(m, dtype=float) for m in rdata.read_rds(path)]
        nodes = np.arange(mats[0].shape[0])
    else:
        raise ValueError(f'unknown input format: {ext}')
    mats = [((m + m.T) > 0).astype(float) for m in mats]
    for m in mats:
        np.fill_diagonal(m, 0)
    return mats, nodes


def aggregate(mats, w):
    """Non-overlapping windows of w snapshots; the weight of an edge is the
    number of snapshots of the window in which the pair interacted."""
    k = len(mats) // w
    return [sum(mats[j * w:(j + 1) * w]) for j in range(k)]


def empty_windows(wins):
    """Indices of windows without any interaction (MajorTrack cannot process
    a time point without groups)."""
    return [j for j, a in enumerate(wins) if not a.any()]


# ---------------------------------------------------------------- detection
def static_communities(adj, seed):
    """Leiden (modularity) on the elements present in the window; absent
    elements (no interaction) get label -1."""
    present = np.flatnonzero(adj.sum(axis=1) > 0)
    lab = -np.ones(adj.shape[0], dtype=int)
    if len(present) == 0:
        return lab
    sub = adj[np.ix_(present, present)]
    g = ig.Graph.Weighted_Adjacency(sub.tolist(), mode='upper', attr='weight', loops=False)
    part = leidenalg.find_partition(g, leidenalg.ModularityVertexPartition,
                                    weights='weight', seed=seed)
    lab[present] = part.membership
    return lab


def run_majortrack(groupings, history, absorption=True):
    """MajorTrack on a series of label vectors (-1 = absent). Returns the DC
    label of every element in every window (-1 = absent)."""
    n = groupings.shape[1]
    clusterings, individuals = [], []
    for lab in groupings:
        clusterings.append([set(np.flatnonzero(lab == c).tolist())
                            for c in np.unique(lab) if c >= 0])
        individuals.append(set(np.flatnonzero(lab >= 0).tolist()))
    mt = MajorTrack(clusterings=clusterings, individuals=individuals, history=history)
    mt.get_group_matchup('fraction')
    mt.get_dcs(absorption=absorption)
    mt.get_community_membership()
    mt.get_individual_membership()
    return np.array([[mt.individual_membership[t].get(i, -1) if i in individuals[t] else -1
                      for i in range(n)] for t in range(len(clusterings))])


# -------------------------------------------------------------- diagnostics
def consistency(labels):
    """C_tot: mean over DCs of the mean Jaccard similarity of their members
    between consecutive windows (elements present in both windows only)."""
    pairs = {}
    for t in range(1, len(labels)):
        res = (labels[t - 1] >= 0) & (labels[t] >= 0)
        for c in (set(labels[t][res]) & set(labels[t - 1][res])):
            a, b = (labels[t - 1] == c) & res, (labels[t] == c) & res
            pairs.setdefault(c, []).append((a & b).sum() / (a | b).sum())
    cmu = [np.mean(v) for v in pairs.values()]
    return np.mean(cmu) if cmu else np.nan


def connections(labels):
    """Mean number of distinct pairs of DCs exchanging elements per transition,
    split into merges, splits and shifts."""
    out = []
    for t in range(1, len(labels)):
        res = (labels[t - 1] >= 0) & (labels[t] >= 0)
        prev, nxt = set(labels[t - 1][labels[t - 1] >= 0]), set(labels[t][labels[t] >= 0])
        edges = {(a, b) for a, b in zip(labels[t - 1][res], labels[t][res]) if a != b}
        merges = sum(a not in nxt for a, b in edges)
        splits = sum(b not in prev and a in nxt for a, b in edges)
        out.append((len(edges), merges, splits, len(edges) - merges - splits))
    return np.array(out).mean(axis=0) if out else np.full(4, np.nan)


def diagnostics(dc):
    con = connections(dc)
    return dict(n_dcs=len(np.unique(dc[dc >= 0])), C_tot=consistency(dc),
                conn_per_transition=con[0], merges_per_transition=con[1],
                splits_per_transition=con[2], shifts_per_transition=con[3])


# --------------------------------------------------------------------- scan
def scan(mats, windows=DEFAULT_WINDOWS, histories=DEFAULT_HISTORIES, out=None,
         absorption=(True, False), replicates=5, verbose=True):
    """Diagnostics for every (window, history, absorption). Static detection is
    repeated with `replicates` Leiden seeds; diagnostics are averaged over the
    replicates (columns *_sd give their standard deviation) and the DC labels
    of the first replicate are stored."""
    rows, labels = [], {}
    for w in windows:
        wins = aggregate(mats, w)
        if len(wins) < MIN_WINDOWS:
            continue
        if empty_windows(wins):
            print(f'  skipped w={w}: windows without interactions '
                  f'{empty_windows(wins)}', flush=True)
            continue
        statics = [np.array([static_communities(a, SEED + 1000 * r + j)
                             for j, a in enumerate(wins)]) for r in range(replicates)]
        n_static = np.mean([len(np.unique(s[s >= 0])) for st in statics for s in st])
        for h in histories:
            if h > len(wins) - 1 and h > 1:
                continue
            for ab in absorption:
                reps = []
                for r, static in enumerate(statics):
                    dc = run_majortrack(static, h, ab)
                    if r == 0:
                        labels[f'w{w}_h{h}_a{int(ab)}'] = dc
                    reps.append(diagnostics(dc))
                reps = pd.DataFrame(reps)
                row = dict(window=w, n_windows=len(wins), history=h, history_steps=h * w,
                           absorption=ab, n_static_groups_mean=n_static)
                row.update(reps.mean().to_dict())
                row.update({f'{k}_sd': v for k, v in reps.std().to_dict().items()})
                rows.append(row)
        if verbose:
            print(f'  scanned w={w}', flush=True)
    df = pd.DataFrame(rows)
    if out:
        os.makedirs(out, exist_ok=True)
        df.to_csv(os.path.join(out, 'scan.csv'), index=False)
        np.savez_compressed(os.path.join(out, 'scan_labels.npz'), **labels)
    return df, labels


# ------------------------------------------------------------------- select
def select(df, tol=0.01, h_max_avg=5):
    """Truth-free choice of window, history and absorption."""
    on = df[df.absorption.astype(bool)]
    avg = on[on.history <= h_max_avg].groupby('window')[['C_tot', 'conn_per_transition']].mean()
    w_c = int(avg.C_tot.idxmax())
    w_conn = int(avg.conn_per_transition.idxmin())
    at_w = on[on.window == w_c].set_index('history').sort_index()
    h_sel = int(at_w.index[at_w.C_tot >= at_w.C_tot.max() - tol].min())
    return dict(window=w_c, history=h_sel, absorption=True,
                window_min_connections=w_conn,
                C_tot=float(at_w.loc[h_sel, 'C_tot']),
                n_dcs=int(at_w.loc[h_sel, 'n_dcs']),
                note=('C_tot maximum and connection minimum agree' if w_c == w_conn else
                      'C_tot maximum and connection minimum differ: inspect both windows'))


def selection_figure(df, sel, path):
    on = df[df.absorption.astype(bool)]
    fig, axes = plt.subplots(1, 4, figsize=(15, 3.4))
    for h, s in on.groupby('history'):
        axes[0].plot(s.window, s.C_tot, 'o-', ms=3, label=f'h={h}')
        axes[1].plot(s.window, s.conn_per_transition, 'o-', ms=3)
        axes[2].plot(s.window, s.n_dcs, 'o-', ms=3)
    for ax, lab in zip(axes[:3], [r'$C_{\mathrm{tot}}$', 'connections / transition',
                                  'number of DCs']):
        ax.axvline(sel['window'], color='k', lw=0.8, ls='--')
        ax.set_xlabel('window size $w$')
        ax.set_ylabel(lab)
    axes[1].axvline(sel['window_min_connections'], color='grey', lw=0.8, ls=':')
    axes[1].set_yscale('log')
    axes[2].set_yscale('log')
    axes[0].legend(fontsize=7, title='history (windows)', title_fontsize=7)
    s = on[on.window == sel['window']].sort_values('history')
    axes[3].plot(s.history, s.C_tot, 'o-', color='k', ms=4)
    axes[3].axvline(sel['history'], color='k', lw=0.8, ls='--')
    axes[3].set_xscale('log')
    axes[3].set_xlabel(f'history $h$ (windows, $w={sel["window"]}$)')
    axes[3].set_ylabel(r'$C_{\mathrm{tot}}$')
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# ---------------------------------------------------------------------- run
def membership_heatmap(ax, labels, order, title, w=1):
    """Element x time heatmap of labels (-1 = absent, white); each window is
    drawn with width w (in snapshots)."""
    lab = np.repeat(labels, w, axis=0)
    codes = {}
    for v in lab.T[order].ravel(order='F'):
        if v >= 0:
            codes.setdefault(v, len(codes))
    img = np.vectorize(lambda v: codes.get(v, -1))(lab[:, order].T).astype(float)
    img[img < 0] = np.nan
    ax.imshow(np.where(np.isnan(img), np.nan, img % 20), aspect='auto',
              interpolation='nearest', cmap=plt.get_cmap('tab20'), vmin=0, vmax=19)
    ax.set_title(f'{title} ({len(codes)} labels)', fontsize=9)
    ax.set_xlabel('time step')


def run(mats, window, history, absorption=True, out=None, nodes=None):
    wins = aggregate(mats, window)
    if empty_windows(wins):
        raise ValueError(f'windows without interactions: {empty_windows(wins)}; '
                         'use a larger window')
    static = np.array([static_communities(a, SEED + j) for j, a in enumerate(wins)])
    dc = run_majortrack(static, history, absorption)
    if out:
        os.makedirs(out, exist_ok=True)
        names = nodes if nodes is not None else np.arange(dc.shape[1])
        pd.DataFrame(dc.T, index=names,
                     columns=[f'window_{k}' for k in range(len(dc))]).to_csv(
            os.path.join(out, f'dc_membership_w-{window}_h-{history}.csv'), index_label='node')
        fig, ax = plt.subplots(figsize=(6, 4))
        order = np.lexsort(dc[::-1])
        membership_heatmap(ax, dc, order, f'DCs, w={window}, h={history}', window)
        ax.set_ylabel('element')
        fig.tight_layout()
        fig.savefig(os.path.join(out, f'dc_heatmap_w-{window}_h-{history}.pdf'))
        plt.close(fig)
    return dc


# ---------------------------------------------------------------------- cli
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('command', choices=['scan', 'select', 'run'])
    p.add_argument('--input', help='edge list (.csv), array (.npy) or list of matrices (.RDS)')
    p.add_argument('--out', required=True, help='output folder')
    p.add_argument('--windows', type=int, nargs='+', default=DEFAULT_WINDOWS)
    p.add_argument('--histories', type=int, nargs='+', default=DEFAULT_HISTORIES)
    p.add_argument('--replicates', type=int, default=5, help='Leiden seeds per window')
    p.add_argument('--tol', type=float, default=0.01, help='C_tot tolerance for h')
    p.add_argument('--window', type=int, help='window for run (default: selection.json)')
    p.add_argument('--history', type=int, help='history for run (default: selection.json)')
    p.add_argument('--no-absorption', action='store_true')
    a = p.parse_args()
    if a.command == 'scan':
        mats, _ = load_snapshots(a.input)
        scan(mats, a.windows, a.histories, a.out, replicates=a.replicates)
    elif a.command == 'select':
        df = pd.read_csv(os.path.join(a.out, 'scan.csv'))
        sel = select(df, a.tol)
        with open(os.path.join(a.out, 'selection.json'), 'w') as f:
            json.dump(sel, f, indent=2)
        selection_figure(df, sel, os.path.join(a.out, 'selection.pdf'))
        print(json.dumps(sel, indent=2))
    else:
        mats, nodes = load_snapshots(a.input)
        sel = {}
        if a.window is None or a.history is None:
            with open(os.path.join(a.out, 'selection.json')) as f:
                sel = json.load(f)
        run(mats, a.window or sel['window'], a.history or sel['history'],
            not a.no_absorption, a.out, nodes)


if __name__ == '__main__':
    main()
