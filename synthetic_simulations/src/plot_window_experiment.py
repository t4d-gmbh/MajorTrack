"""Figures for the window / history / absorption experiment.

Reads data/window_summary.csv (window_experiment.py) and
data/history_Ctot_w1.csv, and writes data/window_effect.pdf,
data/history_calibration.pdf and data/absorption_effect.pdf.

Usage (from synthetic_simulations/):
  python src/plot_window_experiment.py
"""
import os

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use('Agg')
from matplotlib import pyplot as plt  # noqa: E402

DIR_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
PROBLEMS = ['merge', 'merge-split', 'merge-split-drift']
COLORS = {0.2: '#d95f02', 0.5: '#7570b3', 0.8: '#1b9e77'}
HIST_WINDOW_PLOT = 2      # history (in windows) used in the window figure


def load():
    d = pd.read_csv(os.path.join(DIR_DATA, 'window_summary.csv'))
    d['problem'] = d.scenario.str.extract(r'problem_([a-z-]+)_N')
    d['p_in'] = d.scenario.str.extract(r'pin-([0-9.]+)').astype(float)
    return d


def window_figure(d):
    x = d[d.absorption & (d.history == HIST_WINDOW_PLOT)]
    metrics = [('global_NMI_ref', 'global NMI (reference DCs)'),
               ('n_dcs', 'number of DCs'),
               ('conn_per_transition', 'connections / transition'),
               ('C_tot', r'$C_{\mathrm{tot}}$'),
               ('event_F1', 'event detection $F_1$')]
    fig, axes = plt.subplots(len(metrics), 3, figsize=(11, 12), sharex=True)
    for j, prob in enumerate(PROBLEMS):
        for i, (m, lab) in enumerate(metrics):
            ax = axes[i, j]
            for p in sorted(COLORS):
                s = x[(x.problem == prob) & (x.p_in == p) & (x.input == 'leiden')]
                ax.plot(s.window, s[m], 'o-', color=COLORS[p], ms=4,
                        label=rf'Leiden, $p_{{in}}$={p}')
            o = x[(x.problem == prob) & (x.p_in == 0.5) & (x.input == 'oracle')]
            ax.plot(o.window, o[m], 'k--', lw=1, label='oracle partitions')
            if m == 'n_dcs':
                ax.plot(o.window, o.true_n_groups, ':', color='grey', lw=1,
                        label='true groups')
            ax.axvline(10, color='grey', lw=0.6, ls=':')
            if j == 0:
                ax.set_ylabel(lab)
            if i == 0:
                ax.set_title(prob)
            if i == len(metrics) - 1:
                ax.set_xlabel('window size $w$ (time steps)')
            if m in ('n_dcs', 'conn_per_transition') and (x[x.problem == prob][m] > 0).any():
                ax.set_yscale('log')
    axes[1, 0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(DIR_DATA, 'window_effect.pdf'))
    plt.close(fig)


def history_figure():
    h = pd.read_csv(os.path.join(DIR_DATA, 'history_Ctot_w1.csv'))
    h['problem'] = h.sc.str.split(' ').str[0]
    h['p_in'] = h.sc.str.split(' ').str[1].astype(float)
    fig, axes = plt.subplots(2, 3, figsize=(11, 5.5), sharex=True)
    for j, prob in enumerate(PROBLEMS):
        for p in sorted(COLORS):
            s = h[(h.problem == prob) & (h.p_in == p)]
            axes[0, j].plot(s.h, s.C_tot, 'o-', color=COLORS[p], ms=4, label=rf'$p_{{in}}$={p}')
            axes[1, j].plot(s.h, s.refNMI, 'o-', color=COLORS[p], ms=4)
        axes[0, j].set_title(prob)
        axes[1, j].set_xlabel('history $h$ (time steps, $w=1$)')
        for i in range(2):
            axes[i, j].set_xscale('log')
            axes[i, j].axvline(10, color='grey', lw=0.6, ls=':')
    axes[0, 0].set_ylabel(r'$C_{\mathrm{tot}}$')
    axes[1, 0].set_ylabel('global NMI (reference DCs)')
    axes[0, 0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(DIR_DATA, 'history_calibration.pdf'))
    plt.close(fig)


def absorption_figure(d):
    x = d[d.input == 'leiden']
    piv = x.pivot_table(index=['scenario', 'window', 'history'], columns='absorption',
                        values=['n_dcs', 'C_tot', 'conn_per_transition', 'global_NMI_ref'])
    metrics = [('n_dcs', r'$\Delta$ number of DCs'), ('C_tot', r'$\Delta C_{\mathrm{tot}}$'),
               ('conn_per_transition', r'$\Delta$ connections / transition'),
               ('global_NMI_ref', r'$\Delta$ global NMI (ref.)')]
    fig, axes = plt.subplots(1, 4, figsize=(13, 3))
    for ax, (m, lab) in zip(axes, metrics):
        diff = (piv[(m, True)] - piv[(m, False)]).groupby(level='window')
        mean, lo, hi = diff.mean(), diff.quantile(0.1), diff.quantile(0.9)
        ax.fill_between(mean.index, lo, hi, color='#7570b3', alpha=0.25, lw=0)
        ax.plot(mean.index, mean, 'o-', color='#7570b3', ms=4)
        ax.axhline(0, color='k', lw=0.6)
        ax.set_xlabel('window size $w$')
        ax.set_ylabel(lab)
    fig.tight_layout()
    fig.savefig(os.path.join(DIR_DATA, 'absorption_effect.pdf'))
    plt.close(fig)


if __name__ == '__main__':
    d = load()
    window_figure(d)
    history_figure()
    absorption_figure(d)
