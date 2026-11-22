# Synthetic benchmark for MajorTrack

This folder contains a generator of temporal interaction networks with a planted dynamic community structure
(R), and the Python scripts used to evaluate MajorTrack on them: the effect of the history parameter, of the
aggregation of snapshots over time windows, and of absorption, and the validation of the truth-free parameter
selection implemented in `../tuning/majortrack_tune.py`.

## Requirements

- R (>= 4.0) with the packages `gplots` and `rprojroot` (data generation).
- Python with MajorTrack and the packages listed in `environment.yml` (analysis):

```bash
conda env create -f synthetic_simulations/environment.yml
conda activate majortrack-benchmark
pip install -e .            # from the repository root
```

All commands below are run from the `synthetic_simulations/` folder.

## 1. Generation of the synthetic data (R)

`N` elements interact at every time step according to block probability matrices: two elements interact with
probability `p_in` if they belong to the same group, and with probability `p_out` otherwise. The planted
partition changes between consecutive regimes through merges, splits and drifts of elements, and each regime
lasts a fixed number of time steps.

Functions (`src/`):

- `generate_partitions(N, K, subset_sizes)`: random assignment of `N` elements to `K` groups of given sizes.
- `modify_partitions(assignments, operation)`: applies an operation to an assignment vector:
  - merge: `list(type = "merge", groups = c(a, b, ...))` relabels all listed groups with the identifier of the
    first one;
  - split: `list(type = "split", group = g, sizes = c(n1, n2, ...))` divides group `g` into consecutive blocks
    of the given sizes; every block receives a new identifier except the last one, which keeps `g`;
  - drift: `list(type = "drift", groups = c(a, b, ...), size = s)` moves `s` randomly chosen elements from each
    listed group to the next one, cyclically.
- `generate_probability_matrices(...)`: block probability matrix of every regime.
- `simulate_interactions(N, probability_matrices, time_steps)`: Bernoulli draws of the interaction matrices.
- `example_use.R`: minimal example.

Driver: `src/main.R`, with the operations of each problem defined in `src/problem_merge.R`,
`src/problem_merge-split.R` and `src/problem_merge-split-drift.R`. The problem and `p_in` can be edited at the
top of `main.R` or passed on the command line:

```bash
Rscript src/main.R problem_merge-split 0.5
```

Settings of the provided data: `N = 200`, `K = 10` initial groups of sizes 5, 20, 10, 40, 20, 7, 15, 55, 10
and 18, 10 regimes of 10 time steps (`T = 100`), `p_out = 0.05`, `p_in` in {0.2, 0.5, 0.8}, seed `25022025`.
This gives nine scenarios, each in `data/<problem>_N-200_K-10_pin-<p_in>_pbet-0.05/`, with:

- `matrices_interaction_<label>.RDS`: list of the `T` interaction matrices (input of the analysis);
- `assignments_<label>.RDS`: list of the 10 planted partitions (ground truth; regime `m` covers steps
  `10m ... 10m+9`);
- `Err_<label>.pdf`: deviation between the empirical and the planted interaction probabilities.

Heatmaps of the probability and interaction matrices are also produced by `main.R`; they are not tracked in
git.

## 2. Analysis (Python)

| Script | Purpose | Main outputs |
|---|---|---|
| `src/run_majortrack.py` | Leiden communities in every snapshot; MajorTrack for history `h` = 1–50; comparison with the ground truth (per-step and global NMI, number of DCs, spurious identity changes). Optional arguments select scenarios by name. | `data/<scenario>/majortrack_results/`, `data/majortrack_summary.{csv,pdf}` |
| `src/window_experiment.py` | Aggregation of the snapshots over windows of `w` = 1–25 steps; MajorTrack for `h` = 1–10 windows with and without absorption, on Leiden partitions and on the true partitions ("oracle"); accuracy, total consistency `C_tot`, connections between DCs per transition, and label-free event detection. | `data/<scenario>/window_results/`, `data/window_summary.csv` |
| `src/history_consistency.py` | `C_tot` and accuracy as a function of `h` without aggregation (requires `run_majortrack.py`). | `data/history_Ctot_w1.csv` |
| `src/plot_window_experiment.py` | Figures of the window, history and absorption analyses. | `data/window_effect.pdf`, `data/history_calibration.pdf`, `data/absorption_effect.pdf` |
| `src/validate_tuning.py` | Applies the truth-free parameter selection of `../tuning/majortrack_tune.py` to every scenario as if it were raw data, and evaluates the selected parameters with the ground truth. | `data/<scenario>/tuning/`, `data/tuning_validation.{csv,pdf}`, `data/tuning_heatmaps.pdf` |

Reproduction of all results (about 10 minutes in total):

```bash
python src/run_majortrack.py
python src/window_experiment.py
python src/history_consistency.py
python src/plot_window_experiment.py
python src/validate_tuning.py
```

Per-scenario outputs are not tracked in git; the summary files in `data/` are. The scripts accept substrings of
scenario names as arguments (e.g. `python src/run_majortrack.py merge_N-200_K-10_pin-0.5`); the summary files
are then rebuilt from the scenarios analysed so far, so the full sequence above is needed to reproduce the
tracked summaries.

### Evaluation measures

- Mean per-step NMI: normalised mutual information between the inferred and the planted partition at every
  step, averaged over the steps with more than one planted group.
- Global NMI: NMI over all (element, step) pairs, which also requires correct identities over time. Because
  the ground truth keeps the identifier of the first listed group after a merge, whereas MajorTrack keeps the
  identity of the contributor of the majority of the members, the windowed analyses compare with reference DCs
  obtained by running MajorTrack on the true partitions.
- `C_tot`: mean over DCs of the mean Jaccard similarity of their members between consecutive windows.
- Connections per transition: number of distinct pairs of DCs exchanging elements between consecutive
  windows (merges, splits and shifts).
- Event detection: events are connected components of the flows between DCs at a transition; an inferred event
  matches a planted one if it occurs within one window and the member sets involved can be paired one-to-one
  with Jaccard similarity >= 0.5. Precision, recall and F1 are reported.

### Main results

- With the true partitions as input, MajorTrack recovers every planted event; the remaining errors originate
  in the static community detection.
- Aggregating snapshots over windows of 4–10 steps increases the identity-aware accuracy substantially,
  particularly at low `p_in`; windows aligned with the planted regimes perform better than misaligned ones.
- Without ground truth, the window maximising `C_tot` (cross-checked with the window minimising the
  connections per transition) and the history at which `C_tot` saturates give configurations that are better
  than the naive configuration (`w = 1`, `h = 1`) in all scenarios, and within 0.09 NMI of the best
  configuration with `w <= 10`.
- Absorption removes transient DCs, mostly without aggregation, and does not change the accuracy.
