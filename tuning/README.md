# Parameter selection for MajorTrack on raw data

`majortrack_tune.py` selects, without ground truth, the three choices that most affect a MajorTrack analysis of
a temporal network: the size of the window over which snapshots are aggregated, the history parameter and the
absorption option. The criteria were derived and validated on the synthetic benchmark in
`../synthetic_simulations` (see its README).

## Input

- `.csv`: edge list with columns `time,source,target` and an optional `weight` column. With integer times,
  snapshots without edges are kept as empty snapshots so that windows remain aligned in time; with other
  sortable values, each distinct value is a snapshot. Node identifiers are arbitrary.
- `.npy`: array of shape `(T, N, N)` with one interaction matrix per snapshot.
- `.RDS`: R list of `T` interaction matrices (as produced by `synthetic_simulations/src/main.R`).

Interactions are treated as undirected. An element without interactions in a window is absent from that window
(label `-1` in the outputs). Window sizes that produce a window without any interaction are skipped, because
MajorTrack requires at least one group at every time point.

## Usage

Run from the repository root, in an environment with MajorTrack, `python-igraph`, `leidenalg`, `numpy`,
`pandas` and `matplotlib` (see `../synthetic_simulations/environment.yml`).

```bash
# 1. scan windows x histories x absorption (Leiden with 5 seeds per window)
python tuning/majortrack_tune.py scan --input data.csv --out results/
# 2. select the parameters (writes selection.json and the diagnostic figure selection.pdf)
python tuning/majortrack_tune.py select --out results/
# 3. run MajorTrack with the selected parameters (or --window W --history H [--no-absorption])
python tuning/majortrack_tune.py run --input data.csv --out results/
```

Options: `--windows` and `--histories` (candidate values; histories in units of windows), `--replicates`
(Leiden seeds, default 5) and `--tol` (tolerance for the history, default 0.01).

## Selection rules

- Window: the window size that maximises the total consistency `C_tot` (mean over dynamic clusters of the mean
  Jaccard similarity of their members between consecutive windows), averaged over histories `h <= 5` windows
  with absorption. The window that minimises the number of connections between dynamic clusters per transition
  (merges, splits and shifts) is reported as a cross-check; if the two differ by more than one candidate size,
  both should be inspected.
- History: the smallest `h` whose `C_tot`, at the selected window, lies within `tol` of the maximum. `C_tot`
  increases with `h` and saturates, so its saturation point, not its maximum, is used.
- Absorption: on (it removes transient dynamic clusters without changing the accuracy in the benchmark).

## Outputs

- `scan.csv`: diagnostics of every configuration (number of DCs, `C_tot`, connections per transition by type),
  averaged over seeds, with standard deviations; `scan_labels.npz`: DC labels of every configuration (first
  seed).
- `selection.json`, `selection.pdf`: selected parameters and diagnostic figure.
- `dc_membership_w-<w>_h-<h>.csv`: DC label of every element in every window; `dc_heatmap_w-<w>_h-<h>.pdf`.

## Validation

`../synthetic_simulations/src/validate_tuning.py` applies the procedure to the nine synthetic scenarios as if
they were raw data. The selected configurations were better than the naive configuration (`w = 1`, `h = 1`) in
all scenarios, by 0.05–0.29 in global NMI, and within 0.09 (median 0.03) of the best configuration with windows
not longer than the planted regimes (`../synthetic_simulations/data/tuning_validation.{csv,pdf}`).
