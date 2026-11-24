"""Total consistency C_tot and accuracy as a function of the history parameter,
for the DCs obtained without aggregation (w=1) by run_majortrack.py
(data/<scenario>/majortrack_results/dc_membership_history-<h>.csv).

Output: data/history_Ctot_w1.csv

Usage (from synthetic_simulations/):
  python src/history_consistency.py
"""
import glob
import os
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import normalized_mutual_info_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_majortrack import DIR_DATA, HISTORIES, read_scenario  # noqa: E402
from window_experiment import consistency, run_majortrack  # noqa: E402

if __name__ == '__main__':
    rows = []
    for f in sorted(glob.glob(os.path.join(DIR_DATA, 'problem_*'))):
        if not os.path.isdir(os.path.join(f, 'majortrack_results')):
            continue  # scenario not analysed with run_majortrack.py
        label, _, truth = read_scenario(f)
        valid = np.array([len(np.unique(x)) > 1 for x in truth])
        ref = run_majortrack(truth, truth.shape[1], 1, True)
        sc = label.split('_N')[0][len('problem_'):] + ' ' + label.split('pin-')[1][:3]
        for h in HISTORIES:
            dc = pd.read_csv(os.path.join(f, 'majortrack_results',
                                          f'dc_membership_history-{h}.csv'),
                             index_col=0).values.T
            rows.append(dict(sc=sc, h=h, C_tot=consistency(dc)[0],
                             refNMI=normalized_mutual_info_score(ref[valid].ravel(),
                                                                 dc[valid].ravel()),
                             n_dcs=len(np.unique(dc))))
    pd.DataFrame(rows).to_csv(os.path.join(DIR_DATA, 'history_Ctot_w1.csv'), index=False)
