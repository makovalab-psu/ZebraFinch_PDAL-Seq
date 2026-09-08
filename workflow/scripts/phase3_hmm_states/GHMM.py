"""Fit a variational Gaussian HMM to the windowed PDAL-Seq signal.

Adapted from workflow/scripts/07_PhyloHGMP_model/GHMM.py in the PDAL-Seq
manuscript pipeline. Differences:

  * the input carries chr/start/end alongside the feature columns, so the
    feature names are passed explicitly and the coordinates are re-attached
    to the predictions instead of being re-joined by row position later;
  * random_state is fixed, so a re-run reproduces the same segmentation.

Usage:
    python GHMM.py <input.csv> <features> <n_states> <out_states.csv.gz> <out.pkl>
"""

import os
import sys

import joblib
import numpy as np
import pandas as pd
from hmmlearn import vhmm

COORD_COLUMNS = ["chr", "start", "end"]
RANDOM_STATE = 42

def ensure_dir(path):
    """snakemake makes output directories, a bare `python ...` call does not."""
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)



def main():
    if len(sys.argv) != 6:
        sys.exit(__doc__)

    in_csv, features_arg, n_states_arg, out_states, out_pkl = sys.argv[1:6]

    features = features_arg.split(",")
    n_states = int(n_states_arg)

    df = pd.read_csv(in_csv)

    missing = [c for c in COORD_COLUMNS + features if c not in df.columns]
    if missing:
        sys.exit("columns missing from {}: {}".format(in_csv, ", ".join(missing)))

    data_array = df[features].to_numpy(dtype=float)

    if not np.isfinite(data_array).all():
        sys.exit(
            "{} contains non-finite values; the HMM cannot be fit to them".format(in_csv)
        )

    print(
        "fitting {} states to {} windows x {} features".format(
            n_states, data_array.shape[0], data_array.shape[1]
        ),
        flush=True,
    )

    model = vhmm.VariationalGaussianHMM(
        n_components=n_states,
        covariance_type="full",
        n_iter=2000,
        random_state=RANDOM_STATE,
        verbose=True,
    )

    model.fit(data_array)

    hidden_states = model.predict(data_array)

    occupied = len(np.unique(hidden_states))
    print(
        "converged={} monitor_iter={} occupied_states={}/{}".format(
            model.monitor_.converged, model.monitor_.iter, occupied, n_states
        ),
        flush=True,
    )

    out = df[COORD_COLUMNS + features].copy()
    out["hidden_state"] = hidden_states
    ensure_dir(out_states)
    out.to_csv(out_states, index=False)

    ensure_dir(out_pkl)
    joblib.dump(model, out_pkl)


if __name__ == "__main__":
    main()
