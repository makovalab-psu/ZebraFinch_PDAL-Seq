"""Score every fitted HMM and compute BIC, for the scree plot.

Adapted from workflow/scripts/09_Human_HMM/calc_BIC.py in the PDAL-Seq
manuscript pipeline. Differences:

  * the state list is a command-line argument instead of a hardcoded list;
  * an occupied_states column is added. A variational HMM is free to leave
    states empty, so the number of states the model actually uses is itself
    evidence about the optimum: once occupancy stops tracking k, extra states
    are buying nothing.

The parameter count is kept EXACTLY as the manuscript computes it -- means +
full covariances + transitions, with the initial-state probabilities omitted
-- so the BIC values stay comparable with the human analysis. The params
column is written out so BIC can be recomputed under another convention
without refitting anything.

Usage:
    python calc_BIC.py <model_dir> <states> <input.csv> <features> <out.csv>
"""

import os
import sys

import joblib
import numpy as np
import pandas as pd

def ensure_dir(path):
    """snakemake makes output directories, a bare `python ...` call does not."""
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def main():
    if len(sys.argv) != 6:
        sys.exit(__doc__)

    model_dir, states_arg, in_csv, features_arg, out_csv = sys.argv[1:6]

    states = states_arg.split(",")
    features = features_arg.split(",")

    df = pd.read_csv(in_csv)
    data_array = df[features].to_numpy(dtype=float)
    n_observations = data_array.shape[0]
    log_n = float(np.log(n_observations))

    rows = []

    for state in states:
        print("scoring {} states".format(state), flush=True)

        model = joblib.load(os.path.join(model_dir, state, "GHMM.pkl"))

        log_likelihood = float(model.score(data_array))

        n_components = int(model.n_components)
        n_features = int(model.n_features)

        n_means = n_components * n_features
        n_covariances = n_components * (n_features * (n_features + 1)) // 2
        n_transitions = n_components * (n_components - 1)
        n_params = n_means + n_covariances + n_transitions

        bic = -2.0 * log_likelihood + n_params * log_n

        params = pd.read_csv(os.path.join(model_dir, state, "Parameters.csv"))
        occupied = int(
            params.loc[params["feature"] == features[0], "n_windows"].gt(0).sum()
        )

        rows.append(
            {
                "States": int(state),
                "BIC": bic,
                "log_l": log_likelihood,
                "params": n_params,
                "log_obs": log_n,
                "occupied_states": occupied,
                "n_observations": n_observations,
            }
        )

    out = pd.DataFrame(rows).sort_values("States")
    ensure_dir(out_csv)
    out.to_csv(out_csv, index=False)

    best = out.loc[out["BIC"].idxmin(), "States"]
    print("lowest BIC at {} states".format(best), flush=True)

    # Every model saw the same windows, so log likelihood cannot legitimately
    # fall as states are added. Where it does, that fit landed in a worse
    # optimum than a smaller model and its BIC describes the fit, not the
    # data. hmmlearn initialises from k-means, which splits wide clusters and
    # merges narrow ones, so this happens at particular k regardless of seed.
    running_max = out["log_l"].cummax()
    suspect = out.loc[out["log_l"] < running_max - 1e-9, "States"].tolist()
    if suspect:
        print(
            "WARNING: log likelihood dips at {} states -- those fits landed in "
            "a worse optimum than a smaller model, so read their BIC with "
            "suspicion and check plots/phase3/gaussian_fit/".format(
                ", ".join(str(s) for s in suspect)
            ),
            flush=True,
        )


if __name__ == "__main__":
    main()
