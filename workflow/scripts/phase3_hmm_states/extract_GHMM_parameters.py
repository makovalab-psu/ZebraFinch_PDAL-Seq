"""Summarise a fitted HMM and pre-bin its signal for the fit plots.

This does the work of three manuscript scripts in one pass, because they all
need the same two objects (the model and the data) in memory:

  * 07_PhyloHGMP_model/extract_GHMM_parameters.py  -- Gaussian means/variances
  * 07_PhyloHGMP_model/Mean_PDALseq_signal_in_states.R -- empirical mean per state
  * the histogram half of 07_PhyloHGMP_model/plot_to_test_Guassian.R

Pre-binning here is what keeps the plotting rule cheap: Plot_gaussian_fit.R
reads a few hundred rows instead of re-reading ~10^6 windows for each of the
twelve models.

Histogram densities are normalised so they can be compared directly against a
probability density:

  * per state,  count / (windows in that state * bin width) -> dnorm(mu_k, sd_k)
  * "global",   count / (all windows * bin width)           -> sum_k w_k dnorm_k

Usage:
    python extract_GHMM_parameters.py <GHMM.pkl> <GHMM_states.csv.gz> \
        <features> <out_parameters.csv> <out_histogram.csv>
"""

import os
import sys

import joblib
import numpy as np
import pandas as pd

N_BINS = 256

def ensure_dir(path):
    """snakemake makes output directories, a bare `python ...` call does not."""
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)



def main():
    if len(sys.argv) != 6:
        sys.exit(__doc__)

    pkl_path, states_path, features_arg, out_params, out_hist = sys.argv[1:6]

    features = features_arg.split(",")

    model = joblib.load(pkl_path)
    df = pd.read_csv(states_path)

    n_components = int(model.n_components)
    means = np.asarray(model.means_)
    # covariance_type="full" -> (n_components, n_features, n_features).
    variances = np.array([np.diag(np.atleast_2d(cov)) for cov in model.covars_])

    assigned = df["hidden_state"].to_numpy()
    n_total = assigned.shape[0]

    # ---------------------------------------------------------------- params
    rows = []
    for state in range(n_components):
        mask = assigned == state
        n_windows = int(mask.sum())
        for j, feature in enumerate(features):
            if n_windows > 0:
                observed = df.loc[mask, feature].to_numpy(dtype=float)
                empirical_mean = float(observed.mean())
                empirical_sd = float(observed.std(ddof=0))
            else:
                # A variational HMM is free to leave a state empty. That is
                # information, not an error, so the row is kept with NaNs.
                empirical_mean = np.nan
                empirical_sd = np.nan
            rows.append(
                {
                    "n_states": n_components,
                    "state": state,
                    "feature": feature,
                    "gaussian_mean": float(means[state, j]),
                    "gaussian_var": float(variances[state, j]),
                    "weight": n_windows / n_total,
                    "empirical_mean": empirical_mean,
                    "empirical_sd": empirical_sd,
                    "n_windows": n_windows,
                }
            )

    ensure_dir(out_params)
    pd.DataFrame(rows).to_csv(out_params, index=False)

    # ------------------------------------------------------------- histogram
    hist_rows = []
    for feature in features:
        signal = df[feature].to_numpy(dtype=float)
        edges = np.linspace(signal.min(), signal.max(), N_BINS + 1)
        width = float(edges[1] - edges[0])
        mids = (edges[:-1] + edges[1:]) / 2.0

        def add(label, values, denominator):
            counts, _ = np.histogram(values, bins=edges)
            keep = counts > 0
            for mid, count in zip(mids[keep], counts[keep]):
                hist_rows.append(
                    {
                        "n_states": n_components,
                        "feature": feature,
                        "state": label,
                        "bin_mid": float(mid),
                        "bin_width": width,
                        "count": int(count),
                        "density": float(count) / (denominator * width),
                    }
                )

        add("global", signal, float(n_total))

        for state in range(n_components):
            mask = assigned == state
            n_windows = int(mask.sum())
            if n_windows == 0:
                continue
            add(str(state), signal[mask], float(n_windows))

    ensure_dir(out_hist)
    pd.DataFrame(hist_rows).to_csv(out_hist, index=False)


if __name__ == "__main__":
    main()
