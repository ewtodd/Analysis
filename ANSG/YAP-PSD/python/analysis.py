import time
import numpy as np
from analysis_utils.io import load_tree_data
from psd_utils import regress_waveforms, ANALYSIS_CACHE_DIR
from regressors import get_default_regressors
import os

ROOT_FILES_DIR = "../macros/root_files/"
SCALAR_BRANCHES = [
    "pulse_height",
    "trigger_position",
    "long_integral",
    "light_output",
    "charge_comparison",
    "raw_shape_indicator",
    "clean_shape_indicator",
]


def main():
    os.makedirs("plots", exist_ok=True)

    # Check if analysis cache exists — skip expensive ROOT loading if so
    cache_exists = os.path.isdir(ANALYSIS_CACHE_DIR) and all(
        os.path.exists(os.path.join(ANALYSIS_CACHE_DIR, f)) for f in [
            "test_alpha_features.pkl", "test_gamma_features.pkl",
            "test_waveforms.npz", "regressor_names.pkl"
        ])

    if cache_exists:
        print("Analysis cache found — skipping ROOT file loading.")
        alpha_waveforms = gamma_waveforms = None
        alpha_features = gamma_features = None
    else:
        t0 = time.time()
        print("Loading alpha data (Am-241)...")
        alpha_features, alpha_waveforms = load_tree_data(
            ROOT_FILES_DIR + "Am241.root",
            scalar_branches=SCALAR_BRANCHES,
            array_branch="Samples",
        )
        t1 = time.time()
        print(
            f"Alpha events: {len(alpha_features)}, waveform shape: {alpha_waveforms.shape}"
        )
        print(f"Alpha loading took {t1 - t0:.2f}s")

        print("Loading gamma data (Na-22)...")
        gamma_features, gamma_waveforms = load_tree_data(
            ROOT_FILES_DIR + "Na22.root",
            scalar_branches=SCALAR_BRANCHES,
            array_branch="Samples",
        )
        t2 = time.time()
        print(
            f"Gamma events: {len(gamma_features)}, waveform shape: {gamma_waveforms.shape}"
        )
        print(f"Gamma loading took {t2 - t1:.2f}s")

    regressors = get_default_regressors()

    t3 = time.time()
    regress_waveforms(
        (alpha_waveforms, gamma_waveforms),
        (alpha_features, gamma_features),
        regressors=regressors,
    )
    t4 = time.time()
    print(f"Regression took {t4 - t3:.2f}s")


if __name__ == "__main__":
    main()
