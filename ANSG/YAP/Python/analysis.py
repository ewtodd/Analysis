import pandas as pd
import h5py
import numpy as np
from analyze_psd import (
    regress_waveforms,
    compute_and_analyze_psd,
    process_waveforms,
)
import pickle


def main():

    sources = [0, 2, 3]  # Use integers instead of strings
    merged_file = "all_converted_waveforms.h5"

    print("Loading features...")
    features = pd.read_hdf(merged_file, key="features")

    # Check what source_id values actually exist
    print("Unique source_id values:", features["source_id"].unique())
    print("Source_id data type:", features["source_id"].dtype)

    alpha_features = features[features["source_id"] == 0].reset_index(
        drop=True
    )  # Integer 0
    gamma_beta_features = features[features["source_id"] == 2].reset_index(
        drop=True
    )  # Integer 2

    print(f"Alpha features: {len(alpha_features)} events")
    print(f"Gamma/beta features: {len(gamma_beta_features)} events")

    # Load waveforms from the merged file and filter by indices
    print("Loading waveforms...")
    with h5py.File(merged_file, "r") as f:
        all_waveforms = f["raw_waveform"][:]

    # Get indices for each source
    alpha_indices = features[features["source_id"] == 0].index.values
    gamma_beta_indices = features[features["source_id"] == 2].index.values

    # Extract waveforms for each source
    alpha_waveforms_array = all_waveforms[alpha_indices]
    gamma_beta_waveforms_array = all_waveforms[gamma_beta_indices]

    # Convert to DataFrames
    alpha_waveforms = pd.DataFrame(
        alpha_waveforms_array,
        columns=[f"sample_{i}" for i in range(alpha_waveforms_array.shape[1])],
    )
    gamma_beta_waveforms = pd.DataFrame(
        gamma_beta_waveforms_array,
        columns=[f"sample_{i}" for i in range(gamma_beta_waveforms_array.shape[1])],
    )

    print(f"Alpha waveforms: {alpha_waveforms.shape}")
    print(f"Gamma/beta waveforms: {gamma_beta_waveforms.shape}")
    waveforms = (alpha_waveforms, gamma_beta_waveforms)
    features_tuple = (alpha_features, gamma_beta_features)

    (
        (test_alpha_waveforms, test_gamma_beta_waveforms),
        (test_alpha_features, test_gamma_beta_features),
    ) = regress_waveforms(
        waveforms,
        features_tuple,
        process_func=process_waveforms,
        random_state=42,
        model_file="regressor.pkl",
        output_dir="psd_analysis",
    )

    compute_and_analyze_psd(
        (test_alpha_waveforms, test_gamma_beta_waveforms),
        (test_alpha_features, test_gamma_beta_features),
        output_dir="psd_analysis",
    )


if __name__ == "__main__":
    main()
