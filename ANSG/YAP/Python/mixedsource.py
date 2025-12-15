import pandas as pd
import h5py
import numpy as np
import os
import pickle
import ROOT
from sklearn.metrics import roc_curve, auc
from analyze_psd import (
    process_waveforms,
)

# Import your plotting utilities
from PyPlottingUtils import PyPlottingUtils

plot_utils = PyPlottingUtils()
ROOT.gROOT.SetBatch(True)


def plot_classified_spectra_root(features_filtered, output_dir):
    """Plot energy spectra for alpha-like and gamma-like events using strict thresholds"""
    os.makedirs(output_dir, exist_ok=True)

    # Apply strict thresholds
    alpha_like_mask = features_filtered["Regressor_Output"] <= 0.1
    gamma_like_mask = features_filtered["Regressor_Output"] >= 0.9

    # Get energy data
    alpha_like_energies = features_filtered[alpha_like_mask][
        "light_output_keVee"
    ].values
    gamma_like_energies = features_filtered[gamma_like_mask][
        "light_output_keVee"
    ].values
    all_energies = features_filtered["light_output_keVee"].values

    # Create canvas
    canvas = ROOT.TCanvas("c_classified_spectra", "Classified Spectra", 1200, 900)
    plot_utils.ConfigureCanvas(canvas, logy=True)

    # Create histograms
    energy_min, energy_max = 0, 2000
    nbins = 100

    h_all = ROOT.TH1F("h_all", "", nbins, energy_min, energy_max)
    h_alpha = ROOT.TH1F("h_alpha", "", nbins, energy_min, energy_max)
    h_gamma = ROOT.TH1F("h_gamma", "", nbins, energy_min, energy_max)

    # Fill histograms
    for energy in all_energies:
        h_all.Fill(energy)
    for energy in alpha_like_energies:
        h_alpha.Fill(energy)
    for energy in gamma_like_energies:
        h_gamma.Fill(energy)

    # Configure histograms
    plot_utils.ConfigureHistogram(h_all, ROOT.kBlack)
    plot_utils.ConfigureHistogram(h_alpha, ROOT.kRed + 1)
    plot_utils.ConfigureHistogram(h_gamma, ROOT.kBlue + 1)

    h_all.SetLineStyle(2)
    h_all.GetXaxis().SetTitle("Light Output [keVee]")
    h_all.GetYaxis().SetTitle("Counts / 10 keV")
    h_all.SetTitle("")

    # Draw
    max_val = h_all.GetMaximum()
    h_all.SetMaximum(max_val * 1.2)
    h_all.Draw("HIST")
    h_alpha.Draw("HIST SAME")
    h_gamma.Draw("HIST SAME")

    # Legend
    leg = plot_utils.CreateLegend(0.58, 0.33, 0.92, 0.55)
    leg.AddEntry(h_all, f"All Events ", "l")
    leg.AddEntry(h_alpha, f"#alpha-like (#leq 0.02)", "f")
    leg.AddEntry(h_gamma, f"#gamma-ray-like (#geq 0.98)", "f")
    leg.Draw()

    canvas.SaveAs(os.path.join(output_dir, "classified_energy_spectra.pdf"))
    canvas.Close()

    print(f"Classified spectra saved to {output_dir}/classified_energy_spectra.pdf")
    print(f"Alpha-like events (≤0.1): {len(alpha_like_energies)}")
    print(f"Gamma-like events (≥0.9): {len(gamma_like_energies)}")
    print(
        f"Uncertain events (0.1 < x < 0.9): {len(all_energies) - len(alpha_like_energies) - len(gamma_like_energies)}"
    )


def analyze_and_plot_third_source(
    model_file,
    third_source_file,
    third_source_id=3,
    output_dir="third_source_analysis_full",
    forceRECALCULATE=False,
):
    """Load a regressor model, classify a third source file, and plot results using ROOT."""
    os.makedirs(output_dir, exist_ok=True)

    # Load the trained regressor model
    with open(model_file, "rb") as file:
        regressor = pickle.load(file)

    print(f"Loading third source data from {third_source_file}...")

    # Load features
    features = pd.read_hdf(third_source_file, key="features")

    # Filter by source_id
    third_mask = features["source_id"] == third_source_id
    third_features_filtered = features[third_mask].copy()

    print(
        f"Loaded {len(third_features_filtered)} events from third source (source_id={third_source_id})"
    )

    # Check if Regressor_Output already exists in features
    if (
        "Regressor_Output" in third_features_filtered.columns
        and forceRECALCULATE == False
    ):
        print("Regressor_Output already exists in features, using existing values...")
        regressor_output = third_features_filtered["Regressor_Output"].values
    else:
        print("Regressor_Output not found, processing waveforms...")

        # Load waveforms from the merged file and filter by indices
        print("Loading waveforms...")
        with h5py.File(third_source_file, "r") as f:
            all_waveforms = f["raw_waveform"][:]

        # Get indices for the third source
        third_indices = features[features["source_id"] == third_source_id].index.values

        # Extract waveforms for the third source
        third_waveforms_array = all_waveforms[third_indices]

        # Convert to DataFrame with proper column names
        third_waveforms = pd.DataFrame(
            third_waveforms_array,
            columns=[f"sample_{i}" for i in range(third_waveforms_array.shape[1])],
        )

        # Process waveforms for ML prediction
        print("Processing waveforms for ML prediction...")
        processed_waveforms = process_waveforms(third_waveforms)

        # Apply regressor to third source waveforms
        print("Applying ML classifier...")
        regressor_output = regressor.predict(processed_waveforms)

        # Save regressor output back to the HDF5 file
        print("Saving Regressor_Output to features...")

        # Update the full features dataframe with regressor output
        features.loc[third_mask, "Regressor_Output"] = regressor_output

        # Save updated features back to HDF5 file
        try:
            features.to_hdf(third_source_file, key="features", mode="a", complevel=9)
            print(f"Successfully saved Regressor_Output to {third_source_file}")
        except Exception as e:
            print(f"Warning: Could not save Regressor_Output to HDF5 file: {e}")
            print("Continuing with analysis using temporary values...")

    # Update filtered features with regressor output
    third_features_filtered["Regressor_Output"] = regressor_output

    # Simple classification using strict thresholds
    alpha_like_count = np.sum(regressor_output <= 0.02)
    gamma_like_count = np.sum(regressor_output >= 0.98)
    uncertain_count = len(regressor_output) - alpha_like_count - gamma_like_count

    print(f"\n=== Classification Results (Strict Thresholds) ===")
    print(f"Total events: {len(regressor_output)}")
    print(
        f"Alpha-like (≤0.1): {alpha_like_count} ({alpha_like_count/len(regressor_output)*100:.1f}%)"
    )
    print(
        f"Gamma-like (≥0.9): {gamma_like_count} ({gamma_like_count/len(regressor_output)*100:.1f}%)"
    )
    print(
        f"Uncertain (0.1 < x < 0.9): {uncertain_count} ({uncertain_count/len(regressor_output)*100:.1f}%)"
    )

    # Create plots using ROOT
    plot_classified_spectra_root(third_features_filtered, output_dir)

    return third_features_filtered


def main():
    """Main function to analyze the third source"""
    third_results = analyze_and_plot_third_source(
        "regressor.pkl",
        "all_converted_waveforms.h5",
        3,
    )

    print(f"\nAnalysis complete. Results saved to third_source_analysis_full/")


if __name__ == "__main__":
    main()
