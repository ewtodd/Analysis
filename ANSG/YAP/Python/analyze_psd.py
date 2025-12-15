import os
import ROOT
import h5py
import numpy as np
import pandas as pd
import pickle
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, roc_curve, auc
from joblib import parallel_backend, Parallel, delayed
from concurrent.futures import ThreadPoolExecutor
from scipy.stats import norm
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter

# Remove matplotlib imports and replace with your PyPlottingUtils
from PyPlottingUtils import PyPlottingUtils  # Import your class

plot_utils = PyPlottingUtils()
ROOT.gROOT.SetBatch(True)

N_JOBS = 32


def normalize_waveform(waveform):
    """Normalize a waveform by dividing by its maximum value (unless max is zero)."""
    max_val = np.max(waveform)
    if max_val != 0:
        return waveform / max_val
    else:
        return waveform


def process_waveforms(waveform_df, n_jobs=N_JOBS):
    """
    Normalizes each waveform by its maximum value.
    Input:
      waveform_df: DataFrame with shape (n_waveforms, n_samples)
    Output:
      normalized_waveforms: numpy array with shape (n_waveforms, n_samples)
    """
    # Convert DataFrame to numpy array for processing
    waveform_array = waveform_df.values

    with ThreadPoolExecutor(max_workers=n_jobs) as executor:
        normalized = np.array(list(executor.map(normalize_waveform, waveform_array)))
    return normalized


def regress_waveforms(
    waveforms,
    features,
    process_func,
    random_state=42,
    model_file="regressor.pkl",
    output_dir="psd_analysis",
):
    os.makedirs(output_dir, exist_ok=True)

    energy_lower_dict = {
        "0": 500,
        "2": 0,
    }
    energy_upper_dict = {
        "0": 1750,
        "2": 1750,
    }

    alpha_waveforms, gamma_beta_waveforms = waveforms
    alpha_features, gamma_beta_features = features

    # Filter by energy for training ONLY
    energymask_alpha = (
        alpha_features["light_output_keVee"] <= energy_upper_dict["0"]
    ) & (alpha_features["light_output_keVee"] >= energy_lower_dict["0"])

    energymask_gamma_beta = (
        gamma_beta_features["light_output_keVee"] <= energy_upper_dict["2"]
    ) & (gamma_beta_features["light_output_keVee"] >= energy_lower_dict["2"])

    print(alpha_features[energymask_alpha])
    # Apply energy masks to DataFrames using .loc
    alpha_masked_waveforms = alpha_waveforms.loc[energymask_alpha].reset_index(
        drop=True
    )
    gamma_beta_masked_waveforms = gamma_beta_waveforms.loc[
        energymask_gamma_beta
    ].reset_index(drop=True)

    frac = 0.30
    n_alpha_train = int(len(alpha_masked_waveforms) * frac)
    n_gamma_beta_train = int(len(gamma_beta_masked_waveforms) * frac)
    min_samples = min(n_alpha_train, n_gamma_beta_train)
    max_samples = 10000
    if min_samples > max_samples:
        print(f"Min samples exceeds maximum desired number. Using {max_samples}...")
        min_samples = max_samples

    print(f"{frac*100}% of alpha masked: {n_alpha_train}")
    print(f"{frac*100}% of gamma/beta masked: {n_gamma_beta_train}")
    print(f"Samples for balanced training: {min_samples}")

    # Sample training data directly from DataFrames
    np.random.seed(random_state)
    train_alpha_waveforms = alpha_masked_waveforms.sample(
        n=min_samples, random_state=random_state
    )
    train_gamma_beta_waveforms = gamma_beta_masked_waveforms.sample(
        n=min_samples, random_state=random_state
    )

    # Get the original indices of training samples for later exclusion
    alpha_train_original_indices = (
        alpha_waveforms.loc[energymask_alpha].iloc[train_alpha_waveforms.index].index
    )
    gamma_beta_train_original_indices = (
        gamma_beta_waveforms.loc[energymask_gamma_beta]
        .iloc[train_gamma_beta_waveforms.index]
        .index
    )

    # Process training waveforms
    with Parallel(n_jobs=2) as parallel:
        train_results = parallel(
            delayed(process_func)(group)
            for group in [train_alpha_waveforms, train_gamma_beta_waveforms]
        )

    x_train = np.vstack(train_results)
    y_train = np.array([0] * len(train_results[0]) + [1] * len(train_results[1]))

    # Train or load model
    if os.path.exists(model_file):
        print("Loading existing model.")
        with open(model_file, "rb") as file:
            regressor = pickle.load(file)
    else:
        print("Training new model...")
        regressor = RandomForestRegressor(
            n_estimators=250,
            max_depth=30,
            random_state=random_state,
            max_samples=0.632,
            max_features="sqrt",
            n_jobs=-1,
            verbose=1,
        )
        regressor.fit(x_train, y_train)
        with open(model_file, "wb") as file:
            pickle.dump(regressor, file)

    # Get training predictions for histogram
    y_train_pred = regressor.predict(x_train)

    # Create test data by dropping training samples
    test_alpha_waveforms = alpha_waveforms.drop(
        alpha_train_original_indices
    ).reset_index(drop=True)
    test_gamma_beta_waveforms = gamma_beta_waveforms.drop(
        gamma_beta_train_original_indices
    ).reset_index(drop=True)

    test_alpha_features = alpha_features.drop(alpha_train_original_indices).reset_index(
        drop=True
    )
    test_gamma_beta_features = gamma_beta_features.drop(
        gamma_beta_train_original_indices
    ).reset_index(drop=True)

    # Process test waveforms for prediction
    with Parallel(n_jobs=2) as parallel:
        test_results = parallel(
            delayed(process_func)(group)
            for group in [test_alpha_waveforms, test_gamma_beta_waveforms]
        )

    X_test = np.vstack(test_results)
    y_test_pred = regressor.predict(X_test)
    # Plot test scores using ROOT - REPLACE MATPLOTLIB SECTION
    alpha_test_pred = y_test_pred[: len(test_alpha_waveforms)]
    gamma_beta_test_pred = y_test_pred[len(test_alpha_waveforms) :]

    plot_score_histogram_root(
        alpha_test_pred,
        gamma_beta_test_pred,
        "Test Set Scores",
        os.path.join(output_dir, "test_score_histogram.pdf"),
    )
    # Filter test data to 1000-1200 keVee range
    alpha_1000_1200_mask = (test_alpha_features["light_output_keVee"] >= 1000) & (
        test_alpha_features["light_output_keVee"] <= 1200
    )
    gamma_1000_1200_mask = (test_gamma_beta_features["light_output_keVee"] >= 1000) & (
        test_gamma_beta_features["light_output_keVee"] <= 1200
    )

    # Get the filtered features for 1000-1200 keVee
    test_alpha_features_1000_1200 = test_alpha_features[
        alpha_1000_1200_mask
    ].reset_index(drop=True)
    test_gamma_beta_features_1000_1200 = test_gamma_beta_features[
        gamma_1000_1200_mask
    ].reset_index(drop=True)

    # Get predictions for 1000-1200 keVee range
    alpha_1000_1200_pred = y_test_pred[: len(test_alpha_waveforms)][
        alpha_1000_1200_mask
    ]
    gamma_1000_1200_pred = y_test_pred[len(test_alpha_waveforms) :][
        gamma_1000_1200_mask
    ]

    if len(alpha_1000_1200_pred) > 0 and len(gamma_1000_1200_pred) > 0:
        # Add regressor output to the filtered features
        test_alpha_features_1000_1200["Regressor_Output"] = alpha_1000_1200_pred
        test_gamma_beta_features_1000_1200["Regressor_Output"] = gamma_1000_1200_pred

        # Calculate AUC for 1000-1200 keVee range (RF only for printing)
        y_1000_1200_pred = np.concatenate([alpha_1000_1200_pred, gamma_1000_1200_pred])
        y_1000_1200_true = np.array(
            [0] * len(alpha_1000_1200_pred) + [1] * len(gamma_1000_1200_pred)
        )

        from sklearn.metrics import roc_curve, auc

        fpr_1000_1200, tpr_1000_1200, _ = roc_curve(y_1000_1200_true, y_1000_1200_pred)
        auc_1000_1200 = auc(fpr_1000_1200, tpr_1000_1200)

        print(f"\n=== 1000-1200 keVee Test Results ===")
        print(f"Alpha events in range: {len(alpha_1000_1200_pred)}")
        print(f"Gamma/Beta events in range: {len(gamma_1000_1200_pred)}")
        print(f"RF AUC (1000-1200 keVee): {auc_1000_1200:.4f}")

        # Now do the comprehensive analysis with all three methods
        analyze_1000_1200_methods(
            test_alpha_features_1000_1200,
            test_gamma_beta_features_1000_1200,
            output_dir,
        )

        # Individual score histogram for RF in this range
        plot_score_histogram_root(
            alpha_1000_1200_pred,
            gamma_1000_1200_pred,
            "Random Forest Scores - 1000-1200 keVee Range",
            os.path.join(output_dir, "rf_score_histogram_1000_1200_keVee.pdf"),
        )

    else:
        print(f"\n=== 1000-1200 keVee Test Results ===")
        print(
            f"Alpha events in range: {len(alpha_1000_1200_pred) if 'alpha_1000_1200_pred' in locals() else 0}"
        )
        print(
            f"Gamma/Beta events in range: {len(gamma_1000_1200_pred) if 'gamma_1000_1200_pred' in locals() else 0}"
        )
        print("Warning: Insufficient data in 1000-1200 keVee range for AUC calculation")

    # Add regressor output to features
    test_alpha_features["Regressor_Output"] = y_test_pred[: len(test_alpha_waveforms)]
    test_gamma_beta_features["Regressor_Output"] = y_test_pred[
        len(test_alpha_waveforms) :
    ]
    # Add regressor output to features
    test_alpha_features["Regressor_Output"] = y_test_pred[: len(test_alpha_waveforms)]
    test_gamma_beta_features["Regressor_Output"] = y_test_pred[
        len(test_alpha_waveforms) :
    ]
    plot_feature_importance_waveform_with_average_root(
        regressor,  # Your trained model
        (test_alpha_waveforms, test_gamma_beta_waveforms),
        (test_alpha_features, test_gamma_beta_features),
        ["0", "2"],
        output_dir=output_dir,
    )
    return (
        (test_alpha_waveforms, test_gamma_beta_waveforms),
        (test_alpha_features, test_gamma_beta_features),
    )


def compute_and_analyze_psd(test_waveforms, test_features, output_dir):
    """Simplified PSD analysis using pre-computed values for ROC comparison only"""

    test_alpha_waveforms, test_gamma_beta_waveforms = test_waveforms
    test_alpha_features, test_gamma_beta_features = test_features

    # Use pre-computed PSD values - just rename for consistency
    test_alpha_features = test_alpha_features.copy()
    test_gamma_beta_features = test_gamma_beta_features.copy()
    test_alpha_features["charge_comparison_psd"] = test_alpha_features[
        "charge_comparison_psd"
    ]
    test_gamma_beta_features["charge_comparison_psd"] = test_gamma_beta_features[
        "charge_comparison_psd"
    ]
    test_alpha_features["si_psd"] = test_alpha_features["si_psd"]
    test_gamma_beta_features["si_psd"] = test_gamma_beta_features["si_psd"]

    # Only do ROC analysis comparing all three methods
    analyze_all_methods(
        test_alpha_features, test_gamma_beta_features, output_dir=output_dir
    )


def verify_waveform_feature_alignment(waveforms_df, features_df, n_check=5):
    """Verify that waveforms and features are still aligned"""
    print("Checking waveform-feature alignment...")

    for i in range(min(n_check, len(waveforms_df))):
        # Get waveform from DataFrame
        waveform_row = waveforms_df.iloc[i]
        waveform = waveform_row.values

        # Get corresponding feature row
        feature_row = features_df.iloc[i]

        # Check if they make sense together
        wf_peak = np.max(waveform)
        stored_peak = feature_row.get("pulse_height", "N/A")

        print(f"Row {i}:")
        print(f"  Waveform peak: {wf_peak:.3f}")
        print(f"  Stored peak: {stored_peak}")
        print(f"  Match: {'✓' if abs(wf_peak - stored_peak) < 0.001 else '✗'}")


def analyze_all_methods(test_alpha_features, test_gamma_beta_features, output_dir):
    """ROC analysis comparing ML, Charge Comparison, and si_psd PSD"""

    # Filter alphas (500-1750 keVee)
    alpha_mask = (test_alpha_features["light_output_keVee"] >= 500) & (
        test_alpha_features["light_output_keVee"] <= 1750
    )
    test_alpha_features_filtered = test_alpha_features[alpha_mask].reset_index(
        drop=True
    )

    # Create combined dataset
    test_features = pd.concat(
        [test_alpha_features_filtered, test_gamma_beta_features]
    ).reset_index(drop=True)

    y_true = np.array(
        [0] * len(test_alpha_features_filtered) + [1] * len(test_gamma_beta_features)
    )

    # All three methods for comparison
    all_methods = ["Regressor_Output", "charge_comparison_psd", "si_psd"]
    all_method_names = ["Random Forest", "Charge Comparison", "Shape Indicator"]

    # Create ROC curve plot comparing all three
    plot_unified_roc_curves_root(
        test_features,
        y_true,
        all_methods,
        all_method_names,
        output_dir=output_dir,
    )


def plot_feature_importance_waveform_with_average_root(
    regressor,
    test_waveforms,
    test_features,
    sources,
    output_dir="psd_analysis",
):
    """Plot the feature importance and average waveform (both normalized) using ROOT."""
    os.makedirs(output_dir, exist_ok=True)

    # Unpack test data
    test_alpha_waveforms, test_gamma_beta_waveforms = test_waveforms
    test_alpha_features, test_gamma_beta_features = test_features

    # Filter alpha events for training energy range (500-1750 keVee) to match training
    alpha_mask = (test_alpha_features["light_output_keVee"] >= 500) & (
        test_alpha_features["light_output_keVee"] <= 1750
    )
    alpha_waveforms_filtered = test_alpha_waveforms[alpha_mask]

    # Process the filtered waveforms to match ML input format
    alpha_waveforms_processed = process_waveforms(alpha_waveforms_filtered)

    # Calculate average normalized waveform
    avg_waveform = np.mean(alpha_waveforms_processed, axis=0)

    # Extract feature importances
    importances = regressor.feature_importances_

    # Ensure importances match waveform length
    if len(importances) != len(avg_waveform):
        print(
            f"Warning: Feature importance length ({len(importances)}) != waveform length ({len(avg_waveform)})"
        )
        # Pad or truncate as needed
        if len(importances) < len(avg_waveform):
            zero_padded_importances = np.zeros_like(avg_waveform)
            zero_padded_importances[: len(importances)] = importances
        else:
            zero_padded_importances = importances[: len(avg_waveform)]
    else:
        zero_padded_importances = importances

    # Normalize both to [0, 1] for comparison
    avg_waveform_norm = (
        avg_waveform / np.max(avg_waveform)
        if np.max(avg_waveform) > 0
        else avg_waveform
    )
    importances_norm = (
        zero_padded_importances / np.max(zero_padded_importances)
        if np.max(zero_padded_importances) > 0
        else zero_padded_importances
    )

    # Time axis (2 ns sampling)
    x_values = np.arange(len(avg_waveform)) * 2

    # Create ROOT canvas
    canvas = ROOT.TCanvas(
        "c_waveform_importance", "Waveform and Feature Importance", 1600, 1000
    )
    plot_utils.ConfigureCanvas(canvas)

    # Create TGraphs
    graph_waveform = ROOT.TGraph(
        len(x_values), x_values.astype(np.float64), avg_waveform_norm.astype(np.float64)
    )
    graph_importance = ROOT.TGraph(
        len(x_values), x_values.astype(np.float64), importances_norm.astype(np.float64)
    )

    # Configure waveform graph
    graph_waveform.SetLineColor(ROOT.kBlue + 1)
    graph_waveform.SetLineWidth(3)
    graph_waveform.SetTitle("")
    graph_waveform.GetXaxis().SetTitle("Time [ns]")
    graph_waveform.GetYaxis().SetTitle("Normalized Amplitude [a.u.]")
    graph_waveform.GetXaxis().SetRangeUser(0, x_values[-1])
    graph_waveform.GetYaxis().SetRangeUser(0, 1.1)

    # Configure importance graph
    graph_importance.SetLineColor(ROOT.kRed + 1)
    graph_importance.SetLineWidth(3)
    graph_importance.SetLineStyle(2)  # Dashed line

    # Draw graphs
    graph_waveform.Draw("AL")
    graph_importance.Draw("L SAME")

    # Add vertical line at trigger (17 * 2 = 34 ns)
    trigger_line = ROOT.TLine(34, 0, 34, 1.1)
    trigger_line.SetLineColor(ROOT.kGreen + 2)
    trigger_line.SetLineStyle(9)
    trigger_line.SetLineWidth(2)
    trigger_line.Draw()

    # Create legend
    leg = plot_utils.CreateLegend(0.55, 0.6, 0.92, 0.85)
    leg.AddEntry(graph_waveform, "Average #alpha Waveform", "l")
    leg.AddEntry(graph_importance, "Feature Importance", "l")
    leg.AddEntry(trigger_line, "Trigger (34 ns)", "l")
    leg.Draw()

    # Save the plot
    output_path = os.path.join(
        output_dir, "feature_importance_waveform_with_average.pdf"
    )
    canvas.SaveAs(output_path)
    canvas.Close()

    print(f"Feature importance plot saved to {output_path}")


def analyze_auc_and_feature_importance_by_sample_size(
    alpha_waveforms,
    gamma_beta_waveforms,
    alpha_features,
    gamma_beta_features,
    sample_sizes,
    process_func,
    output_dir="psd_analysis",
    energy_cut=(500, 1750),
):
    """
    Trains models with different training sample sizes, calculates AUC, and plots feature importance.
    Applies energy filtering to alpha events like the original code.

    Args:
        alpha_waveforms: DataFrame of alpha waveforms
        gamma_beta_waveforms: DataFrame of gamma/beta waveforms
        alpha_features: DataFrame with alpha features
        gamma_beta_features: DataFrame with gamma/beta features
        sample_sizes: list of integers with training sample sizes
        process_func: preprocessing function for waveforms
        output_dir: directory to save plots
        energy_cut: tuple (min_keV, max_keV) for alpha energy filtering

    Returns:
        auc_list: list of AUC values for each sample size
        feature_importances: list of arrays of feature importances
    """
    importances_all = []
    auc_list = []

    os.makedirs(output_dir, exist_ok=True)

    # Apply energy cut for alpha events (like original code does for training)
    alpha_mask = (alpha_features["light_output_keVee"] >= energy_cut[0]) & (
        alpha_features["light_output_keVee"] <= energy_cut[1]
    )
    alpha_waveforms_filtered = alpha_waveforms.loc[alpha_mask].reset_index(drop=True)
    alpha_features_filtered = alpha_features.loc[alpha_mask].reset_index(drop=True)

    print(
        f"After energy filter ({energy_cut[0]}-{energy_cut[1]} keVee): {len(alpha_waveforms_filtered)} alpha events"
    )
    print(f"Total gamma/beta events: {len(gamma_beta_waveforms)}")

    for n_samples in sample_sizes:
        print(f"Training with {n_samples} samples per class (after energy cut)...")

        # Sample the waveforms and features from filtered alpha data
        alpha_sample = alpha_waveforms_filtered.sample(
            n=min(n_samples, len(alpha_waveforms_filtered)), random_state=42
        ).reset_index(drop=True)
        gamma_sample = gamma_beta_waveforms.sample(
            n=min(n_samples, len(gamma_beta_waveforms)), random_state=42
        ).reset_index(drop=True)

        # Prepare training data
        x_train_alpha = process_func(alpha_sample)
        x_train_gamma = process_func(gamma_sample)
        X_train = np.vstack((x_train_alpha, x_train_gamma))
        y_train = np.array([0] * len(x_train_alpha) + [1] * len(x_train_gamma))

        # Train model
        from sklearn.ensemble import RandomForestRegressor

        model = RandomForestRegressor(
            n_estimators=250,
            max_depth=30,
            max_features="sqrt",
            random_state=42,
            n_jobs=-1,
        )
        model.fit(X_train, y_train)

        # Save feature importances
        importances_all.append(model.feature_importances_)

        # Prepare test data: Use remaining filtered data excluding train samples
        alpha_test = alpha_waveforms_filtered.drop(alpha_sample.index).reset_index(
            drop=True
        )
        gamma_test = gamma_beta_waveforms.drop(gamma_sample.index).reset_index(
            drop=True
        )

        x_test_alpha = process_func(alpha_test)
        x_test_gamma = process_func(gamma_test)
        X_test = np.vstack((x_test_alpha, x_test_gamma))
        y_test_true = np.array([0] * len(x_test_alpha) + [1] * len(x_test_gamma))

        y_test_pred = model.predict(X_test)

        from sklearn.metrics import roc_curve, auc

        fpr, tpr, _ = roc_curve(y_test_true, y_test_pred)
        auc_score = auc(fpr, tpr)
        auc_list.append(auc_score)
        print(f"Sample size: {n_samples}, AUC: {auc_score:.4f}")

    # Plot AUC vs Sample size
    plt.figure(figsize=(30, 20))
    plt.plot(sample_sizes, auc_list, marker="o", linestyle="-", color="blue")
    plt.xlabel("Number of Training Samples per Class")
    plt.ylabel("ROC AUC")
    plt.title(
        f"ROC AUC vs Number of Training Samples\n(Alpha energy: {energy_cut[0]}-{energy_cut[1]} keVee)"
    )
    plt.grid(True)
    plt.xscale("log")
    plt.tight_layout()
    auc_plot_path = os.path.join(output_dir, "auc_vs_samples_filtered.pdf")
    plt.savefig(auc_plot_path, dpi=300)
    plt.close()

    # Plot feature importance comparison
    plt.figure(figsize=(30, 20))
    colors = plt.cm.viridis(np.linspace(0, 1, len(sample_sizes)))
    time_axis = np.arange(len(importances_all[0])) * 2  # Assuming 2 ns sampling

    for i, (imp, n_samples) in enumerate(zip(importances_all, sample_sizes)):
        plt.plot(
            time_axis, imp, label=f"{n_samples} samples", color=colors[i], linewidth=3
        )

    plt.xlabel("Time [ns]")
    plt.ylabel("Feature Importance")
    plt.title(
        f"Feature Importance vs Time for Different Training Sample Sizes\n(Alpha energy: {energy_cut[0]}-{energy_cut[1]} keVee)"
    )
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    fi_plot_path = os.path.join(
        output_dir, "feature_importance_vs_samples_filtered.pdf"
    )
    plt.savefig(fi_plot_path, dpi=300)
    plt.close()

    print(f"Saved AUC plot to {auc_plot_path}")
    plt.close()

    return auc_list, importances_all


def analyze_1000_1200_methods(
    test_alpha_features_1000_1200, test_gamma_beta_features_1000_1200, output_dir
):
    """ROC analysis comparing ML, charge_comparison_psd, and si_psd PSD for 1000-1200 keVee range using ROOT"""

    # Create combined dataset for 1000-1200 keVee
    test_features_1000_1200 = pd.concat(
        [test_alpha_features_1000_1200, test_gamma_beta_features_1000_1200]
    ).reset_index(drop=True)

    y_true_1000_1200 = np.array(
        [0] * len(test_alpha_features_1000_1200)
        + [1] * len(test_gamma_beta_features_1000_1200)
    )

    # All three methods for comparison
    all_methods = ["Regressor_Output", "charge_comparison_psd", "si_psd"]
    all_method_names = ["Random Forest", "Charge Comparison", "Shape Indicator"]

    # Create unified ROC plot for 1000-1200 keVee
    plot_unified_roc_curves_root_1000_1200(
        test_features_1000_1200,
        y_true_1000_1200,
        all_methods,
        all_method_names,
        output_dir=output_dir,
    )

    # Create method comparison histograms for 1000-1200 keVee
    plot_method_histograms_root_1000_1200(
        test_features_1000_1200,
        y_true_1000_1200,
        all_methods,
        all_method_names,
        output_dir=output_dir,
    )


def plot_unified_roc_curves_root_1000_1200(
    test_features, y_true, methods, method_names, output_dir="psd_analysis"
):
    """Plot ROC curves for all methods in 1000-1200 keVee range using ROOT"""
    os.makedirs(output_dir, exist_ok=True)
    colors = [ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2]
    target_fpr = 0.05

    canvas = ROOT.TCanvas("c_unified_roc_1000_1200", "Unified ROC 1000-1200", 1200, 900)
    plot_utils.ConfigureCanvas(canvas)

    roc_graphs = []
    leg = plot_utils.CreateLegend(0.25, 0.2, 0.9, 0.4)

    print(f"\n=== 1000-1200 keVee ROC Analysis ===")

    for i, (method, name) in enumerate(zip(methods, method_names)):
        scores = test_features[method].values
        fpr, tpr, thresholds = roc_curve(y_true, scores)

        # Find threshold and TPR at 5% FPR
        index = np.argmin(np.abs(fpr - target_fpr))
        threshold_at_5pct_fpr = thresholds[index]
        tpr_at_5pct_fpr = tpr[index]
        actual_fpr = fpr[index]

        print(f"{name} (1000-1200 keVee):")
        print(f"  Threshold at {actual_fpr:.3f} FPR: {threshold_at_5pct_fpr:.3f}")
        print(f"  TPR at {actual_fpr:.3f} FPR: {tpr_at_5pct_fpr:.3f}")

        auc_score = auc(fpr, tpr)
        print(f"  AUC: {auc_score:.4f}")
        print()

        # Create TGraph
        roc_graph = ROOT.TGraph(len(fpr), fpr, tpr)
        roc_graph.SetLineColor(colors[i])
        roc_graph.SetLineWidth(3)
        roc_graphs.append(roc_graph)

        # Draw
        if i == 0:
            roc_graph.SetTitle("")
            roc_graph.GetXaxis().SetTitle("False Positive Rate (1 - Specificity)")
            roc_graph.GetYaxis().SetTitle("True Positive Rate (Sensitivity)")
            roc_graph.GetXaxis().SetRangeUser(0, 1)
            roc_graph.GetYaxis().SetRangeUser(0, 1)
            roc_graph.Draw("AL")
        else:
            roc_graph.Draw("L SAME")

        leg.AddEntry(roc_graph, f"{name} (AUC = {auc_score:.3f})", "l")

    # Add diagonal line
    diagonal = ROOT.TLine(0, 0, 1, 1)
    diagonal.SetLineColor(ROOT.kBlack)
    diagonal.SetLineStyle(2)
    diagonal.Draw()

    leg.Draw()

    canvas.SaveAs(
        os.path.join(output_dir, "unified_roc_comparison_1000_1200_keVee.pdf")
    )
    canvas.Close()

    print(
        f"Unified ROC curves (1000-1200 keVee) saved to {output_dir}/unified_roc_comparison_1000_1200_keVee.pdf"
    )


def plot_method_histograms_root_1000_1200(
    test_features, y_true, methods, method_names, output_dir="psd_analysis"
):
    """Plot histograms for all methods in 1000-1200 keVee range using ROOT"""
    os.makedirs(output_dir, exist_ok=True)

    # Create canvas with 3 subpads for all three methods
    canvas = ROOT.TCanvas(
        "c_methods_1000_1200", "Method Histograms 1000-1200 keVee", 2700, 900
    )
    canvas.Divide(3, 1)

    for i, (method, name) in enumerate(zip(methods, method_names)):
        canvas.cd(i + 1)
        plot_utils.ConfigureCanvas(ROOT.gPad, logy=True)

        alpha_scores = test_features[y_true == 0][method].values
        gamma_scores = test_features[y_true == 1][method].values

        # Determine range for all scores
        all_scores = np.concatenate([alpha_scores, gamma_scores])
        score_min = np.min(all_scores)
        score_max = np.max(all_scores)

        # Create histograms
        h_alpha = ROOT.TH1F(f"h_alpha_1000_1200_{i}", "", 50, score_min, score_max)
        h_gamma = ROOT.TH1F(f"h_gamma_1000_1200_{i}", "", 50, score_min, score_max)

        # Fill histograms
        for val in alpha_scores:
            h_alpha.Fill(val)
        for val in gamma_scores:
            h_gamma.Fill(val)

        # Configure
        plot_utils.ConfigureHistogram(h_alpha, ROOT.kRed + 1)
        plot_utils.ConfigureHistogram(h_gamma, ROOT.kBlue + 1)

        h_alpha.GetXaxis().SetTitle(f"{name} Score")
        h_alpha.GetYaxis().SetTitle("Counts")
        h_alpha.SetTitle("")

        # Draw
        max_val = max(h_alpha.GetMaximum(), h_gamma.GetMaximum())
        h_alpha.SetMaximum(max_val * 1.2)
        h_alpha.Draw("HIST")
        h_gamma.Draw("HIST SAME")

        # Legend
        leg = plot_utils.CreateLegend(0.55, 0.75, 0.9, 0.85)
        leg.AddEntry(h_alpha, f"Alpha (n={len(alpha_scores)})", "f")
        leg.AddEntry(h_gamma, f"Gamma/Beta (n={len(gamma_scores)})", "f")
        leg.Draw()

        # Add subplot label and energy range
        plot_utils.AddSubplotLabel(f"({chr(97+i)})")

        # Add energy range text
        energy_text = ROOT.TLatex(0.55, 0.65, "1000-1200 keVee")
        energy_text.SetNDC()
        energy_text.SetTextSize(0.04)
        energy_text.Draw()

    canvas.SaveAs(os.path.join(output_dir, "method_histograms_1000_1200_keVee.pdf"))
    canvas.Close()

    print(
        f"Method histograms (1000-1200 keVee) saved to {output_dir}/method_histograms_1000_1200_keVee.pdf"
    )


def plot_score_histogram_root(alpha_scores, gamma_scores, title, output_path):
    """Plot score histogram using ROOT"""
    canvas = ROOT.TCanvas("c_scores", title, 1200, 900)
    plot_utils.ConfigureCanvas(canvas, logy=True)

    # Determine range
    all_scores = np.concatenate([alpha_scores, gamma_scores])
    score_min = np.min(all_scores)
    score_max = np.max(all_scores)

    # Create histograms
    h_alpha = ROOT.TH1F("h_alpha", "", 75, score_min, score_max)
    h_gamma = ROOT.TH1F("h_gamma", "", 75, score_min, score_max)

    # Fill histograms
    for val in alpha_scores:
        h_alpha.Fill(val)
    for val in gamma_scores:
        h_gamma.Fill(val)

    # Configure
    plot_utils.ConfigureHistogram(h_alpha, ROOT.kRed + 1)
    plot_utils.ConfigureHistogram(h_gamma, ROOT.kGreen + 2)

    h_alpha.GetXaxis().SetTitle("Regressor Output")
    h_alpha.GetYaxis().SetTitle("Counts")
    h_alpha.SetTitle("")

    # Draw
    max_val = max(h_alpha.GetMaximum(), h_gamma.GetMaximum())
    h_alpha.SetMaximum(max_val * 1.2)
    h_alpha.Draw("HIST")
    h_gamma.Draw("HIST SAME")

    # Legend
    leg = plot_utils.CreateLegend(0.65, 0.7, 0.92, 0.85)
    leg.AddEntry(h_alpha, f"Am-241 (#alpha)", "f")
    leg.AddEntry(h_gamma, f"Na-22 (#gamma)", "f")
    leg.Draw()

    canvas.SaveAs(output_path)
    canvas.Close()

def plot_unified_roc_curves_root(
    test_features, y_true, methods, method_names, output_dir="psd_analysis"
):
    """Plot ROC curves for all methods using ROOT"""
    os.makedirs(output_dir, exist_ok=True)
    colors = [ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2]
    target_fpr = 0.05

    canvas = ROOT.TCanvas("c_unified_roc", "Unified ROC", 1200, 900)
    plot_utils.ConfigureCanvas(canvas)

    roc_graphs = []
    # Center the legend at the top - adjust these coordinates
    leg = plot_utils.CreateLegend(0.37, 0.2, 0.92, 0.4)
    leg.SetMargin(0.1)

    for i, (method, name) in enumerate(zip(methods, method_names)):
        scores = test_features[method].values

        # Fix for si_psd: invert scores since higher values correspond to negative class
        if method == "si_psd":
            scores_to_use = -scores  # Invert the scores
        else:
            scores_to_use = scores

        fpr, tpr, thresholds = roc_curve(y_true, scores_to_use)

        # Find threshold and TPR at 5% FPR
        index = np.argmin(np.abs(fpr - target_fpr))
        threshold_at_5pct_fpr = thresholds[index]
        tpr_at_5pct_fpr = tpr[index]
        actual_fpr = fpr[index]

        print(f"{name}:")
        if method == "si_psd":
            # For si_psd, show both the inverted threshold and original threshold
            original_threshold = -threshold_at_5pct_fpr  # Convert back to original scale
            print(f"  Inverted threshold at {actual_fpr:.3f} FPR: {threshold_at_5pct_fpr:.6f}")
            print(f"  Original threshold (lower is better): {original_threshold:.6f}")
        else:
            print(f"  Threshold at {actual_fpr:.3f} FPR: {threshold_at_5pct_fpr:.6f}")
        print(f"  TPR at {actual_fpr:.3f} FPR: {tpr_at_5pct_fpr:.3f}")

        auc_score = auc(fpr, tpr)

        # Create TGraph
        roc_graph = ROOT.TGraph(len(fpr), fpr, tpr)
        roc_graph.SetLineColor(colors[i])
        roc_graph.SetLineWidth(3)
        roc_graphs.append(roc_graph)

        # Draw
        if i == 0:
            roc_graph.SetTitle("")
            roc_graph.GetXaxis().SetTitle("False Positive Rate (1 - Specificity)")
            roc_graph.GetYaxis().SetTitle("True Positive Rate (Sensitivity)")
            roc_graph.GetXaxis().SetRangeUser(0, 1)
            roc_graph.GetYaxis().SetRangeUser(0, 1)
            roc_graph.Draw("AL")
        else:
            roc_graph.Draw("L SAME")

        leg.AddEntry(roc_graph, f"{name} (AUC = {auc_score:.2f})", "l")

    # Add diagonal line
    diagonal = ROOT.TLine(0, 0, 1, 1)
    diagonal.SetLineColor(ROOT.kBlack)
    diagonal.SetLineStyle(2)
    diagonal.Draw()

    leg.Draw()
    canvas.SaveAs(f"{output_dir}/unified_roc_curves.pdf")
    canvas.SaveAs(f"{output_dir}/unified_roc_curves.png")
