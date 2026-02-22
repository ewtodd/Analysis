import os
import time
import ROOT
import numpy as np
import pandas as pd
import pickle
from sklearn.metrics import roc_curve, auc
from joblib import Parallel, delayed, parallel_backend
import analysis_utils

analysis_utils.load_cpp_library()
ROOT.gROOT.SetBatch(True)
ROOT.PlottingUtils.SetStylePreferences(ROOT.PlotSaveFormat.kPNG)

N_JOBS = 32
N_BOOTSTRAP = 250


def _bootstrap_chunk(y_true, scores, indices):
    """Compute AUC for a chunk of bootstrap resamples."""
    aucs = []
    for idx in indices:
        y_b = y_true[idx]
        s_b = scores[idx]
        if len(np.unique(y_b)) < 2:
            continue
        fpr_b, tpr_b, _ = roc_curve(y_b, s_b)
        aucs.append(auc(fpr_b, tpr_b))
    return aucs


def bootstrap_auc(y_true,
                  scores,
                  n_bootstrap=N_BOOTSTRAP,
                  random_state=42,
                  n_jobs=N_JOBS):
    """Estimate AUC and its uncertainty via bootstrap resampling.

    Returns
    -------
    auc_mean : float
    auc_std : float
    """
    rng = np.random.RandomState(random_state)
    n = len(y_true)
    all_idx = rng.randint(0, n, size=(n_bootstrap, n))
    chunks = np.array_split(all_idx, n_jobs)
    with parallel_backend('threading', n_jobs=n_jobs):
        results = Parallel()(delayed(_bootstrap_chunk)(y_true, scores, chunk)
                             for chunk in chunks)
    aucs = np.array([a for chunk_aucs in results for a in chunk_aucs])
    return float(np.mean(aucs)), float(np.std(aucs))


def process_waveforms(waveforms, n_jobs=N_JOBS):
    """Normalize each waveform by its maximum value.

    Parameters
    ----------
    waveforms : numpy.ndarray
        2-D array of shape (n_waveforms, n_samples).
    n_jobs : int
        Unused, kept for interface compatibility.

    Returns
    -------
    numpy.ndarray
        Normalized waveforms.
    """
    maxvals = np.max(waveforms, axis=1, keepdims=True)
    maxvals[maxvals == 0] = 1.0
    return waveforms / maxvals


def column_name(regressor_name):
    """Convert a regressor display name to a DataFrame column name."""
    return regressor_name.replace(" ", "_") + "_Output"


def train_or_load(regressor_cfg, x_train, y_train):
    """Train a regressor or load from cache. Returns (model, train_time_s or None)."""
    model_file = regressor_cfg["file"]
    name = regressor_cfg["name"]

    if os.path.exists(model_file):
        print(f"Loading existing model: {name} ({model_file})")
        with open(model_file, "rb") as f:
            return pickle.load(f), None

    print(f"Training new model: {name}...")
    model = regressor_cfg["model"]
    t0 = time.perf_counter()
    model.fit(x_train, y_train)
    train_time = time.perf_counter() - t0
    with open(model_file, "wb") as f:
        pickle.dump(model, f)
    return model, train_time


def get_training_indices(alpha_features, gamma_features, random_state=42):
    """Reconstruct which event indices were selected for training.

    Mirrors the selection logic in regress_waveforms so that downstream code
    (e.g. proof_of_concept) can exclude training events deterministically.

    Returns
    -------
    alpha_train_idx, gamma_train_idx : numpy arrays of original-index positions
    """
    lo_lower = {"alpha": 375, "gamma": 0}
    lo_upper = {"alpha": 1575, "gamma": 1750}

    lo_mask_alpha = ((alpha_features["light_output"] <= lo_upper["alpha"])
                     & (alpha_features["light_output"] >= lo_lower["alpha"]))
    lo_mask_gamma = ((gamma_features["light_output"] <= lo_upper["gamma"])
                     & (gamma_features["light_output"] >= lo_lower["gamma"]))

    alpha_mask_idx = np.where(lo_mask_alpha.values)[0]
    gamma_mask_idx = np.where(lo_mask_gamma.values)[0]

    min_samples = 5000
    rng = np.random.RandomState(random_state)
    alpha_sel = rng.choice(len(alpha_mask_idx),
                           size=min_samples,
                           replace=False)
    gamma_sel = rng.choice(len(gamma_mask_idx),
                           size=min_samples,
                           replace=False)

    return alpha_mask_idx[alpha_sel], gamma_mask_idx[gamma_sel]


def save_timing_table(timing_records, n_cores):
    """Save training and inference timing data to a LaTeX table."""
    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"  \centering")
    lines.append(r"  \caption{Training and inference time for each regressor"
                 f" ({n_cores} cores)."
                 r"}")
    lines.append(r"  \label{tab:timing}")
    lines.append(r"  \begin{tabular}{lcccc}")
    lines.append(r"    \toprule")
    lines.append(r"    Method & Train [s] & Train/core [s] "
                 r"& Infer [$\mu$s/event] & Infer/core [$\mu$s/event] \\")
    lines.append(r"    \midrule")
    for r in timing_records:
        train = f"{r['train_time']:.3f}" if r.get(
            "train_time") is not None else "--"
        train_pc = (f"{r['train_time_per_core']:.3f}"
                    if r.get("train_time_per_core") is not None else "--")
        infer = (f"{r['infer_per_event_us']:.1f}"
                 if r.get("infer_per_event_us") is not None else "--")
        infer_pc = (f"{r['infer_per_event_per_core_us']:.1f}" if
                    r.get("infer_per_event_per_core_us") is not None else "--")
        lines.append(f"    {r['name']} & {train} & {train_pc}"
                     f" & {infer} & {infer_pc}"
                     r" \\")
    lines.append(r"    \bottomrule")
    lines.append(r"  \end{tabular}")
    lines.append(r"\end{table}")

    table_str = "\n".join(lines)
    output_path = "timing_table.txt"
    with open(output_path, "w") as f:
        f.write(table_str + "\n")
    print(f"\nLaTeX timing table saved to {output_path}")
    print(table_str)


def regress_waveforms(waveforms,
                      features,
                      regressors,
                      process_func=process_waveforms,
                      random_state=42):
    """Train and evaluate ML models using waveforms from ROOT files.

    Parameters
    ----------
    waveforms : tuple of (alpha_waveforms_ndarray, gamma_waveforms_ndarray)
    features : tuple of (alpha_features_df, gamma_features_df)
    regressors : list of dict
        Each dict has keys "name", "model" (unfitted sklearn regressor),
        and "file" (path for caching the trained model).
    """
    lo_lower_dict = {"alpha": 375, "gamma": 0}
    lo_upper_dict = {"alpha": 1575, "gamma": 1750}

    alpha_waveforms, gamma_waveforms = waveforms
    alpha_features, gamma_features = features

    alpha_train_original_indices, gamma_train_original_indices = \
        get_training_indices(alpha_features, gamma_features, random_state)

    train_alpha_wf = alpha_waveforms[alpha_train_original_indices]
    train_gamma_wf = gamma_waveforms[gamma_train_original_indices]

    min_samples = len(alpha_train_original_indices)
    print(f"Samples for balanced training: {min_samples}")

    # Process training waveforms
    with parallel_backend('threading', n_jobs=2):
        train_results = Parallel()(
            delayed(process_func)(group)
            for group in [train_alpha_wf, train_gamma_wf])

    plot_sample_waveforms(train_results,
                          n_samples=5,
                          random_state=random_state)

    x_train = np.vstack(train_results)
    y_train = np.array([0] * len(train_results[0]) +
                       [1] * len(train_results[1]))

    # Train or load each regressor
    n_cores = 32
    trained_models = {}
    timing_records = []
    for cfg in regressors:
        model, train_time = train_or_load(cfg, x_train, y_train)
        trained_models[cfg["name"]] = model
        if train_time is not None:
            print(f"  Training time: {train_time:.3f} s "
                  f"({train_time / n_cores:.3f} s/core)")
            timing_records.append({
                "name": cfg["name"],
                "train_time": train_time,
                "train_time_per_core": train_time / n_cores,
            })

    # Create test data by dropping training samples
    all_alpha_idx = np.arange(len(alpha_waveforms))
    test_alpha_idx = np.setdiff1d(all_alpha_idx, alpha_train_original_indices)
    test_alpha_wf = alpha_waveforms[test_alpha_idx]
    test_alpha_features = alpha_features.iloc[test_alpha_idx].reset_index(
        drop=True)

    all_gamma_idx = np.arange(len(gamma_waveforms))
    test_gamma_idx = np.setdiff1d(all_gamma_idx, gamma_train_original_indices)
    test_gamma_wf = gamma_waveforms[test_gamma_idx]
    test_gamma_features = gamma_features.iloc[test_gamma_idx].reset_index(
        drop=True)

    # Process test waveforms for prediction
    with parallel_backend('threading', n_jobs=2):
        test_results = Parallel()(delayed(process_func)(group)
                                  for group in [test_alpha_wf, test_gamma_wf])

    X_test = np.vstack(test_results)

    # Predict with each regressor and add output columns
    test_alpha_features = test_alpha_features.copy()
    test_gamma_features = test_gamma_features.copy()

    n_test = len(X_test)
    regressor_names = []
    for name, model in trained_models.items():
        col = column_name(name)
        t0 = time.perf_counter()
        y_pred = model.predict(X_test)
        infer_time = time.perf_counter() - t0
        infer_per_event = infer_time / n_test
        print(f"{name}: inference {infer_time:.3f} s total, "
              f"{infer_per_event * 1e6:.1f} us/event, "
              f"{infer_per_event / n_cores * 1e6:.1f} us/event/core")
        # Update timing record or create one (if model was loaded from cache)
        record = next((r for r in timing_records if r["name"] == name), None)
        if record is None:
            timing_records.append({
                "name":
                name,
                "train_time":
                None,
                "train_time_per_core":
                None,
                "infer_time":
                infer_time,
                "infer_per_event_us":
                infer_per_event * 1e6,
                "infer_per_event_per_core_us":
                infer_per_event / n_cores * 1e6,
            })
        else:
            record["infer_time"] = infer_time
            record["infer_per_event_us"] = infer_per_event * 1e6
            record[
                "infer_per_event_per_core_us"] = infer_per_event / n_cores * 1e6
        test_alpha_features[col] = y_pred[:len(test_alpha_wf)]
        test_gamma_features[col] = y_pred[len(test_alpha_wf):]
        regressor_names.append(name)

    # Only save timing table if we actually trained (otherwise we'd overwrite
    # the table with missing training times)
    if any(r.get("train_time") is not None for r in timing_records):
        save_timing_table(timing_records, n_cores)

    # Apply light output filters for plotting score histograms
    alpha_lo_mask = (
        test_alpha_features["light_output"] <= lo_upper_dict["alpha"]) & (
            test_alpha_features["light_output"] >= lo_lower_dict["alpha"])

    gamma_lo_mask = (
        test_gamma_features["light_output"] <= lo_upper_dict["gamma"]) & (
            test_gamma_features["light_output"] >= lo_lower_dict["gamma"])

    test_alpha_features_filtered = test_alpha_features[alpha_lo_mask]
    test_gamma_features_filtered = test_gamma_features[gamma_lo_mask]

    for name in regressor_names:
        col = column_name(name)
        safe_name = name.replace(" ", "_").lower()
        plot_score_histogram(
            test_alpha_features_filtered[col].values,
            test_gamma_features_filtered[col].values,
            f"Test Set Scores ({name})",
            f"test_score_histogram_{safe_name}",
            regressor_name=name,
        )

    alpha_lo_mask_900_1200 = (test_alpha_features["light_output"] <= 1200) & (
        test_alpha_features["light_output"] >= 900)

    gamma_lo_mask_900_1200 = (test_gamma_features["light_output"] <= 1200) & (
        test_gamma_features["light_output"] >= 900)

    test_alpha_features_900_1200 = test_alpha_features[alpha_lo_mask_900_1200]
    test_gamma_features_900_1200 = test_gamma_features[gamma_lo_mask_900_1200]

    for name in regressor_names:
        col = column_name(name)
        safe_name = name.replace(" ", "_").lower()
        plot_score_histogram(
            test_alpha_features_900_1200[col].values,
            test_gamma_features_900_1200[col].values,
            f"Test Set Scores ({name})",
            f"test_score_histogram_900_1200_{safe_name}",
            regressor_name=name,
        )

    analyze_all_methods(test_alpha_features_filtered,
                        test_gamma_features_filtered, regressor_names, "full")

    plot_auc_vs_light_output(test_alpha_features, test_gamma_features,
                             regressor_names)

    return (
        (test_alpha_wf, test_gamma_wf),
        (test_alpha_features, test_gamma_features),
    )


def analyze_all_methods(test_alpha_features, test_gamma_features,
                        regressor_names, output_name):
    """ROC analysis comparing ML regressors, Charge Comparison, and Shape Indicator PSD."""

    test_features = pd.concat([test_alpha_features,
                               test_gamma_features]).reset_index(drop=True)

    y_true = np.array([0] * len(test_alpha_features) +
                      [1] * len(test_gamma_features))

    all_methods = [column_name(n) for n in regressor_names]
    all_method_names = list(regressor_names)

    all_methods += ["charge_comparison", "clean_shape_indicator"]
    all_method_names += ["Charge Comparison", "Shape Indicator"]

    plot_roc_curves(test_features, y_true, all_methods, all_method_names,
                    output_name)


def plot_auc_vs_light_output(test_alpha_features, test_gamma_features,
                             regressor_names):
    """Plot ROC AUC as a function of light output for all classifiers."""
    # Non-uniform bins: wider where alpha events are sparse, narrower where abundant
    bin_edges = ([375, 575, 775, 975, 1075] + list(range(1125, 1575 + 1, 50)))
    lo_bins = [(bin_edges[i], bin_edges[i + 1])
               for i in range(len(bin_edges) - 1)]
    bin_centers = [0.5 * (lo + hi) for lo, hi in lo_bins]

    all_methods = [column_name(n) for n in regressor_names]
    all_method_names = list(regressor_names)
    all_methods += ["charge_comparison", "clean_shape_indicator"]
    all_method_names += ["Charge Comparison", "Shape Indicator"]

    # Compute AUC per bin per method with bootstrap uncertainties
    auc_results = {name: [] for name in all_method_names}
    auc_errors = {name: [] for name in all_method_names}

    for lo, hi in lo_bins:
        alpha_mask = (test_alpha_features["light_output"]
                      >= lo) & (test_alpha_features["light_output"] < hi)
        gamma_mask = (test_gamma_features["light_output"]
                      >= lo) & (test_gamma_features["light_output"] < hi)

        alpha_bin = test_alpha_features[alpha_mask]
        gamma_bin = test_gamma_features[gamma_mask]

        n_alpha = len(alpha_bin)
        n_gamma = len(gamma_bin)
        print(
            f"Light output bin [{lo}, {hi}) keVee: {n_alpha} alpha, {n_gamma} gamma"
        )

        if n_alpha < 5 or n_gamma < 5:
            for name in all_method_names:
                auc_results[name].append(np.nan)
                auc_errors[name].append(np.nan)
            continue

        combined = pd.concat([alpha_bin, gamma_bin]).reset_index(drop=True)
        y_true = np.array([0] * n_alpha + [1] * n_gamma)

        for method, name in zip(all_methods, all_method_names):
            scores = combined[method].values
            if method in ("raw_shape_indicator", "clean_shape_indicator"):
                scores = -scores
            auc_mean, auc_std = bootstrap_auc(y_true, scores)
            auc_results[name].append(auc_mean)
            auc_errors[name].append(auc_std)

    # Plot
    colors = [
        ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kOrange + 1,
        ROOT.kMagenta + 1, ROOT.kCyan + 1, ROOT.kYellow + 2, ROOT.kViolet + 1
    ]

    canvas = ROOT.TCanvas("auc_lo", "", 1200, 600)

    pad_plot = ROOT.TPad("pad_plot", "", 0.0, 0.0, 0.8, 1.0)
    pad_plot.SetLeftMargin(0.15)
    pad_plot.SetRightMargin(0.02)
    pad_plot.Draw()

    pad_leg = ROOT.TPad("pad_leg", "", 0.8, 0.0, 1.0, 1.0)
    pad_leg.SetLeftMargin(0.0)
    pad_leg.SetRightMargin(0.05)
    pad_leg.Draw()

    pad_plot.cd()

    x_arr = np.array(bin_centers, dtype=np.float64)
    graphs = []

    ex_arr = np.array([0.5 * (hi - lo) for lo, hi in lo_bins],
                      dtype=np.float64)

    for i, name in enumerate(all_method_names):
        y_arr = np.array(auc_results[name], dtype=np.float64)
        ey_arr = np.array(auc_errors[name], dtype=np.float64)
        valid = ~np.isnan(y_arr)
        if not np.any(valid):
            continue

        graph = ROOT.TGraphErrors(int(np.sum(valid)), x_arr[valid],
                                  y_arr[valid], ex_arr[valid], ey_arr[valid])
        graph.SetLineColor(colors[i % len(colors)])
        graph.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
        graph.SetMarkerColor(colors[i % len(colors)])
        graph.SetMarkerStyle(20 + i)
        graph.SetMarkerSize(1.2)
        graphs.append(graph)

        if len(graphs) == 1:
            graph.SetTitle("")
            graph.GetXaxis().SetTitle("Light Output [keVee]")
            graph.GetYaxis().SetTitle("ROC AUC")
            graph.GetXaxis().SetRangeUser(400, 2000)
            graph.GetYaxis().SetRangeUser(0.3, 1.05)
            graph.Draw("APL")
        else:
            graph.Draw("PL SAME")

    pad_leg.cd()
    leg = ROOT.PlottingUtils.AddLegend(0.05, 0.95, 0.2, 0.9)
    leg.SetMargin(0.15)
    leg.SetTextSize(0.1)
    for i, name in enumerate(all_method_names):
        if i < len(graphs):
            leg.AddEntry(graphs[i], name, "lp")

    ROOT.PlottingUtils.SaveFigure(canvas, "auc_vs_light_output",
                                  ROOT.PlotSaveOptions.kLINEAR)
    canvas.Close()

    print("AUC vs light output plot saved to auc_vs_light_output")


def plot_score_histogram(alpha_scores,
                         gamma_scores,
                         title,
                         output_path,
                         regressor_name="Regressor"):
    """Plot score histogram using ROOT"""
    canvas = ROOT.PlottingUtils.GetConfiguredCanvas(ROOT.kTRUE)

    all_scores = np.concatenate([alpha_scores, gamma_scores])
    score_min = np.min(all_scores)
    score_max = np.max(all_scores)

    h_alpha = ROOT.TH1F(str(ROOT.PlottingUtils.GetRandomName()), "", 250,
                        score_min, score_max)
    h_gamma = ROOT.TH1F(str(ROOT.PlottingUtils.GetRandomName()), "", 250,
                        score_min, score_max)

    for val in alpha_scores:
        h_alpha.Fill(val)
    for val in gamma_scores:
        h_gamma.Fill(val)

    ROOT.PlottingUtils.ConfigureHistogram(h_alpha, ROOT.kRed + 1)
    ROOT.PlottingUtils.ConfigureHistogram(h_gamma, ROOT.kGreen + 2)

    h_alpha.GetXaxis().SetTitle(f"{regressor_name} Score")
    h_alpha.GetYaxis().SetTitle("Counts")
    h_alpha.SetTitle("")

    max_val = max(h_alpha.GetMaximum(), h_gamma.GetMaximum())
    h_alpha.SetMaximum(max_val * 1.2)
    h_alpha.Draw("HIST")
    h_gamma.Draw("HIST SAME")

    leg = ROOT.PlottingUtils.AddLegend(0.4, 0.6, 0.7, 0.85)
    leg.AddEntry(h_alpha, f"Am-241 (#alpha)", "f")
    leg.AddEntry(h_gamma, f"Na-22 (#gamma)", "f")
    leg.Draw()

    ROOT.PlottingUtils.SaveFigure(
        canvas,
        output_path,
        ROOT.PlotSaveOptions.kLOG,
    )
    canvas.Close()
    h_alpha.Delete()
    h_gamma.Delete()


def plot_roc_curves(test_features, y_true, methods, method_names, output_name):
    """Plot ROC curves for all methods using ROOT"""
    colors = [
        ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kOrange + 1,
        ROOT.kMagenta + 1, ROOT.kCyan + 1, ROOT.kYellow + 2, ROOT.kViolet + 1
    ]

    canvas = ROOT.PlottingUtils.GetConfiguredCanvas()

    roc_graphs = []
    leg = ROOT.PlottingUtils.AddLegend(0.38, 0.88, 0.2, 0.6)
    leg.SetMargin(0.1)

    for i, (method, name) in enumerate(zip(methods, method_names)):
        scores = test_features[method].values

        if method == "raw_shape_indicator" or method == "clean_shape_indicator":
            scores_to_use = -scores
        else:
            scores_to_use = scores

        fpr, tpr, thresholds = roc_curve(y_true, scores_to_use)

        index = np.argmin(np.abs(fpr - 0.05))
        threshold_at_5pct_fpr = thresholds[index]
        tpr_at_5pct_fpr = tpr[index]
        actual_fpr = fpr[index]

        print(f"{name}:")
        if method == "raw_shape_indicator" or method == "clean_shape_indicator":
            original_threshold = -threshold_at_5pct_fpr
            print(
                f"  Inverted threshold at {actual_fpr:.3f} FPR: {threshold_at_5pct_fpr:.6f}"
            )
            print(
                f"  Original threshold (lower is better): {original_threshold:.6f}"
            )
        else:
            print(
                f"  Threshold at {actual_fpr:.3f} FPR: {threshold_at_5pct_fpr:.6f}"
            )
        print(f"  TPR at {actual_fpr:.3f} FPR: {tpr_at_5pct_fpr:.3f}")

        auc_mean, auc_std = bootstrap_auc(y_true, scores_to_use)

        roc_graph = ROOT.TGraph(len(fpr), fpr, tpr)
        roc_graph.SetLineColor(colors[i % len(colors)])
        roc_graph.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
        roc_graphs.append(roc_graph)

        if i == 0:
            roc_graph.SetTitle("")
            roc_graph.GetXaxis().SetTitle(
                "False Positive Rate (1 - Specificity)")
            roc_graph.GetYaxis().SetTitle("True Positive Rate (Sensitivity)")
            roc_graph.GetXaxis().SetRangeUser(0, 1)
            roc_graph.GetYaxis().SetRangeUser(0, 1)
            roc_graph.Draw("AL")
        else:
            roc_graph.Draw("L SAME")

        leg.AddEntry(roc_graph, f"{name} (AUC = {auc_mean:.3f})", "l")

    diagonal = ROOT.TLine(0, 0, 1, 1)
    diagonal.SetLineColor(ROOT.kBlack)
    diagonal.SetLineStyle(2)
    diagonal.Draw()

    leg.Draw()
    ROOT.PlottingUtils.SaveFigure(
        canvas,
        "roc_curves_" + output_name,
        ROOT.PlotSaveOptions.kLINEAR,
    )


def plot_sample_waveforms(waveforms_tuple, n_samples=5, random_state=42):
    """Plot sample waveforms for alpha and gamma classes after normalization.

    Parameters
    ----------
    waveforms_tuple : tuple of (alpha_waveforms, gamma_waveforms)
        Raw waveform arrays before normalization.
    n_samples : int
        Number of sample waveforms to plot per class.
    random_state : int
        Random seed for reproducibility.
    """
    alpha_waveforms, gamma_waveforms = waveforms_tuple

    rng = np.random.RandomState(random_state)
    alpha_indices = rng.choice(len(alpha_waveforms),
                               size=min(n_samples, len(alpha_waveforms)),
                               replace=False)
    gamma_indices = rng.choice(len(gamma_waveforms),
                               size=min(n_samples, len(gamma_waveforms)),
                               replace=False)

    alpha_samples = alpha_waveforms[alpha_indices]
    gamma_samples = gamma_waveforms[gamma_indices]

    n_points = alpha_samples.shape[1]
    x_values = np.arange(n_points) * 2

    canvas_alpha = ROOT.PlottingUtils.GetConfiguredCanvas()

    graphs_alpha = []
    leg_alpha = ROOT.PlottingUtils.AddLegend(0.65, 0.88, 0.65, 0.85)

    colors = [
        ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kOrange + 1,
        ROOT.kMagenta + 1
    ]

    for i, waveform in enumerate(alpha_samples):
        graph = ROOT.TGraph(n_points, x_values.astype(np.float64),
                            waveform.astype(np.float64))
        graph.SetLineColor(colors[i % len(colors)])
        graph.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())

        if i == 0:
            graph.SetTitle("")
            graph.GetXaxis().SetTitle("Time [ns]")
            graph.GetYaxis().SetTitle("Normalized Amplitude [a.u.]")
            graph.GetXaxis().SetRangeUser(0, x_values[-1])
            graph.GetYaxis().SetRangeUser(-0.1, 1.1)
            graph.Draw("AL")
        else:
            graph.Draw("L SAME")

        graphs_alpha.append(graph)
        leg_alpha.AddEntry(graph, f"Am-241 Waveform {i+1}", "l")

    ROOT.PlottingUtils.SaveFigure(canvas_alpha, "sample_waveforms_alpha",
                                  ROOT.PlotSaveOptions.kLINEAR)
    canvas_alpha.Close()

    canvas_gamma = ROOT.PlottingUtils.GetConfiguredCanvas()

    graphs_gamma = []
    leg_gamma = ROOT.PlottingUtils.AddLegend(0.65, 0.88, 0.65, 0.85)

    for i, waveform in enumerate(gamma_samples):
        graph = ROOT.TGraph(n_points, x_values.astype(np.float64),
                            waveform.astype(np.float64))
        graph.SetLineColor(colors[i % len(colors)])
        graph.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())

        if i == 0:
            graph.SetTitle("")
            graph.GetXaxis().SetTitle("Time [ns]")
            graph.GetYaxis().SetTitle("Normalized Amplitude [a.u.]")
            graph.GetXaxis().SetRangeUser(0, x_values[-1])
            graph.GetYaxis().SetRangeUser(-0.1, 1.1)
            graph.Draw("AL")
        else:
            graph.Draw("L SAME")

        graphs_gamma.append(graph)
        leg_gamma.AddEntry(graph, f"Na-22 Waveform {i+1}", "l")

    ROOT.PlottingUtils.SaveFigure(canvas_gamma, "sample_waveforms_gamma",
                                  ROOT.PlotSaveOptions.kLINEAR)
    canvas_gamma.Close()

    print(
        "Sample waveform plots saved: sample_waveforms_alpha, sample_waveforms_gamma"
    )
