"""
Proof-of-Concept: Apply trained ML models and Charge Comparison to all sources

Loads all trained models (from regressors.py config) and applies them to every
available source.  Classification thresholds are calibrated on the pure Am-241
(alpha) and Na-22 (gamma) test data — training events are excluded — so that a
fixed true positive rate for alphas / gammas is retained.  The same
thresholds are then applied to the mixed datasets.  Charge Comparison is
included as a baseline for comparison.

Training was done on:
- Am-241 (alpha)
- Na-22 (gamma)
"""

import numpy as np
import os
import pickle
import ROOT
from analysis_utils.io import load_tree_data
from psd_utils import process_waveforms, column_name, get_training_indices
from regressors import get_default_regressors

import analysis_utils

analysis_utils.load_cpp_library()

ROOT.gROOT.SetBatch(True)
ROOT.PlottingUtils.SetStylePreferences(ROOT.PlotSaveFormat.kPNG)

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

SOURCE_MAP = {
    "Am-241 & Cs-137": "Am241Cs137.root",
    "Am-241 & Co-60": "Am241Co60.root",
}

TRUE_POSITIVE_RATE = 99.9


def _clean_name(source_name):
    return source_name.replace(" ", "_").replace("&", "and").replace("-", "_")


def determine_thresholds(alpha_scores, gamma_scores, tpr):
    """Find thresholds from pure-source test data at a given true positive rate.

    Lower scores are alpha-like, higher scores are gamma-like (matches
    regressor convention alpha=0, gamma=1, and charge_comparison convention
    used in the ROC analysis).

    Returns (alpha_upper, gamma_lower).
    """
    alpha_upper = np.percentile(alpha_scores, tpr)
    gamma_lower = np.percentile(gamma_scores, 100 - tpr)
    return alpha_upper, gamma_lower


def classify_events(scores, alpha_threshold, gamma_threshold):
    """Return (alpha_mask, gamma_mask) boolean arrays.

    Events satisfying both criteria (overlapping thresholds) are uncertain.
    """
    alpha_like = scores <= alpha_threshold
    gamma_like = scores >= gamma_threshold
    both = alpha_like & gamma_like
    return alpha_like & ~both, gamma_like & ~both


def plot_classified_spectra(features, source_name, method_name, alpha_mask,
                            gamma_mask):
    """Plot light-output spectra classified by thresholds.  Returns stats dict."""
    all_lo = features["light_output"].values
    alpha_like_lo = all_lo[alpha_mask]
    gamma_like_lo = all_lo[gamma_mask]

    total = len(all_lo)
    n_alpha = int(np.sum(alpha_mask))
    n_gamma = int(np.sum(gamma_mask))
    n_uncertain = total - n_alpha - n_gamma

    canvas = ROOT.TCanvas(str(ROOT.PlottingUtils.GetRandomName()), "", 1200,
                          600)

    pad_plot = ROOT.TPad("pad_plot", "", 0.0, 0.0, 0.7, 1.0)
    pad_leg = ROOT.TPad("pad_leg", "", 0.7, 0.0, 1.0, 1.0)
    pad_plot.SetLogy()
    pad_plot.SetLeftMargin(0.12)
    pad_plot.SetRightMargin(0.02)
    pad_leg.SetLeftMargin(0.0)
    pad_leg.SetRightMargin(0.05)
    pad_plot.Draw()
    pad_leg.Draw()

    lo_min, lo_max = 0, 2000
    nbins = 200

    h_all = ROOT.TH1F(str(ROOT.PlottingUtils.GetRandomName()), "", nbins,
                      lo_min, lo_max)
    h_alpha = ROOT.TH1F(str(ROOT.PlottingUtils.GetRandomName()), "", nbins,
                        lo_min, lo_max)
    h_gamma = ROOT.TH1F(str(ROOT.PlottingUtils.GetRandomName()), "", nbins,
                        lo_min, lo_max)

    for val in all_lo:
        h_all.Fill(val)
    for val in alpha_like_lo:
        h_alpha.Fill(val)
    for val in gamma_like_lo:
        h_gamma.Fill(val)

    ROOT.PlottingUtils.ConfigureHistogram(h_all, ROOT.kBlack)
    ROOT.PlottingUtils.ConfigureHistogram(h_alpha, ROOT.kRed + 1)
    ROOT.PlottingUtils.ConfigureHistogram(h_gamma, ROOT.kBlue + 1)

    h_all.SetLineStyle(2)
    h_all.GetXaxis().SetTitle("Light Output [keVee]")
    h_all.GetYaxis().SetTitle("Counts / 10 keVee")
    h_all.GetYaxis().SetTitleOffset(1)
    h_all.SetTitle("")

    pad_plot.cd()
    max_val = h_all.GetMaximum()
    h_all.SetMaximum(max_val * 1.2)
    h_all.Draw("HIST")
    h_alpha.Draw("HIST SAME")
    h_gamma.Draw("HIST SAME")
    text = ROOT.PlottingUtils.AddSubplotLabel(method_name, 0.92, 0.84)

    pad_leg.cd()
    leg = ROOT.PlottingUtils.AddLegend(0.05, 0.95, 0.3, 0.7)
    leg.SetMargin(0.22)
    leg.AddEntry(h_all, f"{source_name}", "l")
    leg.AddEntry(h_alpha, "#alpha-like", "f")
    leg.AddEntry(h_gamma, "#gamma-like", "f")
    leg.Draw()

    safe_reg = method_name.replace(" ", "_").lower()
    safe_src = _clean_name(source_name)
    output_path = f"classified_spectra_{safe_src}_{safe_reg}"
    canvas.cd()
    ROOT.PlottingUtils.SaveFigure(canvas, output_path,
                                  ROOT.PlotSaveOptions.kLOG)
    canvas.Close()
    h_all.Delete()
    h_alpha.Delete()
    h_gamma.Delete()

    print(
        f"  [{method_name}] alpha-like: {n_alpha} ({n_alpha/total*100:.1f}%), "
        f"gamma-like: {n_gamma} ({n_gamma/total*100:.1f}%), "
        f"uncertain: {n_uncertain} ({n_uncertain/total*100:.1f}%)")

    return {
        "source": source_name,
        "method": method_name,
        "total": total,
        "n_alpha": n_alpha,
        "n_gamma": n_gamma,
        "n_uncertain": n_uncertain,
    }


def save_classification_table(all_stats,
                              tpr,
                              output_path="classification_table.txt"):
    """Save classification statistics as a LaTeX table."""
    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"  \centering")
    lines.append(r"  \caption{Event classification at "
                 f"{tpr}"
                 r"\% true positive rate for $\alpha$ and $\gamma$ events.}")
    lines.append(r"  \label{tab:classification}")
    lines.append(r"  \begin{tabular}{llrrrr}")
    lines.append(r"    \toprule")
    lines.append(r"    Source & Method & Total"
                 r" & $\alpha$-like [\%]"
                 r" & $\gamma$-like [\%]"
                 r" & Uncertain [\%] \\")
    lines.append(r"    \midrule")

    current_source = None
    for s in all_stats:
        if s["source"] != current_source:
            if current_source is not None:
                lines.append(r"    \midrule")
            current_source = s["source"]
            src_label = s["source"]
        else:
            src_label = ""

        pct_a = s["n_alpha"] / s["total"] * 100
        pct_g = s["n_gamma"] / s["total"] * 100
        pct_u = s["n_uncertain"] / s["total"] * 100
        lines.append(f"    {src_label} & {s['method']}"
                     f" & {s['total']}"
                     f" & {pct_a:.1f}"
                     f" & {pct_g:.1f}"
                     f" & {pct_u:.1f}"
                     r" \\")

    lines.append(r"    \bottomrule")
    lines.append(r"  \end{tabular}")
    lines.append(r"\end{table}")

    table_str = "".join(lines)
    with open(output_path, "w") as f:
        f.write(table_str + "")
    print(f"LaTeX classification table saved to {output_path}")
    print(table_str)


def main():
    os.makedirs("plots/", exist_ok=True)

    regressor_cfgs = get_default_regressors()
    models = {}
    print("Loading trained models...")
    for cfg in regressor_cfgs:
        if os.path.exists(cfg["file"]):
            with open(cfg["file"], "rb") as f:
                models[cfg["name"]] = pickle.load(f)
            print(f"  Loaded: {cfg['name']} ({cfg['file']})")
        else:
            print(
                f"  WARNING: {cfg['file']} not found, skipping {cfg['name']}")

    if not models:
        raise FileNotFoundError(
            "No trained model files found. Run analysis.py first.")

    threshold_cache = "thresholds.npz"
    if os.path.exists(threshold_cache):
        print(f"Loading cached thresholds from {threshold_cache}...")
        data = np.load(threshold_cache, allow_pickle=True)
        thresholds = data["thresholds"].item()
        for name, (at, gt) in thresholds.items():
            print(f"  {name}: alpha <= {at:.4f}, gamma >= {gt:.4f}")
    else:
        print("Loading calibration data...")
        am_features, am_waveforms = load_tree_data(
            ROOT_FILES_DIR + "Am241.root",
            scalar_branches=SCALAR_BRANCHES,
            array_branch="Samples",
        )
        na_features, na_waveforms = load_tree_data(
            ROOT_FILES_DIR + "Na22.root",
            scalar_branches=SCALAR_BRANCHES,
            array_branch="Samples",
        )

        alpha_train_idx, gamma_train_idx = get_training_indices(
            am_features, na_features)

        am_test_mask = np.ones(len(am_features), dtype=bool)
        am_test_mask[alpha_train_idx] = False
        na_test_mask = np.ones(len(na_features), dtype=bool)
        na_test_mask[gamma_train_idx] = False

        am_test_feat = am_features[am_test_mask].reset_index(drop=True)
        na_test_feat = na_features[na_test_mask].reset_index(drop=True)
        am_test_wf = am_waveforms[am_test_mask]
        na_test_wf = na_waveforms[na_test_mask]

        # Apply light output masks matching the training range
        am_lo_mask = ((am_test_feat["light_output"] >= 375)
                      & (am_test_feat["light_output"] <= 1575))
        na_lo_mask = ((na_test_feat["light_output"] >= 0)
                      & (na_test_feat["light_output"] <= 1750))

        am_cal_feat = am_test_feat[am_lo_mask].reset_index(drop=True)
        na_cal_feat = na_test_feat[na_lo_mask].reset_index(drop=True)
        am_cal_processed = process_waveforms(am_test_wf[am_lo_mask.values])
        na_cal_processed = process_waveforms(na_test_wf[na_lo_mask.values])

        print(f"Am-241 calibration events: {len(am_cal_feat)} "
              f"(of {len(am_test_feat)} test)")
        print(f"Na-22  calibration events: {len(na_cal_feat)} "
              f"(of {len(na_test_feat)} test)")

        # Thresholds for ML regressors
        print(f"Determining thresholds at {TRUE_POSITIVE_RATE}% TPR...")
        thresholds = {}
        for name, model in models.items():
            am_scores = model.predict(am_cal_processed)
            na_scores = model.predict(na_cal_processed)
            alpha_thresh, gamma_thresh = determine_thresholds(
                am_scores, na_scores, TRUE_POSITIVE_RATE)
            thresholds[name] = (alpha_thresh, gamma_thresh)
            print(f"  {name}: alpha <= {alpha_thresh:.4f}, "
                  f"gamma >= {gamma_thresh:.4f}")

        # Thresholds for Charge Comparison
        am_cc = am_cal_feat["charge_comparison"].values
        na_cc = na_cal_feat["charge_comparison"].values
        cc_alpha_thresh, cc_gamma_thresh = determine_thresholds(
            am_cc, na_cc, TRUE_POSITIVE_RATE)
        thresholds["Charge Comparison"] = (cc_alpha_thresh, cc_gamma_thresh)
        print(f"  Charge Comparison: alpha <= {cc_alpha_thresh:.4f}, "
              f"gamma >= {cc_gamma_thresh:.4f}")

        np.savez(threshold_cache, thresholds=thresholds)
        print(f"Thresholds cached to {threshold_cache}")

    method_names = list(models.keys()) + ["Charge Comparison"]
    all_stats = []

    for source_name, root_file in SOURCE_MAP.items():
        print(f"Analyzing: {source_name}")
        source_path = ROOT_FILES_DIR + root_file
        if not os.path.exists(source_path):
            print(f"  WARNING: File not found: {source_path}")
            continue

        features, waveforms = load_tree_data(source_path,
                                             scalar_branches=SCALAR_BRANCHES,
                                             array_branch="Samples",
                                             max_events=None)

        if len(features) == 0:
            print(f"  WARNING: No events found")
            continue

        print(f"  Events: {len(features)}")

        processed = process_waveforms(waveforms)

        for name, model in models.items():
            scores = model.predict(processed)
            alpha_thresh, gamma_thresh = thresholds[name]
            alpha_mask, gamma_mask = classify_events(scores, alpha_thresh,
                                                     gamma_thresh)
            stats = plot_classified_spectra(features, source_name, name,
                                            alpha_mask, gamma_mask)
            all_stats.append(stats)

        cc_scores = features["charge_comparison"].values
        alpha_thresh, gamma_thresh = thresholds["Charge Comparison"]
        alpha_mask, gamma_mask = classify_events(cc_scores, alpha_thresh,
                                                 gamma_thresh)
        stats = plot_classified_spectra(features, source_name,
                                        "Charge Comparison", alpha_mask,
                                        gamma_mask)
        all_stats.append(stats)

    save_classification_table(all_stats, TRUE_POSITIVE_RATE)

    print(f"Complete. Analyzed {len(SOURCE_MAP)} sources "
          f"with {len(method_names)} methods.")
    print(f"Plots saved to: plots")


if __name__ == "__main__":
    main()
