import numpy as np
import os
import pickle
import ROOT
from analysis_utils.io import load_tree_data
import pandas as pd
from psd_utils import (process_waveforms, column_name, ANALYSIS_CACHE_DIR)
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

TRUE_POSITIVE_RATE = 99


def _clean_name(source_name):
    return source_name.replace(" ", "_").replace("&", "and").replace("-", "_")


def _determine_thresholds(alpha_scores, gamma_scores, tpr):
    """Find thresholds from pure-source test data at a given true positive rate.

    Lower scores are alpha-like, higher scores are gamma-like (matches
    regressor convention alpha=0, gamma=1).

    At high TPR the thresholds may overlap; events in the overlap zone are
    excluded by _classify_events.

    Returns (alpha_upper, gamma_lower).
    """
    alpha_upper = np.percentile(alpha_scores, tpr)
    gamma_lower = np.percentile(gamma_scores, 100 - tpr)
    return alpha_upper, gamma_lower


def _classify_events(scores, alpha_threshold, gamma_threshold):
    """Return (alpha_mask, gamma_mask) boolean arrays.

    When thresholds overlap, events in the overlap zone satisfy both
    criteria and are excluded from both masks.
    """
    alpha_like = scores <= alpha_threshold
    gamma_like = scores >= gamma_threshold
    both = alpha_like & gamma_like
    return alpha_like & ~both, gamma_like & ~both


def _plot_classified_spectra(features, source_name, method_name, alpha_mask,
                            gamma_mask, alpha_thresh, gamma_thresh):
    """Plot light-output spectra classified by thresholds."""
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
    ROOT.PlottingUtils.ConfigureHistogram(h_alpha, ROOT.kRed + 2)
    ROOT.PlottingUtils.ConfigureHistogram(h_gamma, ROOT.kBlue + 2)

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
    _ = ROOT.PlottingUtils.AddText(method_name, 0.92, 0.84)

    pad_leg.cd()
    leg = ROOT.PlottingUtils.AddLegend(0.05, 0.95, 0.25, 0.77)
    leg.SetMargin(0.22)
    leg.AddEntry(h_all, f"{source_name}", "l")
    effective_alpha = min(alpha_thresh, gamma_thresh)
    effective_gamma = max(alpha_thresh, gamma_thresh)
    efficiency = (n_alpha + n_gamma) / total * 100
    leg.AddEntry(h_alpha, f"#alpha-like (< {effective_alpha:.3f})", "f")
    leg.AddEntry(h_gamma, f"#gamma-like (> {effective_gamma:.3f})", "f")
    leg.AddEntry(ROOT.nullptr, f"Efficiency: {efficiency:.1f}%", "")
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


def _find_crossover_tpr(am_cal_scores, na_cal_scores):
    """Find the TPR at which the alpha and gamma thresholds cross.

    Below this TPR there is a gap between thresholds (dead zone); above it
    the thresholds overlap.  Returns the crossover TPR value.
    """
    tpr_scan = np.linspace(50, 999, 5000)
    for tpr in tpr_scan:
        alpha_upper, gamma_lower = _determine_thresholds(
            am_cal_scores, na_cal_scores, tpr)
        if alpha_upper >= gamma_lower:
            return tpr
    return tpr_scan[-1]


def _plot_efficiency_vs_strictness(method_names, source_scores, am_cal_feat,
                                  na_cal_feat):
    """Plot classification efficiency vs threshold strictness for each method.

    For each method the crossover TPR (where the alpha and gamma thresholds
    just meet) is found from the calibration data.  This is the point of
    maximum efficiency — all events are uniquely classified with no gap and
    no overlap.

    Strictness is then defined as a linear remap of TPR from the crossover
    point (0 %) to 99 % (100 %), so that increasing strictness always
    corresponds to monotonically decreasing efficiency.

    Parameters
    ----------
    method_names : list of str
    source_scores : dict  {source_name: {method_name: scores_array}}
    am_cal_feat, na_cal_feat : DataFrames with calibration scores
    """
    method_columns = {name: column_name(name) for name in method_names}

    max_tpr = 99.9

    crossovers = {}
    for name in method_names:
        col = method_columns[name]
        am_scores = am_cal_feat[col].values
        na_scores = na_cal_feat[col].values
        crossovers[name] = _find_crossover_tpr(am_scores, na_scores)
        print(f"  {name}: crossover TPR = {crossovers[name]:.2f}%")

    for source_name, method_scores in source_scores.items():
        colors = list(ROOT.PlottingUtils.GetDefaultColors())

        canvas = ROOT.TCanvas(str(ROOT.PlottingUtils.GetRandomName()), "",
                              1200, 600)

        pad_plot = ROOT.TPad("pad_plot", "", 0.0, 0.0, 0.68, 1.0)
        pad_plot.SetRightMargin(0.03)
        pad_plot.Draw()

        pad_leg = ROOT.TPad("pad_leg", "", 0.72, 0.0, 1.0, 1.0)
        pad_leg.SetLeftMargin(0.0)
        pad_leg.SetRightMargin(0.05)
        pad_leg.SetBottomMargin(0)
        pad_leg.Draw()

        pad_plot.cd()

        graphs = []
        graph_names = []

        for i, name in enumerate(method_names):
            col = method_columns[name]
            alpha_cal = am_cal_feat[col].values
            gamma_cal = na_cal_feat[col].values

            crossover = crossovers[name]
            tpr_values = np.linspace(crossover, max_tpr, 500)

            strictness = (tpr_values - crossover) / (max_tpr - crossover) * 100

            eff = []
            for tpr in tpr_values:
                alpha_thresh, gamma_thresh = _determine_thresholds(
                    alpha_cal, gamma_cal, tpr)
                scores = method_scores[name]
                alpha_mask, gamma_mask = _classify_events(
                    scores, alpha_thresh, gamma_thresh)
                n_classified = np.sum(alpha_mask) + np.sum(gamma_mask)
                eff.append(n_classified / len(scores) * 100)

            x_arr = np.array(strictness, dtype=np.float64)
            y_arr = np.array(eff, dtype=np.float64)

            graph = ROOT.TGraph(len(x_arr), x_arr, y_arr)
            graph.SetLineColor(colors[i])
            graph.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
            graphs.append(graph)
            graph_names.append(name)

            if len(graphs) == 1:
                graph.SetTitle("")
                graph.GetXaxis().SetTitle("Threshold Strictness [%]")
                graph.GetYaxis().SetTitle("Efficiency [%]")
                graph.GetYaxis().SetTitleOffset(1)
                graph.GetXaxis().SetRangeUser(0, 100)
                graph.GetYaxis().SetRangeUser(30, 101)
                graph.Draw("AL")
            else:
                graph.Draw("L SAME")

        pad_leg.cd()
        leg = ROOT.PlottingUtils.AddLegend(0.0, 0.95, 0.2, 0.85)
        leg.SetMargin(0.15)
        for graph, name in zip(graphs, graph_names):
            leg.AddEntry(graph, name, "l")
        leg.Draw()

        safe_src = _clean_name(source_name)
        output_name = f"efficiency_vs_strictness_{safe_src}"
        ROOT.PlottingUtils.SaveFigure(canvas, output_name,
                                      ROOT.PlotSaveOptions.kLINEAR)
        canvas.Close()

        print(f"Efficiency vs strictness plot saved for {source_name}")


def _save_efficiency_table(records, method_names, source_names, tpr):
    """Save classification efficiency at maximum strictness to a LaTeX table."""
    eff_by = {(r["source"], r["method"]): r["efficiency"] for r in records}

    n_sources = len(source_names)
    col_spec = "l" + "c" * n_sources
    short_sources = [s.replace("&", "+") for s in source_names]

    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"  \centering")
    lines.append(r"  \caption{Classification efficiency [\%] at "
                 f"{tpr}\\% TPR."
                 r"}")
    lines.append(r"  \label{tab:efficiency}")
    lines.append(f"  \\begin{{tabular}}{{{col_spec}}}")
    lines.append(r"    \toprule")
    header = "    Method & " + " & ".join(short_sources) + r" \\"
    lines.append(header)
    lines.append(r"    \midrule")
    for name in method_names:
        vals = [f"{eff_by[(src, name)]:.1f}" for src in source_names]
        lines.append(f"    {name} & " + " & ".join(vals) + r" \\")
    lines.append(r"    \bottomrule")
    lines.append(r"  \end{tabular}")
    lines.append(r"\end{table}")

    table_str = "\n".join(lines)
    output_path = os.path.join(ANALYSIS_CACHE_DIR, "efficiency_table.txt")
    with open(output_path, "w") as f:
        f.write(table_str + "\n")
    print(f"\nLaTeX efficiency table saved to {output_path}")
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

    method_names = list(models.keys())

    alpha_feat_cache = os.path.join(ANALYSIS_CACHE_DIR,
                                    "test_alpha_features.pkl")
    gamma_feat_cache = os.path.join(ANALYSIS_CACHE_DIR,
                                    "test_gamma_features.pkl")
    if not (os.path.exists(alpha_feat_cache)
            and os.path.exists(gamma_feat_cache)):
        raise FileNotFoundError(
            "Analysis cache not found. Run analysis.py first.")

    print("Loading cached test features for threshold calibration...")
    am_cal_feat = pd.read_pickle(alpha_feat_cache)
    na_cal_feat = pd.read_pickle(gamma_feat_cache)

    am_lo_mask = ((am_cal_feat["light_output"] >= 375)
                  & (am_cal_feat["light_output"] <= 1575))
    na_lo_mask = ((na_cal_feat["light_output"] >= 0)
                  & (na_cal_feat["light_output"] <= 1750))
    am_cal_feat = am_cal_feat[am_lo_mask].reset_index(drop=True)
    na_cal_feat = na_cal_feat[na_lo_mask].reset_index(drop=True)

    print(f"Am-241 calibration events: {len(am_cal_feat)}")
    print(f"Na-22  calibration events: {len(na_cal_feat)}")

    print(f"Determining thresholds at {TRUE_POSITIVE_RATE}% TPR...")
    thresholds = {}
    for name in models:
        col = column_name(name)
        am_scores = am_cal_feat[col].values
        na_scores = na_cal_feat[col].values
        alpha_thresh, gamma_thresh = _determine_thresholds(
            am_scores, na_scores, TRUE_POSITIVE_RATE)
        thresholds[name] = (alpha_thresh, gamma_thresh)
        print(f"  {name}: alpha <= {alpha_thresh:.4f}, "
              f"gamma >= {gamma_thresh:.4f}")

    source_scores = {}
    source_features = {}

    for source_name, root_file in SOURCE_MAP.items():
        print(f"Loading: {source_name}")
        source_path = ROOT_FILES_DIR + root_file
        if not os.path.exists(source_path):
            print(f"  WARNING: File not found: {source_path}")
            continue

        features, waveforms = load_tree_data(source_path,
                                             scalar_branches=SCALAR_BRANCHES,
                                             array_branch="Samples",
                                             max_events=int(5e6))

        if len(features) == 0:
            print(f"  WARNING: No events found")
            continue

        print(f"  Events: {len(features)}")

        processed = process_waveforms(waveforms)

        scores = {}
        for name, model in models.items():
            scores[name] = model.predict(processed)
        source_scores[source_name] = scores
        source_features[source_name] = features

    _plot_efficiency_vs_strictness(method_names, source_scores, am_cal_feat,
                                  na_cal_feat)

    efficiency_records = []
    for source_name in source_scores:
        print(f"Classifying: {source_name}")
        features = source_features[source_name]
        for name in method_names:
            scores = source_scores[source_name][name]
            alpha_thresh, gamma_thresh = thresholds[name]
            alpha_mask, gamma_mask = _classify_events(scores, alpha_thresh,
                                                     gamma_thresh)
            n_classified = np.sum(alpha_mask) + np.sum(gamma_mask)
            efficiency = n_classified / len(scores) * 100
            print(f"  [{name}] efficiency: {efficiency:.1f}%")
            efficiency_records.append({
                "source": source_name,
                "method": name,
                "efficiency": efficiency,
            })
            _plot_classified_spectra(features, source_name, name, alpha_mask,
                                    gamma_mask, alpha_thresh, gamma_thresh)

    _save_efficiency_table(efficiency_records, method_names,
                          list(source_scores.keys()), TRUE_POSITIVE_RATE)

    print(f"\nComplete. Analyzed {len(source_scores)} sources "
          f"with {len(method_names)} methods.")
    print(f"Plots saved to: plots")


if __name__ == "__main__":
    main()
