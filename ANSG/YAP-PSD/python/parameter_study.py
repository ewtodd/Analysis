import os
import pickle
import numpy as np
import ROOT
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor
from analysis_utils.io import load_tree_data
from psd_utils import process_waveforms, bootstrap_auc
import analysis_utils

analysis_utils.load_cpp_library()
ROOT.gROOT.SetBatch(True)
ROOT.PlottingUtils.SetStylePreferences(ROOT.PlotSaveFormat.kPNG)

ROOT_FILES_DIR = "../macros/root_files/"
CACHE_DIR = "sweep_cache"

SCALAR_BRANCHES = [
    "pulse_height",
    "trigger_position",
    "long_integral",
    "light_output",
    "charge_comparison",
    "raw_shape_indicator",
    "clean_shape_indicator",
]

N_SEEDS = 5
SEEDS = [42, 123, 256, 789, 1024][:N_SEEDS]
DEFAULT_TRAIN_PER_CLASS = 5000

RF_CONFIG = dict(
    name="Random Forest",
    prefix="rf",
    model_class=RandomForestRegressor,
    color=ROOT.kBlue + 1,
    default_params=dict(
        n_estimators=100,
        max_depth=None,
        max_samples=None,
        max_features=1.0,
        n_jobs=-1,
    ),
    sweeps=[
        dict(
            sweep_name="n_estimators",
            values=[5, 10, 25, 50, 100, 150, 200, 300],
            x_title="Number of Trees",
            param_key="n_estimators",
        ),
        dict(
            sweep_name="n_training_samples",
            values=[500, 1000, 5000, 10000, 25000, 50000, 100000],
            x_title="Training Samples per Class",
            param_key=None,  # special handling
        ),
        dict(
            sweep_name="max_depth",
            values=[3, 5, 10, 20, 30, None],
            x_title="Max Depth (50 = None)",
            param_key="max_depth",
        ),
        dict(
            sweep_name="max_samples",
            values=[0.1, 0.2, 0.3, 0.4, 0.5, 0.632, 0.8, 1.0],
            x_title="Max Samples (Bootstrap Fraction)",
            param_key="max_samples",
        ),
    ],
)

XGB_CONFIG = dict(
    name="XGBoost",
    prefix="xgb",
    model_class=XGBRegressor,
    color=ROOT.kRed + 1,
    default_params=dict(
        n_estimators=100,
        max_depth=6,
        learning_rate=0.3,
        n_jobs=-1,
        verbosity=1,
    ),
    sweeps=[
        dict(
            sweep_name="n_estimators",
            values=[1, 5, 10, 25, 50, 100, 150, 200, 300, 500],
            x_title="Number of Boosting Rounds",
            param_key="n_estimators",
        ),
        dict(
            sweep_name="n_training_samples",
            values=[500, 1000, 5000, 10000, 25000, 50000, 100000],
            x_title="Training Samples per Class",
            param_key=None,
        ),
        dict(
            sweep_name="max_depth",
            values=[1, 2, 3, 5, 7, 10, 15],
            x_title="Max Depth",
            param_key="max_depth",
        ),
        dict(
            sweep_name="learning_rate",
            values=[0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0],
            x_title="Learning Rate",
            param_key="learning_rate",
        ),
    ],
)


def prepare_data(alpha_waveforms,
                 gamma_waveforms,
                 alpha_features,
                 gamma_features,
                 n_train_per_class=10000,
                 random_state=42):
    """Filter by light output, balance classes, split train/test, and process waveforms."""
    lo_lower = {"alpha": 375, "gamma": 0}
    lo_upper = {"alpha": 1575, "gamma": 1750}

    alpha_mask_idx = np.where(
        (alpha_features["light_output"] >= lo_lower["alpha"])
        & (alpha_features["light_output"] <= lo_upper["alpha"]))[0]
    gamma_mask_idx = np.where(
        (gamma_features["light_output"] >= lo_lower["gamma"])
        & (gamma_features["light_output"] <= lo_upper["gamma"]))[0]

    alpha_masked = alpha_waveforms[alpha_mask_idx]
    gamma_masked = gamma_waveforms[gamma_mask_idx]

    min_samples = min(n_train_per_class, len(alpha_masked), len(gamma_masked))
    print(f"Balanced training pool size per class: {min_samples}")

    rng = np.random.RandomState(random_state)
    alpha_train_sel = rng.choice(len(alpha_masked),
                                 size=min_samples,
                                 replace=False)
    gamma_train_sel = rng.choice(len(gamma_masked),
                                 size=min_samples,
                                 replace=False)

    train_alpha_wf = process_waveforms(alpha_masked[alpha_train_sel])
    train_gamma_wf = process_waveforms(gamma_masked[gamma_train_sel])

    x_train = np.vstack([train_alpha_wf, train_gamma_wf])
    y_train = np.array([0] * len(train_alpha_wf) + [1] * len(train_gamma_wf))

    alpha_train_set = set(alpha_train_sel.tolist())
    gamma_train_set = set(gamma_train_sel.tolist())

    alpha_test_sel = [
        i for i in range(len(alpha_masked)) if i not in alpha_train_set
    ]
    gamma_test_sel = [
        i for i in range(len(gamma_masked)) if i not in gamma_train_set
    ]

    test_alpha_wf = process_waveforms(alpha_masked[alpha_test_sel])
    test_gamma_wf = process_waveforms(gamma_masked[gamma_test_sel])

    x_test = np.vstack([test_alpha_wf, test_gamma_wf])
    y_test = np.array([0] * len(test_alpha_wf) + [1] * len(test_gamma_wf))

    print(
        f"Train: {len(x_train)} ({len(train_alpha_wf)} alpha + {len(train_gamma_wf)} gamma)"
    )
    print(
        f"Test:  {len(x_test)} ({len(test_alpha_wf)} alpha + {len(test_gamma_wf)} gamma)"
    )

    return x_train, y_train, x_test, y_test


def plot_sweep(x_values, y_values, y_errors, x_title, output_name, color):
    """Plot a single sweep with colored line + black markers/error bars."""
    canvas = ROOT.PlottingUtils.GetConfiguredCanvas()

    x_arr = np.array(x_values, dtype=np.float64)
    y_arr = np.array(y_values, dtype=np.float64)
    ex_arr = np.zeros(len(x_arr), dtype=np.float64)
    ey_arr = np.array(y_errors, dtype=np.float64)

    graph_line = ROOT.TGraph(len(x_arr), x_arr, y_arr)
    graph_line.SetLineColor(color)
    graph_line.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())

    graph_err = ROOT.TGraphErrors(len(x_arr), x_arr, y_arr, ex_arr, ey_arr)
    graph_err.SetLineColor(ROOT.kBlack)
    graph_err.SetMarkerColor(ROOT.kBlack)
    graph_err.SetMarkerStyle(20)
    graph_err.SetMarkerSize(1.2)

    graph_line.SetTitle("")
    if x_title == "Training Samples per Class":
        canvas.SetLogx(True)
    if "Max Samples" in x_title:
        graph_line.GetYaxis().SetTitleOffset(1.7)
        canvas.SetLeftMargin(0.2)
    if x_title == "Max Depth":
        graph_line.GetYaxis().SetTitleOffset(1.4)

    graph_line.GetXaxis().SetTitle(x_title)
    graph_line.GetYaxis().SetTitle("ROC AUC")

    y_min = min(y_values)
    y_max = max(y_values)
    y_range = y_max - y_min if y_max > y_min else 0.01
    graph_line.GetYaxis().SetRangeUser(y_min - 0.5 * y_range,
                                       y_max + 0.2 * y_range)

    graph_line.Draw("AL")
    graph_err.Draw("P SAME")

    ROOT.PlottingUtils.SaveFigure(canvas, output_name,
                                  ROOT.PlotSaveOptions.kLINEAR)
    canvas.Close()
    print(f"Saved {output_name}")


def run_sweep(config, sweep_name, values, param_key, x_train, y_train, x_test,
              y_test):
    """Run a parameter sweep, caching results to disk."""
    prefix = config["prefix"]
    model_class = config["model_class"]
    default_params = config["default_params"]
    n_per_class = len(x_train) // 2

    cache_file = os.path.join(CACHE_DIR, f"{prefix}_sweep_{sweep_name}.npz")
    if os.path.exists(cache_file):
        print(f"  Loading cached results: {cache_file}")
        data = np.load(cache_file)
        means = data["means"].tolist()
        stds = data["stds"].tolist()
        for val, m, s in zip(values, means, stds):
            print(f"  {sweep_name}={val}    AUC = {m:.4f} +/- {s:.4f}")
        return means, stds

    # For regular sweeps, subsample to DEFAULT_TRAIN_PER_CLASS so we don't
    # train with the entire (potentially huge) pool meant for the
    # n_training_samples sweep.
    if param_key is not None and n_per_class > DEFAULT_TRAIN_PER_CLASS:
        rng = np.random.RandomState(42)
        alpha_sel = rng.choice(n_per_class,
                               size=DEFAULT_TRAIN_PER_CLASS,
                               replace=False)
        gamma_sel = rng.choice(n_per_class,
                               size=DEFAULT_TRAIN_PER_CLASS,
                               replace=False) + n_per_class
        default_idx = np.concatenate([alpha_sel, gamma_sel])
        x_train_sweep = x_train[default_idx]
        y_train_sweep = y_train[default_idx]
        print(f"  Subsampled to {len(x_train_sweep)} training samples "
              f"({DEFAULT_TRAIN_PER_CLASS} per class)")
    else:
        x_train_sweep = x_train
        y_train_sweep = y_train
        print(f"  Using full training pool: {len(x_train_sweep)} samples")

    means, stds = [], []
    for val in values:
        print(f"  {sweep_name}={val}")
        params = dict(default_params, random_state=42)

        if param_key is not None:
            params[param_key] = val
            model = model_class(**params)
            model.fit(x_train_sweep, y_train_sweep)
        else:
            # n_training_samples: subsample from full training pool
            rng = np.random.RandomState(42)
            alpha_sel = rng.choice(n_per_class, size=val, replace=False)
            gamma_sel = rng.choice(n_per_class, size=val,
                                   replace=False) + n_per_class
            idx = np.concatenate([alpha_sel, gamma_sel])
            model = model_class(**params)
            model.fit(x_train[idx], y_train[idx])

        scores = model.predict(x_test)
        auc_mean, auc_std = bootstrap_auc(y_test, scores)
        means.append(auc_mean)
        stds.append(auc_std)
        print(f"    AUC = {auc_mean:.4f} +/- {auc_std:.4f}")

    np.savez(cache_file, means=np.array(means), stds=np.array(stds))
    return means, stds


def plot_feature_importance(model_class, default_params, prefix, color, name,
                            x_train, y_train, avg_waveform):
    """Train or load N_SEEDS models and plot averaged feature importance."""
    # Subsample to DEFAULT_TRAIN_PER_CLASS if pool is larger
    n_per_class = len(x_train) // 2
    if n_per_class > DEFAULT_TRAIN_PER_CLASS:
        rng = np.random.RandomState(0)
        alpha_sel = rng.choice(n_per_class,
                               size=DEFAULT_TRAIN_PER_CLASS,
                               replace=False)
        gamma_sel = rng.choice(n_per_class,
                               size=DEFAULT_TRAIN_PER_CLASS,
                               replace=False) + n_per_class
        idx = np.concatenate([alpha_sel, gamma_sel])
        x_train = x_train[idx]
        y_train = y_train[idx]

    all_importances = []
    for seed in SEEDS:
        model_file = os.path.join(CACHE_DIR,
                                  f"{prefix}_importance_seed{seed}.pkl")
        if os.path.exists(model_file):
            print(f"  Loading cached model: {model_file}")
            with open(model_file, "rb") as fh:
                model = pickle.load(fh)
        else:
            print(f"  Training seed {seed}...")
            params = dict(default_params, random_state=seed)
            model = model_class(**params)
            model.fit(x_train, y_train)
            with open(model_file, "wb") as fh:
                pickle.dump(model, fh)
        all_importances.append(model.feature_importances_)

    importances = np.array(all_importances)
    mean_imp = np.mean(importances, axis=0)
    std_imp = np.std(importances, axis=0)

    wf_max = np.max(avg_waveform)
    avg_wf_norm = avg_waveform / wf_max if wf_max > 0 else avg_waveform

    imp_max = np.max(mean_imp)
    mean_imp_norm = mean_imp / imp_max if imp_max > 0 else mean_imp
    std_imp_norm = std_imp / imp_max if imp_max > 0 else std_imp

    n_points = len(avg_waveform)
    x_values = np.arange(n_points, dtype=np.float64) * 2

    canvas = ROOT.PlottingUtils.GetConfiguredCanvas()

    graph_waveform = ROOT.TGraph(n_points, x_values,
                                 avg_wf_norm.astype(np.float64))
    graph_waveform.SetLineColor(ROOT.kGray + 2)
    graph_waveform.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
    graph_waveform.SetTitle("")
    graph_waveform.GetXaxis().SetTitle("Time [ns]")
    graph_waveform.GetYaxis().SetTitle("Normalized Amplitude [a.u.]")
    graph_waveform.GetXaxis().SetRangeUser(0, x_values[-1])
    graph_waveform.GetYaxis().SetRangeUser(-0.1, 1.1)
    graph_waveform.Draw("AL")

    ex = np.zeros(n_points, dtype=np.float64)
    graph_importance = ROOT.TGraphErrors(n_points, x_values,
                                         mean_imp_norm.astype(np.float64), ex,
                                         std_imp_norm.astype(np.float64))
    graph_importance.SetLineColor(color)
    graph_importance.SetLineWidth(ROOT.PlottingUtils.GetLineWidth())
    graph_importance.SetFillColorAlpha(color, 0.2)
    graph_importance.Draw("L3 SAME")

    leg = ROOT.PlottingUtils.AddLegend(0.42, 0.88, 0.6, 0.85)
    leg.AddEntry(graph_waveform, "Average #alpha Waveform", "l")
    leg.AddEntry(graph_importance, f"{name} Feature Importance", "lf")
    leg.SetMargin(0.1)
    leg.Draw()

    output_file = f"feature_importance_{prefix}"
    ROOT.PlottingUtils.SaveFigure(canvas, output_file,
                                  ROOT.PlotSaveOptions.kLINEAR)
    canvas.Close()
    print(f"Saved {output_file}")


def recommend_value(values, means, default_value=None):
    """Find the recommended parameter value from a sweep.

    If the AUC curve has a clear peak, return the peak value.
    Otherwise, find the plateau onset: the first value where AUC >= 0.99.

    Then compare against the library default: if the default gives higher AUC
    than the plateau onset, prefer the default.
    """
    best_idx = int(np.argmax(means))
    best_auc = means[best_idx]

    # Check if it's a peak (AUC drops on both sides) vs a plateau
    is_peak = (best_idx > 0 and best_idx < len(means) - 1)

    if is_peak:
        rec_val, rec_auc = values[best_idx], best_auc
    else:
        # Plateau: first value that crosses 0.99 AUC
        rec_val, rec_auc = values[best_idx], best_auc
        for i, m in enumerate(means):
            if m >= 0.99:
                rec_val, rec_auc = values[i], m
                break

    # Compare against library default if available
    if default_value is not None and default_value in values:
        default_idx = values.index(default_value)
        default_auc = means[default_idx]
        if default_auc > rec_auc:
            print(
                f"  >> Plateau onset at {rec_val} (AUC={rec_auc:.4f}), "
                f"but default {default_value} is better (AUC={default_auc:.4f})"
            )
            rec_val, rec_auc = default_value, default_auc

    return rec_val, rec_auc


def run_all_sweeps(config, x_train, y_train, x_test, y_test):
    """Run all sweeps for a given model configuration.

    Returns a dict of recommended parameter values (param_key -> value),
    excluding n_training_samples.
    """
    color = config["color"]
    n_per_class = len(x_train) // 2
    recommended = {}

    for sweep in config["sweeps"]:
        sweep_name = sweep["sweep_name"]
        values = sweep["values"]
        x_title = sweep["x_title"]
        param_key = sweep["param_key"]
        output_name = f"{config['prefix']}_auc_vs_{sweep_name}"

        if sweep_name == "n_training_samples":
            values = [v for v in values if v <= n_per_class]

        print(f"{config['name']}: AUC vs {sweep_name}")

        means, stds = run_sweep(config, sweep_name, values, param_key, x_train,
                                y_train, x_test, y_test)

        default_value = config["default_params"].get(
            param_key) if param_key else None
        rec_val, rec_auc = recommend_value(values, means, default_value)
        print(
            f"  >> Recommended {sweep_name} = {rec_val} (AUC = {rec_auc:.4f})")

        if param_key is not None:
            recommended[param_key] = rec_val

        plot_x = values
        if sweep_name == "max_depth" and None in values:
            plot_x = [d if d is not None else 50 for d in values]

        plot_sweep(plot_x, means, stds, x_title, output_name, color)

    return recommended


def main():
    os.makedirs("plots", exist_ok=True)
    os.makedirs(CACHE_DIR, exist_ok=True)

    print("Loading alpha data (Am-241)...")
    alpha_features, alpha_waveforms = load_tree_data(
        ROOT_FILES_DIR + "Am241.root",
        scalar_branches=SCALAR_BRANCHES,
        array_branch="Samples",
    )
    print(
        f"Alpha events: {len(alpha_features)}, waveform shape: {alpha_waveforms.shape}"
    )

    print("Loading gamma data (Na-22)...")
    gamma_features, gamma_waveforms = load_tree_data(
        ROOT_FILES_DIR + "Na22.root",
        scalar_branches=SCALAR_BRANCHES,
        array_branch="Samples",
    )
    print(
        f"Gamma events: {len(gamma_features)}, waveform shape: {gamma_waveforms.shape}"
    )

    # Find the max n_training_samples across all configs so prepare_data
    # allocates enough training data for the full sweep.
    max_train = DEFAULT_TRAIN_PER_CLASS
    for config in [RF_CONFIG, XGB_CONFIG]:
        for sweep in config["sweeps"]:
            if sweep["sweep_name"] == "n_training_samples":
                max_train = max(max_train, max(sweep["values"]))

    print("Preparing train/test data")
    x_train, y_train, x_test, y_test = prepare_data(
        alpha_waveforms,
        gamma_waveforms,
        alpha_features,
        gamma_features,
        n_train_per_class=max_train)

    f = ROOT.TFile.Open(ROOT_FILES_DIR + "Am241.root")
    avg_wf_graph = f.Get("average_waveform")
    avg_waveform = np.array(
        [avg_wf_graph.GetPointY(i) for i in range(avg_wf_graph.GetN())])
    f.Close()

    for config in [RF_CONFIG, XGB_CONFIG]:
        print(f"  {config['name']} Hyperparameter Study")

        recommended = run_all_sweeps(config, x_train, y_train, x_test, y_test)

        optimized_params = dict(config["default_params"])
        optimized_params.update(recommended)
        print(f"{config['name']}: Optimized params: {recommended}")

        print(f"{config['name']}: Feature importance (averaged over seeds)")
        plot_feature_importance(config["model_class"], optimized_params,
                                config["prefix"], config["color"],
                                config["name"], x_train, y_train, avg_waveform)

    print("Done. All plots saved.")


if __name__ == "__main__":
    main()
