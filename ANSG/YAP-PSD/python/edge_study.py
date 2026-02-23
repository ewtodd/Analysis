import numpy as np
import os
from analysis_utils.io import load_tree_data
from psd_utils import regress_waveforms, process_waveforms
from regressors import get_default_regressors
import analysis_utils

analysis_utils.load_cpp_library()
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.PlottingUtils.SetStylePreferences(ROOT.PlotSaveFormat.kPNG)

ROOT_FILES_DIR = "../macros/root_files/"
EDGE_STUDY_CACHE_DIR = "edge_study_cache"

SCALAR_BRANCHES = [
    "pulse_height",
    "trigger_position",
    "long_integral",
    "light_output",
    "charge_comparison",
    "raw_shape_indicator",
    "clean_shape_indicator",
]

SWEEP_CACHE_DIR = "sweep_cache"


def load_peak_sample():
    """Find the peak sample index from the average alpha waveform."""
    wf_cache = os.path.join(SWEEP_CACHE_DIR, "avg_waveform.npy")
    if os.path.exists(wf_cache):
        avg_waveform = np.load(wf_cache)
    else:
        f = ROOT.TFile.Open(ROOT_FILES_DIR + "Am241.root")
        avg_wf_graph = f.Get("average_waveform")
        avg_waveform = np.array(
            [avg_wf_graph.GetPointY(i) for i in range(avg_wf_graph.GetN())])
        f.Close()
    peak = int(np.argmax(avg_waveform))
    print(f"Peak sample from average waveform: {peak}")
    return peak


def make_edge_processors(peak_sample):
    """Return rising/falling edge processing functions for a given peak sample."""

    def process_rising_edge(waveforms):
        """Normalize waveforms and keep only the rising edge (up to peak)."""
        normalized = process_waveforms(waveforms)
        return normalized[:, :peak_sample + 1]

    def process_falling_edge(waveforms):
        """Normalize waveforms and keep only the falling edge (from peak)."""
        normalized = process_waveforms(waveforms)
        return normalized[:, peak_sample:]

    return process_rising_edge, process_falling_edge


def get_regressors_for_variant(variant_name):
    """Return regressor configs with cache paths specific to a variant."""
    regressors = get_default_regressors()
    cache_dir = os.path.join(EDGE_STUDY_CACHE_DIR, variant_name)
    os.makedirs(cache_dir, exist_ok=True)
    for r in regressors:
        r["file"] = os.path.join(cache_dir, os.path.basename(r["file"]))
    return regressors


def main():
    os.makedirs("plots", exist_ok=True)

    peak_sample = load_peak_sample()
    process_rising_edge, process_falling_edge = make_edge_processors(
        peak_sample)

    variants = {
        "rising": process_rising_edge,
        "falling": process_falling_edge,
    }

    all_cached = all(
        os.path.isdir(os.path.join(EDGE_STUDY_CACHE_DIR, v)) and all(
            os.path.exists(os.path.join(EDGE_STUDY_CACHE_DIR, v, f)) for f in [
                "test_alpha_features.pkl", "test_gamma_features.pkl",
                "test_waveforms.npz", "regressor_names.pkl"
            ]) for v in variants)

    if all_cached:
        print("All edge study caches found — skipping ROOT file loading.")
        alpha_waveforms = gamma_waveforms = None
        alpha_features = gamma_features = None
    else:
        print("Loading alpha data (Am-241)...")
        alpha_features, alpha_waveforms = load_tree_data(
            ROOT_FILES_DIR + "Am241.root",
            scalar_branches=SCALAR_BRANCHES,
            array_branch="Samples",
        )
        print(f"Alpha events: {len(alpha_features)}, "
              f"waveform shape: {alpha_waveforms.shape}")

        print("Loading gamma data (Na-22)...")
        gamma_features, gamma_waveforms = load_tree_data(
            ROOT_FILES_DIR + "Na22.root",
            scalar_branches=SCALAR_BRANCHES,
            array_branch="Samples",
        )
        print(f"Gamma events: {len(gamma_features)}, "
              f"waveform shape: {gamma_waveforms.shape}")

    results = {}
    for name, process_func in variants.items():
        print(f"  Running variant: {name}")

        cache_dir = os.path.join(EDGE_STUDY_CACHE_DIR, name)
        regressors = get_regressors_for_variant(name)

        results[name] = regress_waveforms(
            (alpha_waveforms, gamma_waveforms),
            (alpha_features, gamma_features),
            regressors=regressors,
            process_func=process_func,
            cache_dir=cache_dir,
            plot_prefix=f"{name}_",
        )

    print("Edge study complete")


if __name__ == "__main__":
    main()
