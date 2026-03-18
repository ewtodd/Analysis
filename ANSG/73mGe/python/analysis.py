from analysis_utils import load_cpp_library
from analysis_utils.init import set_root_preferences
from analysis_utils.fitting import fit_single_peak, plot_single_peak_fit, fit_double_peak, plot_double_peak_fit
import numpy as np
import constants as C
from data import load_filtered
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = load_cpp_library()

E_AM241 = 59.5409
E_BA133_53 = 53.16
E_BA133_81 = 80.98
E_PB_KA1 = 72.8042
E_PB_KA2 = 74.9694
E_CD114M = 95.9023


def energy_hist(energy, xmin, xmax, title=""):
    nbins = int((xmax - xmin) / C.BIN_WIDTH_KEV)
    hist = ROOT.TH1F(
        "",
        f"{title}; Energy [keV]; Counts / {C.BIN_WIDTH_EV} eV",
        nbins,
        xmin,
        xmax,
    )
    for e in energy:
        if (e > xmin and e < xmax):
            hist.Fill(e)
    return hist


def test_single_peak_fit(filename):
    df = load_filtered(filename)
    energy = df["energykeV"]

    fit_low = 52
    fit_up = 68

    m = fit_single_peak(energy,
                        fit_low,
                        fit_up,
                        E_AM241,
                        cache_path="fit_cache/am241_59keV.pkl")
    hist = energy_hist(energy, C.ZOOMED_XMIN, C.ZOOMED_XMAX)
    plot_single_peak_fit(hist, m, fit_low, fit_up, "test", "am241")


def test_double_peak_fit(filename):
    df = load_filtered(filename)
    energy = df["energykeV"]

    fit_low = 66
    fit_up = 81

    m = fit_double_peak(energy,
                        fit_low,
                        fit_up,
                        E_PB_KA1,
                        E_PB_KA2,
                        cache_path="fit_cache/pb.pkl")
    hist = energy_hist(energy, C.ZOOMED_XMIN, C.ZOOMED_XMAX)
    plot_double_peak_fit(hist, m, fit_low, fit_up, "test", "cdshield_bkg")


def main():
    set_root_preferences()
    test_single_peak_fit(C.POSTREACTOR_AM241_20260115)
    test_single_peak_fit(C.POSTREACTOR_AM241_20260113)
    test_double_peak_fit(C.CDSHIELDBACKGROUND_10PERCENT_20260113)


if __name__ == "__main__":
    main()
