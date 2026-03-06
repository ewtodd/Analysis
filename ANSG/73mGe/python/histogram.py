from analysis_utils import load_cpp_library
from analysis_utils.init import set_root_preferences
from analysis_utils.io import load_tree_data
import numpy as np
import constants as C

ROOT = load_cpp_library()


def add_histogram(filename):
    filepath = f"{C.ROOT_FILES_DIR}/{filename}.root"

    tree_name = "bef_tree_event_summary"
    energy_branch = "totalEnergykeV"

    df = load_tree_data(filepath, tree_name=tree_name)
    energy = df[energy_branch].values

    hist = ROOT.TH1F(
        str(ROOT.PlottingUtils.GetRandomName()),
        f"{filename}; Energy [keV]; Counts / {C.BIN_WIDTH_EV} eV",
        C.HIST_NBINS,
        C.HIST_XMIN,
        C.HIST_XMAX,
    )
    zoomed_hist = ROOT.TH1F(
        str(ROOT.PlottingUtils.GetRandomName()),
        f"{filename}; Energy [keV]; Counts / {C.BIN_WIDTH_EV} eV",
        C.ZOOMED_NBINS,
        C.ZOOMED_XMIN,
        C.ZOOMED_XMAX,
    )
    peak_hist = ROOT.TH1F(
        str(ROOT.PlottingUtils.GetRandomName()),
        f"{filename}; Energy [keV]; Counts / {C.BIN_WIDTH_EV} eV",
        C.PEAK_NBINS,
        C.PEAK_XMIN,
        C.PEAK_XMAX,
    )

    np.vectorize(hist.Fill)(energy)
    zoomed_mask = (energy > C.ZOOMED_XMIN) & (energy < C.ZOOMED_XMAX)
    np.vectorize(zoomed_hist.Fill)(energy[zoomed_mask])
    peak_mask = (energy > C.PEAK_XMIN) & (energy < C.PEAK_XMAX)
    np.vectorize(peak_hist.Fill)(energy[peak_mask])

    print(f"Created histograms for {filename} (branch: {energy_branch})")

    ROOT.PlottingUtils.ConfigureHistogram(hist, ROOT.kP10Violet)
    zoomedCanvas = ROOT.PlottingUtils.GetConfiguredCanvas()
    ROOT.PlottingUtils.ConfigureAndDrawHistogram(zoomed_hist, ROOT.kP10Violet)
    ROOT.PlottingUtils.SaveFigure(zoomedCanvas, f"zoomed_hist_{filename}",
                                  ROOT.PlotSaveOptions.kLOG)
    ROOT.PlottingUtils.ConfigureHistogram(peak_hist, ROOT.kP10Violet)

    hist.GetYaxis().SetTitleOffset(1.2)
    zoomed_hist.GetYaxis().SetTitleOffset(1.2)
    peak_hist.GetYaxis().SetTitleOffset(1.2)

    f = ROOT.TFile(filepath, "UPDATE")
    hist.Write("hist", ROOT.TObject.kOverwrite)
    zoomed_hist.Write("zoomedHist", ROOT.TObject.kOverwrite)
    peak_hist.Write("peakHist", ROOT.TObject.kOverwrite)
    f.Close()

    print(f"Wrote histograms for {filename}")


def main():
    set_root_preferences()
    for dataset in C.ALL_DATASETS:
        add_histogram(dataset)


if __name__ == "__main__":
    main()
