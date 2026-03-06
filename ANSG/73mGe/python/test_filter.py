from analysis_utils import load_cpp_library
from analysis_utils.init import set_root_preferences
from analysis_utils.io import load_tree_data
import numpy as np
import constants as C

ROOT = load_cpp_library()


def run_filter(filenames):
    for filename in filenames:
        filepath = f"root_files/{filename}.root"

        df = load_tree_data(
            f"{C.ROOT_FILES_DIR}/{filename}.root", tree_name="bef_tree"
        )

        energy = df["energykeV"].values
        x = df["xum"].values
        y = df["yum"].values
        z = df["zum"].values

        xy = ROOT.TH2F(
            ROOT.PlottingUtils.GetRandomName(),
            "; X Position [um]; Y Position [um]",
            22,
            -215,
            215,
            22,
            -215,
            215,
        )
        ez = ROOT.TH2F(
            ROOT.PlottingUtils.GetRandomName(),
            "; Energy [keV]; Interaction Z Position [um]",
            C.ZOOMED_NBINS,
            C.ZOOMED_XMIN,
            C.ZOOMED_XMAX,
            200,
            0,
            100,
        )

        included = ROOT.TH1F(
            ROOT.PlottingUtils.GetRandomName(),
            f"{filename}; Energy [keV]; Counts / {C.BIN_WIDTH_EV} eV",
            C.ZOOMED_NBINS,
            C.ZOOMED_XMIN,
            C.ZOOMED_XMAX,
        )
        excluded = ROOT.TH1F(
            ROOT.PlottingUtils.GetRandomName(),
            f"{filename}; Energy [keV]; Counts / {C.BIN_WIDTH_EV} eV",
            C.ZOOMED_NBINS,
            C.ZOOMED_XMIN,
            C.ZOOMED_XMAX,
        )

        # Build exclusion mask
        is_excluded = z < C.FILTER_DEPTH_UM
        if not C.FILTERED:
            is_excluded |= df["nInteractions"].values != 1

        zoomed = (energy > C.ZOOMED_XMIN) & (energy < C.ZOOMED_XMAX)

        for xmin, xmax, ymin, ymax in C.FILTER_REGIONS_EXCLUDE_XY_UM:
            is_excluded |= (
                zoomed & (x >= xmin) & (x <= xmax)
                & (y >= ymin) & (y <= ymax)
            )

        # Fill 2D histograms (zoomed region only)
        zx, zy, ze, zz = x[zoomed], y[zoomed], energy[zoomed], z[zoomed]
        for xi, yi in zip(zx, zy):
            xy.Fill(xi, yi)
        for ei, zi in zip(ze, zz):
            ez.Fill(ei, zi)

        # Fill 1D spectra
        inc_mask = zoomed & ~is_excluded
        exc_mask = zoomed & is_excluded
        np.vectorize(included.Fill)(energy[inc_mask])
        np.vectorize(excluded.Fill)(energy[exc_mask])

        print(f"Created histograms for {filename}")

        # XY position plot with exclusion boxes
        canvas_xy = ROOT.TCanvas("", "", 1200, 800)
        canvas_xy.SetRightMargin(0.2)
        xy.Draw("COLZ")
        ROOT.PlottingUtils.Configure2DHistogram(xy, canvas_xy)
        canvas_xy.SetLogz(False)
        canvas_xy.Update()

        for xmin, xmax, ymin, ymax in C.FILTER_REGIONS_EXCLUDE_XY_UM:
            box = ROOT.TBox(xmin, ymin, xmax, ymax)
            box.SetLineColor(ROOT.kBlack)
            box.SetLineWidth(3)
            box.SetLineStyle(2)
            box.SetFillColorAlpha(ROOT.kBlack, 0.5)
            box.Draw()

        ROOT.PlottingUtils.SaveFigure(canvas_xy, f"{filename}_YvsX",
                                      ROOT.PlotSaveOptions.kLINEAR)

        # Energy vs Z depth plot with depth threshold box
        canvas_ez = ROOT.TCanvas("", "", 1200, 800)
        canvas_ez.SetRightMargin(0.2)
        ez.Draw("COLZ")
        ROOT.PlottingUtils.Configure2DHistogram(ez, canvas_ez)
        canvas_ez.SetLogz(False)
        canvas_ez.Update()

        box = ROOT.TBox(C.ZOOMED_XMIN, 0, C.ZOOMED_XMAX, C.FILTER_DEPTH_UM)
        box.SetLineColor(ROOT.kBlack)
        box.SetLineWidth(3)
        box.SetLineStyle(2)
        box.SetFillColorAlpha(ROOT.kBlack, 0.5)
        box.Draw()

        ROOT.PlottingUtils.SaveFigure(canvas_ez, f"{filename}_ZvsE",
                                      ROOT.PlotSaveOptions.kLINEAR)

        # Accepted vs rejected spectra overlay
        canvas_spectra = ROOT.PlottingUtils.GetConfiguredCanvas(False)

        max_y = max(included.GetMaximum(), excluded.GetMaximum())

        ROOT.PlottingUtils.ConfigureHistogram(included, ROOT.kRed)
        included.GetYaxis().SetRangeUser(0, max_y * 1.1)
        included.SetFillStyle(0)
        included.SetLineWidth(2)
        included.GetYaxis().SetTitleOffset(1.5)
        included.Draw("HIST")

        ROOT.PlottingUtils.ConfigureHistogram(excluded, ROOT.kBlack)
        excluded.SetFillStyle(0)
        excluded.SetLineWidth(2)
        excluded.Draw("HIST SAME")

        if filename == C.CDSHIELDSIGNAL_10PERCENT_20260113:
            line68 = ROOT.TLine(68.752, 0, 68.752, max_y * 1.1)
            line68.SetLineColor(ROOT.kViolet)
            line68.SetLineWidth(2)
            line68.SetLineStyle(2)
            line68.Draw()

            label68 = ROOT.PlottingUtils.AddText(
                "^{73m}Ge 7/2^{+} #rightarrow 9/2^{+} #gamma ray", 0.47, 0.75)
            label68.SetTextColor(ROOT.kViolet)

        leg = ROOT.PlottingUtils.AddLegend(0.72, 0.9, 0.75, 0.88)
        leg.AddEntry(included, "Accepted", "l")
        leg.AddEntry(excluded, "Rejected", "l")
        leg.Draw()

        canvas_spectra.SetLeftMargin(0.2)
        ROOT.PlottingUtils.SaveFigure(canvas_spectra,
                                      f"{filename}_FilterSpectra",
                                      ROOT.PlotSaveOptions.kLINEAR)

        print(f"Wrote plots for {filename}")

        f = ROOT.TFile(filepath, "UPDATE")
        xy.Write("XY", ROOT.TObject.kOverwrite)
        ez.Write("EZ", ROOT.TObject.kOverwrite)
        canvas_spectra.cd()
        canvas_spectra.Write("AcceptedRejectedSpectra",
                             ROOT.TObject.kOverwrite)
        f.Close()


def main():
    set_root_preferences()

    filenames = [
        C.CDSHIELDSIGNAL_10PERCENT_20260113,
        C.POSTREACTOR_AM241_BA133_20260116,
    ]

    run_filter(filenames)


if __name__ == "__main__":
    main()
