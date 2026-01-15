#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TTree.h>

void Stack() {
  InitUtils::SetROOTPreferences();

  TString suffix = Constants::FILTERED ? "_filtered" : "";

  TString CdShieldSignal = "01132026-CdShield-GeSamplesIn-10Percent" + suffix;
  TString CuShieldSignal_01132026 =
      "01132026-CuShield-GeSamplesIn-10Percent" + suffix;
  TString CuShieldSignal_01142026 =
      "01142026-CuShield-GeSamplesIn-10Percent" + suffix;

  Bool_t isFiltered = CdShieldSignal.Contains("_filtered");
  if (isFiltered && !Constants::FILTERED) {
    std::cerr << "ERROR: File names contain '_filtered' but "
                 "Constants::FILTERED is not set!"
              << std::endl;
    return;
  }
  if (!isFiltered && Constants::FILTERED) {
    std::cerr << "ERROR: Constants::FILTERED is set but file names don't "
                 "contain '_filtered'!"
              << std::endl;
    return;
  }
  if ((isFiltered || Constants::FILTERED) && Constants::NORMALIZE_BY_TIME) {
    std::cerr << "Cannot normalize by time on filtered data." << std::endl;
  }

  TFile *file_CdShieldSignal =
      new TFile("root_files/" + CdShieldSignal + ".root", "READ");
  if (!file_CdShieldSignal || file_CdShieldSignal->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << CdShieldSignal << ".root"
              << std::endl;
    return;
  }
  TH1F *zoomedHist_CdShieldSignal =
      static_cast<TH1F *>(file_CdShieldSignal->Get("zoomedHist"));

  TFile *file_CuShieldSignal_01132026 =
      new TFile("root_files/" + CuShieldSignal_01132026 + ".root", "READ");
  if (!file_CuShieldSignal_01132026 ||
      file_CuShieldSignal_01132026->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << CuShieldSignal_01132026 << ".root"
              << std::endl;
    return;
  }
  TH1F *zoomedHist_CuShieldSignal_01132026 =
      static_cast<TH1F *>(file_CuShieldSignal_01132026->Get("zoomedHist"));

  TFile *file_CuShieldSignal_01142026 =
      new TFile("root_files/" + CuShieldSignal_01142026 + ".root", "READ");
  if (!file_CuShieldSignal_01142026 ||
      file_CuShieldSignal_01142026->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << CuShieldSignal_01142026 << ".root"
              << std::endl;
    return;
  }
  TH1F *zoomedHist_CuShieldSignal_01142026 =
      static_cast<TH1F *>(file_CuShieldSignal_01142026->Get("zoomedHist"));

  TH1F *norm_CdShieldSignal = static_cast<TH1F *>(
      zoomedHist_CdShieldSignal->Clone("norm_CdShieldSignal"));
  TH1F *norm_CuShieldSignal_01132026 =
      static_cast<TH1F *>(zoomedHist_CuShieldSignal_01132026->Clone(
          "norm_CuShieldSignal_01132026"));
  TH1F *norm_CuShieldSignal_01142026 =
      static_cast<TH1F *>(zoomedHist_CuShieldSignal_01142026->Clone(
          "norm_CuShieldSignal_01142026"));

  if (norm_CdShieldSignal->Integral() > 0)
    norm_CdShieldSignal->Scale(1.0 / norm_CdShieldSignal->Integral());
  if (norm_CuShieldSignal_01132026->Integral() > 0)
    norm_CuShieldSignal_01132026->Scale(
        1.0 / norm_CuShieldSignal_01132026->Integral());
  if (norm_CuShieldSignal_01142026->Integral() > 0)
    norm_CuShieldSignal_01142026->Scale(
        1.0 / norm_CuShieldSignal_01142026->Integral());

  std::vector<Int_t> colors = PlottingUtils::GetDefaultColors();

  PlottingUtils::ConfigureHistogram(norm_CdShieldSignal, colors[0]);
  PlottingUtils::ConfigureHistogram(norm_CuShieldSignal_01132026, colors[1]);
  PlottingUtils::ConfigureHistogram(norm_CuShieldSignal_01142026, colors[2]);

  norm_CdShieldSignal->SetTitle(Form(
      "; Energy [keV]; Normalized Counts / %d eV", Constants::BIN_WIDTH_EV));

  norm_CdShieldSignal->GetYaxis()->SetRangeUser(0, 0.01);

  TCanvas *canvas = new TCanvas("", "", 1200, 800);
  PlottingUtils::ConfigureCanvas(canvas);

  norm_CdShieldSignal->Draw("HIST");
  norm_CuShieldSignal_01132026->Draw("HIST SAME");
  norm_CuShieldSignal_01142026->Draw("HIST SAME");

  TLegend *legend = new TLegend(0.2, 0.55, 0.4, 0.88);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(norm_CdShieldSignal, "Cd Shield Signal", "l");
  legend->AddEntry(norm_CuShieldSignal_01132026, "Cu Shield Signal (01/13)",
                   "l");
  legend->AddEntry(norm_CuShieldSignal_01142026, "Cu Shield Signal (01/14)",
                   "l");
  legend->Draw();

  PlottingUtils::SaveFigure(canvas, "all_histograms_normalized_zoomed.png",
                            kFALSE);

  std::cout << "All histograms plotted together (normalized by integral)!"
            << std::endl;

  TString outputNameSignal = "Combined_Signal" + suffix;
  TFile *outputFileSignal =
      new TFile("root_files/" + outputNameSignal + ".root", "UPDATE");
  canvas->Write("stackedHist", TObject::kOverwrite);
  outputFileSignal->Close();

  file_CdShieldSignal->Close();
  file_CuShieldSignal_01132026->Close();
  file_CuShieldSignal_01142026->Close();
}
