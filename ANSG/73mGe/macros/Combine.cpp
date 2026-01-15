#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TTree.h>

void Combine() {
  InitUtils::SetROOTPreferences();

  TString suffix = Constants::FILTERED ? "_filtered" : "";

  TString CdShieldSignal = "01132026-CdShield-GeSamplesIn-10Percent" + suffix;
  TString CdShieldBackground =
      "01132026-CdShield-ActiveBackground-10Percent" + suffix;
  TString CuShieldSignal_01132026 =
      "01132026-CuShield-GeSamplesIn-10Percent" + suffix;
  TString CuShieldBackground_01132026 =
      "01132026-CuShield-ActiveBackground-Am241-10Percent" + suffix;
  TString CuShieldSignal_01142026 =
      "01142026-CuShield-GeSamplesIn-10Percent" + suffix;
  TString CuShieldBackground_01142026 =
      "01142026-CuShield-ActiveBackground-10Percent" + suffix;

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
  TH1F *hist_CdShieldSignal =
      static_cast<TH1F *>(file_CdShieldSignal->Get("hist"));
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
  TH1F *hist_CuShieldSignal_01132026 =
      static_cast<TH1F *>(file_CuShieldSignal_01132026->Get("hist"));
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
  TH1F *hist_CuShieldSignal_01142026 =
      static_cast<TH1F *>(file_CuShieldSignal_01142026->Get("hist"));
  TH1F *zoomedHist_CuShieldSignal_01142026 =
      static_cast<TH1F *>(file_CuShieldSignal_01142026->Get("zoomedHist"));

  TFile *file_CdShieldBackground =
      new TFile("root_files/" + CdShieldBackground + ".root", "READ");
  if (!file_CdShieldBackground || file_CdShieldBackground->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << CdShieldBackground << ".root"
              << std::endl;
    return;
  }
  TH1F *hist_CdShieldBackground =
      static_cast<TH1F *>(file_CdShieldBackground->Get("hist"));
  TH1F *zoomedHist_CdShieldBackground =
      static_cast<TH1F *>(file_CdShieldBackground->Get("zoomedHist"));

  TFile *file_CuShieldBackground_01132026 =
      new TFile("root_files/" + CuShieldBackground_01132026 + ".root", "READ");
  if (!file_CuShieldBackground_01132026 ||
      file_CuShieldBackground_01132026->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << CuShieldBackground_01132026 << ".root"
              << std::endl;
    return;
  }
  TH1F *hist_CuShieldBackground_01132026 =
      static_cast<TH1F *>(file_CuShieldBackground_01132026->Get("hist"));
  TH1F *zoomedHist_CuShieldBackground_01132026 =
      static_cast<TH1F *>(file_CuShieldBackground_01132026->Get("zoomedHist"));

  TFile *file_CuShieldBackground_01142026 =
      new TFile("root_files/" + CuShieldBackground_01142026 + ".root", "READ");
  if (!file_CuShieldBackground_01142026 ||
      file_CuShieldBackground_01142026->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << CuShieldBackground_01142026 << ".root"
              << std::endl;
    return;
  }
  TH1F *hist_CuShieldBackground_01142026 =
      static_cast<TH1F *>(file_CuShieldBackground_01142026->Get("hist"));
  TH1F *zoomedHist_CuShieldBackground_01142026 =
      static_cast<TH1F *>(file_CuShieldBackground_01142026->Get("zoomedHist"));

  TString perSecond = Constants::NORMALIZE_BY_TIME ? " / s" : "";

  if (Constants::NORMALIZE_BY_TIME) {
    std::cout << "Running in NORMALIZED mode: performing background subtraction"
              << std::endl;

    TH1F *hist_CdShield_BkgSubtracted = static_cast<TH1F *>(
        hist_CdShieldSignal->Clone("hist_CdShield_BkgSubtracted"));
    hist_CdShield_BkgSubtracted->Add(hist_CdShieldBackground, -1);

    TH1F *zoomedHist_CdShield_BkgSubtracted = static_cast<TH1F *>(
        zoomedHist_CdShieldSignal->Clone("zoomedHist_CdShield_BkgSubtracted"));
    zoomedHist_CdShield_BkgSubtracted->Add(zoomedHist_CdShieldBackground, -1);

    TH1F *hist_CuShield_01132026_BkgSubtracted =
        static_cast<TH1F *>(hist_CuShieldSignal_01132026->Clone(
            "hist_CuShield_01132026_BkgSubtracted"));
    hist_CuShield_01132026_BkgSubtracted->Add(hist_CuShieldBackground_01132026,
                                              -1);

    TH1F *zoomedHist_CuShield_01132026_BkgSubtracted =
        static_cast<TH1F *>(zoomedHist_CuShieldSignal_01132026->Clone(
            "zoomedHist_CuShield_01132026_BkgSubtracted"));
    zoomedHist_CuShield_01132026_BkgSubtracted->Add(
        zoomedHist_CuShieldBackground_01132026, -1);

    TH1F *hist_CuShield_01142026_BkgSubtracted =
        static_cast<TH1F *>(hist_CuShieldSignal_01142026->Clone(
            "hist_CuShield_01142026_BkgSubtracted"));
    hist_CuShield_01142026_BkgSubtracted->Add(hist_CuShieldBackground_01142026,
                                              -1);

    TH1F *zoomedHist_CuShield_01142026_BkgSubtracted =
        static_cast<TH1F *>(zoomedHist_CuShieldSignal_01142026->Clone(
            "zoomedHist_CuShield_01142026_BkgSubtracted"));
    zoomedHist_CuShield_01142026_BkgSubtracted->Add(
        zoomedHist_CuShieldBackground_01142026, -1);

    TH1F *hist_Combined = static_cast<TH1F *>(
        hist_CdShield_BkgSubtracted->Clone("hist_Combined"));
    hist_Combined->Add(hist_CuShield_01132026_BkgSubtracted);
    hist_Combined->Add(hist_CuShield_01142026_BkgSubtracted);
    hist_Combined->SetTitle(
        Form("Background Subtracted; Energy [keV]; Counts / %d eV%s",
             Constants::BIN_WIDTH_EV, perSecond.Data()));

    TH1F *zoomedHist_Combined = static_cast<TH1F *>(
        zoomedHist_CdShield_BkgSubtracted->Clone("zoomedHist_Combined"));
    zoomedHist_Combined->Add(zoomedHist_CuShield_01132026_BkgSubtracted);
    zoomedHist_Combined->Add(zoomedHist_CuShield_01142026_BkgSubtracted);
    zoomedHist_Combined->SetTitle(
        Form("Background Subtracted; Energy [keV]; Counts / %d eV%s",
             Constants::BIN_WIDTH_EV, perSecond.Data()));

    PlottingUtils::ConfigureHistogram(hist_Combined, kBlue);
    TCanvas *canvas = new TCanvas("", "", 1200, 800);
    PlottingUtils::ConfigureCanvas(canvas);
    PlottingUtils::ConfigureAndDrawHistogram(zoomedHist_Combined, kBlue);
    PlottingUtils::SaveFigure(canvas, "background_subtracted_zoomed.png",
                              kFALSE);

    TString outputName = "Combined_BkgSubtracted" + suffix;
    TFile *outputFile =
        new TFile("root_files/" + outputName + ".root", "RECREATE");
    hist_Combined->Write("hist", TObject::kOverwrite);
    zoomedHist_Combined->Write("zoomedHist", TObject::kOverwrite);
    outputFile->Close();

    std::cout << "Background subtraction and combination complete!"
              << std::endl;
    std::cout << "Output saved to: root_files/" << outputName << ".root"
              << std::endl;

  } else {
    std::cout << "Running in NON-NORMALIZED mode: combining without subtraction"
              << std::endl;

    TH1F *hist_CombinedSignal =
        static_cast<TH1F *>(hist_CdShieldSignal->Clone("hist_CombinedSignal"));
    hist_CombinedSignal->Add(hist_CuShieldSignal_01132026);
    hist_CombinedSignal->Add(hist_CuShieldSignal_01142026);
    hist_CombinedSignal->SetTitle(
        Form("Combined Signal; Energy [keV]; Counts / %d eV%s",
             Constants::BIN_WIDTH_EV, perSecond.Data()));

    TH1F *zoomedHist_CombinedSignal = static_cast<TH1F *>(
        zoomedHist_CdShieldSignal->Clone("zoomedHist_CombinedSignal"));
    zoomedHist_CombinedSignal->Add(zoomedHist_CuShieldSignal_01132026);
    zoomedHist_CombinedSignal->Add(zoomedHist_CuShieldSignal_01142026);
    zoomedHist_CombinedSignal->SetTitle(
        Form("Combined Signal; Energy [keV]; Counts / %d eV%s",
             Constants::BIN_WIDTH_EV, perSecond.Data()));

    TH1F *hist_CombinedBackground = static_cast<TH1F *>(
        hist_CdShieldBackground->Clone("hist_CombinedBackground"));
    hist_CombinedBackground->Add(hist_CuShieldBackground_01132026);
    hist_CombinedBackground->Add(hist_CuShieldBackground_01142026);
    hist_CombinedBackground->SetTitle(
        Form("Combined Background; Energy [keV]; Counts / %d eV%s",
             Constants::BIN_WIDTH_EV, perSecond.Data()));

    TH1F *zoomedHist_CombinedBackground = static_cast<TH1F *>(
        zoomedHist_CdShieldBackground->Clone("zoomedHist_CombinedBackground"));
    zoomedHist_CombinedBackground->Add(zoomedHist_CuShieldBackground_01132026);
    zoomedHist_CombinedBackground->Add(zoomedHist_CuShieldBackground_01142026);
    zoomedHist_CombinedBackground->SetTitle(
        Form("Combined Background; Energy [keV]; Counts / %d eV%s",
             Constants::BIN_WIDTH_EV, perSecond.Data()));

    TCanvas *canvas = new TCanvas("", "", 1200, 800);
    PlottingUtils::ConfigureCanvas(canvas);
    PlottingUtils::ConfigureAndDrawHistogram(hist_CombinedSignal, kBlue);
    PlottingUtils::SaveFigure(canvas, "combined_signal.png", kFALSE);
    canvas->Clear();
    canvas->SetLogy(kFALSE);
    PlottingUtils::ConfigureAndDrawHistogram(zoomedHist_CombinedSignal, kBlue);
    PlottingUtils::SaveFigure(canvas, "combined_signal_zoomed.png", kFALSE);

    canvas->Clear();
    canvas->SetLogy(kTRUE);
    PlottingUtils::ConfigureAndDrawHistogram(hist_CombinedBackground, kRed);
    PlottingUtils::SaveFigure(canvas, "combined_background.png", kFALSE);
    canvas->Clear();
    canvas->SetLogy(kFALSE);
    PlottingUtils::ConfigureAndDrawHistogram(zoomedHist_CombinedBackground,
                                             kRed);
    PlottingUtils::SaveFigure(canvas, "combined_background_zoomed.png", kFALSE);

    TString outputNameSignal = "Combined_Signal" + suffix;
    TFile *outputFileSignal =
        new TFile("root_files/" + outputNameSignal + ".root", "RECREATE");
    hist_CombinedSignal->Write("hist", TObject::kOverwrite);
    zoomedHist_CombinedSignal->Write("zoomedHist", TObject::kOverwrite);
    outputFileSignal->Close();

    TString outputNameBackground = "Combined_Background" + suffix;
    TFile *outputFileBackground =
        new TFile("root_files/" + outputNameBackground + ".root", "RECREATE");
    hist_CombinedBackground->Write("hist", TObject::kOverwrite);
    zoomedHist_CombinedBackground->Write("zoomedHist", TObject::kOverwrite);
    outputFileBackground->Close();

    std::cout << "Combination complete!" << std::endl;
    std::cout << "Output saved to: root_files/" << outputNameSignal << ".root"
              << std::endl;
    std::cout << "Output saved to: root_files/" << outputNameBackground
              << ".root" << std::endl;
  }

  file_CdShieldSignal->Close();
  file_CdShieldBackground->Close();
  file_CuShieldSignal_01132026->Close();
  file_CuShieldBackground_01132026->Close();
  file_CuShieldSignal_01142026->Close();
  file_CuShieldBackground_01142026->Close();
}
