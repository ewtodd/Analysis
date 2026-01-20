#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TTree.h>

void SubtractBackground() {
  InitUtils::SetROOTPreferences();
  TString suffix = Constants::FILTERED ? "_filtered" : "";

  TString CdShield_25_Signal =
      "01132026-CdShield-GeSamplesIn-25Percent" + suffix;
  TString CdShield_25_Background =
      "01132026-ActiveBackground-25Percent" + suffix;

  TString CdShield_10_Signal =
      "01132026-CdShield-GeSamplesIn-10Percent" + suffix;
  TString CdShield_10_Background =
      "01132026-CdShield-ActiveBackground-10Percent" + suffix;

  TString CuShield_01132026_Signal =
      "01132026-CuShield-GeSamplesIn-10Percent" + suffix;
  TString CuShield_01132026_Background =
      "01132026-CuShield-ActiveBackground-Am241-10Percent" + suffix;

  TString CuShield_01142026_Signal =
      "01142026-CuShield-GeSamplesIn-10Percent" + suffix;
  TString CuShield_01142026_Background =
      "01142026-CuShield-ActiveBackground-10Percent" + suffix;

  TString NewSetup_Signal = "01152026-NewSetup-GeSamplesIn-5Percent" + suffix;
  TString NewSetup_Background =
      "01152026-NewSetup-ActiveBackground-5Percent" + suffix;

  TString NoShield_GraphiteCastle_Signal =
      "01162026-NoShield-GeSamplesIn-GraphiteCastle-10Percent" + suffix;
  TString NoShield_GraphiteCastle_Background =
      "01162026-NoShield-ActiveBackground-GraphiteCastle-10Percent" + suffix;

  Bool_t isFiltered = CdShield_10_Signal.Contains("_filtered");
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

  TFile *file_CdShield_25_Signal =
      new TFile("root_files/" + CdShield_25_Signal + ".root", "READ");
  if (!file_CdShield_25_Signal || file_CdShield_25_Signal->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << CdShield_25_Signal << ".root"
              << std::endl;
    return;
  }
  TH1F *hist_CdShield_25_Signal =
      static_cast<TH1F *>(file_CdShield_25_Signal->Get("hist"));
  TH1F *zoomedHist_CdShield_25_Signal =
      static_cast<TH1F *>(file_CdShield_25_Signal->Get("zoomedHist"));

  TFile *file_CdShield_25_Background =
      new TFile("root_files/" + CdShield_25_Background + ".root", "READ");
  if (!file_CdShield_25_Background || file_CdShield_25_Background->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << CdShield_25_Background << ".root"
              << std::endl;
    return;
  }
  TH1F *hist_CdShield_25_Background =
      static_cast<TH1F *>(file_CdShield_25_Background->Get("hist"));
  TH1F *zoomedHist_CdShield_25_Background =
      static_cast<TH1F *>(file_CdShield_25_Background->Get("zoomedHist"));

  TFile *file_CdShield_10_Signal =
      new TFile("root_files/" + CdShield_10_Signal + ".root", "READ");
  if (!file_CdShield_10_Signal || file_CdShield_10_Signal->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << CdShield_10_Signal << ".root"
              << std::endl;
    return;
  }
  TH1F *hist_CdShield_10_Signal =
      static_cast<TH1F *>(file_CdShield_10_Signal->Get("hist"));
  TH1F *zoomedHist_CdShield_10_Signal =
      static_cast<TH1F *>(file_CdShield_10_Signal->Get("zoomedHist"));

  TFile *file_CdShield_10_Background =
      new TFile("root_files/" + CdShield_10_Background + ".root", "READ");
  if (!file_CdShield_10_Background || file_CdShield_10_Background->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << CdShield_10_Background << ".root"
              << std::endl;
    return;
  }
  TH1F *hist_CdShield_10_Background =
      static_cast<TH1F *>(file_CdShield_10_Background->Get("hist"));
  TH1F *zoomedHist_CdShield_10_Background =
      static_cast<TH1F *>(file_CdShield_10_Background->Get("zoomedHist"));

  TFile *file_CuShield_01132026_Signal =
      new TFile("root_files/" + CuShield_01132026_Signal + ".root", "READ");
  if (!file_CuShield_01132026_Signal ||
      file_CuShield_01132026_Signal->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << CuShield_01132026_Signal << ".root"
              << std::endl;
    return;
  }
  TH1F *hist_CuShield_01132026_Signal =
      static_cast<TH1F *>(file_CuShield_01132026_Signal->Get("hist"));
  TH1F *zoomedHist_CuShield_01132026_Signal =
      static_cast<TH1F *>(file_CuShield_01132026_Signal->Get("zoomedHist"));

  TFile *file_CuShield_01132026_Background =
      new TFile("root_files/" + CuShield_01132026_Background + ".root", "READ");
  if (!file_CuShield_01132026_Background ||
      file_CuShield_01132026_Background->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << CuShield_01132026_Background
              << ".root" << std::endl;
    return;
  }
  TH1F *hist_CuShield_01132026_Background =
      static_cast<TH1F *>(file_CuShield_01132026_Background->Get("hist"));
  TH1F *zoomedHist_CuShield_01132026_Background = static_cast<TH1F *>(
      file_CuShield_01132026_Background->Get("zoome dHist"));

  TFile *file_CuShield_01142026_Signal =
      new TFile("root_files/" + CuShield_01142026_Signal + ".root", "READ");
  if (!file_CuShield_01142026_Signal ||
      file_CuShield_01142026_Signal->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << CuShield_01142026_Signal << ".root"
              << std::endl;
    return;
  }
  TH1F *hist_CuShield_01142026_Signal =
      static_cast<TH1F *>(file_CuShield_01142026_Signal->Get("hist"));
  TH1F *zoomedHist_CuShield_01142026_Signal =
      static_cast<TH1F *>(file_CuShield_01142026_Signal->Get("zoomedHist"));

  TFile *file_CuShield_01142026_Background =
      new TFile("root_files/" + CuShield_01142026_Background + ".root", "READ");
  if (!file_CuShield_01142026_Background ||
      file_CuShield_01142026_Background->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << CuShield_01142026_Background
              << ".root" << std::endl;
    return;
  }
  TH1F *hist_CuShield_01142026_Background =
      static_cast<TH1F *>(file_CuShield_01142026_Background->Get("hist"));
  TH1F *zoomedHist_CuShield_01142026_Background =
      static_cast<TH1F *>(file_CuShield_01142026_Background->Get("zoomedHist"));

  TFile *file_NewSetup_Signal =
      new TFile("root_files/" + NewSetup_Signal + ".root", "READ");
  if (!file_NewSetup_Signal || file_NewSetup_Signal->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << NewSetup_Signal << ".root"
              << std::endl;
    return;
  }
  TH1F *hist_NewSetup_Signal =
      static_cast<TH1F *>(file_NewSetup_Signal->Get("hist"));
  TH1F *zoomedHist_NewSetup_Signal =
      static_cast<TH1F *>(file_NewSetup_Signal->Get("zoomedHist"));

  TFile *file_NewSetup_Background =
      new TFile("root_files/" + NewSetup_Background + ".root", "READ");
  if (!file_NewSetup_Background || file_NewSetup_Background->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << NewSetup_Background << ".root"
              << std::endl;
    return;
  }
  TH1F *hist_NewSetup_Background =
      static_cast<TH1F *>(file_NewSetup_Background->Get("hist"));
  TH1F *zoomedHist_NewSetup_Background =
      static_cast<TH1F *>(file_NewSetup_Background->Get("zoomedHist"));

  TFile *file_NoShield_GraphiteCastle_Signal = new TFile(
      "root_files/" + NoShield_GraphiteCastle_Signal + ".root", "READ");
  if (!file_NoShield_GraphiteCastle_Signal ||
      file_NoShield_GraphiteCastle_Signal->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << NoShield_GraphiteCastle_Signal
              << ".root" << std::endl;
    return;
  }
  TH1F *hist_NoShield_GraphiteCastle_Signal =
      static_cast<TH1F *>(file_NoShield_GraphiteCastle_Signal->Get("hist"));
  TH1F *zoomedHist_NoShield_GraphiteCastle_Signal = static_cast<TH1F *>(
      file_NoShield_GraphiteCastle_Signal->Get("zoomedHist"));

  TFile *file_NoShield_GraphiteCastle_Background = new TFile(
      "root_files/" + NoShield_GraphiteCastle_Background + ".root", "READ");
  if (!file_NoShield_GraphiteCastle_Background ||
      file_NoShield_GraphiteCastle_Background->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << NoShield_GraphiteCastle_Background
              << ".root" << std::endl;
    return;
  }
  TH1F *hist_NoShield_GraphiteCastle_Background =
      static_cast<TH1F *>(file_NoShield_GraphiteCastle_Background->Get("hist"));
  TH1F *zoomedHist_NoShield_GraphiteCastle_Background = static_cast<TH1F *>(
      file_NoShield_GraphiteCastle_Background->Get("zoomedHist"));

  if (!Constants::NORMALIZE_BY_TIME) {
    std::cerr << "UNSUPPORTED in current version." << std::endl;
  }

  TString perSecond = Constants::NORMALIZE_BY_TIME ? " / s" : "";

  if (Constants::NORMALIZE_BY_TIME) {
    std::cout << "Running in NORMALIZED mode: performing background subtraction"
              << std::endl;

    TH1F *hist_CdShield_25_BkgSubtracted = static_cast<TH1F *>(
        hist_CdShield_25_Signal->Clone("hist_CdShield_25_BkgSubtracted"));
    hist_CdShield_25_BkgSubtracted->Add(hist_CdShield_25_Background, -1);

    TH1F *zoomedHist_CdShield_25_BkgSubtracted =
        static_cast<TH1F *>(zoomedHist_CdShield_25_Signal->Clone(
            "zoomedHist_CdShield_25_BkgSubtracted"));
    zoomedHist_CdShield_25_BkgSubtracted->Add(zoomedHist_CdShield_25_Background,
                                              -1);

    TH1F *hist_CdShield_10_BkgSubtracted = static_cast<TH1F *>(
        hist_CdShield_10_Signal->Clone("hist_CdShield_10_BkgSubtracted"));
    hist_CdShield_10_BkgSubtracted->Add(hist_CdShield_10_Background, -1);

    TH1F *zoomedHist_CdShield_10_BkgSubtracted =
        static_cast<TH1F *>(zoomedHist_CdShield_10_Signal->Clone(
            "zoomedHist_CdShield_10_BkgSubtracted"));
    zoomedHist_CdShield_10_BkgSubtracted->Add(zoomedHist_CdShield_10_Background,
                                              -1);

    TH1F *hist_CuShield_01132026_BkgSubtracted =
        static_cast<TH1F *>(hist_CuShield_01132026_Signal->Clone(
            "hist_CuShield_01132026_BkgSubtracted"));
    hist_CuShield_01132026_BkgSubtracted->Add(hist_CuShield_01132026_Background,
                                              -1);

    TH1F *zoomedHist_CuShield_01132026_BkgSubtracted =
        static_cast<TH1F *>(zoomedHist_CuShield_01132026_Signal->Clone(
            "zoomedHist_CuShield_01132026_BkgSubtracted"));
    zoomedHist_CuShield_01132026_BkgSubtracted->Add(
        zoomedHist_CuShield_01132026_Background, -1);

    TH1F *hist_CuShield_01142026_BkgSubtracted =
        static_cast<TH1F *>(hist_CuShield_01142026_Signal->Clone(
            "hist_CuShield_01142026_BkgSubtracted"));
    hist_CuShield_01142026_BkgSubtracted->Add(hist_CuShield_01142026_Background,
                                              -1);

    TH1F *zoomedHist_CuShield_01142026_BkgSubtracted =
        static_cast<TH1F *>(zoomedHist_CuShield_01142026_Signal->Clone(
            "zoomedHist_CuShield_01142026_BkgSubtracted"));
    zoomedHist_CuShield_01142026_BkgSubtracted->Add(
        zoomedHist_CuShield_01142026_Background, -1);

    TH1F *hist_NewSetup_BkgSubtracted = static_cast<TH1F *>(
        hist_NewSetup_Signal->Clone("hist_NewSetup_BkgSubtracted"));
    hist_NewSetup_BkgSubtracted->Add(hist_NewSetup_Background, -1);

    TH1F *zoomedHist_NewSetup_BkgSubtracted = static_cast<TH1F *>(
        zoomedHist_NewSetup_Signal->Clone("zoomedHist_NewSetup_BkgSubtracted"));
    zoomedHist_NewSetup_BkgSubtracted->Add(zoomedHist_NewSetup_Background, -1);

    TH1F *hist_NoShield_GraphiteCastle_BkgSubtracted =
        static_cast<TH1F *>(hist_NoShield_GraphiteCastle_Signal->Clone(
            "hist_NoShield_GraphiteCastle_BkgSubtracted"));
    hist_NoShield_GraphiteCastle_BkgSubtracted->Add(
        hist_NoShield_GraphiteCastle_Background, -1);

    TH1F *zoomedHist_NoShield_GraphiteCastle_BkgSubtracted =
        static_cast<TH1F *>(zoomedHist_NoShield_GraphiteCastle_Signal->Clone(
            "zoomedHist_NoShield_GraphiteCastle_BkgSubtracted"));
    zoomedHist_NoShield_GraphiteCastle_BkgSubtracted->Add(
        zoomedHist_NoShield_GraphiteCastle_Background, -1);

    TH1F *hist_Combined = static_cast<TH1F *>(
        hist_CdShield_10_BkgSubtracted->Clone("hist_Combined"));
    hist_Combined->Reset();
    hist_Combined->Add(hist_CdShield_25_BkgSubtracted);
    hist_Combined->Add(hist_CdShield_10_BkgSubtracted);
    hist_Combined->Add(hist_CuShield_01132026_BkgSubtracted);
    hist_Combined->Add(hist_CuShield_01142026_BkgSubtracted);
    hist_Combined->Add(hist_NewSetup_BkgSubtracted);
    hist_Combined->Add(hist_NoShield_GraphiteCastle_BkgSubtracted);
    hist_Combined->SetTitle(
        Form("Background Subtracted; Energy [keV]; Counts / %d eV%s",
             Constants::BIN_WIDTH_EV, perSecond.Data()));

    TCanvas *canvasFull = new TCanvas("", "", 1200, 800);
    PlottingUtils::ConfigureCanvas(canvasFull);
    PlottingUtils::ConfigureHistogram(hist_Combined, kP10Violet);
    hist_Combined->Draw("HIST");
    hist_Combined->SetFillStyle(0);
    hist_Combined->SetLineWidth(2);
    PlottingUtils::SaveFigure(canvasFull, "background_subtracted.png", kFALSE);

    TH1F *zoomedHist_Combined = static_cast<TH1F *>(
        zoomedHist_CdShield_10_BkgSubtracted->Clone("zoomedHist_Combined"));
    zoomedHist_Combined->Add(zoomedHist_CdShield_25_BkgSubtracted);
    zoomedHist_Combined->Add(zoomedHist_CuShield_01132026_BkgSubtracted);
    zoomedHist_Combined->Add(zoomedHist_CuShield_01142026_BkgSubtracted);
    zoomedHist_Combined->Add(zoomedHist_NewSetup_BkgSubtracted);
    zoomedHist_Combined->Add(zoomedHist_NoShield_GraphiteCastle_BkgSubtracted);
    zoomedHist_Combined->SetTitle(
        Form("Background Subtracted; Energy [keV]; Counts / %d eV%s",
             Constants::BIN_WIDTH_EV, perSecond.Data()));

    TCanvas *canvas = new TCanvas("", "", 1200, 800);
    PlottingUtils::ConfigureCanvas(canvas);
    PlottingUtils::ConfigureHistogram(zoomedHist_Combined, kP10Violet);
    zoomedHist_Combined->Draw("HIST");
    zoomedHist_Combined->SetFillStyle(0);
    zoomedHist_Combined->SetLineWidth(2);
    PlottingUtils::SaveFigure(canvas, "background_subtracted_zoomed.png",
                              kFALSE);

    TString outputName = "Combined_BkgSubtracted" + suffix;
    TFile *outputFile =
        new TFile("root_files/" + outputName + ".root", "RECREATE");

    hist_Combined->Write("histCombined", TObject::kOverwrite);
    zoomedHist_Combined->Write("zoomedHistCombined", TObject::kOverwrite);
    zoomedHist_CdShield_25_BkgSubtracted->Write(
        "zoomedHist_CdShield_25_BkgSubtracted", TObject::kOverwrite);
    zoomedHist_CdShield_10_BkgSubtracted->Write(
        "zoomedHist_CdShield_10_BkgSubtracted", TObject::kOverwrite);
    zoomedHist_CuShield_01132026_BkgSubtracted->Write(
        "zoomedHist_CuShield_01132026_BkgSubtracted", TObject::kOverwrite);
    zoomedHist_CuShield_01142026_BkgSubtracted->Write(
        "zoomedHist_CuShield_01142026_BkgSubtracted", TObject::kOverwrite);
    zoomedHist_NewSetup_BkgSubtracted->Write(
        "zoomedHist_NewSetup_BkgSubtracted", TObject::kOverwrite);
    zoomedHist_NoShield_GraphiteCastle_BkgSubtracted->Write(
        "zoomedHist_NoShield_GraphiteCastle_BkgSubtracted",
        TObject::kOverwrite);

    outputFile->Close();

    std::cout << "Background subtraction and combination complete!"
              << std::endl;
    std::cout << "Output saved to: root_files/" << outputName << ".root"
              << std::endl;
  } else {
    std::cerr << "Running in NON-NORMALIZED mode: use other macro" << std::endl;
  }

  file_CdShield_25_Signal->Close();
  file_CdShield_25_Background->Close();
  file_CdShield_10_Signal->Close();
  file_CdShield_10_Background->Close();
  file_CuShield_01132026_Signal->Close();
  file_CuShield_01132026_Background->Close();
  file_CuShield_01142026_Signal->Close();
  file_CuShield_01142026_Background->Close();
  file_NewSetup_Signal->Close();
  file_NewSetup_Background->Close();
  file_NoShield_GraphiteCastle_Signal->Close();
  file_NoShield_GraphiteCastle_Background->Close();
}
