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

  if (!Constants::FILTERED)
    std::cerr << "Use other macro for unfiltered data." << std::endl;

  TString CdShield_25_Signal =
      "01132026-CdShield-GeSamplesIn-25Percent" + suffix;

  TString CdShield_10_Signal =
      "01132026-CdShield-GeSamplesIn-10Percent" + suffix;

  TString CuShield_01132026_Signal =
      "01132026-CuShield-GeSamplesIn-10Percent" + suffix;

  TString CuShield_01142026_Signal =
      "01142026-CuShield-GeSamplesIn-10Percent" + suffix;

  TString CuShield_01142026_Signal_90Percent =
      "01142026-CuShield-GeSamplesIn-MovedBack-90Percent" + suffix;

  TString NewSetup_Signal = "01152026-NewSetup-GeSamplesIn-5Percent" + suffix;

  TString NoShield_GraphiteCastle_Signal =
      "01162026-NoShield-GeSamplesIn-GraphiteCastle-10Percent" + suffix;

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

  TFile *file_CuShield_01142026_Signal_90Percent = new TFile(
      "root_files/" + CuShield_01142026_Signal_90Percent + ".root", "READ");
  if (!file_CuShield_01142026_Signal_90Percent ||
      file_CuShield_01142026_Signal_90Percent->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << CuShield_01142026_Signal_90Percent
              << ".root" << std::endl;
    return;
  }
  TH1F *hist_CuShield_01142026_Signal_90Percent =
      static_cast<TH1F *>(file_CuShield_01142026_Signal_90Percent->Get("hist"));
  TH1F *zoomedHist_CuShield_01142026_Signal_90Percent = static_cast<TH1F *>(
      file_CuShield_01142026_Signal_90Percent->Get("zoomedHist"));

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

  if (!Constants::NORMALIZE_BY_TIME) {
    std::cerr << "UNSUPPORTED in current version." << std::endl;
  }

  TString perSecond = Constants::NORMALIZE_BY_TIME ? " / s" : "";

  if (Constants::NORMALIZE_BY_TIME) {
    std::cerr << "Use other macro for normalized by time!" << std::endl;
  }

  TH1F *hist_Combined =
      static_cast<TH1F *>(hist_CdShield_25_Signal->Clone("hist_Combined"));
  hist_Combined->Add(hist_CdShield_10_Signal);
  hist_Combined->Add(hist_CuShield_01132026_Signal);
  hist_Combined->Add(hist_CuShield_01142026_Signal);
  hist_Combined->Add(hist_CuShield_01142026_Signal_90Percent);
  hist_Combined->Add(hist_NewSetup_Signal);
  hist_Combined->Add(hist_NoShield_GraphiteCastle_Signal);
  hist_Combined->SetTitle(Form("; Energy [keV]; Counts / %d eV%s",
                               Constants::BIN_WIDTH_EV, perSecond.Data()));

  TCanvas *canvasFull = new TCanvas("", "", 1200, 800);
  PlottingUtils::ConfigureCanvas(canvasFull);
  PlottingUtils::ConfigureHistogram(hist_Combined, kP10Violet);
  hist_Combined->Draw("HIST");
  hist_Combined->SetFillStyle(0);
  hist_Combined->SetLineWidth(2);
  PlottingUtils::SaveFigure(canvasFull, "combined.png", kFALSE);

  TH1F *zoomedHist_Combined = static_cast<TH1F *>(
      zoomedHist_CdShield_25_Signal->Clone("zoomedHist_Combined"));
  zoomedHist_Combined->Add(zoomedHist_CdShield_10_Signal);
  zoomedHist_Combined->Add(zoomedHist_CuShield_01132026_Signal);
  zoomedHist_Combined->Add(zoomedHist_CuShield_01142026_Signal);
  zoomedHist_Combined->Add(zoomedHist_CuShield_01142026_Signal_90Percent);
  zoomedHist_Combined->Add(zoomedHist_NewSetup_Signal);
  zoomedHist_Combined->Add(zoomedHist_NoShield_GraphiteCastle_Signal);
  zoomedHist_Combined->SetTitle(Form("; Energy [keV]; Counts / %d eV%s",
                                     Constants::BIN_WIDTH_EV,
                                     perSecond.Data()));

  TCanvas *canvas = new TCanvas("", "", 1200, 800);
  PlottingUtils::ConfigureCanvas(canvas);
  PlottingUtils::ConfigureHistogram(zoomedHist_Combined, kP10Violet);
  zoomedHist_Combined->Draw("HIST");
  zoomedHist_Combined->SetFillStyle(0);
  zoomedHist_Combined->SetLineWidth(2);
  PlottingUtils::SaveFigure(canvas, "combined_zoomed.png", kFALSE);

  TString outputName = "Combined" + suffix;
  TFile *outputFile =
      new TFile("root_files/" + outputName + ".root", "RECREATE");

  hist_Combined->Write("histCombined", TObject::kOverwrite);
  zoomedHist_Combined->Write("zoomedHistCombined", TObject::kOverwrite);
  zoomedHist_CdShield_25_Signal->Write("zoomedHist_CdShield_25",
                                       TObject::kOverwrite);
  zoomedHist_CdShield_10_Signal->Write("zoomedHist_CdShield_10",
                                       TObject::kOverwrite);
  zoomedHist_CuShield_01132026_Signal->Write("zoomedHist_CuShield_01132026",
                                             TObject::kOverwrite);
  zoomedHist_CuShield_01142026_Signal->Write("zoomedHist_CuShield_01142026",
                                             TObject::kOverwrite);
  zoomedHist_CuShield_01142026_Signal_90Percent->Write(
      "zoomedHist_CuShield_01142026_90Percent", TObject::kOverwrite);
  zoomedHist_NewSetup_Signal->Write("zoomedHist_NewSetup", TObject::kOverwrite);
  zoomedHist_NoShield_GraphiteCastle_Signal->Write(
      "zoomedHist_NoShield_GraphiteCastle", TObject::kOverwrite);

  outputFile->Close();

  std::cout << "Background subtraction and combination complete!" << std::endl;
  std::cout << "Output saved to: root_files/" << outputName << ".root"
            << std::endl;

  file_CdShield_25_Signal->Close();
  file_CdShield_10_Signal->Close();
  file_CuShield_01132026_Signal->Close();
  file_CuShield_01142026_Signal->Close();
  file_CuShield_01142026_Signal_90Percent->Close();
  file_NewSetup_Signal->Close();
  file_NoShield_GraphiteCastle_Signal->Close();
}
