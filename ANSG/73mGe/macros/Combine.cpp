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
  TString CdShieldSignal_25 =
      "01132026-CdShield-GeSamplesIn-25Percent" + suffix;
  TString CdShieldBackground_25 =
      "01132026-ActiveBackground-25Percent" + suffix;
  TString CdShieldSignal_10 =
      "01132026-CdShield-GeSamplesIn-10Percent" + suffix;
  TString CdShieldBackground_10 =
      "01132026-CdShield-ActiveBackground-10Percent" + suffix;
  TString CuShieldSignal_01132026 =
      "01132026-CuShield-GeSamplesIn-10Percent" + suffix;
  TString CuShieldBackground_01132026 =
      "01132026-CuShield-ActiveBackground-Am241-10Percent" + suffix;
  TString CuShieldSignal_01142026 =
      "01142026-CuShield-GeSamplesIn-10Percent" + suffix;
  TString CuShieldBackground_01142026 =
      "01142026-CuShield-ActiveBackground-10Percent" + suffix;
  TString NewSetupSignal = "01152026-NewSetup-GeSamplesIn-5Percent" + suffix;
  TString NewSetupBackground =
      "01152026-NewSetup-ActiveBackground-5Percent" + suffix;

  Bool_t isFiltered = CdShieldSignal_10.Contains("_filtered");
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

  TFile *file_CdShieldSignal_25 =
      new TFile("root_files/" + CdShieldSignal_25 + ".root", "READ");
  if (!file_CdShieldSignal_25 || file_CdShieldSignal_25->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << CdShieldSignal_25 << ".root"
              << std::endl;
    return;
  }
  TH1F *hist_CdShieldSignal_25 =
      static_cast<TH1F *>(file_CdShieldSignal_25->Get("hist"));
  TH1F *zoomedHist_CdShieldSignal_25 =
      static_cast<TH1F *>(file_CdShieldSignal_25->Get("zoomedHist"));

  TFile *file_CdShieldBackground_25 =
      new TFile("root_files/" + CdShieldBackground_25 + ".root", "READ");
  if (!file_CdShieldBackground_25 || file_CdShieldBackground_25->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << CdShieldBackground_25 << ".root"
              << std::endl;
    return;
  }

  TH1F *hist_CdShieldBackground_25 =
      static_cast<TH1F *>(file_CdShieldBackground_25->Get("hist"));
  TH1F *zoomedHist_CdShieldBackground_25 =
      static_cast<TH1F *>(file_CdShieldBackground_25->Get("zoomedHist"));

  TFile *file_CdShieldSignal_10 =
      new TFile("root_files/" + CdShieldSignal_10 + ".root", "READ");
  if (!file_CdShieldSignal_10 || file_CdShieldSignal_10->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << CdShieldSignal_10 << ".root"
              << std::endl;
    return;
  }

  TH1F *hist_CdShieldSignal_10 =
      static_cast<TH1F *>(file_CdShieldSignal_10->Get("hist"));
  TH1F *zoomedHist_CdShieldSignal_10 =
      static_cast<TH1F *>(file_CdShieldSignal_10->Get("zoomedHist"));
  TFile *file_CdShieldBackground_10 =
      new TFile("root_files/" + CdShieldBackground_10 + ".root", "READ");
  if (!file_CdShieldBackground_10 || file_CdShieldBackground_10->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << CdShieldBackground_10 << ".root"
              << std::endl;
    return;
  }
  TH1F *hist_CdShieldBackground_10 =
      static_cast<TH1F *>(file_CdShieldBackground_10->Get("hist"));
  TH1F *zoomedHist_CdShieldBackground_10 =
      static_cast<TH1F *>(file_CdShieldBackground_10->Get("zoomedHist"));
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
  TFile *file_NewSetupSignal =
      new TFile("root_files/" + NewSetupSignal + ".root", "READ");
  if (!file_NewSetupSignal || file_NewSetupSignal->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << NewSetupSignal << ".root"
              << std::endl;
    return;
  }
  TH1F *hist_NewSetupSignal =
      static_cast<TH1F *>(file_NewSetupSignal->Get("hist"));
  TH1F *zoomedHist_NewSetupSignal =
      static_cast<TH1F *>(file_NewSetupSignal->Get("zoomedHist"));
  TFile *file_NewSetupBackground =
      new TFile("root_files/" + NewSetupBackground + ".root", "READ");
  if (!file_NewSetupBackground || file_NewSetupBackground->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << NewSetupBackground << ".root"
              << std::endl;
    return;
  }
  TH1F *hist_NewSetupBackground =
      static_cast<TH1F *>(file_NewSetupBackground->Get("hist"));
  TH1F *zoomedHist_NewSetupBackground =
      static_cast<TH1F *>(file_NewSetupBackground->Get("zoomedHist"));

  if (!Constants::NORMALIZE_BY_TIME) {
    std::cerr << "UNSUPPORTED in current version." << std::endl;
  }

  TString perSecond = Constants::NORMALIZE_BY_TIME ? " / s" : "";

  if (Constants::NORMALIZE_BY_TIME) {
    std::cout << "Running in NORMALIZED mode: performing background subtraction"
              << std::endl;

    TH1F *hist_CdShield_25_BkgSubtracted = static_cast<TH1F *>(
        hist_CdShieldSignal_25->Clone("hist_CdShield_25_BkgSubtracted"));
    hist_CdShield_25_BkgSubtracted->Add(hist_CdShieldBackground_25, -1);

    TH1F *zoomedHist_CdShield_25_BkgSubtracted =
        static_cast<TH1F *>(zoomedHist_CdShieldSignal_25->Clone(
            "zoomedHist_CdShield_25_BkgSubtracted"));
    zoomedHist_CdShield_25_BkgSubtracted->Add(zoomedHist_CdShieldBackground_25,
                                              -1);

    TH1F *hist_CdShield_10_BkgSubtracted =
        static_cast<TH1F *>(hist_CdShieldSignal_10->Clone("hist_"
                                                          "CdShield_10_"
                                                          "BkgSubtracte"
                                                          "d"));
    hist_CdShield_10_BkgSubtracted->Add(hist_CdShieldBackground_10, -1);

    TH1F *zoomedHist_CdShield_10_BkgSubtracted =
        static_cast<TH1F *>(zoomedHist_CdShieldSignal_10->Clone("zoomedHist_"
                                                                "CdShield_10_"
                                                                "BkgSubtracte"
                                                                "d"));
    zoomedHist_CdShield_10_BkgSubtracted->Add(zoomedHist_CdShieldBackground_10,
                                              -1);

    TH1F *hist_CuShield_01132026_BkgSubtracted =
        static_cast<TH1F *>(hist_CuShieldSignal_01132026->Clone("hist_"
                                                                "CuShield_"
                                                                "01132026_"
                                                                "BkgSubtracte"
                                                                "d"));
    hist_CuShield_01132026_BkgSubtracted->Add(hist_CuShieldBackground_01132026,
                                              -1);

    TH1F *zoomedHist_CuShield_01132026_BkgSubtracted = static_cast<TH1F *>(
        zoomedHist_CuShieldSignal_01132026->Clone("zoomedHist_"
                                                  "CuShield_"
                                                  "01132026_"
                                                  "BkgSubtracte"
                                                  "d"));
    zoomedHist_CuShield_01132026_BkgSubtracted->Add(
        zoomedHist_CuShieldBackground_01132026, -1);

    TH1F *hist_CuShield_01142026_BkgSubtracted =
        static_cast<TH1F *>(hist_CuShieldSignal_01142026->Clone("hist_"
                                                                "CuShield_"
                                                                "01142026_"
                                                                "BkgSubtracte"
                                                                "d"));
    hist_CuShield_01142026_BkgSubtracted->Add(hist_CuShieldBackground_01142026,
                                              -1);

    TH1F *zoomedHist_CuShield_01142026_BkgSubtracted = static_cast<TH1F *>(
        zoomedHist_CuShieldSignal_01142026->Clone("zoomedHist_"
                                                  "CuShield_"
                                                  "01142026_"
                                                  "BkgSubtracte"
                                                  "d"));
    zoomedHist_CuShield_01142026_BkgSubtracted->Add(
        zoomedHist_CuShieldBackground_01142026, -1);

    TH1F *hist_NewSetup_BkgSubtracted =
        static_cast<TH1F *>(hist_NewSetupSignal->Clone("hist_"
                                                       "NewSetup"
                                                       "_BkgSubt"
                                                       "racte"
                                                       "d"));
    hist_NewSetup_BkgSubtracted->Add(hist_NewSetupBackground, -1);

    TH1F *zoomedHist_NewSetup_BkgSubtracted =
        static_cast<TH1F *>(zoomedHist_NewSetupSignal->Clone("zoomedHist_"
                                                             "NewSetup_"
                                                             "BkgSubtracte"
                                                             "d"));
    zoomedHist_NewSetup_BkgSubtracted->Add(zoomedHist_NewSetupBackground, -1);

    TH1F *hist_Combined =
        static_cast<TH1F *>(hist_CdShield_10_BkgSubtracted->Clone("hist_"
                                                                  "Combined"));
    hist_Combined->Add(hist_CuShield_01132026_BkgSubtracted);
    hist_Combined->Add(hist_CuShield_01142026_BkgSubtracted);
    hist_Combined->Add(hist_NewSetup_BkgSubtracted);
    hist_Combined->SetTitle(Form("Background Subtracted; "
                                 "Energy [keV]; Counts / "
                                 "%d eV%s",
                                 Constants::BIN_WIDTH_EV, perSecond.Data()));

    TH1F *zoomedHist_Combined0 = static_cast<TH1F *>(
        zoomedHist_CdShield_25_BkgSubtracted->Clone("zoomedHist_"
                                                    "Combined"));
    zoomedHist_Combined0->Add(zoomedHist_CdShield_10_BkgSubtracted);
    zoomedHist_Combined0->Add(zoomedHist_CuShield_01132026_BkgSubtracted);
    zoomedHist_Combined0->Add(zoomedHist_CuShield_01142026_BkgSubtracted);
    zoomedHist_Combined0->Add(zoomedHist_NewSetup_BkgSubtracted);
    zoomedHist_Combined0->SetTitle(Form("Background "
                                        "Subtracted; Energy "
                                        "[keV]; Counts / %d "
                                        "eV%s",
                                        Constants::BIN_WIDTH_EV,
                                        perSecond.Data()));

    TH1F *zoomedHist_Combined = static_cast<TH1F *>(
        zoomedHist_CdShield_10_BkgSubtracted->Clone("zoomedHist_"
                                                    "Combined"));
    zoomedHist_Combined->Add(zoomedHist_CuShield_01132026_BkgSubtracted);
    zoomedHist_Combined->Add(zoomedHist_CuShield_01142026_BkgSubtracted);
    zoomedHist_Combined->Add(zoomedHist_NewSetup_BkgSubtracted);
    zoomedHist_Combined->SetTitle(Form("Background "
                                       "Subtracted; Energy "
                                       "[keV]; Counts / %d "
                                       "eV%s",
                                       Constants::BIN_WIDTH_EV,
                                       perSecond.Data()));

    TH1F *zoomedHist_Combined2 = static_cast<TH1F *>(
        zoomedHist_CuShield_01132026_BkgSubtracted->Clone("zoomedHist_"
                                                          "Combined"));
    zoomedHist_Combined2->Add(zoomedHist_CuShield_01142026_BkgSubtracted);
    zoomedHist_Combined2->Add(zoomedHist_NewSetup_BkgSubtracted);
    zoomedHist_Combined2->SetTitle(Form("Background "
                                        "Subtracted; Energy "
                                        "[keV]; Counts / %d "
                                        "eV%s",
                                        Constants::BIN_WIDTH_EV,
                                        perSecond.Data()));

    TH1F *zoomedHist_Combined3 = static_cast<TH1F *>(
        zoomedHist_CuShield_01142026_BkgSubtracted->Clone("zoomedHist_"
                                                          "Combined"));
    zoomedHist_Combined3->Add(zoomedHist_NewSetup_BkgSubtracted);
    zoomedHist_Combined3->SetTitle(Form("Background "
                                        "Subtracted; Energy "
                                        "[keV]; Counts / %d "
                                        "eV%s",
                                        Constants::BIN_WIDTH_EV,
                                        perSecond.Data()));

    TCanvas *canvas = new TCanvas("", "", 1200, 800);
    PlottingUtils::ConfigureCanvas(canvas);
    PlottingUtils::ConfigureHistogram(zoomedHist_Combined0, kBlack);
    PlottingUtils::ConfigureHistogram(zoomedHist_Combined, kP10Violet);
    PlottingUtils::ConfigureHistogram(zoomedHist_Combined2, kRed);
    PlottingUtils::ConfigureHistogram(zoomedHist_Combined3, kBlue);
    zoomedHist_Combined0->Draw("HIST");
    zoomedHist_Combined->Draw("HIST SAME");
    zoomedHist_Combined2->Draw("HIST SAME");
    zoomedHist_Combined3->Draw("HIST SAME");
    zoomedHist_Combined0->GetYaxis()->SetRangeUser(-0.006, 0.007);
    PlottingUtils::SaveFigure(canvas,
                              "background_subtracted_"
                              "zoomed.png",
                              kFALSE);

    TString outputName = "Combined_BkgSubtracted" + suffix;
    TFile *outputFile =
        new TFile("root_files/" + outputName + ".root", "RECREATE");
    hist_Combined->Write("histCombined", TObject::kOverwrite);
    zoomedHist_Combined->Write("zoomedHistCombined", TObject::kOverwrite);
    zoomedHist_CdShield_10_BkgSubtracted->Write("zoomedHist_CdShield_"
                                                "10_BkgSubtracted",
                                                TObject::kOverwrite);
    zoomedHist_CuShield_01132026_BkgSubtracted->Write("zoomedHist_CuShield_"
                                                      "01132026_"
                                                      "BkgSubtracted",
                                                      TObject::kOverwrite);
    zoomedHist_CuShield_01142026_BkgSubtracted->Write("zoomedHist_CuShield_"
                                                      "01142026_"
                                                      "BkgSubtracted",
                                                      TObject::kOverwrite);
    zoomedHist_NewSetup_BkgSubtracted->Write("zoomedHist_NewSetup_"
                                             "BkgSubtracted",
                                             TObject::kOverwrite);

    outputFile->Close();

    std::cout << "Background "
                 "subtraction and "
                 "combination complete!"
              << std::endl;
    std::cout << "Output saved "
                 "to: root_files/"
              << outputName << ".root" << std::endl;
  } else {
    std::cerr << "Running in NON-NORMALIZED mode: combining without subtraction"
              << std::endl;
  }
  file_CdShieldSignal_10->Close();
  file_CdShieldBackground_10->Close();
  file_CuShieldSignal_01132026->Close();
  file_CuShieldBackground_01132026->Close();
  file_CuShieldSignal_01142026->Close();
  file_CuShieldBackground_01142026->Close();
  file_NewSetupSignal->Close();
  file_NewSetupBackground->Close();
}
