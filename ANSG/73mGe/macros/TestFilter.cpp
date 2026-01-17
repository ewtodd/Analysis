#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TH1F.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TTree.h>

void Compare(std::vector<TString> filenames) {
  Int_t n_files = filenames.size();

  for (Int_t j = 0; j < n_files; j++) {
    TString filename = filenames.at(j);
    TString filepath = "root_files/" + filename + ".root";
    TFile *file = new TFile(filepath, "UPDATE");
    TTree *tree_with_pos = static_cast<TTree *>(file->Get("bef_tree"));

    Float_t energy = 0;
    Float_t x = 0, y = 0, z = 0;

    tree_with_pos->SetBranchAddress("energykeV", &energy);
    tree_with_pos->SetBranchAddress("xum", &x);
    tree_with_pos->SetBranchAddress("yum", &y);
    tree_with_pos->SetBranchAddress("zum", &z);

    TH2F *XvsE_zoomed =
        new TH2F(PlottingUtils::GetRandomName(),
                 "; Interaction Energy [keV]; Interaction X Position [um]",
                 Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN,
                 Constants::ZOOMED_XMAX, 20, -215, 215);
    TH2F *YvsE_zoomed =
        new TH2F(PlottingUtils::GetRandomName(),
                 "; Interaction Energy [keV]; Interaction Y Position [um]",
                 Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN,
                 Constants::ZOOMED_XMAX, 20, -215, 215);
    TH2F *ZvsE_zoomed =
        new TH2F(PlottingUtils::GetRandomName(),
                 "; Interaction Energy [keV]; Interaction Z Position [um]",
                 Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN,
                 Constants::ZOOMED_XMAX, 200, 0, 100);

    Int_t n_entries = tree_with_pos->GetEntries();

    TParameter<Double_t> *param =
        (TParameter<Double_t> *)file->Get("N42_RealTime_Total");

    if (!param) {
      std::cerr << "WARNING: Could not find N42_RealTime_Total parameter in "
                << filepath << std::endl;
    }
    TString perSecond = Constants::NORMALIZE_BY_TIME && param ? " / s" : "";

    TH1F *zoomedHistSeventyToOneHundred =
        new TH1F(PlottingUtils::GetRandomName(),
                 Form("%s; Energy [keV]; Counts / %d eV%s", filename.Data(),
                      Constants::BIN_WIDTH_EV, perSecond.Data()),
                 Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN,
                 Constants::ZOOMED_XMAX);

    for (Int_t i = 0; i < n_entries; i++) {
      tree_with_pos->GetEntry(i);

      if (energy > Constants::ZOOMED_XMIN && energy < Constants::ZOOMED_XMAX) {
        ZvsE_zoomed->Fill(energy, z);
        zoomedHistSeventyToOneHundred->Fill(energy);
        XvsE_zoomed->Fill(energy, x);
        YvsE_zoomed->Fill(energy, y);
      }
    }

    std::cout << "Created histograms for " << filename << std::endl;

    TCanvas *canvasXvsE_zoomed = new TCanvas("", "", 1200, 800);
    PlottingUtils::ConfigureCanvas(canvasXvsE_zoomed);
    PlottingUtils::ConfigureAndDraw2DHistogram(XvsE_zoomed, canvasXvsE_zoomed);
    TLine *lineX = new TLine(68.75, -150, 68.75, 150);
    lineX->SetLineColor(kRed);
    lineX->SetLineWidth(2);
    lineX->SetLineStyle(2);
    lineX->Draw();
    PlottingUtils::SaveFigure(canvasXvsE_zoomed, filename + "_XvsE_zoomed.png",
                              kFALSE);

    TCanvas *canvasYvsE_zoomed = new TCanvas("", "", 1200, 800);
    PlottingUtils::ConfigureCanvas(canvasYvsE_zoomed);
    PlottingUtils::ConfigureAndDraw2DHistogram(YvsE_zoomed, canvasYvsE_zoomed);
    TLine *lineY = new TLine(68.75, -150, 68.75, 150);
    lineY->SetLineColor(kRed);
    lineY->SetLineWidth(2);
    lineY->SetLineStyle(2);
    lineY->Draw();
    PlottingUtils::SaveFigure(canvasYvsE_zoomed, filename + "_YvsE_zoomed.png",
                              kFALSE);

    TCanvas *canvasZvsE_zoomed = new TCanvas("", "", 1200, 800);
    PlottingUtils::ConfigureCanvas(canvasZvsE_zoomed);
    PlottingUtils::ConfigureAndDraw2DHistogram(ZvsE_zoomed, canvasZvsE_zoomed);
    TLine *lineZ = new TLine(68.75, 0, 68.75, 100);
    lineZ->SetLineColor(kRed);
    lineZ->SetLineWidth(2);
    lineZ->SetLineStyle(2);
    lineZ->Draw();
    PlottingUtils::SaveFigure(canvasZvsE_zoomed, filename + "_ZvsE_zoomed.png",
                              kFALSE);

    TCanvas *canvasSeventyToOneHundred = new TCanvas("", "", 1200, 800);
    PlottingUtils::ConfigureCanvas(canvasSeventyToOneHundred, kFALSE);
    PlottingUtils::ConfigureAndDrawHistogram(zoomedHistSeventyToOneHundred,
                                             kBlue);
    zoomedHistSeventyToOneHundred->SetFillStyle(0);
    zoomedHistSeventyToOneHundred->SetLineWidth(2);
    PlottingUtils::SaveFigure(canvasSeventyToOneHundred,
                              filename + "90_100.png", kFALSE);

    std::cout << "Wrote histograms for " << filename << std::endl;

    ZvsE_zoomed->Write("ZvsE_zoomed", TObject::kOverwrite);

    file->Close();
  }
}

void TestFilter() {
  InitUtils::SetROOTPreferences();

  TString suffix = Constants::FILTERED ? "_filtered" : "";

  // January 12, 2026
  TString PassiveBackground_01122026 = "01122026-PassiveBackground" + suffix;
  TString Calibration_01122026 = "01122026-Calibration" + suffix;

  // January 13, 2026
  TString ActiveBackground_5Percent_01132026 =
      "01132026-ActiveBackground-5Percent" + suffix;
  TString ActiveBackground_25Percent_01132026 =
      "01132026-ActiveBackground-25Percent" + suffix;
  TString ActiveBackground_90Percent_01132026 =
      "01132026-ActiveBackground-90Percent" + suffix;

  TString CdShieldSignal_10Percent =
      "01132026-CdShield-GeSamplesIn-10Percent" + suffix;
  TString CdShieldSignal_25Percent =
      "01132026-CdShield-GeSamplesIn-25Percent" + suffix;
  TString CdShieldBackground =
      "01132026-CdShield-ActiveBackground-10Percent" + suffix;

  TString CuShieldSignal_01132026 =
      "01132026-CuShield-GeSamplesIn-10Percent" + suffix;
  TString CuShieldBackground_01132026 =
      "01132026-CuShield-ActiveBackground-Am241-10Percent" + suffix;

  TString PostReactor_Calibration_01132026 =
      "01132026-PostReactor-Calibration" + suffix;

  // January 14, 2026
  TString CuShieldSignal_01142026 =
      "01142026-CuShield-GeSamplesIn-10Percent" + suffix;
  TString CuShieldBackground_01142026 =
      "01142026-CuShield-ActiveBackground-10Percent" + suffix;
  TString CuShieldSignal_MovedBack_01142026 =
      "01142026-CuShield-GeSamplesIn-MovedBack-90Percent" + suffix;

  // January 15, 2026
  TString NoShieldSignal_01152026 =
      "01152026-NewSetup-GeSamplesIn-5Percent" + suffix;
  TString NoShieldBackground_01152026 =
      "01152026-NewSetup-ActiveBackground-5Percent" + suffix;
  TString PostReactor_Am241_01152026 =
      "01152026-NewSetup-PostReactor-Am241" + suffix;
  TString PostReactor_Ba133_01152026 =
      "01152026-NewSetup-PostReactor-Ba133" + suffix;
  TString ShutterClosed_01152026 =
      "01152026-NewSetup-ShutterClosed-5Percent" + suffix;

  // January 16, 2026
  TString NoShield_GeOnCZT_01162026 =
      "01162026-NoShield-GeOnCZT-0_5Percent" + suffix;
  TString NoShield_ActiveBackground_01162026 =
      "01162026-NoShield-ActiveBackground-0_5Percent" + suffix;
  TString NoShield_GraphiteCastle_Signal_01162026 =
      "01162026-NoShield-GeSamplesIn-GraphiteCastle-10Percent" + suffix;
  TString NoShield_GraphiteCastle_Background_01162026 =
      "01162026-NoShield-ActiveBackground-GraphiteCastle-10Percent" + suffix;
  TString PostReactor_Am241_Ba133_01162026 =
      "01162026-NoShield-PostReactor-Am241-Ba133" + suffix;

  std::vector<TString> filenames;
  filenames.push_back(NoShield_GeOnCZT_01162026);
  filenames.push_back(NoShield_ActiveBackground_01162026);
  filenames.push_back(NoShield_GraphiteCastle_Signal_01162026);
  filenames.push_back(NoShield_GraphiteCastle_Background_01162026);
  filenames.push_back(PostReactor_Am241_Ba133_01162026);
  Compare(filenames);
}
