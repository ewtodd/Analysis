#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TTree.h>

void Filter(std::vector<TString> filenames) {
  Int_t n_files = filenames.size();

  for (Int_t j = 0; j < n_files; j++) {
    TString filename = filenames.at(j);
    TString filepath = "root_files/" + filename + ".root";
    TFile *file = new TFile(filepath, "UPDATE");
    TTree *tree_with_pos = static_cast<TTree *>(file->Get("bef_tree"));

    Float_t energy = 0;
    Float_t x = 0, y = 0, z = 0;

    Int_t nInteractions = 0;

    tree_with_pos->SetBranchAddress("energykeV", &energy);
    tree_with_pos->SetBranchAddress("xum", &x);
    tree_with_pos->SetBranchAddress("yum", &y);
    tree_with_pos->SetBranchAddress("zum", &z);
    if (!Constants::FILTERED)
      tree_with_pos->SetBranchAddress("nInteractions", &nInteractions);

    TH2F *XY = new TH2F(PlottingUtils::GetRandomName(),
                        "; X Position [um]; Y Position [um]", 22, -215, 215, 22,
                        -215, 215);

    TH2F *EZ = new TH2F(PlottingUtils::GetRandomName(),
                        "; Energy [keV]; Interaction Z Position [um]",
                        Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN,
                        Constants::ZOOMED_XMAX, 200, 0, 100);

    Int_t n_entries = tree_with_pos->GetEntries();

    TParameter<Double_t> *param =
        (TParameter<Double_t> *)file->Get("N42_RealTime_Total");

    if (!param) {
      std::cerr << "WARNING: Could not find N42_RealTime_Total parameter in "
                << filepath << std::endl;
    }

    TH1F *includedSpectrum =
        new TH1F(PlottingUtils::GetRandomName(),
                 Form("%s; Energy [keV]; Counts / %d eV", filename.Data(),
                      Constants::BIN_WIDTH_EV),
                 Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN,
                 Constants::ZOOMED_XMAX);

    TH1F *excludedSpectrum =
        new TH1F(PlottingUtils::GetRandomName(),
                 Form("%s; Energy [keV]; Counts / %d eV", filename.Data(),
                      Constants::BIN_WIDTH_EV),
                 Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN,
                 Constants::ZOOMED_XMAX);

    Bool_t in_excluded_region;

    for (Int_t i = 0; i < n_entries; i++) {
      in_excluded_region = kFALSE;

      tree_with_pos->GetEntry(i);

      if (!Constants::FILTERED) {
        if (nInteractions != 1)
          in_excluded_region = kTRUE;
      }

      if (z < Constants::FILTER_DEPTH_UM)
        in_excluded_region = kTRUE;

      if (energy > Constants::ZOOMED_XMIN && energy < Constants::ZOOMED_XMAX) {
        XY->Fill(x, y);
        EZ->Fill(energy, z);

        for (size_t r = 0; r < Constants::FILTER_REGIONS_EXCLUDE_XY_UM.size();
             r++) {
          if (x >= Constants::FILTER_REGIONS_EXCLUDE_XY_UM[r].xmin &&
              x <= Constants::FILTER_REGIONS_EXCLUDE_XY_UM[r].xmax &&
              y >= Constants::FILTER_REGIONS_EXCLUDE_XY_UM[r].ymin &&
              y <= Constants::FILTER_REGIONS_EXCLUDE_XY_UM[r].ymax) {
            in_excluded_region = kTRUE;
            break;
          }
        }
        if (in_excluded_region) {
          excludedSpectrum->Fill(energy);
        } else
          includedSpectrum->Fill(energy);
      }
    }

    std::cout << "Created histograms for " << filename << std::endl;

    TCanvas *canvasXY = new TCanvas("", "", 1200, 800);
    canvasXY->SetRightMargin(0.2);
    XY->Draw("COLZ");
    PlottingUtils::Configure2DHistogram(XY, canvasXY);
    canvasXY->SetLogz(kFALSE);
    canvasXY->Update();

    for (size_t r = 0; r < Constants::FILTER_REGIONS_EXCLUDE_XY_UM.size();
         r++) {
      TBox *box = new TBox(Constants::FILTER_REGIONS_EXCLUDE_XY_UM[r].xmin,
                           Constants::FILTER_REGIONS_EXCLUDE_XY_UM[r].ymin,
                           Constants::FILTER_REGIONS_EXCLUDE_XY_UM[r].xmax,
                           Constants::FILTER_REGIONS_EXCLUDE_XY_UM[r].ymax);
      box->SetLineColor(kBlack);
      box->SetLineWidth(3);
      box->SetLineStyle(2);
      box->SetFillColorAlpha(kBlack, 0.5);
      box->Draw();
    }

    PlottingUtils::SaveFigure(canvasXY, filename + "_YvsX.png", kFALSE);

    TCanvas *canvasEZ = new TCanvas("", "", 1200, 800);
    canvasEZ->SetRightMargin(0.2);
    EZ->Draw("COLZ");
    PlottingUtils::Configure2DHistogram(EZ, canvasEZ);
    canvasEZ->SetLogz(kFALSE);
    canvasEZ->Update();

    TBox *box = new TBox(Constants::ZOOMED_XMIN, 0, Constants::ZOOMED_XMAX,
                         Constants::FILTER_DEPTH_UM);
    box->SetLineColor(kBlack);
    box->SetLineWidth(3);
    box->SetLineStyle(2);
    box->SetFillColorAlpha(kBlack, 0.5);
    box->Draw();

    PlottingUtils::SaveFigure(canvasEZ, filename + "_ZvsE.png", kFALSE);

    TCanvas *canvasRegions = new TCanvas("", "", 1200, 800);
    PlottingUtils::ConfigureCanvas(canvasRegions, kFALSE);

    Double_t maxY = TMath::Max(includedSpectrum->GetMaximum(),
                               excludedSpectrum->GetMaximum());

    PlottingUtils::ConfigureHistogram(includedSpectrum, kRed);
    includedSpectrum->GetYaxis()->SetRangeUser(0, maxY * 1.1);
    includedSpectrum->SetFillStyle(0);
    includedSpectrum->SetLineWidth(2);
    includedSpectrum->GetYaxis()->SetTitleOffset(1.5);
    includedSpectrum->Draw("HIST");

    PlottingUtils::ConfigureHistogram(excludedSpectrum, kBlack);
    excludedSpectrum->SetFillStyle(0);
    excludedSpectrum->SetLineWidth(2);
    excludedSpectrum->Draw("HIST SAME");

    TLine *lineRegional = new TLine(68.75, 0, 68.75, maxY * 1.1);
    lineRegional->SetLineColor(kRed);
    lineRegional->SetLineWidth(2);
    lineRegional->SetLineStyle(2);
    lineRegional->Draw();

    TLegend *leg = new TLegend(0.78, 0.75, 0.9, 0.88);
    leg->AddEntry(includedSpectrum, "Accepted", "l");
    leg->AddEntry(excludedSpectrum, "Rejected", "l");
    leg->Draw();

    canvasRegions->SetLeftMargin(0.2);
    PlottingUtils::SaveFigure(canvasRegions, filename + "_FilterSpectra.png",
                              kFALSE);
    std::cout << "Wrote histograms for " << filename << std::endl;

    XY->Write("XY", TObject::kOverwrite);
    EZ->Write("EZ", TObject::kOverwrite);

    canvasRegions->cd();
    canvasRegions->Write("AcceptedRejectedSpectra", TObject::kOverwrite);

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

  // January 12, 2026
  //  filenames.push_back(PassiveBackground_01122026);
  //  filenames.push_back(Calibration_01122026);
  //
  //  //  // January 13, 2026
  //  filenames.push_back(ActiveBackground_5Percent_01132026);
  //  filenames.push_back(ActiveBackground_25Percent_01132026);
  //  filenames.push_back(ActiveBackground_90Percent_01132026);
  //  filenames.push_back(CdShieldSignal_10Percent);
  filenames.push_back(CdShieldSignal_25Percent);
  //  filenames.push_back(CdShieldBackground);
  //  filenames.push_back(CuShieldSignal_01132026);
  //  filenames.push_back(CuShieldBackground_01132026);
  //  filenames.push_back(PostReactor_Calibration_01132026);
  //
  //  // January 14, 2026
  //  filenames.push_back(CuShieldSignal_01142026);
  //  filenames.push_back(CuShieldBackground_01142026);
  //  filenames.push_back(CuShieldSignal_MovedBack_01142026);
  //
  //  // January 15, 2026
  //  filenames.push_back(NoShieldSignal_01152026);
  //  filenames.push_back(NoShieldBackground_01152026);
  //  filenames.push_back(PostReactor_Am241_01152026);
  //  filenames.push_back(PostReactor_Ba133_01152026);
  //  filenames.push_back(ShutterClosed_01152026);
  //
  //  // January 16, 2026
  //  filenames.push_back(NoShield_GeOnCZT_01162026);
  //  filenames.push_back(NoShield_ActiveBackground_01162026);
  //  filenames.push_back(NoShield_GraphiteCastle_Signal_01162026);
  //  filenames.push_back(NoShield_GraphiteCastle_Background_01162026);
  //  filenames.push_back(PostReactor_Am241_Ba133_01162026);

  Filter(filenames);
}
