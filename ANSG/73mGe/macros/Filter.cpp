#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TROOT.h>
#include <TTree.h>

void MapFiles(std::vector<TString> filenames) {
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
    tree_with_pos->SetBranchAddress("nInteractions", &nInteractions);

    TH2F *XvsE =
        new TH2F(PlottingUtils::GetRandomName(),
                 "; Interaction Energy [keV]; Interaction X Position [um]",
                 Constants::HIST_NBINS, Constants::HIST_XMIN,
                 Constants::HIST_XMAX, 300, -150, 150);
    TH2F *YvsE =
        new TH2F(PlottingUtils::GetRandomName(),
                 "; Interaction Energy [keV]; Interaction Y Position [um]",
                 Constants::HIST_NBINS, Constants::HIST_XMIN,
                 Constants::HIST_XMAX, 300, -150, 150);
    TH2F *ZvsE =
        new TH2F(PlottingUtils::GetRandomName(),
                 "; Interaction Energy [keV]; Interaction Z Position [um]",
                 Constants::HIST_NBINS, Constants::HIST_XMIN,
                 Constants::HIST_XMAX, 200, 0, 100);

    // Zoomed histograms
    TH2F *XvsE_zoomed =
        new TH2F(PlottingUtils::GetRandomName(),
                 "; Interaction Energy [keV]; Interaction X Position [um]",
                 Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN,
                 Constants::ZOOMED_XMAX, 300, -150, 150);
    TH2F *YvsE_zoomed =
        new TH2F(PlottingUtils::GetRandomName(),
                 "; Interaction Energy [keV]; Interaction Y Position [um]",
                 Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN,
                 Constants::ZOOMED_XMAX, 300, -150, 150);
    TH2F *ZvsE_zoomed =
        new TH2F(PlottingUtils::GetRandomName(),
                 "; Interaction Energy [keV]; Interaction Z Position [um]",
                 Constants::ZOOMED_NBINS, Constants::ZOOMED_XMIN,
                 Constants::ZOOMED_XMAX, 200, 0, 100);

    Int_t n_entries = tree_with_pos->GetEntries();
    for (Int_t i = 0; i < n_entries; i++) {
      tree_with_pos->GetEntry(i);
      if (nInteractions != 1)
        continue;

      XvsE->Fill(energy, x);
      YvsE->Fill(energy, y);
      ZvsE->Fill(energy, z);
      if (energy > Constants::ZOOMED_XMIN && energy < Constants::ZOOMED_XMAX) {
        XvsE_zoomed->Fill(energy, x);
        YvsE_zoomed->Fill(energy, y);
        ZvsE_zoomed->Fill(energy, z);
      }
    }

    std::cout << "Created histograms for " << filename << std::endl;

    TCanvas *canvasXvsE = new TCanvas("", "", 1200, 800);
    PlottingUtils::ConfigureCanvas(canvasXvsE);
    PlottingUtils::ConfigureAndDraw2DHistogram(XvsE, canvasXvsE);

    TCanvas *canvasYvsE = new TCanvas("", "", 1200, 800);
    PlottingUtils::ConfigureCanvas(canvasYvsE);
    PlottingUtils::ConfigureAndDraw2DHistogram(YvsE, canvasYvsE);

    TCanvas *canvasZvsE = new TCanvas("", "", 1200, 800);
    PlottingUtils::ConfigureCanvas(canvasZvsE);
    PlottingUtils::ConfigureAndDraw2DHistogram(ZvsE, canvasZvsE);

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

    std::cout << "Wrote histograms for " << filename << std::endl;

    XvsE->Write("XvsE", TObject::kOverwrite);
    YvsE->Write("YvsE", TObject::kOverwrite);
    ZvsE->Write("ZvsE", TObject::kOverwrite);

    XvsE_zoomed->Write("XvsE_zoomed", TObject::kOverwrite);
    YvsE_zoomed->Write("YvsE_zoomed", TObject::kOverwrite);
    ZvsE_zoomed->Write("ZvsE_zoomed", TObject::kOverwrite);

    file->Close();
  }
}

void TripleMapFiles(std::vector<TString> filenames) {
  Int_t n_files = filenames.size();

  for (Int_t j = 0; j < n_files; j++) {
    TString filename = filenames.at(j);
    TString filepath = "root_files/" + filename + ".root";
    TFile *file = new TFile(filepath, "UPDATE");
    TTree *tree_with_pos = static_cast<TTree *>(file->Get("bef_tree"));

    Float_t energykeV = 0;
    Float_t xum = 0, yum = 0, zum = 0;

    tree_with_pos->SetBranchAddress("energykeV", &energykeV);
    tree_with_pos->SetBranchAddress("xum", &xum);
    tree_with_pos->SetBranchAddress("yum", &yum);
    tree_with_pos->SetBranchAddress("zum", &zum);

    TH3F *XYvsE =
        new TH3F(PlottingUtils::GetRandomName(),
                 "; X Position [um]; Y Position [um]; Energy [keV]", 300, -150,
                 150, 300, -150, 150, Constants::HIST_NBINS,
                 Constants::HIST_XMIN, Constants::HIST_XMAX);

    TH3F *XYZ = new TH3F(PlottingUtils::GetRandomName(),
                         "; X Position [um]; Y Position [um]; Z Position [um]",
                         300, -150, 150, 300, -150, 150, 200, 0, 100);

    Int_t n_entries = 1e5; // tree_with_pos->GetEntries();
    for (Int_t i = 0; i < n_entries; i++) {
      tree_with_pos->GetEntry(i);
      XYvsE->Fill(xum, yum, energykeV);
      XYZ->Fill(xum, yum, zum);
    }

    std::cout << "Created histograms for " << filename << std::endl;

    TCanvas *canvasXYvsE = new TCanvas("", "", 1200, 800);
    PlottingUtils::ConfigureCanvas(canvasXYvsE);
    canvasXYvsE->SetRightMargin(0.15);
    gStyle->SetPalette(kRainBow);
    XYvsE->Draw("BOX2Z");
    PlottingUtils::SaveFigure(canvasXYvsE, filename + "_XYvsE.png", kFALSE);

    TCanvas *canvasXYZ = new TCanvas("", "", 1200, 800);
    PlottingUtils::ConfigureCanvas(canvasXYZ);
    canvasXYZ->SetRightMargin(0.15);
    gStyle->SetPalette(kRainBow);
    XYZ->Draw("BOX2Z");
    PlottingUtils::SaveFigure(canvasXYZ, filename + "_XYZ.png", kFALSE);

    std::cout << "Wrote histograms for " << filename << std::endl;

    XYvsE->Write("XYvsE", TObject::kOverwrite);
    XYZ->Write("XYZ", TObject::kOverwrite);
    file->Close();
  }
}

void Map() {
  InitUtils::SetROOTPreferences();
  TString CdShieldSignal = "01132026-CdShield-GeSamplesIn-10Percent";
  TString CdShieldBackground = "01132026-CdShield-ActiveBackground-10Percent";
  TString CuShieldSignal_01132026 = "01132026-CuShield-GeSamplesIn-10Percent";
  TString CuShieldBackground_01132026 =
      "01132026-CuShield-ActiveBackground-Am241-10Percent";
  TString CuShieldSignal_01142026 = "01142026-CuShield-GeSamplesIn-10Percent";
  TString CuShieldBackground_01142026 =
      "01142026-CuShield-ActiveBackground-10Percent";
  TString NoShieldSignal_01152026 = "01152026-NewSetup-GeSamplesIn-5Percent";
  TString NoShieldBackground_01152026 =
      "01152026-NewSetup-ActiveBackground-5Percent";

  std::vector<TString> filenames;
  // filenames.push_back(CdShieldSignal);
  // filenames.push_back(CdShieldBackground);
  filenames.push_back(CuShieldSignal_01132026);
  // filenames.push_back(CuShieldBackground_01132026);
  // filenames.push_back(CuShieldSignal_01142026);
  // filenames.push_back(CuShieldBackground_01142026);
  // filenames.push_back(NoShieldSignal_01152026);
  // filenames.push_back(NoShieldBackground_01152026);

  MapFiles(filenames);
  TripleMapFiles(filenames);
}
