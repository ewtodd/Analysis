#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TTree.h>

void MapBandedEvents(std::vector<TString> filenames) {
  Int_t n_files = filenames.size();

  std::vector<Double_t> bandCenters = {17.25, 36.25,  55.25, 74.25,
                                       93.25, 112.75, 131.75};

  const Double_t bandWidth = 0.35;

  for (Int_t j = 0; j < n_files; j++) {
    TString filename = filenames.at(j);
    TString filepath = "root_files/" + filename + ".root";
    TFile *file = new TFile(filepath, "UPDATE");
    TTree *tree_with_pos = static_cast<TTree *>(file->Get("bef_tree_with_pos"));

    UInt_t energy = 0;
    UInt_t totalEnergy = 0;
    Short_t x = 0, y = 0, z = 0;
    Int_t interaction = 0;
    Int_t nInteractions = 0;
    Char_t crystalNumber = 0;
    Short_t x0 = 0, y0 = 0, z0 = 0;

    tree_with_pos->SetBranchAddress("energy", &energy);
    tree_with_pos->SetBranchAddress("x", &x);
    tree_with_pos->SetBranchAddress("y", &y);
    tree_with_pos->SetBranchAddress("z", &z);
    tree_with_pos->SetBranchAddress("nInteractions", &nInteractions);
    tree_with_pos->SetBranchAddress("interaction", &interaction);

    double binWidth_keV = Constants::BIN_WIDTH_EV / 1000.0;

    int hist_xmin = 0, hist_xmax = 1500;
    int zoom_xmin = 60, zoom_xmax = 80;
    int hist_nbins = (hist_xmax - hist_xmin) / binWidth_keV;
    int zoomed_nbins = (zoom_xmax - zoom_xmin) / binWidth_keV;

    TH1F *EnergySpectrum_BothBands = new TH1F(
        PlottingUtils::GetRandomName(),
        Form("; Energy [keV]; Counts / %d eV", Constants::BIN_WIDTH_EV),
        hist_nbins, hist_xmin, hist_xmax);

    TH1F *EnergySpectrum_BothBands_Zoomed = new TH1F(
        PlottingUtils::GetRandomName(),
        Form("; Energy [keV]; Counts / %d eV", Constants::BIN_WIDTH_EV),
        zoomed_nbins, zoom_xmin, zoom_xmax);

    TH2F *XvsE_XBands =
        new TH2F(PlottingUtils::GetRandomName(),
                 "; Interaction Energy [keV]; Interaction X Position [um]",
                 hist_nbins, hist_xmin, hist_xmax, 600, -150, 150);

    TH1F *EnergySpectrum_BothNotBands = new TH1F(
        PlottingUtils::GetRandomName(),
        Form("; Energy [keV]; Counts / %d eV", Constants::BIN_WIDTH_EV),
        hist_nbins, hist_xmin, hist_xmax);

    TH1F *EnergySpectrum_BothNotBands_Zoomed = new TH1F(
        PlottingUtils::GetRandomName(),
        Form("; Energy [keV]; Counts / %d eV", Constants::BIN_WIDTH_EV),
        zoomed_nbins, zoom_xmin, zoom_xmax);

    TH2F *XvsE_XNotBands =
        new TH2F(PlottingUtils::GetRandomName(),
                 "; Interaction Energy [keV]; Interaction X Position [um]",
                 hist_nbins, hist_xmin, hist_xmax, 600, -150, 150);

    Float_t energykeV = 0.0;
    Float_t xum = 0.0;
    Float_t yum = 0.0;
    Float_t zum = 0.0;

    Int_t n_entries = tree_with_pos->GetEntries();
    for (Int_t i = 0; i < n_entries; i++) {
      tree_with_pos->GetEntry(i);

      if (nInteractions != 1 || interaction != 0)
        continue;

      energykeV = energy / 1000.0;
      xum = x / 10.0;
      yum = y / 10.0;
      zum = z / 10.0;

      Bool_t xInBand = kFALSE;
      for (const auto &center : bandCenters) {
        if (TMath::Abs(xum - center) < bandWidth ||
            TMath::Abs(xum + center) < bandWidth) {
          xInBand = kTRUE;
          break;
        }
      }

      Bool_t yInBand = kFALSE;
      for (const auto &center : bandCenters) {
        if (TMath::Abs(yum - center) < bandWidth ||
            TMath::Abs(yum + center) < bandWidth) {
          yInBand = kTRUE;
          break;
        }
      }

      if (xInBand && yInBand) {
        EnergySpectrum_BothBands->Fill(energykeV);
        EnergySpectrum_BothBands_Zoomed->Fill(energykeV);
      }

      if (xInBand) {
        XvsE_XBands->Fill(energykeV, xum);
      }

      if (!xInBand && !yInBand) {
        EnergySpectrum_BothNotBands->Fill(energykeV);
        EnergySpectrum_BothNotBands_Zoomed->Fill(energykeV);
      }

      if (!xInBand) {
        XvsE_XNotBands->Fill(energykeV, xum);
      }
    }

    std::cout << "Created banded energy spectra for " << filename << std::endl;

    TCanvas *canvasBothBands = new TCanvas("", "", 1200, 800);
    PlottingUtils::ConfigureCanvas(canvasBothBands);
    PlottingUtils::ConfigureAndDrawHistogram(EnergySpectrum_BothBands, kPink);
    PlottingUtils::SaveFigure(
        canvasBothBands, filename + "_EnergySpectrum_BothBands.png", kFALSE);

    TCanvas *canvasBothBands_Zoomed = new TCanvas("", "", 1200, 800);
    PlottingUtils::ConfigureCanvas(canvasBothBands_Zoomed);
    PlottingUtils::ConfigureAndDrawHistogram(EnergySpectrum_BothBands_Zoomed,
                                             kPink);
    PlottingUtils::SaveFigure(canvasBothBands_Zoomed,
                              filename + "_EnergySpectrum_BothBands_Zoomed.png",
                              kFALSE);

    TCanvas *canvasXvsE_XBands = new TCanvas("", "", 1200, 800);
    PlottingUtils::ConfigureCanvas(canvasXvsE_XBands);
    PlottingUtils::ConfigureAndDraw2DHistogram(XvsE_XBands, canvasXvsE_XBands);
    PlottingUtils::SaveFigure(canvasXvsE_XBands, filename + "_XvsE_XBands.png",
                              kFALSE);

    TCanvas *canvasBothNotBands = new TCanvas("", "", 1200, 800);
    PlottingUtils::ConfigureCanvas(canvasBothNotBands);
    PlottingUtils::ConfigureAndDrawHistogram(EnergySpectrum_BothNotBands,
                                             kPink);
    PlottingUtils::SaveFigure(canvasBothNotBands,
                              filename + "_EnergySpectrum_BothNotBands.png",
                              kFALSE);

    TCanvas *canvasBothNotBands_Zoomed = new TCanvas("", "", 1200, 800);
    PlottingUtils::ConfigureCanvas(canvasBothNotBands_Zoomed);
    PlottingUtils::ConfigureAndDrawHistogram(EnergySpectrum_BothNotBands_Zoomed,
                                             kPink);
    PlottingUtils::SaveFigure(
        canvasBothNotBands_Zoomed,
        filename + "_EnergySpectrum_BothNotBands_Zoomed.png", kFALSE);

    TCanvas *canvasXvsE_XNotBands = new TCanvas("", "", 1200, 800);
    PlottingUtils::ConfigureCanvas(canvasXvsE_XNotBands);
    PlottingUtils::ConfigureAndDraw2DHistogram(XvsE_XNotBands,
                                               canvasXvsE_XNotBands);
    PlottingUtils::SaveFigure(canvasXvsE_XNotBands,
                              filename + "_XvsE_XNotBands.png", kFALSE);

    std::cout << "Wrote banded energy spectra for " << filename << std::endl;

    EnergySpectrum_BothBands->Write("EnergySpectrum_BothBands",
                                    TObject::kOverwrite);
    EnergySpectrum_BothBands_Zoomed->Write("EnergySpectrum_BothBands_Zoomed",
                                           TObject::kOverwrite);
    XvsE_XBands->Write("XvsE_XBands", TObject::kOverwrite);
    EnergySpectrum_BothNotBands->Write("EnergySpectrum_BothNotBands",
                                       TObject::kOverwrite);
    EnergySpectrum_BothNotBands_Zoomed->Write(
        "EnergySpectrum_BothNotBands_Zoomed", TObject::kOverwrite);
    XvsE_XNotBands->Write("XvsE_XNotBands", TObject::kOverwrite);

    file->Close();
  }
}

void MapFiles(std::vector<TString> filenames) {
  Int_t n_files = filenames.size();

  for (Int_t j = 0; j < n_files; j++) {
    TString filename = filenames.at(j);
    TString filepath = "root_files/" + filename + ".root";
    TFile *file = new TFile(filepath, "UPDATE");
    TTree *tree_with_pos = static_cast<TTree *>(file->Get("bef_tree_with_pos"));

    UInt_t energy = 0;
    UInt_t totalEnergy = 0;
    UInt_t eventTime = 0;
    Int_t liveTime = 0;
    Short_t x = 0, y = 0, z = 0;
    Int_t interaction = 0;
    Int_t nInteractions = 0;
    Char_t crystalNumber = 0;
    Short_t x0 = 0, y0 = 0, z0 = 0;

    tree_with_pos->SetBranchAddress("energy", &energy);
    tree_with_pos->SetBranchAddress("eventTime", &eventTime);
    tree_with_pos->SetBranchAddress("liveTime", &liveTime);
    tree_with_pos->SetBranchAddress("x", &x);
    tree_with_pos->SetBranchAddress("y", &y);
    tree_with_pos->SetBranchAddress("z", &z);
    tree_with_pos->SetBranchAddress("nInteractions", &nInteractions);
    tree_with_pos->SetBranchAddress("interaction", &interaction);

    double binWidth_keV = Constants::BIN_WIDTH_EV / 1000.0;

    int hist_xmin = 0, hist_xmax = 1500;
    int hist_nbins = (hist_xmax - hist_xmin) / binWidth_keV;

    TH2F *XvsE =
        new TH2F(PlottingUtils::GetRandomName(),
                 "; Interaction Energy [keV]; Interaction X Position [um]",
                 hist_nbins, hist_xmin, hist_xmax, 600, -150, 150);
    TH2F *YvsE =
        new TH2F(PlottingUtils::GetRandomName(),
                 "; Interaction Energy [keV]; Interaction Y Position [um]",
                 hist_nbins, hist_xmin, hist_xmax, 600, -150, 150);
    TH2F *ZvsE =
        new TH2F(PlottingUtils::GetRandomName(),
                 "; Interaction Energy [keV]; Interaction Z Position [um]",
                 hist_nbins, hist_xmin, hist_xmax, 500, 0, 100);

    Float_t energykeV = 0.0;
    Float_t xum = 0.0;
    Float_t yum = 0.0;
    Float_t zum = 0.0;

    Int_t n_entries = tree_with_pos->GetEntries();
    for (Int_t i = 0; i < n_entries; i++) {
      tree_with_pos->GetEntry(i);
      energykeV = energy / 1000.0;
      xum = x / 10.0;
      yum = y / 10.0;
      zum = z / 10.0;

      XvsE->Fill(energykeV, xum);
      YvsE->Fill(energykeV, yum);
      ZvsE->Fill(energykeV, zum);
    }

    std::cout << "Created histograms for " << filename << std::endl;

    TCanvas *canvasXvsE = new TCanvas("", "", 1200, 800);
    PlottingUtils::ConfigureCanvas(canvasXvsE);
    PlottingUtils::ConfigureAndDraw2DHistogram(XvsE, canvasXvsE);
    PlottingUtils::SaveFigure(canvasXvsE, filename + "_XvsE.png", kFALSE);

    TCanvas *canvasYvsE = new TCanvas("", "", 1200, 800);
    PlottingUtils::ConfigureCanvas(canvasYvsE);
    PlottingUtils::ConfigureAndDraw2DHistogram(YvsE, canvasYvsE);
    PlottingUtils::SaveFigure(canvasYvsE, filename + "_YvsE.png", kFALSE);

    TCanvas *canvasZvsE = new TCanvas("", "", 1200, 800);
    PlottingUtils::ConfigureCanvas(canvasZvsE);
    PlottingUtils::ConfigureAndDraw2DHistogram(ZvsE, canvasZvsE);
    PlottingUtils::SaveFigure(canvasZvsE, filename + "_ZvsE.png", kFALSE);

    std::cout << "Wrote histograms for " << filename << std::endl;

    XvsE->Write("XvsE", TObject::kOverwrite);

    YvsE->Write("YvsE", TObject::kOverwrite);
    ZvsE->Write("ZvsE", TObject::kOverwrite);

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

  std::vector<TString> filenames;
  filenames.push_back(CdShieldSignal);
  filenames.push_back(CdShieldBackground);
  filenames.push_back(CuShieldSignal_01132026);
  filenames.push_back(CuShieldBackground_01132026);
  filenames.push_back(CuShieldSignal_01142026);
  filenames.push_back(CuShieldBackground_01142026);

  MapBandedEvents(filenames);
  MapFiles(filenames);
}
