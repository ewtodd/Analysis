#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TH1I.h>

void Plot() {
  InitUtils::SetROOTPreferences(kTRUE);

  int binWidth_eV = 150;
  double binWidth_keV = binWidth_eV / 1000.0;

  int hist_xmin = 0, hist_xmax = 750;  // Range in keV
  int zoom_xmin = 30, zoom_xmax = 100; // Range in keV

  int hist_nbins = (hist_xmax - hist_xmin) / binWidth_keV;
  int zoomed_nbins = (zoom_xmax - zoom_xmin) / binWidth_keV;

  TString filepath = "./root_files/Ba133_TripleCs137_Lead_Gold_Overnight.root";
  TFile *input = new TFile(filepath, "READ");
  TTree *tree = static_cast<TTree *>(input->Get("bef_tree_with_pos"));

  UInt_t totalEnergy;
  tree->SetBranchAddress("totalEnergy", &totalEnergy);

  Int_t n_entries = tree->GetEntries();

  TCanvas *canvas = new TCanvas(PlottingUtils::GetRandomName(),
                                PlottingUtils::GetRandomName(), 1200, 800);

  TH1I *hist = new TH1I(PlottingUtils::GetRandomName(),
                        Form("; Energy [keV]; Counts / %d eV", binWidth_eV),
                        hist_nbins, hist_xmin, hist_xmax);

  TH1I *zoomedHist =
      new TH1I(PlottingUtils::GetRandomName(),
               Form("; Energy [keV]; Counts / %d eV", binWidth_eV),
               zoomed_nbins, zoom_xmin, zoom_xmax);

  Float_t energyKeV;
  for (Int_t i = 0; i < n_entries; i++) {
    tree->GetEntry(i);
    if (i % 50000000 == 0)
      std::cout << i << std::endl;

    energyKeV = totalEnergy / 1000.0;
    hist->Fill(energyKeV);
    if (energyKeV >= zoom_xmin && energyKeV < zoom_xmax)
      zoomedHist->Fill(energyKeV);
  }

  PlottingUtils::ConfigureAndDrawHistogram(hist, kRed);
  gPad->RedrawAxis("g");
  PlottingUtils::SaveFigure(canvas,
                            "Ba133_TripleCs137_Lead_Gold_Overnight.png");

  canvas->Clear();
  canvas->SetLogy(kFALSE);
  PlottingUtils::ConfigureAndDrawHistogram(zoomedHist, kBlue);
  PlottingUtils::SaveFigure(canvas,
                            "Ba133_TripleCs137_Lead_Gold_Overnight_zoomed.png");
  delete canvas;
};
